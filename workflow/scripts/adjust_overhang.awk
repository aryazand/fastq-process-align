# ═══════════════════════════════════════════════════════════════════════════════
# adjust_overhang.awk
# ═══════════════════════════════════════════════════════════════════════════════
#
# PURPOSE:
#   Corrects SAM alignment coordinates for reads mapping to the overhang region
#   of circularized chromosomes. When a genome is circularized for alignment,
#   the first N bases are appended to the 3' end as an overhang to allow reads
#   spanning the circular origin to align properly. This script post-processes
#   alignments to account for this artificial overhang.
#
# BEHAVIOR:
#   1. Reads aligning entirely within the overhang region are "wrapped" back to
#      their true chromosomal coordinates at the start of the sequence
#   2. Reads spanning the junction (normal region + overhang) are flagged as
#      unmapped because the alignment is ambiguous
#   3. Reads in the normal region are passed through unchanged
#   4. Prints per-chromosome statistics to stderr for quality control
#
# USAGE (via gawk):
#   gawk -f adjust_overhang.awk \
#     -v chrom_lengths="chr1:1000,chr2:2000" \
#     input.sam > output.sam

BEGIN {
    FS = OFS = "\t"
    
    # ─── Initialize chromosome lookup tables ─────────────────────────────────
    # Parse the comma-separated chrom_lengths variable (format: "chr:len,chr:len,...")
    # into two associative arrays:
    #   - chrom_length[chrom]: original chromosome length (without overhang)
    #   - is_circular[chrom]: flag indicating this chromosome is circular
    #
    # Example: chrom_lengths="NC_018936.1:1837281"
    #   → chrom_length["NC_018936.1"] = 1837281
    #   → is_circular["NC_018936.1"] = 1
    
    n = split(chrom_lengths, lens, ",")
    for (i = 1; i <= n; i++) {
        split(lens[i], kv, ":")
        chrom_length[kv[1]] = kv[2]
        is_circular[kv[1]]  = 1
    }
}

# ─── FILTER 1: Skip unmapped reads ──────────────────────────────────────────
# Unmapped reads are identified by:
#   - flag & 4: unmapped bit is set (bit 3, decimal value 4)
#   - pos == 0: no position assigned
#   - cigar == "*": no CIGAR string (SAM spec for unmapped)
# These reads pass through unchanged; just count them and move to next record
and($2, 4) || $4 == 0 || $6 == "*" {
    seen_chroms[$3]++
    unmapped[$3]++
    print
    next
}

# ─── FILTER 2: Skip non-circular chromosomes ───────────────────────────────
# Reads on non-circular chromosomes need no correction. Count and pass through.
!($3 in is_circular) {
    seen_chroms[$3]++
    normal[$3]++
    print
    next
}

# ─── PROCESS: Circular chromosome alignments ───────────────────────────────
# Only alignments on circular chromosomes that are mapped reach this block.
# Classify and correct the read based on its position relative to the overhang.
{
    seen_chroms[$3] = 1
    rname = $3
    pos = $4
    cigar = $6
    
    c_len          = chrom_length[rname]  # Original chromosome length
    overhang_start = c_len + 1             # First position of overhang (1-based)

    # ─── Calculate reference span of alignment from CIGAR string ─────────────
    # Only operators that consume reference sequence contribute to length:
    #   M (match)  D (delete)  N (skip)  X (mismatch)  = (match)
    # Note: I (insert), S (soft clip), H (hard clip) do NOT consume reference
    
    ref_len = 0
    tmp = cigar
    while (match(tmp, /[0-9]+[MDNX=]/)) {
        op_str = substr(tmp, RSTART, RLENGTH)
        ref_len += substr(op_str, 1, length(op_str) - 1)
        tmp = substr(tmp, RSTART + RLENGTH)
    }

    # The rightmost reference position covered by this alignment (1-based, inclusive)
    read_end = pos + ref_len - 1

    # ─── CASE 1: Read entirely within overhang ──────────────────────────────
    if (pos >= overhang_start) {
        overhang_counts[rname]++
        $4 = pos - c_len
        
        # If mate is on same chromosome and also in overhang, wrap mate position too
        if ($8 >= overhang_start) { 
            $8 = $8 - c_len 
        }

    # ─── CASE 2: Read spans the junction ────────────────────────────────────
    # This alignment is unreliable so flag as unmapped
    } else if (read_end > c_len) {
        junction[rname]++
        $2 = or($2, 4)
        $4 = 0
        $6 = "*"

    # ─── CASE 3: Read entirely within normal region ─────────────────────────
    # No correction needed
    } else {
        normal[rname]++
    }
    
    print
}

# ─── Generate summary statistics ─────────────────────────────────────────────
END {
    col1 = 16
    col2 = 10

    printf "\n%-*s %*s %*s %*s %*s %*s\n",
        col1, "Chromosome",
        col2, "Normal",
        col2, "Overhang",
        col2, "Junction",
        col2, "Unmapped",
        col2, "Total" > "/dev/stderr"

    divider = ""
    for (i = 1; i <= col1 + (col2 + 1) * 5; i++) divider = divider "-"
    print divider > "/dev/stderr"

    total_normal = 0
    total_overhang = 0
    total_junction = 0
    total_unmapped = 0

    for (chr in seen_chroms) {
        n_val  = (chr in normal)          ? normal[chr]          : 0
        ov_val = (chr in overhang_counts) ? overhang_counts[chr] : 0
        jn_val = (chr in junction)        ? junction[chr]        : 0
        un_val = (chr in unmapped)        ? unmapped[chr]        : 0
        tot    = n_val + ov_val + jn_val + un_val

        label = (chr in is_circular) ? chr "*" : chr

        printf "%-*s %*d %*d %*d %*d %*d\n",
            col1, label,
            col2, n_val,
            col2, ov_val,
            col2, jn_val,
            col2, un_val,
            col2, tot > "/dev/stderr"

        total_normal   += n_val
        total_overhang += ov_val
        total_junction += jn_val
        total_unmapped += un_val
    }

    print divider > "/dev/stderr"
    printf "%-*s %*d %*d %*d %*d %*d\n",
        col1, "TOTAL",
        col2, total_normal,
        col2, total_overhang,
        col2, total_junction,
        col2, total_unmapped,
        col2, total_normal + total_overhang + total_junction + total_unmapped > "/dev/stderr"
    
    printf "\n* = circular chromosome (overhang correction applied)\n\n" > "/dev/stderr"
}