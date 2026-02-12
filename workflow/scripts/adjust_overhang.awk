# adjust_overhang.awk
#
# Corrects SAM coordinates for reads aligning to the overhang region of a
# circularised genome. The overhang is a duplication of the first N bases
# appended to the 3' end of the genome to allow reads spanning the origin
# to align. This script:
#   1. Fixes the @SQ header LN: values to remove the overhang length
#   2. Wraps reads aligning entirely within the overhang back to their
#      true genomic coordinates
#   3. Flags reads spanning the junction as unmapped (ambiguous alignment)
#   4. Prints per-chromosome read classification statistics to stderr
#
# Usage:
#   awk -f adjust_overhang.awk \
#       -v chrom_lengths="chr1:10000,chr2:20000" \
#       -v overhang_len=500                      \
#       input.sam > output.sam

BEGIN {
    FS = OFS = "\t"  # Input and output are tab-delimited
    
    # ── Parse the chrom:length pairs ─────────────────────────────────────────
    # Split the comma-separated chrom:length pairs into the chrom_length
    # associative array, keyed by chromosome name. Simultaneously populate
    # is_circular from the same keys, replacing the old circular_chroms variable.
    # e.g. chrom_lengths="chr1:10000,chr2:20000"
    #   → chrom_length["chr1"] = 10000, is_circular["chr1"] = 1
    #   → chrom_length["chr2"] = 20000, is_circular["chr2"] = 1
    n = split(chrom_lengths, lens, ",")
    for (i = 1; i <= n; i++) {
        # Split each "chrom:length" pair on ":" into a temporary key-value array
        split(lens[i], kv, ":")
        chrom_length[kv[1]] = kv[2]
        is_circular[kv[1]]  = 1
    }
}

# ── @SQ header lines ─────────────────────────────────────────────────────────
# For circular chromosomes, correct the LN: tag to reflect the original genome
# length without the overhang. All other header lines are passed through unchanged.
/^@SQ/ {
    # Loop over all fields to find the SN: tag, since SAM does not guarantee
    # field order within a header line
    sn = ""
    for (i = 1; i <= NF; i++) {
        if ($i ~ /^SN:/) { sn = substr($i, 4) }  # substr skips the "SN:" prefix
    }

    # Only modify LN: for chromosomes listed as circular
    if (sn in is_circular) {
        out = ""
        for (i = 1; i <= NF; i++) {
            if ($i ~ /^LN:/) {
                # Replace LN: with the original genome length (without overhang)
                out = out (i > 1 ? "\t" : "") "LN:" chrom_length[sn]
            } else {
                # All other tags (SN:, M5:, etc.) are kept as-is
                out = out (i > 1 ? "\t" : "") $i
            }
        }
        print out
    } else {
        print $0
    }
    next  # Skip to next line — don't process header lines as alignments
}

# ── All other header lines (@HD, @RG, @PG etc.) ──────────────────────────────
# Pass through unchanged
/^@/ { print $0; next }

# ── Alignment lines ───────────────────────────────────────────────────────────
{
    # Extract key SAM fields by column position
    flag  = $2   # Bitwise flag encoding read properties (mapped/unmapped, strand etc.)
    rname = $3   # Reference chromosome name
    pos   = $4   # 1-based leftmost mapping position (0 if unmapped)
    cigar = $6   # CIGAR string describing the alignment

    # Track every chromosome name seen in the file for the summary table
    seen_chroms[rname] = 1

    # ── Unmapped reads ────────────────────────────────────────────────────────
    # A read is unmapped if bit 3 of the flag is set (decimal 4), or if POS
    # is 0, or if the CIGAR string is absent. Pass through unchanged.
    if (and(flag, 4) || pos == 0 || cigar == "*") {
        unmapped[rname]++
        print $0
        next
    }

    # ── Non-circular chromosomes ──────────────────────────────────────────────
    # Reads on non-circular chromosomes need no correction. Pass through unchanged.
    if (!(rname in is_circular)) {
        normal[rname]++
        print $0
        next
    }

    # ── Circular chromosome — classify and correct the read ───────────────────
    genome_len     = chrom_length[rname]  # Original genome length (without overhang)
    overhang_start = genome_len + 1       # First position of the overhang region (1-based)

    # Compute the reference-consuming length of the alignment from the CIGAR string.
    # Only operators that consume the reference (M, D, N, X, =) contribute to
    # the length — insertions (I), soft clips (S), and hard clips (H) do not.
    ref_len = 0
    tmp = cigar
    while (match(tmp, /[0-9]+[MDNX=]/)) {
        op_str = substr(tmp, RSTART, RLENGTH)
        # Add the numeric part of the operator (everything except the final letter)
        ref_len += substr(op_str, 1, length(op_str) - 1)
        # Consume the matched operator from the front of the CIGAR string
        tmp = substr(tmp, RSTART + RLENGTH)
    }

    # The rightmost reference position this read covers (1-based, inclusive)
    read_end = pos + ref_len - 1

    if (pos >= overhang_start) {
        # ── Case 1: Read lies entirely within the overhang ────────────────────
        # Subtract the genome length to wrap the position back to its true
        # coordinate at the start of the genome
        overhang_counts[rname]++
        $4 = pos - genome_len

        # For paired-end reads, also fix the mate position (PNEXT, field 8)
        # if the mate also falls within the overhang region
        if ($8 >= overhang_start) { $8 = $8 - genome_len }

        print $0

    } else if (read_end > genome_len) {
        # ── Case 2: Read spans the junction ───────────────────────────────────
        # The read starts in the normal region but ends in the overhang.
        # This alignment is unreliable so we flag the read as unmapped:
        #   $2: set the unmapped bit (bit 3, decimal 4) using bitwise OR
        #   $4: set POS to 0 (the SAM spec value for unmapped reads)
        #   $6: set CIGAR to * (the SAM spec value for absent CIGAR strings)
        junction[rname]++
        $2 = or(flag, 4)
        $4 = 0
        $6 = "*"

        print $0

    } else {
        # ── Case 3: Read lies entirely within the normal region ───────────────
        # No correction needed
        normal[rname]++
        print $0
    }
}

# ── Summary statistics ────────────────────────────────────────────────────────
# Print a per-chromosome breakdown of read classifications to stderr so it
# doesn't interfere with the SAM output written to stdout
END {
    col1 = 16  # Width of the chromosome name column
    col2 = 10  # Width of each count column

    # Print the table header and a divider line to stderr
    printf "%-16s %10s %10s %10s %10s %10s\n",
        "Chromosome",
        "Normal",
        "Overhang",
        "Junction",
        "Unmapped",
        "Total" > "/dev/stderr"

    divider = ""
    for (i = 1; i <= col1 + (col2 + 1) * 5; i++) divider = divider "-"
    print divider > "/dev/stderr"

    total_normal = 0; total_overhang = 0
    total_junction = 0; total_unmapped = 0

    # Sort chromosome names alphabetically using gawk's asorti()
    n = asorti(seen_chroms, sorted_chroms)

    for (i = 1; i <= n; i++) {
        chr = sorted_chroms[i]

        # Use 0 as default if a chromosome has no reads in a given category
        n_val  = (chr in normal)          ? normal[chr]          : 0
        ov_val = (chr in overhang_counts) ? overhang_counts[chr] : 0
        jn_val = (chr in junction)        ? junction[chr]        : 0
        un_val = (chr in unmapped)        ? unmapped[chr]        : 0
        tot    = n_val + ov_val + jn_val + un_val

        # Mark circular chromosomes with an asterisk in the output table
        label = (chr in is_circular) ? chr "*" : chr

        printf "%-16s %10d %10d %10d %10d %10d\n",
            label,
            n_val,
            ov_val,
            jn_val,
            un_val,
            tot > "/dev/stderr"

        total_normal   += n_val
        total_overhang += ov_val
        total_junction += jn_val
        total_unmapped += un_val
    }

    # Print the totals row
    print divider > "/dev/stderr"
    printf "%-16s %10d %10d %10d %10d %10d\n",
        "TOTAL",
        total_normal,
        total_overhang,
        total_junction,
        total_unmapped,
        total_normal + total_overhang + total_junction + total_unmapped > "/dev/stderr"
    printf "\n* = circular chromosome (overhang correction applied)\n\n" > "/dev/stderr"
}