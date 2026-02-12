BEGIN {

    FS = OFS = "\t"
    
    # Parse circular chromosome names into a lookup set
    n = split(circular_chroms, chroms, ",")
    for (i = 1; i <= n; i++) {
        is_circular[chroms[i]] = 1
    }

    # Parse per-chromosome genome lengths into a lookup table
    n = split(chromosome_lengths, lens, ",")
    for (i = 1; i <= n; i++) {
        split(lens[i], kv, ":")
        chrom_length[kv[1]] = kv[2]
    }
}

# ── Header lines ─────────────────────────────────────────────────────────────
/^@/ {
    if ($0 ~ /^@SQ/) {
        sn = ""
        for (i = 1; i <= NF; i++) {
            if ($i ~ /^SN:/) { sn = substr($i, 4) }
        }
        # Only replace LN if we have a valid length for this chrom
        if (sn in is_circular && (sn in chrom_length) && chrom_length[sn] > 0) {
            ln_found = 0
            for (i = 1; i <= NF; i++) {
                if ($i ~ /^LN:/) {
                    $i = "LN:" chrom_length[sn]
                    ln_found = 1
                }
            }
            if (!ln_found) {
                # append LN if missing
                NF = NF + 1
                $NF = "LN:" chrom_length[sn]
            }
            print $0
        } else {
            print $0
        }
    } else {
        print $0
    }
    next
}

# ── Alignment lines ───────────────────────────────────────────────────────────
{
    flag  = $2
    rname = $3
    pos   = $4
    cigar = $6

    # Leave unmapped reads or non-circular chromosomes untouched
    if ((and(flag, 4)) || !(rname in is_circular) || pos == 0 || cigar == "*") {
        print $0
        next
    }

    genome_len    = chrom_length[rname]
    overhang_start = genome_len + 1

    # Compute reference-consuming length from CIGAR
    ref_len = 0
    tmp = cigar
    while (match(tmp, /[0-9]+[MDNX=]/)) {
        op_str = substr(tmp, RSTART, RLENGTH)
        ref_len += substr(op_str, 1, length(op_str) - 1)
        tmp = substr(tmp, RSTART + RLENGTH)
    }
    read_end = pos + ref_len - 1

    if (pos >= overhang_start) {
        # Entirely in overhang → wrap position back
        $4 = pos - genome_len

        # Fix mate position if also in overhang (paired-end)
        if ($7 == rname || $7 == "=") {
            if ($8 >= overhang_start) { $8 = $8 - genome_len }
        }

    } else if (read_end > genome_len) {
        # Spans junction → flag as unmapped
        $2  = or(flag, 4)
        $4  = 0
        $6  = "*"
    }

    print $0
}