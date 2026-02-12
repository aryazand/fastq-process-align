# ─────────────────────────────────────────────────────────────────────────────
# adjust_bam_header.awk
# ─────────────────────────────────────────────────────────────────────────────
# Extract and modify SAM header lines to correct chromosome lengths.
# Removes overhang length from @SQ LN: tags for circular chromosomes.

BEGIN { 
    FS = OFS = "\t"
    
    # Parse chromosome:length pairs into a lookup table
    n = split(chrom_lengths, lens, ",")
    for (i = 1; i <= n; i++) {
        split(lens[i], kv, ":")
        chrom_length[kv[1]] = kv[2]
    }
}

# Process @SQ header lines
/^@SQ/ {
    # Extract chromosome name from the SN: tag
    sn = ""
    for (i = 1; i <= NF; i++) {
        if ($i ~ /^SN:/) sn = substr($i, 4)
    }
    
    # Replace LN: tag with original length (no overhang) if known
    if (sn in chrom_length) {
        for (i = 1; i <= NF; i++) {
            if ($i ~ /^LN:/) $i = "LN:" chrom_length[sn]
        }
    }
}

# Print all lines (modified or unchanged)
{ print }