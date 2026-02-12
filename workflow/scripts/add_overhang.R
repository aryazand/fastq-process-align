library(Biostrings)

genome <- readDNAStringSet(snakemake@input[[1]])
circular_chrom_info <- snakemake@params[["overhang"]]

chrom_names <- unlist(circular_chrom_info$circular_chromosomes)
overhang_length <- circular_chrom_info$overhang_length

for (i in seq_along(chrom_names)) {
    chr_name <- chrom_names[i]
    n <- overhang_length

    # Using regex to match chromosomes that start with the specified name
    
    for(j in grep(paste0("^", chr_name), names(genome))) {
        chr <- genome[[j]]
        chr_overhang <- subseq(chr, 1, n)
        chr_with_overhang <- xscat(chr, chr_overhang)
        genome[[j]] <- chr_with_overhang

    }
}

writeXStringSet(genome, snakemake@output[[1]])