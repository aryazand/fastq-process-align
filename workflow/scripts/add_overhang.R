library(Biostrings)

genome <- readDNAStringSet(snakemake@input[[1]])
shift_amounts <- snakemake@params[["overhang"]]

for (i in seq_along(shift_amounts)) {
    chr_name <- names(shift_amounts)[i]
    n <- shift_amounts[[i]]
    
    print(n)

    # Using regex to match chromosomes that start with the specified name
    
    for(j in grep(paste0("^", chr_name), names(genome))) {
        chr <- genome[[j]]
        chr_overhang <- subseq(chr, 1, n)
        chr_with_overhang <- xscat(chr, chr_overhang)
        genome[[j]] <- chr_with_overhang

    }
}

writeXStringSet(genome, snakemake@output[[1]])