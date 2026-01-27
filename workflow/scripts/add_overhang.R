library(Biostrings)

genome <- readDNAStringSet(snakemake@input[[1]])
shift_amounts <- snakemake@params[["overhang"]][[1]]

for (i in seq_along(shift_amounts)) {
    chr_name <- names(shift_amounts)[i]
    n <- shift_amounts[[i]]

    
    chr <- genome[[chr_name]]
    chr_overhang <- subseq(chr, 1, n)
    chr_with_overhang <- xscat(chr, chr_overhang)
    genome[[chr_name]] <- chr_with_overhang
}

writeXStringSet(genome, snakemake@output[[1]])