txdb <- txdbmaker::makeTxDbFromUCSC(genome = "sacCer3", tablename = "sgdGene") #import gtf annotation from UCSC
GenomeInfoDb::seqlevels(txdb) <- GenomeInfoDb::seqlevels(txdb)[1:16] #use only annotations on canonical chromosomes

prom_rng <- GenomicFeatures::promoters(txdb, upstream = 200, downstream = 50, use.names = TRUE) 
prom_rng <- GenomicRanges::trim(prom_rng) 
prom_seq <- Biostrings::getSeq(x = BSgenome.Scerevisiae.UCSC.sacCer3::BSgenome.Scerevisiae.UCSC.sacCer3, 
                               prom_rng) #promoter sequences
opts <- list()
opts[["collection"]] <- "CORE"
opts[["tax_group"]] <- "fungi"


J2022 <- TFBSTools::getMatrixSet(JASPAR2022::JASPAR2022, opts) #core Jaspar 2022 profiles for fungi

J2022_PSBG <- PscanR::ps_build_bg(prom_seq, J2022, BPPARAM = BiocParallel::MulticoreParam(12)) #Build Pscan Background

PscanR::ps_write_bg_to_file(J2022_PSBG, "J2022_sacCer3_200u_50d_UCSC.psbg.txt")