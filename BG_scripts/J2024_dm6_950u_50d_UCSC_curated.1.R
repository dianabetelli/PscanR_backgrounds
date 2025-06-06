txdb <- txdbmaker::makeTxDbFromUCSC(genome="dm6", tablename="ncbiRefSeqCurated") #import gtf annotation from UCSC
GenomeInfoDb::seqlevels(txdb) <- GenomeInfoDb::seqlevels(txdb)[1:7] #use only annotations on canonical chromosomes

prom_rng <- GenomicFeatures::promoters(txdb, upstream = 950, downstream = 50, use.names = TRUE)
prom_rng <- GenomicRanges::trim(prom_rng)
prom_seq <- Biostrings::getSeq(x = BSgenome.Dmelanogaster.UCSC.dm6::BSgenome.Dmelanogaster.UCSC.dm6, prom_rng) #promoter sequences

opts <- list()
opts[["collection"]] <- "CORE"
opts[["tax_group"]] <- "insects"

httr::set_config(httr::config(ssl_verifypeer = 0L))
JASPAR2024 <- JASPAR2024::JASPAR2024()
JASPARConnect <- RSQLite::dbConnect(RSQLite::SQLite(), JASPAR2024::db(JASPAR2024))
J2024 <- TFBSTools::getMatrixSet(JASPARConnect, opts) #core Jaspar 2024 profiles for insects

J2024_PSBG <- PscanR::ps_build_bg(prom_seq, J2024, BPPARAM = BiocParallel::MulticoreParam(12)) #Build Pscan Background

PscanR::ps_write_bg_to_file(J2024_PSBG, "J2024_dm6_950u_50d_UCSC.psbg.txt")