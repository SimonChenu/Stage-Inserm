require(gaston)
require(OmicKriging)
x1 <- read.bed.matrix("EOAD") #EOAD etant la bed matrice (avec .fam et .bim) apres QC

#MODELE GCTA
system("./ldak5.1.linux --calc-kins-direct grm --bfile EOAD --ignore-weights YES --power -1")
details <- read.table('grm.grm.details', header = T, sep ='')
write.table(details, 'grm.grm.N.bin', quote = FALSE)
grm_initial <- read_GRMBin(prefix = '/home/simon/grm.grm') 
eigK <- eigen(grm_initial)
eigK$values [eigK$values < 0] <- 0
pc <- sweep(eigK$vectors, 2, sqrt(eigK$values), '*')
covar_pheno(10)
system("./ldak5.1.linux --adjust-grm pcgc.covar --grm grm --covar QC.covar")
system("./ldak5.1.linux --pcgc resultats --pheno QC.pheno  --grm pcgc.covar --covar QC.covar --prevalence .02")

#MODELE LDAK
gc()
system("./ldak5.1.linux --cut-weights sections --bfile EOAD")
system("./ldak5.1.linux --calc-weights-all sections --bfile EOAD")
system("./ldak5.1.linux --calc-kins-direct grm --bfile EOAD --weights sections/weights.all --power -0.25")
details <- read.table('grm.grm.details', header = T, sep ='')
write.table(details, 'grm.grm.N.bin', quote = FALSE)
grm_initial <- read_GRMBin(prefix = '/home/simon/grm.grm')
eigK <- eigen(grm_initial)
eigK$values [eigK$values < 0] <- 0
pc <- sweep(eigK$vectors, 2, sqrt(eigK$values), '*')
covar_pheno(10)
system("./ldak5.1.linux --adjust-grm pcgc.covar --grm grm --covar QC.covar")
system("./ldak5.1.linux --pcgc resultats --pheno QC.pheno  --grm pcgc.covar --covar QC.covar --prevalence .02")

#MODELE LDAK-Thin
Thin <- read.table("sections/thin.in")
Thin[,2] <- 1
write.table(Thin, 'Thin.in', quote = FALSE, row.names=FALSE, col.names = FALSE)
system("./ldak5.1.linux --calc-kins-direct grm --bfile EOAD --weights Thin.in --power -0.25")
details <- read.table('grm.grm.details', header = T, sep ='')
write.table(details, 'grm.grm.N.bin', quote = FALSE)
grm_initial <- read_GRMBin(prefix = '/home/simon/grm.grm')
eigK <- eigen(grm_initial)
eigK$values [eigK$values < 0] <- 0
pc <- sweep(eigK$vectors, 2, sqrt(eigK$values), '*')
covar_pheno(10)
system("./ldak5.1.linux --adjust-grm pcgc.covar --grm grm --covar QC.covar")
system("./ldak5.1.linux --pcgc resultats --pheno QC.pheno  --grm pcgc.covar --covar QC.covar --prevalence .02")
