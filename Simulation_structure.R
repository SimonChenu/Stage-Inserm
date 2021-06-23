require(gaston)
require(OmicKriging)
#seuil de liabilite pour K=1%
seuil <- qnorm(0.99, mean = 0, sd=1)
#refresh memoire vive
gc()

#stockage valeurs
PCGC <- rep(NA,20)
REML <- rep(NA,20)

PCGC_extreme <- rep(NA,20)
REML_extreme <- rep(NA,20)


#Equilibre de structure cas et controles

for (a in 1:20){
  x1 <- generation(0.5,0.5,0)
  gc()
  x1 <- set.hwe(x1)
  x1 <- select.snps(x1, x1@snps$hwe > 0.001)
  gc()
  write.bed.matrix(x1, 'simu')
  
  system("./ldak5.1.linux --calc-kins-direct grm --bfile simu --ignore-weights YES --power -1")
  details <- read.table('grm.grm.details', header = T, sep ='')
  write.table(details, 'grm.grm.N.bin', quote = FALSE)
  grm_initial <- read_GRMBin(prefix = '/home/simon/grm.grm')
  eigK <- eigen(grm_initial)
  eigK$values [eigK$values < 0] <- 0
  pc <- sweep(eigK$vectors, 2, sqrt(eigK$values), '*')
  covar_pheno(10)
  system("./ldak5.1.linux --adjust-grm pcgc.covar --grm grm --covar QC.covar")
  system("./ldak5.1.linux --pcgc resultats --pheno QC.pheno --covar QC.covar --grm pcgc.covar --prevalence .01")
  her_liab <- read.table("/home/simon/resultats.pcgc")
  her_liab <- as.matrix(her_liab)
  PCGC[a] <- as.numeric(her_liab[23,2])
  system("./ldak5.1.linux --reml resultats --pheno QC.pheno --grm grm --covar QC.covar --prevalence .01")
  her_liab <- read.table("/home/simon/resultats.reml.liab")
  her_liab <- as.matrix(her_liab)
  REML[a] <- as.numeric(her_liab[26,2])
} 
  


#Fort desequilibre de structure dans les controles
for (i in 1:20){
  
  x1 <- generation(0.5,0.5,0.5)
  gc()
  x1 <- set.hwe(x1)
  x1 <- select.snps(x1, x1@snps$hwe > 0.001)
  gc()
  write.bed.matrix(x1, 'simu')
  
  system("./ldak5.1.linux --calc-kins-direct grm --bfile simu --ignore-weights YES --power -1")
  details <- read.table('grm.grm.details', header = T, sep ='')
  write.table(details, 'grm.grm.N.bin', quote = FALSE)
  grm_initial <- read_GRMBin(prefix = '/home/simon/grm.grm')
  eigK <- eigen(grm_initial)
  eigK$values [eigK$values < 0] <- 0
  pc <- sweep(eigK$vectors, 2, sqrt(eigK$values), '*')
  covar_pheno(10)
  system("./ldak5.1.linux --adjust-grm pcgc.covar --grm grm --covar QC.covar")
  system("./ldak5.1.linux --pcgc resultats --pheno QC.pheno --covar QC.covar --grm pcgc.covar --prevalence .01")
  her_liab <- read.table("/home/simon/resultats.pcgc")
  her_liab <- as.matrix(her_liab)
  PCGC_extreme[i] <- as.numeric(her_liab[23,2])
  system("./ldak5.1.linux --reml resultats --pheno QC.pheno --grm grm --covar QC.covar --prevalence .01")
  her_liab <- read.table("/home/simon/resultats.reml.liab")
  her_liab <- as.matrix(her_liab)
  REML_extreme[i] <- as.numeric(her_liab[26,2])
}
