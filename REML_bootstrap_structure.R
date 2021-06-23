set.seed(sample(1:100e5,1))
require(gaston)
require(OmicKriging)
#seuil de liabilite pour K = 1%
seuil <- qnorm(0.99, mean = 0, sd=1)

x1 <- generation(0.8,0.5,0)
gc()
x1 <- set.hwe(x1)
x1 <- select.snps(x1, x1@snps$hwe > 0.001)
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
PCGC <- as.numeric(her_liab[23,2])
system("./ldak5.1.linux --reml resultats --pheno QC.pheno --grm grm --covar QC.covar --prevalence .01")
her_liab <- read.table("/home/simon/resultats.reml.liab")
her_liab <- as.matrix(her_liab)
REML <- as.numeric(her_liab[26,2])

#BOOTSTRAP 100 ITERATIONS
stockage <- rep(NA,100)
for (i in 1:100) {
  gc()
  liste1 <- x1@ped$id[x1@ped$id & x1@ped$pheno == 1]
  liste2 <- x1@ped$id[x1@ped$id & x1@ped$pheno == 0]
  #permet de garantir un Ke proche de 0.875
  liste2 <- sample(liste2, length(liste1)*125/875)
  liste <- c(liste1,liste2)
  x2 <- select.inds(x1, x1@ped$id %in% liste == TRUE)
  write.bed.matrix(x2, 'boot')
  system("./ldak5.1.linux --calc-kins-direct reml --bfile boot --ignore-weights YES --power -1")
  details_reml <- read.table('reml.grm.details', header = T, sep ='')
  write.table(details_reml, 'reml.grm.N.bin', quote = FALSE)
  grm_reml <- read_GRMBin(prefix = '/home/simon/reml.grm')
  eigK <- eigen(grm_reml)
  eigK$values [eigK$values < 0] <- 0
  pc <- sweep(eigK$vectors, 2, sqrt(eigK$values), '*')
  covar_pheno_bootstrap(10)
  system("./ldak5.1.linux --reml resultats --pheno QC.pheno --grm reml --covar QC.covar --prevalence .01")
  her_liab <- read.table("/home/simon/resultats.reml.liab")
  her_liab <- as.matrix(her_liab)
  stockage[i] <- as.numeric(her_liab[26,2])
  #100 iterations minimum
  cat(
    'Moyenne apres', i, 'iteration(s) :', mean(stockage[1:i])

  )
}

REML_bootstrap <- mean(stockage[1:100])
var(stockage[1:100])
