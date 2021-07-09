#Code pour aboutir aux graphiques de h2 en fonction des PC en covariables

require(gaston)
require(OmicKriging)
#Calcul de la GRM par modele GCTA
x1 <- generation(0.5,0.5,0)
write.bed.matrix(x1, 'simu')
system("./ldak5.1.linux --calc-kins-direct grm --bfile simu --ignore-weights YES --power -1")
details <- read.table('grm.grm.details', header = T, sep ='')
write.table(details, 'grm.grm.N.bin', quote = FALSE)
grm_initial <- read_GRMBin(prefix = '/home/simon/grm.grm')
eigK <- eigen(grm_initial)
eigK$values [eigK$values < 0] <- 0
pc <- sweep(eigK$vectors, 2, sqrt(eigK$values), '*')

#Creation de la matrice de stockage des resultats pour 100 PC
ec <- matrix(nrow=303, ncol=4)
ec <- as.data.frame(ec)
colnames(ec) <- c('Methode','PC','h2', 'sd')
ec$Methode[1:101] <- 'PCGC_LDAK'
ec$Methode[102:202] <- 'REML'
ec$Methode[203:303] <- 'PCGC_AUTEURS'
ec$PC[1:101] <- seq(0,100,1)
ec$PC[102:202] <- seq(0,100,1)
ec$PC[203:303] <- seq(0,100,1)
matrice_pheno <- matrix(nrow=length(x1@ped[["famid"]]), ncol=3)
matrice_pheno[,1] <- x1@ped[["famid"]]
matrice_pheno[,2] <- x1@ped[["id"]]
matrice_pheno[,3] <- x1@ped[["pheno"]]
write.table(matrice_pheno, file='QC.pheno', quote=FALSE, row.names=FALSE, col.names = FALSE)
# CALCUL PCGC NO COVARIABLE
system("./ldak5.1.linux --pcgc resultats --pheno QC.pheno --grm grm --prevalence .02")
her_liab <- read.table("/home/simon/resultats.pcgc")
her_liab <- as.matrix(her_liab)
ec$h2[1] <- as.numeric(her_liab[17,2])
ec$sd[1] <- as.numeric(her_liab[18,1])
#REML NO COVARIABLE
system("./ldak5.1.linux --reml resultats --pheno QC.pheno --grm grm --prevalence .02")
her_liab <- read.table("/home/simon/resultats.reml.liab")
her_liab <- as.matrix(her_liab)
ec$h2[102] <- as.numeric(her_liab[26,2])
ec$sd[102] <- as.numeric(her_liab[27,1])
#PCGC LDAK 0 PC COVARIABLE
ec$h2[203] <- ec$h2[1]
ec$sd[203] <- ec$sd[1]
for (i in 1:100) {
  PC <- matrix(nrow=length(x1@ped[["famid"]]), ncol=i+2)
  PC[,1] <- x1@ped[["famid"]]
  PC[,2] <- x1@ped[["id"]]
  for (j in 1:i){
    PC[,j+2] <- pc[,j]
  }
  write.table(PC, file='QC.covar', quote=FALSE, row.names=FALSE, col.names = FALSE)
  #PCGC LDAK
  system("./ldak5.1.linux --adjust-grm pcgc.covar --grm grm --covar QC.covar")
  system("./ldak5.1.linux --pcgc resultats --pheno QC.pheno --covar QC.covar --grm pcgc.covar --prevalence .02")
  her_liab <- read.table("/home/simon/resultats.pcgc")
  her_liab <- as.matrix(her_liab)
  ec$h2[i+1] <- as.numeric(her_liab[17,2])
  ec$sd[i+1] <- as.numeric(her_liab[18,1])
  #REML
  system("./ldak5.1.linux --reml resultats --pheno QC.pheno --covar QC.covar --grm grm --prevalence .02")
  her_liab <- read.table("/home/simon/resultats.reml.liab")
  her_liab <- as.matrix(her_liab)
  ec$h2[i+102] <- as.numeric(her_liab[26,2])
  ec$sd[i+102] <- as.numeric(her_liab[27,1])
  #PCGC AUTEURS
  her_liab <- read.table("/home/simon/resultats.pcgc.marginal")
  her_liab <- as.matrix(her_liab)
  ec$h2[i+203] <- as.numeric(her_liab[17,2])
  ec$sd[i+203] <- as.numeric(her_liab[18,1])
}

require(ggplot2)
ec$Methode <- factor(ec$Methode, levels=c("PCGC_LDAK","PCGC_AUTEURS", "REML"))
ec$Methode <- as.factor(ec$Methode)

p <- ggplot(ec,aes(x=PC,y=h2,group=Methode,color=Methode))+
  geom_line()+
  geom_point()+
  geom_errorbar(aes(ymin=h2-sd,ymax=h2+sd),width=3,position=position_dodge(0.05))+
  labs(title=" ",x="PC en covariables", y='h²')+
  scale_color_manual(values=c("#1A1AB3","#53B4E9","#CC1A4D99"))
print(p)

p+geom_hline(yintercept=0.5, color='red')

