#Une serie de fonctions pour la partie simulatons sur structure

#creation des fichiers .covar et .pheno avec, format plink, n = nombre de PC
covar_pheno <- function(n){
  PC <- matrix(nrow=length(x1@ped[["famid"]]), ncol=n+2)
  PC[,1] <- x1@ped[["famid"]]
  PC[,2] <- x1@ped[["id"]]
  for (j in 1:n){
    PC[,j+2] <- pc[,j]
  }
  write.table(PC, file='QC.covar', quote=FALSE, row.names=FALSE, col.names = FALSE)
  matrice_pheno <- matrix(nrow=length(x1@ped[["famid"]]), ncol=3)
  matrice_pheno[,1] <- x1@ped[["famid"]]
  matrice_pheno[,2] <- x1@ped[["id"]]
  matrice_pheno[,3] <- x1@ped[["pheno"]]
  write.table(matrice_pheno, file='QC.pheno', quote=FALSE, row.names=FALSE, col.names = FALSE)
}

#idem pour le bootstrap de REML avec structure de population
covar_pheno_bootstrap <- function(n){
  PC <- matrix(nrow=length(x2@ped[["famid"]]), ncol=n+2)
  PC[,1] <- x2@ped[["famid"]]
  PC[,2] <- x2@ped[["id"]]
  for (j in 1:n){
    PC[,j+2] <- pc[,j]
  }
  write.table(PC, file='QC.covar', quote=FALSE, row.names=FALSE, col.names = FALSE)
  matrice_pheno <- matrix(nrow=length(x2@ped[["famid"]]), ncol=3)
  matrice_pheno[,1] <- x2@ped[["famid"]]
  matrice_pheno[,2] <- x2@ped[["id"]]
  matrice_pheno[,3] <- x2@ped[["pheno"]]
  write.table(matrice_pheno, file='QC.pheno', quote=FALSE, row.names=FALSE, col.names = FALSE)
}

#fonction generant la bed matrice avec structure de population
generation <- function(h2,Ke,score_IBS){
  u <- rnorm(444601, sd=sqrt(h2/444601))
  #1
  x <- read.bed.matrix('structure01')
  gc()
  standardize(x) <- "p"
  G <- x %*% u
  E <- rnorm(length(x@ped[["famid"]]), sd= sqrt(1-h2))
  Liab <- G+E
  x@ped[,35] <- G
  x@ped[,36] <- E
  
  for(i in 1:length(Liab)){
    ifelse(Liab[,1][i] >= seuil, Liab[,1][i] <- 1, Liab[,1][i] <- 0)
  }
  x@ped[["pheno"]] <- Liab[,1]

  liste1 <- x@ped$id[x@ped$id & x@ped$pheno == 1]
  liste2 <- x@ped$id[x@ped$id & x@ped$pheno == 0 & x@ped$IBS >= score_IBS]
  liste2 <- sample(liste2, length(liste1)*(1-Ke)/Ke)
  liste <- c(liste1,liste2)
  x1 <- select.inds(x, x@ped$id %in% liste == TRUE)
  gc()
  
  #2
  x <- read.bed.matrix('structure02')
  gc()
  standardize(x) <- "p"
  G <- x %*% u
  E <- rnorm(length(x@ped[["famid"]]), sd= sqrt(1-h2))
  Liab <- G+E
  x@ped[,35] <- G
  x@ped[,36] <- E
  
  for(i in 1:length(Liab)){
    ifelse(Liab[,1][i] >= seuil, Liab[,1][i] <- 1, Liab[,1][i] <- 0)
  }
  x@ped[["pheno"]] <- Liab[,1]
  
  liste1 <- x@ped$id[x@ped$id & x@ped$pheno == 1]
  liste2 <- x@ped$id[x@ped$id & x@ped$pheno == 0 & x@ped$IBS >= score_IBS]
  liste2 <- sample(liste2, length(liste1)*(1-Ke)/Ke)
  liste <- c(liste1,liste2)
  x2 <- select.inds(x, x@ped$id %in% liste == TRUE)
  gc()
  
  #3
  x <- read.bed.matrix('structure03')
  gc()
  standardize(x) <- "p"
  G <- x %*% u
  E <- rnorm(length(x@ped[["famid"]]), sd= sqrt(1-h2))
  Liab <- G+E
  x@ped[,35] <- G
  x@ped[,36] <- E
  
  for(i in 1:length(Liab)){
    ifelse(Liab[,1][i] >= seuil, Liab[,1][i] <- 1, Liab[,1][i] <- 0)
  }
  x@ped[["pheno"]] <- Liab[,1]
  
  liste1 <- x@ped$id[x@ped$id & x@ped$pheno == 1]
  liste2 <- x@ped$id[x@ped$id & x@ped$pheno == 0 & x@ped$IBS >= score_IBS]
  liste2 <- sample(liste2, length(liste1)*(1-Ke)/Ke)
  liste <- c(liste1,liste2)
  x3 <- select.inds(x, x@ped$id %in% liste == TRUE)
  gc()
  
  #4
  x <- read.bed.matrix('structure04')
  gc()
  standardize(x) <- "p"
  G <- x %*% u
  E <- rnorm(length(x@ped[["famid"]]), sd= sqrt(1-h2))
  Liab <- G+E
  x@ped[,35] <- G
  x@ped[,36] <- E
  
  for(i in 1:length(Liab)){
    ifelse(Liab[,1][i] >= seuil, Liab[,1][i] <- 1, Liab[,1][i] <- 0)
  }
  x@ped[["pheno"]] <- Liab[,1]
  
  liste1 <- x@ped$id[x@ped$id & x@ped$pheno == 1]
  liste2 <- x@ped$id[x@ped$id & x@ped$pheno == 0 & x@ped$IBS >= score_IBS]
  liste2 <- sample(liste2, length(liste1)*(1-Ke)/Ke)
  liste <- c(liste1,liste2)
  x4 <- select.inds(x, x@ped$id %in% liste == TRUE)
  gc()
  
  #5
  x <- read.bed.matrix('structure05')
  gc()
  standardize(x) <- "p"
  G <- x %*% u
  E <- rnorm(length(x@ped[["famid"]]), sd= sqrt(1-h2))
  Liab <- G+E
  x@ped[,35] <- G
  x@ped[,36] <- E
  
  for(i in 1:length(Liab)){
    ifelse(Liab[,1][i] >= seuil, Liab[,1][i] <- 1, Liab[,1][i] <- 0)
  }
  x@ped[["pheno"]] <- Liab[,1]
  
  liste1 <- x@ped$id[x@ped$id & x@ped$pheno == 1]
  liste2 <- x@ped$id[x@ped$id & (x@ped$pheno == 0 & x@ped$IBS >= score_IBS)]
  liste2 <- sample(liste2, length(liste1)*(1-Ke)/Ke)
  liste <- c(liste1,liste2)
  x5 <- select.inds(x, x@ped$id %in% liste == TRUE)
  gc()
  
  #6
  x <- read.bed.matrix('structure06')
  gc()
  standardize(x) <- "p"
  G <- x %*% u
  E <- rnorm(length(x@ped[["famid"]]), sd= sqrt(1-h2))
  Liab <- G+E
  x@ped[,35] <- G
  x@ped[,36] <- E
  
  for(i in 1:length(Liab)){
    ifelse(Liab[,1][i] >= seuil, Liab[,1][i] <- 1, Liab[,1][i] <- 0)
  }
  x@ped[["pheno"]] <- Liab[,1]
  
  liste1 <- x@ped$id[x@ped$id & x@ped$pheno == 1]
  liste2 <- x@ped$id[x@ped$id & x@ped$pheno == 0 & x@ped$IBS >= score_IBS]
  liste2 <- sample(liste2, length(liste1)*(1-Ke)/Ke)
  liste <- c(liste1,liste2)
  x6 <- select.inds(x, x@ped$id %in% liste == TRUE)
  gc()
  
  #7
  x <- read.bed.matrix('structure07')
  gc()
  standardize(x) <- "p"
  G <- x %*% u
  E <- rnorm(length(x@ped[["famid"]]), sd= sqrt(1-h2))
  Liab <- G+E
  x@ped[,35] <- G
  x@ped[,36] <- E
  
  for(i in 1:length(Liab)){
    ifelse(Liab[,1][i] >= seuil, Liab[,1][i] <- 1, Liab[,1][i] <- 0)
  }
  x@ped[["pheno"]] <- Liab[,1]
  
  liste1 <- x@ped$id[x@ped$id & x@ped$pheno == 1]
  liste2 <- x@ped$id[x@ped$id & x@ped$pheno == 0 & x@ped$IBS >= score_IBS]
  liste2 <- sample(liste2, length(liste1)*(1-Ke)/Ke)
  liste <- c(liste1,liste2)
  x7 <- select.inds(x, x@ped$id %in% liste == TRUE)
  gc()

  
  #8
  x <- read.bed.matrix('structure08')
  gc()
  standardize(x) <- "p"
  G <- x %*% u
  E <- rnorm(length(x@ped[["famid"]]), sd= sqrt(1-h2))
  Liab <- G+E
  x@ped[,35] <- G
  x@ped[,36] <- E
  
  for(i in 1:length(Liab)){
    ifelse(Liab[,1][i] >= seuil, Liab[,1][i] <- 1, Liab[,1][i] <- 0)
  }
  x@ped[["pheno"]] <- Liab[,1]
  
  liste1 <- x@ped$id[x@ped$id & x@ped$pheno == 1]
  liste2 <- x@ped$id[x@ped$id & x@ped$pheno == 0 & x@ped$IBS >= score_IBS]
  liste2 <- sample(liste2, length(liste1)*(1-Ke)/Ke)
  liste <- c(liste1,liste2)
  x8 <- select.inds(x, x@ped$id %in% liste == TRUE)
  gc()
  
  #9
  x <- read.bed.matrix('structure09')
  gc()
  standardize(x) <- "p"
  G <- x %*% u
  E <- rnorm(length(x@ped[["famid"]]), sd= sqrt(1-h2))
  Liab <- G+E
  x@ped[,35] <- G
  x@ped[,36] <- E
  
  for(i in 1:length(Liab)){
    ifelse(Liab[,1][i] >= seuil, Liab[,1][i] <- 1, Liab[,1][i] <- 0)
  }
  x@ped[["pheno"]] <- Liab[,1]
  
  liste1 <- x@ped$id[x@ped$id & x@ped$pheno == 1]
  liste2 <- x@ped$id[x@ped$id & x@ped$pheno == 0 & x@ped$IBS >= score_IBS]
  liste2 <- sample(liste2, length(liste1)*(1-Ke)/Ke)
  liste <- c(liste1,liste2)
  x9 <- select.inds(x, x@ped$id %in% liste == TRUE)
  gc()
  
  #10
  x <- read.bed.matrix('structure10')
  gc()
  standardize(x) <- "p"
  G <- x %*% u
  E <- rnorm(length(x@ped[["famid"]]), sd= sqrt(1-h2))
  Liab <- G+E
  x@ped[,35] <- G
  x@ped[,36] <- E
  
  for(i in 1:length(Liab)){
    ifelse(Liab[,1][i] >= seuil, Liab[,1][i] <- 1, Liab[,1][i] <- 0)
  }
  x@ped[["pheno"]] <- Liab[,1]
  liste1 <- x@ped$id[x@ped$id & x@ped$pheno == 1]
  liste2 <- x@ped$id[x@ped$id & x@ped$pheno == 0 & x@ped$IBS >= score_IBS]
  liste2 <- sample(liste2, length(liste1)*(1-Ke)/Ke)
  liste <- c(liste1,liste2)
  x10 <- select.inds(x, x@ped$id %in% liste == TRUE)
  gc()
  
  x <- rbind(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10)
  x@ped$id <- seq(1,length(x@ped$id),1)
  x@ped$famid <- seq(1,length(x@ped$famid),1)
  standardize(x) <- 'p'
  gc()
  return(x)
}

