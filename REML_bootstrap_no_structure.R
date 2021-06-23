require(gaston)
require(genio)

#h2_sim = heritabilite simulee
#p = nb de variants
#n = nb d'individus
#K = prevalence maladie dans la population
#Ke = prevalence dans l'etude
#n_simu = nombre de simulations
#iterations = nombre d'iterations pour bootstrap

simulation_bootstrap <- function(h2_sim,p,n,K,Ke,n_simu,iterations) {
  stockage_heritabilite_PCGC <- rep(NA,n_simu)
  stockage_heritabilite_GCTA <- rep(NA,n_simu)
  stockage_heritabilite_GCTA_bootstrap <- rep(NA,n_simu)
  stockage_correlation <- rep(NA,n_simu)
  for (b in 1:n_simu){
    stockage_G <- rep(NA,n)
    stockage_E <- rep(NA,n)
    cat(
      'Simulation', b, 'sur', n_simu
      )
    #Simulation des frequences alleliques
    f <- runif(p,0.05,0.5)
    # Simulation du vecteur d'effets des snps u, et determination de la variance de G 
    u <- rnorm(p, 0, sqrt(h2_sim/p))
    g <- rep(NA,10000)
    for (i in 1:10000){
      genotypes <- rbinom (p, 2, f)
      g[i] <- sum(genotypes*u)
    }
    var_g <- var(g)
    
    #Determination du seuil t de liabilite adequat
    cat("
    Etablissement du seuil t adequat  ...
      ")
    l <- rep(NA,2000)
    stockage_t <- rep(NA,100)
    for (i in 1:100) {
      for (j in 1:2000){
        genotypes <- rbinom (p,2,f)
        g <- sum(genotypes*u)
        e <- rnorm(1, 0, sqrt(var_g*(1/h2_sim-1)))
        l[j] <- e+g
      }
      l <- sort(l,decreasing=TRUE)
      stockage_t[i] <- l[2000*K]
    }
    t <- mean(stockage_t)
    matrice_stockage <- matrix(nrow=p, ncol=n)
    save_pheno <- rep(NA, n)
    n_etude <- 0
    cat("
        Inclusion de", n, 'individus ...
        '
    )
    while(n_etude < n){
      Z <- rbinom (p, 2, f)
      g <- sum(Z*u)
      e <- rnorm(1, 0, sqrt(var_g*(1/h2_sim-1)))
      l <- g + e
      ifelse (l > t, y <- 2, y <- 1) 
      if (y==2){
        matrice_stockage[,n_etude+1] <- Z
        save_pheno[n_etude+1] <- y
        n_etude <- n_etude+1
        stockage_G[n_etude] <- g
        stockage_E[n_etude] <- e
      } 
      else {
        s <- as.numeric(rbinom (1, 1, (K*(1-Ke)) / (Ke*(1-K))))
        if (s==1){
          matrice_stockage[,n_etude+1] <- Z
          save_pheno[n_etude+1] <- y
          n_etude <- n_etude+1
          stockage_G[n_etude] <- g
          stockage_E[n_etude] <- e
        }
      }
      stockage_correlation[b] <- cor(stockage_G,stockage_E,method="pearson")
    }
    pheno <- matrix(nrow=n,ncol=3)
    pheno[,1] <- 1:n
    pheno[,2] <- 1:n
    pheno[,3] <- save_pheno
    #creation du .bed .fam et .bim
    write_plink(file = "/home/simon/bootstrap", X=matrice_stockage)
    write.table(pheno, file="simu.pheno", row.names=FALSE, col.names=FALSE)
    
    #PCGC
    system("./ldak5.1.linux --calc-kins-direct grm --bfile bootstrap --ignore-weights YES --power -1")
    system("./ldak5.1.linux --pcgc resultats --pheno simu.pheno --grm grm --prevalence .01")
    her_liab <- read.table("/home/simon/resultats.pcgc")
    her_liab <- as.matrix(her_liab)
    stockage_heritabilite_PCGC[b] <- as.numeric(her_liab[23,2])
    
    #REML
    system("./ldak5.1.linux --reml resultats --pheno simu.pheno --grm grm --prevalence .01")
    her_liab <- read.table("/home/simon/resultats.reml.liab")
    her_liab <- as.matrix(her_liab)
    stockage_heritabilite_GCTA[b] <- as.numeric(her_liab[26,2])
    #REML bootstrap
    x1 <- read.bed.matrix('bootstrap')
    x1@ped$pheno <- pheno[,3]
    stockage <- rep(NA,iterations)
    for (z in 1:iterations) {
      liste1 <- x1@ped$id[x1@ped$id & x1@ped$pheno == 2]
      liste2 <- x1@ped$id[x1@ped$id & x1@ped$pheno == 1]
      #garantir Ke = 0.875
      liste2 <- sample(liste2, length(liste1)*125/875)
      liste <- c(liste1,liste2)
      x2 <- select.inds(x1, x1@ped$id %in% liste == TRUE)
      write.bed.matrix(x2, 'bootstrap2')
      matrice_pheno <- matrix(nrow=length(x2@ped[["famid"]]), ncol=3)
      matrice_pheno[,1] <- x2@ped[["famid"]]
      matrice_pheno[,2] <- x2@ped[["id"]]
      matrice_pheno[,3] <- x2@ped[["pheno"]]
      write.table(matrice_pheno, file='simu.pheno', quote=FALSE, row.names=FALSE, col.names = FALSE)
      system("./ldak5.1.linux --calc-kins-direct reml --bfile bootstrap2 --ignore-weights YES --power -1")
      system("./ldak5.1.linux --reml resultats --pheno simu.pheno --grm reml --prevalence .01")
      her_liab <- read.table("/home/simon/resultats.reml.liab")
      her_liab <- as.matrix(her_liab)
      stockage[z] <- as.numeric(her_liab[26,2])
    }
    stockage_heritabilite_GCTA_bootstrap[b] <- mean(stockage)
  }
  x <- matrix(nrow = n_simu, ncol = 4)
  x <- as.data.frame(x)
  colnames(x) <- c('REML', 'REML_bootstrap', 'PCGC', 'Correlation GxE')
  x$REML <- stockage_heritabilite_GCTA
  x$PCGC <- stockage_heritabilite_PCGC
  x$`Correlation GxE` <- stockage_correlation
  x$REML_bootstrap <- stockage_heritabilite_GCTA_bootstrap
  return(x)
}

x <- simulation_bootstrap(0.5,5000,2000,0.01,0.5,1,100)
