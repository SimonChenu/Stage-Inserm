require(gaston)
require(genio)
#h2_sim = heritabilite simulee
#p = nb de variants
#n = nb d'individus
#K = prevalence maladie dans la population
#Ke = prevalence dans l'etude
#n_simu = nombre de simulations


simulation_norm <- function(h2_sim,p,n,K,Ke,n_simu) {
  stockage_heritabilite_PCGC <- rep(NA,n_simu)
  stockage_heritabilite_REML <- rep(NA,n_simu)
  stockage_correlation <- rep(NA,n_simu)
  t <- qnorm(1-K, mean = 0, sd=1) #seuil de liabilite
  for (b in 1:n_simu){
    stockage_G <- rep(NA,n)
    stockage_E <- rep(NA,n)
    cat(
      'Simulation', b, 'sur', n_simu
    )
    f <- runif(p,0.05,0.5) #Simulation des frequences alleliques
    u <- rnorm(p, 0, sqrt(h2_sim/p)) #Simulation des effets des variants
    E <- 2*f #Pour la normalisation
    SD <- sqrt(2*f*(1-f)) #Pour la normalisation
    matrice_stockage <- matrix(nrow=p, ncol=n) #Matrice de stockage des genotypes
    save_pheno <- rep(NA, n) #vecteur de sauvegarde des phenotypes
    n_etude <- 0
    cat("
        Inclusion de", n, 'individus ...
        '
    )
    while(n_etude < n){
      Z <- rbinom (p, 2, f) #tirage des genotypes
      Z_norm <- (Z - E)/SD #normalisation des genotypes
      g <- sum(Z_norm*u) #determination de G pour l'individu
      e <- rnorm(1, 0, sqrt(1-h2_sim)) #determination de E pour l'individu
      l <- g + e
      ifelse (l > t, y <- 2, y <- 1) 
      if (y==2){ #sauvegarde si l'individu est un cas
        matrice_stockage[,n_etude+1] <- Z
        save_pheno[n_etude+1] <- y
        n_etude <- n_etude+1
        stockage_G[n_etude] <- g
        stockage_E[n_etude] <- e
      } 
      else { #tirage si controle
        s <- as.numeric(rbinom (1, 1, (K*(1-Ke)) / (Ke*(1-K))))
        if (s==1){
          matrice_stockage[,n_etude+1] <- Z
          save_pheno[n_etude+1] <- y
          n_etude <- n_etude+1
          stockage_G[n_etude] <- g
          stockage_E[n_etude] <- e
        }
      }
      stockage_correlation[b] <- cor(stockage_G,stockage_E,method="pearson") #stockage des correlations
    }
    pheno <- matrix(nrow=n,ncol=3) #fichier pheno format plink pour le logiciel LDAK
    pheno[,1] <- 1:n
    pheno[,2] <- 1:n
    pheno[,3] <- save_pheno
    #creation du .bed .fam et .bim
    write_plink(file = "/home/simon/simu", X=matrice_stockage)
    write.table(pheno, file="simu.pheno", row.names=FALSE, col.names=FALSE)
    
    #PCGC
    system("./ldak5.1.linux --calc-kins-direct grm --bfile simu --ignore-weights YES --power -1") #Calcul de la GRM
    system("./ldak5.1.linux --pcgc resultats --pheno simu.pheno --grm grm --prevalence .01") # /!\ Changer la prevalence si changement du K de la fonction
    her_liab <- read.table("/home/simon/resultats.pcgc") #Lecture du fichier de resultat cree par LDAK
    her_liab <- as.matrix(her_liab)
    stockage_heritabilite_PCGC[b] <- as.numeric(her_liab[23,2]) #Sauvegarde de la valeur d'interet
    
    #REML
    system("./ldak5.1.linux --reml resultats --pheno simu.pheno --grm grm --prevalence .01")
    her_liab <- read.table("/home/simon/resultats.reml.liab")
    her_liab <- as.matrix(her_liab)
    stockage_heritabilite_REML[b] <- as.numeric(her_liab[26,2])
  }
  x <- matrix(nrow = n_simu, ncol = 3)
  x <- as.data.frame(x)
  colnames(x) <- c('REML', 'PCGC', 'Correlation GxE')
  x$REML <- stockage_heritabilite_REML
  x$PCGC <- stockage_heritabilite_PCGC
  x$`Correlation GxE` <- stockage_correlation
  return(x)
}

x <- simulation_norm(0.5,5000,2000,0.01,0.5,10)

