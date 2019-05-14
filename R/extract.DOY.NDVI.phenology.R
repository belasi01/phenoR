#' Extrait des événements phénologiques d'un cycle annuel de NDVI
#' à partir de données extraites de
#' Google Earth Engine générées avec la fonction
#' ui.Chart.image.doySeriesByYear.
#'
#' La fonction trouve les dates (jour de l'année) de
#' plusieurs événements phénologiques
#' en fonction des seuils qui sont passés à la
#' fonction comme paramètres:
#' le jour du début de la saison de croissance (C);
#' le jour où le NDVI atteint son maximum (M);
#' le jour marquant le début de la sénescence (S);
#' le jour marquant le début de l'hiver (W);
#' la durée de la période de croissance (DC = M - C);
#' la durée de la période de la saison de production  (DP = S - C);
#' la durée de la période de sénescence (DS = W - S).
#'
#'
#' @param filen est le nom du fichier CVS telechargé de GEE et produit par
#' ui.Chart.image.doySeriesByYear. Le fichier peut contenir plusieurs annees.
#' @param SpringBloom.Threshold fraction de la valeur maximum de NDVI que l'on
#' considere en début de saison pour identifier le début de la croissance.
#' (valeur par défault est 0.25)
#' @param Maximum.Threshold fraction de la valeur maximum de NDVI que l'on
#' considère en début de saison pour identifier la date où la végétation a
#' atteint sa maturité.
#'  (valeur par défault est 0.9)
#' @param Scenescence.Threshold fraction de la valeur maximum de NDVI que l'on
#' considère à l'automne quand la végétation commence sa sénescence.
#'  (valeur par défault est 0.85)
#' @param Winter.Threshold fraction de la valeur maximum de NDVI que l'on
#' considère après la sénescence.
#'  (valeur par défault est 0.5)
#'
#'
#' @return Retourne un tableau avec le jour de l'année de chaque événements.
#' De plus, une figure en format png est créée dans le répertoire de travail pour
#' chaque année contenue dans le fichier CSV.
#'
#' @author Simon Bélanger
#'
#' @export
#'
extract.DOY.NDVI.phenology <- function(filen,
                                  SpringBloom.Threshold=0.25,
                                  Maximum.Threshold=0.9,
                                  Scenescence.Threshold=0.85,
                                  Winter.Threshold=0.5) {
  ndvi <- read.csv(filen, header = T)

  nyear=dim(ndvi)[2]-1
  cyears=colnames(ndvi)[2:(nyear+1)]
  fyears=as.numeric(str_sub(cyears, 2,5))

  C<- rep(NA, nyear)
  M<- rep(NA, nyear)
  S<- rep(NA, nyear)
  W<- rep(NA, nyear)


  for(i in 1:nyear) {

    ## Normalisation de la série par le maximum atteint
    df = ndvi[,c(1,(i+1))]
    names(df) <- c("doy", "year")

    ndvi_n <-  (df$year - min(df$year,na.rm=T))/
                (max(df$year,na.rm=T) - min(df$year,na.rm=T))

    df = data.frame(doy=ndvi$doy, ndvi_n)
    names(df) <- c("doy", "year")
    smoothed <- loess(year ~ doy, data=df, span = 0.1)
    x  = predict(smoothed, 1:365)

    delta = (x[2:365] - x[1:364])
    delta = c(delta,NA)


    C[i] <- which(x > SpringBloom.Threshold & delta > 0)[1]
    M[i] <- which(x > Maximum.Threshold)[1]
    tmp <- which(x > Scenescence.Threshold)
    S[i] <- tmp[length(tmp)]
    tmp <- which(x > Winter.Threshold)
    W[i] <- tmp[length(tmp)]

    png(file=paste("ndvi_",cyears[i],".png",sep=""),res=300,units="in", width=7, height=5)

    plot(df,
         xlab="Jour de l'annee",
         ylab="NDVI normalisé",
         pch=19,
         col="grey",
         main=paste("Série Temporelle: NDVI", fyears[i]))
    lines(1:C[i], x[1:C[i]], col="blue", lwd=5)
    lines(C[i]:M[i], x[C[i]:M[i]], col="lightgreen", lwd=5)
    lines(M[i]:S[i], x[M[i]:S[i]], col="darkgreen", lwd=5)
    lines(S[i]:W[i], x[S[i]:W[i]], col="darkgoldenrod2", lwd=5)
    lines(W[i]:365, x[W[i]:365], col="blue", lwd=5)
    text(0,1, paste("Debut de la croissance:", C[i]),pos=4)
    text(0,0.9, paste("Duree de la croissance:", M[i]-C[i]),pos=4)
    text(0,0.8, paste("Debut de la scenescence:", S[i]),pos=4)
    text(0,0.7, paste("Duree de la saison de production:", S[i]-C[i]),pos=4)
    text(0,0.6, paste("Duree de la scenescence:", W[i]-S[i]),pos=4)
    dev.off()
  }
  DC  <- (M - C)
  DP  <- (S - C)
  DS  <- (W - S)

  return(data.frame(year=fyears,C,M,S,W,DC,DP,DS))

}

