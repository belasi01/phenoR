#' Extrait des événements phénologiques d'un cycle annuel de NDSI
#' à partir de données extraites de
#' Google Earth Engine générées avec la fonction
#' ui.Chart.image.doySeriesByYear.
#'
#' La fonction trouve les dates (jour de l'année) de plusieurs événements phénologiques :
#' le jour du début de la fonte de la neige (Melt.Onset);
#' le jour de la fonte totale de la neige (Snow.Free);
#' le jour des premières neige d'automne (Fall.Snow);
#' la durée de la fonte (DF = Snow.Free - Melt.Onset);
#' la durée de la saison sans neige (DSSN= Fall.Snow - Snow.Free).
#'
#' @param filen est le nom du fichier CVS telechargé de GEE et produit par
#' ui.Chart.image.doySeriesByYear. Le fichier peut contenir plusieurs années.
#' @param Melt.Onset.Threshold fraction de la valeur maximum de NDSI que l'on
#' considere en début de saison pour identifier le début de la fonte de la neige.
#' (valeur par défault est 0.75)
#' @param Snow.Free.Threshold fraction de la valeur maximum de NDSI que l'on
#' considere pour identifier la fonte complète de la neige.
#' (valeur par défault est 0.15)
#' @param Fall.Snow.Thresold fraction de la valeur maximum de NDSI que l'on
#' considere à la fin de l'été pour identifier l'apparition de la neige.
#' (valeur par défault est 0.15)
#'
#' @return Retourne un tableau avec le jour de l'année de chaque evenements.
#' De plus, une figure en format png est creee dans le répertoire de travail pour
#' chaque année contenue dans le fichier CSV.
#'
#' @author Simon Bélanger
#'
#' @export
#'
extract.DOY.NDSI.phenology <- function(filen,
                                       Melt.Onset.Threshold=0.75,
                                       Snow.Free.Threshold=0.15,
                                       Fall.Snow.Thresold=0.15) {
  NDSI <- read.csv(filen, header = T)

  nyear=dim(NDSI)[2]-1
  cyears=colnames(NDSI)[2:(nyear+1)]
  fyears=as.numeric(str_sub(cyears, 2,5))

  Snow.Free<- rep(NA, nyear)
  Melt.Onset<- rep(NA, nyear)
  Fall.Snow<- rep(NA, nyear)
  W<- rep(NA, nyear)


  for(i in 1:nyear) {

    ## Normalisation de la série par le maximum atteint
    df = NDSI[,c(1,(i+1))]
    names(df) <- c("doy", "year")

    NDSI_n <-  (df$year - min(df$year,na.rm=T))/
      (max(df$year,na.rm=T) - min(df$year,na.rm=T))

    df = data.frame(doy=NDSI$doy, NDSI_n)
    names(df) <- c("doy", "year")
    smoothed <- loess(year ~ doy, data=df, span = 0.1)
    x  = predict(smoothed, 1:365)

    delta = (x[2:365] - x[1:364])
    delta = c(delta,NA)


    Melt.Onset[i] <- which(x < Melt.Onset.Threshold & delta < 0)[1]
    Snow.Free[i] <- which(x < Snow.Free.Threshold)[1]
    tmp <- which(x < Fall.Snow.Thresold)
    Fall.Snow[i] <- tmp[length(tmp)]

    png(file=paste("ndsi_",cyears[i],".png",sep=""),res=300,units="in", width=7, height=5)

    plot(df,
         xlab="Jour de l'annee",
         ylab="NDSI normalisé",
         pch=19,
         col="grey",
         main=paste("Série Temporelle: NDSI", fyears[i]))
    lines(1:Melt.Onset[i], x[1:Melt.Onset[i]], col="blue", lwd=5)
    lines(Melt.Onset[i]:Snow.Free[i], x[Melt.Onset[i]:Snow.Free[i]], col="lightblue", lwd=5)
    lines(Snow.Free[i]:Fall.Snow[i], x[Snow.Free[i]:Fall.Snow[i]], col="darkgreen", lwd=5)
    lines(Fall.Snow[i]:365, x[Fall.Snow[i]:365], col="blue", lwd=5)
    text(0,0.3, paste("Debut de la fonte:", Melt.Onset[i]),pos=4)
    text(0,0.2, paste("Duree de la fonte:", Snow.Free[i]-Melt.Onset[i]),pos=4)
    text(0,0.1, paste("Première neige d'automne:", Fall.Snow[i]),pos=4)
    text(0,0., paste("Duree de la saison sans neige:", Fall.Snow[i]-Snow.Free[i]),pos=4)
    dev.off()
  }
  DF  <- Snow.Free-Melt.Onset
  DSSN  <- Fall.Snow-Snow.Free

  return(data.frame(year=fyears,Melt.Onset, Snow.Free, Fall.Snow, DF, DSSN))

}

