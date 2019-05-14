#' Extrait le jour de l'année du début du réchauffement
#' des eaux de surface, du pic de SST
#' et sa valeur à partir de donnees extraites de
#' Google Earth Engine générées avec la fonction
#' ui.Chart.image.doySeriesByYear.
#'
#' @param filen est le nom du fichier CVS telechargé de GEE et produit par
#' ui.Chart.image.doySeriesByYear. Le fichier peut contenir plusieurs annees.
#'
#' @return Retourne un tableau avec le jour de l'année ou le pic de Chla
#' est atteint et la valeur du pic.
#' De plus, une figure en format png est crée dans le répertoire de travail pour
#' chaque annee contenu dans le fichier CSV.
#'
#' @author Simon Bélanger
#'
#' @export
extract.DOY.peak.SST.warming.onset <- function(filen) {

  sst <- read.csv(filen, header = T)

  nyear=dim(sst)[2]-1
  cyears=colnames(sst)[2:(nyear+1)]
  DOY.peak.sst<- rep(NA, nyear)
  peak.sst<- rep(NA, nyear)
  warming.onset<- rep(NA, nyear)


  for(i in 1:nyear) {
    df = sst[,c(1,(i+1))]
    names(df) <- c("doy", "year")
    smoothed <- loess(year ~ doy, data=df, span = 0.25)
    x  = predict(smoothed, 1:365)

    DOY.peak.sst[i] <- which.max(x)
    peak.sst[i] <- df$year[DOY.peak.sst[i]]

    ## Get the date when the warming begin in the season.
    # compute the difference in temperature between 2 consecutive days
    delta = (x[2:365] - x[1:364])
    warming.onset[i] <- which(delta > 0)[1]



    png(file=paste("sst_",cyears[i],".png",sep=""),res=300,units="in", width=7, height=5)
    #plot(smoothed)
    #df %>%
    myplot<-ggplot(aes(x = doy), data=df) +
      geom_point(aes(y = year, color = cyears[i])) +
      geom_smooth(aes(y = year, color = cyears[i]), span=0.25, data=df) +
      geom_point(aes(x = DOY.peak.sst[i], y = x[DOY.peak.sst[i]], color="SST peak"), size=3) +
      geom_point(aes(x = warming.onset[i], y = x[warming.onset[i]], color="Warming onset"), size=3) +
      theme_bw() +
      xlab("Jour de l'année") +
      ylab("Température de surface de la mer (°C)") +
      theme(legend.title = element_blank())
    print(myplot)
    dev.off()
  }

  return(data.frame(year=as.numeric(str_sub(cyears, 2,5)),warming.onset,peak.sst, DOY.peak.sst))
}
