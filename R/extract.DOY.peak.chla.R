#' Extrait le jour de l'année du pic de Chlorophylle-a
#' et sa valeur à partir de données extraites de
#' Google Earth Engine générées avec la fonction
#' ui.Chart.image.doySeriesByYear.
#'
#' @param filen est le nom du fichier CVS telechargé de GEE et produit par
#' ui.Chart.image.doySeriesByYear. Le fichier peut contenir plusieurs années.
#' @param start est le jour de l'année correspondant au debut de la période qui nous interesse
#' @param end est le jour de l'année correspondant à la fin de la période qui nous interesse
#'
#' @return Retourne un tableau avec le jour de l'année ou le pic de Chla
#' est atteint et la valeur du pic.
#' De plus, une figure en format *.png est crée dans le repertoire de travail pour
#' chaque année contenu dans le fichier CSV.
#'
#' @author Simon Bélanger
#'
#' @export
extract.DOY.peak.chla <- function(filen,
                                  start=30,
                                  end=180) {

  chl <- read.csv(filen, header = T)

  nyear=dim(chl)[2]-1
  cyears=colnames(chl)[2:(nyear+1)]
  DOY.peak.chla <- rep(NA, nyear)
  peak.chla <- rep(NA, nyear)


  for(i in 1:nyear) {
    df = chl[,c(1,(i+1))]
    names(df) <- c("doy", "year")
    smoothed_cha <- loess(log10(year) ~ doy, data=df[start:end,], span = 0.25)
    peak.chla[i] <- 10^max(smoothed_cha$fitted, na.rm = T)
    DOY.peak.chla[i] <- smoothed_cha$x[which(smoothed_cha$fitted == max(smoothed_cha$fitted, na.rm = T))]

    png(file=paste("Chla_",cyears[i],".png",sep=""),res=300,units="in", width=7, height=5)
    myplot<-ggplot(aes(x = doy), data=df) +
      geom_point(aes(y = year, color = cyears[i])) +
      geom_smooth(aes(y = year, color = cyears[i]), span=0.25, data=df) +
      theme_bw() +
      xlab("Jour de l'année") +
      scale_y_log10() +
      ylab("Concentration de chlorophylle") +
      theme(legend.title = element_blank())
    print(myplot)
    dev.off()
  }

  return(data.frame(year=as.numeric(str_sub(cyears, 2,5)),peak.chla, DOY.peak.chla))
}
