#' Extrait le jour de l'année du début du bloom printanier,
#' du jour de l'atteinte du pic de Chlorophylle-a
#' et sa valeur à partir de données extraites de
#' Google Earth Engine générées avec la fonction
#' ui.Chart.image.doySeriesByYear.
#'
#' @param filen est le nom du fichier CVS telechargé de GEE et produit par
#' ui.Chart.image.doySeriesByYear. Le fichier peut contenir plusieurs années.
#' @param start est le jour de l'année correspondant au debut de la période qui nous interesse
#' (valeur par défault est 32 (1 fevrier)).
#' @param end est le jour de l'année correspondant à la fin de la période qui nous interesse
#' (valeur par défault est 182 (1 juillet))
#' @param growth.rate.threshold est le taux de croissance à partir duquel on considère
#' le début du bloom. (valeur par défault est 0.01 mg Chl/m^3/jour)
#'
#' @return Retourne un tableau avec le jour de l'année marquant le début du bloom (bloom.onset),
#' le jour où le pic de Chla est atteint (DOY.peak.chla), la valeur du pic (peak.chla). et la durée de la
#' croissance du bloom (DC).
#' De plus, une figure en format *.png est crée dans le repertoire de travail pour
#' chaque année contenu dans le fichier CSV.
#'
#' @author Simon Bélanger
#'
#' @export
extract.DOY.chla.phenology <- function(filen,
                                  start=32,
                                  end=182,
                                  growth.rate.threshold = 0.01,
                                  weekly=FALSE) {

  chl <- read.csv(filen, header = T)
  semaine <- week(as.Date(chl$doy-1, origin="1970-01-01"))

  nyear=dim(chl)[2]-1
  cyears=colnames(chl)[2:(nyear+1)]
  DOY.peak.chla <- rep(NA, nyear)
  peak.chla <- rep(NA, nyear)
  bloom.onset <- rep(NA, nyear)


  for(i in 1:nyear) {

    # Try to bin the data on a weekly basis
    if (weekly) {
      df = data.frame(semaine, chla=chl[,(i+1)])
      df <- df %>% group_by(semaine) %>% summarise(Chla.med = median(chla, na.rm=T))

      # # find the spring bloom onset, maximum Chla reach and the date
      smoothed <- loess(log10(Chla.med) ~ semaine, data=df, span = 0.25)
      x  = 10^predict(smoothed, 1:52)
      delta = (x[2:52] - x[1:51])
      bloom.onset[i] <- which(delta > growth.rate.threshold*7)[1]
      peak.chla[i] <- max(x[floor(start/7):round(end/7)], na.rm=T)
      DOY.peak.chla[i] <- which(x == max(x[floor(start/7):round(end/7)], na.rm=T))
      #
      png(file=paste("Chla_weekly_",cyears[i],".png",sep=""),res=300,units="in", width=7, height=5)
      myplot<-ggplot(aes(x = semaine), data=df) +
        geom_point(aes(y = Chla.med, color = cyears[i])) +
        geom_smooth(aes(y = Chla.med, color = cyears[i]), span=0.25, data=df) +
        geom_point(aes(x = DOY.peak.chla[i], y = x[DOY.peak.chla[i]], color="Chla peak"), size=3) +
        geom_point(aes(x = bloom.onset[i], y = x[bloom.onset[i]], color="Bloom onset"), size=3) +
        theme_bw() +
        xlab("Jour de l'année") +
        scale_y_log10() +
        ylab("Concentration de chlorophylle") +
        theme(legend.title = element_blank())
      print(myplot)
      dev.off()

      bloom.onset[i]   <- (bloom.onset[i]  *7) - 3 # convert week into DOY
      DOY.peak.chla[i] <- (DOY.peak.chla[i]*7) - 3


    } else {
      df = chl[,c(1,(i+1))]
      names(df) <- c("doy", "year")

      # find the spring bloom onset, maximum Chla reach and the date
      smoothed <- loess(log10(year) ~ doy, data=df, span = 0.25)
      x  = 10^predict(smoothed, 1:365)
      delta = (x[2:365] - x[1:364])
      bloom.onset[i] <- which(delta > growth.rate.threshold)[1]
      peak.chla[i] <- max(x[start:end], na.rm=T)
      DOY.peak.chla[i] <- which(x == max(x[start:end], na.rm=T))


      # plot the raw data and the smothed data using the loess function.
      png(file=paste("Chla_",cyears[i],".png",sep=""),res=300,units="in", width=7, height=5)
      myplot<-ggplot(aes(x = doy), data=df) +
        geom_point(aes(y = year, color = cyears[i])) +
        geom_smooth(aes(y = year, color = cyears[i]), span=0.25, data=df) +
        geom_point(aes(x = DOY.peak.chla[i], y = x[DOY.peak.chla[i]], color="Chla peak"), size=3) +
        geom_point(aes(x = bloom.onset[i], y = x[bloom.onset[i]], color="Bloom onset"), size=3) +
        theme_bw() +
        xlab("Jour de l'année") +
        scale_y_log10() +
        ylab("Concentration de chlorophylle") +
        theme(legend.title = element_blank())
      print(myplot)
      dev.off()
    }


  }

  DC <- (DOY.peak.chla - bloom.onset)
  return(data.frame(year=as.numeric(str_sub(cyears, 2,5)), bloom.onset, peak.chla, DOY.peak.chla, DC))
}
