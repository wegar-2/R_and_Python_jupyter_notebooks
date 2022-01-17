library(okeanos.astro)

lData <- vector(mode = "list", length = 10)
names(lData) <- 2012:2021
for (iIterYear in names(lData)) {
  lData[[iIterYear]] <-
    okeanos.astro::dtDownloadGeomagneticIndicesObsTimeSeriesForYear(
      iYear = as.integer(iIterYear)
    )
}
dtGeomagData <- data.table::rbindlist(l = lData, use.names = TRUE)
data.table::fwrite(x = dtGeomagData, file = "GFZ_geomagnetic_data.csv")
