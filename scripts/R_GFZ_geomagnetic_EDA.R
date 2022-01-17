library(data.table)
library(tidyverse)
library(ggplot2)
library(stringr)
library(zoo)
library(lubridate)
library(scales)
library(broom)
source("scripts/gamma_fit_utils.R")


# --------- 1. load the data set ---------
dtGeomagneticData <- data.table::fread(file = file.path("data", "GFZ_geomagnetic_data.csv"))
# set key
print(data.table::key(dtGeomagneticData))

# --------- 2. drop redundant columns columns ---------
print(colnames(dtGeomagneticData))
cColsToDrop <- c(
  "days_since_origin", "days_since_origin_threehourly_interval"
)
dtGeomagneticData <- dtGeomagneticData[, (cColsToDrop) := NULL]
print(colnames(dtGeomagneticData))

# --------- 3. prepare date column ---------
# 3.1. pad months and days with zeros
print(lapply(X = dtGeomagneticData, FUN = typeof))
dtGeomagneticData[, date_month := stringr::str_pad(string = date_month, pad = "0", width = 2, side = "left")]
dtGeomagneticData[, date_day := stringr::str_pad(string = date_day, pad = "0", width = 2, side = "left")]
dtGeomagneticData[, date_year := as.character(date_year)]
# 3.2. concatenate date columns into a date string
dtGeomagneticData[, date_isodate := paste(date_year, date_month, date_day, sep = "-") %>%
                    as.Date(format = "%Y-%m-%d")]
# View(dtGeomagneticData)
# 3.3. drop redundant columns
cColsToDrop <- c("date_year", "date_month", "date_day")
dtGeomagneticData[, (cColsToDrop) := NULL]

# --------- 4. analyze obs times data ---------
table(dtGeomagneticData[, c("which_threehourly_interval", "specific_observation_time")])
# hence, column "which_threehourly_interval" can be dropped
dtGeomagneticData[, which_threehourly_interval := NULL]

# --------- 5. check if all observations are definitive ---------
# According to the description here: https://www-app3.gfz-potsdam.de/kp_index/Kp_ap_Ap_SN_F107_since_1932.txt
# any value of is_definitive column mean data is OK
table(dtGeomagneticData[, is_definitive])
print(nrow(dtGeomagneticData))
dtGeomagneticData <- dtGeomagneticData[dtGeomagneticData$is_definitive >= 1, ]
print(nrow(dtGeomagneticData))

# --------- 6. calculate quarterly descriptive statistics ---------
dtGeomagneticData[, quarter_date := zoo::as.yearqtr(date_isodate)]
dtSummaryStats <- dtGeomagneticData[
  , list(
    avg_Kp = mean(Kp), min_Kp = min(Kp), max_Kp = max(Kp),
    median_Kp = median(Kp),
    q_025_Kp = quantile(x = Kp, probs = 0.25),
    q_075_Kp = quantile(x = Kp, probs = 0.75),
    avg_ap = mean(ap), min_ap = min(ap), max_ap = max(ap),
    median_ap = median(ap),
    q_025_ap = quantile(x = ap, probs = 0.25),
    q_075_ap = quantile(x = ap, probs = 0.75)
  ),
  by = "quarter_date"
]
data.table::setkey(x = dtSummaryStats, "quarter_date")

# --------- 7. plot average and median Kp over time ---------
dtSummaryStatsPlot <- dtSummaryStats[, c("quarter_date", "avg_Kp", "median_Kp")]
dtSummaryStatsPlot <- data.table::melt.data.table(
  data = dtSummaryStatsPlot, id.vars = "quarter_date",
  measure.vars = c("avg_Kp", "median_Kp"),
  variable.name = "which_stat",
  value.name = "stat_value"
)
ggplot2::ggplot(data = dtSummaryStatsPlot) +
  ggplot2::geom_line(
    mapping = ggplot2::aes(x = quarter_date, y = stat_value, color = which_stat),
    size = 1
  ) +
  ggplot2::ggtitle(label = "Average & median of Kp index by quarter",
                   subtitle = "based on GFZ data for years 2012-2021") +
  ggplot2::xlab("Which quarter") +
  ggplot2::ylab("Average / median of Kp index") +
  ggplot2::theme_bw() +
  theme(
    plot.title = ggplot2::element_text(hjust = 0.5),
    plot.subtitle = ggplot2::element_text(hjust = 0.5),
    legend.position = "bottom"
  ) +
  ggplot2::scale_color_discrete(
    name = "Color of median & average lines"
  )


# --------- 8. scatter plot of ap against Kp index ---------
# Kp index is derived from ap using the rule described here: https://www.ngdc.noaa.gov/stp/GEOMAG/kp_ap.html
ggplotIndicesScatterPlot <- ggplot2::ggplot(
  data = dtGeomagneticData, mapping = ggplot2::aes(x = Kp, y = ap)
) + ggplot2::geom_point() + ggplot2::theme_bw() + ggplot2::xlab("Kp index") + ggplot2::ylab("ap index") +
  ggplot2::ggtitle(label = "Scatter plot of ap index vs Kp index",
                   subtitle = "Based on GFZ data for years 2012-2021") +
  ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),
                 plot.subtitle = ggplot2::element_text(hjust = 0.5)) +
  ggplot2::scale_x_continuous(breaks = seq(0, 9, 1), minor_breaks = seq(0, 9, 1/3)) +
  ggplot2::scale_y_continuous(breaks = seq(0, 300, 50), minor_breaks = seq(0, 300, 25))
ggplotIndicesScatterPlot

# --------- 9. histogram of the Kp index ---------
ggplotKpRootHistogram <- ggplot2::ggplot(
  data = dtGeomagneticData, mapping = aes(x = Kp)
) + geom_histogram(
  bins = length(unique(dtGeomagneticData$Kp)),
  color = "black", fill = "red", alpha = 0.2
) + theme_bw() +
  scale_x_continuous(
    breaks = seq(0, 9, 1),
    minor_breaks = seq(0, 9, 1/3),
    labels = scales::label_number(accuracy = 0.01)
  )
# 9.1. obs counts on the OY axis
ggplotKpAbsHistogram <- ggplotKpRootHistogram + scale_y_continuous(
  breaks = seq(0, 4000, 500),
  minor_breaks = seq(0, 4000, 250)
) + ggtitle(label = "Histogram of the Kp index",
            subtitle = "counts by Kp index value; GFZ data for years 2012-2021") +
  theme(
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5)
  ) + ylab(label = "Count of observations in a bin") +
  xlab(label = "Kp index value")
ggplotKpAbsHistogram
# 9.3. obs shares on the OY axis
ggplotKpRelHistogram <- ggplotKpRootHistogram +
  ggplot2::aes(y = ..count.. / sum(..count..)) +
  scale_y_continuous(
    breaks = seq(0, 0.15, 0.02), minor_breaks = seq(0, 0.15, 0.01),
    labels = scales::percent
  )+ ggtitle(label = "Histogram of the Kp index",
             subtitle = "shares in sample by Kp index value; GFZ data for years 2012-2021") +
  theme(
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5)
  ) + ylab(label = "Share of observations in a bin") +
  xlab(label = "Kp index value")
ggplotKpRelHistogram

# --------- 10. estimate gamma distribution parameters for Kp index ---------
# use helper function to fit the non-shifted gamma distribution
dGrid <- seq(min(dtGeomagneticData$Kp), max(dtGeomagneticData$Kp), length.out = 1000)
# 10.1. method of moments fit with zeros removed
lResult <- lFitGammaDistribution(dVec = dtGeomagneticData$Kp, cEstMeth = "MoM",
                                 bRemoveZerosMoM = TRUE)
dAlphaShapeHatMomNoZeros <- lResult$alpha_shape_hat
dBetaRateHatMomNoZeros <- lResult$beta_rate_hat
dGammaMomNoZerosDensity <- dgamma(x = dGrid, shape = dAlphaShapeHatMomNoZeros,
                                  rate = dBetaRateHatMomNoZeros)
# 10.2. method of moments fit with zeros NOT removed
lResult <- lFitGammaDistribution(dVec = dtGeomagneticData$Kp, cEstMeth = "MoM",
                                 bRemoveZerosMoM = FALSE)
dAlphaShapeHatMomZeros <- lResult$alpha_shape_hat
dBetaRateHatMomZeros <- lResult$beta_rate_hat
dGammaMomZerosDensity <- dgamma(x = dGrid, shape = dAlphaShapeHatMomZeros,
                                  rate = dBetaRateHatMomZeros)
# 10.3. maximum likelihood estimates
# 10.3.1. find parameters estimate
lResult <- lFitGammaDistribution(dVec = dtGeomagneticData$Kp, cEstMeth = "ML")
# 10.3.2. plot shape parameter estimate agains the nonlinear euqation solver's
#         starting points
dtEstimatesTable <- lResult$estimates_table
plotShapeParamEstimVsStartingPoint <- ggplot2::ggplot(
  data = dtEstimatesTable, mapping = aes(x = starting_point, y = alpha_shape_hat)
) + geom_line() + theme_bw() +
  ggtitle(
    label = "ML shape param estimate against starting point of solver"
  ) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_x_continuous(breaks = seq(0, 10, 2),
                     minor_breaks = seq(0, 20, 1))
plotShapeParamEstimVsStartingPoint
# load the parameters - given low std of estimates pick the first one - as good
# as any other
dAlphaShapeHatMaxLik <- dtEstimatesTable[1, alpha_shape_hat]
dBetaRateHatMaxLik <- dtEstimatesTable[1, beta_rate_hat]
dGammaMaxLikDensity <- dgamma(x = dGrid, shape = dAlphaShapeHatMaxLik,
                                rate = dBetaRateHatMaxLik)
# 10.4. prepare table with values of estimates
dtEstimatesSummary <- data.table::data.table(
  "Estimation method" = c("MoM, zeros not removed", "MoM, zeros removed", "ML"),
  "Shape parameter estimate" = c(dAlphaShapeHatMomZeros, dAlphaShapeHatMomNoZeros, dAlphaShapeHatMaxLik),
  "Rate parameter estimate" = c(dBetaRateHatMomZeros, dBetaRateHatMomNoZeros, dBetaRateHatMaxLik)
)
dtEstimatesSummary
# 10.5. make plot of densities fitted using various approaches to estimation of parameters
dtDensitiesData <- data.table::data.table(
  grid_point = dGrid,
  mom_no_zeros_density = dGammaMomNoZerosDensity,
  mom_zeros_density = dGammaMomZerosDensity,
  ml_density = dGammaMaxLikDensity
)
dtDensitiesData <- data.table::melt.data.table(
  data = dtDensitiesData,
  id.vars = "grid_point",
  measure.vars = c("mom_no_zeros_density", "mom_zeros_density", "ml_density"),
  variable.name = "which_density",
  value.name = "density_function_value"
)
View(dtDensitiesData)
ggplotComparativeDensitiesPlot <- ggplot2::ggplot(
  data = dtDensitiesData,
  mapping = aes(x = grid_point, y = density_function_value,
                color = which_density)
) + geom_line(size=1.1) + theme_bw() + theme(legend.position = "bottom") +
  ggtitle(label = "Density functions of Kp index for various estimation methods",
          subtitle = "Fitted for 3H freq data, using GFZ data for years 2012-2021") +
  scale_x_continuous(breaks = seq(0, 9, 1), minor_breaks = seq(0, 9, 0.5)) +
  scale_color_discrete(
    name = "Which method: ",
    labels = c("MoM, zero values not removed", "MoM, zero values removed", "MaxLik")
  ) + theme(plot.title = element_text(hjust=0.5), plot.subtitle = element_text(hjust = 0.5)) +
  xlab(label = "Kp index value") + ylab(label = "fitted gamma density function value")
ggplotComparativeDensitiesPlot

# --------- 10. fit non-linear regression to the ap vs Kp relationship ---------
# 10.1. fitting quadratic relationship: y = a*x^2 + b*x + c
dtRegData <- dtGeomagneticData[, c("Kp", "ap")]
objQuadraticRegModel <- lm(formula = ap ~ Kp + I(Kp^2), data = dtRegData)
# predict(nonlinRegModel)
# 10.2. fitting cubic relationship: y = a*x^2 + b*x + c
objCubicRegModel <- lm(formula = ap ~ Kp + I(Kp^2) + I(Kp^3), data = dtRegData)
# 10.3. compare the results of fitting both regressions
# 10.3.1. quadratic model
broom::glance(objQuadraticRegModel)
broom::augment(objQuadraticRegModel)
broom::tidy(objQuadraticRegModel)
# 10.3.2. cubic model fit
broom::glance(objCubicRegModel)
broom::augment(objCubicRegModel)
broom::tidy(objCubicRegModel)
# 10.4. make plot of the index with fitted curves:
dtQuadraticFit <- broom::augment(x = objQuadraticRegModel)
dtCubicFit <- broom::augment(x = objCubicRegModel)
dtFittedData <- data.table::data.table(
  Kp = dtQuadraticFit$Kp, ap = dtQuadraticFit$ap,
  quadratic_fit = dtQuadraticFit$.fitted,
  cubic_fit = dtCubicFit$.fitted
) %>% dplyr::mutate(ap = as.double(ap)) %>% as.data.table()

ggplot(data = dtFittedData) +
  geom_point(mapping = aes(x = Kp, y = ap, color = "blue"),  size = 2, show.legend = T) +
  geom_line(mapping = aes(x = Kp, y = quadratic_fit, color="grey"), size = 0.75, show.legend = T) +
  geom_line(mapping = aes(x = Kp, y = cubic_fit, color="black"), size = 0.75, show.legend = T) +
  theme_bw() +
  scale_x_continuous(breaks = seq(0, 9, 1), minor_breaks = seq(0, 9, 0.25)) +
  scale_y_continuous(breaks = seq(0, 300, 50), minor_breaks = seq(0, 300, 12.5)) +
  scale_color_manual(
    name = "Legend:", guide = "legend",
    values = c("blue", "grey", "black"),
    labels = c("data points", "quadratic fit", "cubic fit")
  ) + ggtitle(
    label = "Regression of ap index vs Kp index",
    subtitle = "based on GFZ data 2012-2021"
  ) + theme(
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5)
  ) + xlab(label = "Kp index") +
  ylab(label = "ap index, predicted values of quadratic, cubic regressions")
