library(nleqslv)


lFitGammaDistribution <- function(dVec, cEstMeth = "ML", bRemoveZerosMoM = TRUE) {

  # cEstMeth - decide which estimation method to use: "Maximum Likelihood" or
  # "Method of Moments"
  cEstMeth <- match.arg(arg = cEstMeth, choices = c("MoM", "ML"), several.ok = FALSE)

  # check if all observations are nonnegative doubles
  assertive::is_double(dVec)
  assertive::assert_all_are_non_negative(x = dVec)
  # if there is any zero value obs and the Max Likelihood method is used - remove the zero values
  bCheckZeros <- dVec == 0
  if (any(bCheckZeros) & cEstMeth == "ML") {
    warning("Zero values have been found in the dVec argument of the dFitGammaDistribution function",
            " and the maximum likelihood estimation method has been indicated: ",
            " removing the zeros from the vector of observations! ")
    dVec <- dVec[!bCheckZeros]
  }
  # if there is any zero value obs and MoM is picked and user wants to remove zeros - remove them
  if (any(bCheckZeros) & cEstMeth == "MoM" & bRemoveZerosMoM) {
    warning("Zero values have been found in the dVec argument of the dFitGammaDistribution function",
            " and the method of moments estimation method has been indicated and user wants to remove zeros: ",
            " removing the zeros from the vector of observations! ")
    dVec <- dVec[!bCheckZeros]
  }
  # if dVec is empty after removal of zeros - throw an error
  if (length(dVec) == 0) {
    stop("There are no observations left after the removal of the zeros in the ",
         "vector dVec --- halting function execution with error")
  }

  # ============================================================================
  # =====================       METHOD OF MOMENTS       ========================
  if (cEstMeth == "MoM") {
    # calculate the values of estimators
    dBetaRateHat <- mean(dVec)/var(dVec)
    dAlphaShapeHat <- mean(dVec)*dBetaRateHat
    # prepare function output
    lResult <- list(
      estimation_method = "MoM",
      alpha_shape_hat = dAlphaShapeHat,
      beta_rate_hat = dBetaRateHat
    )
  }
  # ============================================================================


  # ============================================================================
  # =======================    MAXIMUM LIKELOHOOD   ============================
  if (cEstMeth == "ML") {
    # calculate average of observations
    avg_X = mean(dVec)
    # calculate average of logarithms of observations' values
    avg_ln_X = mean(log(dVec))
    # define the nonlinear equation to be solved
    funcParamEq <- function(x, avg_X, avg_ln_X) {
      log(x) - log(avg_X) + avg_ln_X - digamma(x)
    }
    # iterate over the grids of the starting points and solve the nonlinear equation
    # for the shape parameter of gamma distribution for all of them
    dParamStartPointsGrid <- seq(0.1*avg_X, 10*avg_X, 0.1*avg_X)
    dtEstimates <- data.table::data.table(
      starting_point = dParamStartPointsGrid,
      alpha_shape_hat = vector(mode = "double", length = length(dParamStartPointsGrid)),
      beta_rate_hat = vector(mode = "double", length = length(dParamStartPointsGrid)),
      number_of_newton_steps = vector(mode = "double", length = length(dParamStartPointsGrid)),
      nleslv_termination_code = vector(mode = "double", length = length(dParamStartPointsGrid)),
      check_function_value_at_solution = vector(mode = "double", length = length(dParamStartPointsGrid))
    )
    for (k in 1:length(dParamStartPointsGrid)) {
      cIterPoint <- dParamStartPointsGrid[[k]]
      # find parameters estimates
      lSovlerResult <- nleqslv::nleqslv(
        x = cIterPoint, fn = funcParamEq,
        avg_X = avg_X, avg_ln_X = avg_ln_X,
        method = "Newton",
        control = list(maxit = 1000)
      )
      iNewtonStepsNumber <- lSovlerResult$iter
      dAlphaShapeHat <- lSovlerResult$x
      dBetaRateHat <- dAlphaShapeHat / avg_X
      iTerminationCode <- lSovlerResult$termcd
      dFunValueCheck <- funcParamEq(x = dAlphaShapeHat, avg_X = avg_X, avg_ln_X = avg_ln_X)
      # save parameters estimates
      dtEstimates[k, alpha_shape_hat := dAlphaShapeHat]
      dtEstimates[k, beta_rate_hat := dBetaRateHat]
      dtEstimates[k, number_of_newton_steps := iNewtonStepsNumber]
      dtEstimates[k, nleslv_termination_code := iTerminationCode]
      dtEstimates[k, check_function_value_at_solution := dFunValueCheck]
    }
    # prepare function output
    lResult <- list(
      estimation_method = "ML",
      estimates_table = dtEstimates
    )
  }
  # ============================================================================

  return(lResult)
}

# dVec <- rgamma(n = 150, shape = 10, scale = 5)
# res <- lFitGammaDistribution(dVec = dVec, cEstMeth = "ML")
# View(res$estimates_table)
