#'
#'@export
ASIRclass <- setClass(
  "ASIR",
  slots = c(
    Model = "ANY",
    MCMC = "ANY",
    Samples = "ANY"
  )
)
#'
#'@export
ASIR <- function(newR,
                N,
                t.step = 1,
                Frequency = TRUE,
                ChangePoint,
                TotalInfections){
  tempCode <- nimbleCode({
    # Set priors
    wGamma <- (DGamma*detectedInfections + UGamma*hiddenInfections)/(detectedInfections + hiddenInfections)
    Betas[1] ~ T(dnorm(mean = R0Means[1]*(wGamma*Pop^(Frequency == FALSE)), sd = R0SDs[1]), 0, Inf)
    Betas[2] ~ T(dnorm(mean = R0Means[2]*(wGamma*Pop^(Frequency == FALSE)), sd = R0SDs[2]), 0, Inf)
    UGamma ~ dgamma(shape = UGammaShape, rate = UGammaRate)
    DGamma ~ dgamma(shape = DGammaShape, rate = DGammaRate)
    # likelihood
    #S <- integer(length = Timepe)
    S[1] <- Pop - 1
    I[1] <- 1
    for(i in 1:(ChangePoint-1)){
      newI[i] ~ dbinom(size = S[i],
                       prob =  probGen(I[i]*Betas[1^(i < ChangePoint) + 1^(i>=ChangePoint)]*t.step/(Pop^Frequency)))
      newDR[i] ~ dbinom(size = I[i], prob =  probGen(DGamma*t.step))
      newUR[i] ~ dbinom(size = I[i] - newDR[i], prob =  probGen(UGamma*t.step))
      S[i+1] <- S[i] - newI[i]
      I[i+1] <- I[i] + newI[i] - newDR[i] - newUR[i]
    }
    for(i in ChangePoint:TimePeriod){
      newI[i] ~ dbinom(size = S[i],
                       prob =  probGen(I[i]*Betas[2]*t.step/(Pop^Frequency)))
      newDR[i] ~ dbinom(size = I[i], prob =  probGen(DGamma*t.step))
      newUR[i] ~ dbinom(size = I[i] - newDR[i], prob =  probGen(UGamma*t.step))
      S[i+1] <- S[i] - newI[i]
      I[i+1] <- I[i] + newI[i] - newDR[i] - newUR[i]
    }
  })
  hiddenInfections <- TotalInfections - sum(newR)
  return(ASIRclass(
    Model = compileNimble(
      nimbleModel(
        code = tempCode,
        constants = list(TimePeriod = length(newR),
                         ChangePoint = ChangePoint,
                         hiddenInfections = hiddenInfections,
                         detectedInfections = sum(newR)),
        data = list(newDR = newR,
                    t.step = t.step,
                    Pop = N,
                    Frequency = Frequency,
                    #priors
                    UGammaShape = 1,
                    UGammaRate = 1,
                    DGammaShape = 1,
                    DGammaRate = 1,
                    R0Means = rep(1, 2),
                    R0SDs = rep(1, 2)
                    ),
        inits = list(Betas = rep(1, 2),
                     UGamma = 1,
                     DGamma = 1,
                     newI = c(TotalInfections-1, rep(0, length(newR)-1)),
                     newUR = c(rep(0, length(newR)-1), hiddenInfections)
        ),
        calculate = FALSE
      )
    ),
    MCMC = NA,
    Samples = NA
  )
  )
}
#'
#'@export
initialValues.ASIR <- function(epiModel, hyperParameters){
  epiModel@Model$Betas <- hyperParameters$`Initial Values`$Betas
  epiModel@Model$UGamma <- hyperParameters$`Initial Values`$UGamma
  epiModel@Model$DGamma <- hyperParameters$`Initial Values`$DGamma
  epiModel@Model$UGammaShape <- hyperParameters$Priors$RecoveryRate$Shape
  epiModel@Model$UGammaRate <- hyperParameters$Priors$RecoveryRate$Rate
  epiModel@Model$DGammaShape <- hyperParameters$Priors$DetectionRate$Shape
  epiModel@Model$DGammaRate <- hyperParameters$Priors$DetectionRate$Rate
  epiModel@Model$R0Means <- hyperParameters$Priors$R0$Means
  epiModel@Model$R0SDs <- hyperParameters$Priors$R0$SDs

  mcmc <- configureMCMC(epiModel@Model, nodes = NULL)
  mcmc$addSampler(target = "newI",
                  type = sampler,
                  control = list(
                    TMax = 20,
                    DeltaMax = 20,
                    R = 1
                  ))
  mcmc$addSampler(target = "newUR",
                  type = sampler,
                  control = list(
                    TMax = 20,
                    DeltaMax = 20,
                    R = 1
                  ))
  mcmc <- buildMCMC(
    mcmc
  )
  mcmc <- compileNimble(mcmc, project = epiModel@Model, resetFunctions = TRUE)
  mcmc$run(hyperParameters$`Initial Values`$Runs)
  return(
    epiModel
  )
}
#'
#'@export
buildMCMCInternal.ASIR <- function(epiModel, hyperParameters){
  output <- configureMCMC(epiModel@Model, nodes = NULL)
  output$addSampler(target = c('Betas[1]', 'Betas[2]'),
                    type = sampler_RW_block,
                    control = hyperParameters[["RWM"]])
  output$addSampler(target = c('UGamma', 'DGamma'),
                    type = sampler_RW_block,
                    control = hyperParameters[["RWM"]])
  output$addSampler(target = "newI",
                    type = sampler,
                    control = hyperParameters[["N+-Delta"]])
  output$addSampler(target = "newUR",
                    type = sampler,
                    control = hyperParameters[["N+-Delta"]])
  output$addMonitors('Betas')
  output$addMonitors(c('UGamma', 'DGamma'))
  output <- buildMCMC(
    output
  )
  output <- compileNimble(output, project = epiModel@Model, resetFunctions = TRUE)
  return(output)
}