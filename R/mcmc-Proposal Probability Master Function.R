#' Master function for find the log proposals for each proposal function.
#' Sends the models and revelant parameters to the relevant function for the
#' proposal method given.
#' @export
logProposal <- function(epiModel1, epiModel2, proposal, hyperParameters){
  if(identical(proposal, independenceSampler)){
    return(logIndependenceSampler(epiModel2, hyperParameters$`Independence Sampler`))
  }
  else if(identical(proposal, randomWalk)){
    return(logRandomWalk(epiModel1, epiModel2))
  }
  else if(identical(proposal, stepProp)){
    return(logStep(epiModel1, epiModel2))
  }
  else if(identical(proposal, hybridProp)){
    return(logHybrid(epiModel1, epiModel2))
  }
}
