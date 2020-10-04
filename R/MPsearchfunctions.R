######################################################################
##  Copyright 2018-2020 Carl F. Falk
##
##  This program is free software: you can redistribute it and/or
##    modify it under the terms of the GNU General Public License as
##    published by the Free Software Foundation, either version 3 of
##    the License, or (at your option) any later version.
##
##    This program is distributed in the hope that it will be useful,
##    but WITHOUT ANY WARRANTY; without even the implied warranty of
##    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##    GNU General Public License for more details.
##    <http://www.gnu.org/licenses/>
##

# Modifies currently selected "k" for each item to produce a neighboring model
#
# x - is a matrix that has indicators that correspond to what value of k is currently selected for each item
#   The matrix is of dimensions (kmax+1) X ni. Where kmax is the maximum k allowed (user-specified) and ni is the number of items
#   An example of such a matrix will be given in the MPexamples.R file. In retrospect, maybe I shouldn't have written it this way.
#   Couldn't I just have a vector that indicates the value for "k"? Oh well.
# items - the most number of items to change at a time. e.g., 2 -> 1-2 items can be changed
# step - selects the most amount of change in k allowed (e.g., step 1 -> k=0 to k=1 allowed, but not k=0 to k=2)
#   Really only step=1 is probably reasonable. Larger steps might yield estimation problems.
computeNeighbor<-function(x,items=1,step=NULL){
  kmax<-nrow(x)-1 # detect max k allowed based on input matrix
  ni<-ncol(x) # detect
  items<-sample(1:items,1) # how many items to change
  i<-sample(1:ni,items) # randomly select those items

  # loop over selected items (which might be only 1 item if items == 1)
  # and change k for those items
  for(indx in i){
    k<-which(x[,indx]==1)
    newk<-1:(kmax+1)
    newk<-newk[-k]
    if(!is.null(step)){
      newk<-newk[newk >= k-step  & newk <= k+step]
    }
    if(length(newk)>1){
      newk<-sample(newk,1)
    } else {
      newk<-newk
    }
    x[,indx]<-0
    x[newk,indx]<-1
  }
  return(x)
}

# Computes "energy" of a given model, which amounts to fitting the model and extracting information criterion.
#
# dat - N X ni data matrix in format usable by OpenMx
# k.mat - matrix that represents possible values of "k" for each item
# aic - logical value whether to compute AIC (TRUE) or BIC (FALSE)
# itemtype - character vector indicating the type of items
# priors - logical value whether prior distributions are used on alpha and tau
# startimat - starting item parameter matrix (more lengthly description given above for fitMP)
# pvar - prior variance for alpha and tau
# taumean - prior mean for tau
# qpoints - number of quadrature points
# qwidth - width of quadrature grid; 5 means -5 to +5
#
# Returned value of AIC or BIC also as item parameter matrix "imat" as an attribute
energy<-function(dat,k.mat, type="aic", itemtype=NULL, priors=TRUE,startimat=NULL,
                 qpoints=101,qwidth=5,...){
  N<-nrow(dat)
  k<-vector("numeric")
  for(j in 1:ncol(k.mat)){
    k<-c(k,which(k.mat[,j]==1)-1)
  }
  if(type=="sic"){
    se<-TRUE
  } else {
    se<-FALSE
  }
  fitmodel<-fitMP(dat,k=k,fit=TRUE,itemtype=itemtype,priors=priors,
                             startimat=startimat,
                             qpoints=qpoints,qwidth=qwidth,se=se,...)

  #fitfunc<-fitmodel$output$fit-fitmodel$output$algebras$Model.fitfunction

  en<-getIC(fitmodel,type,N)

  attr(en,"imat")<-fitmodel$itemModel@matrices$item
  attr(en,"mod")<-fitmodel
  return(en)
}

# compute transition probability for model (e) vs. candidate model (eprime)
transition.prob<-function(e,eprime,temp,random=FALSE){

  if(eprime<e | random){
    tprob<-1
  } else {
    tprob<-exp(-(eprime-e)/temp)
  }
  return(tprob)
}

#' Extract information criterion or other info from a fitted model
#' @param x Fitted \code{mxModel}, e.g., from \code{\link{fitMP}}.
#' @param type String indicating information criterion to extract. See details.
#' @param N Sample size (used in BIC computations). Could be auto-detected, but not done yet.
#' @param usefitfunc Logical value. Toggles how to obtain fit function (log-likelihood). At some point how these are stored may have changed.
#' @param infotype If \code{type} is set to \code{"sic"}, this determines how the information matrix is computed, if not already available from the fitted model.
#'   May require model re-fitting.
#' @details Use of Bayesian priors for MP models sometimes complicates computation of some information criterion as typically done by some popular software packages,
#'   and as done by Mislevy (1986). Typically computation of AIC, BIC, and the log-likelihood is done by plugging in parameter estimates (e.g., based on the posterior mode)
#'   into the equation for the marginal log-likelihood, instead of computing the value of the log-posterior. As the latter is typically done by OpenMx, the former is
#'   done by this function.
#'
#'   Currently supported are AIC ("aic"), BIC ("bic"), log-likelihood ("ll"), number of parameters ("np"),
#'   and stochastic information criterion ("sic").
#'
#' @references Mislevy, R.J. (1986) Bayes modal estimation in item response models. Psychometrika 51, 177–195. \url{https://doi.org/10.1007/BF02293979}
#' @export
getIC<-function(x,type=c("aic","bic","sic","ll","np"),N,usefitfunc=FALSE,
                infotype="oakes1999"){

  type <- match.arg(type)

  fitfunc<-x$output$fit-x$output$algebras$Model.fitfunction
  if(type=="aic"){
    IC<-2*length(x$output$estimate) + x$output$algebras$itemModel.fitfunction
    if(usefitfunc){
      IC<-2*length(x$output$estimate) + fitfunc
    }
  } else if (type=="bic"){
    # Try BIC?
    IC<-log(N)*length(x$output$estimate) + x$output$algebras$itemModel.fitfunction
    if(usefitfunc){
      IC<-log(N)*length(x$output$estimate) + fitfunc
    }
  } else if (type=="sic"){
    # SIC requires information matrix
    # what if standard errors not available?
    if(is.null(x$output$hessian)){
      newcompute<-computeSeq<-mxComputeSequence(list(
        mxComputeEM('itemModel.expectation',
                    'scores',
                    mxComputeNewtonRaphson(maxIter=500L,tolerance=1e-9),
                    maxIter=2000L,tolerance=1e-7,
                    information=infotype,
                    infoArgs=list(fitfunction='fitfunction')),
        mxComputeReportDeriv(),
        mxComputeStandardError()))
      x$compute<-newcompute
      x<-mxRun(x)
      fitfunc<-x$output$fit-x$output$algebras$Model.fitfunction
    }

    penalty<-try(log(det(x$output$hessian/2)))
    IC<-fitfunc + penalty
    if(is.nan(IC)){
      IC<-.Machine$double.xmax
    }
  } else if (type=="ll"){
    if(usefitfunc){
      IC<- fitfunc
    } else {
      IC<-x$output$algebras$itemModel.fitfunction
    }
  } else if (type=="np"){
    IC<-length(x$output$estimate)
  }
  return(IC)
}

#' Main function that performs simulated annealing (SA)
#'
#' @param dat The data, typically in a format prepared by \code{\link[OpenMx]{mxFactor}}.
#' @param k.mat Starting matrix that represents k for each item (e.g., as generated by \code{\link{newkmat}}).
#' @param itermax Number of iterations for SA algorithm.
#' @param inittemp Staring temperature for SA algorithm.
#' @param type Which information criterion to use as the objective function to minimize. Anything supported by \code{\link{getIC}} is currently possible and "aic" or "bic" are typical choices.
#' @param random Logical value indicating whether to just automatically accept neighboring models (meaning not actually do SA).
#' @param startimat Starting item parameter matrix, if any. Used only for the first iteration.
#' @param step Maximum allowed change in k for candidate models. e.g., step 1 -> k=0 to k=1 allowed, but not k=0 to k=2. Smaller steps less likely to have estimation problems.
#' @param items How many items' k should be perturbed for computing a neighboring model?
#' @param temptype Method for temperature schedule decreases (see \code{\link{newTemp}} function).
#' @param ... Additional arguments to be passed to \code{\link{fitMP}}.
#' @details In monotonic polynomial (MP) models, each item may have a different polynomial order. This function performs simulated annealing to select polynomial order.
#'   Details of the exact implementation are given in Falk (2019). In brief, at a given iteration of the algorithm, a candidate model is generated by
#'   perturbing the polynomial order for a certain number of items. The candidate model is fit to the data and model fit is compared to the current model.
#'   If fit improves, the candidate model is accepted and will be the current model in the next iteration. If the candidate model's fit is worse,
#'   it may still be accepted with some non-zero probability that depends on the discrepancy in fit and the current "temperature" of the algorithm.
#'   High temperatures will result in a higher probability of accepting a poorer candidate model. In the algorithm, temperature declines towards zero over time, yet depends
#'   on the temperature schedule and number of iterations(see \code{\link{newTemp}}).
#'
#'   This function assumes that all items on the test are candidates for higher order polynomials; functionality to restrict the search space is possible, but not yet implemented.
#'   The only functionality for search space restriction is the highest k to consider for each item, as determined by the size of the matrix
#'   in \code{k.mat}. Note also that the starting polynomial order is also determined by this matrix.
#'
#' @return A list with the following elements
#' @slot k.mat A matrix of the same type as its input that encodes the polynomial order for each item.
#' @slot beste A recording of the best objective function encountered from the best fitting model.
#' @slot bestk A vector encoding the values of k for the best model.
#' @slot bestmod The best fitted \code{mxModel}.
#' @slot bestimat The item parameter matrix from the best fitted model (useful if one wants to do more iterations or re-fit the same model).
#'
#' @importFrom stats runif
#' @export
simAnneal<-function(dat,k.mat,itermax=1000,inittemp=10,type="aic",
                     random=FALSE,startimat=NULL,step=1,items=1,
                     temptype="straight",...){

  # compute "energy" for first state
  e<-energy(dat,k.mat,type=type,startimat=startimat,...)

  beste<-c(e)
  bestk<-k.mat
  bestmod<-attr(e,"mod")
  bestimat<-attr(e,"imat")

  for(i in 1:itermax){
    # compute "temperature"
    temp<-newTemp(inittemp, itermax, i, type=temptype)

    # compute neighbor
    k.mat.neighbor<-computeNeighbor(k.mat, items=items, step=step)

    # obtain "energy" of candidate
    startimat<-attr(e,"imat")
    eprime<-energy(dat,k.mat.neighbor,type=type,startimat=startimat,...)

    # compute transition probability
    tprob<-transition.prob(e,eprime,temp,random)

    # whether to accept
    r<-runif(1,0,1)
    if(r<tprob){
      k.mat<-k.mat.neighbor
      e<-eprime
    } else {
      k.mat<-k.mat
      e<-e
    }

    # Test whether this is best observed
    if(e<beste){
      beste<-e
      bestk<-k.mat
      bestmod<-attr(e,"mod")
      bestimat<-attr(e,"imat")
    }
  }

  # obtain best model, and return
  ret<-list(k.mat=k.mat,beste=c(beste),bestk=bestk,bestmod=bestmod,bestimat=bestimat)
  return(ret)
}

#' Computes algorithm "temperature" based on initial temperature, max number of iterations, current iteration, and last temperature
#'
#' @param inittemp Initial temperature.
#' @param itermax Maximum number of iterations.
#' @param i Current iteration
#' @param lasttemp Last temperature (i.e., when i == itermax )
#' @param type Type of temperature schedule
#   "straight" - used in initial simulations for IMPS paper (Falk, 2019)
#   "straight2" - used by Stander & Silverman (1994)
#   "logarithmic" - best performing logarithmic schedule by Stander & Silverman (1994)
#     oops, lasttemp can't be 0 if logarithmic schedule is used.
newTemp<-function(inittemp, itermax, i, lasttemp=.00001, type=c("straight","straight2","logarithmic")){
  type <- match.arg(type)
  if(type=="straight"){
    # Used for IMPS paper
    temp<-inittemp*((itermax-i)/itermax)
  } else if (type=="straight2"){
    # Stander & Silverman, 1994
    # What shows up in the paper from 1994 (just slightly different, depending on how iterations are numbered, otherwise the same)
    temp<-(lasttemp-inittemp)/(itermax-1)*(i -1)+inittemp
  } else if (type=="logarithmic"){
    if(lasttemp<=0){
      stop("Pick a value for lasttemp that's a tiny bit above zero for logarithmic temperature schedule")
    }
    # Also from Stander & Silverman (1994)
    temp<-(lasttemp*inittemp*(log(itermax+1)-log(2)))/(lasttemp*log(itermax+1)-inittemp*log(2)+(inittemp-lasttemp)*log(i+1))
  }
  return(temp)
}

# might not be tested yet
stepwise<-function(dat,kmax,type="aic",itemtype=NULL,priors=TRUE,randstart=FALSE,startimat=NULL,pvar=500,
                    taumean=-10,...){
  fitted.models<-0
  ni<-ncol(dat)
  k.mat<-matrix(0,nrow=kmax+1,ncol=ni)
  k.mat[1,]<-1 # initialize at k=0 for all items

  # compute energy for first state
  e<-energy(dat,k.mat,type=type,itemtype=itemtype,randstart=randstart,startimat=startimat,pvar=pvar,taumean=taumean,...)
  startimat<-attr(e,"imat")
  beste<-e
  bestk<-k.mat

  flag<-TRUE
  while(flag){
    aic.tmp<-vector("numeric")
    k.tmp<-list()
    imat.tmp<-list()
    indx<-1
    for(j in 1:ni){
      k.mat<-bestk
      # detect current k
      k<-which(bestk[,j]==1)

      # increase by 1
      knew<-k+1

      # fit model and save
      if(knew<=(kmax+1)){
        k.mat[,j]<-0
        k.mat[knew,j]<-1
        e<-energy(dat,k.mat,type=type,itemtype=itemtype,randstart=randstart,startimat=startimat,pvar=pvar,taumean=taumean,...)
        fitted.models<-fitted.models+1
        aic.tmp<-c(aic.tmp,e)
        k.tmp[[indx]]<-k.mat
        imat.tmp[[indx]]<-attr(e,"imat")
        indx<-indx+1
      }
    }

    # If one has better AIC than current, it's the new state
    if(length(aic.tmp)>0){
      best.indx<-which.min(aic.tmp)
      # Is the best criterion better than the current?
      if(aic.tmp[best.indx]<beste){
        beste<-aic.tmp[best.indx]
        startimat<-imat.tmp[[best.indx]]
        attr(beste,"imat")<-startimat
        bestk<-k.tmp[[best.indx]]
      } else {
        flag<-FALSE
      }
    } else {
      # if not, break out of loop
      flag<-FALSE
    }
  }

  ret<-list(beste=beste,bestk=bestk,fitted.models=fitted.models)
  return(ret)
}
