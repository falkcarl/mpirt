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
compute.neighbor<-function(x,items=1,step=NULL){
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
    #print(paste0("k =",k))
    #print(paste0("newk =",newk))
    x[,indx]<-0
    x[newk,indx]<-1
  }
  return(x)
}

# Computes "energy" of a given model, which here I have hard-coded as either AIC or BIC
#
# dat - N X ni data matrix
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
                #pvar=1000,taumean=-10,qpoints=101,qwidth=5){
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
                             #startimat=startimat,pvar=pvar,taumean=taumean,
                             startimat=startimat,
                             qpoints=qpoints,qwidth=qwidth,se=se,...)

  fitfunc<-fitmodel$output$fit-fitmodel$output$algebras$Model.fitfunction

  en<-getIC(fitmodel,type,N)

  #if(aic){
    #en<-2*length(fitmodel$output$estimate) + fitmodel$output$algebras$itemModel.fitfunction
  #  en<-2*length(fitmodel$output$estimate) + fitfunc
  #} else {
    # Try BIC?
    #en<-log(N)*length(fitmodel$output$estimate) + fitmodel$output$algebras$itemModel.fitfunction
  #  en<-log(N)*length(fitmodel$output$estimate) + fitfunc
  #}
  attr(en,"imat")<-fitmodel$itemModel@matrices$item
  attr(en,"mod")<-fitmodel
  return(en)
}

transition.prob<-function(e,eprime,temp,random=FALSE){

  if(eprime<e | random){
    tprob<-1
  } else {
    tprob<-exp(-(eprime-e)/temp)
  }
  return(tprob)
}

#' To extract an information criterion from a fitted model
#' @param fitmodel blah
#' @param type blah
#' @param N blah
#' @param usefitfunc blah
#' @param infotype blah
#' @export
getIC<-function(fitmodel,type="aic",N,usefitfunc=FALSE,
                infotype="oakes1999"){
  fitfunc<-fitmodel$output$fit-fitmodel$output$algebras$Model.fitfunction
  if(type=="aic"){
    IC<-2*length(fitmodel$output$estimate) + fitmodel$output$algebras$itemModel.fitfunction
    if(usefitfunc){
      IC<-2*length(fitmodel$output$estimate) + fitfunc
    }
  } else if (type=="bic"){
    # Try BIC?
    IC<-log(N)*length(fitmodel$output$estimate) + fitmodel$output$algebras$itemModel.fitfunction
    if(usefitfunc){
      IC<-log(N)*length(fitmodel$output$estimate) + fitfunc
    }
  } else if (type=="sic"){
    # SIC requires information matrix
    #browser()
    #IC<- fitfunc + log(det(fitmodel$output$hessian/2)) # implicitly multiplying by 2 already (b/c 1/2 times the thing on the right)
    #log(det(N*fitmodel$output$hessian/2)) # would result in larger penalty term, but does it matter?
    #IC<- fitfunc + log(det(N*fitmodel$output$hessian/2)) # or do we need to solve?
    #IC<- fitfunc + log(det(N*solve(fitmodel$output$hessian/2/N)))
    #browser()
    #log(det(N*solve(solve(fitmodel$output$hessian/2/N))))
    #IC<- fitfunc + log(det(fitmodel$output$hessian/2)) # looks correct? not as severe as BIC, but... what do we get for fit?
    #IC<- fitfunc + log(det(fitmodel$output$hessian/2*N))
    #penalty<-try(log(det(N*solve(solve(fitmodel$output$hessian/2/N)))))
    #print(penalty)

    # what if standard errors not available?
    if(is.null(fitmodel$output$hessian)){
      newcompute<-computeSeq<-mxComputeSequence(list(
        mxComputeEM('itemModel.expectation',
                    'scores',
                    mxComputeNewtonRaphson(maxIter=500L,tolerance=1e-9),
                    maxIter=2000L,tolerance=1e-7,
                    information=infotype,
                    infoArgs=list(fitfunction='fitfunction')),
        mxComputeReportDeriv(),
        mxComputeStandardError()))
      fitmodel$compute<-newcompute
      fitmodel<-mxRun(fitmodel)
      fitfunc<-fitmodel$output$fit-fitmodel$output$algebras$Model.fitfunction
    }

    penalty<-try(log(det(fitmodel$output$hessian/2)))
    #print(penalty)
    #browser()
    IC<-fitfunc + penalty
    if(is.nan(IC)){
      IC<-.Machine$double.xmax
    }
    #IC<- fitfunc + log(det(N*))
  } else if (type=="ll"){
    if(usefitfunc){
      IC<- fitfunc
    } else {
      IC<-fitmodel$output$algebras$itemModel.fitfunction
    }
  } else if (type=="np"){
    IC<-length(fitmodel$output$estimate)
  }
  #print(2*length(fitmodel$output$estimate))
  #print(log(N)*length(fitmodel$output$estimate))
  #print(log(det(N*solve(solve(fitmodel$output$hessian/2/N)))))
  #print(IC)
  return(IC)
}

#' Main function that performs simulated annealing
#'
#' @param dat - data matrix of item responses
#' @param k.mat - matrix that represents k for each item
#' @param itermax - maximum number of iterations
#' @param inittemp - staring temperature for SA
#' @param type - blah
#' @param random - logical value indicating whether to just automatically accept neighboring models (meaning not actually do SA)
#' @param priors - logical value indicating whether to use prior distributions on alpha and tau parameters
#' @param startimat - starting item parameter matrix, if any. For example, say we start with a model with k=1 for all items
#' @param step - allowed change in k for candidate models (see documentation for compute.neighbor function)
#' @param items - number of items for perturbing k (see documentation for compute.neighbor function)
#' @param pvar - prior variance
#' @param taumean - mean of prior distribution for tau parameters; currently all parameters will have the same prior mean
#' @param itemtype - character vector that has length same as number of items; indicates what type of model to use for each item
#' @param qpoints - number of quadrature points
#' @param qwidth - width for quadrature grid
#' @param temptype - method for temperature schedule decreases (see newtemp function)
#' @param ... blah
#' @importFrom stats runif
#' @export
sim.anneal<-function(dat,k.mat,itermax=1000,inittemp=10,type="aic",
                     #random=FALSE,priors=TRUE,startimat=NULL,step=1,items=1,pvar=1000,taumean=-10,
                     random=FALSE,priors=TRUE,startimat=NULL,step=1,items=1,pvar=1000,taumean=-10,
                     itemtype=NULL,qpoints=101,qwidth=5,
                     temptype="straight",...){

  # compute "energy" for first state
  #e<-energy(dat,k.mat,type=type,itemtype=itemtype,priors=priors,pvar=pvar,startimat=startimat,taumean=taumean,
  #          qpoints=qpoints,qwidth=qwidth)
  e<-energy(dat,k.mat,type=type,itemtype=itemtype,priors=priors,pvar=pvar,startimat=startimat,taumean=taumean,
            qpoints=qpoints,qwidth=qwidth,...)

  beste<-c(e)
  bestk<-k.mat
  bestmod<-attr(e,"mod")
  bestimat<-attr(e,"imat")

  for(i in 1:itermax){
    # compute "temperature"
    temp<-newtemp(inittemp, itermax, i, type=temptype)
    #temp<-inittemp*((itermax-i)/itermax)

    # compute neighbor
    k.mat.neighbor<-compute.neighbor(k.mat, items=items, step=step)

    # obtain "energy" of candidate
    startimat<-attr(e,"imat")
    #eprime<-energy(dat,k.mat.neighbor,type=type,itemtype=itemtype,priors=priors,pvar=pvar,startimat=startimat,taumean=taumean,
    #               qpoints=qpoints,qwidth=qwidth)
    eprime<-energy(dat,k.mat.neighbor,type=type,itemtype=itemtype,priors=priors,pvar=pvar,startimat=startimat,taumean=taumean,
                   qpoints=qpoints,qwidth=qwidth,...)

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

# Computes algorithm "temperature" based on initial temperature, max number of iterations, current iteration, and last temperature
#
# inittemp - initial temperature (e.g., 5)
# itermax - maximum number of iterations
# i - current iteration
# lasttemp - last temperature (i.e., when i == itermax )
# type - type of temperature schedule
#   "straight" - used in initial simulations for IMPS paper
#   "straight2" - used by Stander & Silverman (1994)
#   "logarithmic" - best performing logarithmic schedule by Stander & Silverman (1994)
#     oops, lasttemp can't be 0 if logarithmic schedule is used.
newtemp<-function(inittemp, itermax, i, lasttemp=.00001, type="straight"){
  if(type=="straight"){
    # what I used
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

  # compute AIC for first state
  #e<-energy(dat,k.mat,type=type,itemtype=itemtype,randstart=randstart,startimat=startimat,pvar=pvar,taumean=taumean)
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
      # Is the best AIC better than the current?
      if(aic.tmp[best.indx]<beste){
        beste<-aic.tmp[best.indx]
        #startimat<-attr(beste,"imat")
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
