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
###########################################################################################

#' Wrapper for fitting monotonic polynomial (MP) models with OpenMx.
#'
#' @param dat The data, typically in a format prepared by \code{\link[OpenMx]{mxFactor}}.
#' @param k A vector of integers (greater than or equal to 0) that controls the order of the polynomial for each item. Polynomial order is equal to 2*k+1.
#' @param fit Logical value. If true, actually fit the model. Otherwise, just return everything prepared for OpenMx (useful for inspecting set-up prior to running).
#' @param itemtype Character vector of same length as number of items that may allow different types of item on a single test. These currently support the following
#'   monotonic polynomial models from \code{\link[rpf]{rpf}}.
#'   1. Logistic function of monotonic polynomial ("lmp"). Dichotomous items, no asymptote. Uses \code{\link[rpf]{rpf.lmp}}).
#'   2. Monotonic Polynomial Graded Response Model ("grmp"). Based on graded response model. Uses \code{\link[rpf]{rpf.grmp}})
#'   3. Monotonic Polynomial Generalized Partial Credit Model ("gpcmp"). Based on generalized partial credit model model. Uses \code{\link[rpf]{rpf.gpcmp}})
#'   If an option is not specified, it tries to auto-detect dichotomous vs. polytomous items.
#' @param ncat Number of categories per item; if NULL, tries to auto-detect.
#' @param priors Logical value. If true, use prior distributions for alpha and tau parameters.
#' @param randstart Logical value. If true, tries to use random starting values (experimental).
#' @param startimat Accepts an item parameter matrix for custom starting values. Custom starting values are possible, but currently I only support custom starting values based on an item parameter matrix in the format
#'   of OpenMx. This is useful if we have already fitted a model and wish to increase/decrease polynomial order for just an item or two.
#'   We can use estimates from the initially fitted model for all other items, and any corresponding item parameters.
#'   To see an example of what such a matrix looks like, fit a model with this function, then extract item paramter matrix to see its format.
#' @param pvar Single numeric value indicating prior variance for alpha and tau parameters. Uses normal prior. Used only if \code{pvartau} or \code{pvaralpha} are NULL.
#' @param pvartau Single numeric value indicating prior variance for tau parameters. Uses normal prior.
#' @param pvaralpha Single numeric value indicating prior variance for tau parameters. Uses normal prior.
#' @param taumean Single numeric value indicating prior mean for tau parameters. -10 apparently works well in practice and I believe this is what Falk & Cai (2016, Psychometrika) used.
#' @param alphamean Single numeric value indicating prior mean for alpha parameters. Defaults to 0.
#' @param qpoints Integer indicating number of quadrature points - passed to \code{\link[OpenMx]{mxExpectationBA81}}.
#' @param qwidth Defines limits of quadrature grid. I think "5" here means -5 to 5. See \code{\link[OpenMx]{mxExpectationBA81}} documentation.
#' @param se Logical value. If TRUE, try to compute standard errors for item parameters.
#' @param infotype String passed as \code{information} argument to \code{\link[OpenMx]{mxComputeEM}} to determine how to compute information matrix (for standard errors).
#' @param semMethod If "mr1991" is chosen for \code{infotype}, then supplemented EM is used. This argument then takes a string that determines
#'  which variant of S-EM is performed. e.g., "mr" = as applied by Cai (2008) to IFA models, "tian" = as specified by Tian, Cai, Thissen, & Xin (2013), "agile" = Joshua Pritikin's method for S-EM.
#' @param ... Not used yet, but some arguments may later be passed directly to some OpenMx functions.
#' @details Setting up monotonic polynomial models in OpenMx can be a bit of a pain. This wrapper function
#'  attempts to make it easier. In a single line, it will set up and estimate a single group, unidimensional item response model with three choices
#'  for item models. Note that the models for each item are extensions of the two-parameter logistic, generalized partial credit, and graded response models.
#'  Setting \code{k} equal to 0 will essentially be equivalent to these simpler models, provided that all items are keyed the same direction (the GRMP can relax this assumption).
#'
#'  Inspiration for starting values is generally borrowed from the \code{MonoPoly} package and Elphinstone parameterization.
#'  The estimation procedure is generally that used by Falk & Cai (2016) and Falk (2020) in that normal priors for alpha and tau parameters typically help to stabilize estimation.
#'  Many estimation options are otherwise still hard-coded, but are those that we have found to be useful in past research.
#'  Note finally that estimating models with \code{k} greater than zero is not recommended until a simpler model is fit, as immediately fitting high order polynomial models without good starting values may lead to estimation difficulty.
#'
#'  Attempts at further generalization of this code are welcome.
#' @references Falk, C. F., & Cai, L. (2016). Maximum marginal likelihood
#' estimation of a monotonic polynomial generalized partial credit model with
#' applications to multiple group analysis. \emph{Psychometrika, 81}, 434-460.
#' \url{http://dx.doi.org/10.1007/s11336-014-9428-7}
#'
#' Falk, C. F. (2020). The monotonic polynomial
#' graded response model: Implementation and a comparative study. \emph{Applied Psychological
#' Measurement, 44}, 465-481. \url{https://doi.org/10.1177/0146621620909897}
#' @examples
#' \donttest{
#' # TODO examples here
#' }
#' @return Returns an object of class \code{\link[OpenMx]{MxModel-class}}. Convenience functions for extracting more from this object are available.
#' @seealso Anything else?
#' @export
#' @importFrom rpf rpf.lmp rpf.grmp rpf.gpcmp rpf.rparam
#' @importFrom ifaTools univariatePrior
#' @import OpenMx
#' @importFrom stats rnorm
fitMP<-function(dat,k=rep(0,ncol(dat)),fit=TRUE,itemtype=NULL,ncat=NULL,
                           priors=FALSE,randstart=FALSE, startimat=NULL,
                           pvar=500,pvartau=NULL,pvaralpha=NULL,taumean=-10,alphamean=0,
                           qpoints=101,qwidth=5,
                           se=FALSE,infotype="oakes1999",semMethod=NULL,...){

  # Process some input regarding priors
  pmalpha<-alphamean # prior mean for alpha
  pmtau<-taumean # prior mean for tau
  if(is.null(pvartau)){
    pvartau<-pvar
  }
  if(is.null(pvaralpha)){
    pvaralpha<-pvar
  }

  # Determine number of items and number of categories per item
  ni<-ncol(dat) # number of items
  if(is.null(ncat)){
    ncat<-vector("numeric")
    for(j in 1:ni){
      ncat<-c(ncat,sum(!is.na(unique(dat[,j]))))
    }
  }

  # Set up type of item
  if(is.null(itemtype)){
    itemtype<-vector("character")
    for(j in 1:ni){
      if(ncat[j]>2){
        itemtype<-c(itemtype,"grmp")
      } else {
        itemtype<-c(itemtype,"lmp")
      }
    }
  }

  npmax<-max(ncat) + 2*max(k) # max number of parameters per item
  kmax<-max(k)
  datlabs<-colnames(dat)
  nk<-unique(k)

  # create list of items - each slot has an object created from rpf
  spec<-list()
  for(j in 1:ni){
    if(itemtype[j]=="lmp"){
      spec[[j]]<-rpf.lmp(k[j])
    } else if (itemtype[j]=="grmp"){
      spec[[j]]<-rpf.grmp(ncat[j],k[j])
    } else if (itemtype[j]=="gpcmp"){
      spec[[j]]<-rpf.gpcmp(ncat[j],k[j])
    }
  }

  names(spec)<-datlabs

  # set up starting values; most of these are hard-coded; I suppose I could have added arguments to do these in a custom way
  alphastart<-0
  omegastart<- -.5
  lambdastart<- exp(-.5)
  startingValues <- mxSimplify2Array(lapply(spec, rpf.rparam))
  for(j in 1:ni){
    # polytomous items - intercept
    if(ncat[j]>2){
      xistart<-seq(1.5,-1.5,length.out=ncat[j]-1)
    } else {
      # dichotomous items intercept
      xistart<-0
    }

    ## Murray 2013 starting values for tau; actually from the MonoPoly package and Elphinstone parameterization
    alphataustart<-vector("numeric")
    if(kmax>0){
      taustart<-seq(.1,1,length=kmax)
      taustart<-log(taustart)
      for(q in 1:kmax){
        alphataustart<-c(alphataustart,0,taustart[q])
      }
    }

    # combine all starting values into a vector
    if(itemtype[j]!="grmp"){
      startvals<-c(omegastart,xistart,alphataustart)
    } else {
      startvals<-c(lambdastart,xistart,alphataustart)
    }

    # over-write randomly generated values from rpf.rparam and use these instead?
    if(randstart){
      startingValues[which(!is.na(startingValues[,j])),j]<-startvals[which(!is.na(startingValues[,j]))]+ rnorm(length(startvals[!is.na(startingValues[,j])]),0,.5)
    } else {
      startingValues[which(!is.na(startingValues[,j])),j]<-startvals[which(!is.na(startingValues[,j]))]#+rnorm(length(startvals),0,.1) # add a tiny bit of noise, otherwise
    }
  }

  rownames(startingValues) <- dimnames <- paste0('p', 1:nrow(startingValues))

  # create item parameter matrix
  imat <- mxMatrix(name='item', values=startingValues, free=!is.na(startingValues),
                   dimnames=list(dimnames,datlabs))

  # tweak parameter labels - useful for priors
  for(j in 1:ni){
    if(itemtype[j]!="grmp"){
      imat$labels[1,datlabs[j]]<-paste('omega',j,sep="")
    } else {
      imat$labels[1,datlabs[j]]<-paste('lambda',j,sep="")
    }
    for(i in 1:(ncat[j]-1)){
      imat$labels[i+1,datlabs[j]]<-paste(paste0("xi",i,"_"),j,sep="")
    }
    if(k[j]>0){
      for(i in 1:k[j]){
        imat$labels[ncat[j]+2*(i-1)+1,datlabs[j]]<-paste(paste("alpha",i,"_",sep=""),j,sep="")
        imat$labels[ncat[j]+2*(i-1)+2,datlabs[j]]<-paste(paste("tau",i,"_",sep=""),j,sep="")
      }
    }
  }

  # If a starting item parameter matrix is given as input, clobber relevant values in the item parameter matrix
  if(!is.null(startimat)){
    indx<-match(startimat$labels,imat$labels)
    imat$values[indx[!is.na(indx)]]<-startimat$values[!is.na(indx)]
  }

  # Create item model
  itemModel <- mxModel(model="itemModel", imat,
                       mxData(observed=dat, type="raw"),
                       mxExpectationBA81(spec,qwidth=qwidth,qpoints=qpoints),
                       mxFitFunctionML())

  # Compute sequence - here is where we may tweak whether standard errors are computed
  if(se){
    computeSeq<-mxComputeSequence(list(
      mxComputeEM('itemModel.expectation',
                  'scores',
                  mxComputeNewtonRaphson(maxIter=500L,tolerance=1e-9),
                  maxIter=2000L,tolerance=1e-7,
                  information=infotype, # mr1991 is S-EM
                  infoArgs=list(fitfunction='fitfunction',semMethod=semMethod)),
      mxComputeReportDeriv(),
      mxComputeStandardError()
    )) }else {
      computeSeq<-mxComputeSequence(list(
        mxComputeEM('itemModel.expectation',
                    'scores',
                    mxComputeNewtonRaphson(maxIter=500L,tolerance=1e-9),
                    maxIter=2000L,tolerance=1e-7)))
  }

  # If the largest k is greater than 0, set up (optional) priors and then set up the mxModel
  # Note that some estimation options are hard-coded. These seem to work well, though I suppose
  # some options could be passed along as arguments.
  if(kmax>0 & priors){

    ## priors
    priorAlphaLabels <- grep("alpha",imat$labels,value=TRUE)
    priorAlphaMode <- rep(NA, length(priorAlphaLabels))
    priorAlphaMode[1:length(priorAlphaLabels)] <- pmalpha
    priorAlphaModel <- univariatePrior('logit-norm',
                                       priorAlphaLabels, priorAlphaMode,
                                       strength=sqrt(pvaralpha),
                                       name="priorAlpha")

    priorTauLabels <- grep("tau",imat$labels,value=TRUE)
    priorTauMode <- rep(NA, length(priorTauLabels))
    priorTauMode[1:length(priorTauLabels)] <- pmtau
    priorTauModel <- univariatePrior('logit-norm',
                                     priorTauLabels, priorTauMode,
                                     strength=sqrt(pvartau),
                                     name="priorTau")
    fitfunc<-mxFitFunctionMultigroup(groups=c('itemModel.fitfunction',
                                              'priorAlpha.fitfunction',
                                              'priorTau.fitfunction'))

    Model <- mxModel(model="Model", itemModel, priorAlphaModel, priorTauModel,
                     fitfunc,computeSeq)
  } else {
    # otherwise set up model w/o priors
    fitfunc<-mxFitFunctionMultigroup(groups=c('itemModel.fitfunction'))
    Model <- mxModel(model="Model", itemModel,
                     fitfunc,computeSeq)
  }

  # fit model?
  if(fit){
    Model<-mxRun(Model)
  }

  return(Model)

}


###############################################################################################
##  Helper functions for extracting stuff from fitted OpenMx models
## This is annoying, I have to do these so often... Maybe I should write functions

#' This gets all item parameters (pars), number of categories (ncat), and k for an item
#' and returns all of these values in a list
#'
#' @param x - fitted model
#' @param j - item index
#' @param mp - whether the item in question is an MP model (otherwise, don't compute k)
#' @export
extractItemMx<-function(x,j,mp=TRUE){
  pars<-x$itemModel$item$values[,j]
  #ncat<-length(summary(x$itemModel$data$observed[,j]))
  #ncat<-length(unique(x$itemModel$data$observed[,j]))
  ncat<-sum(!is.na(unique(x$itemModel$data$observed[,j])))
  spec<-x@submodels$itemModel@expectation@ItemSpec[[j]]
  out<-list(pars=pars[!is.na(pars)],ncat=ncat,spec=spec)
  if(mp){
    k<-computek(x,j)
    out$k<-k
  }
  return(out)
}

#' Compute k (used in specifying polynomial order) for an item
#' This seems to work at auto-detecting k when the item is an MP item
#' I suppose if this is used on an item that isn't an MP item, it'll yield some value,
#' but that value sure won't be k?
#'
#' @param x - fitted mxModel or matrix that contains a "1" in some row indicating the value of k for each item
#' @param j - item index
#' @export
computek<-function(x,j){
  # If an mx Model, we can guess k from the number of parameters and categories
  if(class(x)=="MxModel"){
    pars<-x$itemModel$item$values[,j]
    npar<-length(pars[!is.na(pars)])
    #ncat<-length(summary(x$itemModel$data$observed[,j]))
    ncat<-sum(!is.na(unique(x$itemModel$data$observed[,j])))
    k<-(npar-ncat)/2
  } else if (class(x)=="matrix"){
    # if x is a matrix, determine k based on where there's a "1"
    k<-which(x[,j]==1)-1
  }
  return(k)
}

#' computes k for all items and returns k as a vector
#' @param x blah
#' @param ni blah
#' @export
computekrec<-function(x,ni){
  k<-vector("numeric")
  for(j in 1:ni){
    k<-c(k,computek(x,j))
  }
  return(k)
}

#' Extracts T(theta), or response functions, for an item from a fitted mxModel
#' and for all categories - one category per column
#'
#' @param x - fitted mxModel
#' @param j - item index
#' @param itemtype - type of item (lmp, grmp, gpcmp)
#' @param theta - grid for theta
#' @export

#' @importFrom rpf rpf.prob
traceProbMx<-function(x,j,itemtype="",theta=seq(-5,5,.1)){
  i<-extractItemMx(x,j)
  #if(itemtype=="lmp"){
  #  item<-rpf.lmp(i$k)
  #} else if (itemtype=="grmp"){
  #  item<-rpf.grmp(i$ncat,i$k)
  #} else if (itemtype=="gpcmp"){
  #  item<-rpf.gpcmp(i$ncat,i$k)
  #}
  item<-i$spec
  P<-t(rpf.prob(item,i$par,theta))
  return(P)
}

#' Computes item information
#' @param x - fitted mxModel
#' @param j - item index
#' @param theta - grid for theta
#' @importFrom rpf rpf.info
#' @export
infoMx<-function(x,j,theta=seq(-5,5,.1)){
  i<-extractItemMx(x,j)
  item<-i$spec
  info<-sapply(theta, function(th){
    rpf.info(item,i$par,th)
  })
  return(info)
}

#' Compute TCC
#' @param x blah
#' @param items blah
#' @param theta blah
#' @param scale blah
#' @param scaletest blah
#' @importFrom rpf rpf.prob
TCCMx<-function(x,items,theta=seq(-5,5,.1),scale=FALSE,scaletest=FALSE){
  P<-rep(0,length(theta))
  ni<-length(items)
  idx<-1
  for(j in items){
    i<-extractItemMx(x,j)
    item<-i$spec

    catmat<-matrix(0:(i$ncat-1),length(theta),i$ncat, byrow=TRUE)
    tmpP<-t(rpf.prob(item,i$par,theta))
    tmpP<-rowSums(catmat*tmpP)
    if(scale){
      tmpP<-tmpP/(i$ncat-1)
    }
    P<-P+tmpP
    idx<-idx+1
  }
  if(scaletest){
    P<-P/ni
  }
  return(P)
}

#' test information function
#' @param x blah
#' @param items blah
#' @param theta blah
#' @export
TIFMx<-function(x,items,theta=seq(-5,5,.1)){
  info<-rep(0,length(theta))
  ni<-length(items)
  idx<-1
  for(j in items){
    info<-info+infoMx(x,j,theta)
    idx<-idx+1
  }
  return(info)
}

#' Plots traceline
#' @param x - is an fitted mxModel
#' @param j - item index
#' @param itemtype - is the type of item; must be specified here
#' @param th grid
#' @param col color
#' @param ... blah
#' @importFrom graphics plot lines
#' @export
plotMx<-function(x,j,itemtype="",th=seq(-5,5,.1),col="black",...){
  P<-traceProbMx(x,j,itemtype,th)
  ncat<-ncol(P)
  plot(th,P[,1],type="l",ylim=c(0,1),col=col,...)
  for(j in 2:ncat){
    lines(th,P[,j],col=col)
  }
}

#' plots of info
#' @param x blah
#' @param j blah
#' @param th blah
#' @param ... blah
#' @importFrom graphics plot
#' @export
plotInfoMx<-function(x,j,th=seq(-5,5,.1),...){
  info<-infoMx(x,j,th)
  plot(th,info,type="l",...)
}

#' Function to create matrix that represents "k" for each item
#' @param kstart - starting value for k for all items
#' @param kmax - maximum k
#' @param ni - number of items
#' @export
newkmat<-function(kstart,kmax,ni){
  k.mat<-matrix(0,nrow=kmax+1,ncol=ni)
  k.mat[kstart+1,]<-1
  return(k.mat)
}

#' Compute EAP scores
#' @param x blah
#' @param items blah
#' @param dat blah
#' @export
#' @importFrom rpf EAPscores
EAPscoresMx<-function(x,items,dat){

  spec<-list()
  idx<-1
  for(j in 1:length(items)){
    i<-extractItemMx(x,items[j])
    #if(itemtype[j]=="lmp"){
    #  item<-rpf.lmp(i$k)
    #} else if (itemtype[j]=="grmp"){
    #  item<-rpf.grmp(i$ncat,i$k)
    #} else if (itemtype[j]=="gpcmp"){
    #  item<-rpf.gpcmp(i$ncat,i$k)
    #}
    spec[j]<-list(i$spec)
  }
  pars<-x$itemModel$item$values[,items,drop=FALSE]
  grp<-list(spec=spec,param=pars,data=dat)
  scores<-EAPscores(grp)
  return(scores)
}
