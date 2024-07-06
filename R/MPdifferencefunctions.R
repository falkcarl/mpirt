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

# RIMSD, RIMSE, and effect size measures to examine discrepancies among response functions

#' Root Integrated Mean Square Difference
#'
#' @param T1 Vector of values corresponding to traceline or information function 1.
#' @param T2 Vector of values corresponding to traceline or information function 2.
#' @param w Vector of weights; usually these should sum to 1.
#' @details Computes root integrated mean square difference (RIMSD) between \code{T1} and \code{T2}, with weights for integration
#' from \code{w}. The difference is squared, a weighted sum is computed based on \code{w} and then the square root is taken.
#'
#' @export
RIMSD<-function(T1,T2,w){
  #out<-sqrt(sum((T1-T2)^2*w)/sum(w))
  #w<-w/sum(w)
  out<-sqrt(sum((T1-T2)^2*w))
  out
}

#' Integrated Absolute Difference
#'
#' @param T1 Vector of values corresponding to traceline or information function 1.
#' @param T2 Vector of values corresponding to traceline or information function 2.
#' @param w Vector of weights; usually these should sum to 1.
#' @details Computes integrated absolute difference (IAD) between \code{T1} and \code{T2}, with weights for integration
#' from \code{w}. Note, this function is essentially the same as AAD from Bolt (2002), wABC from Edelen, Stuck, and Chandra (2015), and uDTF from Chalmers, Counsell, and Flora (2016).
#'
#' @export
IAD<- function(T1,T2,w){
  wABC(T1,T2,w)
}

# Edelen et al; also in Hansen et al
# similar to uDTF
wABC<-function(T1,T2,w){
  #w<-w/sum(w)
  out<-sum(abs(T1-T2)*w)
  out
}

#' Integrated Signed Difference
#'
#' @param T1 Vector of values corresponding to traceline or information function 1.
#' @param T2 Vector of values corresponding to traceline or information function 2.
#' @param w Vector of weights; usually these should sum to 1.
#' @details Computes integrated signed difference (ISD) between \code{T1} and \code{T2}, with weights for integration
#' from \code{w}. Note, this function is essentially the same as sDTF from Chalmers, Counsell, and Flora (2016).
#'
#' @export
# integrated signed difference
ISD<-function(T1,T2,w){
  sDTF(T1,T2,w)
}

# Chalmers et al
sDTF<-function(T1,T2,w){
  #w<-w/sum(w)
  out<-sum((T1-T2)*w)
  out
}

# some attempt at dEAP function; not yet used
dEAP<-function(mod1,mod2,j1,j2,cats){
  catmat<-as.data.frame(matrix(cats,length(cats),1, byrow=TRUE))
  catmat<-mxFactor(catmat,levels=0:length(cats))
  colnames(catmat)<-colnames(mod1$itemModel$item$values[,j1,drop=FALSE])
  EAP1<-EAPscoresMP(mod1,j1,catmat)
  EAP2<-EAPscoresMP(mod2,j2,catmat)

  EAP1[,1]-EAP2[,1] # in Edelen et al; In Hansen et al, these are further weighted by observed probabilities of response to each category
}


#' Difference in tracelines or TCC
#'
#' @param mod1 Fitted \code{mxModel}, e.g., from \code{\link{fitMP}}.
#' @param mod2 Fitted \code{mxModel}, e.g., from \code{\link{fitMP}}.
#' @param j1 Item index (or vector of indices in the case of TCC) for the first model.
#' @param j2 Item index (or vector of indices in the case of TCC) for the second model.
#' @param method Which type of discrepancy to compute? RIMSD, ISD, or IAD.
#' @param scale logical value to determine whether to scale each expected score value to 0-1. See \code{\link{TCC}}.
#' @param scaletest logical value to determine whether to scale the overall test by dividing by the number of items. See \code{\link{TCC}}.
#' @param theta Grid over which to compute the traceline or TCC.
#' @param w Vector of weights; usually these should sum to 1.
#' @details This is a convenience function which allows computation of RIMSD, ISD, or IAD directly from two fitted models
#'   and will do so on either single items from each or for TCCs if a vector of items from each is selected. Note that RIMSD
#'   is usually scaled so that the resulting value is on a metric between 0 and 1. Variants of ISD and IAD (wABC an sDTF) in previous
#'   papers were dependent on the number of score categories. However, Chalmers et al. mentions a variant that can serve as an effect
#'   size if further scaled.
#' @export
#'
#' @importFrom stats dnorm
tracedif<-function(mod1,mod2,j1,j2,
                  method=c("RIMSD","ISD","IAD"),scale=TRUE,scaletest=TRUE,
                  theta=seq(-5,5,.1),w=dnorm(theta)/sum(dnorm(theta))){

  #w<-w/sum(w)
  method <- match.arg(method)

  ni<-length(j1)

  ES1<-TCC(mod1,j1,theta,scale,scaletest)
  ES2<-TCC(mod2,j2,theta,scale,scaletest)

  if(method=="RIMSD"){
    out<-RIMSD(ES1,ES2,w)
  } else if (method=="ISD"){
    out<-ISD(ES1,ES2,w)
  } else if (method=="IAD"){
    out<-IAD(ES1,ES2,w)
  }
  out
}

#' Difference in information functions or TIF
#'
#' @param mod1 Fitted \code{mxModel}, e.g., from \code{\link{fitMP}}.
#' @param mod2 Fitted \code{mxModel}, e.g., from \code{\link{fitMP}}.
#' @param j1 Item index (or vector of indices in the case of TIF) for the first model.
#' @param j2 Item index (or vector of indices in the case of TIF) for the second model.
#' @param method Which type of discrepancy to compute? RIMSD, ISD, or IAD.
#' @param theta Grid over which to compute the information function or TIF.
#' @param w Vector of weights; usually these should sum to 1.
#' @details This is a convenience function which allows computation of RIMSD, ISD, or IAD for item or test information (TIF) directly from two fitted models.
#' @export
#'
#' @importFrom stats dnorm
#' @importFrom stats dnorm
infodif<-function(mod1,mod2,j1,j2,
                    method=c("RIMSD","IAD","ISD"),
                    theta=seq(-5,5,.1),w=dnorm(theta)/sum(dnorm(theta))){
  #w<-w/sum(w)
  method<-match.arg(method)
  ni<-length(j1)

  EIF1<-TIF(mod1,j1,theta)
  EIF2<-TIF(mod2,j2,theta)

  if(method=="RIMSD"){
    out<-RIMSD(EIF1,EIF2,w)
  } else if (method=="ISD"){
    out<-ISD(EIF1,EIF2,w)
  } else if (method=="IAD"){
    out<-IAD(EIF1,EIF2,w)
  }
  out
}

#' Difference in tracelines for just a pair of items
#'
#' @param item1 Item specification. e.g., from \code{rpf}.
#' @param item2 Item specification. e.g., from \code{rpf}.
#' @param par1 Item parameters for item 1.
#' @param par2 Item parameters for item 2.
#' @param method Which type of discrepancy to compute? RIMSD, ISD, or IAD.
#' @param scale logical value to determine whether to scale each expected score value to 0-1.
#' @param theta Grid over which to compute the traceline.
#' @param w Vector of weights; usually these should sum to 1.
#' @details This is a convenience function which allows computation of RIMSD, ISD, or IAD from two item specifications and their parameters.
#' @export
#'
#' @importFrom stats dnorm
tracedifitem<-function(item1,item2,par1,par2,
                       method=c("RIMSD","ISD","IAD"),scale=TRUE,
                       theta=seq(-5,5,.1),
                       w=dnorm(theta)/sum(dnorm(theta))){

  method<- match.arg(method)

  P1 <- rpf.prob(item1,par1,theta)
  P2 <- rpf.prob(item2,par2,theta)
  ncat<-nrow(P1)

  catmat<-matrix(0:(ncat-1),length(theta),ncat, byrow=TRUE)

  ES1<-rowSums(catmat*t(P1))
  ES2<-rowSums(catmat*t(P2))
  if(scale){
    ES1<-ES1/(ncat-1)
    ES2<-ES2/(ncat-1)
  }

  if(method=="RIMSD"){
    out<-RIMSD(ES1,ES2,w)
  } else if (method=="ISD"){
    out<-ISD(ES1,ES2,w)
  } else if (method=="IAD"){
    out<-IAD(ES1,ES2,w)
  }
  out

}

#' @importFrom rpf rpf.prob
bolt2002.nll<- function(x,grm.item,true.item,trupars,theta,w){

  P<-rpf.prob(grm.item,x,theta)
  Ptrue<-rpf.prob(true.item,trupars,theta)

  tmp<-colSums(Ptrue*log(P)+(1-Ptrue)*log(1-P))
  fun<-sum(tmp*w)

  return(-fun)
  # TO DO: implement derivatives

}

#' @importFrom stats nlminb
bolt2002 <- function(item,npitem,startvals,nppars,theta,w){
  w<-w/sum(w)
  est<-nlminb(startvals, bolt2002.nll,
              grm.item=item,
              true.item=npitem,
              trupars=nppars,
              theta=theta,
              w=w)
  # could do some post-processing to make this easier to use
  return(est)
}

#' Fit graded response model to MP model based on procedure by Bolt (2002)
#'
#' @param x Fitted \code{mxModel}, e.g., from \code{\link{simAnneal}}.
#' @param items Vector of indices indicating which items to perform procedure for.
#' @param theta Grid across theta over which to perform estimation for Bolt (2002)'s approach.
#' @param w Weights corresponding to \code{theta} for performing estimation for Bolt (2002)'s approach.
#' @details This function goes item by item for the model in \code{x} and uses the procedure by Bolt (2002) to fit
#'   a graded response model as best it can to each item's (semi-parametric) response functions.
#' @return A list with the following elements
#' @slot conv Convergence code for each item returned from \code{nlminb}.
#' @slot gr.items A list of item specifications from the graded model.
#' @slot pars A list of item parameters from the fitted graded model for each item.
#' @slot mp.items A list of item specifications from the MP model.
#' @slot mp.pars A list of item parameters from the MP model for each item.
#'
#' @export
#' @importFrom stats dnorm
#' @importFrom rpf rpf.grm
bolt2002Model <- function(x, items, theta=seq(-4,4,length.out=500), w=dnorm(theta)/sum(dnorm(theta))){

  pars<-list()
  conv<-list()
  gr.items<-list()
  mp.items<-list()
  mp.pars<-list()
  for(j in 1:length(items)){

    mp.i<-extractItem(x,items[j])
    mp.par<-mp.i$pars
    mp.item<-mp.i$spec
    ncat<-mp.item$spec[2]

    initpar<-c(1,seq(2,-2,length.out=ncat-1))
    grm.item<-rpf.grm(ncat)
    est<-bolt2002(grm.item,mp.item,initpar,mp.par,theta,w)

    pars[[j]]<-est$par
    conv[[j]]<-est$conv
    gr.items[[j]]<-grm.item
    mp.items[[j]]<-mp.item
    mp.pars[[j]]<-mp.par

  }
  return(list(conv=conv,gr.items=gr.items,pars=pars,mp.items=mp.items,mp.pars=mp.pars))
}

#' @importFrom rpf rpf.sample
bootWellsBolt<-function(i, boltmod, N, ni, mask, kmax=3, theta, w, ...){

  tmp<-rpf.sample(N,items=boltmod$gr.items,params=boltmod$pars) # we want to generate under the GRM

  tmp.miss<-tmp
  tmp.miss[mask]<-NA # make missing data in same pattern as actual missing data, if applicable

  kmat<-newkmat(0,kmax,ni) #k = 3 -> 7th order polynomial

  #test<-simAnneal(tmp.miss,kmat,itemtype=rep("grmp",ni),itermax=500,
  #                 inittemp=5,type="aic",pvar=500,taumean=-1,temptype="logarithmic",
  #                 items=1)

  test<-simAnneal(tmp.miss,kmat,itemtype=rep("grmp",ni),...)

  res.tmp<-bolt2002Model(test$bestmod, 1:ni, theta=theta, w=w)

  rmsd.wb.tmp<-sapply(1:ni, function(x){
    tracedifitem(res.tmp$gr.items[[x]],res.tmp$mp.items[[x]],res.tmp$pars[[x]],res.tmp$mp.pars[[x]],
                 method = "RIMSD", theta=theta, w=w)
  })

  rmsd.wb.tmp
}


#' Wells and Bolt procedure for identifying item misfit.
#'
#' @param MPmod Fitted \code{mxModel}, e.g., from \code{\link{simAnneal}}.
#' @param dat Original data.
#' @param nrep Number of replications for parametric bootstrap.
#' @param kmax Max for integer that controls polynomial order.
#' @param theta Grid across theta over which to perform estimation for Bolt (2002)'s approach.
#' @param w Weights corresponding to \code{esttheta} for performing estimation for Bolt (2002)'s approach.
#' @param parallel Method for parallel processing, if any. \code{"lapply"}, which does not actually do parallel processing or \code{"furrr"}, which
#'   uses the \code{furrr} package with the number of processing cores specified by \code{ncores}.
#' @param ncores Number of processing cores to use for parallel processing.
#' @param seed Integer for setting seed when doing multicore processing.
#' @param p.adjust.method Optional argument passed to \code{p.adjust} for making adjustments to p-values (e.g., Bonferroni, Bejamini-Hochberg).
#' @param ... Arguments passed to \code{simAnneal}.
#' @details THIS FUNCTION IS EXPERIMENTAL: additional functionality will be added, which may change the function
#'  signature and capabilities of the function. The function is also very COMPUTATIONALLY INTENSIVE (i.e., slow).
#'
#'  Procedure in Wells and Bolt (orig Douglas and Cohen) to identify item misfit. This requires that a non- or semi-parametric
#'  model is first fit to the data. So far, MP-based item models are supported and must be fit by a separate function
#'  (e.g., \code{\link{simAnneal}}). The procedure from Bolt (2002) is followed to obtain a parametric model for each item
#'  that bests fits each non- or semi-parametrically estimated response function. This typically involves fitting the
#'  parametric model item-by-item to the estimated response function from the non- or semi-parametric model.
#'  Currently the graded response model is used as the parametric model for this purpose, optimization or fitting is done
#'  over a grid for theta with weights from a standard normal distribution (as this is typical for calibration), and \code{nlminb}
#'  is used for model fitting with numerical derivatives (the form of the log-likelihood for fitting is given by Bolt, 2002).
#'  A discrepancy between the non/semi-parametric model and this newly fitted parametric model can then be determined,
#'  such as RIMSD (functionality for IAD or ISD may be forthcoming). This value represents the best that the parametric model
#'  can get to fitting the non/semi-parametric function, or in other words, how much the non/semi-parametric response function
#'  differs from this parametric model. However, the sampling distribution for this discrepancy is unknown.
#'
#'  To then obtain the sampling distribution for the discrepancy (currently RIMSD) and obtain a p-value, what resembles a parametric
#'  bootstrap is performed:
#'  1. Generate M datasets under the parametric model estimated using the procedure from Bolt (2002) as described above. Here, we
#'   generate data with characteristics (sample size and pattern of missing data) that are identical to the original dataset.
#'  2. Fit non- or semi-parametric approach to each simulated dataset. Here, we use the \code{\link{simAnneal}} with some defaults chosen
#'   to be computationally efficient and currently only the graded version of the MP model is possible (support may change in future
#'   versions of this function). The defaults currently are the following, and cannot yet be changed: \code{itermax=500,
#'   inittemp=5,type="aic",pvar=500,taumean=-1,temptype="logarithmic",items=1}.
#'  3. For each estimated model, the procedure by Bolt (2002) is again used to obtain a discrepancy value (RIMSD).
#'  4. Since computation of RIMSD here is essentially under the null hypothesis that the true response function is the parametric model,
#'   then the obtained p-value for each RIMSD can be obtained (i.e., how far in the tail is our observed value?).
#' @return A list with the following elements
#' @slot boltModel A list that essentially contains the information from the Bolt (2002) procedure. See \code{\link{bolt2002Model}}
#' @slot dif A vector that provides the discrepancy between the non/semi-parametric model from the Bolt (2002) procedure.
#' @slot pval A vector that provides p-values for the discrepancy measure.
#' @slot bootResults A matrix that contains the results of the discrepancy measure from the parametric bootstrap. Replications are rows, columns are items.
#' @examples
#' \donttest{
#'
#' # For now, just load something from mirt
#' #library(mirt)
#' data(Science)
#'
#' dat <- mxFactor(Science,levels=1:4)
#' safit <- simAnneal(dat, k.mat=newkmat(0,2,4),
#'                    itermax = 4*6,
#'                    inittemp = 5,
#'                    type = "aic",
#'                    step = 1,
#'                    items = 1,
#'                    temptype = "logarithmic",
#'                    itemtype=rep("grmp",4))
#'
#' samod <- safit$bestmod # best model according to SA
#'
#' getkrec(samod, 4) # value of k for each item from best model
#'
#' # Estimation settings similar to SA, but fewer iterations
#' # If generating under the graded model, fewer iterations should be necessary anyway
#' # 100 replications also probably not enough to get very accurate p-values
#' WB <- WellsBolt(samod, dat, nrep=100, kmax=2, seq(-4,4,length.out=81),
#'                parallel="furrr", ncores=2,
#'                itermax = 12,
#'                inittemp = 5,
#'                type = "aic",
#'                step = 1,
#'                items = 1,
#'                temptype = "logarithmic"
#'                )
#'
#' WB$dif # observed RMSD values
#' WB$pval # p-values
#'
#' }
#'
#' @export
#' @importFrom stats dnorm ecdf p.adjust
#' @importFrom future multisession multicore sequential plan
#' @importFrom furrr future_map furrr_options
#' @importFrom parallelly supportsMulticore
WellsBolt<-function(MPmod, dat, nrep=500, kmax=3,
                    theta=seq(-4,4,length.out=500), w=dnorm(theta)/sum(dnorm(theta)),
                    parallel=c("furrr","lapply"),ncores=2,seed=5224L,
                    p.adjust.method=NULL,
                    ...){

  parallel <- match.arg(parallel)

  ni <- ncol(dat)
  N <- nrow(dat)
  mask <- is.na(dat)

  boltmod<-bolt2002Model(MPmod, 1:ni, theta=theta, w=w)
  rmsd.wb<-sapply(1:ni, function(x){
    tracedifitem(boltmod$gr.items[[x]],boltmod$mp.items[[x]],boltmod$pars[[x]],boltmod$mp.pars[[x]],
                 method = "RIMSD", theta=theta, w=w)
  })

  if(parallel=="lapply"){
    set.seed(seed)
    stop("lapply not yet tested")
  } else if (parallel=="furrr"){
    #require(furrr)
    if(supportsMulticore()){
      plan(multicore, workers = ncores)
    } else {
      plan(multisession, workers = ncores)
    }

    rmsd.wb.boot<-future_map(1:nrep, bootWellsBolt,
                             boltmod = boltmod, N=N, ni=ni, mask=mask, kmax=kmax, theta=theta, w=w,
                             .progress=TRUE,.options=furrr_options(seed=seed), ...=...)
    plan(sequential)
  }

  rmsd.result<-do.call("rbind",rmsd.wb.boot)

  ecdf.rmsd.boot<-apply(rmsd.result,2,ecdf)
  pval.rmsd<-vector("numeric")
  for(j in 1:ni){
    pval.rmsd<-c(pval.rmsd,1-ecdf.rmsd.boot[[j]](rmsd.wb[j]))
  }

  if(!is.null(p.adjust.method)){
    pval.rmsd<-p.adjust(pval.rmsd, method=p.adjust.method)
  }

  out <- list (boltModel = boltmod,
               dif = rmsd.wb,
               pval = pval.rmsd,
               bootResults = rmsd.result)

  return(out)
}





