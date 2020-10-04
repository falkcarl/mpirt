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

# only RIMSE or RIMSD is usually scaled so that the resulting value is on a metric between 0 and 1
# Otherwise, wABC and sDTF are dependent on the number of score categories, but Chalmers mentions an effect size that's scaled
#' @importFrom stats dnorm
tracedifmod<-function(mod1,mod2,j1,j2,
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
bolt2002 <- function(item,npitem,startvals,nppars,theta,w,...){
  w<-w/sum(w)
  est<-nlminb(startvals, bolt2002.nll,
              grm.item=item,
              true.item=npitem,
              trupars=nppars,
              theta=theta,
              w=w,...)
  # could do some post-processing to make this easier to use
  return(est)
}

# only designed for GRM at the moment
#' @importFrom stats dnorm
#' @importFrom rpf rpf.grm
bolt2002.model <- function(x, items, theta=seq(-4,4,length.out=500), w=dnorm(theta)/sum(dnorm(theta)),...){

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
    est<-bolt2002(grm.item,mp.item,initpar,mp.par,theta,w,...)

    pars[[j]]<-est$par
    conv[[j]]<-est$conv
    gr.items[[j]]<-grm.item
    mp.items[[j]]<-mp.item
    mp.pars[[j]]<-mp.par

  }
  return(list(pars=pars,conv=conv,gr.items=gr.items,mp.items=mp.items,mp.pars=mp.pars))
}

# will only do 1 pair of items
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



