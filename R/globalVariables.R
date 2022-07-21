# Create environment ----
pkg.env <- new.env()


#Supported models ----

pkg.env$models=list(
  
  ## a model
  a = StMoMo::StMoMo(link="log",
                       staticAgeFun = TRUE,
                       periodAgeFun = NULL,
                       cohortAgeFun = NULL),
  ## p model
  
  p = StMoMo::StMoMo(link="log",
                        staticAgeFun = FALSE,
                        periodAgeFun = c("1"),
                        cohortAgeFun = NULL),
  ## c model
  c = StMoMo::StMoMo(link="log",
                        staticAgeFun = FALSE,
                        periodAgeFun = NULL,
                        cohortAgeFun = c("1")),
  ## ac model
  
  ac = StMoMo::StMoMo(link="log",
                     staticAgeFun = TRUE,
                     periodAgeFun = NULL,
                     cohortAgeFun = c("1")),
  
  ## ap model
  
  ap = StMoMo::StMoMo(link="log",
                     staticAgeFun = TRUE,
                     periodAgeFun = c("1"),
                     cohortAgeFun = NULL),
  
  
  ## pc model
  
  pc = StMoMo::StMoMo(link="log",
                     staticAgeFun = FALSE,
                     periodAgeFun = c("1"),
                     cohortAgeFun = c("1")),
  
  ## apc model
  
  apc = StMoMo::apc(),
  
  
  ## cbd model
  
  cbd = StMoMo::cbd(link = c("log")),
  
  ## m6 model
  
  m6 = StMoMo::m6(link = c("log")),
  
  ## m7 model
  
  m7 = StMoMo::m7(link = c("log"))
  
  
  )



# Utils ----
## Transform the upper run-off triangle to the life-table representation.
pkg.env$t2c <- function(x){
  "
  Transform the upper run-off triangle to the life-table representation.
  "
  I= dim(x)[1]
  J= dim(x)[2]
  
  mx=matrix(NA,nrow=I,ncol=J)
  for(i in 1:(I)){
    for(j in 1:(J)){
      if(i+j<=J+1){
        mx[j,(i+j-1)]=x[i,j]
      }
    }
  }
  return(mx)
}

## Rotate a matrix
pkg.env$rotate <- function(x){ 
  "
  It rotates a matrix by 90 degrees.
  "
  
  
  t(apply(x, 2, rev))}

pkg.env$c2t <- function(x){
  I= dim(x)[1]
  J= dim(x)[2]
  
  mx=matrix(NA,nrow=I,ncol=J)
  for(i in 1:(I)){
    for(j in 1:(J)){
      if(i+j<=J+1){
        mx[i,j]=x[j,(i+j-1)]
      }
    }
  }
  return(mx)
}


pkg.env$t2c.full.square <- function(x){
    "
  Function to transform a full run-off triangle into a square.
  
  This function takes a run-off triangle as input.
  
  It returns a square.
  "
    
    I= dim(x)[1]
    J= dim(x)[2]
    
    mx = matrix(NA, nrow = I, ncol = 2 * J)
    for (i in 1:(I)) {
      for (j in 1:(J)) {
        mx[j,(i+j-1)]=x[i,j]
        
      }
    }
    return(mx)
  }

pkg.env$fit.lc.nr <- function(data.T,
                              iter.max=1e+04,
                              tolerance.max=1e-06){
  "
  Fits the LC model via a Netwon-Rhapson from scratch implementation to the data.
  
  "
  
  data.O <- data.T$occurrance
  data.E <- data.T$exposure
  
  wxc <- matrix(1,
                nrow = dim(data.O)[1],
                ncol=dim(data.O)[2])
  wxc[is.na(data.E)]<-0
  
  
  ## parameters initial values
  ax = runif(dim(data.O)[1])
  bx = rep(1/dim(data.O)[1],dim(data.O)[1])
  kt = rep(0,dim(data.O)[2])
  
  #create matrices to compute the linear predictor
  ax.mx = matrix(rep(ax,dim(data.O)[2]),
                 byrow = F,
                 ncol=dim(data.O)[2])
  
  bx.mx = matrix(rep(bx,
                     dim(data.O)[2]),
                 byrow = F,
                 ncol=dim(data.O)[2])
  
  kt.mx = matrix(rep(kt,
                     dim(data.O)[1]),
                 byrow = T,
                 nrow=dim(data.O)[1])
  
  
  mu.mx = exp(ax.mx+bx.mx*kt.mx)
  
  dhat = data.E*mu.mx
  
  dxt0=1000
  niter=1000
  
  deltaax=100
  deltabx=100
  deltakt=100
  
  deltadxt=100
  
  kt0=(abs(kt))
  ax0=(abs(ax))
  bx0=(abs(bx))
  
  citer=0
  
  dxt.v=NULL
  
  deltakt.v=NULL
  deltaax.v=NULL
  deltabx.v=NULL
  
  deltadxt.v=NULL
  
  kt1.v=NULL
  
  converged=TRUE
  
  while((abs(deltaax)>tolerance.max|abs(deltabx)>tolerance.max|abs(deltakt)>tolerance.max)&citer<iter.max){
    #update alpha
    ax.num = apply((data.O-dhat)*wxc, 1, sum, na.rm=T)
    ax.den = apply(dhat*wxc, 1, sum, na.rm=T)
    ax = ax - ax.num/(-ax.den)
    
    deltaax=sum(abs(ax)-ax0)
    ax0=abs(ax)
    deltaax.v=c(deltaax.v,deltaax)
    
    ## compute the new estimates
    ax.mx = matrix(rep(ax,dim(data.O)[2]),
                   byrow = F,
                   ncol=dim(data.O)[2])
    mu.mx = exp(ax.mx+bx.mx*kt.mx)
    
    dhat = data.E*mu.mx
    
    #update gc
    kt.num = apply((data.O-dhat)*wxc*bx.mx, 2, sum, na.rm=T)
    kt.den = apply(dhat*wxc*(bx.mx^2), 2, sum, na.rm=T)
    kt = kt - kt.num/(-kt.den)
    kt[1]=0
    
    deltakt=sum(abs(kt)-kt0)
    kt0=(abs(kt))
    deltakt.v=c(deltakt.v,deltakt)
    
    kt1.v=c(kt1.v,kt[2])
    
    ## compute the new estimates
    
    kt.mx = matrix(rep(kt,
                       dim(data.O)[1]),
                   byrow = T,
                   nrow=dim(data.O)[1])
    
    mu.mx = exp(ax.mx+bx.mx*kt.mx)
    
    dhat = data.E*mu.mx
    
    #update beta
    bx.num = apply((data.O-dhat)*wxc*kt.mx, 1, sum, na.rm=T)
    bx.den = apply(dhat*wxc*(kt.mx^2), 1, sum, na.rm=T)
    bx = bx - bx.num/(-bx.den)
    bx=bx/sum(bx)
    
    deltabx=sum(abs(bx)-bx0)
    bx0=abs(bx)
    deltabx.v=c(deltabx.v,deltabx)
    
    ## compute the new estimates
    bx.mx = matrix(rep(bx,
                       dim(data.O)[2]),
                   byrow = F,
                   ncol=dim(data.O)[2])
    mu.mx = exp(ax.mx+bx.mx*kt.mx)
    
    dhat = data.E*mu.mx
    
    dxt1= sum(2*wxc*(data.O*log(data.O/dhat)-(data.O-dhat)),na.rm = T)
    deltadxt= dxt1-dxt0
    
    deltadxt.v=c(deltadxt.v,deltadxt)
    dxt.v=c(dxt.v,dxt1)
    
    dxt0=dxt1
    citer=citer+1
  }
  
  if((abs(deltaax)>tolerance.max|abs(deltabx)>tolerance.max|abs(deltakt)>tolerance.max)){
    warning('The Newton-Rhapson algorithm for the lee-carter may have not converged')
    converged=FALSE
  }
  
  return(list(ax=ax,
              bx=bx,
              kt=kt,
              citer=citer,
              data.T=data.T,
              converged=converged,
              deltaax=deltaax,
              deltabx=deltabx,
              deltakt=deltakt,
              dxt0=dxt0,
              deltadxt=deltadxt,
              converged=converged
              ))
  
}


pkg.env$fit.aac.nr <- function(data.T,
                              iter.max=1e+04,
                              tolerance.max=1e-06,
                              ax.start=NULL,
                              gc.start=NULL){
  "
  Fits the AAC model via a Netwon-Rhapson from scratch implementation to the data.
  
  "
  data.O <- pkg.env$rotate(pkg.env$c2t(data.T$occurrance))
  data.E <- pkg.env$rotate(pkg.env$c2t(data.T$exposure))
  
  wxc <- matrix(1,
                nrow = dim(data.O)[1],
                ncol=dim(data.O)[2])
  wxc[is.na(data.E)]<-0
  
  
  ## parameters initial values
  ax = ax.start
  bx = rep(1/dim(data.O)[1],dim(data.O)[1])
  gc = gc.start
  
  #create matrices to compute the linear predictor
  ax.mx = matrix(rep(ax,dim(data.O)[2]),
                 byrow = F,
                 ncol=dim(data.O)[2])
  
  bx.mx = matrix(rep(bx,
                     dim(data.O)[2]),
                 byrow = F,
                 ncol=dim(data.O)[2])
  
  gc.mx = matrix(rep(gc,
                     dim(data.O)[1]),
                 byrow = T,
                 nrow=dim(data.O)[1])
  
  
  mu.mx = exp(ax.mx+bx.mx*gc.mx)
  
  dhat = data.E*mu.mx
  
  dxt0=1000
  niter=1000
  
  deltaax=100
  deltabx=100
  deltagc=100
  
  deltadxt=100
  
  gc0=(abs(gc))
  ax0=(abs(ax))
  bx0=(abs(bx))
  
  citer=0
  
  dxt.v=NULL
  
  deltagc.v=NULL
  deltaax.v=NULL
  deltabx.v=NULL
  
  deltadxt.v=NULL
  
  gc1.v=NULL
  
  converged=TRUE
  
  while((abs(deltaax)>tolerance.max|abs(deltabx)>tolerance.max|abs(deltagc)>tolerance.max)&citer<iter.max){
    #update alpha
    ax.num = apply((data.O-dhat)*wxc, 1, sum, na.rm=T)
    ax.den = apply(dhat*wxc, 1, sum, na.rm=T)
    ax = ax - ax.num/(-ax.den)
    
    deltaax=sum(abs(ax)-ax0)
    ax0=abs(ax)
    deltaax.v=c(deltaax.v,deltaax)
    
    ## compute the new estimates
    ax.mx = matrix(rep(ax,dim(data.O)[2]),
                   byrow = F,
                   ncol=dim(data.O)[2])
    mu.mx = exp(ax.mx+bx.mx*gc.mx)
    
    dhat = data.E*mu.mx
    
    #update gc
    gc.num = apply((data.O-dhat)*wxc*bx.mx, 2, sum, na.rm=T)
    gc.den = apply(dhat*wxc*(bx.mx^2), 2, sum, na.rm=T)
    gc = gc - gc.num/(-gc.den)
    gc[1]=0
    
    deltagc=sum(abs(gc)-gc0)
    gc0=(abs(gc))
    deltagc.v=c(deltagc.v,deltagc)
    
    # gc1.v=c(gc1.v,gc[2])
    
    ## compute the new estimates
    
    gc.mx = matrix(rep(gc,
                       dim(data.O)[1]),
                   byrow = T,
                   nrow=dim(data.O)[1])
    
    mu.mx = exp(ax.mx+bx.mx*gc.mx)
    
    dhat = data.E*mu.mx
    
    #update beta
    bx.num = apply((data.O-dhat)*wxc*gc.mx, 1, sum, na.rm=T)
    bx.den = apply(dhat*wxc*(gc.mx^2), 1, sum, na.rm=T)
    bx = bx - bx.num/(-bx.den)
    bx=bx/sum(bx)
    
    deltabx=sum(abs(bx)-bx0)
    bx0=abs(bx)
    deltabx.v=c(deltabx.v,deltabx)
    
    ## compute the new estimates
    bx.mx = matrix(rep(bx,
                       dim(data.O)[2]),
                   byrow = F,
                   ncol=dim(data.O)[2])
    mu.mx = exp(ax.mx+bx.mx*gc.mx)
    
    dhat = data.E*mu.mx
    
    dxt1= sum(2*wxc*(data.O*log(data.O/dhat)-(data.O-dhat)),na.rm = T)
    deltadxt= dxt1-dxt0
    
    deltadxt.v=c(deltadxt.v,deltadxt)
    dxt.v=c(dxt.v,dxt1)
    
    dxt0=dxt1
    citer=citer+1
  }
  
  if((abs(deltaax)>tolerance.max|abs(deltabx)>tolerance.max|abs(deltagc)>tolerance.max)){
    warning('The Newton-Rhapson algorithm for the lee-carter may have not converged')
    converged=FALSE
  }
  
  return(list(ax=ax,
              bx=bx,
              gc=rev(gc), #reverse it as you read the data in the opposite direction
              citer=citer,
              converged=converged,
              deltaax=deltaax,
              deltabx=deltabx,
              deltagc=deltagc,
              dxt0=dxt0,
              deltadxt=deltadxt,
              converged=converged
  ))
  
}






