###
#
# calibration evaluation under a simple ground truth
#
# prior: theta~gamma(a,b)
# model: poisson(theta)
# target: rate
#
# Add log file to document all settings
###
library("here")
library("imabc")
library("epitools") #pois.exact
outDir<-file.path(here("OutputTmp"))

# Known True Parameter theta0
set.seed(20240507)
shape0<-0.2
rate0<-1
theta0<-shape0/rate0 #rgamma(n=1,rate=rate0,scale=scale0)

# target from gamma(theta0)
nDat<-100
x<-rpois(n=nDat,lambda=theta0)
print(tgt<-mean(x))
# sample size from posterior
nPost<-10000

### Set up for IMABC
# Define Priors
priors<-define_priors(rateP=add_prior(dist_base_name="gamma", shape=shape0,rate=rate0,min=0,max=Inf),nDatT1=add_prior(nDat))
# Define Target Functions
fn1<-function(rateP,nDatT1){
  mean(rpois(n=nDatT1,lambda=rateP))
}
fn <- function(rateP,nDatT1) {
  res <- c()
  res["T1"] <- fn1(rateP,nDatT1)
  return(res)
}
CI<-as.numeric(pois.exact(sum(x),pt=length(x),conf.level=0.99)[,c("lower","upper")])
targets<-define_targets(
  G1 = group_targets(
    T1=add_target(
      target=tgt,
      # ToDo: revisit this
      starting_range = c(0.25,4)*CI,
      stopping_range = c(1,1)*CI,
    )
  )
)
target_fun <- define_target_function(
  targets, priors, FUN = fn, use_seed = FALSE
)
calibration_results <-tryCatch( imabc(priors = priors,
                                      targets = targets,
                                      target_fun = target_fun,
                                      #seed = NULL,
                                      N_start = 300,
                                      N_centers = 5, # why is this 1? makes sense to have more?
                                      Center_n = 500,
                                      N_cov_points = 50,
                                      N_post = 200,
                                      verbose=FALSE
), error=function(e){
  message(paste("Error i tgt:",tgt,":\n"),e)
  e
})
if(sum(class(calibration_results)=="error")>0){
  nCatch<-nCatch+1
}else{
  postDist<-calibration_results$good_parm_draws[,c("rateP","sample_wt","iter")]
  postSmpl<-sample(x=postDist$rateP,size=nPost,replace=T,prob=postDist$sample_wt)
} 

postDist_df = as.data.frame(postDist)
write_csv(postDist_df, "C:/Users/prave/Downloads/IMABC_Loop_Results.csv")

weighted.mean(postDist$rateP, w=postDist$sample_wt) #weighted mean of rateP

hist(postSmpl)
postCDF<-ecdf(postSmpl)

hist(postDist$rateP)

## create violin plot of posterior distributions per true lambda
## true lambda = 0.2, 0.3, 0.4, 0.5
# read in aggregated results
violin_df = read_csv("IMABC_Loop_Results_2.csv")

ggplot(violin_df, aes(x = lambda_true, y = rateP, group = lambda_true)) +
  geom_violin(aes(fill = lambda_estimate), scale = "width")

#known posterior
shapeP<-shape0+sum(x)
rateP<-rate0+nDat #why did we add nDat here?

maxRate<-2
qDF<-seq(from=0,to=maxRate,by=0.02)
cDF_0<-pgamma(q=qDF,shape=shape0,rate=rate0) #cdf of gamma prior
dDF_0<-dgamma(x=qDF,shape=shape0,rate=rate0) # pdf of gamma prior

cDF_P<-pgamma(q=qDF,shape=shapeP,rate=rateP) # cdf of gamma posterior
dDF_P<-dgamma(x=qDF,shape=shapeP,rate=rateP)

maxD<-3 #max(setdiff(dDF_0,Inf),setdiff(dDF_P,Inf))

png(file.path("C:/Users/Manpreet Dhankhar/Documents","rateTgt_Params2grant.png"))
#png(file.path(outDir,"rateTgt_Params2.png"))
plot(x=qDF,y=cDF_0,type='n',ylab="Density",xlab="Poisson rate",lty=2,xlim=c(0,maxRate),ylim=c(0,maxD),xaxs="i",yaxs="i",main="Ground Truth (GT) Experiment")
#png(file.path(outDir,"rateTgt_zoom_Params2.png"))
#   plot(x=qDF,y=dDF_0,type='n',ylab="Target: Rate",xlab="Poisson rate",lty=2,xlim=theta0+c(-1,1)/5,ylim=c(0,maxD),xaxs="i",yaxs="i",main="Ground Truth (GT) Experiment")

# prior
lines(x=qDF,y=cDF_0,lty=2,col=1,lwd=2)
#lines(x=qDF,y=dDF_0,lty=2,col=1,lwd=2)
#posterior CDF from sample
lines(postCDF,lty=2,col=4,lwd=1,pch=19,cex=0.5)
#lines(density(postSmpl,adjust=5),lty=2,col=2,lwd=1,pch=19,cex=0.5)

# posterior - known closed form
lines(x=qDF,y=cDF_P,lty=3,col=1,lwd=2)
#lines(x=qDF,y=dDF_P,lty=3,col=1,lwd=2)
# #target
# abline(v=theta0,col=3,lwd=2)   
# abline(h=tgt,col=3,lwd=2)  

legend("bottomright",
       c(#paste0(
         "Prior CDF", #: Gamma(shape=",shape0,"; rate=",rate0,")"),
         "GT posterior CDF",
         "IMABC estimated posterior CDF"),
       #paste0("Target estimated from GT, n=",nDat),
       #paste0("GT rate: ",theta0)),
       lty=c(2,3,1,1,1),col=c(1,1,4,3,3),bg="white",cex=0.7)
dev.off()

ggpairs(postDist)




