library(plyr)
library(dplyr)
library(reshape2)
library(tidyr)
library(viridisLite)
library(viridis)
library(patchwork)
library(ggplot2)
library(deSolve)
library(pracma)				# contains sigmoid function
library(TrenchR)
library(rvmethod) #gaussian function
library(ggpubr)
library(scico)

#read in data and output
temps.all<- read.csv("processed_data/TempTimeSeries.csv")
PerfDat<- read.csv("processed_data/PerformanceData.csv")
out<- read.csv("out/out_fec.csv")
out.dr<- read.csv("out/out_dr.csv")

#performance metric
pms<- c("dr", "sur", "long", "fec")
pm.ind<- 4

#rename treatment in expt 3
t3<- c("22_0","22_5","22_9","22_13")
t3.lab<- c("00","05","09","13")
temps.all$treatment[temps.all$expt==3]<- t3.lab[match(temps.all$treatment[temps.all$expt==3],t3)]
PerfDat$treatment[PerfDat$expt==3]<- t3.lab[match(PerfDat$treatment[PerfDat$expt==3],t3)]

#drop field treatments
#PerfDat <- PerfDat[-which(PerfDat$population=="field"),]
#====================
#Constant rate TPCs
#English grain aphid, Sitobion avenae
#Zhao et al 2013
#https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/1365-2656.12196

#https://journals.biologists.com/jeb/article/218/14/2289/14375/Daily-temperature-extremes-play-an-important-role

#development rate
#use Ma et al 2015
dr= function(T, Tmax=32.947, a=0.137, b=4.514, c=7.267){ 
  d=exp(a*T)-exp(b-(Tmax-T)/c)
  d[d<0]<- 0
  return(d)
}

#survival
sur= function(T, a=345.51, b=0.35) {
  sur=1-T/(a-b*T^2)
  sur[T>30]<- 0
  return(sur)
}

long= function(T, a=32.73, b= -0.91) 
{
  lon= a + b*T
  lon[lon<0]<- 0
  return(lon)
}

#fedundity TPC (nymphs/adult)
fec= function(T, a= -69.1, b=12.49, c= -0.34){
  fec=a +b*T +c*T^2
  fec[fec<0]<- 0
  return(fec)
}

#time series estimates
perf.nodamage<- function(pm, series,scale)  {
  perf=NA
  for(i in 1:length(series)){
    if(pm==1) perf= dr(series[i])
    if(pm==2) perf= sur(series[i])
    if(pm==3) perf= long(series[i])
    if(pm==4) perf= fec(series[i])
    
    if(i==1) perf.all=perf*scale
    if(i>1) perf.all=c(perf.all, perf*scale)
  }
  return(perf.all)
}

#===================
#Survival curves
#from TDT analysis.Rmd

library(car)
library(MASS)

aphid<-read.table("/Users/lbuckley/Library/CloudStorage/GoogleDrive-lbuckley@uw.edu/My\ Drive/Buckley/Work/AphidThermalTol/data/aphidTDT/Thermal\ death\ time\ analysis/surv.txt", header = T)
aphid<-aphid%>%mutate(dur_h=dur/60)%>%within(.,spec<-factor(spec,levels = c("MD","SA","RP")),age<-factor(age,levels=c("2","6","12")))
aphid_high<-aphid%>%filter(temp>0)
aphid_low<-aphid%>%filter(temp<=0)

## Survival curve fitting to logistic model
#fit logistic regression model to the survival data and estimate the goodness of fit using McFadden’s pseudo-R squared

model1<-function(data){glm(cbind(surv,total-surv) ~ dur,data, family = binomial(link = logit))}
species = unique(as.character(aphid$spec))
ages = unique(aphid$age)
temps = unique(aphid$temp)
TB <- data.frame(species=character(),ages=integer(),temp=integer(),slope=numeric(),intercept=numeric(),p_value=numeric(),rsquare=numeric())
for(i in 1:length(species)){
  for(j in 1:length(ages)){
    for(k in 1:length(temps)){
      tmp = aphid%>%filter(spec == species[i]&age == ages[j]&temp == temps[k])
      if(dim(tmp)[1]!=0){
        fit = model1(tmp)
        fit1 = Anova(fit)
        p_value = summary(fit)$coef[2,4]
        slope=summary(fit)$coef[2,1]
        intercept=summary(fit)$coef[1,1]
        rsquare = (fit$null.deviance/-2 - fit$deviance/-2)/(fit$null.deviance/-2)
        fit.model=data.frame(species[i],ages[j],temps[k],slope,intercept,p_value,rsquare)
        names(fit.model)=c("species","ages","temp","slope","intercept","p_value","rsquare")
        TB=rbind(TB,fit.model)}
      
    }
  }
}

#------
##estimate Lt50 at each temperature from survival curve
###Lt50 at high temperatures

aphid_high<-aphid%>%filter(temp>0)

model2<-function(data){dose.p(glm(cbind(surv,total-surv) ~ dur,data, family = binomial(link=logit)),p=0.5)}

species = unique(as.character(aphid_high$spec))
ages = unique(aphid_high$age)
temps=unique(aphid_high$temp)
lt50_high <- data.frame(spec=character(),ages=integer(),temps=integer(),lt50=numeric(),SE=numeric())
for(i in 1:length(species)){
  for(j in 1:length(ages)){
    for(k in 1:length(temps)){
      
      data = aphid_high%>%filter(spec == species[i]&age == ages[j]&temp == temps[k])
      if(dim(data)[1]!=0){
        lt50<-model2(data)
        LT50_high=data.frame(species[i],ages[j],temps[k],round(lt50[1],2),round(unname(attributes(lt50)$SE[,1]),2))
        names(LT50_high)=c("spec","age","temp","lt50","SE")
        lt50_high=rbind(LT50_high,lt50_high)
      }
      
    }
  }
}

###Lt50 at low temperatures
aphid_low<-aphid%>%filter(temp<=0)
species = unique(as.character(aphid_low$spec))
ages = unique(aphid_low$age)
temps=unique(aphid_low$temp)
lt50_low <- data.frame(spec=character(),ages=integer(),temps=integer(),lt50=numeric(),SE=numeric())
for(i in 1:length(species)){
  for(j in 1:length(ages)){
    for(k in 1:length(temps)){
      
      data = aphid_low%>%filter(spec == species[i]&age == ages[j]&temp == temps[k])
      if(dim(data)[1]!=0){
        lt50<-model2(data)
        LT50_low=data.frame(species[i],ages[j],temps[k],round(lt50[1],2),round(unname(attributes(lt50)$SE[,1]),2))
        names(LT50_low)=c("spec","age","temp","lt50","SE")
        lt50_low=rbind(LT50_low,lt50_low)
      }
      
    }
  }
}

##TDT model
### TDT model at high temperatures

species = unique(as.character(lt50_high$spec))
ages = unique(lt50_high$age)
tdt_high<-data.frame(species=character(),age=integer(),slope=numeric(),z=numeric(),ctm=numeric(),ctm1h=numeric(),rsquare=numeric(),pvalue=numeric(),intercept=numeric(),type = character())

for(i in 1:length(species)){
  for(j in 1:length(ages)){
    
    tmp = lt50_high%>%filter(spec == species[i]&age == ages[j])
    fit = lm(log10(lt50/60)~temp,data=tmp)###ctmax 1h
    fit1 = lm(log(lt50,10)~temp,data=tmp)
    intercept = summary(fit1)$coefficient[1,1]
    slope = summary(fit1)$coefficient[2,1]
    ctm1h = -summary(fit)$coefficient[1,1]/summary(fit)$coefficient[2,1]
    ctmax=-summary(fit1)$coefficient[1,1]/summary(fit1)$coefficient[2,1]
    rsquare = summary(fit1)$r.squared
    pvalue = summary(fit1)$coefficient[2,4]
    z = -1/slope
    type = c("warm")
    p=data.frame(species[i],ages[j],slope,round(z,2),round(ctmax,2),round(ctm1h,2),rsquare,pvalue,intercept,type)
    names(p)=c("species","age","slope","z","ctm","ctm1h","rsquare","pvalue","intercept","type")
    tdt_high<-rbind(tdt_high,p)
  }
}

### TDT model at low temperatures

tdt_low<-data.frame(species=character(),age=integer(),slope=numeric(),z=numeric(),ctm=numeric(),ctm1h=numeric(),rsquare=numeric(),pvalue=numeric(),intercept=numeric(),type = character())

for(i in 1:length(species)){
  for(j in 1:length(ages)){
    
    tmp = lt50_low%>%filter(spec == species[i]&age == ages[j])
    fit = lm(log10(lt50/60)~temp,data=tmp)###ctmax 1h
    fit1 = lm(log10(lt50)~temp,data=tmp)
    intercept = summary(fit1)$coefficient[1,1]
    slope = summary(fit1)$coefficient[2,1]
    ctm1h =-summary(fit)$coefficient[1,1]/summary(fit)$coefficient[2,1]
    ctmax=-summary(fit1)$coefficient[1,1]/summary(fit1)$coefficient[2,1]
    rsquare = summary(fit1)$r.squared
    pvalue = summary(fit1)$coefficient[2,4]
    z = -1/slope
    type = c("cold")
    p=data.frame(species[i],ages[j],slope,round(z,2),round(ctmax,2),round(ctm1h,2),rsquare,pvalue,intercept,type)
    names(p)=c("species","age","slope","z","ctm","ctm1h","rsquare","pvalue","intercept","type")
    tdt_low<-rbind(tdt_low,p)
  }
}

#write.csv(rbind(tdt_high,tdt_low),"tdt.csv")


#===================
#Expt
#Expt 1: vary min
#Expt 2: vary max
#Expt 3: vary variance (expt 3, mild means)
#Expt 4: vary means and variance (expt 3, high means)
#Expt 5: vary duration of heatwave
#Expt 6: vary length of heatwave and timing 
#Expt 7: vary length of heatwave and timing 

