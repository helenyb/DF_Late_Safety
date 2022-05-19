#POMM
#simulations
#stage 1

#packages
library(doParallel)
registerDoParallel(cores=24)

#source code
source("POMM_function_stage2_comp-v7.R")


##define scenarios
#different dose levels MTD (linear around MTD)
s1<-c(0.30, 0.40, 0.50 ,0.60, 0.70, 0.80)
s2<-c(0.20, 0.30, 0.40, 0.50, 0.60, 0.70)
s3<-c(0.10, 0.20, 0.30, 0.40, 0.50, 0.60)
s4<-c(0.05, 0.10, 0.20, 0.30, 0.40, 0.50)
s5<-c(0.05, 0.10, 0.15, 0.20, 0.30 ,0.40)
s6<-c(0.02, 0.05, 0.10, 0.15, 0.20, 0.30)

#non-linear around MTD
s7<-c(0.15, 0.20, 0.25, 0.30, 0.45, 0.60)
s8<-c(0.05, 0.15, 0.30, 0.35, 0.40, 0.45)

#no dose has exactly target
s10<-c(0.05, 0.15, 0.25, 0.35, 0.45, 0.55)
s11<-c(0.15, 0.20, 0.35, 0.40, 0.45, 0.50)
s12<-c(0.05, 0.10, 0.15, 0.20, 0.25, 0.40)

#all unsafe (but closer to MTD than D)
s9<-c(0.40, 0.45, 0.50, 0.55, 0.60, 0.65)

#defined scenarios
s13<-c(0.06, 0.07, 0.08, 0.09, 0.11, 0.12) #A
s14<-c(0.10, 0.14, 0.21, 0.30, 0.46, 0.58) #B
s15<-c(0.16, 0.30, 0.50, 0.70, 0.89, 0.95) #C
s16<-c(0.55, 0.91, 0.99, 1.00, 1.00, 1.00) #D
s17<-c(0.05, 0.05, 0.05, 0.80, 0.80, 0.80) #E

scen_list<-lapply(c(1:17),function(x) get(paste(c("s",x),collapse="")))


##define fixed inputs
co_size<-3
ncohorts<-10
doses<-c(1.5,2.5,3.5,4.5,6,7)
target<-c(0.3,0.391)
ncycles<-3
cycle_dec<-1/3
nsims<-5000
ndoses<-length(doses)

##define hyper parameter values
priorDLT<-c(0.15,.2,.25,.3,.35,.4)



#number of pseudo patients per dose
n0<-2

grade_ratios<-c(0.2,0.3,0.4,0.5,0.6,0.7)




POMM_prior_data<-dousseau_prior_func2(prior_doses=priorDLT, num0=n0, graderatio_vec=grade_ratios, ncycles=ncycles, cycle_dec=cycle_dec) 





#execute calibration

for(sc in 1:length(scen_list)){
  
  
  cal_time <- system.time({
    assign(paste(c("POMM_simulation_s",sc),collapse = ""), foreach(i=1:nsims, combine = list) %dopar% {
      ##function
      POMM_stage2(seed=i,tru_cyc1=scen_list[[sc]],co_size=co_size,ncohorts=ncohorts,target=target,doses=doses,prior_pseudodata=POMM_prior_data,
                  nTTP_array=nTTP_array,ncycles=ncycles,cycle_dec=cycle_dec,
                  hard.safety.rule = 0, safety.stopping.low.unsafe= F, safety.stopping.high.toosafe = F,precision.stopping = F)
      
    })
  })
  save.image("POMM_simulations_stage1full.RData")
}

