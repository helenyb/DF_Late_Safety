##R_mTPI2
##simulations
#STAGE 1

#packages
library(doParallel)
registerDoParallel(cores=24)

#source code
source("RTPI_function_stage2_comp-v5.R")

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
ndoses<-length(doses)
target<-0.391
ncycles<-3
cycle_dec<-1/3
nsims<-5000


eps1<-0.0391
  eps2<-0.1173

#execute simulations


for(sc in 1:length(scen_list)){
  
  
  cal_time <- system.time({
    assign(paste(c("R_mTPI2_simulation_s",sc),collapse = ""), foreach(i=1:nsims, combine = list) %dopar% {
      ##function
      R_TPI_stage2(seed=i,tru_cyc1=scen_list[[sc]],co_size=co_size,ncohorts=ncohorts,
                   eps1 =eps1,  eps2=eps2, 
                   target=target,doses=doses,nTTP_array=nTTP_array,ncycles=ncycles,cycle_dec=cycle_dec,
                   hard.safety.rule = 0, safety.stopping.low.unsafe= F, safety.stopping.high.toosafe = F)
      
    })
  })
  save.image("R_mTPI2_simulations_stage1full.RData")
}

