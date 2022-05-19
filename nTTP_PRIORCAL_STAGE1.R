##prior calibration nTTP

#packages
library(doParallel)
registerDoParallel(cores=24)

#source code
source("nTTP_function_stage2_comp-v9.R")


##define scenarios
P_s1<-c(0.30, 0.40, 0.45 ,0.50, 0.55, 0.60)
P_s2<-c(0.05, 0.07, 0.10, 0.15, 0.2, 0.30)
P_s3<-c(0.10, 0.20, 0.30,  0.40,  0.50, 0.60)
P_s4<-c(0.15, 0.20, 0.25, 0.30,  0.35,  0.40)
P_s5<-c(0.40, 0.45, 0.5, 0.55, 0.6, 0.65)
P_s6<-c(0.07, 0.09, 0.11, 0.13, 0.15, 0.17)


prior_cal_scen_list<-list(P_s1,P_s2,P_s3,P_s4,P_s5,P_s6)

##define fixed inputs
co_size<-3
ncohorts<-10
doses<-c(1.5,2.5,3.5,4.5,6,7)
target<-0.289
ncycles<-3
cycle_dec<-1/3



##define grid for parameter values
doseparmean<-c(0.05,0.5)
doseparvar<-c(10,100,1000)
intparmean<-c(0,0.1)
intparvar<-c(10,100,1000)
cycleparmean<-c(0)
cycleparvar<-c(10,100,1000)


nTTP_prior_par_list<-list()
nTTP_prior_index_list<-list()
index<-1

for(p1 in 1: length(doseparmean)){
  for(p2 in 1: length(doseparvar)){
    for(p3 in 1: length(intparmean)){
      for(p4 in 1: length(intparvar)){
        for(p5 in 1: length(cycleparmean)){
          for(p6 in 1: length(cycleparvar)){
            nTTP_prior_index_list[[index]]<-c(p1,p2,p3,p4,p5,p6)
            nTTP_prior_par_list[[index]]<-c(doseparmean[p1],
                                            doseparvar[p2],
                                            intparmean[p3],
                                            intparvar[p4],
                                            cycleparmean[p5],
                                            cycleparvar[p6])
            index<-index+1
            
          }
        }
      }
    }
    
  }
}

#execute calibration

for(pc in 1:length(nTTP_prior_par_list)){
  
  for(sc in 1:length(prior_cal_scen_list)){
    
    
    cal_time <- system.time({
      assign(paste(c("nTTP_priorcal_s",sc,"_p",pc),collapse = ""), foreach(i=1:1000, combine = list) %dopar% {
        ##function
        nTTP_stage2(seed=i,tru_cyc1=prior_cal_scen_list[[sc]],co_size=co_size,ncohorts=ncohorts,target=target,doses=doses,nTTP_array=nTTP_array,ncycles=ncycles,prior_vec=nTTP_prior_par_list[[pc]],cycle_dec=cycle_dec,
                    hard.safety.rule = 0, safety.stopping.low.unsafe= F, safety.stopping.high.toosafe = F,precision.stopping = F)
        
      })
    })
    save.image("nTTP_priorcal_stage1.RData")
  }
}
