##prior calibration outline

#packages
library(doParallel)
registerDoParallel(cores=24)

#source code
source("POMM_function_stage2_comp-v7.R")


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
target<-c(0.3,0.391)
ncycles<-3
cycle_dec<-1/3


##define grid for parameter values
dosesmat<-matrix(c(0.05,.1,.15,.2,.25,.3,
                   0.05,.1,.2,.3,.5,.7,
                   0.15,.2,.25,.3,.35,.4,
                   0.15,.2,.3,.4,.5,.7),nrow=4,byrow=T)



#number of pseudo patients per dose
n0<-c(2,3,4,5,6)

grade_ratio_mat<-matrix(c(0,0.2,0.4,0.6,0.8,1,
                          0.05,0.1,0.15,0.2,0.25,0.3,
                          0.1,0.2,0.3,0.5,0.7,0.9,
                          0.2,0.3,0.4,0.5,0.6,0.7),nrow=4,byrow=T)



POMM_prior_par_list<-list()
POMM_prior_index_list<-list()
POMM_prior_data_list<-list()
index<-1

for(p1 in 1: nrow(dosesmat)){
  for(p2 in 1:length(n0)){
    for(p3 in 1:nrow(grade_ratio_mat)){
      
      POMM_prior_index_list[[index]]<-c(p1,p2,p3)
      POMM_prior_par_list[[index]]<-c(dosesmat[p1,],
                                      n0[p2],
                                      grade_ratio_mat[p3,])
      
      POMM_prior_data_list[[index]]<-dousseau_prior_func2(prior_doses=dosesmat[p1,], num0=n0[p2], graderatio_vec=grade_ratio_mat[p3,], ncycles=ncycles, cycle_dec=cycle_dec) 
      
      index<-index+1
      
      
    }
  }
}




#execute calibration


for(pc in 1:length(POMM_prior_par_list)){
  
  
  for(sc in 1:length(prior_cal_scen_list)){
    
    
    cal_time <- system.time({
      assign(paste(c("POMM_priorcal_s",sc,"_p",pc),collapse = ""), foreach(i=1:1000, combine = list) %dopar% {
        ##function
        POMM_stage2(seed=i,tru_cyc1=prior_cal_scen_list[[sc]],co_size=co_size,ncohorts=ncohorts,target=target,doses=doses,prior_pseudodata=POMM_prior_data_list[[pc]],
                    nTTP_array=nTTP_array,ncycles=ncycles,cycle_dec=cycle_dec,
                    hard.safety.rule = 0, safety.stopping.low.unsafe= F, safety.stopping.high.toosafe = F,precision.stopping = F)
        
      })
    })
    save.image("POMM_priorcal_stage1.RData")
  }
}
