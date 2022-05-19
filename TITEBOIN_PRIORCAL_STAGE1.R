##prior calibration outline
#TITEBOIN 
#STAGE 1

#packages
library(doParallel)
registerDoParallel(cores=24)

#source code
source("TITEBOIN_function_stage2_comp-v5.R")


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
target<-0.391
ncycles<-3
cycle_dec<-1/3
ndoses<-length(doses)


##define grid for parameter values
#a_vec<-seq(0.6,0.9,0.1)
#b_vec<-seq(1.1,1.4,0.1)
##define grid for parameter values
a_vec<-c(0.8,0.9)
b_vec<-c(1.3,1.4)

prior_ss<-c(1,2)
prior_me<-c(0.1,0.5*0.391,0.391,0.5)



TITEBOIN_prior_par_list<-list()

index<-1
for(p1 in 1: length(a_vec)){
  for(p2 in 1:length(b_vec)){
    for(p3 in 1:length(prior_ss)){
      for(p4 in 1:length(prior_me)){
        TITEBOIN_prior_par_list[[index]]<-c(a_vec[p1],
                                            b_vec[p2],prior_ss[p3]*c(prior_me[p4],1-prior_me[p4]))
        
        index<-index+1
        
        
      }
    }
    
  }
}




#execute calibration


for(pc in 1:length(TITEBOIN_prior_par_list)){
  
  
  for(sc in 1:length(prior_cal_scen_list)){
    
    
    cal_time <- system.time({
      assign(paste(c("TITEBOIN_priorcal_s",sc,"_p",pc),collapse = ""), foreach(i=1:1000, combine = list) %dopar% {
        ##function
        TITEBOIN_stage2(seed=i,tru_cyc1=prior_cal_scen_list[[sc]],co_size=co_size,ncohorts=ncohorts,
                        a=TITEBOIN_prior_par_list[[pc]][1],  b=TITEBOIN_prior_par_list[[pc]][2], alpha=TITEBOIN_prior_par_list[[pc]][3], beta=TITEBOIN_prior_par_list[[pc]][4],
                        # a=TITEBOIN_prior_par_list[[pc]][1],  b=TITEBOIN_prior_par_list[[pc]][2], alpha=0.5*target, beta=1-0.5*target,
                        target=target,doses=doses,nTTP_array=nTTP_array,ncycles=ncycles,cycle_dec=cycle_dec,
                        hard.safety.rule = 0, safety.stopping.low.unsafe= F, safety.stopping.high.toosafe = F)
        
      })
    })
   # save.image("TITEBOIN_priorcal_stage1.RData")
    save.image("TITEBOIN_priorcal_stage1add.RData")
  }
}
