#####################################
##    RMD for comparison  ######
#####################################
#STAGE 2: hard safety stopping included (with options)
#STAGE 2: safety stopping for low and high doses
#STAGE 2: no dose skipping (k-fold rise on experimented doses)
#v8: max admissible in output
#v9: precision stopping only once 9 patients in trial
library(rjags)


#full nTTP data gen
source("data_generation_functions-V2.R")

##pre-amble functions
#number of cycle observations per patient
n4Gibbs<-function(vec_in){
  if(sum(!is.na(vec_in))!=0){
    if(max(vec_in[!is.na(vec_in)])==1){
      return( min(which.max(vec_in[!is.na(vec_in)])))
    }else{
      return(sum(!is.na(vec_in)))
    }
  }else{
    return(sum(!is.na(vec_in)))
  }
  
}

#gibbs sampler function
gibbs_sampler<-function(data_list,iter,model_string){
  model1.spec<-textConnection(model_string)
  jags <- jags.model(model1.spec,data =data_list,n.chains=1,n.adapt=iter,quiet=T)
  update(jags, iter,progress.bar="none")
  tt<-jags.samples(jags,c('beta1','beta02','gamma'),iter,progress.bar="none")
  return(tt)
}

loss_func<-function(beta,doses){
  beta[1]+beta[2]*doses +beta[3]
}



##INPUT:
#seed=seed for reproducability
#tru_cyc1=true underlying probability of DLT in cycle 1
#co_size=cohort size
#ncohorts=maximum number of cohorts in study
#target=target toxicity (all cycles)
#doses= real doses
#prior_vec= vector of hyperparameters for prior (dose_par:mean, variance, intercept:mean, variance, cycle_par:mean,variance)
#nTTP_array= 5x5x5 array of nTTP values for each of the 125 combinations of type and grade of toxicity
#ncycles=number of cycles
#cycle_dec=decrease in P(DLT) by cycle
#sufficient.information==enforce stopping for sufficient information? no more than 9 patients per dose. T=enforce stopping for sufficient information (default)
#hard.safety.rule=enforce hard safety rule based on Beta(1.1)? 
#options: 85= 2/3,3/6,4/9. 90=2/3,4/6,5/9, 95=3/3,4/6,5/9, 0=no hard safety enforced
#safety.stopping.low.unsafe= (T=stop when P(p1>0.3)>0.8 (cycle 1))
#safety.stopping.high.toosafe= (T= stop when P(pJ>0.3)>0.8 (cycle 1))
#kfold.skipping: Is a kfold dose skipping rule implemented?
#kfold: the k-fold rise thsat is allowable in dose skipping
#precision.stopping: Is a precision stopping rule implemented?



##OUTPUT:
#list with elements corresponding to:
#dose.rec : dose recommendation
#num.pat: number of patients
#dose.ass : number of cohorts per dose (vector)
#stop.code: stopping reason (4=sufficient information,3= max patients reached)
#num.DLT: number of DLTS (total)
#grade.mat: 5xnumdoses matrix. rows for grades 0,1,2,3,4. cols for doses. max grade observed by each patient.
#duration: total trial duration (until all recruited patients are fully observed)
#max_admissible: max admissible dose at the end of the trial (has hard safety eliminated any?)

nTTP_stage2<-function(seed,tru_cyc1,co_size,ncohorts,target,doses,nTTP_array,ncycles,prior_vec,cycle_dec,
                      sufficient.information=T,hard.safety.rule=95,safety.stopping.low.unsafe=T,safety.stopping.high.toosafe=T,
                      kfold.skipping=T,kfold=2,precision.stopping=T){
  set.seed(seed)
  rands<-runif(co_size*ncohorts)
  stop_vec<-rep(0,7)
  try({
    
    ynTTP<-matrix(NA,nrow=ncohorts*co_size,ncol=ncycles)
    yDLT<-matrix(NA,nrow=ncohorts*co_size,ncol=ncycles)
    ydose<-matrix(NA,nrow=ncohorts*co_size,ncol=ncycles)
    ycycle<-matrix(NA,nrow=ncohorts*co_size,ncol=ncycles)
    
    model1.string <-"
model {
for (i in 1:m){
gamma[i] ~ dnorm(0,taugamma)
for(j in 1:n[i]){
eta[i,j] <- beta02[1] +beta1*x[i,j] +beta02[2]*tee[i,j] +gamma[i] 
ynTTP[i,j] ~ dnorm(eta[i,j],taueps)

}
}
beta1 ~ dnorm(u, 1/nu)
beta02[1:2] ~ dmnorm(prioralphabeta,priorBbeta)
taugamma ~ dgamma(agamma,bgamma)
taueps ~ dgamma(aeps,beps)


}
"
stop<-0
nextdose<-c()
nextdose[1]<-1


ndoses<-length(doses)

#define all doses as admissible before any are dropped for safety
max_admissible<-ndoses
if(hard.safety.rule==85){
  hard.safety<-T
  hard.safety.mat<-matrix(c(2,3,3,6,4,9),nrow=2)
}else if(hard.safety.rule==90){
  hard.safety<-T
  hard.safety.mat<-matrix(c(2,3,4,6,5,9),nrow=2)
}else if(hard.safety.rule==95){
  hard.safety<-T
  hard.safety.mat<-matrix(c(3,3,4,6,5,9),nrow=2)
}else{
  hard.safety<-F
}

dosevec<-c()
patient_ID1<-1
current.time<-0


initial<-1

while(stop==0){
  # time<-time+1
  
  #DATA GENERATION
  
  if(current.time==0){
    #first cohort
    all.data<-multiple_patients_nTTP(patient_ID1,tru_cyc1,ncycles ,nextdose[current.time+1],entry_time=current.time,nTTP_array = nTTP_array,num_patients = co_size,rands=rands,cycle_dec=cycle_dec)
    current.time<-current.time+1
    patient_ID1<-max(all.data$patient_ID)+1
    current.data<-all.data[all.data$time_of<=current.time,]
    dosevec[current.time]<-nextdose[1]
    ynTTP[c(1:co_size),1]<-current.data$nTTP
    yDLT[c(1:co_size),1]<-current.data$DLT
    ydose[c(1:3),1]<-1
    ycycle[c(1:3),1]<-1
    
    #initial phase before model fit
    if(sum(current.data$DLT)==0){
      nextdose[2]<-nextdose[1]+1
    } else{
      nextdose[2]<-nextdose[1]
    }
    
  }else{
    
    
    if(current.time<ncohorts){
      #subsequent cohorts
      
      
      all.data<-rbind(all.data,multiple_patients_nTTP(patient_ID1,tru_cyc1,ncycles ,nextdose[current.time+1],entry_time=current.time,nTTP_array = nTTP_array,num_patients = co_size,rands=rands,cycle_dec=cycle_dec))
      current.time<-current.time+1
      patient_ID1<-max(all.data$patient_ID)+1
      current.data<-all.data[all.data$time_of<=current.time,]
      dosevec[current.time]<-nextdose[current.time]
      
      
      #format for rjags
      for(obser in 1:nrow(current.data)){
        ynTTP[current.data$patient_ID[obser],current.data$cycle_num[obser]]<-current.data$nTTP[obser]
        yDLT[current.data$patient_ID[obser],current.data$cycle_num[obser]]<-current.data$DLT[obser]
        ydose[current.data$patient_ID[obser],current.data$cycle_num[obser]]<-doses[current.data$dose_level[obser]]
        ycycle[current.data$patient_ID[obser],current.data$cycle_num[obser]]<-current.data$cycle_num[obser]
        
      }
      
      
      
      
      
      n<-apply(yDLT,1,n4Gibbs)
      n<-n[n>0]
      
      currentdata<-list(ynTTP=ynTTP,tee=ycycle,x=ydose,m=length(n),
                        n=n,
                        u=prior_vec[1],
                        nu=prior_vec[2],
                        prioralphabeta=c(prior_vec[3],prior_vec[5]),
                        priorBbeta=c(1/prior_vec[4],1/prior_vec[6])*diag(2),
                        agamma=0.001,
                        bgamma=0.001,
                        aeps=0.001,
                        beps=0.001)
      
      gibbs_out<-gibbs_sampler(data_list=currentdata,iter=10000,model_string=model1.string)
      
      
      #as described in paper, only include output that has beta1>0. If no ouput satisfies, run again (extremely unlikely!)
      while(length(gibbs_out$beta1[gibbs_out$beta1>0])==0){
        gibbs_out<-gibbs_sampler(data_list=currentdata,iter=10000,model_string=model1.string)
      }
      
      
      
      nextdose[current.time+1]<-which.min(abs(target-loss_func(beta=c(rowMeans(gibbs_out$beta02[,gibbs_out$beta1>0,])[1],
                                                                      mean(gibbs_out$beta1[gibbs_out$beta1>0]),
                                                                      rowMeans(gibbs_out$beta02[,gibbs_out$beta1>0,])[2]),
                                                               doses)))
      #k fold dose increase on EXPERIMENTED DOSES
      if(kfold.skipping==T){
        
        prevdoseval<-doses[max(dosevec)]
        nextdoseval<-doses[nextdose[current.time+1]]
        
        if(nextdoseval>(kfold*prevdoseval)){
          maxdoseval<-kfold*prevdoseval
          nextdose[current.time+1]<-max(which(doses<=maxdoseval))
        }
      }
      
      
      #number of patietnts per dose level
      npats_doses<-co_size*tabulate(dosevec,nbins = ndoses)
      
      #safety stopping
      if((safety.stopping.low.unsafe==T)&(npats_doses[1]>0)){
        
        posterior1<- gibbs_out$beta02[1,gibbs_out$beta1>0,] + doses[1]*gibbs_out$beta1[gibbs_out$beta1>0] + gibbs_out$beta02[2,gibbs_out$beta1>0,]
        
        cyc1_0.3g<-  mean(posterior1>target)
        if(cyc1_0.3g>0.8){
          stop<-6
          stop_vec[6]<-1
          nextdose[current.time+1]<-NA
          #  break
        }
      }
      
      if((safety.stopping.high.toosafe==T)&(npats_doses[ndoses]>0)){
        
        posteriorJ<-gibbs_out$beta02[1,gibbs_out$beta1>0,] + doses[ndoses]*gibbs_out$beta1[gibbs_out$beta1>0] + gibbs_out$beta02[2,gibbs_out$beta1>0,]
        
        
        cycJ_0.3l<-mean(posteriorJ<target)
        if(cycJ_0.3l>0.8){
          stop<-7
          stop_vec[7]<-1
          nextdose[current.time+1]<-NA
          #  break
        }
      }
      
      #number of DLTs per dose level
      nDLTs_doses<-sapply(c(1:ndoses),function(x) sum(current.data$DLT[(current.data$dose_level==x)&(current.data$patient_ID>0)]))
      nDLTs_doses_c1<-sapply(c(1:ndoses),function(x) sum(current.data$DLT[(current.data$dose_level==x)&(current.data$patient>0)&(current.data$cycle==1)]))
      
      ##hard safety
      if(hard.safety==T){
        
        explored<-c(1:ndoses)[npats_doses>0]
        for(do in explored){
          
          if(nDLTs_doses_c1[do]>=hard.safety.mat[1, which(hard.safety.mat[2,]==npats_doses[do])]){
            
            max_admissible<-min(c(max_admissible,do-1))
          }
        }
        if(!is.na(nextdose[current.time+1])){
          if(nextdose[current.time+1]>max_admissible){
            nextdose[current.time+1]<-max_admissible
          }
        }
        if(max_admissible==0){
          stop<-5
          stop_vec[5]<-1
          #   break
        }
        
        
      }
      
      
      
      
      ##sufficient information 
      if((sufficient.information==T)&(!is.na(nextdose[current.time+1]))){
        if(npats_doses[nextdose[current.time+1]]>=9){
          stop<-4
          stop_vec[4]<-1
          #   break
        }
      }
      
      
      ##precision stopping
      
      if((precision.stopping==T)&(sum(npats_doses)>=9)){
        MTD_vector<-(target - gibbs_out$beta02[1,gibbs_out$beta1>0,]-gibbs_out$beta02[2,gibbs_out$beta1>0,])/gibbs_out$beta1[gibbs_out$beta1>0]
        CV<-mad(MTD_vector, center = median(MTD_vector), constant = 1.4826, 
                na.rm = FALSE, low = FALSE, high = FALSE)/median(MTD_vector)
        # browser()
        if(CV<0.3){
          stop_vec[2]<-1
          stop<-2
          
        }
      }
      
      
    }else{
      #max patients
      stop<-3
      stop_vec[3]<-1
      #   break
    }
    
  }
  
  
}
  
  
  
  
  })
  
  #follow up for all patients
  current.data<-all.data
  current.time<-max(current.data$time_of)
  
  
  #grade matrix out
  grade.matrix.out<-matrix(0,ncol=ndoses,nrow=5)
  for (pat in 1:(max(current.data$patient))){
    pat_max_grade<-max(current.data[current.data$patient_ID==pat,]$max_grade)
    pat_dose<-current.data[current.data$patient_ID==pat,]$dose_level[1]
    grade.matrix.out[pat_max_grade+1,pat_dose]<- grade.matrix.out[pat_max_grade+1,pat_dose]+1
  }
  
  if(any(stop_vec[c(5:7)]>=1)){
    dose.rec<-NA
  }else{
    dose.rec<-nextdose[length(nextdose)]
  }
  
  # browser()
  #output list
  output<-list(
    dose.rec=dose.rec,#
    num.pat=max(current.data$patient) ,#: number of patients
    dose.ass=tabulate(dosevec,nbins = ndoses) ,# : number of cohorts per dose (vector)
    stop.code=stop_vec, #: stopping reason
    num.DLT= sum(current.data$DLT[current.data$patient_ID>0]),#: number of DLTS (total)
    grade.mat=grade.matrix.out ,#: 5xnumdoses matrix. rows for grades 0,1,2,3,4. cols for doses.
    duration= current.time,#: total trial duration (until all recruited patients are fully observed)
    max_admissible=max_admissible # max admissible dose at the end of the trial (has hard safety eliminated any?)
    
  )
  
  
  return(output)
  
  
}





