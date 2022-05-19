
#####################################
####    POMM for comparison    ######
#####################################
#STAGE 2: hard safety stopping included (with options)
#STAGE 2: safety stopping for low and high doses
#STAGE 2: no dose skipping (k-fold rise on experimented doses)

#v5: if we have seen at least one DLT, use GLM (as opposed to needed BOTH DLT and non-DLT)
#v6: max admissable in output
#v7: precision stopping only once 9 patients in trial

library(ordinal)
library(mvtnorm)
source("data_generation_functions-V2.R")
##pre-amble functions

expit<-function (x) {
  return(exp(x)/(1 + exp(x)))
}

logit<-function(p){
  return(log(p/(1-p)))
}
#generate prior pseudo-data
#INPUT
#prior_doses=prob(DLT) per dose
#num0=number of pseudo patients per dose
#graderatio_vec=P(grade 2|no DLT)
#ncycles=number of cycles
#cycle_dec=decrease in P(DLT) by cycle

#OUTPUT
#dataframe of pseudo data
dousseau_prior_func2<-function (prior_doses, num0, graderatio_vec, ncycles, cycle_dec) 
{
  pseudo_data<-c()
  for(i in 1:length(prior_doses)){
    probs_cycles<-prior_doses[i] * (cycle_dec)^c(0:(ncycles - 1))
    pseudo_pat<-matrix(0, nrow=3,ncol=ncycles)
    #grade 1
    for(cycs in 1:ncycles){
      #grade 1 or 0
      pseudo_pat[1,cycs]<-prod((1-probs_cycles)[c(1:cycs)])*num0*(1-graderatio_vec[i])
      pseudo_data1<-c(-i,1,cycs,i,-3,i-4,0,NA,NA,pseudo_pat[1,cycs])
      #grade 2
      pseudo_pat[2,cycs]<-prod((1-probs_cycles)[c(1:cycs)])*num0*(graderatio_vec[i])
      pseudo_data2<-c(-i,2,cycs,i,-3,i-4,0,NA,NA,pseudo_pat[2,cycs])
      #grade 3
      pseudo_pat[3,cycs]<-prod(c(1,1-probs_cycles)[c(1:cycs)])*num0*probs_cycles[cycs]
      pseudo_data3<-c(-i,3,cycs,i,-3,i-4,1,NA,NA,pseudo_pat[3,cycs])
      pseudo_data<-rbind(pseudo_data,pseudo_data1,pseudo_data2,pseudo_data3)
    }
  }
  
  pseudo_data<-data.frame(pseudo_data)
  names(pseudo_data) <- c("patient", "Yobs", "cycle", "dose_level", 
                          "entry_time", "timeof", "DLT", "nTTP","max_grade", "weight")
  return(pseudo_data)
}


#probability of DLT in cycle 1 based on parameter estimation
P3<-function(parameter_est,dose){
  1-expit(parameter_est[2]-parameter_est[3]*dose)
  
}

#probability of DLT in all cycles based on vector of conditional probabilities of DLT for each cycle
cond_cyc_func<-function(cond_vec){
  length_c<-length(cond_vec)
  p3<-c()
  p3[1]<-cond_vec[1]
  for ( i in 2:length_c){
    p3[i]<-prod(1-cond_vec[1:(i-1)])*cond_vec[i]
  }
  return(sum(p3))
}


#addition of dose value to patient data frame
dose_level2dose_val_frame<-function(doses,frame){
  
  values<-sapply(frame$dose_level, function(x) doses[x])
  out<-cbind(frame,values)
  names(out)<-c(names(frame),"doseval")
  return(out)
}


##precision stopping

#preamble function
MTD_confunc<-function(dose_value,parameters,target){
  p3_vec<-1-expit(parameters[1]-parameters[2]*dose_value-parameters[3]*c(1:ncycles))
  return(cond_cyc_func(p3_vec)-target)
}
##INPUT:
#seed=seed for reproducability
#tru_cyc1=true underlying probability of DLT in cycle 1
#co_size=cohort size
#ncohorts=maximum number of cohorts in study
#target=target toxicity vector (cycle 1, all cycles)
#doses= real doses
#prior_pseudodata=prior pseudo dataframe
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
#stop.code: stopping reason (1=model fit failure,2=safety, all DLTs in first part,4=sufficient information,3= max patients reached)
#num.DLT: number of DLTS (total)
#grade.mat: 5xnumdoses matrix. rows for grades 0,1,2,3,4. cols for doses. max grade observed by each patient.
#duration: total trial duration (until all recruited patients are fully observed)
#max_admissible: max admissible dose at the end of the trial (has hard safety eliminated any?)


POMM_stage2<-function(seed,tru_cyc1,co_size,ncohorts ,target,doses,prior_pseudodata,nTTP_array,ncycles,cycle_dec,
                      sufficient.information=T,hard.safety.rule=95,safety.stopping.low.unsafe=T,safety.stopping.high.toosafe=T,
                      kfold.skipping=T,kfold=2,precision.stopping=T){
  set.seed(seed)
  rands<-runif(co_size*ncohorts)
  stop_vec<-rep(0,7)
  if(length(target)!=2){
    stop("target must be defined for first and all cycles")
  }
  
  ndoses<-length(doses)
  
  #define all doses as admissable before any are dropped for safety
  max_admissable<-ndoses
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
  
  #browser()
  #initialize
  dosevec<-c()
  patient_ID1<-1
  current.time<-0
  nextdose<-1
  recdose<-NA
  stop<-0
  initial<-1
  ass_doses<-c()
  
  while((current.time<=ncohorts)&(stop==0)){
    #DATA GENERATION
  #  browser()
    if(current.time==0){
      #first cohort
      all.data<-multiple_patients_nTTP(patient_ID1,tru_cyc1,ncycles ,nextdose,entry_time=current.time,nTTP_array = nTTP_array,num_patients = co_size,rands=rands,cycle_dec=cycle_dec)
      all.data<-cbind(all.data,1)
      names(all.data)<-c("patient" ,"Yobs","cycle" , "dose_level", "entry_time", "timeof",   
                         "DLT","nTTP","max_grade","weight")
      all.data<-rbind(prior_pseudodata,all.data)
      all.data<-dose_level2dose_val_frame(doses=doses,frame=all.data)
      current.time<-current.time+1
      patient_ID1<-max(all.data$patient)+1
      current.data<-all.data[all.data$timeof<=current.time,]
      dosevec[current.time]<-nextdose
    }else if(current.time<ncohorts){
      #subsequent cohorts
      new.data<-cbind(multiple_patients_nTTP(patient_ID1,tru_cyc1,ncycles ,nextdose,entry_time=current.time,nTTP_array = nTTP_array,num_patients = co_size,rands=rands,cycle_dec=cycle_dec),1)
      names(new.data)<-c("patient" ,"Yobs","cycle" , "dose_level", "entry_time", "timeof",   
                         "DLT","nTTP","max_grade","weight")
      new.data<-dose_level2dose_val_frame(doses=doses,frame=new.data)
      all.data<-rbind(all.data,new.data)
      current.time<-current.time+1
      patient_ID1<-max(all.data$patient)+1
      current.data<-all.data[all.data$timeof<=current.time,]
      dosevec[current.time]<-nextdose
    }else{
      #no new cohorts, but time moves on
      current.time<-current.time+1
      current.data<-all.data[all.data$timeof<=current.time,]
    }
    
    
    #assigned doses
    ass_doses[current.time]<-nextdose
    
    
    if((current.time*co_size)<=15){
      #we only look at first cycles for the first 15 patients
      
      currentdata_c1<-current.data[current.data$cycle==1,]
      Y1bin<-as.numeric(currentdata_c1$Yobs==3)
      currentdata_c1<-cbind(currentdata_c1,Y1bin)
      # if we have seen a DLT, use glm
      if((max(currentdata_c1$Yobs[currentdata_c1$patient>0])==3)){
        
        
        try(  model_CRM<-glm(currentdata_c1$Y1bin ~ currentdata_c1$doseval, family="binomial",weights=currentdata_c1$weight),silent = T)
        nextdose<-which.min(abs(target[1]- expit(model_CRM$coefficients[1] +model_CRM$coefficients[2]*doses)))
        
        
        #k fold dose increase on EXPERIMENTED DOSES
        if(kfold.skipping==T){
          
          prevdoseval<-doses[max(dosevec)]
          nextdoseval<-doses[nextdose]
          
          if(nextdoseval>(kfold*prevdoseval)){
            maxdoseval<-kfold*prevdoseval
            nextdose<-max(which(doses<=maxdoseval))
          }
        }
        
        #number of patients per dose (currently)
        npats_doses<-co_size*tabulate(dosevec,nbins = ndoses)
        
        #safety stopping
        if((safety.stopping.low.unsafe==T)&(npats_doses[1]>0)){
          #extract Hessian from model
          
          varmat<-vcov(model_CRM)
          #1st derivative of linear part wrt parameters
          dev_vec<-c(1,doses[1])
          #mean of linear part of model
          mu1<-sum(dev_vec*model_CRM$coefficients)
          #variance of linear part of model
          sigma1<-dev_vec%*%varmat%*%dev_vec
          
          #linear part assumed normality - pi is logit transformed
          cyc1_0.3g<-pnorm(log(0.3/(1-0.3)),mean=mu1,sd=sqrt(sigma1),lower.tail = F)
          # interest in upper tail of pi
          if(cyc1_0.3g>0.8){
            
            stop<-6
            stop_vec[6]<-1
            nextdose<-NA
        #    break
          }
        }
        
        if((safety.stopping.high.toosafe==T)&(npats_doses[ndoses]>0)){
          #extract Hessian from model
          varmat<-vcov(model_CRM)
          #1st derivative of linear part wrt parameters
          dev_vec<-c(1,doses[ndoses])
          #mean of linear part of model
          muJ<-sum(dev_vec*model_CRM$coefficients)
          #variance of linear part of model
          sigmaJ<-dev_vec%*%varmat%*%dev_vec
          # interest in lower tail of pi, 
          cycJ_0.3l<-pnorm(log(0.3/(1-0.3)),mean=muJ,sd=sqrt(sigmaJ),lower.tail=T)
          if(cycJ_0.3l>0.8){
            stop<-7
            stop_vec[7]<-1
            nextdose<-NA
        #    break
          }
        }
        
        #number of DLTs per dose level
        nDLTs_doses<-sapply(c(1:ndoses),function(x) sum(current.data$DLT[(current.data$dose_level==x)&(current.data$patient>0)]))
        
        
        ##hard safety
        if(hard.safety==T){
          
          explored<-c(1:ndoses)[npats_doses>0]
          for(do in explored){
           # browser()
            if(nDLTs_doses[do]>=hard.safety.mat[1, which(hard.safety.mat[2,]==npats_doses[do])]){
              
              max_admissable<-min(c(max_admissable,do-1))
            }
          }
          if(!is.na(nextdose)){
          if(nextdose>max_admissable){
            nextdose<-max_admissable
          }
          }
          if(max_admissable==0){
            stop<-5
            stop_vec[5]<-1
            nextdose<-NA
            recdose<-NA
       #     break
          }
          
          
        }
        
        
      # browser()
        ##sufficient information 
        if((sufficient.information==T)&(!is.na(nextdose))){
          if(npats_doses[nextdose]>=9){
            stop<-4
            stop_vec[4]<-1
            recdose<-nextdose
        #    break
          }
        }
        
        #glm precision stopping
        if((precision.stopping==T)&(sum(npats_doses)>=9)){
          varmat<-vcov(model_CRM)
          
          par_dist<-rmvnorm(10000,mean=model_CRM$coefficients,sigma=varmat)
          
          MTD_vector<-(logit(target[1])-par_dist[,1])/par_dist[,2]
          
          CV<-mad(MTD_vector, center = median(MTD_vector), constant = 1.4826, 
                  na.rm = FALSE, low = FALSE, high = FALSE)/median(MTD_vector)
          
          if(CV<0.3){
            
            stop<-2
            stop_vec[2]<-1
            recdose<-nextdose
        #    break
          }
        }
        
        
      } else { 
      #  browser()
        
        if(sum(currentdata_c1$Y1bin[currentdata_c1$patient>0])==0){ #if we have seen no DLTs in cycle 1, next dose up
          nextdose<-min(c(nextdose+1,length(doses)))
        #  browser()
        }else{ #otherwise all DLTs stop for safety
          stop<-5
          stop_vec[5]<-1
          stop_vec[1]<-1
          nextdose<-NA
          recdose<-NA
        }
        
        
      }
      
    }  else if(current.time<ncohorts){ #applying POMM from patient 16 onwards
      #model 1
      
      try(model1_t<-clmm(as.factor(Yobs) ~ doseval + cycle + (1|patient),data=current.data,weights=current.data$weight,Hess = T),silent=T)
      
      if(exists("model1_t")){
        
        parms<-c(model1_t$coefficients[1:2],model1_t$coefficients[3:4],(model1_t$ST$patient)^2)
        if(any(is.na(parms))){ #if model1 has warnings, stop
          stop<-1
          stop_vec[1]<-1
          recdose<-NA
          break
        }
        
        #if all random effects are 0, refit the model without the random effects
        #(otherwise we get into problems with the variance-covariance matrix)
        if(all(model1_t$ranef==0)){
          #browser()
          try(model1_t<-clm(as.factor(Yobs) ~ doseval + cycle ,data=current.data,weights=current.data$weight),silent=T)
          parms<-c(model1_t$coefficients[1:2],model1_t$coefficients[3:4],0)
          if(any(is.na(parms))){ #if model1 has warnings, stop
            stop<-1
            stop_vec[1]<-1
            recdose<-NA
            break
          }
          
          
        }else if(model1_t$Hessian[5,5]==0){
          #browser()
          try(model1_t<-clm(as.factor(Yobs) ~ doseval + cycle ,data=current.data,weights=current.data$weight),silent=T)
          parms<-c(model1_t$coefficients[1:2],model1_t$coefficients[3:4],0)
          if(any(is.na(parms))){ #if model1 has warnings, stop
            stop<-1
            stop_vec[1]<-1
            recdose<-NA
            break
          }
          
        }else if(any(diag(solve(model1_t$Hessian))<=0)){
          try(model1_t<-clm(as.factor(Yobs) ~ doseval + cycle ,data=current.data,weights=current.data$weight),silent=T)
          parms<-c(model1_t$coefficients[1:2],model1_t$coefficients[3:4],0)
          if(any(is.na(parms))){ #if model1 has warnings, stop
            stop<-1
            stop_vec[1]<-1
            recdose<-NA
            break
          }
        }
        #calculating the P(DLT) for each dose
        p3_doses<-c()
        for (do in 1:ndoses){
          p3_vec<-1-expit(parms[2]-parms[3]*doses[do]-parms[4]*c(1:ncycles))
          p3_doses[do]<- cond_cyc_func(p3_vec)
        }
        
        #choose dose with P(DLT) closest to target
        nextdose<-which.min(abs(target[2]- p3_doses))
        
        
        #k fold dose increase on EXPERIMENTED DOSES
        if(kfold.skipping==T){
          
          prevdoseval<-doses[max(dosevec)]
          nextdoseval<-doses[nextdose]
          
          if(nextdoseval>(kfold*prevdoseval)){
            maxdoseval<-kfold*prevdoseval
            nextdose<-max(which(doses<=maxdoseval))
          }
        }
        
        #number of patietnts per dose level
        npats_doses<-co_size*tabulate(dosevec,nbins = ndoses)
        
        #safety stopping
        if((safety.stopping.low.unsafe==T)&(npats_doses[1]>0)){
          # browser()
          #adding the extra column and row of 0 ensures that if the clm was fitted, the variance matrix is of the correct dims
          # (if clmm was used, adding the 0s have no effect anyway)
          # browser()
          varmat<-cbind(rbind(solve(model1_t$Hessian),0),0)[c(2:5),c(2:5)]
          
          #1st derivative of linear part wrt parameters
          dev_vec<-c(1,-doses[1],-1,1)
          #mean of linear part of model
          mu1<-sum(dev_vec*c(parms[c(2:4)],0))
          #variance of linear part of model
          sigma1<-dev_vec%*%varmat%*%dev_vec

          #linear part assumed normality - pi is logit transformed
          cyc1_0.3g<-pnorm(log(0.7/(1-0.7)),mean=mu1,sd=sqrt(sigma1))
          # interest in upper tail of 1-pi(t), equivalent to taking lower tail of pi(1-t)
          if(cyc1_0.3g>0.8){
            stop<-6
            stop_vec[6]<-1
            nextdose<-NA
           # break
          }
        }
        
        if((safety.stopping.high.toosafe==T)&(npats_doses[ndoses]>0)){
          #adding the extra column and row of 0 ensures that if the clm was fitted, the variance matrix is of the correct dims
          # (if clmm was used, adding the 0s have no effect anyway)
          varmat<-cbind(rbind(solve(model1_t$Hessian),0),0)[c(2:5),c(2:5)]
          #1st derivative of linear part wrt parameters
          dev_vec<-c(1,-doses[ndoses],-1,1)
          #mean of linear part of model
          muJ<-sum(dev_vec*c(parms[c(2:4)],0))
          #variance of linear part of model
          sigmaJ<-dev_vec%*%varmat%*%dev_vec
          # interest in lower tail of 1-pi(t), equivalent to taking upper tail of pi(1-t))
          cycJ_0.3l<-pnorm(log(0.7/(1-0.7)),mean=muJ,sd=sqrt(sigmaJ),lower.tail=FALSE)
          if(cycJ_0.3l>0.8){
            stop<-7
            stop_vec[7]<-1
            nextdose<-NA
         #   break
          }
        }
        
        #number of DLTs per dose level
        nDLTs_doses<-sapply(c(1:ndoses),function(x) sum(current.data$DLT[(current.data$dose_level==x)&(current.data$patient>0)]))
        nDLTs_doses_c1<-sapply(c(1:ndoses),function(x) sum(current.data$DLT[(current.data$dose_level==x)&(current.data$patient>0)&(current.data$cycle==1)]))
        
        
        ##hard safety
        if(hard.safety==T){
          
          explored<-c(1:ndoses)[npats_doses>0]
          for(do in explored){
            
            if(nDLTs_doses_c1[do]>=hard.safety.mat[1, which(hard.safety.mat[2,]==npats_doses[do])]){
              
              max_admissable<-min(c(max_admissable,do-1))
            }
          }
          if(!is.na(nextdose)){
          if(nextdose>max_admissable){
            nextdose<-max_admissable
          }
          }
          if(max_admissable==0){
            stop<-5
            stop_vec[5]<-1
            recdose<-NA
            nextdose<-NA
         #   break
          }
          
          
        }
        
        
        ##sufficient information 
        if((sufficient.information==T)&(!is.na(nextdose))){
          if(npats_doses[nextdose]>=9){
            stop<-4
            stop_vec[4]<-1
            recdose<-nextdose
        #    break
          }
        }
        
        #pomm precision stopping
        if((precision.stopping==T)&(sum(npats_doses)>=9)){
          #adding the extra column and row of 0 ensures that if the clm was fitted, the variance matrix is of the correct dims
          # (if clmm was used, adding the 0s have no effect anyway)
          varmat<-solve(model1_t$Hessian)[c(2:4),c(2:4)]
          
          par_dist<-rmvnorm(10000,mean=parms[c(2:4)],sigma=varmat)
          #browser()
          MTD_vector<-unlist(apply(par_dist,1, function(x) uniroot(MTD_confunc,c(-100,100)/abs(x[2]),x,target[2])$root))
          CV<-mad(MTD_vector, center = median(MTD_vector), constant = 1.4826, 
                  na.rm = FALSE, low = FALSE, high = FALSE)/median(MTD_vector)
          
          if(CV<0.3){
            
            stop<-2
            stop_vec[2]<-1
            recdose<-nextdose
          #  break
          }
        }
        
      }else{ # stop for model fitting failure
        stop<-1
        recdose<-NA 
        break
      }
      
    }else { #final cohort
      current.data<-all.data
      
      try( model1_t<-clmm(as.factor(Yobs) ~ doseval + cycle + (1|patient),data=current.data,weights=current.data$weight),silent=T)
      
      if(exists("model1_t")){
        
        parms<-c(model1_t$coefficients[1:2],model1_t$coefficients[3:4],(model1_t$ST$patient)^2)
        if(any(is.na(parms))){ #if model1 has warnings
          
          stop<-1
          stop_vec[1]<-1
          recdose<-NA
          break
        }
        
        
        
        #calculating the P(DLT) for each dose
        p3_doses<-c()
        for (do in doses){
          p3_vec<-1-expit(parms[2]-parms[3]*do-parms[4]*c(1:ncycles))
          p3_doses[doses[do]]<- cond_cyc_func(p3_vec)
        }
        
        #choose dose with P(DLT) closest to target
        recdose<-which.min(abs(target[2]- p3_doses))
        
       
        
        
        
      }else{ #stop for model fitting failure
        stop<-1
        stop_vec[1]<-1
        recdose<-NA
        break
      }
      
      
      
    }
    
    
    #remove model stored
    if(exists("model1_t")){
      rm(model1_t)
    }
    
  } #ends while (current.time<=ncohorts)
  
  
  
  
  
  #follow up for all patients
  current.data<-all.data
  current.time<-max(current.data$timeof)
  
  
  #grade matrix out
  grade.matrix.out<-matrix(0,ncol=ndoses,nrow=5)
  for (pat in 1:(max(current.data$patient))){
    pat_max_grade<-max(current.data[current.data$patient==pat,]$max_grade)
    pat_dose<-current.data[current.data$patient==pat,]$dose_level[1]
    grade.matrix.out[pat_max_grade+1,pat_dose]<- grade.matrix.out[pat_max_grade+1,pat_dose]+1
  }
  
  #to allign stop code with other methods
  if(stop==0){
    stop<-3
    stop_vec[3]<-1
  }
  #output list
  output<-list(
    dose.rec=recdose ,#
    num.pat=max(current.data$patient) ,#: number of patients
    dose.ass=tabulate(dosevec,nbins = ndoses) ,# : number of cohorts per dose (vector)
    stop.code=stop_vec, #: stopping reason 
    num.DLT= sum(current.data$DLT[current.data$patient>0]),#: number of DLTS (total)
    grade.mat=grade.matrix.out ,#: 5xnumdoses matrix. rows for grades 0,1,2,3,4. cols for doses.
    duration= current.time,#: total trial duration (until all recruited patients are fully observed)
    max_admissable=max_admissable # max admissable dose at the end of the trial (has hard safety eliminated any?)
  
    )
  

  return(output)
  
}





