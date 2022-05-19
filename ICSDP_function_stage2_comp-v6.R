#####################################
##    ICSDP for comparison  #########
#####################################
#STAGE 2: hard safety stopping included (with options)
#STAGE 2: safety stopping for low and high doses
#STAGE 2: no dose skipping (k-fold rise on experimented doses)
#v5: max admissable in output
#v6: precision stopping only once 9 patients in trial

library(maxLik)
source("data_generation_functions-V2.R")

##pre-amble functions
#likelihood functions
Lik_func<-function(parms,rmat,qmat,rprior,qprior,doses){
  
  log_lik<-matrix(nrow=length(doses),ncol=3)
  for(j in 1:(length(doses))){
    for(l in 1:3){
      log_lik[j,l]<-log(
        ((1-exp(-exp(parms[l]+parms[4]*log(doses[j]))))^(rmat[j,l]+rprior[j,l]))*
          ((exp(-exp(parms[l]+parms[4]*log(doses[j]))))^(qmat[j,l]+qprior[j,l]))
      )
    }
  }
  return(sum(log_lik))
}

Lik_func_cycles<-function(parms,rmat,qmat,rprior,qprior,doses){
  
  log_lik<-matrix(nrow=length(doses),ncol=(length(parms)-1))
  for(j in 1:(length(doses))){
    for(l in 1:(length(parms)-1)){
      log_lik[j,l]<-log(
        ((1-exp(-exp(parms[l]+parms[length(parms)]*log(doses[j]))))^(rmat[j,l]+rprior[j,l]))*
          ((exp(-exp(parms[l]+parms[length(parms)]*log(doses[j]))))^(qmat[j,l]+qprior[j,l]))
      )
    }
  }
  return(sum(log_lik))
}

pi_j_function<-function(dose,estimates){
  out_vec<-c()
  
  
  for(i in 1:(length(estimates)-1)){
    out_vec[i]<-1-exp(-exp(estimates[i]+estimates[length(estimates)]*log(dose)))
  }
  return(out_vec)
}

p_j_hat<-function(dose,estimates){
  pi_vec<-pi_j_function(dose,estimates)
  p_vec<-c()
  
  
  
  p_vec[1]<-pi_vec[1]
  for(i in 2:(length(estimates)-1)){
    p_vec[i]<-prod(1-pi_vec[1:(i-1)])*pi_vec[i]
    
  }
  return(sum(p_vec))
  
}

gain_func<-function(delta,pjhat){
  
  gval<-1/(delta-pjhat)^2
  
  return(which.max(gval))
}

Rfunc<-function(j,l,estimates,rpost,qpost){
  pi_j<-pi_j_function(doses[j],estimates) # pi vector of length 3 (l)
  return( (rpost[j,l]*((log(1-pi_j[l]))^2) + rpost[j,l]*pi_j[l]*log(1-pi_j[l])*(1-log(1-pi_j[l])) -
             (rpost[j,l]+qpost[j,l])*((pi_j[l])^2)*log(1-pi_j[l]))/((pi_j[l])^2))
}

sub_Rfunc<-function(ncycles,estimates,rpost,qpost){
  R_submat<-matrix(0,nrow=ndoses,ncol=ncycles)
  for (i in 1:ndoses){
    for (j in 1:ncycles){
      R_submat[i,j]<-Rfunc(i,j,estimates,rpost,qpost)
    }
  }
  return(colSums(R_submat))
}

dose_level2dose_val_frame<-function(doses,frame){
  
  values<-sapply(frame$dose_level, function(x) doses[x])
  out<-cbind(frame,values)
  names(out)<-c(names(frame),"doseval")
  return(out)
}


#prior pseudo data generation
#INPUT
#pi1=P(DLT) on lowest dose
#piJ= P(DLT) on highest dose
#num0=number of pseudo patients per end dose (only pseudo patients on lowest and highest dose)
#ncycles=number of cycles
#cycle_dec=decrease in P(DLT) per cycle
#ndoses= number of doses
#rorq= "r" or "q" depending on whether output required is "r" for DLT or "q" for non DLT

#OUTPUT
#matrix of pseudo patient information
ICSDP_prior_func<-function(pi1,piJ,num0,ncycles,cycle_dec,ndoses,rorq){
  
  dk_prior_r1<-c()
  dk_prior_rJ<-c()
  
  dk_prior_n1<-c()
  dk_prior_nJ<-c()
  
  dk_priors1<-pi1*(cycle_dec)^c(0:(ncycles-1))
  dk_priorsJ<-piJ*(cycle_dec)^c(0:(ncycles-1))
  
  dk_prior_n1[1]<-dk_prior_nJ[1]<-num0
  
  
  for(i in 1:ncycles){
    dk_prior_r1[i]<-dk_prior_n1[i]*dk_priors1[i]
    dk_prior_n1[i+1]<-dk_prior_n1[i]-dk_prior_r1[i]
    
    dk_prior_rJ[i]<-dk_prior_nJ[i]*dk_priorsJ[i]
    dk_prior_nJ[i+1]<-dk_prior_nJ[i]-dk_prior_rJ[i]
    
  }
  
  dk_prior_n1<-dk_prior_n1[c(1:ncycles)]
  dk_prior_nJ<-dk_prior_nJ[c(1:ncycles)]
  
  rprior<-matrix(c(dk_prior_r1,rep(rep(0,ncycles),ndoses-2),dk_prior_rJ),byrow=T,ncol=ncycles,nrow=ndoses)
  nprior<-matrix(c(dk_prior_n1,rep(rep(0,ncycles),ndoses-2),dk_prior_nJ),byrow=T,ncol=ncycles,nrow=ndoses)
  qprior<-nprior-rprior
  
  
  if(rorq=="r"){
    out<-rprior
  }else if(rorq=="q"){
    out<-qprior
  }
  #browser()
  return(out)
}



##INPUT:
#seed=seed for reproducability
#tru_cyc1=true underlying probability of DLT in cycle 1
#co_size=cohort size
#ncohorts=maximum number of cohorts in study
#target=target toxicity (all cycles)
#doses= real doses
#rprior=pseudo patient r matrix
#qprior=pseudo patient q matrix
#nTTP_array= 5x5x5 array of nTTP values for each of the 125 combinations of type and grade of toxicity
#ncycles=number of cycles
#cycle_dec=decrease in P(DLT) by cycle
#dose.skipping=enforce dose skipping rule? T=no dose skipping (default)
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
#MTD.est: estimate of MTD on cts scale
#max_admissible: max admissible dose at the end of the trial (has hard safety eliminated any?)


ICSDP_stage2<-function(seed,tru_cyc1,co_size,ncohorts ,target,doses,rprior,qprior,nTTP_array,ncycles,cycle_dec,
                       sufficient.information=T,hard.safety.rule=95,safety.stopping.low.unsafe=T,safety.stopping.high.toosafe=T,
                       kfold.skipping=T,kfold=2,precision.stopping=T){
  set.seed(seed)
  rands<-runif(co_size*ncohorts)
  stop_vec<-rep(0,7)
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
  
  out_mat<-matrix(rep(0,((ncohorts+ncycles))*4),nrow=(ncohorts+ncycles),ncol=4)
  dosevec<-c()
  patient_ID1<-1
  current.time<-0
  
  nextdose<-1
  stop<-0
  
  
  while(stop==0){
    #these need to be reset each time!
    rmat<-matrix(0,nrow=ndoses,ncol=ncycles)
    qmat<-matrix(0,nrow=ndoses,ncol=ncycles)
    if(current.time==0){
      #first cohort
      all.data<-multiple_patients_nTTP(patient_ID1,tru_cyc1,ncycles ,nextdose,entry_time=current.time,nTTP_array = nTTP_array,num_patients = co_size,rands=rands,cycle_dec=cycle_dec)
      current.time<-current.time+1
      all.data<-dose_level2dose_val_frame(doses=doses,frame=all.data)
      patient_ID1<-max(all.data$patient_ID)+1
      current.data<-all.data[all.data$time_of<=current.time,]
      dosevec[current.time]<-nextdose
    }else if(current.time<ncohorts){
      #subsequent cohorts
      all.data<-rbind(all.data,dose_level2dose_val_frame(doses=doses,frame=multiple_patients_nTTP(patient_ID1,tru_cyc1,ncycles ,nextdose,entry_time=current.time,nTTP_array = nTTP_array,num_patients = co_size,rands=rands,cycle_dec=cycle_dec)))
      current.time<-current.time+1
      patient_ID1<-max(all.data$patient_ID)+1
      current.data<-all.data[all.data$time_of<=current.time,]
      dosevec[current.time]<-nextdose
    }else{
      current.time<-current.time+1
      current.data<-all.data[all.data$time_of<=current.time,]
    }
    
    #rmat and qmat for ICSDP
    for(p in 1:(nrow(current.data))){
      rmat[current.data$dose_level[p],current.data$cycle_num[p]]<-rmat[current.data$dose_level[p],current.data$cycle_num[p]] +current.data$DLT[p]
      qmat[current.data$dose_level[p],current.data$cycle_num[p]]<-qmat[current.data$dose_level[p],current.data$cycle_num[p]] +1-current.data$DLT[p]
      
    }
    
    
    
    
    ##MODEL FITTING
    #choosing the next dose
    
    max_log_lik<- maxLik(Lik_func_cycles,start=c(rep(-10,ncycles),1.8),rmat=rmat,qmat=qmat,rprior=rprior,qprior=qprior,doses=doses)
    
    
    pj_hat<- sapply(doses,p_j_hat,max_log_lik$estimate) #calculaes the pjhat for all doses
    
    
    
    ##DOSE RECOMMENDATION
    nextdose<-gain_func(target,pj_hat)
    

    #k fold dose increase on EXPERIMENTED DOSES
    if(kfold.skipping==T){
      
      prevdoseval<-doses[max(dosevec)]
      nextdoseval<-doses[nextdose]
      
      if(nextdoseval>(kfold*prevdoseval)){
        maxdoseval<-kfold*prevdoseval
        nextdose<-max(which(doses<=maxdoseval))
      }
    }
    
    
    
    ##STOPPING RULES
    
    #accuracy
    logTD_est<-log(log(1-target)/(-sum(exp(max_log_lik$estimate[c(1:ncycles)]))))/max_log_lik$estimate[ncycles+1]
    
    ##NEW FISHER MAT FOR ALL CYCLES
    Rjl<-matrix(rep(0,ndoses*ncycles),nrow=ndoses,ncol=ncycles)
    for(j in 1:ndoses){
      for(l in 1:ncycles){
        Rjl[j,l]<-Rfunc(j,l,max_log_lik$estimate,rpost=rmat+rprior,qpost=qmat+qprior)
      }
    }
    
    Fisher_mat<-matrix(0,nrow=ncycles+1,ncol=ncycles+1)
    diags<-sub_Rfunc(ncycles,max_log_lik$estimate,rpost=rmat+rprior,qpost=qmat+qprior)
    diag(Fisher_mat)[c(1:ncycles)]<-diags
    
    for(cycnum in 1:ncycles){
      Fisher_mat[cycnum,(ncycles+1)]<-Fisher_mat[(ncycles+1),cycnum]<-sum(log(doses)*Rjl[,cycnum])
    }
    
    Fisher_mat[(ncycles+1),(ncycles+1)]<-sum(((log(doses))^2)*rowSums(Rjl))
    
    logTD_grad<--c(exp(max_log_lik$estimate[c(1:ncycles)])/(max_log_lik$estimate[ncycles+1]*sum(exp(max_log_lik$estimate[c(1:ncycles)]))), logTD_est/max_log_lik$estimate[ncycles+1])
    
    
    logTD_var<-t(logTD_grad)%*%solve(Fisher_mat,tol=1e-22)%*%logTD_grad
    
    CI<-exp(c(logTD_est-1.96*sqrt(logTD_var),logTD_est+1.96*sqrt(logTD_var)))
    CIratio<-CI[2]/CI[1]
    
    #stopping rule on CI
    # if((CIratio<accuracy)&(stop!=1)){
    #   # browser()
    #   stop<-1
    #   stoptime<-current.time
    #   
    # }
    # 
    
    
    #number of patietnts per dose level
    npats_doses<-co_size*tabulate(dosevec,nbins = ndoses)
    
    #safety stopping
    if((safety.stopping.low.unsafe==T)&(npats_doses[1]>0)){
      cyc1_pr_q<-qprior[1,1]
      cyc1_pr_r<-rprior[1,1]
      cyc1_obs_q<-qmat[1,1]
      cyc1_obs_r<-rmat[1,1]
      
      cyc1_0.3g<-pbeta(0.3,1+cyc1_pr_r+cyc1_obs_r,1+cyc1_pr_q+cyc1_obs_q,lower.tail = FALSE)
      if(cyc1_0.3g>0.8){
        stop<-6
        stop_vec[6]<-1
        nextdose<-NA
       # break
      }
    }
    
    if((safety.stopping.high.toosafe==T)&(npats_doses[ndoses]>0)){
      cycJ_pr_q<-qprior[ndoses,1]
      cycJ_pr_r<-rprior[ndoses,1]
      cycJ_obs_q<-qmat[ndoses,1]
      cycJ_obs_r<-rmat[ndoses,1]
      
      cycJ_0.3l<-pbeta(0.3,1+cycJ_pr_r+cycJ_obs_r,1+cycJ_pr_q+cycJ_obs_q,lower.tail = TRUE)
      if(cycJ_0.3l>0.8){
        stop<-7
        stop_vec[7]<-1
        nextdose<-NA
       # break
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
      #  break
      }
  
      
    }
    
    ##sufficient information 
    if((sufficient.information==T)&(!is.na(nextdose))){
      if(npats_doses[nextdose]>=9){
        stop<-4
        stop_vec[4]<-1
      #  break
      }
    }
    
    #max cohort
    if(max(current.data$patient_ID)==(ncohorts*co_size)){
      stop<-3
      stop_vec[3]<-1
    #  break
    }
    
    ##precision stopping
    if((precision.stopping==T)&(sum(npats_doses)>=9)){
      
      MTD_vector<-exp(rnorm(10000,mean=logTD_est,sd=sqrt(logTD_var)))
      if(any(is.infinite(MTD_vector))){
        CV<-1
      }else{
        CV<-mad(MTD_vector, center = median(MTD_vector), constant = 1.4826, 
                na.rm = FALSE, low = FALSE, high = FALSE)/median(MTD_vector)
      }
      if(CV<0.3){
        
        stop<-2
        stop_vec[2]<-1
        #break
      }
    }
    
  }
  
  
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
  
  
  #output list
  output<-list(
    dose.rec=nextdose ,#
    num.pat=max(current.data$patient) ,#: number of patients
    dose.ass=tabulate(dosevec,nbins = ndoses) ,# : number of cohorts per dose (vector)
    stop.code=stop_vec, #: stopping reason 
    num.DLT= sum(current.data$DLT[current.data$patient_ID>0]),#: number of DLTS (total)
    grade.mat=grade.matrix.out ,#: 5xnumdoses matrix. rows for grades 0,1,2,3,4. cols for doses.
    duration= current.time,#: total trial duration (until all recruited patients are fully observed)
    MTD.est=exp(logTD_est), # estimate of MTD on cts scale
    max_admissable=max_admissable # max admissable dose at the end of the trial (has hard safety eliminated any?)
  )
  
  
  return(output)
}

