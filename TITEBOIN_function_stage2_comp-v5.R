#####################################
##    TITE-BOIN for comparison  ######
#####################################
library(Iso)
#STAGE 2: hard safety stopping included
#STAGE 2: safety stopping for low and high doses
#STAGE 2: no dose skipping k-fold
#v5: max admissable in output


source("data_generation_functions-V2.R")
##pre-amble functions
current_patient_data_frame<-function(current_time,patient.dataframe,follow_up){
  #input the whole data and convert to format for TITECRM
  #each row is a patient, with entries c("entry.time","dose.level","DLT","DLT.time") 
  num_patients<-max(patient.dataframe$patient_ID)
  
  patient_dataframe<-matrix(NA,nrow=num_patients,ncol=4)
  for (i in 1:num_patients){
    patient_data_ind<-patient.dataframe[patient.dataframe$patient_ID==i,]
    patient_DLT<-max(patient_data_ind$DLT)
    patient_entry.time<-patient_data_ind$entry_time[1]
    patient_dose.level<-patient_data_ind$dose_level[1]
    if(patient_DLT==1){
      patient_DLT.time<-min(which(patient_data_ind$DLT==1))+patient_entry.time
    }else{
      patient_DLT.time<-NA
    }
    patient_dataframe[i,]<-c(patient_entry.time,patient_dose.level,patient_DLT,patient_DLT.time)
  }
  
  patient_dataframe<-data.frame(patient_dataframe)
  names(patient_dataframe)<-c("entry.time","dose.level","DLT","DLT.time") 
  follow_up_time<-function(x) min(follow_up,x)
  patient_follow_up<- sapply(current_time-patient_dataframe$entry.time,follow_up_time)
  patient_weights<-patient_follow_up/follow_up
  
  #browser()
  #only DLT if we have seen it at the current time
  current.DLT<-patient_dataframe$DLT
  current.DLT[which(current.DLT==1)]<-patient_dataframe$DLT.time[which(current.DLT==1)]<=current_time
  
  #only DLT time observed id current DLT is true
  current.DLT.time<-patient_dataframe$DLT.time
  current.DLT.time[which(current.DLT==0)]<-NA
  
  patient_weights[which(current.DLT==1)]<-1
  entry.time<-patient_dataframe$entry.time
  dose.level<-patient_dataframe$dose.level
  return(data.frame(entry.time,dose.level,current.DLT,current.DLT.time,patient_weights,patient_follow_up))
  
}


patient_data_frame<-function(patient.dataframe){
  #input the whole data and convert to format for TITEBOIN
  #each row is a patient, with entries c("entry.time","dose.level","DLT","DLT.time") 
  num_patients<-max(patient.dataframe$patient_ID)
  #browser()
  patient_dataframe<-matrix(NA,nrow=num_patients,ncol=4)
  for (i in 1:num_patients){
    patient_data_ind<-patient.dataframe[patient.dataframe$patient_ID==i,]
    patient_DLT<-max(patient_data_ind$DLT)
    patient_entry.time<-patient_data_ind$entry_time[1]
    patient_dose.level<-patient_data_ind$dose_level[1]
    if(patient_DLT==1){
      patient_DLT.time<-min(which(patient_data_ind$DLT==1))+patient_entry.time
    }else{
      patient_DLT.time<-NA
    }
    patient_dataframe[i,]<-c(patient_entry.time,patient_dose.level,patient_DLT,patient_DLT.time)
  }
  
  patient_dataframe<-data.frame(patient_dataframe)
  names(patient_dataframe)<-c("entry.time","dose.level","DLT","DLT.time") 
  return(patient_dataframe)
  
}


TITE_BOIN_boundary<-function(c,n,s,phi,alpha,beta,a,b){
  
  phi1<-a*phi
  phi2<-b*phi
  
  r<-n-c
  p_tilda<-(s+alpha)/(r+alpha+beta)
  
  lambda_e<-log((1-phi1)/(1-phi))/log((phi*(1-phi1))/(phi1*(1-phi)))
  lambda_d<-log((1-phi)/(1-phi2))/log((phi2*(1-phi))/(phi*(1-phi2)))
  
  if((s/n)>phi){
    pi_e<-Inf
  }else{
    pi_e<-(c-((1-p_tilda)/p_tilda)*(n*lambda_e-s))*as.numeric((s/n)<=phi)
  }
  
  
  if((s/n)<phi){
    pi_d<--Inf
  }else{
    pi_d<-(c-((1-p_tilda)/p_tilda)*(n*lambda_d-s))*as.numeric((s/n)>=phi)
    
  }
  # browser()
  return(c(pi_d,pi_e))
}


#function to figure out the STFT
STFT_dose<-function(dose_level,current_data_frame,total_followup){
  dose_data<-current_data_frame[current_data_frame$dose.level==dose_level,]
  
  STFT<-sum(total_followup-dose_data$patient_followup) #total patient time left
  n<-nrow(dose_data) #number on dose
  c<-sum(dose_data$patient_followup<total_followup) #pending patients
  s<-sum(dose_data$current.DLT)
  return(c(STFT,c,n,s))
}
##INPUT:
#seed=seed for reproducability
#tru_cyc1=true underlying probability of DLT in cycle 1
#co_size=cohort size
#ncohorts=maximum number of cohorts in study
#target=target toxicity (all cycles)
#doses=scaled doses
#a=lower multiplier to target for interval 
#b=upper multiplier to target for interval 
#alpha= prior par 1 (Beta(par 1, par 2))
#beta=prior par 2 (Beta(par 1, par 2))
#nTTP_array= 5x5x5 array of nTTP values for each of the 125 combinations of type and grade of toxicity
#ncycles=number of cycles
#cycle_dec=decrease in P(DLT) by cycle
#dose.skipping=enforce dose skipping rule? T=no dose skipping (default)
#sufficient.information==enforce stopping for sufficient information? no more than 9 patients per dose. T=enforce stopping for sufficient information (default)
#hard.safety=enforce hard safety rule based on Beta(1.1)? T= enforce hard safety (default)
#safety.stopping.low.unsafe= (T=stop when P(p1>0.3)>0.8 (cycle 1))
#safety.stopping.high.toosafe= (T= stop when P(pJ>0.3)>0.8 (cycle 1))
#kfold.skipping: Is a kfold dose skipping rule implemented?
#kfold: the k-fold rise thsat is allowable in dose skipping


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



TITEBOIN_stage2<-function(seed,tru_cyc1,co_size,ncohorts ,target,doses,a,b,alpha,beta,nTTP_array,ncycles,cycle_dec,
                          sufficient.information=T,hard.safety.rule=95,safety.stopping.low.unsafe=T,safety.stopping.high.toosafe=T,
                          kfold.skipping=T,kfold=2){
  set.seed(seed)
  rands<-runif(co_size*ncohorts)
  stop_vec<-rep(0,7)
  ndoses<-length(doses)
  out_mat<-matrix(rep(0,((ncohorts+ncycles))*4),nrow=(ncohorts+ncycles),ncol=4)
  dosevec<-c()
  patient_ID1<-1
  current.time<-0
  
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
  
  
  dose_rec<-NA
  nextdose<-1
  stop<-0

  
  initial<-1
  
  while(stop==0){
    
    
    #DATA GENERATION
    
    if(current.time==0){
      #first cohort
      all.data<-multiple_patients_nTTP(patient_ID1,tru_cyc1,ncycles ,nextdose,entry_time=current.time,nTTP_array = nTTP_array,num_patients = co_size,rands=rands,cycle_dec=cycle_dec)
      current.time<-current.time+1
      patient_ID1<-max(all.data$patient_ID)+1
      current.data<-all.data[all.data$time_of<=current.time,]
      dosevec[current.time]<-nextdose
    }else{
      #subsequent cohorts
      all.data<-rbind(all.data,multiple_patients_nTTP(patient_ID1,tru_cyc1,ncycles ,nextdose,entry_time=current.time,nTTP_array = nTTP_array,num_patients = co_size,rands=rands,cycle_dec=cycle_dec))
      current.time<-current.time+1
      patient_ID1<-max(all.data$patient_ID)+1
      current.data<-all.data[all.data$time_of<=current.time,]
      dosevec[current.time]<-nextdose
    }
    patient.data<-patient_data_frame(all.data)
    
    current.patient.data<-current_patient_data_frame(current_time=current.time,patient.dataframe=all.data,follow_up=ncycles)
    
    
    STFT_current<-STFT_dose(dose_level=nextdose,current_data_frame=current.patient.data,total_followup=ncycles)
    #browser()
    #calculate boundaries for escalation/de-escalation
    current_boundaries<- TITE_BOIN_boundary(c=STFT_current[2],n=STFT_current[3],s=STFT_current[4],phi=target,alpha=alpha,beta=beta,a=a,b=b)
    
    #escalate/de-escalate
    if(STFT_current[1]>current_boundaries[2]){
      nextdose<-min(nextdose+1,ndoses)
    }else if(STFT_current[1]<current_boundaries[1]){
      nextdose<-max(nextdose-1,1)
    }
    
    #k fold dose increase on EXPERIMENTED DOSES
    if(!is.na(nextdose)){
    if(kfold.skipping==T){
      
      prevdoseval<-doses[max(dosevec,na.rm=T)]
      nextdoseval<-doses[nextdose]
     # browser()
      if(nextdoseval>(kfold*prevdoseval)){
        maxdoseval<-kfold*prevdoseval
        nextdose<-max(which(doses<=maxdoseval))
      }
    }
    }
    #stopping based on beta-binomial
    npats_doses<-co_size*tabulate(dosevec,nbins = ndoses)
    
    
    #number of DLTs per dose level
    nDLTs_doses<-sapply(c(1:ndoses),function(x) sum(current.data$DLT[(current.data$dose_level==x)&(current.data$patient_ID>0)]))
    nDLTs_doses_c1<-sapply(c(1:ndoses),function(x) sum(current.data$DLT[(current.data$dose_level==x)&(current.data$patient>0)&(current.data$cycle==1)]))
    
    #safety stopping
    if((safety.stopping.low.unsafe==T)&(npats_doses[1]>0)){
      
      
      cyc1_0.3g<-pbeta(0.3,alpha+nDLTs_doses_c1[1],alpha+beta+npats_doses[1]-nDLTs_doses_c1[1],lower.tail = FALSE)
      if(cyc1_0.3g>0.8){
        stop<-6
        stop_vec[6]<-1
        nextdose<-NA
      #  break
      }
    }
    
    if((safety.stopping.high.toosafe==T)&(npats_doses[ndoses]>0)){
      
      
      cycJ_0.3l<-pbeta(0.3,alpha+nDLTs_doses_c1[ndoses],alpha+beta+npats_doses[ndoses]-nDLTs_doses_c1[ndoses],lower.tail = TRUE)
      if(cycJ_0.3l>0.8){
        stop<-7
        stop_vec[7]<-1
        nextdose<-NA
      #  break
      }
    }
    
    
    
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
        dose_rec<-NA
     #   break
      }
      
      
    }
    
    ##sufficient information 
    if((sufficient.information==T)&(!is.na(nextdose))){
      if(npats_doses[nextdose]>=9){
        stop<-4
        stop_vec[4]<-1
        dose_rec<-nextdose
     #   break
      }
    }
    ##precision stopping
    
    # if(precision.stopping==T){
    #   MTD_vector<-
    #     CV<-mad(MTD_vector, center = median(MTD_vector), constant = 1.4826, 
    #           na.rm = FALSE, low = FALSE, high = FALSE)/median(MTD_vector)
    #   # browser()
    #   if(CV<0.3){
    #     stop_vec[2]<-1
    #     stop<-2
    #     
    #   }
    # }
    
    
    if(nrow(patient.data)==(ncohorts*co_size)){
      #follow up for all patients
      current.data<-all.data
      current.time<-max(current.data$time_of)
      
      post_est<-sapply(c(1:ndoses),function(x) (nDLTs_doses[x]+alpha)/(npats_doses[x]+alpha+beta))
      iso_out<-pava(doses,post_est)
      #choose MTD
      dose_rec<-max(which(abs(iso_out-target)==min(abs(iso_out-target))))
      stop<-3
      stop_vec[3]<-1
    }
    
  }
  
  
  
  #grade matrix out
  grade.matrix.out<-matrix(0,ncol=ndoses,nrow=5)
  for (pat in 1:(max(current.data$patient))){
    pat_max_grade<-max(current.data[current.data$patient_ID==pat,]$max_grade)
    pat_dose<-current.data[current.data$patient_ID==pat,]$dose_level[1]
    grade.matrix.out[pat_max_grade+1,pat_dose]<- grade.matrix.out[pat_max_grade+1,pat_dose]+1
  }
  
  
  #output list
  output<-list(
    dose.rec=dose_rec ,#dose recommendation
    num.pat=max(current.data$patient) ,#: number of patients
    dose.ass=tabulate(dosevec,nbins = ndoses) ,# : number of cohorts per dose (vector)
    stop.code=stop_vec, #: stopping reason 
    num.DLT= sum(current.data$DLT[current.data$patient_ID>0]),#: number of DLTS (total)
    grade.mat=grade.matrix.out ,#: 5xnumdoses matrix. rows for grades 0,1,2,3,4. cols for doses.
    duration= current.time,#: total trial duration (until all recruited patients are fully observed)
    max_admissable=max_admissable # max admissable dose at the end of the trial (has hard safety eliminated any?)
    
  )
  
  
  return(output)
}

