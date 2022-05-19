###############################
##    DATA GENERATION        ##
###############################

##weights matrix
weights.mat<-matrix(c(0,0.5,0.75,1,1.5,
                      0,0.5,0.75,1,1.5,
                      0,0,0,0.5,1),nrow=3,ncol=5,byrow=5)



#gives the nTTP value for each of the 125 (5x5x5) combinations of event
nTTP_array<-array(dim=c(5,5,5))
for(i in 1:5){
  for(j in 1:5){
    for(k in 1:5){
      nTTP_array[i,j,k]<-sqrt((weights.mat[1,i])^2+(weights.mat[2,j])^2+(weights.mat[3,k])^2)/2.5
    }
  }
}



index_125_04<-array(NA,c(5,5,5))
index_125<-array(NA,c(5,5,5))

for(i in 1:5){
  for(j in 1:5){
    for(k in 1:5){
      index_125_04[i,j,k]<-max(c(i,j,k))-1
    }
  }
}
grades<-c()
for(grade in 0:4){
  grades[grade+1]<-sum(index_125_04==grade)
}

index_125[1,1,1]<-1
for(grade in 1:4){
  
  indices<-which(index_125_04==grade,arr.ind = T)
  grades_vec<-c(1:(grades[grade+1]))+sum(grades[1:(grade)])
  for(j in 1:nrow(indices)){
    index_125[indices[j,1], indices[j,2],indices[j,3]]<-grades_vec[j]
  }
}

##in order to generate data to correspond to benchmark, we must vectorize the array
vectorize_array<-function(array_in,INDEX=index_125){
  dims<-c(5,5,5)
  vector_out<-c()
  for(i in 1:dims[1]){
    for(j in 1:dims[2]){
      for(k in 1:dims[3]){
        vector_out[INDEX[i,j,k]]<-array_in[i,j,k]
      }
    }
  }
  
  return(vector_out)
}


#take the index of the vector, gives the indices of the array
index_array<-function(index_in,INDEX=index_125){
  vector_out<-as.vector(which(INDEX==index_in,arr.ind = T))
  return(vector_out)
}

#takes the grade, gives the category for POMM method
grade2cat<-function(grade){
  if(grade==5){
    out<-3
  }else{
    out<-max(grade-1,1)
  }
  return(out)
}


cond_cyc_func<-function(cond_vec){
  length_c<-length(cond_vec)
  p3<-c()
  p3[1]<-cond_vec[1]
  for ( i in 2:length_c){
    p3[i]<-prod(1-cond_vec[1:(i-1)])*cond_vec[i]
  }
  return(sum(p3))
}

##generates the individual array of probabilities of each combination
ind_array_S<-function (tru_cyc1, cycle_num, nextdose,cycle_dec) 
{
  #decreasing for subsequent cycles
  tru_cyc_dose<-tru_cyc1[nextdose]*(cycle_dec^(cycle_num-1))
  #this is for the simple distribution of grades - can be changed depending on requirements
   revcumsum5<-c(1-0.5*tru_cyc_dose*c(0:2),1-0.5*tru_cyc_dose*c(4,5))
  revcumsum5[revcumsum5<0]<-0
  cumsum5 <- rev(revcumsum5)
  div_vec <- c()
  div_vec[1] <- 1
  array_out <- array(NA, dim = c(5, 5, 5))
  array_out[1, 1, 1] <- cumsum5[1]
  for (i in 2:5) {
    div_vec[i] <- (i^3) - sum(div_vec[c(1:(i - 1))])
    array_out[i, , ] <- array_out[, i, ] <- array_out[, , 
                                                      i] <- (cumsum5[i] - cumsum5[i - 1])/div_vec[i]
  }
  return(array_out)
}

## PATIENT DATA FRAMES ##
#Each row is a new cycle. Each patient can have from 1 to 'ncycles' rows, depending on response.
#patient_ID : unique index of patient
#Y_obs : category for POMM method
#cycle_num : cycle number
#dose_level : dose level (not value)
#entry_time : time of entry 
#time_of : time of this cycle's observation in the trial
#DLT : 1=DLT in this cycle , 0=no DLT in this cycle
#nTTP : nTTP value for this cycle
#max_grade : maximum toxicity grade observed in this cycle for this patient


##for direction comparison:

#generates a single patient's COMPLETE outcome as a data-frame
single_patient_nTTP<-function(patient,tru_cyc1,ncycles,thenextdose,entry_time,nTTP_array,rands,cycle_dec){
  
  Y_obs<-c()
  nTTP_obs<-c()
  nTTP_obs_ind<-c()
  #offset conditional vector of DLT probs
  cond_vec_off<-c(0,tru_cyc1[thenextdose]*(cycle_dec^(c(0:(ncycles-2)))))
  
  #reduce the cutoffs for subsequent cycles
  DLT_cutoffs<-c(1,sapply(c(2:ncycles), function(x) 1-cond_cyc_func(cond_vec_off[1:x])))
  
  rand<-rands[patient]
  
  for(i in 1:ncycles){
    individual_array<-ind_array_S(tru_cyc1,cycle_num=i,thenextdose,cycle_dec)
    individual_vector<-vectorize_array(individual_array)
    individual_vector_cum<-cumsum(individual_vector)
    individual_vector_cum<-individual_vector_cum*DLT_cutoffs[i]
    #nTTP_obs is for nTTP
    nTTP_obs_cyc_ind<-min(which(rand<individual_vector_cum))
    nTTP_obs_cyc<-nTTP_array[index_array(nTTP_obs_cyc_ind)[1],index_array(nTTP_obs_cyc_ind)[2],index_array(nTTP_obs_cyc_ind)[3]]
    #Y_obs is for POMM
    Y_obs_cyc<-grade2cat(max(index_array(nTTP_obs_cyc_ind)))
    Y_obs<-c(Y_obs,Y_obs_cyc) #vector of observations for one patient on all cycles
    nTTP_obs<-c(nTTP_obs,nTTP_obs_cyc)
    nTTP_obs_ind<-c(nTTP_obs_ind,nTTP_obs_cyc_ind)
    if(max(Y_obs)==3){
      break
    }
  }
  
  if(max(Y_obs)==3){
    Y_obs<-Y_obs[c(1:(min(which(Y_obs==3))))]
    nTTP_obs<-nTTP_obs[c(1:(min(which(Y_obs==3))))]
    nTTP_obs_ind<-nTTP_obs_ind[c(1:(min(which(Y_obs==3))))]
  }
  out_data<-data.frame(patient,Y_obs,c(1:(length(Y_obs))),thenextdose,entry_time,entry_time+c(1:(length(Y_obs))),
                       as.numeric(Y_obs==3),nTTP_obs,unlist(lapply(nTTP_obs_ind,function(x) max(index_array(x))))-1)
                                                            
  
  # out_data<-data.frame(patient,Y_obs,c(1:(length(Y_obs))),thenextdose,entry_time,entry_time+c(1:(length(Y_obs))),
  #                      as.numeric(Y_obs==3),nTTP_obs,max(unlist(lapply(nTTP_obs_ind,index_array)))-1)
  names(out_data)<-c("patient_ID","Y_obs","cycle_num","dose_level","entry_time","time_of","DLT","nTTP","max_grade")
  return(out_data)
}

#generates multiple patients' COMPLETE outcome as a data-frame
multiple_patients_nTTP<-function(patient_ID1,tru_cyc1,ncycles,thenextdose,entry_time,nTTP_array,num_patients,rands,cycle_dec){
  out_data<-single_patient_nTTP(patient_ID1,tru_cyc1,ncycles,thenextdose,entry_time,nTTP_array,rands,cycle_dec)
  for(i in 2:num_patients){
    out_data<-rbind(out_data,single_patient_nTTP(patient_ID1+i-1,tru_cyc1,ncycles,thenextdose,entry_time,nTTP_array,rands,cycle_dec))
  }
  return(out_data)
}


