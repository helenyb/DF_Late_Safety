##prior calibrations

load("TITE_mTPI2_priorcal_stage1.RData")
method<-"TITE_mTPI2"
prior_cal_scen_list<-list(P_s1,P_s2,P_s3,P_s4)
#number of prior parameter combinations
num_par<-length(get(paste(c(method,"_prior_par_list"),collapse = "")))

#which dose is correct in each scenario
correct_list<-unlist(lapply(prior_cal_scen_list, function(x) which.min(abs(0.3-x))))

#empty results matrix
result_matrix<-matrix(NA, nrow=num_par,ncol=length(prior_cal_scen_list))      

#fill results matrix
for (scen in 1:length(prior_cal_scen_list)){
  if(all(prior_cal_scen_list[[scen]]>0.3)){
    #safety stopping rules for all dose unsafe, safety code 5 or 6  
    for(parcomb in 1:num_par){
      #count the number of simulations recommending correct dose
      result<-get(paste(c(method,"_priorcal_s",scen,"_p",parcomb),collapse = ""))
      stop_list<-(lapply(result,'[[','stop.code'))
      stop_matrix<-do.call(rbind,stop_list)
      result_matrix[parcomb,scen]<-sum(rowSums(stop_matrix[,c(5,6)],na.rm=T)>0)
    }
  }else{
    
    for(parcomb in 1:num_par){
      #count the number of simulations recommending correct dose
      result<-get(paste(c(method,"_priorcal_s",scen,"_p",parcomb),collapse = ""))
      result_matrix[parcomb,scen]<-sum(unlist(lapply(result,"[[",'dose.rec'))==correct_list[scen],na.rm=T)
    }
  }
}


#proportion
result_matrix<-result_matrix/length(result)
#geometric mean
geomean_vec<-exp(rowMeans(log(result_matrix)))
#choice of parameters (index)
geomean_choice<-which.max(geomean_vec)
#choice of parameters (parameters)
par_choice<- get(paste(c(method,"_prior_par_list"),collapse = ""))[[geomean_choice]]
