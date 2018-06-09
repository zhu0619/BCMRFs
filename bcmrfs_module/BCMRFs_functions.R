#----------------------------------------------------------------------#  
# This script contains the BCMRFs functions
# Author : Lu Zhu
# Email : lzhu@techfak.uni-bielefel.de
#----------------------------------------------------------------------#  

BCMRFs_single_imbalance = function( wA, L, Lname, titer, unknowns, knowns)
{
  #---------------------------------------------------------#
  # INPUT:
  # wA: The weighted symmetric adjacency matrix of the network with a dialog of zero.
  # L: the labellings vector with 0, 1 and -1
  # burnin, niter: MCMC iterations 
  #
  #OUTPUT:
  # res : a list of Probability, SCL labelings, Parameters and Likelihood
  #---------------------------------------------------------#
  
  #---------------------------------------------------------
  # Set some internal parameters	
  #---------------------------------------------------------
  
  # report the parameters for the proposal distribution
  para_recoder= c();
  # report the likelihood
  cpc_recorder = c();
  # report the labelings
  label_recorder = L;
  # imbalance coefficient of binary class  
  Ci = rep(1,length(L))
  # record  
  initZlength = 1000; # NB OF SAMPLER FOR z
  Z.schd = seq(from = 1, to = titer, by = 1);	
  
  # Degree of each protein in the PPI network
  dwA = rowSums(wA);
  names(dwA) = rownames(wA);
  
  #------------------------------------------------------------------------------------------#
  # Initialize MRF parameters by logistic regression : linear model on all SCL known proteins!
  #------------------------------------------------------------------------------------------#
  
  # nb of interacting partner in this location
  K1 = as.vector(wA[knowns,knowns] %*% L[knowns]);
  # nb of interacting partnerS
  NS = as.vector(dwA[knowns]);  
  # nb of interacting partner NOT in this location
  K0 = NS-K1; 
  #  SCL of known proteins
  Lk = as.vector(L[knowns]);
  
  
  regtable = as.data.frame(cbind(Lk, K1, K0));
  colnames(regtable) = c("SCL", "K1", "K0");
  
  if(sum(K1)!=0){
    regtable.fit = brglm(regtable$SCL ~ regtable$K1 + regtable$K0,
                         family=binomial(logit), 
                         method = "brglm.fit",nIter = 5000,data = regtable );
  }else{
    regtable.fit = brglm(regtable$SCL ~ regtable$K0,
                         family=binomial(logit), 
                         method = "brglm.fit",nIter = 5000,data = regtable );
  }
  
  if(any(is.na(regtable.fit$coefficients)))
    warning(Lname, ": brglm has NA coefficients - ",
            paste0(names(regtable.fit$coefficients)[is.na(regtable.fit$coefficients)], collapse=", "),
            call.=FALSE
    )
  
  if (!regtable.fit$converged) {
    warning(Lname, ": brglm fitting did not converge", call.=FALSE)
  }
  
  
  Z = rmvnorm(mean = regtable.fit$coefficients,sigma=vcov(regtable.fit), n=initZlength);
  para_names = colnames(Z)
  
  MRFparams = as.vector(regtable.fit$coefficients);
  
  #--------------------------------------------------------------------------------#
  #Initialization of Labelling using the model above 
  #--------------------------------------------------------------------------------#
  
  # a matrix with unknown proteins with their interacting partners(only SCL known proteins)
  wAuk = wA[unknowns,knowns]; 
  K1uk = as.vector(wAuk %*% L[knowns]);
  NSuk = dwA[unknowns];
  K0uk = NSuk-K1uk;
  
  v = MRFparams[1] + K1uk*MRFparams[2] + K0uk*MRFparams[3];
  
  P = (1/(1 + exp(-v)))  
  L[unknowns[P >= pcut]] = 1;
  L[unknowns[P  < pcut]] = 0;
  
  # imblance level of SCL labelings
  Lc_single = (length(L)-sum(L))/sum(L)
  
  Ci[which(L==1)] = Lc_single
  Ci[which(L==0)] = 1
  
  #--------------------------------------------------------------------------------#
  # MRF simulations
  #--------------------------------------------------------------------------------#
  
  # sigma for parameter sampling
  e.sigma = matrix(ncol= ncol(Z), nrow= ncol(Z));
  e.sigma[] = 0;
  diag(e.sigma) = 0.0001;	
  
  # a matrix with unknown proteins with theirinteracting partners 
  # (all, including initialized unknown proteins)
  wAuk = wA[unknowns,];
  K1uk = as.vector(wAuk %*% L);
  NSuk = dwA[unknowns];
  K0uk = NSuk-K1uk;
  
  # a vector store all the probablity values during the simulation
  probs = vector(mode = "numeric", length = length(L));  
  names(probs) = names(L);
  probs[]  = 0;
  counter  = 0;
  it = 0
  temp= 10000
  for(t in 1:titer)
  {
    # cat("Cycle:",t,"with Temperature:",temp,"\n")
    para_recoder_x = c()
    cpc_recorder_x = c()
    
    #--------------------------------------------------#
    # update the labeling
    #--------------------------------------------------#
    # estimate the unknown protein with ALL the interacting  proteins
    K1uk = as.vector(wAuk %*% L); 
    K0uk = NSuk-K1uk;
    
    v =  MRFparams[1] + K1uk*MRFparams[2] + K0uk*MRFparams[3];
    
    P = (1/(1 + exp(-v))) 
    L[unknowns[P >= pcut]] = 1;
    L[unknowns[P < pcut]] = 0;
    Lc_single = (length(L)-sum(L))/sum(L)
    
    # update the imbalance coefficient
    Ci[which(L==1)] = Lc_single
    Ci[which(L==0)] = 1
    
    # times of sample parameters without optimization
    ntt = 1 
    for(tt in 1: ntt){
      #--------------------------------------------------#
      # Update MRFparams. Propose a candidate		
      #--------------------------------------------------#
      s = sample(1:nrow(Z),2, replace=F);
      e = rmvnorm(n=1, mean = c(rep(0,ncol(Z))), sigma = e.sigma);
      # nb of parameters in the model 
      d = length(MRFparams)
      gamma_s=2.38/sqrt(2*d)
      gamma = runif(min = gamma_s/2, max = gamma_s,n=1);
      
      MRFparamsP = as.vector(MRFparams + (gamma*(Z[s[1],] - Z[s[2],])) + e);
      
      K1 = as.vector (wA %*% L);
      NS = as.vector(rowSums(wA));
      K0 = NS - K1;
      
      vc = as.vector(MRFparams[1]  + K1*MRFparams[2] +  K0*MRFparams[3] );
      vp = as.vector(MRFparamsP[1] + K1*MRFparamsP[2] + K0*MRFparamsP[3] );		
      
      # PLF for current parameters
      # cp is the conditional probablity for xi=1 given it's neighbors
      cpc = 1/(1 + exp(-vc));			    
      cpc.1 = cpc;
      # update the information of proteins which have no evidence located in this SCL.
      cpc[L == 0] = 1 - cpc[L == 0];	
      cpc[cpc < 1e-12] = 1e-12;
      cpc = `^`(cpc,Ci)
      logc = log(cpc);
      
      #PLF for proposed parameters
      cpp = 1/(1 + exp(-vp));			
      cpp.1 = cpp;			
      cpp[L == 0] = 1 - cpp[L == 0];
      cpp[cpp < 1e-12] = 1e-12;
      cpp = `^`(cpp,Ci)
      logp = log(cpp);
      
      if(acceptanceProbability(sum(logc),sum(logp),temp)>runif(1))
      {
        # accept the proposal and update the parameters
        MRFparams = MRFparamsP
        # record parameters and PL
        it=it+1
        para_recoder_x = rbind(para_recoder_x,c(MRFparams,sum(logp),it));  
        cpc_recorder_x = rbind(cpc_recorder_x, cpp.1) 
      }
      
    }
    #--------------------------------------------------#
    if(length(para_recoder_x)!=0){
      tt_max_PLF = max(para_recoder_x[,4])
      
      best_iter_x = max(which(para_recoder_x[,4]==tt_max_PLF))
      best_x = para_recoder_x[best_iter_x,]
      
      cpc_recorder = rbind(cpc_recorder, cpc_recorder_x[best_iter_x,]) 
      para_recoder = rbind(para_recoder,c(best_x[1:4],t))
      
      MRFparams = best_x[1:3]
      if(is.element(t,Z.schd)){Z = rbind(Z, MRFparams);
      }
      label_recorder = rbind(label_recorder, L) 
    }
    
    # cooling
    if(temp>1){
      temp=temp*(1-0.08)
    }
  }
  
  colnames(para_recoder) = c(para_names,"log_PLF","iteration");
  paranames = para_names
  
  #------------------------------------------------------------------
  # prepare output
  #------------------------------------------------------------------  
  
  # take the average of the last iterations of simulated annealing
  iter_record = 20
  last_iter= seq(nrow(para_recoder)-iter_record,nrow(para_recoder))
  bestPLF = round(mean(para_recoder[last_iter,"log_PLF"]),3)
  bestProbs = colMeans(cpc_recorder[last_iter,])
  bestMRFparams  = round(colMeans(para_recoder[last_iter,para_names]),3) 
  
  # calibrate the final probabilities
  P = calibrate(bestProbs)
  
  names(P) = names(L)
  pcut = 0.5  
  L[unknowns[P[unknowns] >= pcut]] = 1;
  L[unknowns[P[unknowns] < pcut]] = 0;
  
  res = list(P,L,bestMRFparams,bestPLF)
  names(res) = c('Probability','SCL_label','Parameters','Likelihood')
  return(res)
}


BCMRFs_adj_imbalance = function(wA, L, Adjs,Lname, titer,L_all ,unknowns,knowns)
{
  #---------------------------------------------------------#
  # INPUT:
  # wA: The weighted symmetric adjacency matrix of the network with a dialog of zero.
  # L: the labellings vector with 0, 1 and -1
  # burnin, niter: MCMC iterations 
  #
  #OUTPUT:
  # res : a list of Probability, SCL labelings, Parameters and Likelihood
  #---------------------------------------------------------#
  
  #---------------------------------------------------------
  # Set some internal parameters	
  #---------------------------------------------------------
  
  # report the parameters for the proposal distribution
  para_recoder= c();
  # report the likelihood
  cpc_recorder = c();
  # report the labelings
  label_recorder = L;
  # imbalance coefficient of binary class  
  Ci = rep(1,length(L))
  # record  
  initZlength = 1000; # NB OF SAMPLER FOR z
  Z.schd = seq(from = 1, to = titer, by = 1);	
  
  # Degree of each protein in the PPI network
  dwA = rowSums(wA);
  names(dwA) = rownames(wA);
  
  
  #------------------------------------------------------------------------------------------#
  # Initialize MRF parameters by logistic regression : linear model on all SCL (initialized)!
  #------------------------------------------------------------------------------------------#
  
  # nb of interacting partner in this location
  K1 = as.vector(wA[knowns,knowns] %*% L[knowns]);
  # nb of interacting partnerS
  NS = as.vector(dwA[knowns]);  
  # nb of interacting partner NOT in this location
  K0 = as.vector(NS-K1); 
  #  SCL of known proteins
  Lk = as.vector(L[knowns]); 
  
  if(length(Adjs)>1){
    Adjs = Adjs[colSums(L_all[knowns,Adjs])>0]
  }else{
    Adjs = Adjs[sum(L_all[knowns,Adjs])>0]
  }
  
  Adjs_2= c()
  Kh=c()
  names=c()
  for(h in 1:length(Adjs)){
    Kadj= as.vector(wA[knowns,knowns] %*% L_all[knowns,Adjs[h]])
    if(sum(Kadj)>0){
      Kh=cbind(Kh,Kadj)
      names=c(names,paste(c("Kh_",h),collapse = ""))
      Adjs_2 = c(Adjs_2,Adjs[h])
    }
  }
  
  Adjs = Adjs_2
  
  regtable = cbind(Lk, K1, K0, Kh);
  regtable = as.data.frame(regtable)
  colnames(regtable) = c("SCL", "K1", "K0",names)
  
  model_command =" regtable.fit = brglm(regtable$SCL ~ regtable$K1 + regtable$K0 "
  
  for(h in names){
    model_command =  paste(c(model_command,"+ regtable$",h),collapse = "")
  }
  
  
  model_command = paste(model_command,",family=binomial(link = 'logit'),method = 'brglm.fit')")
  
  eval(parse(text=model_command))
  
  coef = regtable.fit$coefficient;
  na = which(is.na(coef)==T)
  if(any(na)){
    coef = regtable.fit$coefficient[-na]
    names = names[-na]
    Adjs = Adjs [-(na-3)]
  }
  
  sigma = vcov(regtable.fit)
  Z = rmvnorm(mean = coef,sigma=sigma, n=initZlength);
  MRFparams = as.vector(coef);
  para_names =  colnames(Z)
  
  #--------------------------------------------------------------------------------#
  # MRF analysis
  #--------------------------------------------------------------------------------#
  ### sigma for parameter sampling
  e.sigma = matrix(ncol= ncol(Z), nrow= ncol(Z));
  e.sigma[] = 0;
  diag(e.sigma) = 0.0001;	
  
  # a matrix with unknown proteins with theirinteracting partners
  # (all, including initialized unknown proteins)
  wAuk = wA[unknowns,];
  K1uk = as.vector(wAuk %*% L);
  NSuk = dwA[unknowns];
  K0uk = NSuk-K1uk;
  
  # a vector store all the probablity values during the simulation
  probs = vector(mode = "numeric", length = length(L));  
  names(probs) = names(L);
  probs[]  = 0;
  counter  = 0;
  
  Khuk=c()
  names=c()
  for(h in 1:length(Adjs)){
    Kadjuk= as.vector(wAuk %*% L_all[,Adjs[h]])
    Khuk=cbind(Khuk,Kadjuk)
    names=c(names,paste(c("Kh_",h,"uk"),collapse = ""))
  }
  
  it = 0
  temp= 1000
  coolrate = 0.04
  
  for(t in 1:titer)
  {
    #print(t)
    para_recoder_x = c()
    cpc_recorder_x = c()
    
    #--------------------------------------------------#
    # update the labeling
    #--------------------------------------------------#
    # estimate the unknown protein with ALL the interacting  proteins
    K1uk = as.vector(wAuk %*% L); 
    K0uk = NSuk-K1uk;
    # Kt1uk = as.vector(tAuk %*% L);  
    
    V_command = " v= MRFparams[1] + K1uk*MRFparams[2] + K0uk*MRFparams[3]" 
    
    for(h in 1:length(names)){
      V_command =  paste(c(V_command," + Khuk[,",h,"]*MRFparams[",(h+3),"]"),collapse = "")
    }
    
    
    eval(parse(text=V_command))
    
    P = (1/(1 + exp(-v))) 
    L[unknowns[P >= pcut]] = 1;
    L[unknowns[P  < pcut]] = 0;	
    
    Lc_single = (length(L)-sum(L))/sum(L)
    
    Ci[which(L==1)] = Lc_single
    Ci[which(L==0)] = 1
    
    for(tt in 1: 5){
      #--------------------------------------------------#
      # Update MRFparams. Propose a candidate		
      #--------------------------------------------------#
      s = sample(1:nrow(Z),2, replace=F);
      
      e = rmvnorm(n=1, mean = c(rep(0,ncol(Z))), sigma = e.sigma);
      
      d = length(MRFparams)# nb of parameters in the model 
      gamma_s=2.38/sqrt(2*d)
      gamma = runif(min = gamma_s/2, max = gamma_s,n=1);
      
      MRFparamsP = as.vector(MRFparams + (gamma*(Z[s[1],] - Z[s[2],])) + e);
      
      K1 = as.vector (wA %*% L);
      NS = as.vector(rowSums(wA));
      K0 = NS - K1;
      
      Kh=c()
      names=c()
      for(h in 1:length(Adjs)){
        Kadj= as.vector(wA %*% L_all[,Adjs[h]])
        Kh=cbind(Kh,Kadj)
        names=c(names,paste(c("Kh_",h),collapse = ""))
      }
      
      vc_command = " vc= MRFparams[1] + K1*MRFparams[2] + K0*MRFparams[3]" 
      vp_command = " vp= MRFparamsP[1] + K1*MRFparamsP[2] + K0*MRFparamsP[3]" 
      
      for(h in 1:length(names)){
        vc_command =  paste(c(vc_command," + Kh[,",h,"]*MRFparams[",(h+3),"]"),collapse = "")
        vp_command =  paste(c(vp_command," + Kh[,",h,"]*MRFparamsP[",(h+3),"]"),collapse = "")
      }
      
      eval(parse(text=vc_command))
      eval(parse(text=vp_command))
      
      
      #PLF for current parameters
      # cp is the conditional probablity for xi=1 given it's neighbors
      cpc = 1/(1 + exp(-vc));			 
      cpc.1 = cpc;			
      # update the information of proteins which have no evidence located in this SCL.
      cpc[L == 0] = 1 - cpc[L == 0];
      cpc[cpc < 1e-12] = 1e-12;
      cpc = `^`(cpc,Ci)
      logc = log(cpc);
      
      
      #PLF for proposed parameters
      cpp = 1/(1 + exp(-vp));			
      cpp.1 = cpp;			
      cpp[L == 0] = 1 - cpp[L == 0];
      cpp[cpp < 1e-12] = 1e-12;
      cpp = `^`(cpp,Ci)
      logp = log(cpp);
      
      if(temp>=0)
      {
        #print(paste("temperature is ",temp,"\n "))
        if(acceptanceProbability(sum(logc),sum(logp),temp)>runif(1) )
        {
          MRFparams = MRFparamsP
          #print("update parameter!");
          it=it+1
          para_recoder_x = rbind(para_recoder_x,c(MRFparams,sum(logc),it));  
          cpc_recorder_x = rbind(cpc_recorder_x, cpp.1) 
        }
      }
      
    }
    #-----------------------#
    if(!is.null(para_recoder_x)){
      colnames(para_recoder_x) = c(para_names,'PLF','iter')
      best_iter_x = max(which.max(para_recoder_x[,'PLF']))
      best_x = para_recoder_x[best_iter_x,]
      
      cpc_recorder = rbind(cpc_recorder, cpc_recorder_x[best_iter_x,]) 
      para_recoder = rbind(para_recoder,c(best_x[c(para_names,'PLF')],t))
      
      MRFparams = best_x[para_names]
      if(is.element(t,Z.schd)){Z = rbind(Z, MRFparams);}
    }
    
    
    # cooling
    if(temp>1){
      temp=temp*(1-coolrate)
    }
  }
  
  #-----------------------#
  # take the parameter set which the maximal log PLFs
  #-----------------------#
  colnames(para_recoder) = c(para_names,"log_PLF","iteration");
  paranames = para_names
  
  #------------------------------------------------------------------
  # prepare output
  #------------------------------------------------------------------  
  
  # take the average of the last iterations of simulated annealing
  iter_record = 20
  last_iter= seq(nrow(para_recoder)-iter_record,nrow(para_recoder))
  bestPLF = round(mean(para_recoder[last_iter,"log_PLF"]),3)
  bestProbs = colMeans(cpc_recorder[last_iter,])
  bestMRFparams  = round(colMeans(para_recoder[last_iter,para_names]),3) 
  
  
  P = calibrate(bestProbs)
  names(P) = names(L)
  
  pcut = 0.5  
  L[unknowns[P[unknowns] >= pcut]] = 1;
  L[unknowns[P[unknowns] < pcut]] = 0;
  
  res = list(P,L,bestMRFparams,bestPLF)
  names(res) = c('Probability','SCL_label','Parameters','Likelihood')
  return(res)
}



BMRF_adj_feature_imbalance = function(wA, L, Adjs,Lname,titer,feature,L_all,unknowns,knowns)
{
  #---------------------------------------------------------#
  # INPUT:
  # wA: The weighted symmetric adjacency matrix of the network with a dialog of zero.
  # L: the labellings vector with 0, 1 and -1
  # burnin, niter: MCMC iterations 
  #
  #OUTPUT:
  # res : a list of Probability, SCL labelings, Parameters and Likelihood
  #---------------------------------------------------------#
  
  #---------------------------------------------------------
  # Set some internal parameters	
  #---------------------------------------------------------
  
  # report the parameters for the proposal distribution
  para_recoder= c();
  # report the likelihood
  cpc_recorder = c();
  # report the labelings
  label_recorder = L;
  # imbalance coefficient of binary class  
  Ci = rep(1,length(L))
  # record  
  initZlength = 1000; # NB OF SAMPLER FOR z
  Z.schd = seq(from = 1, to = titer, by = 1);	
  
  # Degree of each protein in the PPI network
  dwA = rowSums(wA);
  names(dwA) = rownames(wA);
  
  feat = feature
  #------------------------------------------------------------------------------------------#
  # Initialize MRF parameters by logistic regression : linear model on all SCL (initialized)!
  #------------------------------------------------------------------------------------------#
  
  # for one SCL L[,i]
  # nb of interacting partner in this location for the tissue
  K1 = as.vector(wA[knowns,knowns] %*% L[knowns]);
  # nb of interacting partnerS
  NS = as.vector(dwA[knowns]);
  # nb of interacting partner NOT in this location
  K0 = as.vector(NS-K1); 
  #  SCL of known proteins
  Lk = as.vector(L[knowns]); 
  
  fk = feature[knowns]
  nbf = ncol(data.frame(fk))
  
  if(!is.null(Adjs)){
    if(length(Adjs) >1){
      Adjs = Adjs[colSums(L_all[knowns,Adjs])>0]
    }else{
      Adjs = Adjs[sum(L_all[knowns,Adjs])>0]
    }
    
    Kh=c()
    names=c()
    Adjs_2= c()
    for(h in 1:length(Adjs)){
      Kadj= as.vector(wA[knowns,knowns] %*% L_all[knowns,Adjs[h]])
      if(sum(Kadj)>0){
        Kh=cbind(Kh,Kadj)
        names=c(names,paste(c("Kh_",h),collapse = ""))
        Adjs_2 = c(Adjs_2,Adjs[h])
      }
    }
    
    Adjs = Adjs_2
    
    regtable = cbind(Lk, K1, K0, Kh, fk);
  }else{
    regtable = cbind(Lk, K1, K0, fk);
  }
  regtable = as.data.frame(regtable)
  colnames(regtable) = c("SCL", "K1", "K0",names, "features")
  
  model_command =" regtable.fit = brglm(regtable$SCL ~ regtable$K1 + regtable$K0 "
  
  if(!is.null(Adjs)){
    for(h in names){
      model_command =  paste(c(model_command,"+ regtable$",h),collapse = "")
    }
  }
  
  
  model_command =  paste(c(model_command," + regtable$features"),collapse = "")
  
  
  model_command = paste(model_command,",family=binomial(link = 'logit'),method = 'brglm.fit')")
  
  eval(parse(text=model_command))
  
  coef = regtable.fit$coefficient;
  na = which(is.na(coef)==T)
  if(any(na)){
    coef = regtable.fit$coefficient[-na]
    names = names[-na]
    Adjs = Adjs [-(na-3)]
    if(2 %in% na){
      return(NULL)
    }
    
  }
  
  sigma = vcov(regtable.fit)
  Z = rmvnorm(mean = coef,sigma=sigma, n=initZlength);
  MRFparams = as.vector(coef);
  para_names =  colnames(Z)
  
  #--------------------------------------------------------------------------------#
  # MRF simulations
  #--------------------------------------------------------------------------------#
  ### sigma for parameter sampling
  e.sigma = matrix(ncol= ncol(Z), nrow= ncol(Z));
  e.sigma[] = 0;
  diag(e.sigma) = 0.0001;	
  # a matrix with unknown proteins with theirinteracting partners
  #all, including initialized unknown proteins)
  wAuk = wA[unknowns,];
  K1uk = as.vector(wAuk %*% L);
  NSuk = dwA[unknowns];
  K0uk = NSuk-K1uk;
  fuk = feature[unknowns]
  
  # a vector store all the probablity values during the simulation
  probs = vector(mode = "numeric", length = length(L));  
  names(probs) = names(L);
  probs[]  = 0;
  counter  = 0;
  
  names_uk=c()
  if(!is.null(Adjs)){
    Khuk=c()
    for(h in 1:length(Adjs)){
      Kadjuk= as.vector(wAuk %*% L_all[,Adjs[h]])
      Khuk=cbind(Khuk,Kadjuk)
      names_uk=c(names_uk,paste(c("Kh_",h,"uk"),collapse = ""))
    }
  }
  
  it = 0
  temp= 10000
  coolrate = 0.04
  
  for(t in 1:titer)
  {
    #print(t)
    para_recoder_x = c()
    cpc_recorder_x = c()
    
    #--------------------------------------------------#
    # update the labeling
    #--------------------------------------------------#
    # estimate the unknown protein with ALL the interacting  proteins
    K1uk = as.vector(wAuk %*% L); 
    K0uk = NSuk-K1uk;
    # Kt1uk = as.vector(tAuk %*% L);  
    
    V_command = " v= MRFparams[1] + K1uk*MRFparams[2] + K0uk*MRFparams[3]" 
    
    if(!is.null(Adjs)){
      for(h in 1:length(names_uk)){
        V_command =  paste(c(V_command," + Khuk[,",h,"]*MRFparams[",(h+3),"]"),collapse = "")
      }
    }else{
      h=0
    }
    
    V_command =  paste(c(V_command,"+ fuk*MRFparams[",(h+4),"]"),collapse = "")
    
    
    eval(parse(text=V_command))
    
    P = (1/(1 + exp(-v))) 
    L[unknowns[P >= pcut]] = 1;
    L[unknowns[P  < pcut]] = 0;	
    
    Lc_single = (length(L)-sum(L))/sum(L)
    
    Ci[which(L==1)] = Lc_single
    Ci[which(L==0)] = 1
    
    # times of sample parameters without optimization
    ntt = 1 
    for(tt in 1: ntt){
      #--------------------------------------------------#
      # Update MRFparams. Propose a candidate		
      #--------------------------------------------------#
      s = sample(1:nrow(Z),2, replace=F);
      
      e = rmvnorm(n=1, mean = c(rep(0,ncol(Z))), sigma = e.sigma);
      # nb of parameters in the model 
      d = length(MRFparams)
      gamma_s=2.38/sqrt(2*d)
      gamma = runif(min = gamma_s/2, max = gamma_s,n=1);
      
      MRFparamsP = as.vector(MRFparams + (gamma*(Z[s[1],] - Z[s[2],])) + e);
      
      K1 = as.vector (wA %*% L);
      NS = as.vector(rowSums(wA));
      K0 = NS - K1;
      
      Kh=c()
      names=c()
      if(!is.null(Adjs)){
        for(h in 1:length(Adjs)){
          Kadj= as.vector(wA %*% L_all[,Adjs[h]])
          Kh=cbind(Kh,Kadj)
          names=c(names,paste(c("Kh_",h),collapse = ""))
        }
      }
      vc_command = " vc= MRFparams[1] + K1*MRFparams[2] + K0*MRFparams[3]" 
      vp_command = " vp= MRFparamsP[1] + K1*MRFparamsP[2] + K0*MRFparamsP[3]" 
      
      if(!is.null(Adjs)){
        for(h in 1:length(names)){
          vc_command =  paste(c(vc_command," + Kh[,",h,"]*MRFparams[",(h+3),"]"),collapse = "")
          vp_command =  paste(c(vp_command," + Kh[,",h,"]*MRFparamsP[",(h+3),"]"),collapse = "")
        }
      }else{h=0}
      
      
      vc_command =  paste(c(vc_command,"+ feat*MRFparams[",(h+4),"]"),collapse = "")
      vp_command =  paste(c(vp_command,"+ feat*MRFparamsP[",(h+4),"]"),collapse = "")
      
      eval(parse(text=vc_command))
      eval(parse(text=vp_command))
      
      
      #PLF for current parameters
      # cp is the conditional probablity for xi=1 given it's neighbors
      cpc = 1/(1 + exp(-vc));			 
      cpc.1 = cpc;			
      # update the information of proteins which have no evidence located in this SCL.
      cpc[L == 0] = 1 - cpc[L == 0];	
      cpc[cpc < 1e-12] = 1e-12;
      cpc = `^`(cpc,Ci)
      logc = log(cpc);
      
      
      #PLF for proposed parameters
      cpp = 1/(1 + exp(-vp));			
      cpp.1 = cpp;			
      cpp[L == 0] = 1 - cpp[L == 0];
      cpp[cpp < 1e-12] = 1e-12;
      cpp = `^`(cpp,Ci)
      logp = log(cpp);
      
      #decide if accept proposed parameters or not 
      if(temp>=0)
      {
        #print(paste("temperature is ",temp,"\n "))
        if(acceptanceProbability(sum(logc),sum(logp),temp)>runif(1) )
        {
          MRFparams = MRFparamsP
          #print("update parameter!");
          it=it+1
          para_recoder_x = rbind(para_recoder_x,c(MRFparams,sum(logc),it));  
          cpc_recorder_x = rbind(cpc_recorder_x, cpp.1) 
        }
      }
      
    }
    #-----------------------------------------------------------------
    if(!is.null(para_recoder_x)){
      colnames(para_recoder_x) = c(para_names,'PLF','iter')
      best_iter_x = max(which.max(para_recoder_x[,'PLF']))
      best_x = para_recoder_x[best_iter_x,]
      
      cpc_recorder = rbind(cpc_recorder, cpc_recorder_x[best_iter_x,]) 
      para_recoder = rbind(para_recoder,c(best_x[c(para_names,'PLF')],t))
      MRFparams = best_x[para_names]
      if(is.element(t,Z.schd)){Z = rbind(Z, MRFparams);}
    }
    
    
    # cooling
    if(temp>1){
      temp=temp*(1-coolrate)
    }
  }
  
  colnames(para_recoder) = c(para_names,"log_PLF","iteration");
  paranames = para_names
  
  #------------------------------------------------------------------
  # prepare output
  #------------------------------------------------------------------  
  
  # take the average of the last iterations of simulated annealing
  iter_record = 20
  last_iter= seq(nrow(para_recoder)-iter_record,nrow(para_recoder))
  bestPLF = round(mean(para_recoder[last_iter,"log_PLF"]),3)
  bestProbs = colMeans(cpc_recorder[last_iter,])
  bestMRFparams  = round(colMeans(para_recoder[last_iter,para_names]),3) 
  
  
  P = calibrate(bestProbs)
  names(P) = names(L)
  
  pcut = 0.5  
  L[unknowns[P[unknowns] >= pcut]] = 1;
  L[unknowns[P[unknowns] < pcut]] = 0;
  
  res = list(P,L,bestMRFparams,bestPLF)
  names(res) = c('Probability','SCL_label','Parameters','Likelihood')
  return(res)
}


