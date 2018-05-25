# basic functions
adjacentMatrixFile = function(file){
  adjMatrix= as.matrix(read.csv(file,sep=";",row.names=1))
  return(adjMatrix) 
}


loadProteinSCL = function(FILE, A, S)
{
  data =  unique(read.table(FILE, sep = S,header = T));
  proteins = as.character(rownames(A));
  #print("proteins")
  #print(proteins)
  
  ids = 1:length(proteins)
  names(ids) = proteins;
  p = as.character(data[,1]);
  
  # if the annotated proteins are more than the protein in the network
  data =data[which(is.na(ids[p])==F),]
  
  ann = as.character(data[,2]);
  ann = unique(ann);
  ann = sort(ann);
  
  ida = 1:length(ann);
  names(ida) = ann;    
  
  p = as.character(data[,1]);
  a = as.character(data[,2]);
  
  iv = ids[p];
  jv = ida[a];
  
  L = sparseMatrix(i = iv, j = jv , x = 1, dims = c(length(proteins), length(ann)));	
  L@x[] = 1;
  
  rownames(L) = proteins;
  colnames(L) = ann;
  print("loadAnn... ok");	
  ## L is a binary matrix of proteins and their annotation
  return(as.matrix(L));
}

rebalance = function(L,prob){
  
  Lc_single = (length(L)-sum(L))/sum(L)
  
  Ci[which(L==1)] = Lc_single
  Ci[which(L==0)] = 1
}


acceptanceProbability = function(PLFc, PLFp, temp)
{
  if (PLFp >= PLFc) {
    return(1)
  }else{
    r=exp((PLFp-PLFc)/(temp))
    #print(c(temp,PLFp-PLFc,r,runif(1)))
    return(r)
  }
}

calibrate = function(P)
{
  #Calibrate the input probabilities
  P[P >= 1 - (1e-10)] = 1- (1e-10);
  P[P <= (1e-10)] =  (1e-10);
  Ppriors = mean(P);
  logitP = log(P/(1-P))
  logitPpriors = log(Ppriors/(1-Ppriors));
  a = 2;
  P2 = logitPpriors + a*(logitP - logitPpriors);
  P2 = exp(P2)/(1+exp(P2));
  P = P2;
  return(P);
  
}

GiantCluster = function(gh)
{
  g=as.undirected(gh)
  clu <- components(g)
  cluster = subgraph(g, groups(clu)$'1')
  return(cluster)
}
