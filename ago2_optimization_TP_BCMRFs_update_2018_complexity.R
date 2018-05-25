
  pkg_to_require = c("ROCR","Matrix","brglm","mvtnorm","glmnet","plyr","igraph",'utiml','mygene')
  pkg_to_install = setdiff(pkg_to_require, rownames(installed.packages()))
if (length(pkg_to_install) > 0) {
  for(pkg_name in pkg_to_install){
    print(paste0('Installing package ',pkg_name))
    install.packages(pkg_name)
  }
}
for(pkg_name in pkg_to_require){
  print(paste0('requring package ',pkg_name))
  require(pkg_name,character.only = TRUE)
}

  
source("bcmrfs_module/basic_functions.R")
# library(ROCR);
# library(Matrix);
# library(brglm);
# library(mvtnorm); 
# library(glmnet);
# library(plyr);
# library(igraph);
# library(utiml)
# library(mldr)
# library(mygene)
# #library(GO.db);
# 
# 
# 
# source("bmrf_functions.R");
# source("BMRF_single.R");
# source("BMRF_single_feature.R");
# source("BMRF_adj_2.R");
# source("BMRF_adj_known.R");
# source("BMRF_adj_known_feature.R");
# source("bmrf_adj_feature.R");
# source("evaluation_sheet.R")
# source('BMRF_single.R')
# source('BMRF_single_imbalance.R')
# source('BMRF_single_imbalance_K0.R')
# source('BMRF_single_tissue_nv.R')
# source('BMRF_single_tissue_ct.R')
# source('BMRF_adj_known_feature_imbalance_tissue.R')
# # source('BMRF_adj_feature_imbalance_tissue.R')
# source('BMRF_adj_feature_imbalance.R')
# source('BMRF_adj_feature_imbalance_K0.R')
# source('BMRF_adj_imbalance_tissue.R')
# source('BMRF_adj_imbalance.R')
# 
# source('BMRF_adj_imbalance_tissue_2.R')
#   

#-----------------------------------------------#
# Preparation
#-----------------------------------------------#

## network preparation
##-----------------------------#
net_T_file = "test_dataset/test.gml"
net_T = read_graph(net_T_file, format = "graphml")
# option: select the largest subnetwork in case of unconnected graph
#A_cluster = GiantCluster(net_T)

A_ori = get.adjacency(net_T, type="both",attr=NULL, names=TRUE, sparse=TRUE)
# option: weighted network
# wA_ori = get.adjacency(net_T, type="both",attr='weight', names=TRUE, sparse=TRUE)[A_cluster,A_cluster]

# remove the self-interacting PPI
diag(A_ori)=0
#diag(wA_ori)=0

# load the spatial adjacent matrix of subcellular compartments
# user can customize the matrix for their needs
adjMatrix = adjacentMatrixFile("test_dataset/test_adjMatrix.csv")


## features
##-----------------------------#
 # features = read.csv("all_features_GO.csv",sep=";",row.names = 1,check.names=FALSE)
# features_ori = unique(read.csv("260917_no_GO_feature_table_entrezID.csv",sep=";",header = T))
features = read.csv("test_dataset/test_feature.txt",sep="\t",header = T,row.names = 1)
# head(features_ori)
# feat_colums = colnames(features_ori)
# ver = c("Lysosome","Peroxisome","Endosome","Lipid.droplet",'Vesicle')
# ver_col = feat_colums[feat_colums %in% ver]
# if(length(ver_col)>1){ 
#   features = cbind(features_ori[,feat_colums[!feat_colums %in% ver]],apply(features_ori[,feat_colums[feat_colums %in% ver]], 1, sum))
# }else if(length(ver_col)==1){
#   features = cbind(features_ori[,feat_colums[!feat_colums %in% ver]],features_ori[,ver_col])
# }
# colnames(features)[ncol(features)]= 'Vesicle' 
# feat_colums  = colnames(features)

# Vesicle = rowSums(features_ori[,c("Lysosome","Peroxisome")])
# feat_colums = colnames(features_ori)
# features = features_ori[,!feat_colums %in% c("Lysosome","Peroxisome")]
# features = cbind(features,Vesicle)
# features = unique(features[,!feat_colums %in% 'Entrez'])
# features = features[complete.cases(features),]
# head(features)



# dup = features[duplicated(features[,1]),1]
# for(i in dup){
#   features = features[features[,'ACC']!=i,]
# }
# rownames(features) = features[,1]
# features = features[,-1]
# colnames(features) = gsub('.',' ',colnames(features),fix=T)
# dim(features)


# intersection of proteins which have feature data and proteins in the tp-phy graph
prot = row.names(features)
A_prot = row.names(A_ori)
common = prot[prot %in% A_prot]
uncommon = A_prot[!A_prot %in% prot]
features = features[common,]
if(length(uncommon)>0){
  f = matrix( rep(0,length(uncommon)*ncol(features)),nrow=length(uncommon),ncol = ncol(features),dimnames = list(uncommon,colnames(features)))
  features = rbind(features,f)
}
features = features[,colSums(features) != 0  ]


A = as.matrix(A_ori)
# wA = as.matrix(wA_ori)



# A = A_ori[common,common]
# wA = wA_ori[common,common]
# dim(A)
# rowSums = Matrix::rowSums
# colSums = Matrix::colSums

# scls = c("Cytoskeleton", "Cytosol","Endoplasmic.Reticulum","Golgi.apparatus"
         # ,"Mitochondrion","Nucleus","Plasma.membrane","Extracellular.space", "Cell.junction")
# SCL data
# exp =  '../SCL_data/CMPT_hpa_scl_acc_2018.txt'
# exp = 'all_ago2_related_prots.txt_scl.txt'
# data_exp = read.csv(exp,';',header = T)
# ago2_scls = unique(data_exp[,2])
# data_exp = data_exp[!grepl("complex",data_exp[,2]),]
# 
# # write.csv(ago2_scls[!grepl("complex",ago2_scls) & !ago2_scls%in%rownames(data_exp_map)],file = 'ago2_scls_3.csv',quote=F,row.names = F)
# 
# exp_map= 'ago_scl_map.csv'
# data_exp_map = read.csv(exp_map,';',header = T)
# # unique(data_exp[,2]) [!unique(data_exp[,2]) %in% row.names(data_exp_map)]
# exp_mapped = unique(merge(data_exp_map ,data_exp, by.y = 'SCL',by.x='scls')[,c('ACC','map')])

exp = 'test_dataset/test_scl.txt'
# write.table(na.omit(exp_mapped),file = exp,qpouote=F,row.names = F,sep= ';')


L_real = loadProteinSCL(exp , A, ';')
# L_real = as.matrix(L_real)
scls = colnames(L_real)
print('Investigated SCLs are ...')
print(scls)
# ver = c("Lysosome","Peroxisome","Endosome","Lipid.droplet",'Vesicle')
# ver_col = scls[scls %in% ver]
# if(length(ver_col)>1){ 
  # L_real = cbind(L_real[,scls[!scls %in% ver]],rowSums(L_real[,scls[scls %in% ver]]))
# }else if(length(ver_col)==1){
  # L_real = cbind(L_real[,scls[!scls %in% ver]],L_real[,ver_col])
# }
# colnames(L_real)[ncol(L_real)]= 'Vesicle' 
# scls = colnames(L_real)
# scls_p = scls[3,4] 

# scls_f = colnames(features)

# L_real[,L_real['Q9UKV8',]==1]


# Nucleoplasm Nucleus P-body Polysome Cytoplasm Cytosol Extracellular exosome,Membrane
# kn_nb = rownames(L_real)
# all_nb[all_nb %in% kn_nb]
# L_real[nb,]

# par(xpd=FALSE)
# pie(colSums(L_real[nb,]))
# box()
# grid()

# scl_off = scls[!scls %in% scl_on]
# L_off = sparseMatrix(i = {}, j = {} , x = 1, dims = c(nrow(L_real), length(scl_off)))
# colnames(L_off) = scl_off
# rownames(L_off) = rownames(L_real)
# L_real = cbind(L_real,L_off)

# L_real = L_real[common,scls]
# features = features[t_prot_name,]
# scls = scls[!scls %in% c('Myofibril','Sarcomere')]
# L_real = L_real[,scls]
# scls


# overall sta
u = rowSums(as.matrix(L_real));
names(u)=rownames(L_real)

unlabeled = names(which(u == 0));
# unlabeled = c("27161",unlabeled) 
labeled = names(which(u > 0));
L= L_real
# scls = colnames(L)


cat("Nb of unlabeled proteins:",length(unlabeled),"\n")
#print(unknowns)
cat("Nb of labeled proteins",length(labeled),"\n")
if(length(unlabeled)>0){
  L[unlabeled,] = -1
}

# Lc = colSums(L_real[labeled,])/(nrow(L_real[labeled,])-colSums(L_real[labeled,]))
Lc = (nrow(L_real[labeled,])-colSums(L_real[labeled,]))/colSums(L_real[labeled,])


L_all = L


niter = 250;
burnin = 50
pcut = 0.5

# nb_mask = sample()
# masked = "Q9UKV8"
# seed = labeled[!labeled %in% masked]
# scls =scls[colSums(L_real[seed,])>0]
# L_all = L_all[,scls]
# L_all[L_all>1] =1
# L_real[L_real>1] =1
features = features[,colnames(features)%in%scls]


L_all_prob = L_real[,scls]
auc_changes = c()
auc = c()
L_all_plf = L_real[1,]
cols = rainbow(length(scls))
names(cols) = scls

adjMatrix = adjMatrix[scls,scls]


ini_record=c()
for(Lname in scls){
# for(Lname in c("Vesicle")){
  tryCatch({
  # Lname = scls[3]
  cat("SCL:",Lname,"\n")
  iteration = 0
  #---------------------
  
  L_all[masked,Lname] = -1

  # U is the proteins has no annotation
  cat("Nb of masked proteins:",length(masked),"\n")
  cat("Nb of seed proteins",length(seed),"\n")
  
  #-----------------------------#
  #  basic
  #-----------------------------#  
  Lsignle_ori  = L_real[,Lname]
  Lsingle = L_all[,Lname];
  Lc_single = Lc[Lname]
  L_seed =L_real[seed,Lname]
  type=Lsingle
  type[seed]="s"
  type[masked]="m"
  type[unlabeled]="u"
  
  unknown = unique(c(unlabeled,masked))
  #unknown_all[[Lname]] = unknown
  cat("Nb of unknown proteins:",length(unknown),"\n")
  
  
  #-----------------------------#     
  #  Initiation 
  #-----------------------------#
  cat("ppi!\n")
  #pdf(paste(Lname,"ppi.pdf"))
  # 1/(1 + exp(-vc))
  # cat(1  )
  L_initialized = try(BMRF_single_imbalance(wA, L = Lsingle, Lname=Lname, burnin = burnin, niter = niter,vis=FALSE,0.5,masked,L_seed,Lsignle_ori), silent=FALSE);
  if(length(L_initialized)<5){
    L_initialized = try(BMRF_single_imbalance_K0(wA, L = Lsingle, Lname=Lname, burnin = burnin, niter = niter,vis=FALSE,0.5,masked,L_seed,Lsignle_ori), silent=FALSE);
  }
  
  L_all[,Lname] = L_initialized[[2]]
  L_all_prob[,Lname]  =  L_initialized[[1]]
  L_all_plf[Lname]  =  L_initialized[[5]]
  # auc is f measure here
  # auc=c(auc,L_initialized[[3]])
  #   
  # L_write = L_real[,Lname]
  # L_write[unlabeled]=-1
  
  # towrite = cbind(L_initialized[[1]],L_initialized[[2]],L_write,type,L_initialized[[3]])
  # colnames(towrite) = c("prob","predicted","real",Lname,"F_score")
  # write.table(towrite,file = paste(c(dirname,"/",iteration,"_",Lname,"_fold_",CV_index,"_tissue_ppi_output.txt"),collapse=""))
}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  }


L_all_initialized = L_all
colnames(L_all_initialized ) = scls
L_all_prob_initialized = L_all_prob
colnames(L_all_prob_initialized ) = scls
# names(cols) = scls
# names(auc) = scls

# plot(NA, xlim=c(0,50),ylim=c(0,1) , main ="Collective MRFs",ylab="F1 score",xlab = "Iterations")
# points(rep(0,length(auc)),auc,pch=16,col = cols)
# legend("topright",scls , pch=16,col=cols,cex=0.6)

# initiation end
ago2_scl[[CV_index]]=L_all[masked,]
# }

# for(ago2 in ago2_scl){
#   write(ago2,file = 'ago2_scl_out.txt',append = T)
# }
#-----------------------------------------------#

#-----------------------------------------------#

L_all_iterative = list()
L_all_iterative[[1]] = L_all_initialized
L_all_plf_iterative = c(0,L_all_plf)
L_all_prob_iterative = list()
L_all_prob_iterative[[1]] = L_all_prob_initialized
# auc_changes_iterative = list()
# auc_changes_iterative[[1]] = auc
print("Initialized!")

# auc_changes = auc
# auc_m = auc
# names(auc_m) = scls
# names(auc_changes) = scls
temp= 1000
coolrate = 0.2

L_all = L_all_initialized
colnames(L_all_initialized ) = scls
L_all_prob_initialized = L_all_prob
colnames(L_all_prob_initialized ) = scls
names(cols) = scls
# names(auc) = scls
nit=50
for(iter in 1:nit){
  cat("iteration:",iter,"\n")
  
  # L_all is initialized here!
  iteration = iter
  
  for(Lname in scls){
    result = NULL
    tryCatch({
    # then = Sys.time()
    cat("SCL:",Lname,"\n")
    Lsingle = L_all[,Lname];
    Lsignle_ori  = L_real[,Lname]
    # Lc_single = Lc[Lname]
    Adjs=names(which(adjMatrix[scls,Lname]==1))
    if(length(Adjs)>1){
      Adjs = Adjs[colSums(L_all[seed,Adjs])>0]
    }else{
      Adjs = Adjs[sum(L_all[seed,Adjs])>0]
    }
    if(Lname == "Cytoplasm") {
      Adjs = c()
    }
    #----------------------------#
    # + Adjacent SCL + features
    #----------------------------#
    #cat("ppi+adj!\n")
    if(Lname %in% colnames(features)){
      cat("ppi+adj+feat!\n")
      # then = Sys.time()
      feature = features[names(Lsingle),Lname]
      names(feature) = names(Lsingle)

      # BMRF_adj_result= try(BMRF_adj(wA=wA,Adjs = Adjs,knowns=seed,unknowns=unknown,L_all=L_all,L = Lsingle ,Lname=Lname,
      #                               para=NULL,sigma=NULL,burnin = burnin, niter = niter,flag=0,vis=FALSE), silent=FALSE);
      # BMRF_adj_feature_result = try(BMRF_adj_feature(wA=wA, L = Lsingle, Adjs,Lname=Lname, burnin = burnin, niter = niter,vis=FALSE,pcut,feature,seed,unknowns=unknown,masked,L_masked =L_real[masked,Lname],Lsignle_ori), silent=FALSE);
      
      
      # BMRF_adj_known_feature_imbalance_tissue_result = try(BMRF_adj_known_feature_imbalance_tissue(wA=A,tA, L = Lsingle, Adjs,Lname, burnin, niter ,vis=FALSE,pcut,feature,masked,L_all,L_seed,Lsignle_ori,Lc_single,knowns = seed, unknowns = c(masked, unlabeled)), silent=FALSE);
      
      # BMRF_adj_feature_imbalance_tissue_result = try(BMRF_adj_feature_imbalance_tissue(wA,tA, L = Lsingle, Adjs,Lname, burnin, niter ,vis=FALSE,pcut,feature,masked,L_all,L_seed,Lsignle_ori,Lc_single,knowns = seed, unknowns = c(masked, unlabeled)), silent=FALSE);
      # result = BMRF_adj_feature_imbalance_tissue_result
      BMRF_adj_feature_imbalance_tissue_result = try(BMRF_adj_feature_imbalance(wA, L = Lsingle, Adjs,Lname, burnin, niter ,vis=FALSE,pcut,feature,masked,L_all,L_seed,Lsignle_ori,knowns = seed, unknowns = c(masked, unlabeled)), silent=FALSE);
      result = BMRF_adj_feature_imbalance_tissue_result
      if(length(result)<5){
        BMRF_adj_feature_imbalance_tissue_result = try(BMRF_adj_feature_imbalance_K0(wA, L = Lsingle, Adjs,Lname, burnin, niter ,vis=FALSE,pcut,feature,masked,L_all,L_seed,Lsignle_ori,knowns = seed, unknowns = c(masked, unlabeled)), silent=FALSE);
        result = BMRF_adj_feature_imbalance_tissue_result
      }
    }else{
      BMRF_adj_imbalance_tissue_result = try(BMRF_adj_imbalance_tissue_2(wA, L = Lsingle, Adjs,Lname, burnin, niter ,vis=FALSE,pcut,masked,L_all,L_seed,Lsignle_ori,knowns = seed, unknowns = c(masked, unlabeled)), silent=FALSE);
      # BMRF_adj_imbalance_tissue_result = try(BMRF_adj_imbalance(wA, L = Lsingle, Adjs,Lname, burnin, niter ,vis=FALSE,pcut,masked,L_all,L_seed,Lsignle_ori,Lc_single,knowns = seed, unknowns = c(masked, unlabeled)), silent=FALSE);
      result = BMRF_adj_imbalance_tissue_result
    }
    
    if(!is.null(result)){
      logc = L_all_plf[Lname]
      logp = result[[5]]
      
      # current = auc_m[Lname]
      # proposed = BMRF_adj_feature_result[[3]]
      #if(!is.na(proposed) & proposed > current) # strict
      # if(sum(logp)> sum(logc)) # strict
      acpt = acceptanceProbability(sum(logc),sum(logp),temp)
      rand = runif(1)
      print(paste(c('sum(logc)','sum(logp)','temp','acpt','runif(1)'),collapse = ' - '))
      print(paste(round(c(sum(logc),sum(logp),temp,acpt,rand),3),collapse = ' - '))
      if( acpt > rand)
        # if(sum(logc) < sum(logp))
      {
        print("accept!")
        L_all_prob[,Lname] = result[[1]]
        L_all[,Lname] = result[[2]]
        # auc_m[Lname] = result[[3]]
        L_all_plf[Lname] = result[[5]]
      }
    }
    # points(iter, auc_m[Lname],pch=16,col = cols[Lname])
    # now = Sys.time()
    # L_write = L_real[,Lname]
    # L_write[unlabeled]=-1
    #towrite = cbind(BMRF_adj_feature_result[[1]],BMRF_adj_feature_result[[2]],L_write,type, BMRF_adj_feature_result[[3]])
    # towrite = cbind( L_all_prob[,Lname], L_all[,Lname],L_write,type, auc_m[Lname])
    
    # colnames(towrite) = c("prob","predicted","real",Lname,"F1")
    # write.table(towrite,file = paste(c(dirname,"/",iteration,"_",Lname,"_fold_",CV_index,"_ppi_adj_output.txt"),collapse=""))
    
    # print(auc_m[Lname])
    # print(now - then)
    
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  }
  # auc_changes = rbind(auc_changes, auc_m)
  # L_all_iterative[iter] = L_all
  print(iter)
  print(L_all[masked,])
  # L_all_prob_iterative[iter] = L_all_prob
  # L_all_plf_iterative = rbind(L_all_plf_iterative,c(iter,L_all_plf))
  
  ago2_scl_2[[CV_index]][[iter]]=L_all[masked,]
  if(temp>1){
    temp=temp*(1-coolrate)
  }
}

df <- data.frame(matrix(unlist(ago2_scl_2[[CV_index]]), nrow=nit, byrow=T))
df = rbind(ago2_scl[[CV_index]],df) 
colnames(df) = scls
write.table(df,file = paste0(CV_index,'_result_cyto.txt'),quote = F)
}
  print(ago2_scl)
  print(ago2_scl_2)


# par(mfrow=c(1,1))
# 
# # plot(NA, xlim=c(0,50),ylim=c(0,1) , main ="Collective MRFs",ylab="F1 score",xlab = "Iterations")
# # points(rep(0,length(auc)),auc,pch=16,col = cols)
# # legend("topright",scls , pch=16,col=cols,cex=0.6)
# 
# # initiation end
# 
# pdf(paste0(dirname,'/energy_f1.pdf'))
# par(mfrow=c(1,1))
# 
# 
# #}
# # change of PLF over interations
# # plf_changes = t(sapply(L_all_plf_iterative, function(x){colSums(log(x))}))
# plf_changes = L_all_plf_iterative[,-1]
# matplot(-plf_changes,type='l',col=cols,ylim=c(min(-plf_changes),max(-plf_changes)),lty=c(1:11),lwd = 2
#         ,ylab = 'Energy(-log(PLF))',xlab="iteration")
# legend("topright",legend = scls,col=cols,lwd = 2,lty=c(1:11),bty = "n",cex=.5)
# abline(v = 2,lty='dotted')
# grid()
# 
# # changes of  F score over iterations
# nbcol = ncol(auc_changes)
# cols = rainbow(nbcol)
# rownames(auc_changes) = c('ini',as.character(1:50))
# # plot(NA, xlim=c(1,nrow(auc_changes)),ylim=c(0,1) , main ="Collective MRFs",ylab="F1",xlab = "Iterations")
# # for(i in 1 : length(scls)){
# #   lines(auc_changes[,i]/100,col = cols[i],lwd =2)
# # }
# # legend("bottomright",legend = scls,col=cols,lwd =2,lty=1,bty = "n")
# matplot(auc_changes,type='l',col=cols,ylim=c(0,1),lty=c(1:11),ylab="F1 score", xlab="iteration",lwd=2)
# abline(v = 2,lty='dotted')
# abline(v = 20, col='grey',lwd=2)
# grid()
# legend("topright",legend = scls,col=cols,lwd = 2,lty=c(1:11),bty = "n",cex=.5)
# legend("bottomright","* The F1 scores are calculated by applying the generic SCLs ground truth data on TSP PPI network",bty = "n",cex=.7)
# dev.off()
