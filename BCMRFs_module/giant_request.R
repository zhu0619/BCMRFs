#---------------------------------------------------------#  
# This script is used for request protien interaction data 
# from Mentha RESTful API 
# Author : Lu Zhu
# Email : lzhu@techfak.uni-bielefel.de
#---------------------------------------------------------#  



#---------------------------------------------------------#  
# Load dependencies
#---------------------------------------------------------#
pkg_to_require = c("httr","jsonlite","igraph","mygene")
pkg_to_install = setdiff(pkg_to_require, rownames(installed.packages()))
# install packages
if (length(pkg_to_install) > 0) {
  for(pkg_name in pkg_to_install){
    print(paste0('Installing package ',pkg_name))
    install.packages(pkg_name)
  }
}

# require packages
for(pkg_name in pkg_to_require){
  print(paste0('requring package ',pkg_name))
  require(pkg_name,character.only = TRUE)
}

#---------------------------------------------------------#  

#---------------------------------------------------------#

HOME_giant = "http://giant-api.princeton.edu/networks/"
giant_api <- function(tissue,prots) {
  url <- modify_url(paste0(HOME_giant,tissue,'/',prots))
  print(url)
  return(fromJSON(url))
}

symbol_acc = function(net){
  ids <- unique(c(net[,"source"],net[,"target"]))
  result = queryMany(ids,scopes = 'symbol',fields='uniprot',species="human")
  map = cbind(attributes(result)$listData$query,attributes(result)$listData$uniprot.Swiss.Prot)
  colnames(map) = c('Symbol','ACC')
  rownames(map) = map[,'Symbol']
  col_a <- lapply(map[net[,1],"ACC"], function(x) ifelse(is.null(x), NA, x))
  col_b <- lapply(map[net[,2],"ACC"], function(x) ifelse(is.null(x), NA, x))
  mapped_net =cbind(unlist(col_a), unlist(col_b))
  mapped_net= cbind(mapped_net,net[,'weight'])
  rownames(mapped_net)=NULL
  colnames(mapped_net) = c("source","target","weight")
  return(data.frame(na.omit(mapped_net)))
}


ts_ppi= function(prot_list_file,tissue){
  prots <- levels(unlist(read.table(prot_list_file,header = F,sep="\n"),use.names = FALSE))
  resp <- giant_api(tissue,paste(prots,collapse = "+"))
  mapped_net = symbol_acc(resp$edges)
  phy_ppi_acc <- graph.edgelist(as.matrix(mapped_net[,c("source","target")]))
  E(phy_ppi_acc )$weight=as.numeric(mapped_net[,"weight"])
  return(phy_ppi_acc)
}

prot_list_file <- "gene_symbol.txt"  # test protein file

tissue = "brain"

ts_ppi_graph = ts_ppi(prot_list_file,tissue)
mentha_ppi(prot_list_file)

