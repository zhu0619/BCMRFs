#---------------------------------------------------------#  
# This script is used for request protien interaction data 
# from Mentha RESTful API 
# Author : Lu Zhu
# Email : lzhu@techfak.uni-bielefel.de
#---------------------------------------------------------#  



#---------------------------------------------------------#  
# Load dependencies
#---------------------------------------------------------#
pkg_to_require = c("httr","jsonlite","igraph")
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
HOME_mentha = "http://mentha.uniroma2.it:8080/server/"


github_api <- function(prots) {
  url <- modify_url(paste0(HOME_mentha,'getInteractions?org=9096&ids=',prots))
  print(url)
  GET(url)
}

mentha_ppi= function(prot_list_file){
  prots = readChar(prot_list_file,file.info(prot_list_file)$size)
  resp <- github_api(prots)
  phy_ppi_acc <- read.csv(text=content(resp),sep=';',header=F)
  phy_ppi <- graph.data.frame(phy_ppi_acc, directed=FALSE)
  return(phy_ppi)
}

# prot_list_file <- "prots.txt"  # test protein file
# 
# mentha_ppi(prot_list_file)

