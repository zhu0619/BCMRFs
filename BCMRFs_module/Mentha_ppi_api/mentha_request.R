#---------------------------------------------------------#  
# This script is used for request protien interaction data 
# from Mentha RESTful API 
# Author : Lu Zhu
# Email : lzhu@techfak.uni-bielefel.de
#---------------------------------------------------------#  



#---------------------------------------------------------#  
# Load dependencies
#---------------------------------------------------------#
pkg_to_require = c("httr","jsonlite")
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
# Input: list of protein (UniProt ID)
# Output: physical protein interaction networks 
# Reliablity score: 0.2 (by default) (Optinal)
#---------------------------------------------------------#





