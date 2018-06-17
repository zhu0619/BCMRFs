
# This script show an example of SCL prediction using BCMRFs algorithm
# Author : Lu Zhu
# Email : lzhu@techfak.uni-bielefel.de
#---------------------------------------------------------#
source('bcmrf_test_optimization.R')


#---------------------------------------------------------#
# General SCL prediction using only physical PPI networks
file_test = "test_datasets/phy_ppi_test.gml"
try(optimization(file_test,'test_results/generic'),silent=T)


#---------------------------------------------------------#

# Tissue-specific SCL prediction using tissue-specific PPI networks
file_ts_test = "test_datasets/ts_phy_ppi_test.gml"
try(optimization(file_ts_test,'test_results/ts'),silent=T)
