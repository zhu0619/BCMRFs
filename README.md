# BCMRFs
Bayesian collective Markov Random Fields

Advanced biotechnology makes it possible to access a multitude of heterogeneous proteomic, interactomic, genomic, and functional annotation data. One challenge in computational biology is to inte- grate these data to enable automated prediction of the Subcellular Localizations (SCL) of human proteins. For proteins that have mul- tiple biological roles, their correct in silico assignment to different SCL can be considered as an imbalanced multi-label classification problem. In this study, we developed a Bayesian Collective Markov Random Fields (BCMRFs) model for multi-SCL prediction ofhuman proteins. Given a set of unknown proteins and their corresponding protein-protein interaction (PPI) network, the SCLs of each pro- tein can be inferred by the SCLs of its interacting partners. To do so, we integrate PPIs, the adjacency of SCLs and protein features, and perform transductive learning on the re-balanced dataset. Our experimental results show that the spatial adjacency of the SCLs improves multi-SCL prediction, especially for the SCLs with few annotated instances. Our approach outperforms the state-of-art PPI- based and feature-based multi-SCL prediction method for human proteins.

The detail of the method can be found in our publication

Zhu,L. and Ester,M. (2017) Bayesian Collective Markov Random Fields for Subcellular Localization Prediction of Human Proteins. In Proceedings of the 8th ACM International Conference on Bioinformatics, Computational Biology,and Health Informatics - ACM-BCB ’17. ACM Press, Boston, MA, USA, pp. 321–329.  http://dl.acm.org/citation.cfm?doid=3107411.3107412


## 

BCMRFs_module folder includes the source codes of BCMRFs

## BCMRFs for tissue-specific SCL prediction
An example of generic SCL prediction and an example of tissue-specific SCL prediction can be found in R script 'test.R'
