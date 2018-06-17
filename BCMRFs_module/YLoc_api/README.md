yloc_customized.py
######################################################################

Usage: python yloc.py <fastafile.fasta> <model_name> <origin> <use of homology(Yes/No)> <output(Advanced/Simple)>

Available models: YLoc-LowRes, YLow-HighRes, YLoc+
Available origins: Animals, Fungi, Plants

Example:
--------

-- With GO
python yloc_costumized.py test.fasta YLoc-HighRes Animals Yes Advanced


-- no GO 
python yloc_costumized.py test.fasta YLoc-HighRes Animals No Advanced


feature_extrac.py
######################################################################
    This script is for collecting all the feature-based values from the result text files from YLoc into a table.

Usage: pyhton feature_extract.py <output>

pyhton feature_extract.py file_out.txt