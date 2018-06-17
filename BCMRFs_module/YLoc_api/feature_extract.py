import re;
import os;
import sys;
'''
    This script is for collecting all the feature-based probablity values from the result text files from Yloc into a table.
    Author: Lu Zhu
    Date: 2016.11.04
    Usage：
    pyhton feature_extract.py file_out.txt
'''
import mygene
mg = mygene.MyGeneInfo()

if len(sys.argv) > 1:
    fout = sys.argv[1]
    cwd = os.getcwd()
    dirlist = os.listdir(cwd)
    all_loc=[];
    all_features = [];
    all_acc = []
    for files in dirlist :
        if 'result_' in files:
            with open(files,'r') as infile:
                features = dict();
                copy = False
                for line in infile:
                    if "Prediction of sequence" in line:
                        #print line;
                        seqID = line.split('|',)
                        #print seqID[1];
                        features['ID'] = seqID[1]; 
                        all_acc.append(seqID[1])
                    if line.strip() == "Location (Probability)":
                        copy = True
                    elif line.strip() == "Detail 2 - Attributes":
                        copy = False
                    elif copy:
                        #print line;
                        loc = re.findall('[a-zA-Z]+ *[a-zA-Z]*',line)[0].rstrip();
                        if loc not in all_loc:
                            all_loc.append(loc);
                            
                        prob = re.findall('[0-9]+.[0-9]+',line)[0].rstrip();
                        features[loc] = prob
            infile.close();
            all_features.append(features.copy());
    
    maps = mg.querymany(list(set(all_acc)),scopes='uniprot',field= 'entrezgene', species='human',returnall=True)['out']

    with open(fout, 'w') as outfile:
        outfile.write('ACC;Entrez;');
        for loc in all_loc:
            outfile.write(loc+';');
        outfile.write('\n');
        # feature values
        for prot in all_features:
            if 'ID' in prot:
                outfile.write(prot['ID']+';');
                entrez=str([ i['entrezgene'] for i in maps if i['query']== prot['ID']][0])
                outfile.write(entrez+';');
                for loc in all_loc:
                    if loc in prot:
                        outfile.write(prot[loc]);
                    outfile.write(';');
                outfile.write('\n');
    outfile.close();
else:
    print "Error! No exporting file name!"
