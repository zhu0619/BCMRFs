import re;
import os;
import errno;
import sys;
import YLocHTTPclient_lu;
from YLocHTTPclient_lu import Error;


# gets an open file stream as input and tests whether this file is in fasta format or not
# returns a list of dictionaries with the id and the sequence (keys="id" "sequence")
def __parse_fasta_file(data):
    try:
        file = data
        tmp = file.tell()
    except AttributeError:
        file = open(data, 'r')
    id=""
    sequence=""
    proteins=[]
    line = file.readline()
    if line[0] != ">": raise Error("wrong format in fasta file %s!" % data)
    id = re.sub("^>","",line)
    id = re.findall("^[\S]+", id)
    if len(id) == 0:
        id = "Seq 1"
    else:
        id = id[0]
    i=0
    proteins.append({'id':id,'sequence':""})
    while 1:
        line = file.readline()
        if i==0 and not line and not sequence: raise Error("wrong format in fasta file %s!" % data)
        if i>0 and not sequence and not line: raise Error("wrong format in fasta file %s!" % data)
        if not line: break
        if line[0] == ">":
            if sequence == "": raise Error("wrong format in fasta file %s!" % data)
            proteins[i]['sequence']=sequence
            sequence = ""
            i=i+1
            id = re.sub("^>","",line)
            id = re.findall("^[\S]+", id)
            if len(id) == 0:
                id = "Seq %s" %(i+1)
            else:
                id = id[0]
            proteins.append({'id':id,'sequence':""})
        else:
            if (re.findall("[A-Za-z]+",line)):
                sequence += re.findall("[A-Za-z]+",line)[0]
    proteins[i]['sequence']=sequence
    return proteins
	
def printResult(result_element):
	# check simple or advanced output
	if len(result_element)==5:
		#simple output
		(name, location, probability, confidence_score, reasoning) = result_element;
		print "Prediction of sequence \""+name+"\":";
		print "Location: "+location+" ("+str(probability)+"%)";
		print "Confidence score: "+str(confidence_score);
		print "Reasoning:";
		print reasoning;
		print "";
	else:
		#advanced outout
		(name, location, predicted_prob, probabilities, confidence_score, most_similar, reasoning, attributes) = result_element;
		print "Prediction of sequence \""+name+"\":";
		print "Location: "+location+" ("+str(predicted_prob)+"%)";
		print "Confidence score: "+str(confidence_score);
		print "Reasoning:";
		print reasoning;
		print "Detail 1 - Probability Distribution";
		print "\tLocation (Probability)";
		for elem in probabilities.keys():
			print "\t"+elem+" ("+str(probabilities[elem])+"%)";
		print "Detail 2 - Attributes";
		print "\tDiscrimination score\t Attribute";
		for elem in attributes.keys():
			print "\t"+str(attributes[elem])+" \t"+elem;
		print "";

def writeResult(result_element):
    #print result_element;
    # check simple or advanced output
    if len(result_element)==5:
        #simple output
        (name, location, probability, confidence_score, reasoning) = result_element;
        # open file
        fout = file("result_"+name+".txt","a+")
        fout.write("Sequence:"+name+"\n");
        fout.write("Location:"+location+" ("+str(probability)+"%)\n");
        fout.write( "Confidence score:"+str(confidence_score)+"\n");
        fout.write( "Reasoning:"+ reasoning + "\n");
        fout.close()
    else:
    #advanced outout
        (name, location, predicted_prob, probabilities, confidence_score, most_similar, reasoning, attributes) = result_element;
        # open file
        fout = file("result_"+name+".txt","a+")
        fout.write("Sequence:"+name+"\n");
        fout.write( "Location: "+location+" ("+str(predicted_prob)+"%)\n");
        fout.write( "Confidence score:"+str(confidence_score)+"\n");
        fout.write( "Reasoning:"+ reasoning + "\n");
        fout.write( "Detail 1 - Probability Distribution\n");
        fout.write( "\tLocation (Probability)\n");
        for elem in probabilities.keys():
            fout.write( "\t"+elem+" ("+str(probabilities[elem])+"%)\n");
        fout.write( "Detail 2 - Attributes\n");
        fout.write( "\tDiscrimination score\t Attribute\n");
        for elem in attributes.keys():
            fout.write( "\t"+str(attributes[elem])+" \t"+elem+ "\n");
        fout.write("");
        fout.close()

def make_sure_path_exists(path):
        try:
                os.makedirs(path)
        except OSError as exception:
                if exception.errno != errno.EEXIST:
                        raise
# MAIN
try:
    # check number of arguments
	if len(sys.argv) ==6:
		#check arguments
		if not sys.argv[2] in ["YLoc-LowRes","YLoc-HighRes","YLoc+"]:
			raise Error("Model with name "+sys.argv[2]+" is not available!");
		if not sys.argv[3] in ["Animals", "Fungi", "Plants"]:
			raise Error("Origin "+sys.argv[3]+" is not available!");
		if not sys.argv[4] in ["Yes", "No"]:
			raise Error("Use of GO terms should be indicated with either 'Yes' or 'No'");
		if not sys.argv[5] in ["Advanced", "Simple"]:
			raise Error("Output options are either 'Advanced' or 'Simple'");
		else:
			advanced = False;
			if sys.argv[5] == "Advanced":
				advanced = True;
		# parse FASTA
		sequence_list = __parse_fasta_file(sys.argv[1]);
                folder = sys.argv[1]+'_folder'
                make_sure_path_exists(folder)
                print folder;
                # call YLoc
		YLocHTTPclient_lu.predict3_bis(folder,sequence_list, sys.argv[2], sys.argv[3], sys.argv[4], advanced);        
		# print results
		print "YLoc - Prediction Output ";
		print "-------------------------";
            #for elem in result:
            #writeResult(elem);
            #for elem in result:
            #writeResult(elem);
	else:
		# print help / manual output
		print "YLoc - Interpretable Subcellular Localization Prediction";
		print "--------------------------------------------------------";
		print "Usage: python yloc.py <fastafile.fasta> <model_name> <origin> <use of GO terms(Yes/No)> <output(Advanced/Simple)>";
		print " ";
		print "Available models: YLoc-LowRes, YLow-HighRes, YLoc+";
		print "Available origins: Animals, Fungi, Plants";
		print "";
except Error, e:
	print "An ERROR occured: "+str(e.value);
	print "Use no arguments to get help!"
