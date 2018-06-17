import urllib2;
import urllib;
import socket;
import re;
import time;
import os;
from os.path import exists;
import multiprocessing;

# YLoc HTTP client, 2009 Sebastian Briesemeister
# The script uses the http protocol, so make sure the connection
# is not blocked.

# global definition of YLoc webpage
yloc_url = "http://abi.inf.uni-tuebingen.de/Services/YLoc/webloc.cgi";

# Error class
class Error(Exception):
	def __init__(self, value):
		self.value = value;
	def __str__(self):
		return repr(self.value);

# calls YLoc webservice using the given protein sequences, model name, protein origin, and the use of homology
# the function returns a prediction id
def __startPrediction(aasequence_string, model_name, origin, homology):

	# set up a dictionary with the post form variables
	yloc_input = { "plain_sequence" : str(aasequence_string),
		"model" : str(model_name),
		"origin" : str(origin),
		"goterms" : str(homology),
		"page" : "predict"};
		
	
	# call YLoc via http and retrieve webpage
	data = urllib.urlencode(yloc_input);
	req = urllib2.Request(yloc_url, data);
	response = urllib2.urlopen(req);
	webpage = response.read();	
		
	# extract prediction id
	re_res = re.findall("Query ID: .{32}", webpage);
	
	if len(re_res) == 1:
		# return prediction ID
		return re_res[0][10:]
	else:
		# if not found, raise error
		print webpage
		raise Error("ERROR: Did not receive prediction ID from YLoc webservice!");

# returns page retrieved with given post input
def __getResultPage(yloc_input):
	
	# call YLoc via http and retrieve webpage
	data = urllib.urlencode(yloc_input);
	req = urllib2.Request(yloc_url, data);
	response = urllib2.urlopen(req);
	webpage = response.read();		
		
	re_res = re.findall("<h2>Prediction in progress</h2>", webpage);
	
	if len(re_res) == 1:
		# return empty string
		return ""
	else:
		return webpage;
	
# retrieves detailed result page from YLoc server for the given protein nr in the prediction
# it returns a tupel of (dictionary of location proabilities, Swiss-Prot AC of most similar protein in SP 42.0, dictionary of attributes with discrimination score)
def __parseDetailedResultPage(id, nr):
	#get webpage
	yloc_input = { "id" : str(id), "detailed" : str(nr)};
	webpage = __getResultPage(yloc_input);
	
	
	# get probability distribution
	probabilities = {};
	re_res = re.findall("<tr><td class=amatrix.{1,80}>.{1,30}</td><td class=amatrix.{1,90}>.{3,10}</td></tr>", webpage);
	for elem in re_res:
		re_res2 = re.findall(">.{2,20}</td><", elem);
		locs = len(re_res2)/2;
		for i in range(locs):
			probabilities[re_res2[i*2][1:-6]] = re_res2[i*2+1][1:-7];

	
	#get most similar protein
	most_similar ="";	
	re_res = re.findall("The most similar protein from Swiss-Prot 42\.0 to the query is <a href='http://www\.expasy\.org/uniprot/.{6}", webpage);
	if len(re_res) == 2:
		most_similar = re_res[0][-6:];

	#get attributes with discrimination score
	attributes = {};
	re_res = re.findall("<tr><td class='amatrix.'.{0,100}>.{1,100}</td>\n<td class='amatrix.'.{0,100}>\n<a class='helptext'.{1,130} href=#G>.{1,130}\n", webpage);
	re_res2 = re.findall("</span>\n</a></td>\n<td class='amatrix.'.{0,100}>\n.{2,8}</td>", webpage);
	
	for i in range(len(re_res)):
		quantity = re.findall(">.{1,30}</td>\n",re_res[i])[0][1:-6];
		attr = re.findall("href=#G>.{1,130}\n",re_res[i])[0][8:-1];
		score = re.findall(">\n.{2,8}</td>$",re_res2[i])[0][2:-5];
		attributes[quantity+" "+attr] = score;
		
	return (probabilities, most_similar, attributes);
	
	
# calls YLoc webservice using a prediction id
# the function will return a list of predictions
# simple mode:
#	each prediction is a tupel of (sequence name, predicted location,  probability, confidence score, reasoning)
# advanced mode:
#	each prediction is a tupel of (sequence name, dictionary of locations with probability, confidence score, SwissProt AC of most similar protein, reasoning,  dictionary of features with feature name and discrimination score)
def __getResults(id, advanced=False):
	
	# output list
	output = [];
	
	yloc_input = { "id" : str(id)};
	
	webpage = __getResultPage(yloc_input);
	while len(webpage)==0:
		time.sleep(2);
		webpage = __getResultPage(yloc_input);
	
	# extract number of predictions and basic informations
	re_res = re.findall("id="+str(id)+"&detailed=.{1,2}>\n.{1,120}\n</a></td>", webpage);
	
	if len(re_res) % 4 != 0:
			# if results do not match
			raise Error("ERROR: Could not parese YLoc output!");
	
	no_predictions = len(re_res) / 4;
	sequence_names = [];
	predicted_locations = [];
	probabilities = [];
	confidence_scores = [];
	
	for i in range(len(re_res)):
		if i%4 == 2:
			#exctract probability
			re_res[i] = re.findall("\n.{1,10} %", re_res[i])[0];
			re_res[i] = re_res[i][1:-1];
			probabilities.append(re_res[i]);
		elif i%4 == 3:
			# extract confidence score
			re_res[i] = re.findall("\(.{1,10}\)", re_res[i])[0];
			re_res[i] = re_res[i][1:-1];
			confidence_scores.append(re_res[i]);
		else:
			# extract value
			re_res[i] = re.findall("\n.{1,120}\n", re_res[i])[0];
			re_res[i] = re_res[i][1:-1];
			
			if i%4 == 0:
				sequence_names.append(re_res[i]);
			elif i%4 == 1:
				predicted_locations.append(re_res[i]);	
	
	#get reasoning
	re_res = re.findall("<h3>Why\?</h3><p>.{300,1500}</b></p>",webpage);
				
	if len(re_res) != no_predictions:
		# if results do not match
		raise Error("ERROR: Could not parese YLoc output!");

	reasonings = [];
	for i in range(len(re_res)):
		#delete html tags from reasoning
		reasonings.append( re_res[i][16:-8] );
		reasonings[i] = re.sub("<.{1,20}>","", reasonings[i]);
	
	if advanced:
		for i in range(no_predictions):
			(probability_distr, most_similar, attributes) = __parseDetailedResultPage(id, i);
			output.append( (sequence_names[i], predicted_locations[i], probabilities[i], probability_distr, confidence_scores[i], most_similar, reasonings[i], attributes)  );	
	else:
		for i in range(no_predictions):
			output.append( (sequence_names[i], predicted_locations[i], probabilities[i], confidence_scores[i], reasonings[i]) );
	
	return output;
		
		



def predict2(folder,sequence_list, model_name, origin, homology, advanced=False):
    # split sequences in packs of 20 since the webserver allows at most 20 protein sequences
    i = 0;
    #results = [];
    result = [];
    while i < len(sequence_list):
        fasta_string = ">"+str(sequence_list[i]['id'])+"\n"+sequence_list[i]['sequence']+"\n";
        #check file
        fname= str(sequence_list[i]['id']);
        print fname;
        fname = fname.replace("|", "--");
        fname = folder + "/result_"+fname+".txt";
        i += 1;       
        if check_file(fname)== False:
            #print(fasta_string);
            id = __startPrediction(fasta_string, model_name, origin, homology);
            result = __getResults(id, advanced)
            if len(result) > 0:
            #results += __getResults(id, advanced)
            #print(result);
                result = result[0]
                try:
                    writeResult2(result,folder);
                    print "Done!"
                except Error, e:
                    print "An ERROR occured: "+str(e.value);
                    print "Use no arguments to get help!"
            else:
                fout = open(folder+"/No_result_list.txt","a");
                fout.write(fasta_string+"\n");
                fout.close();
                print "No result!";
        else:
            print  fname+" exists!"
 


def predict3(folder,sequence_list, seq_range ,model_name, origin, homology, advanced=False):
    # split sequences in packs of 20 since the webserver allows at most 20 protein sequences
    i = 0;
    #results = [];
    result = [];
    for i in seq_range:
        fasta_string = ">"+str(sequence_list[i]['id'])+"\n"+sequence_list[i]['sequence']+"\n";
        #check file
        fname= str(sequence_list[i]['id']);
        print fname;
        fname = fname.replace("|", "--");
        fname = folder + "/result_"+fname+".txt";
        i += 1;       
        if check_file(fname)== False:
            #print(fasta_string);
            id = __startPrediction(fasta_string, model_name, origin, homology);
            result = __getResults(id, advanced)
            if len(result) > 0:
            #results += __getResults(id, advanced)
            #print(result);
                result = result[0]
                try:
                    writeResult2(result,folder);
                    print "Done!"
                except Error, e:
                    print "An ERROR occured: "+str(e.value);
                    print "Use no arguments to get help!"
            else:
                fout = open(folder+"/No_result_list.txt","a");
                fout.write(fasta_string+"\n");
                fout.close();
                print "No result!";
        else:
            print  fname+" exists!"


def predict3_bis(folder,sequence_list, model_name, origin, homology, advanced=False):
    # split sequences in packs of 20 since the webserver allows at most 20 protein sequences
 	nb_thread = 4
	seq_range_1 = range(len(sequence_list)/nb_thread)	
        p_1 = multiprocessing.Process(target=predict3, args=(folder,sequence_list,seq_range_1 ,model_name, origin, homology, False,))
        seq_range_2 = range(seq_range_1[-1]+1,len(sequence_list)*2/nb_thread)
        p_2 = multiprocessing.Process(target=predict3, args=(folder,sequence_list,seq_range_2,model_name, origin, homology, False,))
        seq_range_3 = range(seq_range_2[-1]+1,len(sequence_list)*3/nb_thread)
        p_3 = multiprocessing.Process(target=predict3, args=(folder,sequence_list, seq_range_3,model_name, origin, homology, False,))
        seq_range_4 = range(seq_range_3[-1]+1,len(sequence_list))
        p_4 = multiprocessing.Process(target=predict3, args=(folder,sequence_list,seq_range_4 ,model_name, origin, homology, False,))
        p_1.start()
        # p_2.start()
        # p_3.start()
        # p_4.start()
        p_1.join()
        # p_2.join()
        # p_3.join()
        # p_4.join()



# def predict3_b(i,folder,sequence_list, model_name, origin, homology, advanced=False):
# 	fasta_string = ">"+str(sequence_list[i]['id'])+"\n"+sequence_list[i]['sequence']+"\n";
# 	#check file
# 	fname= str(sequence_list[i]['id']);
# 	print fname;
# 	fname = fname.replace("|", "--");
# 	fname = folder + "/result_"+fname+".txt";
# 	i += 1;       
# 	if check_file(fname)== False:
# 	    #print(fasta_string);
# 	    id = __startPrediction(fasta_string, model_name, origin, homology);
# 	    result = __getResults(id, advanced)
# 	    if len(result) > 0:
# 	    #results += __getResults(id, advanced)
# 	    #print(result);
# 	        result = result[0]
# 	        try:
# 	            writeResult2(result,folder);
# 	            print "Done!"
# 	        except Error, e:
# 	            print "An ERROR occured: "+str(e.value);
# 	            print "Use no arguments to get help!"
# 	    else:
# 	        fout = open(folder+"/No_result_list.txt","a");
# 	        fout.write(fasta_string+"\n");
# 	        fout.close();
# 	        print "No result!";
# 	else:
# 	    print  fname+" exists!"







def writeResult2(result_element,folder):
    #print result_element;
    # check simple or advanced output
    if len(result_element)==5:
        #simple output
        (name, location, probability, confidence_score, reasoning) = result_element;
        # open file
        fname = name.replace("|", "--");
        fname = folder+"/result_"+fname+".txt";
        fout = open(fname,"w");
        fout.write("Sequence:"+name+"\n");
        fout.write("Location:"+location+" ("+str(probability)+"%)\n");
        fout.write( "Confidence score:"+str(confidence_score)+"\n");
        fout.write( "Reasoning:"+ reasoning + "\n");
        fout.close();
    else:
        #advanced outout
        (name, location, predicted_prob, probabilities, confidence_score, most_similar, reasoning, attributes) = result_element;
        # open file
        fname = name.replace("|", "--");
#        fname = "result_"+fname+".txt";
        fname = folder+ "/result_"+fname+".txt";
        fout = open(fname,"w");
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
        fout.close();


def check_file(fname):
    cwd = os.getcwd();
    return exists(cwd+"\\"+fname)
