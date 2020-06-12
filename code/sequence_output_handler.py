import os    
import tempfile
#required for webgraphics/apache
os.environ['MPLCONFIGDIR'] = tempfile.mkdtemp()
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys
import numpy as np
import math
import re
import shlex

def build_single_graphs(timestamp, seq_pred_file):
    #feature names for the dev graph
    feat_names = ["KmR","DmE","len","A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y","KpR","DpE","PmN","PpN","aro","pI","mem","chr","fld","dis","ent","bet"]
    
    #line cleaning function and then split on comma
    def arrayLine(line):
        line = line.replace("\n","")
        line = line.replace("\r","")
        line = line.replace(" ","")
        return line.split(",")
    
    #open the sequence prediction file and extract features from jim's perl output
    with open(seq_pred_file,"r") as seq_file:
        for line in seq_file:
            if line.startswith("SEQUENCE DEVIATIONS"):
                clean_line = arrayLine(line)
                deviations = clean_line[2:]
            if line.startswith("SEQUENCE PREDICTIONS"):
                clean_line = arrayLine(line)
                scaledSol = clean_line[3]
                popSol = clean_line[4]
            if line.startswith("SEQUENCE PROFILE charge,>"):
                clean_line = arrayLine(line)
                charge_window = clean_line[2:]
            if line.startswith("SEQUENCE PROFILE Uversky"):
                clean_line = arrayLine(line)
                fold_window = clean_line[2:]    
    
    #MAKE SOLUBILITY BAR GRAPH
    #get the solubility values extracted above
    sols = np.array([float(popSol),float(scaledSol)])
    sol_names = np.array(["PopAvrSol","QuerySol"])
    
    plt.figure(figsize=(3, 5))
    plt.bar("PopAvrSol",float(popSol),color="#454654")
    plt.bar("QuerySol",float(scaledSol),color="#8587a1")
    plt.box(on=None)
    plt.yticks(fontsize=15)
    plt.xticks([0,1], sol_names, rotation=45,ha="center", fontsize=15)
    plt.ylim(0,1)
    plt.title("Solubility",fontsize=20)
    plt.ylabel("Calculated value",fontsize=20)
    plt.axhline(0,c="black")
    plt.tight_layout()
    plt.tick_params(top=False, bottom=False, left=False, right=False, labelleft=True, labelbottom=True)
    plt.savefig("/tmp/run-"+timestamp+"/Solub-bargraph.png")
    plt.savefig("/tmp/run-"+timestamp+"/Solub-bargraph.svg")
    
    #MAKE DEVIATIONS GRAPH
    deviations = np.array([float(i) for i in deviations])
    #extract devations above and below 0
    devup = np.array([i >= 0 for i in deviations])
    devdown = np.array([i < 0 for i in deviations])
    
    len_graph = np.arange(len(feat_names))
    
    plt.figure(figsize=(20, 3))
    plt.bar(len_graph[devup],deviations[devup],color="#ffe80e")
    plt.bar(len_graph[devdown],deviations[devdown],color="#a6c600")
    plt.box(on=None)
    plt.axhline(0,c="black")
    plt.xticks(range(0, len(feat_names), 1),rotation=45)
    plt.xlim(-1,35)
    plt.xticks(len_graph, feat_names,fontsize=15,ha="center")
    plt.yticks(fontsize=15)
    plt.title("Deviations from population average",fontsize=18)
    plt.tight_layout()
    plt.tick_params(top=False, bottom=False, left=False, right=False, labelleft=True, labelbottom=True)
    plt.savefig("/tmp/run-"+timestamp+"/Solub-devgraph.svg")
    plt.savefig("/tmp/run-"+timestamp+"/Solub-devgraph.png")
    
    #MAKE CHARGE GRAPH
    charge_window = np.array([float(i) for i in charge_window])
    chargeup = np.array([i >= 0 for i in charge_window])
    chargedown = np.array([i < 0 for i in charge_window])
    
    len_graph = np.arange(len(charge_window))
    len_feat_names = np.arange(1,len(charge_window)+1)
    
    plt.figure(figsize=(20, 3))
    plt.bar(len_graph[chargeup],charge_window[chargeup],color="#3892d0")
    plt.bar(len_graph[chargedown],charge_window[chargedown],color="#cc0000")
    plt.box(on=None)
    plt.axhline(0,c="black")
    plt.xticks(range(0, len(len_feat_names), 1),rotation=45)
    n_ten = int(len(charge_window)/10)
    n_int = math.ceil(n_ten/10)*10
    plt.xticks(np.arange(0, len(charge_window)+1, n_int),fontsize=15)
    plt.yticks(fontsize=15)
    plt.xlim(-1,len(charge_window))
    plt.title("Windowed charge score per amino acid",fontsize=18)
    plt.tight_layout()
    plt.tick_params(top=False, bottom=False, left=False, right=False, labelleft=True, labelbottom=True)
    plt.savefig("/tmp/run-"+timestamp+"/Solub-chargegraph.svg")
    plt.savefig("/tmp/run-"+timestamp+"/Solub-chargegraph.png")
    
    #MAKE FOLD GRAPH
    
    fold_window = np.array([float(i) for i in fold_window])
    foldup = np.array([i >= 0 for i in fold_window])
    folddown = np.array([i < 0 for i in fold_window])
    
    len_graph = np.arange(len(fold_window))
    len_feat_names = np.arange(1,len(fold_window)+1)
    
    plt.figure(figsize=(20, 3))
    plt.bar(len_graph[foldup],fold_window[foldup],color="#ffa838")
    plt.bar(len_graph[folddown],fold_window[folddown],color="#484cd7")
    plt.box(on=None)
    plt.axhline(0,c="black")
    plt.xticks(range(0, len(len_feat_names), 1),rotation=45)
    n_ten = int(len(fold_window)/10)
    n_int = math.ceil(n_ten/10)*10
    plt.xticks(np.arange(0, len(fold_window)+1, n_int),fontsize=15)
    plt.yticks(fontsize=15)
    plt.xlim(-1,len(fold_window))
    plt.title("Windowed fold propensity per amino acid",fontsize=18)
    plt.tight_layout()
    plt.tick_params(top=False, bottom=False, left=False, right=False, labelleft=True, labelbottom=True)
    plt.savefig("/tmp/run-"+timestamp+"/Solub-foldgraph.svg")
    plt.savefig("/tmp/run-"+timestamp+"/Solub-foldgraph.png")
    
    # save files
    with open("/tmp/run-"+timestamp+"/predicted_deviations.csv","w") as deviationsCSV:
        for index, deviationValue in enumerate(deviations):
            deviationsCSV.write(feat_names[index])
            deviationsCSV.write(",")
            deviationsCSV.write(str(deviationValue))
            deviationsCSV.write("\n")
    
    with open("/tmp/run-"+timestamp+"/proteinsol_prediction.csv","w") as solubilityCSV:
        for index, solubilityValue in enumerate(sols):
            solubilityCSV.write(sol_names[index])
            solubilityCSV.write(",")
            solubilityCSV.write(str(solubilityValue))
            solubilityCSV.write("\n")
    
    with open("/tmp/run-"+timestamp+"/charge_window.csv","w") as chargewindowCSV:
        for index, chargeValue in enumerate(charge_window):
            chargewindowCSV.write(str(index))
            chargewindowCSV.write(",")
            chargewindowCSV.write(str(chargeValue))
            chargewindowCSV.write("\n")
    
    with open("/tmp/run-"+timestamp+"/fold_window.csv","w") as foldwindowCSV:
        for index, foldValue in enumerate(fold_window):
            foldwindowCSV.write(str(index))
            foldwindowCSV.write(",")
            foldwindowCSV.write(str(foldValue))
            foldwindowCSV.write("\n")

def make_single_results(timestamp):
    #location of the results
    resultsdir = "/var/www/html/results/run-"+timestamp

    #get scaled sol and pI
    with open ("/tmp/run-"+timestamp+"/all_predictions.csv") as pred_file:
        for line in pred_file:
            if line.startswith(">"):
                array_line = line.split(",")
                idname = array_line[0]
                scaledSol = array_line[2]
                pI = array_line[4]

    #get TM present
    TMpresent = False
    with open ("/tmp/run-"+timestamp+"/blah.txt") as blah_file:
        for line in blah_file:
            if line.startswith("TM REGION PREDICTED"):
                TMpresent = True
                line = " ".join(line.split())
                line_array = line.split(" ")
                TMval = line_array[10]

    html_header_text = """
    <!DOCTYPE html>
    <html lang='en'>
    
    <head>
      <title>Protein-sol solubility prediction</title>
      <meta charset='utf-8'>
      <meta name='viewport' content='width=device-width, user-scalable=no, minimum-scale=1.0, maximum-scale=1.0'>
      <link rel="icon" type="image/png" href="/resources/favicon-32x32.png" sizes="32x32" />
      <link rel="icon" type="image/png" href="/resources/favicon-16x16.png" sizes="16x16" />
      <link rel='stylesheet' href='/css/stylesheet.css?ver=1'>
      <link rel='stylesheet' href='/css/sequence.css?ver=1'>
    </head>
    
    <body>
    <div id='headerbar'> 
            <a class='headerlogo' href='/'>protein-sol</a>
            <a class='manunilogo' href='http://www.manchester.ac.uk'>UoM</a>
            <a class='bioprologo' href='http://www.biopronetuk.org'>BioProNET</a>
    </div>
    <?php include '/var/www/templates/navigation_links.php';?>
    <div id='content'>
    <div class='solub-grid'>
    <div class='protein'>
    
    <h5 style="display:inline">Protein:</h5> <p style="display:inline">"""+idname+"""</p><br>
    <h5 style="display:inline">Predicted scaled solubility:</h5> <p style="display:inline">"""+scaledSol+"""</p><br>
    <h5 style="display:inline">pI:</h5> <p style="display:inline">"""+pI+"""</p><br>
    <br>
    <p style="display:inline">Job ID:</p> <p style="display:inline"> """+timestamp+"""</p>
    <p style="font-size:0.8em">These results are shareable via the unique url for 7 days from creation date </p>
    
    """

    if TMpresent:
        TM_text = """
	<br>
	<p class='warning'>POSSIBLE TM REGION PREDICTED</p>
	<p class='warning'>Kyle-Doolitle hydropathy value of """+TMval+""" compared to threshold of 1.6</p>
	<p class='warning'>Protein-sol solubility prediction invalid for membrane proteins</p>
	"""
        html_header_text = html_header_text + TM_text

    html_footer_text = """
    <br>
    </div>
    <div class='container'>
    <div class='row'>
    <div class='cell-left'>
    <a href="Solub-bargraph.png"><img id='solub-bar-graph' src='Solub-bargraph.svg'/></a>
    </div>
    <div class='cell-right'>
    <div class='prediction-header'>
    <p>The scaled solubility value (QuerySol) is the predicted solubility. The population average for the experimental dataset (PopAvrSol) is 0.45, and therefore any scaled solubility value greater than 0.45 is predicted to have a higher solubility than the average soluble <em>E.coli</em> protein from the experimental solubility dataset <a href='https://doi.org/10.1073/pnas.0811922106'> Niwa <em>et al</em> 2009</a>, and any protein with a lower scaled solubility value is predicted to be less soluble.</p><br>
    <p>The protein-sol sequence algorithm calculated 35 sequence features. This includes the composition of the standard 20 amino acids and sequence length (<code>len</code>), as well as the following features which are calculated over a sliding 21 amino acid window. </p>
    <p>There are 7 amino acid composite scores:<br><code>KmR</code> = K minus R, <code>DmE</code>=DminusE, <code>KpR</code> = K plus R, <code>DpE</code> = D+E, <code>PmN</code> = K+R-D-E, <code>PpN</code> = K+R+D+E, <code>aro</code> = F + W + Y</code><p>
    <p>We then calculate a further 7 sequence features:<br>
    <code>fld</code> = folding propensity <a href='https://doi.org/10.1002/1097-0134(20001115)41:3%3C415::AID-PROT130%3E3.0.CO;2-7'>Uversky <em>et al</em> 2000<a>, <code>dis</code> = disorder propensity <a href='https://doi.org/10.1093/nar/gkg519'>Linding <em>et al</em> 2003</a>, <code>bet</code> = beta strand propensities <a href='https://doi.org/10.1016/j.bbrc.2006.01.159'>Costantini <em>et al</em> 2006</a>, <code>mem</code> = Kyte-Doolittle hydropathy <a href='https://doi.org/10.1016/0022-2836(82)90515-0'>Kyte and Doolittle 1982</a>, <code>pI</code>, <code>ent</code> = sequence entropy, <code>abs</code> = absolute charge at pH 7.
    <p>Further information is available in the <a href='https://doi.org/10.1093/bioinformatics/btx345'>paper</a>.
    
    <br><br>
    </div>
    <br>
    <br>
    </div>
    </div>
    </div>
    <div class='container'>
    <div class='row'>
    <a href="Solub-devgraph.png"><img style='width:100%;' src='Solub-devgraph.svg'></a>
    <a href="Solub-chargegraph.png"><img style='width:100%;' src='Solub-chargegraph.svg' ></a>
    <a href="Solub-foldgraph.png"><img style='width:100%;' src='Solub-foldgraph.svg'></a>
    </div>
    </div
    <div class='cell-results'>
    <div style='text-align:center;'>
    <button style='margin:auto;display:block;' id='predictionButton'>Download Prediction</button>
    </div>
    <br>
        <br>
        <h5>Citation</h5>
        <span style="line-height:1.5em">
            <p>Hebditch M, Carballo-Amador M.A., Charonis S, Curtis R, Warwicker J<br>
            <a href="https://doi.org/10.1093/bioinformatics/btx345">
            Protein-Sol: a web tool for predicting protein solubility from sequence.</a><br>
    Bioinformatics (2017)</p></span>
        </span>
        <br>
    <br>
    </div
    </div>
    </div>

<script type='text/javascript'>
    var timeStamp = '"""+timestamp+"""';
    document.getElementById('predictionButton').onclick = function () {
        location.href = '/code/download_file.php?app=solubility&dirname=run&timestamp='+timeStamp+'&idname='+timeStamp+'&file=seq_prediction.txt';
    };
</script>
"""

    #location to write results file
    with open("results.html","w") as html_file:
        html_file.write(html_header_text)
        html_file.write(html_footer_text)

def make_multiple_results(timestamp, count):

    results = {}
    with open ("/tmp/run-"+timestamp+"/all_predictions.csv") as pred_file:
        for line in pred_file:
            if line.startswith(">"):
                array_line = line.split(",")
                idname = array_line[0]
                scaledSol = array_line[2]
                results[idname] = float(scaledSol)


    most_sol_id = max(results, key=results.get)
    most_sol_val = results[most_sol_id]

    least_sol_val = min(results.values())

    html_header_text = """
    <!DOCTYPE html>
    <html lang='en'>
    
    <head>
      <title>Protein-sol solubility prediction</title>
      <meta charset='utf-8'>
      <meta name='viewport' content='width=device-width, user-scalable=no, minimum-scale=1.0, maximum-scale=1.0'>
      <link rel="icon" type="image/png" href="/resources/favicon-32x32.png" sizes="32x32" />
      <link rel="icon" type="image/png" href="/resources/favicon-16x16.png" sizes="16x16" />
      <link rel='stylesheet' href='/css/stylesheet.css?ver=1'>
      <link rel='stylesheet' href='/css/sequence.css?ver=1'>
    </head>
    
    <body>
    <div id='headerbar'> 
            <a class='headerlogo' href='/'>protein-sol</a>
            <a class='manunilogo' href='http://www.manchester.ac.uk'>UoM</a>
            <a class='bioprologo' href='http://www.biopronetuk.org'>BioProNET</a>
    </div>
    <?php include '/var/www/templates/navigation_links.php';?>
    <div id='content'>
    <h4>"""+str(count-1)+""" sequences submitted</h1>
    <br>
    <p>The predicted solubilities are in the range of """+str(least_sol_val)+"""-"""+str(most_sol_val)+"""</p>
    <br>
    <p>"""+most_sol_id+""" was the most soluble with a predicted score of """+str(most_sol_val)+"""</p>
    <br><br>
    """

    html_footer_text = """
    <div style='text-align:center;'>
    <button style='margin:auto;display:block;' id='predictionButton'>Download Prediction</button>
    </div>
    <br>
        <br>
        <h5>Citation</h5>
        <span style="line-height:1.5em">
            <p>Hebditch M, Carballo-Amador M.A., Charonis S, Curtis R, Warwicker J<br>
            <a href="https://doi.org/10.1093/bioinformatics/btx345">
            Protein-Sol: a web tool for predicting protein solubility from sequence.</a><br>
    Bioinformatics (2017)</p></span>
        </span>
        <br>
    <br>
    </div
    </div>
    </div>

<script type='text/javascript'>
    var timeStamp = '"""+timestamp+"""';
    document.getElementById('predictionButton').onclick = function () {
        location.href = '/code/download_file.php?app=solubility&dirname=run&timestamp='+timeStamp+'&idname='+timeStamp+'&file=results.zip';
    };
</script>
"""



    #location to write results file
    with open("results.html","w") as html_file:
        html_file.write(html_header_text)
        html_file.write(html_footer_text)


def handle_single(timestamp, seq_pred_file):

    build_single_graphs(timestamp, seq_pred_file)
    make_single_results(timestamp)

def handle_multiple(timestamp, count):
    make_multiple_results(timestamp, count)

#get timestamp
timestamp = sys.argv[1]
#get naame of prediction file
seq_pred_file = sys.argv[2]

filedir = "/tmp/run-"+timestamp+"/"

count = 0
with open(filedir+"all_predictions.csv","r") as predfile:
    for line in predfile:
        count+=1

if count == 2:
    handle_single(timestamp, seq_pred_file)
else:
    handle_multiple(timestamp, count)
