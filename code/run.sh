#!/bin/bash

fileName=$1
timeStamp=$2
proteinName=$3

utime=$(date +%s)

#go to the temp dir where Jim's code is saved
cd /tmp/run-$timeStamp

#run jim's code
bash /var/www/html/code/solubility_code/run_model.sh $fileName

#delete un-needed files from Jim's code
rm /tmp/run-$timeStamp/sequenceIn-"$timeStamp"_ORIGINAL

#parseoutput file
python3 /var/www/html/code/parse_solubility.py --input seq_prediction.txt

#determine number of sequences
numberSequences=$(grep -c ">" predictions.csv)

#run python code to handle output
python3 /var/www/html/code/sequence_output_handler.py $timeStamp "/tmp/run-$timeStamp/seq_prediction.txt"

zip -jr /tmp/run-$timeStamp/results.zip /tmp/run-$timeStamp/{*_*.txt,*_*.csv,*.png,sequenceIn*}
