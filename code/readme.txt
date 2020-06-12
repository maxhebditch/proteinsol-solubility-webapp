#sequenceprediction.php
Main code that calls the other functions
1) Handles input from user
2) Error checks it
3) Calls perl code to run algorithm
4) Calls python code to make graphs
5) Calls bash code to make directories and move files etc
6) Calls bash to generate results html

#run.sh
Bash code that 
1) enters directory
2) calls Jim's perl code
3) calls python code to make graphs

#solubility_code
Directory containing Jim's code to run the algorithm
Contains run_model.sh which is the bash script wrapper

#sequence_output_handler.py
python code that extracts data from Jim's code and makes the graphs

#results-build.sh
Bash code that generates results html

#env
The python virtual environment

#requirements.txt
List of python packages required to run the software
