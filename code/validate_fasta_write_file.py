from tidyfasta.common.process import ProcessFasta
import argparse
import re
import os

parser = argparse.ArgumentParser()
parser.add_argument("--input", type=str, help="Input file name")
parser.add_argument("--min", type=int, required=False, help="Minimum number of AA")
parser.add_argument("--max", type=int, required=False, help="Maximum number of AA")
args = parser.parse_args()

input_file = args.input
min_aa = args.min
max_aa = args.max

def write_error(message, input_file):

    outputfile = str(input_file + "-ERROR")
    with open(outputfile, "w") as output:
        output.write(message+"\n")


try:
    split_path = os.path.split(os.path.abspath(input_file))
    candidate_output_file = split_path[0] + "/tidied-" + split_path[1]

    PF = ProcessFasta(input_file, strict=True, single=False)

    fasta_array = PF.validated_array

    autocount = 0


    final_array = []

    for fasta in fasta_array:
        if min_aa or max_aa:
            if min_aa and len(fasta.sequence) < min_aa:
                write_error("Sequence "+str(fasta.id)+" has "+str(len(fasta.sequence))+", the calculation requires a sequence of at least "+str(min_aa)+" amino acids", input_file)
            if max_aa and len(fasta.sequence) > max_aa:
                write_error("Sequence "+str(fasta.id)+" has "+str(len(fasta.sequence))+", the calculation requires a sequence of at most "+str(max_aa)+" amino acids", input_file)

        if fasta.id.startswith("> sequence"):
            fasta.id = ">Protein-sol-"+str(autocount)
            autocount+=1

        if fasta.id.startswith("> "):
            fix = fasta.id.replace("> ",">")
            fasta.id = fix

        final_array.append(fasta)


    index = 0
    with open(candidate_output_file, "w") as output:
        for object in final_array:
            index += 1
            output.write(object.id + "\n")
            if index == len(final_array):
                output.write(object.sequence + "\n")
            else:
                output.write(object.sequence + "\n\n")


    os.rename(candidate_output_file,input_file)


except Exception as e:
    error_msg = str(e)
    if re.match("^Output location",error_msg):
        write_error("Sequence not saved correctly", input_file)
    elif re.match("^Path.+doesn't exist",error_msg):
        write_error("Sequence not saved correctly", input_file)
    elif re.match("^No data in file",error_msg):
        write_error("Sequence not saved correctly", input_file)
    elif re.match("^Combined array not generated",error_msg):
        write_error("Sequence not saved correctly", input_file)
    elif re.match("^Cleaned array not generated",error_msg):
        write_error("Sequence not saved correctly", input_file)
    elif re.match("^Named array issue",error_msg):
        write_error("Issue processing ID and sequence line", input_file)
    elif re.match("^Unpaired ID and sequence",error_msg):
        write_error("Issue processing ID and sequence line", input_file)
    elif re.match("^Object array failed",error_msg):
        write_error("Issue processing ID and sequence line", input_file)
    elif re.match("^Issue identifying ID",error_msg):
        write_error("Issue processing ID", input_file)
    elif re.match("^ID without sequence",error_msg):
        write_error("ID only with no sequence", input_file)
    elif re.match("^No valid AA in input",error_msg):
        write_error("Non-standard AA identified, please check your sequence only contains the 20 canonical amino acid characters", input_file)
    elif re.match("^Non canonical",error_msg):
        write_error("Non-standard AA identified, please check your sequence only contains the 20 canonical amino acid characters", input_file)
    elif re.match("^Bad characters in ID",error_msg):
        write_error("Non-allowed characters in ID line, please check your ID only contains alphanumeric characters, >, _ or -", input_file)
    elif re.match("^More than 1 sequence",error_msg):
        write_error("More than one sequence submitted", input_file)
    else:
        write_error("Error processing sequence", input_file)
        print(e)

