import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--input", type=str, help="Input file name")
args = parser.parse_args()

input_file = args.input

class Results:
    def __init__(self, filename, header_text):
        self.filename = filename
        self.header_text = header_text
        self.dict = {}

def parse_line(line,obj):
    line = line.replace("\n","")
    line = line.replace(" ","")
    line = line.split(",")
    line = line[1:]

    obj.dict[line[0]] = line[1:]

def write_obj(obj):
    with open(obj.filename,"w") as output_file:
        output_file.write(obj.header_text)
        for protein, predictions in obj.dict.items():
            output_file.write(protein)
            output_file.write(",")
            for idx, value in enumerate(predictions):
                output_file.write(value)
                if idx < len(predictions)-1:
                    output_file.write(",")

            output_file.write("\n")

with open (input_file,"r") as results:

    Predictions = Results("all_predictions.csv",
                          "ID,percentage_sol,scaled_sol,population_sol,pI\n")

    Deviations = Results("all_deviations.csv",
                         "ID, KmR, DmE, len, A, C, D, E, F, G, H, I, K, L, M, N, P, Q, R, S, T,V, W, Y, KpR, DpE, PmN, "
                         "PpN," " aro, pI, mem, chr, fld, dis, ent, bet\n")

    KD = Results("all_KyteDoolittle.csv",
                 "KyteDoolittle sequence profile\n")

    FI = Results("all_foldIndex.csv",
                 "FoldIndex sequence profile\n")

    Entropy = Results("all_Entropy.csv",
                 "entropy sequence profile\n")

    Charge = Results("all_charge.csv",
                      "Charge sequence profile\n")

    for line in results:
        if line.startswith("SEQUENCE PREDICTIONS"):
            parse_line(line,Predictions)
        elif line.startswith("SEQUENCE DEVIATIONS"):
            parse_line(line,Deviations)
        elif line.startswith("SEQUENCE PROFILE"):
            if "KyteDoolittle" in line:
                parse_line(line,KD)
            elif "Uversky" in line:
                parse_line(line,FI)
            elif "entropy" in line:
                parse_line(line,Entropy)
            elif "charge" in line:
                parse_line(line,Charge)

for obj in [Predictions,Deviations,KD,FI,Entropy,Charge]:
    write_obj(obj)

