import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", type=str, help="input FASTA")
parser.add_argument("-a", "--acc", action="store_true", help="Remove duplicate accessions")
parser.add_argument("-o", "--output", type=str, help="output name for deduplicate FASTA")
args = parser.parse_args()

seq_dict ={}
with open(args.input, "r") as ifile:
    seq = ""
    lines = ifile.readlines()
    for line in lines:
        if line != "" and line[0] == ">":
            if seq != "":
                seq_dict[header] = seq + "\n"
            if args.acc:
                header = line.split(" ")[0].split("|")[0]
            else:
                header = line
            seq = ""
        elif line != "" and line != "\n":
            seq += line.strip()
    if seq != "":
        seq_dict[header] = seq + "\n"
buffer = ""
for i in seq_dict:
    if args.acc:
        buffer += F"{i}\n{seq_dict[i]}"
    else:
        buffer += F"{i}{seq_dict[i]}"
with open(args.output, "w") as ofile:
    ofile.write(buffer)
