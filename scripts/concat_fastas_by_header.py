import glob
# import argparse
#
# parser = argparse.ArgumentParser()
# parser.add_argument("-i", "--input", type=str, help="input FASTAs, separated by commas")
# parser.add_argument("-a", "--acc", action="store_true", help="Remove duplicate accessions")
# parser.add_argument("-o", "--output", type=str, help="output name for deduplicate FASTA")
# args = parser.parse_args()
fastas = glob.glob("../data/reference_data/class*.fasta")


seq_dict ={}
headers = []
for fasta in fastas:
    reg = fasta.split(".dedup")[0].split("seqs.")[1]
    seq_dict[reg] = {}
    with open(fasta, "r") as ifile:
        header_count = 0
        seq = ""
        lines = ifile.readlines()
        for line in lines:
            if line != "" and line[0] == ">":
                if seq != "":
                    seq_dict[reg][header] = seq.replace("N", "")
                header = line
                seq = ""
            elif line != "" and line != "\n":
                seq += line.strip()
        if seq != "":
            seq_dict[reg][header] = seq.replace("N", "")

header_set = set(seq_dict["ITS1"].keys())
for reg in ["5_8S", "ITS2"]:
    header_set &= set(seq_dict[reg].keys())
header_list = list(header_set)

buffer = ""
for header in header_list:
    print(header[1:].strip())
    buffer += header
    for reg in ["ITS1", "5_8S", "ITS2"]:
        buffer += F"{seq_dict[reg][header]}"
    buffer += "\n"
with open("../data/reference_data/tree_ITS.fasta", "w") as ofile:
    ofile.write(buffer)
