import subprocess, json
import pandas as pd


study_data = pd.read_excel("../data/metadata/biosamples_from_nonphyllo_cleaned.xlsx")

srr_dict = {}
print(study_data)

bioprojs = list(set(study_data.BioProject))
for bioproj in bioprojs:
        print(bioproj)
        cmd1 = F'esearch -db sra -query "{bioproj}[bioproject]" | efetch -mode xml | xtract -pattern EXPERIMENT_PACKAGE -element EXTERNAL_ID SRAFile@filename EXPERIMENT@alias'
        # cmd2 = F'esearch -db sra -query "{bioproj}[bioproject]" | efetch -mode xml | xtract -pattern EXPERIMENT_PACKAGE -element EXTERNAL_ID SAMPLE/TITLE'
        try:
            out1 = subprocess.run(cmd1, capture_output = True, shell = True).stdout.decode("utf-8")
        except:
            out1 = ""
        # try:
        #     out2 = subprocess.run(cmd2, capture_output = True, shell = True).stdout.decode("utf-8")
        # except:
        #     out2 = ""
        if "\t" in out1:
            res = out1.split("\n")
            print(res[0].split("\t"))
        srr_dict[bioproj] = out1
        # title_dict[bioproj] = out2

with open("../data/metadata/bioproject_query_sra_files_nonphyllo.json", "w") as fp:
    json.dump(srr_dict, fp)
