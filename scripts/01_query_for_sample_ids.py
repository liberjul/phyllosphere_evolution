import subprocess, json
import pandas as pd


study_data = pd.read_excel("../data/metadata/Study_metadata_nonphyllo.xlsx")

tag_dict = {}
title_dict = {}
for i in range(len(study_data)):
        bioproj = study_data.BioProject[i]
        print(bioproj)
        # esearch -db sra -query "PRJNA480528[bioproject]" | efetch -mode xml | xtract -pattern EXPERIMENT_PACKAGE -element EXTERNAL_ID TAG VALUE
        cmd1 = F'esearch -db sra -query "{bioproj}[bioproject]" | efetch -mode xml | xtract -pattern EXPERIMENT_PACKAGE -element EXTERNAL_ID TAG VALUE'
        cmd2 = F'esearch -db sra -query "{bioproj}[bioproject]" | efetch -mode xml | xtract -pattern EXPERIMENT_PACKAGE -element EXTERNAL_ID SAMPLE/TITLE'
        try:
            out1 = subprocess.run(cmd1, capture_output = True, shell = True).stdout.decode("utf-8")
        except:
            out1 = ""
        try:
            out2 = subprocess.run(cmd2, capture_output = True, shell = True).stdout.decode("utf-8")
        except:
            out2 = ""
        if "\t" in out1:
            res = out1.split("\n")
            print(res[0].split("\t"))
        tag_dict[bioproj] = out1
        title_dict[bioproj] = out2

with open("../data/metadata/bioproject_query_samples_nonphyllo.json", "w") as fp:
    json.dump(tag_dict, fp)

with open("../data/metadata/bioproject_query_sample_titles_nonphyllo.json", "w") as fp:
    json.dump(title_dict, fp)
