import json, copy
import pandas as pd

with open("../data/metadata/bioproject_query_sra_files_nonphyllo.json", "r") as fp:
    srr_dict = json.load(fp)

study_metadata = pd.read_excel("../data/metadata/Study_metadata_nonphyllo.xlsx").set_index("BioProject").to_dict(orient="index")
cleaned_sam_data = pd.read_excel("../data/metadata/biosamples_from_nonphyllo_cleaned.xlsx").set_index("BioSample").to_dict(orient="index")

bac_terms = ["785F-1064R", "16S", "V4", "341F-785R", "799F", "bac", "Bac"]
fun_terms = ["ITS", "fun", "Fun", "Ftrad"]
oom_terms = ["Otrad"]

srr_to_sam = {}
for bioproj in srr_dict:
    ext_met = study_metadata[bioproj]["Extraction_method"]
    res = srr_dict[bioproj]
    records = res.split("\n")
    for rec in records:
        tab_spl = rec.split("\t")
        # print(tab_spl)
        sam_acc = [x for x in tab_spl if "SAM" in x]
        if len(sam_acc) == 0:
            print(tab_spl)
        elif sam_acc[0] not in cleaned_sam_data:
            print(F"{sam_acc[0]} is not in cleaned list")
        else:
            sam_acc = sam_acc[0]
            srr_acc = tab_spl[-2].split(".")[0] #[x for x in tab_spl if "RR" in x]
            alias = tab_spl[-1]
            files = [x for x in tab_spl if "." in x]
            if len(files) == 1:
                layout = "single"
                fwd = files[0]
                rev = ""
            else:
                layout = "paired"
                fwd, rev = files[0], files[1]
            if ";" not in study_metadata[bioproj]["Target"]:
                target = study_metadata[bioproj]["Target"]
                primers = study_metadata[bioproj]["Primers"]
                region = study_metadata[bioproj]["Region"]
            elif pd.notnull(cleaned_sam_data[sam_acc]["Target"]):
                target = cleaned_sam_data[sam_acc]["Target"]
                primer_spl = study_metadata[bioproj]["Primers"].split(";")
                region_spl = study_metadata[bioproj]["Region"].split(";")
                if target == "Bacteria":
                    primers = primer_spl[0]
                    region = region_spl[0]
                elif target == "Fungi":
                    primers = primer_spl[1]
                    region = region_spl[1]
                elif target == "Oomycota":
                    primers = primer_spl[2]
                    region = region_spl[2]
            else:
                target = ""
                primers = ""
                for term in bac_terms:
                    if term in files[0]:
                        target = "Bacteria"
                for term in fun_terms:
                    if term in files[0]:
                        target = "Fungi"
                for term in oom_terms:
                    if term in files[0]:
                        target = "Oomycota"
                if target == "":
                    for term in bac_terms:
                        if term in alias:
                            target = "Bacteria"
                    for term in fun_terms:
                        if term in alias:
                            target = "Fungi"
                    for term in oom_terms:
                        if term in alias:
                            target = "Oomycota"
                    if "PPK" in alias and "-B" in alias:
                        target = "Bacteria"
                    if "PPK" in alias and "-F" in alias:
                        target = "Fungi"
                if target != "":
                    primer_spl = study_metadata[bioproj]["Primers"].split(";")
                    region_spl = study_metadata[bioproj]["Region"].split(";")
                    if target == "Bacteria":
                        primers = primer_spl[0]
                        region = region_spl[0]
                    elif target == "Fungi":
                        primers = primer_spl[1]
                        region = region_spl[1]
                    elif target == "Oomycota":
                        primers = primer_spl[2]
                        region = region_spl[2]

            srr_to_sam[srr_acc] = {"BioSample" : sam_acc, "BioProject" : bioproj, "Layout" : layout, "Target" : target, "Primer_set" : primers,
                                    "Region" : region, "Extraction_method" : ext_met, "Fwd_file" : fwd, "Rev_file" : rev}


with open("../data/metadata/bioproject_query_srr_to_sam_nonphyllo.json", "w") as fp:
    json.dump(srr_to_sam, fp)

with open("../data/metadata/bioproject_query_srr_to_sam_nonphyllo.json", "r") as fp:
    srr_to_sam = json.load(fp)




srr_data_dict = {}
for srr_acc in srr_to_sam:
    sam_acc = srr_to_sam[srr_acc]["BioSample"]
    sam_acc = srr_to_sam[srr_acc]["BioSample"]
    srr_data_dict[srr_acc] = cleaned_sam_data[sam_acc].copy()
    srr_data_dict[srr_acc]["BioSample"] = sam_acc
    for field in ["Target", "Layout", "Primer_set", "Region", "Extraction_method", "Fwd_file", "Rev_file"]:
        srr_data_dict[srr_acc][field] = srr_to_sam[srr_acc][field]

srr_data = pd.DataFrame.from_dict(srr_data_dict, orient= "index")

cols_to_drop = []
for i in srr_data.columns:
    if srr_data[i].isnull().sum() == len(srr_data[i]):
        cols_to_drop.append(i)
for i in cols_to_drop:
    srr_data.drop(i, axis=1, inplace=True)

srr_data.to_excel("../data/metadata/srr_accessions_with_metadata_nonphyllo.xlsx")
srr_data.to_csv("../data/metadata/srr_accessions_with_metadata_nonphyllo.txt", sep = "\t")

regs = ["ITS1", "ITS1-ITS2", "ITS2"]
targets = ["Fungi"]
layouts = ["single", "paired"]

buffer = ""
for r in regs:
    for t in targets:
        for layout in layouts:
            # print((srr_data['Target'] == t) & (srr_data['Region'] == r) & (srr_data['Layout'] == layout))
            sub_df = srr_data.loc[(srr_data['Target'] == t) & (srr_data['Region'] == r) & (srr_data['Layout'] == layout)]
            print(F"{len(sub_df)} records have target = {t}, region = {r}, and layout = {layout}.")
            if len(sub_df) > 0:
                with open(F"../data/metadata/srr_list_{t}_{r}_{layout}.txt", "w") as ofile:
                    ofile.write("\n".join([x.split(".")[0] for x in sub_df.index]))
                buffer += F"{t} {r} {layout}\n"
with open("../data/metadata/srr_targets_regions.txt", "w") as ofile:
    ofile.write(buffer)
