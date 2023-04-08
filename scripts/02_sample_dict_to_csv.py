import json
import pandas as pd

with open("../data/metadata/bioproject_query_samples_nonphyllo.json", "r") as fp:
    tag_dict = json.load(fp)
with open("../data/metadata/bioproject_query_sample_titles_nonphyllo.json", "r") as fp:
    title_dict = json.load(fp)

parsed_dict = {}
for bioproj in title_dict:
    res1 = tag_dict[bioproj]
    res2 = title_dict[bioproj]
    for biosamp, title in zip(res1.split("\n"), res2.split("\n")):
        tab_spl = biosamp.split("\t")
        if len(tab_spl) > 1:
            attr = [x for x in tab_spl if ("SAM" not in x and "PRJ" not in x)]
            sample_name = [x for x in tab_spl if "SAM" in x]
            if len(sample_name) == 0:
                print(tab_spl)
            else:
                sample_name = sample_name[0]
                bioproject = [x for x in tab_spl if "PRJ" in x][0]
                parsed_dict[sample_name] = {"bioproject" : bioproject}
                if "SAMN" not in title.split("\t")[-1]:
                    parsed_dict[sample_name]["Title"] = title.split("\t")[-1]
                for x,y in zip(attr[:int(len(attr)/2)], attr[int(len(attr)/2):]):
                    # print(sample_name, x, y)
                    # if x not in parsed_dict[tab_spl[1]].keys():
                    parsed_dict[sample_name][x] = y

with open("../data/metadata/bioproject_query_samples_parsed_nonphyllo.json", "w") as fp:
    json.dump(parsed_dict, fp)

with open("../data/metadata/bioproject_query_samples_parsed_nonphyllo.json", "r") as fp:
    parsed_dict = json.load(fp)

sample_data = pd.DataFrame.from_dict(parsed_dict, orient= "index")
old_data = pd.read_excel("../data/metadata/biosamples_from_nonphyllo_studies.xlsx", sheet_name=0, index_col=0)
print(sample_data.index)
print(old_data.index)

combined_data = pd.concat([sample_data, old_data])
combined_data.to_excel("../data/metadata/biosamples_from_studies_nonphyllo.xlsx")
combined_data.to_csv("../data/metadata/biosamples_from_studies_nonphyllo.txt", sep = "\t")
