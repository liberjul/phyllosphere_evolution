import pandas as pd
from Bio.Seq import Seq
import subprocess, os, glob

primer_dict = {
"338F-806R" : ["ACTCCTACGGGAGGCAGCA", "GGACTACNVGGGTWTCTAAT"],
"341F-785R" : ["CCTACGGGNGGCWGCAG", "GACTACHVGGGTATCTAATCC"],
"515F-806R" : ["GTGYCAGCMGCCGCGGTAA", "GGACTACNVGGGTWTCTAAT"],
"341F-805R" : ["CCTACGGGNGGCWGCAG", "GACTACHVGGGTATCTAATCC"],
"5.8SFun-ITS4Fun" : ["AACTTTYRRCAAYGGATCWCT", "AGCCTCCGCTTATTGATATGCTTAART"],
"501F-706R" : ["CAGCCTACGGGNGGCWGCAG", "ACTACHVGGGTATCTAATCC"],
"515F-926R" : ["GTGYCAGCMGCCGCGGTAA", "CCGYCAATTYMTTTRAGTTT"],
"785F-1064R" : ["GGATTAGATACCC", "CGACRRCCATGCANCACCT"],
"799F-1115R" : ["AACMGGATTAGATACCCKG", "AGGGTTGCGCTCGTTG"],
"799F-1193R" : ["AACMGGATTAGATACCCKG", "ACGTCATCCCCACCTTCC"],
"799F-1492R" : ["AACMGGATTAGATACCCKG", "TACCTTGTTACGACTT"],
"926F-1392R" : ["AAACTYAAAKGAATTGACGG", "ACGGGCGGTGTGTRC"],
"B341F-B806R" : ["CCTACGGGAGGCAGCAG", "GGACTACHVGGGTWTCTAAT"],
"B799F-B1194R" : ["AACMGGATTAGATACCCKG", "ACGTCATCCCCACCTTCC"],
"BAC799F-BAC1193R" : ["AACMGGATTAGATACCCKG", "ACGTCATCCCCACCTTCC"],
"fITS7-ITS4" : ["GTGAATCATCGAATCTTTG" , "TCCTCCGCTTATTGATATGC"],
"fITS7-ITS4R" : ["GTGAATCATCGAATCTTTG", "TCCTCCGCTTATTGATATGC"],
"Illumina V1-V3" : ["AGAGTTTGATCCTGGCTCAG", "ATTACCGCGGCTGCTGG"],
"ITS-1F-F-ITS1-1F-R" : ["CTTGGTCATTTAGAGGAAGTAA", "GCTGCGTTCTTCATCGATGC"],
"ITS1F-ITS2" : ["CTTGGTCATTTAGAGGAAGTAA", "GCTGCGTTCTTCATCGATGC"],
"ITS1F-ITS2R" : ["CTTGGTCATTTAGAGGAAGTAA", "GCTGCGTTCTTCATCGATGC"],
"ITS1f-ITS4" : ["TTGGTCATTTAGAGGAAGTA", "TCCTCCGCTTATTGATATGC"],
"ITS1F-KYO1-ITS2-KYO1" : ["CTHGGTCATTTAGAGGAASTAA", "TAGAGGAAGTAAAAGTCGTAA"],
"ITS1O-5.8s-O-R" : ["CGGAAGGATCATTACCAC", "AGCCTAGACATCCACTGCTG"],
"ITS2-F-ITS2-R" : ["GCATCGATGAAGAACGC", "CCTCCGCTTATTGATATGC"],
"ITS7-ITS4" : ["GTGAATCATCGAATCTTTG", "TCCTCCGCTTATTGATATGC"],
"TS1F12-5.8S" : ["GAACCWGCGGARGGATCA", "CGCTGCGTTCTTCATCG"],
"ITS3-ITS4" : ["GCATCGATGAAGAACGCAGC", "TCCTCCGCTTATTGATATGC"],
"ITS3_KYO2-ITS4" : ["TAGAGGAAGTAAAAGTCGTAA", "TCCTCCGCTTATTGATATGC"]
}

srr_dat = pd.read_csv("../data/metadata/srr_accessions_with_metadata_nonphyllo.txt", sep = "\t",  index_col=0)

srr_files = glob.glob("../data/SRA_files/*/SRR[0-9]*.fastq.gz")
srr_files = [os.path.basename(x).split(".")[0] for x in srr_files]
srrs_incl = set([x.split("_")[0] for x in srr_files])

print(srrs_incl)

# subprocess.run("conda init bash; conda activate qiime2-2023.2", shell=True)
for i in range(len(srr_dat)):
    if srr_dat.index[i] in srrs_incl and srr_dat["Target"][i] == "Fungi":
        primer_set, layout = srr_dat["Primer_set"][i], srr_dat["Layout"][i]
        if "SRR749" in srr_dat.index[i]:
            layout = "single"
        fwd_primer_seq, rev_primer_seq = primer_dict[primer_set]
        fwd_rc, rev_rc = str(Seq(fwd_primer_seq).reverse_complement()), str(Seq(rev_primer_seq).reverse_complement())
        if layout == "paired":
            inpath = F'../data/SRA_files/uncut_seqs/{srr_dat.index[i]}'
            outpath = F'../data/SRA_files/{srr_dat["Target"][i]}_{srr_dat["Region"][i]}_paired/JL_{srr_dat.index[i]}_L001_R'
            if len(glob.glob(F'../data/SRA_files/{srr_dat["Target"][i]}_{srr_dat["Region"][i]}_paired/JL_...._{srr_dat.index[i]}_L001_R._001.fastq.gz')) > 0:
                pass
            elif os.path.exists(inpath + "_1.sampled.fastq.gz"):
                cmd = F"cutadapt -a {fwd_primer_seq}...{rev_rc} -A {rev_primer_seq}...{fwd_rc} -e .1 -j 24 -m 150 -o {outpath}1_001.cutadapt.fastq.gz -p {outpath}2_001.cutadapt.fastq.gz {inpath}_1.sampled.fastq.gz {inpath}_2.sampled.fastq.gz"
                subprocess.run(cmd, shell = True)
            elif os.path.exists(inpath + "_1.fastq.gz") and os.path.exists(inpath + "_2.fastq.gz") and not os.path.exists(F"{outpath}1_001.cutadapt.fastq.gz"):
                cmd = F"cutadapt -a {fwd_primer_seq}...{rev_rc} -A {rev_primer_seq}...{fwd_rc} -e .1 -j 24 -m 150 -o {outpath}1_001.cutadapt.fastq.gz -p {outpath}2_001.cutadapt.fastq.gz {inpath}_1.fastq.gz {inpath}_2.fastq.gz"
                subprocess.run(cmd, shell = True)
        elif layout == "single":
            inpath = F'../data/SRA_files/uncut_seqs/{srr_dat.index[i]}'
            outpath = F'../data/SRA_files/{srr_dat["Target"][i]}_{srr_dat["Region"][i]}_single/JL_{srr_dat.index[i]}_L001_R'
            if os.path.exists(inpath + ".fastq.gz") and not os.path.exists(F"{outpath}1_001.cutadapt.fastq.gz"):
                cmd = F"cutadapt -a {fwd_primer_seq}...{rev_rc} -e .1 -j 24 -m 150 -o {outpath}1_001.cutadapt.fastq.gz {inpath}.fastq.gz"
                subprocess.run(cmd, shell = True)
        else:
            print(F"No layout provided for {srr_dat.index[i]}")
