#!/usr/bin/python3

# -*- coding: utf-8 -*-
"""
Created on Fri May 13 16:39:44 2022

@author: Roberto

Convert BLAST results (genomic) to .bed format for IGV

Usage:

1) Download the BLAST results as Hit Table (CSV)

2) Run the script with the name of the file as argument

3) The resulting bed file (same file name, ending .bed) can be loaded in IGV

"""

import pandas as pd
import sys

def main(args):
    
    data = pd.read_csv(args[1])
    
    # NCBI Hit Table column layout as of 13/05/2022
    data.columns = ["query", "subject", "pc_identity", "align_lengh", "n_mismatches", "n_gaps", 
                    "query_start", "query_end", "subject_start", "subject_end", "evalue", "bitscore"]
    
    # Saving only "NC" contigs to the .bed
    data = data[data["subject"].apply(lambda x: x.startswith("NC"))]
    
    # name = "Query start-end, identity%"
    data["name"] = data["query"] + " " + data["query_start"].astype(str) + "-" + \
            data["query_end"].astype(str) + ", " + data["pc_identity"].astype(str) + "%"
    data["chrom"] = "chr" + data["subject"].apply(lambda x: str(int(x.split(".")[0][7:])))
    data["chromStart"] = data[["subject_start", "subject_end"]].min(axis=1)
    data["chromEnd"] = data[["subject_start", "subject_end"]].max(axis=1)
    data["score"] = data["pc_identity"] * 100
    data["score"] = data["score"].astype(int)
    data["strand"] = data["subject_start"] > data["subject_end"]
    data["strand"] = data["strand"].apply(lambda x: ["-", "+"][x])
    data["thickStart"] = data["chromStart"]
    data["thickEnd"] = data["chromEnd"]
    # Column order conforms to the official UCSC bed format
    data = data[["chrom", "chromStart", "chromEnd", "name", "score", \
                 "strand", "thickStart", "thickEnd"]]
    print(data)
    data.to_csv(args[1].split(".")[0]+".bed", sep="\t", index=False, header=False)
        

if __name__ == "__main__":
    if len(sys.argv) > 0:
        main(sys.argv)
    else:
        print("Usage: blastigv.py <My-Alignment-HitTable.csv>")
