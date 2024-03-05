import sys
import pandas as pd
import numpy as np
fileName = sys.argv[1]
df = pd.read_table(fileName, comment='#',index_col=False, header=None)
fileName=fileName.replace(".vcf", "")
df.columns=["CHROM", "POS","ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT","NORMAL","TUMOR"]
df = df[df['FILTER'].str.contains('PASS')]
columns_to_drop = ["QUAL", "INFO", "FORMAT", "NORMAL"]
df.drop(columns=columns_to_drop, axis=1, inplace=True)
df[['DP','DP2','TAR','TIR','TOR','DP50','FDP50','SUBDP50', 'BCN50']] = df.TUMOR.str.split(":", expand = True)
columns_to_drop = ["TUMOR", "DP2", "TOR", "DP50", "FDP50", "SUBDP50", "BCN50"]
df.drop(columns=columns_to_drop, axis=1, inplace=True)
df.insert(7, "RD1", "")
df.insert(8, "RD2", "")
df[['RD1','RD2']]=df.TAR.str.split(",", expand=True)
df.insert(10, "AD1", "")
df.insert(11, "AD2", "")
df[['AD1','AD2']]=df.TIR.str.split(",", expand=True)
columns_to_drop = ["TAR", "RD2", "TIR", "AD2"]
df.drop(columns=columns_to_drop, axis=1, inplace=True)
df["RD1"]=pd.to_numeric(df["RD1"])
df["AD1"]=pd.to_numeric(df["AD1"])
df["DP"]=(((df["RD1"]) + (df["AD1"]))) 
df.insert(9, "AF", "")
df["AF"]=((df["AD1"]) / (df["DP"])) 
df["AF"]=(df["AF"]*100) 
df.insert(6, "DP1", "DP")
df.insert(8, "RD2", "RD")
df.insert(10, "AD2", "AD")
df.insert(12, "AF2", "AF")
df = df.astype(str)
df["DP3"] = df[["DP1", "DP"]].apply("=".join, axis=1)
df["RD3"] = df[["RD2", "RD1"]].apply("=".join, axis=1)
df["AD3"] = df[["AD2", "AD1"]].apply("=".join, axis=1)
df["AF3"] = df[["AF2", "AF"]].apply("=".join, axis=1)
columns_to_drop = ["DP1", "DP", "RD2", "RD1", "AD2", "AD1", "AF2", "AF"]
df.drop(columns=columns_to_drop, axis=1, inplace=True)
df["FORMAT"]=df[["DP3", "RD3", "AD3", "AF3"]].apply(";".join, axis=1)
columns_to_drop = ["DP3", "AD3", "RD3", "AF3"]
df.drop(columns=columns_to_drop, axis=1, inplace=True)
df.insert(6, "GT", "GT")
df.insert(7, "GT1", ".")
df.insert(8, "GT2", "")
df["GT2"] = df[["GT", "GT1"]].apply("=".join, axis=1)
df["FORMAT1"]=df[["GT2","FORMAT"]].apply(";".join, axis=1)
columns_to_drop = ["GT", "GT1", "GT2", "FORMAT"]
df.drop(columns=columns_to_drop, axis=1, inplace=True)
df.insert(6, "INFO", ".")
df.to_csv("_".join([fileName,'rf.vcf']), header=None, index=None, sep='\t')
