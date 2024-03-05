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
df[['DP','FDP','SDP','SUBDP','AU','CU','GU','TU']] = df.TUMOR.str.split(":", expand = True)
columns_to_drop = ["TUMOR", "FDP", "SDP", "SUBDP"]
df.drop(columns=columns_to_drop, axis=1, inplace=True)
df.insert(7, "A1", "")
df.insert(8, "A2", "")
df[['A1','A2']]=df.AU.str.split(",", expand=True)
df=df.drop("A2", axis=1)
df=df.drop("AU", axis=1)
df.insert(8, "C1", "")
df.insert(9, "C2", "")
df[['C1','C2']]=df.CU.str.split(",", expand=True)
df=df.drop("C2", axis=1)
df=df.drop("CU", axis=1)
df.insert(9, "G1", "")
df.insert(10, "G2", "")
df[['G1','G2']]=df.GU.str.split(",", expand=True)
df=df.drop("G2", axis=1)
df=df.drop("GU", axis=1)
df.insert(10, "T1", "")
df.insert(11, "T2", "")
df[['T1','T2']]=df.TU.str.split(",", expand=True)
df=df.drop("T2", axis=1)
df=df.drop("TU", axis=1)
conditions = [(df['REF'] == 'A'), (df['REF'] == 'C'),(df['REF'] == 'G'),(df['REF'] == 'T')]
rd=[df["A1"],df["C1"],df["G1"],df["T1"]]
df['RD'] = np.select(conditions, rd)    
conditions_1 = [(df['ALT'] == 'A'), (df['ALT'] == 'C'),(df['ALT'] == 'G'),(df['ALT'] == 'T')]
ad=[df["A1"],df["C1"],df["G1"],df["T1"]]
df['AD'] = np.select(conditions_1, ad) 
df["RD"]=pd.to_numeric(df["RD"])
df["AD"]=pd.to_numeric(df["AD"])
df["DP"]=(((df["RD"]) + (df["AD"]))) 
df["AF"]=((df["AD"]) / (df["DP"])) 
df["AF"]=(df["AF"]*100)
columns_to_drop = ["A1", "C1", "G1", "T1"]
df.drop(columns=columns_to_drop, axis=1, inplace=True)
df.insert(6, "DP1", "DP")
df.insert(8, "RD1", "RD")
df.insert(10, "AD1", "AD")
df.insert(12, "AF1", "AF")
df = df.astype(str)
df["DP2"] = df[["DP1", "DP"]].apply("=".join, axis=1)
df["RD2"] = df[["RD1", "RD"]].apply("=".join, axis=1)
df["AD2"] = df[["AD1", "AD"]].apply("=".join, axis=1)
df["AF2"] = df[["AF1", "AF"]].apply("=".join, axis=1)
df["FORMAT"]=df[["DP2", "RD2", "AD2", "AF2"]].apply(";".join, axis=1)
columns_to_drop = ["DP1", "DP", "RD1", "RD", "AD1", "AD", "AF1", "AF", "DP2", "RD2", "AD2", "AF2"]
df.drop(columns=columns_to_drop, axis=1, inplace=True)
df.insert(6, "GT", ".")
df.insert(7, "GT1", "GT")
df.insert(8, "GT2", ".")
df["GT3"]=df[["GT1","GT2"]].apply("=".join, axis=1)
df["FORMAT1"]=df[["GT3","FORMAT"]].apply(";".join, axis=1)
columns_to_drop = ["GT", "GT1", "GT2", "FORMAT", "GT3"]
df.drop(columns=columns_to_drop, axis=1, inplace=True)
df.to_csv("_".join([fileName,'rf.vcf']), header=None, index=None, sep='\t')
