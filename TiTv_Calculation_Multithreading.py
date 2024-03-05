import sys
import pandas as pd
import numpy as np
from multiprocessing import Pool
fileName = sys.argv[1]
df = pd.read_table(fileName, comment='#',index_col=False, header=None)
fileName=fileName.replace(".vcf", "")
df.columns=["CHROM", "POS","ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT","NORMAL","TUMOR"]
df = df[df['FILTER'].str.contains('PASS')]

def calculate_ti_tv(row):
    ti = 1 if ((row['REF'], row['ALT']) in [('A', 'G'), ('G', 'A'), ('C', 'T'), ('T', 'C')]) else 0
    tv = 1 if ((row['REF'], row['ALT']) in [('A', 'C'), ('C', 'A'), ('G', 'T'), ('T', 'G'),
                                             ('A', 'T'), ('T', 'A'), ('G', 'C'), ('C', 'G')]) else 0
    return ti, tv

cores = 30
with Pool(cores) as pool:
    results = pool.map(calculate_ti_tv, df.to_dict(orient='records'))

df['Ti'], df['Tv'] = zip(*results)
ti_tv_ratio = df['Ti'].sum() / df['Tv'].sum()
print(f"Ti/Tv ratio: {ti_tv_ratio:.4f}")
