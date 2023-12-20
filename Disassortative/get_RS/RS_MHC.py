#!/home/wangxuedong/app/miniconda3/envs/mhc/bin/python
import pysam
from pysam import VariantFile as vcf
import operator
from math import log2
import pandas as pd
from pandas import DataFrame as dataframe
import matplotlib.pyplot as plt
import numpy as np

from scipy.spatial.distance import pdist, squareform
import scipy
import  os
import os.path
import matplotlib.colors as mcolors
from scipy import stats
import csv
import requests
import json
from xml.etree import ElementTree as ET
from statsmodels.stats.multitest import multipletests
import re

df=pd.read_csv("/data2/wangxuedong/mhc_test_data/first_review/a_rsform/csvs/rs_tables4.csv")

df=df[df['dbSNP'].isnull()]

df_split= df[df['isMHC'] == 1]
positions=df_split["Position of SNP"].values.tolist()
chrs=df_split["Position of gene.1"].values.tolist()
allrecords=[]
for i in range(len(chrs)):
    arr=[]
    ch=chrs[i].split(":")[0]
    arr.append(ch)
    arr.append(positions[i])
    allrecords.append(arr)
allrecords = list(set(tuple(x) for x in allrecords))
allrecords = [list(x) for x in allrecords]
allrecords=[a[1] for a in allrecords]
sorted_pos=sorted(allrecords)

sorted_pos=sorted_pos

allrecords_mhc=[]
vcf_in_mhc=vcf("/data2/wangxuedong/mhc_test_data/mhcdataset/1kgp.29720000-33130000.vcf","r")
for mhc_pos in sorted_pos:
    for rec in vcf_in_mhc.fetch():
        if rec.pos==mhc_pos:
            onerecord_mhc=["chromosome"]
            onerecord_mhc+=["6",rec.pos,rec.ref,rec.alts[0],1]
            allrecords_mhc.append(onerecord_mhc)
        else:
            continue
        
        

def generatedf(columns,allrecords):
    dictforDF=dict()
    for i in range(len(columns)):
        midarr=[]
        for record in allrecords:           
            midarr.append(record[i])
        dictforDF[columns[i]]=midarr
    dataframe(dictforDF).to_csv("/data2/wangxuedong/mhc_test_data/first_review/a_rsform/rs_dis_res/dis_mhc_p0.csv")
    return dataframe(dictforDF)


generatedf(columns=["chrom","chr","pos","ref","alts[0]","not know"],allrecords=allrecords_mhc)


