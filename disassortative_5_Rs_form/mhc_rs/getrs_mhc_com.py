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

df=pd.read_csv("/data2/wangxuedong/mhc_test_data/python/getRS/df_sorted_0513_top100_csv3.csv")
df_mhc= df[df['isMHC'] == 1]
positions=df_mhc["Position of SNP"]
allpositions=list(set(positions))
sorted_pos=sorted(allpositions)

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
    dataframe(dictforDF).to_csv("/data2/wangxuedong/mhc_test_data/getRS_0616/forms/COM_MHC.csv")
    return dataframe(dictforDF)


generatedf(columns=["chrom","chr","pos","ref","alts[0]","not know"],allrecords=allrecords_mhc)