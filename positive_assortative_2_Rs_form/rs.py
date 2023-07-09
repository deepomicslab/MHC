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
import pickle

split_race_region={'ACB': {'split_chr1/xbb': [[123992935, 124938490]],
  'split_chr1/xba': [[121579777, 123992934]],
  'split_chr13/xaa': [[16000187, 18181463]],
  'split_chr2/xav': [[91720904, 94530414]],
  'split_chr3/xau': [[90103657, 93851087]]},
  'All':{'split_chr13/xaa': [[16000187, 18181463]],
 'split_chr1/xba': [[121579777, 123992934]],
 'split_chr1/xbb': [[123992935, 124938490]],
 'split_chr2/xav': [[91720904, 94530414]],
 'split_chr3/xau': [[90362505, 91477491]]},
 'ASW': {'split_chr1/xbb': [[123992935, 124938490]],
  'split_chr1/xba': [[121635735, 123992934]],
  'split_chr13/xaa': [[16040220, 18181463]],
  'split_chr16/xaj': [[35510197, 46630361]],
  'split_chr2/xav': [[91720904, 94530414]],
  'split_chr3/xau': [[90103657, 93851087]]},
 'BEB': {'split_chr1/xbb': [[123992935, 124938490]],
  'split_chr1/xba': [[121579777, 123992934]],
  'split_chr13/xaa': [[16165134, 17417031]],
  'split_chr16/xaj': [[35510197, 46692526]],
  'split_chr2/xav': [[91992563, 94195617]]},
 'CDX': {},
 'CEU': {'split_chr1/xbb': [[123992935, 124920289]],
  'split_chr1/xba': [[121579777, 123992934]],
  'split_chr13/xaa': [[16446009, 18181463]],
  'split_chr2/xav': [[92076199, 94195617]]},
 'CHB': {},
 'CHS': {'split_chr1/xba': [[123102017, 123978193]]},
 'CLM': {'split_chr1/xbb': [[123992935, 124938490]],
  'split_chr1/xba': [[121579777, 123992934]],
  'split_chr13/xaa': [[16446009, 18181463]],
  'split_chr16/xaj': [[35510197, 46630361]],
  'split_chr16/xai': [[32030198, 35510168]],
  'split_chr2/xav': [[91992563, 94654240]]},
 'ESN': {'split_chr1/xbb': [[123992935, 124938490]],
  'split_chr1/xba': [[121635735, 123992934]],
  'split_chr13/xaa': [[16000187, 18181463]],
  'split_chr2/xav': [[91720904, 94654240]],
  'split_chr3/xau': [[90103657, 93851087]]},
 'FIN': {'split_chr1/xbb': [[124048603, 124877211]],
  'split_chr1/xba': [[121635735, 123992934]],
  'split_chr2/xav': [[91992563, 94195617]]},
 'GBR': {'split_chr1/xbb': [[123992935, 124877211]],
  'split_chr1/xba': [[122678788, 123992934]]},
 'GIH': {'split_chr1/xbb': [[123992935, 124877211]],
  'split_chr1/xba': [[121635735, 123992934]],
  'split_chr16/xaj': [[35510197, 46692526]],
  'split_chr2/xav': [[91951194, 94195617]]},
 'GWD': {'split_chr1/xbb': [[123992935, 124938490]],
  'split_chr1/xba': [[121579777, 123992934]],
  'split_chr13/xaa': [[16000187, 18181463]],
  'split_chr16/xaj': [[35510197, 46630361]],
  'split_chr2/xav': [[91720904, 94195617]],
  'split_chr3/xau': [[90103657, 93851087]]},
 'IBS': {'split_chr1/xbb': [[123992935, 124938490]],
  'split_chr1/xba': [[121579777, 123992934]],
  'split_chr13/xaa': [[16446009, 18024589]],
  'split_chr16/xaj': [[35510197, 46630361]],
  'split_chr2/xav': [[91992563, 94654240]]},
 'ITU': {'split_chr1/xbb': [[123992935, 124877211]],
  'split_chr1/xba': [[121635735, 123992934]],
  'split_chr16/xaj': [[35510197, 46692526]],
  'split_chr2/xav': [[91992563, 94195617]]},
 'JPT': {},
 'KHV': {},
 'LWK': {'split_chr1/xbb': [[123992935, 124938490]],
  'split_chr1/xba': [[121579777, 123992934]],
  'split_chr13/xaa': [[16040220, 18181463]],
  'split_chr2/xav': [[91720904, 94654240]],
  'split_chr3/xau': [[90103657, 93851087]]},
 'MSL': {'split_chr1/xbb': [[123992935, 124938490]],
  'split_chr1/xba': [[121635735, 123992934]],
  'split_chr13/xaa': [[16040220, 18181463]],
  'split_chr2/xav': [[91720904, 94195617]],
  'split_chr3/xau': [[90103657, 93851087]]},
 'MXL': {'split_chr1/xbb': [[123992935, 124920289]],
  'split_chr1/xba': [[121635735, 123992934]],
  'split_chr16/xaj': [[35510197, 46692526]],
  'split_chr2/xav': [[91720904, 94735664]]},
 'PEL': {'split_chr1/xbb': [[123992935, 124782972]],
  'split_chr1/xba': [[122885884, 123992934]]},
 'PJL': {'split_chr1/xbb': [[123992935, 124938490]],
  'split_chr1/xba': [[121579777, 123992934]],
  'split_chr13/xaa': [[16165134, 18024589]],
  'split_chr16/xaj': [[35510197, 46692526]],
  'split_chr2/xav': [[91768349, 94195617]]},
 'PUR': {'split_chr1/xbb': [[123992935, 124938490]],
  'split_chr1/xba': [[121579777, 123992934]],
  'split_chr13/xaa': [[16152209, 18181463]],
  'split_chr16/xaj': [[35510197, 46630361]],
  'split_chr2/xav': [[91720904, 94809756]]},
 'STU': {'split_chr1/xbb': [[123992935, 124920289]],
  'split_chr1/xba': [[121635735, 123992934]],
  'split_chr16/xaj': [[35510197, 46692526]],
  'split_chr16/xai': [[33059491, 35510168]],
  'split_chr2/xav': [[91992563, 94195617]]},
 'TSI': {'split_chr1/xbb': [[123992935, 124920289]],
  'split_chr1/xba': [[121579777, 123992934]],
  'split_chr2/xav': [[91992563, 94195617]]},
 'YRI': {'split_chr1/xbb': [[123992935, 124938490]],
  'split_chr1/xba': [[121635735, 123992934]],
  'split_chr13/xaa': [[16000187, 18181463]],
  'split_chr2/xav': [[91720904, 94530414]],
  'split_chr3/xau': [[90103657, 93851087]]}}



result_dict=dict()
for race, region_ in split_race_region.items():
    for split,regions in region_.items():
        for region in regions:
            result_dict[tuple(region)]=split

df=pd.read_csv("/data2/wangxuedong/mhc_test_data/similarity_region_pythonfiles/csvfiles_95/csv3_top100_with_corr_frompy.csv")
df_split= df[df['isMHC'] == 0]



def getwhichsplit(resultdict,position,chr):
    allsplit=[]
    for digitregion,split_ in resultdict.items():
        if position>=digitregion[0] and position<=digitregion[1]:
            s=split_.split("/")[0].split("_")[-1]
            print(s)
            if chr==s:
                allsplit.append(split_)
        else:
            continue
    return list(set(allsplit))

positions=df_split["Position of SNP"].values.tolist()
chrs=df_split["Position of gene"].values.tolist()
allrecords=[]
for i in range(len(chrs)):
    arr=[]
    ch=chrs[i].split(":")[0]
    arr.append(ch)
    arr.append(positions[i])
    allrecords.append(arr)
allrecords = list(set(tuple(x) for x in allrecords))
allrecords = [list(x) for x in allrecords]


def extract_numbers(input_string):
    return re.findall(r'\d+', input_string)

outputs=[]

for arr in allrecords:
    whichsplit=getwhichsplit(result_dict,arr[1],arr[0])
    if len(whichsplit)==1:
        filepath="/data2/wangxuedong/mhc_test_data/"+whichsplit[0].split("/")[0]+"/"+whichsplit[0].split("/")[1]+".vcf.gz"
        print(filepath)
        vcf_in=vcf(filepath)
        for rec in vcf_in.fetch():
            if rec.pos==arr[1]:
                onerecord=["chromosome"]
                onerecord+=[extract_numbers(arr[0])[0],rec.pos,rec.ref,rec.alts[0],1]
                outputs.append(onerecord)
    else:
        print("with2split",arr)

def generatedf(columns,allrecords):
    dictforDF=dict()
    for i in range(len(columns)):
        midarr=[]
        for record in allrecords:           
            midarr.append(record[i])
        dictforDF[columns[i]]=midarr
    dataframe(dictforDF).to_csv("/data2/wangxuedong/mhc_test_data/similarity_region_pythonfiles/rsform_95/split_part0_rs.csv")
    return dataframe(dictforDF)


df2=generatedf(columns=["chrom","chr","pos","ref","alts[0]","not know"],allrecords=outputs)