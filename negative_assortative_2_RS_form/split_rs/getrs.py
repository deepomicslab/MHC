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

split_race_region={'ACB': {},
 'ASW': {},
 'BEB': {'split_chr3/xav': [[98046342, 98273693]]},
 'CDX': {'split_chr13/xai': [[55102733, 55443491]]},
 'CEU': {'split_chr17/xaj': [[45670177, 46418024]]},
 'CHB': {},
 'CHS': {'split_chr3/xav': [[98046342, 98273693]]},
 'CLM': {'split_chr17/xaj': [[45836662, 46224960]],
  'split_chr3/xav': [[98046342, 98273693]]},
 'ESN': {'split_chr5/xai': [[34264948, 34562837]]},
 'FIN': {},
 'GBR': {'split_chr1/xbt': [[226978066, 227352188]],
  'split_chr17/xaj': [[45670177, 46418024]],
  'split_chr3/xav': [[98046342, 98273693]]},
 'GIH': {'split_chr3/xav': [[98046342, 98276757]]},
 'GWD': {},
 'IBS': {'split_chr17/xaj': [[45670177, 46418024]],
  'split_chr3/xav': [[98046342, 98273693]]},
 'ITU': {'split_chr17/xak': [[52575391, 52889013]]},
 'JPT': {},
 'KHV': {},
 'LWK': {},
 'MSL': {'split_chr5/xai': [[34264948, 34562837]]},
 'MXL': {'split_chr1/xbt': [[227016487, 227352188]],
  'split_chr3/xav': [[98046342, 98273693]]},
 'PEL': {'split_chr11/xaf': [[23168664, 23564494]],
  'split_chr3/xav': [[98046342, 98276757]]},
 'PJL': {},
 'PUR': {'split_chr17/xaj': [[45836662, 46418024]],
  'split_chr3/xav': [[98046342, 98276757]]},
 'STU': {'split_chr13/xai': [[55102733, 55443491]],
  'split_chr3/xav': [[98046342, 98273693]]},
 'TSI': {'split_chr17/xaj': [[45670177, 46418024]]},
 'YRI': {}}



result_dict=dict()
for race, region_ in split_race_region.items():
    for split,regions in region_.items():
        for region in regions:
            result_dict[tuple(region)]=split

df=pd.read_csv("/data2/wangxuedong/mhc_test_data/similarity_region_pythonfiles/csvfiles/csv3_top100_with_corr.csv")
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
    dataframe(dictforDF).to_csv("/data2/wangxuedong/mhc_test_data/similar_region_Z_RS/form/split_part_rs.csv")
    return dataframe(dictforDF)


df2=generatedf(columns=["chrom","chr","pos","ref","alts[0]","not know"],allrecords=outputs)