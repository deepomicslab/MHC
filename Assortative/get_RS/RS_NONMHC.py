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

df=pd.read_csv("/data2/wangxuedong/mhc_test_data/first_review/a_rsform/csvs/rs_tables6.csv")

df_split= df[df['isMHC'] == 0]
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
allrecords=sorted(allrecords)

def getwhichsplit(position,chr):
    allsplit=[]
    allsplits_filepaths=[]
    for file in os.listdir('/data2/wangxuedong/mhc_test_data/first_review/similar_region_csvfiles/'):
        if file!="posmhc.csv":
            if chr ==file.split(".")[0].split("_")[3]:
                allsplits_filepaths.append(file)
    allsplits_filepaths=['/data2/wangxuedong/mhc_test_data/first_review/similar_region_csvfiles/'+i for i in allsplits_filepaths]
    #print(allsplits_filepaths)
    if(position>=29720403) &(position<=33129982):
        allsplit.append("mhc")
        return allsplit
    for split_path in allsplits_filepaths:
        df=pd.read_csv(split_path)
        poses=sorted(df["Unnamed: 0"].to_list())
        #print(poses)

        if (position>=poses[0]) & (position<=poses[-1]):
            print('yeah')
            allsplit.append(chr+"/"+split_path.split("/")[-1].split(".")[0].split("_")[-1])
    return allsplit

allrecords=allrecords


def extract_numbers(input_string):
    return re.findall(r'\d+', input_string)

outputs=[]

for arr in allrecords:
    whichsplit=getwhichsplit(arr[1],arr[0])
    if len(whichsplit)==1:
        filepath="/data2/wangxuedong/mhc_test_data/split_"+whichsplit[0].split("/")[0]+"/"+whichsplit[0].split("/")[1]+".vcf.gz"
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
    dataframe(dictforDF).to_csv("/data2/wangxuedong/mhc_test_data/first_review/a_rsform/rs_ass_res/ASS_SPLIT_PART1.csv")
    return dataframe(dictforDF)


df2=generatedf(columns=["chrom","chr","pos","ref","alts[0]","not know"],allrecords=outputs)
