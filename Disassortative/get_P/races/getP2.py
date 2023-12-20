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
import sys

fileName = '/data2/wangxuedong/mhc_test_data/1000_population.tsv'
lines = []
categories = []
samplenames=[]

with open(fileName, 'r') as f:
    text = f.read()

lines = text.split('\n')
for l in lines:
    samplenames.append(l.split('\t')[0])
    categories.append(l.split('\t')[-1])
dic_population_category=dict()


for i in range(len(samplenames)):
    dic_population_category[samplenames[i]]=categories[i]

def calPA_and_PB_PAB_P(dic):#AA aa Aa Aa     A:4  a:4 
    result=dict()
    for race,GTs in dic.items():
        count_01=0
        count_10=0
        count_0=0
        count_1=0
        list1=[x[0] for x in GTs]
        list2=[x[1] for x in GTs]
        for m in list1:
            if m==0:
                count_0+=1
            else:
                count_1+=1
        for n in list2:
            if n==0:
                count_0+=1
            else:
                count_1+=1
        p0=(count_0+4)/(2*len(GTs)+8)
        p1=(count_1+4)/(2*len(GTs)+8)
        # p0=0.2+0.8*p0
        # p1=0.2+0.8*p1

        for GT in GTs:
            if(operator.eq(GT,(0,1))):
                count_01+=1
            if(operator.eq(GT,(1,0))):
                count_10+=1
                
        num=(count_01+count_10+2)/(len(GTs)+4)
        # if num==0:
        #     cal=0
        # else:
        cal=num*log2(num/(p0*p1))
        result[race]=cal
    return result

def generatedfandcsv(columns,dic,filename):#filename is split_chr1_xaa_NEU.csv

    positions=list(dic.keys())
    df=pd.DataFrame(index=positions,columns=columns)
    for i in range(len(positions)):
        a=dic[positions[i]]
        for j in range(len(columns)):
            df.iloc[i][j]=a[columns[j]]
    #print("df",df)
    df.to_csv("/data2/wangxuedong/mhc_test_data/first_review/race/"+filename+".csv")
    return df



allcategories=sorted(list(set(categories)))
def runall(filepath):

    
    vcf_in=vcf(filepath)
    dic_GT_different_races=dict()
    resultdic=dict()
    samplelist=list((vcf_in.header.samples))
    for rec in vcf_in.fetch():      
        for i in range(len(allcategories)):
            dic_GT_different_races[allcategories[i]]=[]                 
        for samplename in samplelist:
            if len(rec.samples[samplename]['GT'])==2:
                dic_GT_different_races[dic_population_category[samplename]].append(rec.samples[samplename]['GT'])
        #print("dic_GT_different_Race",dic_GT_different_races)
        middle_dict=calPA_and_PB_PAB_P(dic=dic_GT_different_races)
        #print("middle dic",middle_dict)
        resultdic[rec.pos]=middle_dict
    print("resultdic",resultdic)
    filenameprefix=filepath.split(".")[0].split("/")[-2]+"_"+filepath.split(".")[0].split("/")[-1]
    generatedfandcsv(columns=allcategories,dic=resultdic,filename=filenameprefix)

#generate all path
foldernames=[]
for i in range(1,23):
    foldernames.append('split_chr'+str(i))

patharray=[]
for foldername in foldernames:
    for curDir, dirs, files in os.walk(top="/data2/wangxuedong/mhc_test_data/"+foldername+"/"):
        for file in files:
            if file.endswith(".vcf.gz"):
                path=os.path.join(curDir,file)
                patharray.append(path)

sysarg=sys.argv[1]
index=int(sysarg)-1

runall(patharray[index])