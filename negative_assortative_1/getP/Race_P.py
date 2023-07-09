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

def calPA_and_PB_PAB_P(dic):
    result=dict()
    for race,GTs in dic.items():
        count_00=0
        count_11=0
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
        p0=count_0/(2*len(GTs))
        p1=count_1/(2*len(GTs))
        p0=0.2+0.8*p0
        p1=0.2+0.8*p1

        for GT in GTs:
            if(operator.eq(GT,(0,0))):
                count_00+=1
            if(operator.eq(GT,(1,1))):
                count_11+=1
        p00=count_00/(len(GTs))
        p11=count_11/(len(GTs))

        if p00==0:
            cal0=0
        else:
            cal0=p00*log2(p00/(p0*p0))
        if p11==0:
            cal1=0
        else:
            cal1=p11*log2(p11/(p1*p1))
  
        result[race]=cal0+cal1
    return result

def generatedfandcsv(columns,dic,filename):#filename is split_chr1_xaa_NEU

    positions=list(dic.keys())
    df=pd.DataFrame(index=positions,columns=columns)
    for i in range(len(positions)):
        a=dic[positions[i]]
        for j in range(len(columns)):
            df.iloc[i][j]=a[columns[j]]
    print("df",df)
    df.to_csv("/data2/wangxuedong/mhc_test_data/similar_region_race/"+filename+".csv")
    return df

def check_file_not_exists(filepath):
    folder_path="/data2/wangxuedong/mhc_test_data/similar_region_race"
    file_name=filepath.split("/")[4]+"_"+filepath.split("/")[5].split(".")[0]+".csv"
    file_path = os.path.join(folder_path, file_name)
    if os.path.exists(file_path):
        
        return False
    else:
        return True


allcategories=sorted(list(set(categories)))
def runall(filepath):

    if check_file_not_exists(filepath=filepath):

        vcf_in=vcf(filepath)
        dic_GT_different_races=dict()
        resultdic=dict()
        samplelist=list((vcf_in.header.samples))
        for rec in vcf_in.fetch():      
            for i in range(len(allcategories)):
                dic_GT_different_races[allcategories[i]]=[]                 
            for samplename in samplelist:
                dic_GT_different_races[dic_population_category[samplename]].append(rec.samples[samplename]['GT'])
            #print("dic_GT_different_Race",dic_GT_different_races)
            middle_dict=calPA_and_PB_PAB_P(dic=dic_GT_different_races)
            #print("middle dic",middle_dict)
            resultdic[rec.pos]=middle_dict
        print("resultdic",resultdic)
        filenameprefix=filepath.split(".")[0].split("/")[-2]+"_"+filepath.split(".")[0].split("/")[-1]
        generatedfandcsv(columns=allcategories,dic=resultdic,filename=filenameprefix)
    else:
        print("file already exists")

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


patharray = sorted(patharray)
patharray=patharray


sysarg=sys.argv[1]
index=int(sysarg)-1

runall(patharray[index])