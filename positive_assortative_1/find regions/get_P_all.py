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
import sys

class CalProb:
    def __init__(self,vcf_in):
        self.vcf_in=vcf_in
        self.result=dict()

    def calculateP(self):
        for rec in self.vcf_in.fetch():
            count_0=0
            count_00=0
            count_1=0
            count_11=0
            GTs=[]
            samplelist=list((self.vcf_in.header.samples))
            # print(samplelist)
            for samplename in samplelist:
                t1=rec.samples[samplename]['GT']
                if len(t1)==2:
                    GTs.append(t1)
                if(operator.eq(t1,(0,0))):
                    count_00+=1
                if(operator.eq(t1,(1,1))):
                    count_11+=1
                    
            #print(GTs)
            list1=[x[0] for x in GTs]
            list2=[x[1] for x in GTs]
            #p0 and p1
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
            p0=count_0/(2*len(samplelist))
            p1=count_1/(2*len(samplelist))
            p0=0.2+0.8*p0
            p1=0.2+0.8*p1
            arr=[]#arr[p0,p1,p00,p11,p]
            arr.append(p0)
            arr.append(p1)
            #p00
            p00=count_00/(len(samplelist))
            p11=count_11/(len(samplelist))
            arr.append(p00)
            arr.append(p11)
            if p00==0:
                p00=0
            else:
                p00=p00*log2(p00/(p0*p0))
            if p11==0:
                p11=0
            else:
                p11=p11*log2(p11/(p1*p1))
            finalp=p00+p11
            arr.append(finalp)
            self.result[rec.pos]=arr

    def getresult(self):
        self.calculateP()
        return self.result
    
    def df_pos_prob(self,pos_prob_csvname):
        positions=list(self.result.keys())
        column_names=['P0','P1','P00','P11','P']
        df=dataframe(index=positions,columns=column_names)
        for i in range(len(positions)):
            a=self.result[positions[i]]
            df.iloc[i][0]=a[0]
            df.iloc[i][1]=a[1]
            df.iloc[i][2]=a[2]
            df.iloc[i][3]=a[3]
            df.iloc[i][4]=a[4]
        df.to_csv("/data2/wangxuedong/mhc_test_data/similar_region_csvfiles/"+pos_prob_csvname+".csv")
        return df
#/data2/wangxuedong/mhc_test_data/split_chr1/xaa.vcf.gz
def runall(filepath):
    fileprefix=filepath.split("/")[5].split(".")[0]
    print("fileprefix",fileprefix)
    foldername=filepath.split("/")[4]
    print("foldername",foldername)
    pos_prob_name="pos_prob_"+foldername+"_"+fileprefix
    
    vcf_in=vcf(filepath)
    testclass=CalProb(vcf_in=vcf_in)
    testclass.calculateP()
    df_pos_prob=testclass.df_pos_prob(pos_prob_csvname=pos_prob_name)


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


sysarg=sys.argv[1]
index=int(sysarg)-1

runall(patharray[index])
