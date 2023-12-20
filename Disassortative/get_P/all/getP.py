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

class CalProb:
    def __init__(self,vcf_in):
        self.vcf_in=vcf_in
        self.result=dict()
    
    def calPaAndPb(self):#AA aa Aa Aa     A:4  a:4 
        for rec in self.vcf_in.fetch():
            count_0=0
            count_1=0
            GTs=[]
            samplelist=list((self.vcf_in.header.samples))
            for samplename in samplelist:
                t1=rec.samples[samplename]['GT']
                if len(t1)==2:
                    GTs.append(t1)
                
            #print(GTs)
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
            p0=(count_0+4)/(2*len(samplelist)+8)
            p1=(count_1+4)/(2*len(samplelist)+8)
            # p0=0.2+0.8*p0
            # p1=0.2+0.8*p1
            arr=[]
            arr.append(p0)
            arr.append(p1)
            self.result[rec.pos]=arr
    def calPAB(self):     
        for rec in self.vcf_in.fetch():
            count_01=0
            count_10=0
            samplelist=list((self.vcf_in.header.samples))
            for samplename in samplelist:
                t1=rec.samples[samplename]['GT']
                if(operator.eq(t1,(0,1))):
                    count_01+=1
                if(operator.eq(t1,(1,0))):
                    count_10+=1
            num=(count_01+count_10+2)/(len(samplelist)+4)
            self.result[rec.pos].append(num)

    def calP(self):
        positions=list(self.result.keys())
        for position in positions:
            arr=self.result[position]
            # if arr[2]==0:
            #     cal=0
            # else:
                # if arr[0]*arr[1]==0:
                #     cal=arr[2]*log2(arr[2]/(1e-10))
             
            cal=arr[2]*log2(arr[2]/(arr[0]*arr[1]))
            
            self.result[position].append(cal)
    def calculateP(self):
        self.calPaAndPb()
        self.calPAB()
        self.calP()
        
    def __getResult__(self):
        return self.result
    
    def df_pos_prob(self,pos_prob_csvname):
        positions=list(self.result.keys())
        column_names=['PA','PB','PAB','P']
        df=dataframe(index=positions,columns=column_names)
        for i in range(len(positions)):
            a=self.result[positions[i]]
            df.iloc[i][0]=a[0]
            df.iloc[i][1]=a[1]
            df.iloc[i][2]=a[2]
            df.iloc[i][3]=a[3]
        df.to_csv("/data2/wangxuedong/mhc_test_data/first_review/csvfiles/"+pos_prob_csvname+".csv")
        return df
    
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