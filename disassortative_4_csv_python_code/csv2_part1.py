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


df_csv2=pd.read_csv("/data2/wangxuedong/mhc_test_data/python/newrequirements/filtered_csv2.csv").values.tolist()


all_mean_arr=-0.019077723390066942
all_std_arr=0.006128439002438663
mean_dict_different_race={'ACB': -0.014011151662888553,
 'ASW': -0.012693412503283,
 'BEB': -0.0038608773522230016,
 'CDX': -0.000870975201046489,
 'CEU': -0.002626550761399808,
 'CHB': -0.000840617782511786,
 'CHS': -0.0016030210248706835,
 'CLM': -0.006497201413442441,
 'ESN': -0.012407888605124987,
 'FIN': -0.0013888123211964092,
 'GBR': -0.0016421960521901954,
 'GIH': -0.0029655045615120147,
 'GWD': -0.013383054011719266,
 'IBS': -0.0036274568161741316,
 'ITU': -0.003697573961398361,
 'JPT': -0.0004422965646677968,
 'KHV': -0.0014824412579579188,
 'LWK': -0.013214685221562708,
 'MSL': -0.012902619299438558,
 'MXL': -0.004064132377929528,
 'PEL': -0.0032051883341473507,
 'PJL': -0.004661546996727815,
 'PUR': -0.008252100717713692,
 'STU': -0.003961436794723387,
 'TSI': -0.002358177545570631,
 'YRI': -0.012978184004079407}
std_dict_different_race={'ACB': 0.007679014864719199,
 'ASW': 0.008041527656473924,
 'BEB': 0.008061181937627181,
 'CDX': 0.008574772155886829,
 'CEU': 0.007805712926154353,
 'CHB': 0.00858396197930932,
 'CHS': 0.008210661522162982,
 'CLM': 0.007855389049385354,
 'ESN': 0.007528169549036397,
 'FIN': 0.008233373830804186,
 'GBR': 0.008280482586380011,
 'GIH': 0.008318428884896655,
 'GWD': 0.00741808675422558,
 'IBS': 0.007946168606506706,
 'ITU': 0.007944085877932789,
 'JPT': 0.008520600208015637,
 'KHV': 0.008457471872578115,
 'LWK': 0.00775427047797345,
 'MSL': 0.00788405338804666,
 'MXL': 0.008244211463289736,
 'PEL': 0.008220014590454098,
 'PJL': 0.007728364636381623,
 'PUR': 0.007789028974214059,
 'STU': 0.007989081367465312,
 'TSI': 0.008171487042712589,
 'YRI': 0.007574843363674686}

allrecords=[]
#in csv2
for record in df_csv2:
    print("record",record)
    print("record6",record)
    if type(record[6])!=float:
        if record[1]=="All":
            if record[2]==0:#no mhc csvfiles                        
                position_of_gene=record[6]
                chr=position_of_gene.split(":")[0]
                start_position=int(position_of_gene.split(":")[1].split("-")[0])

                end_position=int(position_of_gene.split(":")[1].split("-")[1])
                filenames=os.listdir("/data2/wangxuedong/mhc_test_data/csvfiles/")
                for filename in filenames:
                    if filename.startswith("pos_prob_split_"+chr+"_"):
                        df_split=pd.read_csv("/data2/wangxuedong/mhc_test_data/csvfiles/"+filename)  

                        allrigion_P_top_1000=df_split[(df_split["Unnamed: 0"]>=start_position) & (df_split["Unnamed: 0"]<=end_position)]["P"].to_list()
                        if allrigion_P_top_1000:
                            mean_100_P=np.mean(allrigion_P_top_1000)
                            std_100_P=np.std(allrigion_P_top_1000, ddof=1)
                            pvalue=round((1-stats.norm.cdf(mean_100_P, all_mean_arr, all_std_arr)),15)
                            onerecord=[]
                            onerecord+=record[1:9]
                            onerecord+=[mean_100_P,std_100_P,pvalue]
                            allrecords.append(onerecord)
                                

                        else:
                            continue

            else:
                position_of_gene=record[6]
                start_position=int(position_of_gene.split(":")[1].split("-")[0])
                end_position=int(position_of_gene.split(":")[1].split("-")[1])
                df_mhc=pd.read_csv("/data2/wangxuedong/mhc_test_data/csvfiles/posmhc.csv")
                allrigion_P_top_1000=df_mhc[(df_mhc["Unnamed: 0"]>=start_position) & (df_mhc["Unnamed: 0"]<=end_position)]["P"].to_list()
                if allrigion_P_top_1000:
                    mean_100_P=np.mean(allrigion_P_top_1000)
                    std_100_P=np.std(allrigion_P_top_1000, ddof=1)
                    pvalue=round((1-stats.norm.cdf(mean_100_P, all_mean_arr, all_std_arr)),15)
                    onerecord=[]
                    onerecord+=record[1:9]
                    onerecord+=[mean_100_P,std_100_P,pvalue]
                    allrecords.append(onerecord)
                else:
                    continue

        else:#race!=all
            if record[2]==0:#no mhc
                position_of_gene=record[6]
                chr=position_of_gene.split(":")[0]
                start_position=int(position_of_gene.split(":")[1].split("-")[0])

                end_position=int(position_of_gene.split(":")[1].split("-")[1])
                filenames=os.listdir("/data2/wangxuedong/mhc_test_data/race/")
                for filename in filenames:
                    if filename.startswith("split_"+chr+"_"):
                        df_split=pd.read_csv("/data2/wangxuedong/mhc_test_data/race/"+filename)  
                    
                        allrigion_P_top_1000=df_split[(df_split["Unnamed: 0"]>=start_position) & (df_split["Unnamed: 0"]<=end_position)][record[1]].to_list()
                        if allrigion_P_top_1000:
                            mean_100_P=np.mean(allrigion_P_top_1000)
                            std_100_P=np.std(allrigion_P_top_1000, ddof=1)
                            pvalue=round((1-stats.norm.cdf(mean_100_P,mean_dict_different_race[record[1]],std_dict_different_race[record[1]])),15)
                                
                            onerecord=[]
                            onerecord+=record[1:9]
                            onerecord+=[mean_100_P,std_100_P,pvalue]
                            allrecords.append(onerecord)
                        else:
                            continue
                    else:
                        continue
            else:#ismhc different race
                position_of_gene=record[6]
                start_position=int(position_of_gene.split(":")[1].split("-")[0])
                end_position=int(position_of_gene.split(":")[1].split("-")[1])
                df_mhc=pd.read_csv("/data2/wangxuedong/mhc_test_data/race/mhc.csv")
                
                allrigion_P_top_1000=df_mhc[(df_mhc["Unnamed: 0"]>=start_position) & (df_mhc["Unnamed: 0"]<=end_position)][record[1]].to_list()
                if allrigion_P_top_1000:                
                    mean_100_P=np.mean(allrigion_P_top_1000)
                    std_100_P=np.std(allrigion_P_top_1000, ddof=1)
                    pvalue=round((1-stats.norm.cdf(mean_100_P,mean_dict_different_race[record[1]],std_dict_different_race[record[1]])),15)              
                    onerecord=[]
                    onerecord+=record[1:9]
                    onerecord+=[mean_100_P,std_100_P,pvalue]
                    allrecords.append(onerecord)
                else:
                    continue
    else:
        print("record",record)
        print("record6",record)
        onerecord=[]
        onerecord+=record[1:9]
        onerecord+=["null","null","null"]
        allrecords.append(onerecord)


#generatedf
def generatedf(columns,allrecords):
    dictforDF=dict()
    for i in range(len(columns)):
        midarr=[]
        for record in allrecords:           
            midarr.append(record[i])
        dictforDF[columns[i]]=midarr
    dataframe(dictforDF).to_csv("csv2_new_05_11_without_correct_pvalue.csv")
    return dataframe(dictforDF)


dfnew=generatedf(columns=["Race","isMHC","Position of disassortative mating region","Gene ID","Gene name","Position of gene","is_complement","Gene description","Average of P","Std of P","P-value"],allrecords=allrecords)


pvalue_list=dfnew["P-value"].tolist()
corr_pvals = multipletests(pvalue_list, method='fdr_bh')[1]
dfnew = dfnew.assign(corr_pvals=corr_pvals)
dfnew.to_csv("/data2/wangxuedong/mhc_test_data/python/newrequirements/csv2_again_xu_05_11.csv")