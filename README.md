#   					MHC

### Disassortative mating Regions:

#### For All: ('disassortative_1_all_find_region' folder)

##### The following steps are used to get complementary regions for 'All'.

- ***get_P_region.ipynb***: This file is used to get the probability from the formula

  input: vcf files 

  output: csv files with probability

- ***get_region_first_part.ipynb***: 

  Take the average of every 1000 P values and save it in the file, named with 'csvfile_bin_1000_AveragewithoutX' folder.   

  Then, use the file in 'csvfile_bin_1000_AveragewithoutX' to get the pvalue, if the p-value is greater than 0.95, mark it as 1, otherwise, mark it as 0. Put the above files into 'csvfiles_compare_with_0.01pvaluewithoutX' folder

- ***get_region_second_part.ipynb***:

  From the 'csvfiles_compare_with_0.01pvaluewithoutX' folder, get the number of consecutive '1's in each file, and save the position of the starting '1' and the corresponding length of each group in 'csvfile_start_lengthwithoutX' folder. 

  Get the size of threshold by analyzing the length of continuous '1'.

- ***get_region_third_part.ipynb***:

  Obtain disassortative mating regions according to threshold.

#### For Race :('disassortative_2_race_find_region' folder)

##### The following steps are used to get complementary regions for different races.

- ***get_P_race.py***: This file is used to get the probability from the formula

  input: vcf files

  output: csv files with probability

- ***get_position_race.ipynb***

  Get disassortative mating regions for different race, whose method is the same as 'All'.

### Csvfiles and Figures:(disassortative mating regions)

#### Figures:（‘disassortative_3_figs_python_code’ folder）

- ***'FIGURE 1_MHC.ipynb':*** For the FIGURE 1 about disassortative mating mhc regions
- ***'FIGURE 1_NONMHC.ipynb':*** For the FIGURE 1 about disassortative mating non-mhc regions
- ***'FIGURE 3.ipynb':*** For the FIGURE 3
- ***'FIGURE 4.ipynb':*** For the FIGURE 4
- ***'FIGURE 5.ipynb':*** For the FIGURE 5

#### Csvfiles :('disassortative_4_csv_python_code' folder)

- ***'Table S1.ipynb':***  For Table S1 in supplementary material
- ***'Table S3_part0.ipynb', 'Table S3_part1.ipynb':*** For TableS3 in supplementary material
- ***'Table S4.ipynb':*** For TableS4 in supplementary material



### 'disassortative_5_Rs_form' folder:

Code in this folder is used to generate the form that can get #rs value for the disassortative mating regions position.

#### 1.'mhc_rs' folder:

- ***'getrs_mhc_com.py':*** get rs form for disassortative mating mhc regions position

#### 2.'split_rs' folder:

- ***'getrs.py':*** get rs form for disassortative mating non-mhc regions position 

  


### 'Positive_assortative_1' folder:

This folder is about 1) finding regions, 2) figures python code 3) csvfiles python code in positive-assortative mating regions.

#### 1). 'figures' folder:

- ***'Figure 2.ipynb'***: For FIGURE 2 


#### 2). 'csvfiles' folder:

- ***'Table S2.ipynb'***: For Table S2 in supplementary material
- ***'Table S5_part0.ipynb', 'Table S5_part1.py'***: For Table S5 in supplementary material
- ***'Table S6.py'***: For Table S6 in supplementary material



#### 'Positive_assortative_2_Rs_form' folder:

Code in this folder is used to generate the form that can get #rs value for the positive-assortative mating regions positions.

- ***'rs.py'***: get rs form for positive-assortative mating regions position 



 















