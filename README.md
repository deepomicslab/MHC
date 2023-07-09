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

#### Csvfiles and Figures:(disassortative mating regions)

##### Figures:（‘disassortative_3_figs_python_code’ folder）

- ***fig1_except_mhc.ipynb":*** For the pic1 about disassortative mating non-mhc regions
- ***fig1_mhc.ipynb:*** For the pic1 about disassortative mating mhc regions
- ***fig2.ipynb:*** For the pic2 in disassortative mating regions
- ***fig3.ipynb:*** For the pic3 in disassortative mating regions
- ***fig4.ipynb:*** For the pic4 in disassortative mating regions

##### Csvfiles :('disassortative_4_csv_python_code' folder)

- ***csv1_mhc.ipynb:***  For TableS1 in disassortative mating regions
- ***csv2_part0,ipynb, csv2_part1.py, csv2_part2.ipynb:*** For TableS2 in disassortative mating regions
- ***csv3_part0.py, csv3_part1.ipynb, csv3_part2.ipynb***: For TableS3 in disassortative mating regions

- ***form_main_body.ipynb:*** For the form in the main body



### 'disassortative_5_Rs_form' folder:

Code in this folder is used to generate the form that can get #rs value for the disassortative mating regions position.

#### 1.'mhc_rs' folder:

- ***'getrs_mhc_com.py'***: get rs form for disassortative mating mhc regions position

#### 2.'split_rs' folder:

- ***'getrs.py'***: get rs form for disassortative mating non-mhc regions position 

  

### 'negative-assortative_1' folder:

This folder is about 1) finding regions, 2) figures python code 3) csvfiles python code in negative-assortative mating regions.

#### 1.'getP' folder:

- ***'All_P.py':***  Get P according to the formula for 'all'.
- ***'Race_P.py':*** Get P according to the formula for each race.

#### 2. ‘getRegion’ folder:

- ***'get_region_for_all.ipynb':***get negative assortative mating regions for 'all'
- ***'get_region_for_race.ipynb':*** get negative assortative mating regions for each race

#### 3. 'figures' folder:

- ***'pic1_1.ipynb'***: For the pic1 about negative-assortative mating regions 
- ***'pic2.ipynb'***: For the pic2 about negative-assortative mating regions and Average Probability 
- ***'pic3.ipynb'***: For the pic3 about genes included in negative-assortative mating regions 
- ***'pic4.ipynb'***: For the pic4 about comparison between negative-assortative mating non-MHC regions  and negative-assortative mating MHC regions

#### 4.'csvfiles' folder:

- ***'csv1.ipynb'***: TableS1 for negative-assortative mating regions 
- ***'csv2_first_part.ipynb', 'csv2_second_part.py'***: TableS2 for negative-assortative mating regions 
- ***'csv3.py'***: TableS3 for negative-assortative mating regions



### 'negative_assortative_2_RS_form' folder:

Code in this folder is used to generate the form that can get #rs value for the negative assortative mating regions positions.

#### 1. 'mhc_rs' folder:

- ***'rs_mhc_part.py':***  get rs form for negative-assortative mating mhc regions

#### 2. 'split_rs' folder:

- ***'getrs.py'***: get rs form for negative-assortative mating non-mhc regions

### 'Positive_assortative_1' folder:

This folder is about 1) finding regions, 2) figures python code 3) csvfiles python code in positive-assortative mating regions.

#### 1. 'figures' folder:

- ***'pic1.ipynb'***: pic1 for Positive-assortative mating MHC regions
- ***'pic1_part2.ipynb'***: pic1 for Positive-assortative mating non-MHC regions

### 2. 'csvfiles' folder:

- ***'csv1.ipynb'***: TableS1 for Positive-assortative mating regions
- ***'csv2.ipynb', 'csv2_second_part.py'***: TableS2 for for Positive-assortative mating regions
- ***'csv3.py'***: TableS3 for Positive-assortative mating regions



### 'Positive_assortative_2_Rs_form' folder:

Code in this folder is used to generate the form that can get #rs value for the positive-assortative mating regions positions.

- ***'rs.py'***: get rs form for positive-assortative mating regions position 

### 'TableS4' folder:

#### 1.'dis_vs_neg_repeat_region' folder: 

- ***'dis_vs_neg.ipynb':***Find the repeat regions for mhc and non-mhc between disassortative mating regions and negative-assortative mating regions.

- ***'mhc.ipynb'***: Find the genes in the above repeat region in the mhc part.
- ***'non-mhc.ipynb'***: Find the genes in the above repeat region in the non-mhc part.

#### 2. 'other_repeat_region' folder:

- ***'dis_vs_pos.ipynb'***:  Find the repeat region between dissortative mating complementary region and positive assortative mating region
- ***'neg_vs_pos.ipynb'***: Find the repeat region between negative assortative mating region and positive assortative mating region



 















