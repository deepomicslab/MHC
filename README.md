#   					MHC

### Complementary Region:

#### For All: ('All_find_region' folder)

##### If you want to run the code, please proceed in the following order.

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

  Obtain complementary regions according to threshold.

#### For Race :

##### If you want to run the code, please proceed in the following order.

- ***get_P_race.py***: This file is used to get the probability from the formula

  input: vcf files

  output: csv files with probability

- ***get_position_race.ipynb***

  Get complementary region for differenct race, whose method is the same as 'All'.

#### Csvfiles and Figures:

##### Csvfiles in complementary region:('csv_python_code')

- ***csv1_mhc.ipynb:***  For TableS1 in complementary region
- ***csv2_part0,ipynb, csv2_part1.py, csv2_part2.ipynb:*** For TableS2 in complementary region
- ***csv3_part0.py, csv3_part1.ipynb, csv3_part2.ipynb***: For TableS3 in complmentary region

- ***form_main_body.ipynb:*** For the form in the main body

##### Figures in complementary region:（‘figs_python_code’）

- ***fig1_except_mhc.ipynb":*** For the pic1 about complementary non-mhc region
- ***fig1_mhc.ipynb:*** For the pic1 about complementary mhc region
- ***fig2.ipynb:*** For the pic2 in complementary region
- ***fig3.ipynb:*** For the pic3 in complementary region
- ***fig4.ipynb:*** For the pic4 in complementary region

### 'complementary_Region_Rs_form' folder:

Code in this folder is used to generate the form that can get #rs value for the dissorsative complementary region position.

#### 1.'mhc_rs' folder:

- ***'getrs_mhc_com.py'***: get rs form for dissorsative complementary mhc region position

#### 2.'split_rs' folder:

- ***'getrs.py'***: get rs form for dissorsative complementary non-mhc region position 

  

### 'Negative-assortative mating region' folder:

#### 1.'getP' folder:

- ***'All_P.py':***  Get P according to the formula for 'all'.
- ***'Race_P.py':*** Get P according to the formula for each race.

#### 2. ‘getRegion’ folder:

- ***'get_region_for_all.ipynb':***get negative similar region for 'all'
- ***'get_region_for_race.ipynb':*** get negative similar region for each race

#### 3. 'figures' folder:

- ***'pic1_1.ipynb'***: For the pic1 about Negative-assortative mating regions 
- ***'pic2.ipynb'***: For the pic2 about Negative-assortative mating regions and Average Probability 
- ***'pic3.ipynb'***: For the pic3 about genes included in Negative-assortative mating regions 
- ***'pic4.ipynb'***: For the pic4 about comparison between Negative-assortative mating non-MHC regions  and Negative-assortative mating MHC regions

#### 4.'Csvfiles' folder:

- ***'csv1.ipynb'***: TableS1 for Negative-assortative mating regions 
- ***'csv2_first_part.ipynb', 'csv2_second_part.py'***: TableS2 for Negative-assortative mating regions 
- ***'csv3.py'***: TableS3 for Negative-assortative mating regions



### 'Negative_assortative_mating_region_RS' folder:

Code in this folder is used to generate the form that can get #rs value for the negative similar region position.

#### 1. 'mhc_rs' folder:

- ***'rs_mhc_part.py':***  get rs form for negative-assortative mating mhc region

#### 2. 'split_rs' folder:

- ***'getrs.py'***: get rs form for negative-assortative mating non-mhc region

### 'Positive similar region' folder:

#### 1. 'figures' folder:

- ***'pic1.ipynb'***: pic1 for Positive-assortative mating MHC region 
- ***'pic1_part2.ipynb'***: pic1 for Positive-assortative mating non-MHC region

### 2. 'csvfiles' folder:

- ***'csv1.ipynb'***: TableS1 for Positive similar regions
- ***'csv2.ipynb', 'csv2_second_part.py'***: TableS2 for for Positive-assortative mating regions
- ***'csv3.py'***: TableS3 for Positive-assortative mating regions



### 'Positive_assortative_mating_region_rs' folder:

Code in this folder is used to generate the form that can get #rs value for the positive assortative mating region position.

- ***'rs.py'***: get rs form for positive similar region position 

### 'TableS4' folder:

#### 1.'repeat_region.ipynb': 

​      Find the repeat region in the mhc part and non-mhc part between 'All' and 'Race'

####2. 'mhc' folder:

- ***'mhc.ipynb'***: find the genes in the repeat region in the mhc part

#### 3. 'split' folder:

- ***'split.ipynb'***: find the genes in the repeat region in the non-mhc part

#### 4. 'others'  folder:

- ***'discom_vs_pos.ipynb'***:  Find the repeat region between dissortative mating complementary region and positive assortative mating region
- ***'neg_vs_pos.ipynb'***: Find the repeat region between negative assortative mating region and positive assortative mating region



 

##### 













