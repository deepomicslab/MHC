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

* ***get_P_race.py***: This file is used to get the probability from the formula

  input: vcf files

  output: csv files with probability

* ***get_position_race.ipynb***

  Get complementary region for differenct race, whose method is the same as 'All'.

  

  
