# filter.py python script filters 2 replicates indepandently by user parameters of minimum frequency, minimum Base Count and minimum Coverged value
# Input for script
- Folder name containing 2 .tsv files named replicate1 and replicate 2
# Output: Results directory containing:
-replicate1/2.csv - containing all data per replicate after filtering minimum frequency, minimum Base Count and minimum Coverged value
-merged.csv - containing all data for join replicates (only common mutations) after filtering minimum frequency, minimum Base Count and minimum Coverged value
-Output1 - containing mutation name and frequency for replicate1
-Output2 - containing mutation name and frequency for replicate2
-Output3 - containing mutation name and frequency for merged replicates after more filtering on both replicates. Also contains if there is a critical delta between 2 replicates

# filtered_2_BN_input.py python script recives 2 filtered samples of two time points and create exact and approx format input files for Bottle Neck algorithem
# Input for script
- Folder name containing 2 .tsv files named replicate1 and replicate2
# Output: Results directory containing:
-replicate1/2.csv - containing all data per replicate after filtering minimum frequency, minimum Base Count and minimum Coverged value
-merged.csv - containing all data for join replicates (only common mutations) after filtering minimum frequency, minimum Base Count and minimum Coverged value
-Output1 - containing mutation name and frequency for replicate1
-Output2 - containing mutation name and frequency for replicate2
-Output3 - containing mutation name and frequency for merged replicates after more filtering on both replicates. Also contains if there is a critical delta between 2 replicates