# Information
- main_script.py python script filters 2 replicates indepandently by hardcoded changeable parameters (FREQ, BASECOUNT, COVERAGE)

# Input for script
- Script will ask where are you running the script for path variables

# Output:
# Results directory organized by Patient name and time point named after script running date and hour:
- replicate_1/2_freq_base_count_coverage.csv - containing all data per replicate after filtering minimum frequency, minimum Base Count and minimum Coverged value
- merged_freq_base_count_coverage.csv - containing all data for join replicates (only common mutations) after filtering minimum frequency, minimum Base Count and minimum Coverged value
- mut_freq_1/2_freq_base_count_coverage - containing mutation name and frequency for replicate1/2
- usecase_freq_base_count_coverage - containing mutation name and inforamtion about use case and CriticalDelta
- Output3 - containing mutation name and frequency for merged replicates after more filtering on both replicates