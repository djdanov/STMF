# STMF
Short-term mortality fluctuations  

To run the scripts you need the following data:  

1. STMF input file - see examples in the STMF folder. 
The most recent version available on the HMD website (https://www.mortality.org/Public/STMF/Inputs/STMFinput.zip)  

2. HMD life tables and HMD Lexis files - see examples in the HMD folder. Data for all countries are available at https://www.mortality.org/hmd/zip/all_hmd/hmd_statistics.zip
Please note that life tables should presented as csv files with first first six columns as follows: Year, Year, Age, Population exposures, Deaths, death rates.  

To run the script:  

1. Check that folders \output and \checks exist  

2. Source expForecast.R and stmf.R  

3. Run stmf(country_code) to get a coutry-specific STMF file (check the folder \output). Here the country_code is the HMD country code used in the STMF.
For example, stmf("DNK"). The whole list of codesis available here: https://www.mortality.org/Public/STMF/Outputs/stmf.xlsx  

4. Run  
a = stmf(country_code)  
stmf.graphs(a, country_code)  
to get data quality checks. Result will be saved in country_codechecks.pdf    


