************************************************************
Workflows for measuring epigenetic age of mouse RRBS samples
************************************************************

Author: Csaba Kerepesi, Brigham and Women's Hospital and Harvard Medical School, Boston, MA, USA
Info/request: ckerepesi@bwh.harvard.edu

Requirements (versions used by us in parenthesis): 

   -- Linux (Ubuntu 16.04.6 LTS) 
   -- Python 3.7.4 
   -- Python packages: pandas (0.25.1), numpy (1.17.2)

Usage:

  -- Running command example: python AppMouseGenomicClocks.py FT.csv Meta.csv
  -- FT.csv should be a transposed feature table (rows: CpG sites, columns: samples)
  -- Meta.csv should be a ',' separated metadata file (rows: samples, columns: attributes). The list of samples should be the same as in FT.csv
  -- Missing values are filled by the mean of the mean values of the clock sites of the samples. It may not practical in some cases. In that case, you can also use pre-filled feature tables.
  -- The script expects methylation values between 0 and 1 and provide preidcted age in months. 
  
Please download and output the following files to the 'ClockData' directory:

    -- elife-40675-supp3-v2.xlsx (Meer et al. 2018. Elife, Supplementary file 3, https://elifesciences.org/articles/40675/figures#supp3)
    -- TrainingData_Babraham_Reizel_Cannon.txt (Stubbs et al. 2017. Genome Biology, https://github.com/EpigenomeClock/MouseEpigeneticClock/blob/master/TrainingMatrix/TrainingData_Babraham_Reizel_Cannon.txt      
    -- Supplemental_Table_S2.xlsx (Wang and Lemos 2019. Genome Research, Supplementary Table S2 https://genome.cshlp.org/content/suppl/2019/02/07/gr.241745.118.DC1/Supplemental_Table_S2.xlsx)
