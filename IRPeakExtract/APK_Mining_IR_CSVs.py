#IMPORT dependencies
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
import glob

#GET all the csv files in the working folder with glob (later: change to opening a 'select directory' prompt or drag/drop)
path = r'C:\Users\taekw\Desktop\PythonScripts\IRPeakExtract\CSVs\220202_cbPDMS_Cure_Extent_Tracking'
csv_files = glob.glob(os.path.join(path, "*.csv"))

WN_Normal = 5316 #line 5316 -> WN 2926
WN_SiH = 3669 #line 3669 -> WN 2168
WN_SiO = 420 #idk what WN or line corresponds to SiO yet


#Goal 1: LOOP through CSVs one by one, finding max number and appending it as a column in a new CSV
#Goal 2: find Signal_WN_normal, divide signal_WN_interest, use worked up numbers instead
#Goal 3: take WN's as inputs for normal and interest
#Goal 4: take arbitrary number of WN's as interest inputs 

#Create a new CSV to add data to (Later: prompt for name of new CSV)

# loop over the list of csv files
for f in csv_files:
    #WRITE the filename as a column header (currently just prints the filename)
    print('File Name:', f.split("\\")[-1])

    #WRITE the CSV into temp dataframe
    df = pd.read_csv(f)
    
    #EXTRACT the signal for a specified wavenumber and append it in new CSV (normal ~ line_5316 for WN_2926, Si-H ~ line_3669 for WN_2168)
    df_all = pd.concat([df_first_shift, df_second_shift, df_third_shift])
    print(df_all) 

    # print the content
    print('Content:')
    display(df)
    print()

pivot = df_all.groupby(['Shift']).mean()
shift_productivity = pivot.loc[:,"Production Run Time (Min)":"Products Produced (Units)"]

#TROUBLESHOOTING GRAPH
#shift_productivity.plot(kind='bar')
#plt.show()

#OUTPUT new CSV
df_all.to_csv("output.csv")