from math import floor
from matplotlib import pyplot as plt
import pandas as pd
import numpy as np
from scipy import stats
import glob, os

#nate's beat me to a lot, seems to use the Mining_IR script as a method in the PU script

#next idea is make dataframe with just the wavenumber of the first csv then specifically add the second column from each csv starting from the top

WN_Normal = 1310 #excel 1311 -> WN 2926
WN_Normal2 = 960 #excel 961 --> WN 2250, use to check if ratio between normalizing WN changes
WN_SiH = 917 #excel 918 -> WN 2168
WN_SiOSi = 336 #excel 337 --> WN 1047
WN_SiC = 204 #excel 205 --> WN 792

WN_List = [WN_Normal, WN_Normal2, WN_SiH, WN_SiOSi, WN_SiC]

path = r'C:\Users\taekw\Desktop\PythonScripts\IRPeakExtract\CSVs\220202_cbPDMS_Cure_Extent_Tracking'
os.chdir(path)

#! import wavenumber and intensity from csv as dataframe
# 1) make dataframe with first csv in directory
# 2) iteratively append temporary dataframes with just the second column from each proceeding CSV using glob to rename the columns

all_filenames = [i for i in glob.glob('*.{}.csv')]

df_IR_total = pd.concat(map(pd.read_csv, glob.glob(os.path.join(path, '*.csv'))), axis=1, ignore_index = True)
df_IR_total.drop(df_IR_total.columns[[2, 4, 6, 8, 10, 12, 14, 16, 18]], axis = 1, inplace = True)
df_IR_peaks = df_IR_total.iloc[WN_List]

df_IR_peaks = df_IR_peaks.rename(columns={0:'Wavenumber', 3: '0A_1'})
df_IR_peaks = df_IR_peaks.rename(index={1310:'Normal (2926)', 960: 'Normal 2 (2250)', 917 : 'SiH (2168)', 336 : 'SiOSi (1047)', 204 : 'SiC (792)'})

print(all_filenames)
print(df_IR_peaks)
#print(df_IR_total)

#! find normalization intensity at (wavenumber)

#! append intensities from other CSVs stepwise, extracting column header from filename

#! write new dataframe with extracted intensities at 2950-2960 cm-1 (CH3's C-H stretch) and 1020-1074 cm-1 (Si-O-Si asymmetric and symmetric stretches) (Joe did 1597 and 2162)

#! plot extracted intensities with trend line and inset graph of the normalized IR plots

#! write both new dataframes into a new CSV