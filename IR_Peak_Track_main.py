from math import floor
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from scipy import stats
from scipy.integrate import simps
import glob, os

readpath = r'C:\Users\taekw\Desktop\1_PythonScripts\IRPeakExtract\CSVs\220218_cbPDMS_50ppm_LaserTracking_10A'
writepath = r'C:\Users\taekw\Desktop\1_PythonScripts\IRPeakExtract\CSVs\Output'
output = '220317_5e7 cbPDMS_laser 10A oven 100C_abstest.csv'

os.chdir(readpath) #change working directory to folder with CSVs of interest
filelist = sorted(glob.glob('*.csv')) #make list of names of csvs from readpath directory  

#band = 'SiC'
WN_low = 2100 #2100
WN_high = 2200 #2200

colorlist = ['#ffbaba', '#fc6d6d', '#f10f0f', '#ca1212', '#800c0c']
colorlength = len(colorlist)

df_IR_total = pd.read_csv(filelist[0], names = ['wavenumber','temp'])
df_IR_total = df_IR_total.drop(['temp'], axis = 1)
WN_array = np.array(df_IR_total['wavenumber'])
df_IR_total.set_index('wavenumber', append = False, inplace = True)
# can check here later if temp column is abs or trans
df_IR_total_abs = df_IR_total

def make_integration_arrays(lower_bound, upper_bound, wave_array):
    index = 0
    for item in wave_array:
        if item < lower_bound:
            index += 1
    index2 = 0
    for item in wave_array:
        if item < upper_bound:
            index2 += 1
    return index, index2

index_low, index_high = make_integration_arrays(WN_low, WN_high, WN_array)

areaarray = []
namearray = []
ax = plt.gca()
colorcount = 0

for file in filelist:     #append signal columns from each csv in directory to the intial dataframe stepwise in loop
    columnname = file[0:-4] #extract column name (slicing off '.csv' from filename)
    df_temp = pd.read_csv(file, names = ['wavenumber', columnname]) #temp dataframe with csv
    #df_temp.set_index('wavenumber', append = False, inplace = True)
    #df_temp.plot(kind = 'scatter', x = 'wavenumber', y = columnname, s = 2, c = colorlist[colorcount], ax = ax)
    df_IR_total = pd.concat([df_IR_total, df_temp[columnname]], axis = 1) #add current signals column to total dataframe
    
    # transmittance-to-absorbance conversion
    df_temp[columnname] /= 100
    df_temp[columnname] += 1e-10
    df_temp[columnname] = np.log10(df_temp[columnname]) * -1
    #df_temp.plot(kind = 'scatter', x = 'wavenumber', y = columnname, s = 2, c = colorlist[colorcount], ax = ax)    
    # /transmittance-to-absorbance conversion
    
    # peak integration & baseline subtraction
    abs_array = np.array(df_temp[columnname])
    area = simps(abs_array[index_low:index_high], WN_array[index_low:index_high])
    m = (abs_array[index_low] - abs_array[index_high])/(WN_array[index_low] - WN_array[index_high])
    b = abs_array[index_low] - m*WN_array[index_low]
    baseline_y = np.array(m*WN_array[index_low:index_high] + b)
    baseline_area = simps(baseline_y, WN_array[index_low:index_high])
    areaarray.append(area - baseline_area)
    namearray.append(columnname)
    # /peak integration & baseline subtraction
    
    colorcount += 1
    if colorcount == colorlength:
        colorcount = 0

#print(df_IR_total.iloc[500:510])
#print(areaarray)
#print(namearray)
plt.bar(namearray, areaarray) #BAR GRAPH INTEGRATION 
plt.show()