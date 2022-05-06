from math import floor
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from scipy import stats
from scipy.integrate import simps
import glob, os

# assign file paths/names
readpath = r'C:\Users\taekw\Desktop\1_PythonScripts\IRPeakExtract\CSVs\220504_808nm_5A_5e-5cbPDMS'
writepath = r'C:\Users\taekw\Desktop\1_PythonScripts\IRPeakExtract\CSVs\Output\220504_808nm_5A_5e-5cbPDMS.csv'
output = '220504_808nm_5A_5e-5cbPDMS.csv'
willExport = True

os.chdir(readpath) #change working directory to folder with CSVs of interest
filelist = sorted(glob.glob('*.csv')) #make list of names of csvs from readpath directory  

# normalization and baseline correction wavenumbers
WN_normal_CH = 2960
WN_normal_CH_high = 3000
WN_baseline_low = 3400 
WN_baseline_high = 3600

# choose wavenumber range to integrate over for bar graph
WN_group = 1

if WN_group == 1: # (Si-O-Si): more cure --> HIGHER signal
    WN_low = 715
    WN_high = 830
    groupname = 'Si-O-Si'
elif WN_group == 2: # (Si-O-Si): more cure --> HIGHER signal
    WN_low = 940
    WN_high = 1230
    groupname = 'Si-O-Si'
elif WN_group == 3: # silane (Si-H): more cure --> LOWER signal
    WN_low = 2290
    WN_high = 2390
    groupname = 'Si-H'
elif WN_group == 4: # (-CH2- and -CH3): more cure --> SAME signal (for checking if internal standards are changing)
    WN_low = 2900
    WN_high = 2970
    groupname = 'CH2 + CH3'
elif WN_group == 5: # vinyl (CH=CH2) asym.: more cure --> LOWER signal (often so small, gets lost to baseline subtraction)
    WN_low = 3060
    WN_high = 3080
    groupname = 'vinyl'
    
# select type of plot to show (true = scatter, false = bar)
plotScatter = False

# matplot manual formatting/scaling values
manual_y = False
ymin = -50
ymax = 65
xmin = 400
xmax = 3500
width = 0.8 # primary bar plot bar width
plotSize = 1 # size of dots in scatterplot

# assign control number (loop will average over this number of initial areas for control %change-from value)
control_number = 10

# assign colors for differentiating overlapping spectra in scatterplots (currently cycles across custom contrast gradient)
colorlist = ['#cd6155', '#ec7063', '#af7ac5', '#a569bd', '#5499c7', '#5dade2', '#48c9b0', '#45b39d', '#52be80', '#58d68d', '#f4d03f', '#f5b041'] #'eb984e', '#dc7633', '#f0f3f4', '#cacfd2', '#aab7b8', '#99a3a4', '#5d6d7e', '#566573']#['#d55e00', '#cc79a7', '#0072b2', '#f0e442', '#009e73', '#24ff24','#000000","#004949","#009292","#ff6db6","#ffb6db","#490092","#006ddb","#b66dff","#6db6ff","#b6dbff","#920000","#924900","#db6d00","#24ff24","#ffff6d]
colorlength = len(colorlist)

# prepare dataframes/arrays/variables for loop
df_IR_total = pd.read_csv(filelist[0], names = ['wavenumber','temp']) # read first CSV in designated folder and write wavenumbers to primer dataframe and an array
isAbs = False # check if CSVs are for absorbance or transmittance signals before primer dataframe drops signal column
if df_IR_total['temp'].max() < 20:
    isAbs = True
df_IR_total = df_IR_total.drop(['temp'], axis = 1)
WN_array = np.array(df_IR_total['wavenumber'])
df_IR_total.set_index('wavenumber', append = False, inplace = True)
df_IR_total_abs = df_IR_total
namearray = []
areaarray = []
changearray = []
area_control = 0
colorcount = colorlength # start at max value for modulo loop count
ax = plt.gca() # define constant axis for iterable additions

# find indices that correspond to integration bounding and normalization wavenumbers so they can be referenced directly
index_low = 0
index_high = 0
index_normal_CH = 0
index_normal_CH_high = 0
index_baseline_low = 0
index_baseline_high = 0
for item in WN_array:
    if item < WN_low:
        index_low += 1
    if item < WN_high:
        index_high += 1
    if item < WN_normal_CH:
        index_normal_CH += 1
    if item < WN_normal_CH_high:
        index_normal_CH_high += 1
    if item < WN_baseline_low:
        index_baseline_low += 1
    if item < WN_baseline_high:
        index_baseline_high += 1

# append signal columns from each CSV in directory to the intial dataframe stepwise in a loop
for file in filelist:     
    columnname = file[0:-4] # extract column name (slicing off '.csv' from filename)
    df_temp = pd.read_csv(file, names = ['wavenumber', columnname]) # temp dataframe with csv
    normal_temp = df_temp.iloc[index_normal_CH][columnname] # signal value for current csv at normalization index
    baseline_temp = df_temp.iloc[index_baseline_low:index_baseline_high][columnname].mean(axis = 0) # baseline correction value to subtract
    #baseline_temp = df_temp.iloc[index_baseline_low][columnname] # baseline correction value to subtract (legacy - keeping in case average method above breaks)

    # transmittance-to-absorbance conversion (on condition that it is not already in absorbance)
    if isAbs == False:
        df_temp[columnname] /= 100
        df_temp[columnname] += 1e-10 #bandaid for 'log of 0' errors
        df_temp[columnname] = np.log10(df_temp[columnname]) * -1
        
    # baseline subtraction (breaks sometimes and makes values negative?)
    df_temp[columnname] -= baseline_temp
    
    # normalization
    df_temp[columnname] /= normal_temp 
    
    # add current signals column to total dataframe
    df_IR_total = pd.concat([df_IR_total, df_temp[columnname]], axis = 1)
    
    # superimpose corrected IR data with cycling color to existing axis for scatterplot
    if plotScatter == True:
        df_temp.plot(kind = 'scatter', x = 'wavenumber', y = columnname, s = plotSize, c = colorlist[colorcount % colorlength], label = columnname, ax = ax)
    colorcount += 1 # out of conditional because also serves as loop counter below
    
    # peak integration, baseline subtraction, and array collection of values for bar graph (ty Nate and Sarah)
    if plotScatter == False:
        abs_array = np.array(df_temp[columnname])
        area = simps(abs_array[index_low:index_high], WN_array[index_low:index_high])
        m = (abs_array[index_low] - abs_array[index_high])/(WN_array[index_low] - WN_array[index_high])
        b = abs_array[index_low] - m*WN_array[index_low]
        baseline_y = np.array(m*WN_array[index_low:index_high] + b)
        baseline_area = simps(baseline_y, WN_array[index_low:index_high])
        namearray.append(columnname)
        areaarray.append(area - baseline_area)

        
        # find % change values to populate second bar graph (starting with checking if first time through loop for control value)
        if colorcount <= colorlength + control_number:
            area_control += (area - baseline_area)
            changearray.append(0)
        else:    
            changearray.append(100 * ((area - baseline_area) - area_control)/area_control)

        if colorcount == colorlength + control_number:
            area_control /= control_number
# create scatterplot
if plotScatter == True:
    if manual_y == True:
        plt.ylim((ymin, ymax))
    plt.xlim((xmin, xmax))
    x_label = ('wavenumbers /cm-1')
    plt.xlabel(x_label)
    y_label = 'absorbance (normalized)'
    plt.ylabel(y_label)
    plt.title("{} IR spectra from {}-{} cm-1".format(os.path.basename(os.path.normpath(readpath)), xmin, xmax)) # uses folder name as plot title

# create bar graph
if plotScatter == False:
    if manual_y == True:
        plt.ylim((ymin, ymax))
    plt.bar(namearray, areaarray, width = width, color = 'b', label = 'Integrated Signals')
    plt.bar(namearray, changearray, width = 0.5 * width, color = 'r', alpha = 0.8, label = 'percent change')
    y_label = 'integrated signal' #('integrated signal from {WN_low} cm-1 to {WN_high} cm-1')
    plt.ylabel(y_label)
    plt.xticks(rotation = 45)
    plt.title("{} integrated over {}-{} cm-1 ({})".format(os.path.basename(os.path.normpath(readpath)), WN_low, WN_high, groupname)) # uses folder name as plot title

# global plot formatting and printing
plt.legend()
plt.show()

# export (arrays to DF to CSV)
if willExport == True:
    df_export = pd.DataFrame({"sample (cor./norm.)" : namearray, "area" : areaarray, "pct change" : changearray})
    df_export.to_csv(writepath, index = False)
