from math import floor
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from scipy import stats
from scipy.integrate import simps
import glob, os

# assign file paths/names
readpath = r'C:\Users\taekw\Desktop\1_PythonScripts\IRPeakExtract\CSVs\220413_808nm_5A_5e-5cbPDMS'
writepath = r'C:\Users\taekw\Desktop\1_PythonScripts\IRPeakExtract\CSVs\Output'
output = '220317_5e7 cbPDMS_laser 10A oven 100C_abstest.csv'
willExport = False

os.chdir(readpath) #change working directory to folder with CSVs of interest
filelist = sorted(glob.glob('*.csv')) #make list of names of csvs from readpath directory  

# assign normalization and baseline correction wavenumbers
WN_normal_CH = 2960
WN_baseline_low = 3200
WN_baseline_high = 3400

# select type of plot to show (true = scatter, false = bar)
plotScatter = False
#plotScatter = True
plotSize = 2

# wavenumber range to integrate over (shown in bar graph, see info.txt for functional groups)
WN_low = 940
WN_high = 1230

# matplot manual axes scaling values (scatter plot)
ymin = -0.5
ymax = 3.5
xmin = 700
xmax = 3000

# assign colors for differentiating overlapping spectra in scatterplots (currently alternates between 6 color-blind-friendly colors)
colorlist = ['#d55e00', '#cc79a7', '#0072b2', '#f0e442', '#009e73', '#24ff24']#, '#000000","#004949","#009292","#ff6db6","#ffb6db","#490092","#006ddb","#b66dff","#6db6ff","#b6dbff","#920000","#924900","#db6d00","#24ff24","#ffff6d]
colorlength = len(colorlist)

# prepare dataframes/arrays for loop
df_IR_total = pd.read_csv(filelist[0], names = ['wavenumber','temp']) # read first CSV in designated folder and write wavenumbers to primer dataframe and an array
isAbs = False # check if CSVs are for absorbance or transmittance signals before primer dataframe drops signal column
if df_IR_total['temp'].max() < 20:
    isAbs = True
df_IR_total = df_IR_total.drop(['temp'], axis = 1)
WN_array = np.array(df_IR_total['wavenumber'])
df_IR_total.set_index('wavenumber', append = False, inplace = True)
df_IR_total_abs = df_IR_total
areaarray = []
namearray = []
realcolorarray = []
colorcount = 0
ax = plt.gca() # plot axis for iterable additions

# find indices that correspond to integration bounding and normalization wavenumbers (stupidly inefficient brute force)
index_low = 0
for item in WN_array:
    if item < WN_low:
        index_low += 1
   
index_high = 0
for item in WN_array:
    if item < WN_high:
        index_high += 1
        
index_normal_CH = 0
for item in WN_array:
    if item < WN_normal_CH:
        index_normal_CH += 1
        
index_baseline_low = 0
for item in WN_array:
    if item < WN_baseline_low:
        index_baseline_low += 1
        
index_baseline_high = 0
for item in WN_array:
    if item < WN_baseline_high:
        index_baseline_high += 1

# append signal columns from each CSV in directory to the intial dataframe stepwise in a loop
for file in filelist:     
    columnname = file[0:-4] # extract column name (slicing off '.csv' from filename)
    df_temp = pd.read_csv(file, names = ['wavenumber', columnname]) # temp dataframe with csv
    normal_temp = df_temp.iloc[index_normal_CH][columnname] # signal value for current csv at normalization index
    baseline_temp = df_temp.iloc[index_baseline_low][columnname] #baseline correction value to subtract RETURN TO TO ACTAULLY AVERAGE VALUES FOR BASELINE
    
    # transmittance-to-absorbance conversion (on condition that it is not already in absorbance)
    if isAbs == False:
        df_temp[columnname] /= 100
        df_temp[columnname] += 1e-10 #bandaid for 'log of 0' errors
        df_temp[columnname] = np.log10(df_temp[columnname]) * -1
    
    # baseline subtraction
    df_temp[columnname] -= baseline_temp
    
    # normalization
    df_temp[columnname] /= normal_temp 
    
    # add current signals column to total dataframe
    df_IR_total = pd.concat([df_IR_total, df_temp[columnname]], axis = 1)
    
    if plotScatter == True:
        df_temp.plot(kind = 'scatter', x = 'wavenumber', y = columnname, s = plotSize, c = colorlist[colorcount], label = columnname, ax = ax)
        realcolorarray.append(colorlist[colorcount])
        colorcount += 1
        if colorcount == colorlength:
            colorcount = 0 # cycles back to beginning of color array after final entry
    
    namearray.append(columnname)
    
    # peak integration & baseline subtraction for bar graph (ty Nate and Sarah)
    if plotScatter == False:
        abs_array = np.array(df_temp[columnname])
        area = simps(abs_array[index_low:index_high], WN_array[index_low:index_high])
        m = (abs_array[index_low] - abs_array[index_high])/(WN_array[index_low] - WN_array[index_high])
        b = abs_array[index_low] - m*WN_array[index_low]
        baseline_y = np.array(m*WN_array[index_low:index_high] + b)
        baseline_area = simps(baseline_y, WN_array[index_low:index_high])
        areaarray.append(area - baseline_area)

# create scatterplot
if plotScatter == True:
    plt.ylim((ymin, ymax))
    plt.xlim((xmin, xmax))
    x_label = ('wavenumbers /cm-1')
    plt.xlabel(x_label)
    y_label = 'absorbance (normalized)'
    plt.ylabel(y_label)
    plt.legend()

# create bar graph
if plotScatter == False:
    #areaarray -= areaarray[0] # bandaid for too broad a y scale when all values are similar
    #plt.ylim((230, 340))
    plt.bar(namearray, areaarray)
    y_label = 'integrated signal' #('integrated signal from {WN_low} cm-1 to {WN_high} cm-1')
    plt.ylabel(y_label)
    plt.xticks(rotation = 45)

# global plot formatting and printing
plt.title("{} integrated over {}-{} cm-1".format(os.path.basename(os.path.normpath(readpath)), WN_low, WN_high)) # uses folder name as plot title
plt.show()

# export (arrays to DF to CSV)
if willExport == True:
    df_export = pd.DataFrame({"sample (cor./norm.)" : namearray, "area" : areaarray})
    df_export.to_csv("output.csv", index = False)
