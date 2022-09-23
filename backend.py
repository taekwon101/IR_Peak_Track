from distutils.command.build_scripts import first_line_re
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from scipy.integrate import simps
import glob, os

readpath = r'CSVs\220921_apk_PDMSCompare_uncured-7A30s808nm'
writepath = r'CSVs\Output'
plotScatter = True
plotBar = False
control_number = 5

def main(readpath, writepath, plotScatter, plotBar, control_number):
        
    # assign file paths/names
    colortext = open('textfiles/colorlist_lear.txt') #alternative: colorlist.txt
    output = '220513_cbPDMS_laserrots_oven.csv'
    willExport = False

    # change working directory to folder with CSVs of interest
    os.chdir(readpath) 

    # make list of names of csvs from readpath directory 
    filelist = sorted(glob.glob('*.csv'))  

    # normalization and baseline correction wavenumbers
    WN_normal_CH_low = 2957
    WN_normal_CH_high = 2967
    WN_baseline_low = 3400 
    WN_baseline_high = 3600

    # choose wavenumber range of interest (to integrate over for bar graph)
    WN_group = 0
    WN_low = [715, 940, 2290, 2900, 3060]
    WN_high = [830, 1230, 2390, 2970, 3080]
    groupname = ['Si-O-Si (?)', 'Si-O-Si (?)', 'Si-H', 'CH3', 'vinyl']
    note = ['?', '?', 'more cure = lower signal', 'more cure = same signal (troubleshooting)', 'more cure = lower signal']
        
    # select type of plot to show (true = scatter, false = bar)
    plotScatter = True
    plotBar = False
    plotBox = False
    plotkind = 'line' #line, bar, barh, hist, box, scatter, etc. for plotScatter plot

    # matplot formatting/scaling values
    xmin = 600
    xmax = 4000
    width = 0.8 # primary bar plot bar width
    lwidth = 2

    # assign colors for differentiating overlapping spectra in scatterplots (currently cycles across custom contrast gradient)
    colorlist = colortext.read().split()
    colorlength = len(colorlist)

    # prepare dataframes/arrays/variables for loop
    df_IR_total = pd.read_csv(filelist[0], skiprows = 2, names = ['wavenumber','temp']) # read first CSV in designated folder and write wavenumbers to primer dataframe and an array
    if df_IR_total['temp'].max() < 20: # check if CSVs are for absorbance or transmittance signals before primer dataframe drops signal column
        isAbs = True
    else:
        isAbs = False

    df_IR_total = df_IR_total.drop(['temp'], axis = 1)
    WN_array = np.array(df_IR_total['wavenumber'])
    df_IR_total.set_index('wavenumber', append = False, inplace = True)
    df_IR_total_abs = df_IR_total

    namearray = []
    areaarray = []
    changearray = []
    firstarray = []
    listcount = 0
    area_control = 0
    colorcount = colorlength # starts at max value for modulo loop count
    ax = plt.gca() # define constant axis for iterable additions

    # find indices that correspond to integration bounding and normalization wavenumbers so they can be referenced directly
    index_low = 0
    index_high = 0
    index_normal_CH_low = 0
    index_normal_CH_high = 0
    index_baseline_low = 0
    index_baseline_high = 0

    for WN in WN_array:
        if WN < WN_low[WN_group]:
            index_low += 1
        if WN < WN_high[WN_group]:
            index_high += 1
        if WN < WN_normal_CH_low:
            index_normal_CH_low += 1
        if WN < WN_normal_CH_high:
            index_normal_CH_high += 1
        if WN < WN_baseline_low:
            index_baseline_low += 1
        if WN < WN_baseline_high:
            index_baseline_high += 1

    # append signal columns from each CSV in directory to the intial dataframe stepwise in a loop
    for file in filelist:     
        columnname = file[0:-4] # extract column name (slicing off '.csv' from filename)
        df_temp = pd.read_csv(file, skiprows = 2, names = ['wavenumber', columnname]) # temp dataframe with csv

        # transmittance-to-absorbance conversion (on condition that it is not already in absorbance)
        if isAbs == False:
            df_temp[columnname] /= 100
            df_temp[columnname] += 1e-10 #bandaid for 'log of 0' errors
            df_temp[columnname] = np.log10(df_temp[columnname]) * -1
            
        # baseline subtraction
        baseline_temp = df_temp.iloc[index_baseline_low:index_baseline_high][columnname].mean(axis = 0) # baseline correction value to subtract
        df_temp[columnname] -= baseline_temp

        # normalization
        normal_temp = df_temp.iloc[index_normal_CH_low:index_normal_CH_high][columnname].mean(axis = 0) # signal value for current csv at normalization index
        df_temp[columnname] /= normal_temp 
        
        # add current signals column to total dataframe
        df_IR_total = pd.concat([df_IR_total, df_temp[columnname]], axis = 1)
        
        # superimpose corrected IR data with cycling color to existing axis for scatterplot
        if plotScatter == True:
            df_temp.plot(kind = plotkind, x = 'wavenumber', y = columnname, linewidth = lwidth, c = colorlist[colorcount % colorlength], label = columnname, ax = ax) # old: s = plotSize, c = colorlist[colorcount % colorlength]
        colorcount += 1 # out of conditional because also serves as loop counter below
        
        # peak integration, baseline subtraction, and array collection of values for bar graph (ty Nate and Sarah)
        if plotBar == True:
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
                pctchange = 100 * ((area - baseline_area) - area_control)/area_control
                changearray.append(pctchange)
            
            # turn the area_control value from a sum into an average by dividing by # of elements
            if colorcount == colorlength + control_number:
                area_control /= control_number

        # stores startpoints for repeat trials
        if plotBox == True:
            listcount += 1
            if columnname[0:2] == '01':
                firstarray.append(listcount)

    # create scatterplot
    if plotScatter == True:
        plt.xlim((xmin, xmax))
        x_label = ('wavenumbers /cm-1')
        plt.xlabel(x_label)
        y_label = 'absorbance (normalized)'
        plt.ylabel(y_label)
        plt.title("{} IR spectra from {}-{} cm-1".format(os.path.basename(os.path.normpath(readpath)), xmin, xmax)) # uses folder name as plot title

    #df_IR_total.plot(kind = plotkind, y = df_IR_total.columns, colormap = 'tab20') # old: y = columnname, s = plotSize, c = colorlist[colorcount % colorlength]

    # create bar graph
    if plotBar == True:
        print(changearray)
        plt.bar(namearray, areaarray, width = width, color = 'b', label = 'Integrated Signals')
        plt.bar(namearray, changearray, width = 0.5 * width, color = 'r', alpha = 0.8, label = 'percent change')
        y_label = 'integrated signal' #('integrated signal from {WN_low} cm-1 to {WN_high} cm-1')
        plt.ylabel(y_label)
        plt.xticks(rotation = 45, ha = 'right')
        plt.title("{} integrated over {}-{} cm-1 ({})".format(os.path.basename(os.path.normpath(readpath)), WN_low[WN_group], WN_high[WN_group], groupname[WN_group])) # uses folder name as plot title

    # if plotBox == True:
    #     boxdata = []
    #     for first in firstarray:
    #         boxdata[first] = first
    #         for 



    # global plot formatting and printing
    plt.legend()
    plt.show()

    # export (arrays to DF to CSV)
    if willExport == True:
        df_export = pd.DataFrame({"sample (cor./norm.)" : namearray, "area" : areaarray, "pct change" : changearray})
        df_export.to_csv(writepath, index = False)

main(readpath, writepath, plotScatter, plotBar, control_number);