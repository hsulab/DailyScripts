# -*- coding: utf-8 -*-
"""
Created on Thu Oct 19 21:55:56 2017
@author: Alexander
"""
#!/usr/bin/env python3
#import libs
import sys
import csv
import matplotlib.pyplot as plt
#plot_style
plt.style.use('ggplot')
#datas
input_file = sys.argv[1]
my_columns = ['suf','ab','ts','fs']
# data_form
# suf ab ts fs
# TOTEN ~~~~
# this is methane TOTEN is calculated by PBE+D3
methane_TOTEN=float(-24.07015819)
# get data
with open(input_file, 'r', newline='') as csv_in_file:
    filereader = csv.reader(csv_in_file, delimiter=',')
    header = next(filereader, None)
    my_columns_index = []
    for index_value in  range(len(header)):
        if header[index_value] in my_columns:
            my_columns_index.append(index_value)
    plot_data = []
    for row_list in filereader:
        for index_value in my_columns_index:
            plot_data.append(row_list[index_value])
    plot_data = [float(i) for i in plot_data]
    plot_data[0] = plot_data[0] + methane_TOTEN
    max_data = max(plot_data)
    plot_data = [i - max_data for i in plot_data]    
#start pic
fig = plt.figure()
ax1 = fig.add_subplot(1,1,1)
ax1.plot(plot_data, marker=r'*', color=u'red', linestyle='-',\
         label='State')
for i in range(len(plot_data)):
    plt.text(i+0.05, plot_data[i]+0.05, '%.4f' % plot_data[i], ha='center', va= 'bottom',fontsize=7)
    plt.text(i+0.05, plot_data[i]-0.1, '%s' % my_columns[i], ha='right', va= 'bottom',fontsize=10)
ax1.xaxis.set_ticks_position('bottom')
ax1.yaxis.set_ticks_position('left')
pic_title = str(input_file).split('.')[0]
ax1.set_title(pic_title)
plt.xlabel('Reaction Coordinate')
plt.ylabel('Energy / eV')
plt.legend(loc='best')
plt.savefig(pic_title, dpi=400, bbox_inches='tight')
#plt.show()
