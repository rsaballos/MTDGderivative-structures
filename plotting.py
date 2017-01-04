#!/usr/bin/python

"""
Created on Mon Sep 26 09:33:45 2016

@author: nenian charles
"""

import re
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

def plot():
    plt.ioff()
    pp = PdfPages('ordering_stats.pdf')

    f = open('ORDERING_OUTPUT.txt', 'r')

    SG_tally = []
    crystal_class_tally = []

    lines = f.readlines()
    for indx, line in enumerate(lines):
        if re.match('space_group #', line):
            SG_tally.append([int(line.split()[2]),int(line.split()[3])])
        if re.match("Anion ordering stats by crystal system.", line):
            for x in range(1,7):
                crystal_class_tally.append([lines[indx+x].split()[0],int(lines[indx+x].split()[1])])

    # Format space group data for plotting            
    SG_tally = np.array(sorted(SG_tally, key=lambda x: x[0]))
    xpos_SG = np.arange(len(SG_tally))
    
    # Format crystal class data for plotting
    crystal_class_tally = np.array(crystal_class_tally)
    count_crys_class = np.array([float(x) for x in crystal_class_tally[:,1]])

    #Compute the percent of each 
    percent = 100.*count_crys_class/count_crys_class.sum()
    
    
    #Begin plotting Ordering data
    fig = plt.figure()
    fig.set_size_inches(8.5, 11) #Output and 8.5x11 pdf
    
    plt.subplot(2,1,1)
    
    plt.bar(xpos_SG, SG_tally[:,1], align='center', alpha=0.5)
    plt.xticks(xpos_SG, SG_tally[:,0], rotation=-90)
    plt.title('Ordering stats')
    plt.xlabel('Space group number')
    plt.ylabel('Count')
    
    plt.subplot(2, 1, 2)
    colors = ['yellowgreen', 'gold', 'lightskyblue', 'lightcoral', 'peru', 'grey', 'darkblue']
    patches, texts = plt.pie(count_crys_class, colors=colors, shadow=True, startangle=90)
    labels = ['{0} - {1:1.2f} %'.format(i,j) for i,j in zip(crystal_class_tally[:,0], percent)]

    patches, labels, dummy =  zip(*sorted(zip(patches, labels, count_crys_class),
                                          key=lambda x: x[2],
                                          reverse=True))
    plt.legend(patches, labels,loc='best',fontsize=8)
    plt.axis('equal')
    plt.tight_layout()
    plt.savefig(pp, format='pdf', pad_inches=1)
    plt.close(fig)
    pp.close()
    
plot()