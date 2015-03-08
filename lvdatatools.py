#better ideas....

#all of these are meant to work with the CNT current measurements. assume you are always measuring current versus 
#some other parameters. 

#there are really only two types of data... 2D data (current versus something) and 3D data (current
#versus something as a function of something else)

#split into a few things... grab data and relevant information, plot data (return axes instance like Dan suggested?)

from __future__ import division
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import cm
from scipy import signal
import os, struct, time

def get_swp_info(filename):
    """ get data shape from vi file """
    point_bytes = 8.0 #bytes per point
    points = os.path.getsize(filename)/point_bytes
    f = open(filename, 'rb')
    f.seek(0)
    columns = struct.unpack('>d', f.read(8))[0]
    rows = points/columns
    f.close()
    return (int(rows), int(columns))
    
def get_column_data(filename, column_num, shape):
    """ get data from one column without reading the whole file """
    point_bytes = 8.0
    row_bytes = shape[1]*point_bytes
    column = np.zeros(shape[0])
    f = open(filename, 'rb')
    for i in range(shape[0]):
        f.seek(int(i*row_bytes + column_num*point_bytes), 0) #could be faster to seek from current position not start
        column[i] = struct.unpack('>d', f.read(8))[0]
    f.close()
    return column
    
def vi_gateswp_build(filelist):
    """ take a list of files created by the vi_gateswp VI and compile
        them into a matrix of current values that can be plotted quickly """
    
    #this is so fucking slow it's unbelievable
    
    start_time = time.time()
    
    shape = get_swp_info(filelist[0]) #get shape from first file
    bias = get_column_data(filelist[0], 1, shape) #get bias from first file
    
    gate = np.zeros(len(filelist)) 
    current = np.zeros((len(bias),len(gate)))
    
    print 'after setup = {0}'.format(time.time()-start_time)
    
    for i, filename in enumerate(filelist):   
        current[:,i] = get_column_data(filename, 2, shape)
        m = re.search('_(n?)(\d+).bin', filename)
        gate[i] =  (1-2*len(m.group(1)))*float(m.group(2))/1000 #make it negative if i find an 'n'
        print 'after {0} = {1}'.format(i, time.time()-start_time)
    return pd.DataFrame(current, index = bias, columns=gate).sort(axis=1)




############# new and more useful stuff below here ##################


def get_txt(filename):
    """ creates a string array out of .txt data autosaved
        by LabView"""
    if not filename.endswith('.txt'):
        filename = filename[:-4]+'.txt'
    if os.path.isfile(filename):
        f = open(filename)
    elif os.path.isfile(filename[:-9]+'.txt'):
        f = open(filename[:-9]+'.txt') #if it is a section of a 3d sweep
    else:
        f = open(filename[:-10]+'.txt') #if it is a negative section of a 3d sweep
    f.seek(0)
    header = []
    while True:
        line = f.readline()
        if len(line) == 0:
            break
        else:
            line = line.rstrip('\n\r')
            line = line.split(': ')
            if len(line)==1:
                line.append('')
            header.append(line)
    return header
    
def output_format(header):
    """ takes the header file and determines what type of measurement
        was performed. Returns a string array corresponding to column headers."""
    for i,line in enumerate(header):
        if line[0] == 'Binary Output Format':
        	names = line[1].split(', ')
        	if names[-1] == 'I(t)':
        		return names[:-1]
        	else:
        		return names

def sweep_variables(header):
    """ find the two variables swept in the measurement. this fails hard
        if something is labelled wrong in the txt file. """
    fast = output_format(header)[1]
    slow = ''
    for line in header:
        if ' Sweep' in line[0]:
            for word in line[0].split():
                if len(word)>1 and word!='Sweep' and word.lower()!=fast.lower():
                    slow=word
    return fast, slow
    
def split_sweeps(df, col_name):
	""" get indicies to split multiple sweeps and plot separately. """
	# this could really be much better
	# should have different options for different sweep types
	#	right now it just returns zeros of the first derivative
    ind = np.where(np.diff(df[col_name])==0)[0]
    return np.concatenate(([0],ind,[len(df[col_name])]), axis=0)
    
def is_2d(filename):
    """ uses header info to try to determine if the data is 2d or 3d """   
    if os.path.isfile(filename[:-4]+'.txt'):
        header = get_txt(filename)
        fast, slow = sweep_variables(header)
        if slow:
            return 3 #has two sweep axes
        else: 
            return 2 #only has one sweep axis
            
    else:
        return 1 #if it doesn't have it's own .txt file it must be a 2d cut of a 3d measurement 
                 #or totally screwed up


def get_data_2d(filename):
    """ returns a pandas array with the 2d data from the binary file
        tagged with relevant information from the header file """
    header = get_txt(filename)
    data = np.fromfile(filename,'>d')
    data = data.reshape((-1,data[0]))
    data = data.byteswap().newbyteorder()
    col_names = output_format(header)
    if data[0,0] <= 6.0:
        return pd.DataFrame(data, columns=col_names)
    else:
        col_names.extend(['i'+str(n) for n in range(int(data[0,0]-6))])
        return pd.DataFrame(data, columns=col_names)

def get_data_3d(filename):
    """ returns a pandas array with the 3d data from the binary file
        tagged with the relevant information from the header file """
    header  = get_txt(filename)
    fast_scan, slow_scan = sweep_variables(header)
    data = np.fromfile(filename,'>d')
    data = data.reshape((-1, data[0]))
    data = data.byteswap().newbyteorder()
    data = data.transpose() #20254, 20472
    df = pd.DataFrame(data[1:, 1:], index=data[1:, 0], columns=data[0, 1:])
    df.columns.name = slow_scan
    df.index.name = fast_scan
    return df

def df_extent(df):
    """ get the extent for imshow images using dataframs """
    return [df.columns.values.min(), df.columns.values.max(),
            df.index.values.min(), df.index.values.max()]
    
def replace_ticks(ax, xlabels = None, ylabels = None):
    """ replace existing labels with labels generated elsewhere. 
        should work for pretty much any axes instance. """
    if xlabels is not None:
        xmin, xmax = ax.get_xlim()
        ax.set_xticks(np.linspace(xmin, xmax, len(xlabels)))
        ax.set_xticklabels(xlabels)
    if ylabels is not None:
        ymin, ymax = ax.get_ylim()
        ax.set_yticks(np.linspace(ymin, ymax, len(ylabels)))
        ax.set_yticklabels(ylabels)
    
def plot_simple_2d(df):
    """ a fast way of recreating the labview plots so I can tell whats what """
    fig = plt.figure(figsize=(8,6))
    ax = fig.add_subplot(1,1,1)
    xlabel = df.columns.values[1]
    ylabel = df.columns.values[2]
    ax.plot(df.iloc[:,1], df.iloc[:,2], linewidth=2)
    ax.set_xlabel(r'${0}$'.format(xlabel), fontsize=14)
    ax.set_ylabel(r'${0}$'.format(ylabel), fontsize=14)
    return fig, ax #this isn't quite right/works well enough for now

def plot_simple_3d(df):
    """ this is a little sloppy, but it gets the point across quickly """
    df = df.groupby(df.index).mean()
    fig = plt.figure(figsize=(12,6))
    ax = fig.add_subplot(1,1,1)
    limits = [df.columns.values.min(), df.columns.values.max(),df.index[0],df.index[-1]]
    im = ax.imshow(df*1e9, extent = limits, 
               cmap = plt.cm.seismic,
               origin = 'lower', aspect = 'auto', interpolation = 'None')
    ax.set_xlabel(r'${0}$'.format(df.columns.name), fontsize=16)
    ax.set_ylabel(r'${0}$'.format(df.index.name), fontsize=16)
    cb = plt.colorbar(im)
    cb.set_label(r'$I (nA)$', fontsize=16)
    return fig, ax, cb

def quick_plot_any(filename):
    """ for use in a for-loop when trolling through a whole folder of data.
        trying to catch exceptions as they arise. """
    dimension = is_2d(filename)
    if dimension==2:
        try:
            df = get_data_2d(filename)
            fig, ax = plot_simple_2d(df)
        except(AttributeError, ValueError, IndexError): 
            return 'Fail -- 2d. A real problem.'
    elif dimension==3:
        try:
            df = get_data_3d(filename)
            fig, ax, cb = plot_simple_3d(df)
        except(AttributeError, ValueError, IndexError): 
            return 'Fail -- 3d. A real problem.'
    else:
        return 'Ignored. No .txt file to match.'
    ax.set_title(filename[:-4])
    fig.savefig(filename[:-4]+'_quick-plot.png')
    plt.close()
    return 'Pass -- {0}d. Good job, buddy.'.format(dimension)
