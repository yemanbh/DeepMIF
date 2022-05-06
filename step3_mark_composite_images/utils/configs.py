# -*- coding: utf-8 -*-
"""
Created on 29/05/2020

@author: yhagos
"""
# +++++++++++++++++++++++++++++++++++++++ Analysis section variables ++++++++++++++++++++++++++++++++


hue = 'Class'  # 'Class' 'Component'

using_co_expression = True
scale = 1.0
on_composite = True
composite_legend = False
marker_size = 16
scatter_plot = True
scatter_legend = False
scatter_plot_save_format = '.jpg'
on_composite_save_format = '.jpg'

def setFigureObjectProperties():
	from matplotlib import rc, rcParams, pyplot as plt
	import matplotlib as mpl
	mpl.rcParams['axes.spines.right'] = False
	mpl.rcParams['axes.spines.top'] = False
	rcParams['pdf.fonttype'] = 42
	font = {'family' : 'Arial',
			'size':7,
	        'weight' : 'bold'
	        }
	axes= { 'labelsize':7,
	        'labelweight' : 'bold'
	        }
	rc('font', **font)  # pass in the font dict as kwargs
	plt.rc('axes', **axes)
