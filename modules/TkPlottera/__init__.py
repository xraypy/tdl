#/usr/bin/python
    
import os
import sys
import time
import types
import Tkinter
import matplotlib
matplotlib.rcParams['numerix'] ='numpy'

import matplotlib.numerix as numpy

matplotlib.use('TkAgg')


import PlotFrame, PlotPanel
PlotFrame = PlotFrame.PlotFrame
PlotPanel = PlotPanel.PlotPanel

