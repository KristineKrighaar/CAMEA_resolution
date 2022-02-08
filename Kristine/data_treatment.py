import numpy as np
import matplotlib.pyplot as plt
from iminuit import Minuit
from scipy import stats
import os
import subprocess

#import seaborn as sns
import sys
import scipy
#sys.path.append('External_Functions')
#from ExternalFunctions import Chi2Regression
#from ExternalFunctions import nice_string_output, add_text_to_ax 

from sympy import *
import random
import math

import matplotlib as mpl
mpl.rcParams['figure.figsize']   = (10,6)
mpl.rcParams['font.size']        = 15 # standard er 45
mpl.rcParams['lines.color']      = 'r'
mpl.rcParams['lines.markersize'] = 15
plt.rcParams['figure.constrained_layout.use'] = True



calib1 = np.arange(0,8,1,)

pathname_old_mid = '~/Dropbox/Masters_2/CAMEA_Resolution_project/Resolution_intrument/Small_stat_mid_detector/elastic_old_mono_'

for i in calib1:
    path = str(pathname_old_mid+str(i)+'calib1/')
    to_folder_command = 'cd '+ path
    print(to_folder_command)
    subprocess.call(to_folder_command, shell=True)
    
    long_file = 'res_mon_mid_'+str(i)+'_list.ki_x.ki_y.ki_z.kf_x.kf_y.kf_z.x.y.z.p_i.p_f'
    mc_ras_run = 'mcresplot '+long_file
    print(mc_ras_run)
    subprocess.call(mc_ras_run, shell=True)

