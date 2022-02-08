# -*- coding: utf-8 -*-
"""
Created on Wed Sep  8 15:53:34 2021

@author: lass_j
"""

from MJONIR.Geometry import Instrument
import numpy as np


folder = r'C:/Users/lass_j/Documents/McStasCAMEATestData/phonon_test_6_5mev_2'

calibrationFiles = ["C:/Users/lass_j/Documents/CAMEA2021/Normalization_1.calib",
                    "C:/Users/lass_j/Documents/CAMEA2021/Normalization_3.calib",
                    "C:/Users/lass_j/Documents/CAMEA2021/Normalization_5.calib",
                    "C:/Users/lass_j/Documents/CAMEA2021/Normalization_8.calib"]


Instrument.convertToHDF(r'C:\Users\lass_j\Documents\McStasCAMEATestData\test.hdf',
             'Noget_tm',
             sample='phonon',
             fname = folder,
             CalibrationFile=calibrationFiles,
             plane_normal=np.array([0,0,1]))
