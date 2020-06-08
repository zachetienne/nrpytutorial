# Piecewise Polytrope Dictionary
# Author: Leo Werneck
# Source: Read et al. (2009) PRD 79,124032 (https://arxiv.org/pdf/0812.2163.pdf)
from collections import namedtuple

# Initialize the named tuple EOS_Read_et_al_parameters
EOS_Read_et_al_tuple = namedtuple('Read_et_al_EOS_Parameters', 'log_of_p4 Gamma4 Gamma5 Gamma6')

# Initialize the dictionary ShockAvoidLapseDict
EOS_Read_et_al_dict = {}

# PAL6 EOS parameters
EOS_Read_et_al_dict['PAL6']   = ( EOS_Read_et_al_tuple( 34.380, 2.227, 2.189, 2.159 ) )

# SLy EOS parameters
EOS_Read_et_al_dict['SLy']    = ( EOS_Read_et_al_tuple( 34.384, 3.005, 2.988, 2.851 ) )

# APR1 EOS parameters
EOS_Read_et_al_dict['APR1']   = ( EOS_Read_et_al_tuple( 33.943, 2.442, 3.256, 2.908 ) )

# APR2 EOS parameters
EOS_Read_et_al_dict['APR2']   = ( EOS_Read_et_al_tuple( 34.126, 2.643, 3.014, 2.945 ) )

# APR3 EOS parameters
EOS_Read_et_al_dict['APR3']   = ( EOS_Read_et_al_tuple( 34.392, 3.166, 3.573, 3.281 ) )

# APR4 EOS parameters
EOS_Read_et_al_dict['APR4']   = ( EOS_Read_et_al_tuple( 34.269, 2.830, 3.445, 3.348 ) )

# FPS EOS parameters
EOS_Read_et_al_dict['FPS']    = ( EOS_Read_et_al_tuple( 34.283, 2.985, 2.863, 2.600 ) )

# WFF1 EOS parameters
EOS_Read_et_al_dict['WFF1']   = ( EOS_Read_et_al_tuple( 34.031, 2.519, 3.791, 3.660 ) )

# WFF2 EOS parameters
EOS_Read_et_al_dict['WFF2']   = ( EOS_Read_et_al_tuple( 34.233, 2.888, 3.475, 3.517 ) )

# WFF3 EOS parameters
EOS_Read_et_al_dict['WFF3']   = ( EOS_Read_et_al_tuple( 34.283, 3.329, 2.952, 2.589 ) )

# BBB2 EOS parameters
EOS_Read_et_al_dict['BBB2']   = ( EOS_Read_et_al_tuple( 34.331, 3.418, 2.835, 2.832 ) )

# BPAL12 EOS parameters
EOS_Read_et_al_dict['BPAL12'] = ( EOS_Read_et_al_tuple( 34.358, 2.209, 2.201, 2.176 ) )

# ENG EOS parameters
EOS_Read_et_al_dict['ENG']    = ( EOS_Read_et_al_tuple( 34.437, 3.514, 3.130, 3.168 ) )

# MPA1 EOS parameters
EOS_Read_et_al_dict['MPA1']   = ( EOS_Read_et_al_tuple( 34.495, 3.446, 3.572, 2.887 ) )

# MS1 EOS parameters
EOS_Read_et_al_dict['MS1']    = ( EOS_Read_et_al_tuple( 34.858, 3.224, 3.033, 1.325 ) )

# MS2 EOS parameters
EOS_Read_et_al_dict['MS2']    = ( EOS_Read_et_al_tuple( 34.605, 2.447, 2.184, 1.855 ) )

# MS1b EOS parameters
EOS_Read_et_al_dict['MS1b']   = ( EOS_Read_et_al_tuple( 34.855, 3.456, 3.011, 1.425 ) )

# PS EOS parameters
EOS_Read_et_al_dict['PS']     = ( EOS_Read_et_al_tuple( 34.671, 2.216, 1.640, 2.365 ) )

# GS1 EOS parameters
EOS_Read_et_al_dict['GS1']    = ( EOS_Read_et_al_tuple( 34.504, 2.350, 1.267, 2.421 ) )

# GS2 EOS parameters
EOS_Read_et_al_dict['GS2']    = ( EOS_Read_et_al_tuple( 34.642, 2.519, 1.571, 2.314 ) )

# BGN1H1 EOS parameters
EOS_Read_et_al_dict['BGN1H1'] = ( EOS_Read_et_al_tuple( 34.623, 3.258, 1.472, 2.464 ) )

# BNH3 EOS parameters
EOS_Read_et_al_dict['GNH3']   = ( EOS_Read_et_al_tuple( 34.648, 2.664, 2.194, 2.304 ) )

# H1 EOS parameters
EOS_Read_et_al_dict['H1']     = ( EOS_Read_et_al_tuple( 34.564, 2.595, 1.845, 1.897 ) )

# H2 EOS parameters
EOS_Read_et_al_dict['H2']     = ( EOS_Read_et_al_tuple( 34.617, 2.775, 1.855, 1.858 ) )

# H3 EOS parameters
EOS_Read_et_al_dict['H3']     = ( EOS_Read_et_al_tuple( 34.646, 2.787, 1.951, 1.901 ) )

# H4 EOS parameters
EOS_Read_et_al_dict['H4']     = ( EOS_Read_et_al_tuple( 34.669, 2.909, 2.246, 2.144 ) )

# H5 EOS parameters
EOS_Read_et_al_dict['H5']     = ( EOS_Read_et_al_tuple( 34.609, 2.793, 1.974, 1.915 ) )

# H6 EOS parameters
EOS_Read_et_al_dict['H6']     = ( EOS_Read_et_al_tuple( 34.593, 2.637, 2.121, 2.064 ) )

# H7 EOS parameters
EOS_Read_et_al_dict['H7']     = ( EOS_Read_et_al_tuple( 34.559, 2.621, 2.048, 2.006 ) )

# PCL2 EOS parameters
EOS_Read_et_al_dict['PCL2']   = ( EOS_Read_et_al_tuple( 34.507, 2.554, 1.880, 1.977 ) )

# ALF1 EOS parameters
EOS_Read_et_al_dict['ALF1']   = ( EOS_Read_et_al_tuple( 34.055, 2.013, 3.389, 2.033 ) )

# ALF2 EOS parameters
EOS_Read_et_al_dict['ALF2']   = ( EOS_Read_et_al_tuple( 34.616, 4.070, 2.411, 1.890 ) )

# ALF3 EOS parameters
EOS_Read_et_al_dict['ALF3']   = ( EOS_Read_et_al_tuple( 34.283, 2.883, 2.653, 1.952 ) )

# ALF4 EOS parameters
EOS_Read_et_al_dict['ALF4']   = ( EOS_Read_et_al_tuple( 34.314, 3.009, 3.438, 1.803 ) )
