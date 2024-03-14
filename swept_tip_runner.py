import shutil
import os
import warnings
import numpy as np
import scipy.io as spio

import sharpy.sharpy_main

from wing_generator import swept_tip_goland
import helper_functions.case_data_extract as case_data_extract
import helper_functions.biot_savart as biot_savart

### Disable warnings for more readable output
warnings.filterwarnings("ignore")    

### Remove old case files and outputs
try:
        shutil.rmtree('' + str(os.path.dirname(os.path.realpath(__file__))) + '/output/')
        shutil.rmtree('' + str(os.path.dirname(os.path.realpath(__file__))) + '/cases/')
except: ()

case_output = dict()

flow =  ['BeamLoader', 
        'AerogridLoader',
        # 'StaticCoupled',
        # 'Modal',
        # 'AerogridPlot',
        # 'BeamPlot',
        'DynamicCoupled',
        'AeroForcesCalculator',
        # 'LinearAssembler',
        # 'AsymptoticStability',
        ]

# Case parameters
ang_h = np.deg2rad(30)
pos_frac_h = 0.7
u_inf = 50
c_ref = 1.8288
ar = 7  

# Generate unique case name for each case
case_name = 'swept_tip_case'

# Generate wing and run SHARPy
wing = swept_tip_goland(case_name, flow,
                        disc_mode = 5,
                        c_ref = c_ref,
                        ang_h = ang_h,
                        pos_frac_h = pos_frac_h,
                        u_inf = u_inf,
                        n_surf = 1,
                        b_ref = c_ref*ar,
                        gust_intensity = 0.5,
                        gust_length = 0.2 * u_inf,
                        gust_offset = 0.1 * u_inf,
                        physical_time = 1.5,
                        M = 10,
                        N = 40,
                        Mstar_fact = 4,
                        wake_cfl1 = False,
                        write_screen = 'on')
        
case_data_curr = sharpy.sharpy_main.main(['', wing.route + wing.case_name + '.sharpy'])