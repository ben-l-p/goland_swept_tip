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
        'StaticCoupled',
        # 'Modal',
        'AerogridPlot',
        'BeamPlot',
        # 'DynamicCoupled',
        # 'AeroForcesCalculator',
        # 'LinearAssembler',
        # 'AsymptoticStability',
        ]

# Case parameters
ang_h_vect = list(np.deg2rad(np.linspace(-45, 45, 7)))
pos_frac_h_vect = list(np.linspace(0.5, 0.95, 10))


u_inf = 50
c_ref = 1.8288
ar = 7  

data_out = np.empty([len(ang_h_vect), len(pos_frac_h_vect)], dtype=dict)

for i_ang, _ in enumerate(ang_h_vect):
        for i_pos, _ in enumerate(pos_frac_h_vect):
                data_out[i_ang, i_pos] = dict()

for i_ang, ang_h in enumerate(ang_h_vect):
        for i_pos, pos_frac_h in enumerate(pos_frac_h_vect):
                # Generate unique case name for each case
                case_name = f'swept_tip_ang_{np.rad2deg(ang_h):.2f}_pos_{pos_frac_h:.2f}'
                case_name = case_name.replace('.', '_')
                case_name = case_name.replace('-', 'n') 
                print(case_name) 
                print(f'Case {i_pos + i_ang*len(pos_frac_h_vect) + 1} of {len(ang_h_vect)*len(pos_frac_h_vect)}\n')

                # Skip multiple straight wing cases
                if np.abs(ang_h) < 1e-5 and i_pos >= 1:
                        print("Skipping duplicate straight wing case")
                        data_out[i_ang, i_pos] = data_out[i_ang, 0]
                        continue

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
                                        physical_time = 1.3,
                                        M = 12,
                                        N = 60,
                                        sigma_1 = 1.0,
                                        Mstar_fact = 3.5,
                                        wake_cfl1 = False,
                                        write_screen = 'on')
        
                case_data = sharpy.sharpy_main.main(['', wing.route + wing.case_name + '.sharpy'])
                data_out[i_ang, i_pos] = case_data_extract.case_data_extract(wing, case_data)
                spio.savemat('case_data.mat', {'data': data_out})