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
sweep = np.deg2rad(30)
u_inf = 50       
c_ref = 1.8288
ar = 10   
pnts = np.array([[3.0, 0.0, 0.0], [7.5, 7.5, 0.0],\
                [11.0, 14.0, 0.0], [15.0, 14.0, 0.0], \
                [11.0, 15.0, 0.0], [15.0, 15.0, 0.0],\
                [11.0, 16.0, 0.0], [15.0, 16.0, 0.0]]).T

# Dictionary of final converged paramters and case data   
case_data_final = dict()
wing_final = dict()
induced_v = dict()

### Discretisation study
### Loop this code, varying parameters each run
for i_disc in range(2):
        # Generate unique case name for each case
        case_name = ['aligned_conv', 'misaligned_conv'][i_disc]

        # Generate wing and run SHARPy
        wing = swept_tip_goland(case_name, flow,
                                c_ref = [c_ref/np.cos(sweep), c_ref][i_disc],
                                sweep_beam = sweep,
                                sweep_panel = [0, -sweep][i_disc],
                                u_inf = u_inf,
                                n_surf = 1,
                                b_ref = c_ref*ar,
                                gust_intensity = 0.5,
                                gust_length = 0.2 * u_inf,
                                gust_offset = 0.1 * u_inf,
                                physical_time = 1.5,
                                M = [18, 12][i_disc],
                                N = 66,
                                Mstar_fact = 4.2,
                                sym = False,
                                wake_cfl1 = False,
                                write_screen = 'on')
                
        case_data_curr = sharpy.sharpy_main.main(['', wing.route + wing.case_name + '.sharpy'])
                
        # biot_savart.plot_points(case_data_curr, pnts, 0)

        case_data_final.update({['aligned', 'misaligned'][i_disc] : case_data_curr})
        wing_final.update({['aligned', 'misaligned'][i_disc] : wing})
        induced_v.update({['aligned', 'misaligned'][i_disc] : biot_savart.biot_savart_points(case_data_curr, pnts)})
### End of loop
                                
# Extract key data from system and output to matlab
case_data_out = dict()
case_data_out.update({'aligned': case_data_extract.case_data_extract(\
        wing_final['aligned'], case_data_final['aligned'])})
case_data_out.update({'misaligned': case_data_extract.case_data_extract(\
        wing_final['misaligned'], case_data_final['misaligned'])})
case_data_out['aligned'].update({'induced_v': induced_v['aligned']})
case_data_out['misaligned'].update({'induced_v': induced_v['misaligned']})
case_data_out.update({'induced_v_pnts': pnts})
spio.savemat('convergence_out.mat', case_data_out)