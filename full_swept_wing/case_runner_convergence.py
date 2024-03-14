import shutil
import os
import warnings
import numpy as np
import scipy.io as spio

import sharpy.sharpy_main

from wing_generator import swept_tip_goland
import helper_functions.convergence_utils as convergence
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
disc_final = {'aligned': dict(), 'misaligned': dict()}
case_data_final = dict()
wing_final = dict()

# Discretisation study parameters, in M, N, M* order
disc_init = [3, 8, 20]
disc_mult = [1.4, 1.5, 1.5]
disc_order = ['Msf', 'M', 'N']
max_iter = 20                           # Maximum number of convergence iterations to run

### Discretisation study
### Loop this code, varying parameters each run
for i_disc in range(2):

        case_data_curr = None                   # Case data for current case
        case_data_conv = None                   # Case data for converged case                            
        n_conv = 0                              # Number of sequential runs which are converged
        is_final = False                        # Flag for running the final case

        # Generate the initial discretisation parameters
        disc_dict = convergence.disc_generate(disc_init, disc_mult, disc_order)

        for i_case in range(max_iter):
                # Select parameters to use for the case
                params = convergence.param_to_increase(disc_dict, disc_order, i_case, is_final)

                # Generate unique case name for each case
                case_name = ['aligned', 'misaligned'][i_disc] + '_M_{}_N_{}_Mstar_{}'.format(\
                        params['M'], params['N'], "{:.1f}".format(params['Msf']))
                case_name = case_name.replace('.', '_')
                print('\n' + case_name)

                # Copy previous case data
                case_data_prev = case_data_curr

                
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
                                        M = params['M'],
                                        N = params['N'],
                                        Mstar_fact = params['Msf'],
                                        sym = False,
                                        wake_cfl1 = False,
                                        write_screen = 'off')
                
                case_data_curr = sharpy.sharpy_main.main(['', wing.route + wing.case_name + '.sharpy'])
                
                # biot_savart.plot_points(case_data_curr, pnts, 0)

                # Check for initial or final case
                if i_case == 0:
                        continue
                elif is_final or i_case + 1 == max_iter:
                        case_data_final.update({['aligned', 'misaligned'][i_disc] : case_data_curr})
                        wing_final.update({['aligned', 'misaligned'][i_disc] : wing})
                        disc_final.update({['aligned', 'misaligned'][i_disc] : params})
                        break

                # Test for convergence
                changed_param = convergence.curr_param(disc_order, i_case)
                if convergence.is_converged(case_data_prev, case_data_curr, pnts, False):
                        disc_dict[changed_param]['conv'] = True
                else:
                        disc_dict = convergence.disc_update(disc_dict, changed_param)
                        disc_dict[changed_param]['conv'] = False

                # Check if three converged cases in a row are run
                for i_param in disc_order:
                        if disc_dict[i_param]['conv'] == False:
                                n_conv = -1
                
                # Update counts if converged
                n_conv += 1
                if n_conv == 3:
                        is_final = True         

### End of loop
                                
### Case comparision
# Will always return True, prints values to screen
convergence.is_converged(case_data_final['aligned'], case_data_final['misaligned'], pnts, True)

# Extract key data from system and output to matlab
case_data_out = dict()
case_data_out.update({'aligned': case_data_extract.case_data_extract(\
        wing_final['aligned'], case_data_final['aligned'])})
case_data_out.update({'misaligned': case_data_extract.case_data_extract(\
        wing_final['misaligned'], case_data_final['misaligned'])})
case_data_out.update({'disc_final': disc_final})
spio.savemat('convergence_out.mat', case_data_out)