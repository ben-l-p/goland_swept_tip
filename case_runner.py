import numpy as np
import matplotlib.pyplot as plt
import os
import fnmatch
import scipy.io as spio
import pandas

import flying_wings as wings                    #Goland wing defintion
import sharpy.sharpy_main 

import warnings
warnings.filterwarnings("ignore")               #Disable warnings to make output more readable

### Setup parameters
# ang_h = np.deg2rad(np.linspace(-50, 50, 11))      # Tip sweep angles in radians. Base case of (-50, 50, 11)
# pos_frac_h = np.linspace(0.5, 0.95, 10)           # Fraction of total straight span between root and hinge. Base case of (0.5, 0.95, 10)

ang_h = np.deg2rad(33)      # Tip sweep angles in radians
pos_frac_h = 0.69           # Fraction of total straight span between root and hinge
n_surfaces = 2              # Number of wings (1 or 2)

u_inf = 80.                 # Velocity for static analysis
alpha_deg = 3.              # Define angle of attack for static aeroelastic analsis
rho = 1.02                  # Air density
M = 16                      # Number of chordwise panels
N = 20*n_surfaces           # Number of spanwise panels
M_star_fact = 10            # Length of the wake in chords.
c_ref = 1.8288              # Goland wing reference chord
num_modes =  8              # Number of vibration modes retained in the structural model.

asym_v_min = 50             # Minimum velocity for stability analysis (m/s)
asym_v_max = 250            # Maximum velocity for stability analysis (m/s)
asym_v_num = 101            # Number of velocities to use for stability analysis

freq_min = 1e-2             # Minimum frequency for frequency analysis (rad/s (w), reduced (k))
freq_max = 1e2              # Maximum frequency for frequency analysis (rad/s (w), reduced (k))
freq_num = 150              # Number of frequencies to use for frequency analysis

physical_time = 0.6         # Simulation runtime for dynamic coupled (s)

gust_intensity = 0.5
gust_length = 0.2 * u_inf
gust_offset = 0.05 * u_inf

dt = c_ref / M / u_inf
n_tstep = int(np.round(physical_time / dt))

flow =  ['BeamLoader', 
        'AerogridLoader',
        'StaticCoupled',
        'Modal',
        'AerogridPlot',
        'BeamPlot',
        'LinearAssembler',
        'AsymptoticStability',
        'FrequencyResponse',
        # 'DynamicCoupled',
        'AeroForcesCalculator',
        ]

### Prevent linearising problem if only a single wing (throws errors due to asymmetry)
try:
    if n_surfaces == 1:     
        flow.remove('LinearAssembler')
        flow.remove('AsyptoticStability')
        flow.remove('FrequencyResponse')
except: 
    pass

### Allows for running a single case
try:
    n_angs = len(ang_h)
except:
    n_angs = 1

try:
    n_pos_h = len(pos_frac_h)
except:
    n_pos_h = 1

### Output array allocation
v_flutter = np.zeros([n_angs, n_pos_h])
v_flutter_filt = np.zeros([n_angs, n_pos_h])
v_div = np.zeros([n_angs, n_pos_h])
v_max = np.zeros([n_angs, n_pos_h])
ang_flutter = np.zeros([n_angs, n_pos_h])
pos_flutter = np.zeros([n_angs, n_pos_h])
unst_modes = np.zeros([n_angs, n_pos_h, asym_v_num])

n_force_steps = 1
if 'DynamicCoupled' in flow:
    n_force_steps = n_tstep

moments_g_s = np.zeros([n_angs, n_pos_h, n_force_steps, 3])
moments_a_s = np.zeros([n_angs, n_pos_h, n_force_steps, 3])
forces_g_s = np.zeros([n_angs, n_pos_h, n_force_steps, 3])
forces_a_s = np.zeros([n_angs, n_pos_h, n_force_steps, 3])
moments_g_u = np.zeros([n_angs, n_pos_h, n_force_steps, 3])
moments_a_u = np.zeros([n_angs, n_pos_h, n_force_steps, 3])
forces_g_u= np.zeros([n_angs, n_pos_h, n_force_steps, 3])

forces_a_u = np.zeros([n_angs, n_pos_h, n_force_steps, 2])  #No Z term

### Universal case setup
rom_settings = dict()
rom_settings['algorithm'] = 'mimo_rational_arnoldi'                     # Reduction algorithm
rom_settings['r'] = 6                                                   # Krylov subspace order
frequency_continuous_k = np.array([0.])                                 # Interpolation point in the complex plane with reduced frequency units
frequency_continuous_w = 2 * u_inf * frequency_continuous_k / c_ref
rom_settings['frequency'] = frequency_continuous_w

case_name = 'goland_cs'
case_nlin_info = 'M%dN%dMs%d_nmodes%d' % (M, N, M_star_fact, num_modes)
case_rom_info = 'rom_MIMORA_r%d_sig%04d_%04dj' % (rom_settings['r'], frequency_continuous_k[-1].real * 100,
                                                  frequency_continuous_k[-1].imag * 100)
route_test_dir = os.path.dirname(os.path.realpath(__file__))

for i in range(n_angs):                 # Loop through sweep angles
    for j in range(n_pos_h):        # Loop through sweep positions

        try:
            i_ang = ang_h[i]
        except:
            i_ang = ang_h

        try:
            j_pos_h = pos_frac_h[j]
        except:
            j_pos_h = pos_frac_h

        print("\nCase ", j+i*n_pos_h+1, " of ", n_angs*n_pos_h)
        print("Angle: ", np.rad2deg(i_ang), "\nPosition: ", j_pos_h)

        ang_flutter[i, j] = i_ang
        pos_flutter[i, j] = j_pos_h

        case_name = 'goland_ang_{}_pos_{}'.format("{:.2f}".format(np.rad2deg(i_ang)), "{:.2f}".format(j_pos_h))
        case_name = case_name.replace('.', '_')
        case_name = case_name.replace('-', 'n')  
        
        ### Generate Goland wing object
        ws = wings.GolandControlSurface(M=M,
                                        N=N,
                                        Mstar_fact=M_star_fact,
                                        u_inf=u_inf,
                                        alpha=alpha_deg,
                                        n_control_surfaces=0,
                                        cs_deflection=[],
                                        rho=rho,
                                        sweep=0,
                                        ang_h=i_ang,
                                        pos_frac_h=j_pos_h,
                                        physical_time=physical_time,
                                        n_surfaces=n_surfaces,
                                        route=route_test_dir + '/cases',
                                        case_name=case_name)

        ws.clean_test_files()
        ws.update_derived_params()
        ws.set_default_config_dict()

        ws.generate_aero_file()
        ws.generate_fem_file()

        ### SHARPy case settings
        ws.config['SHARPy'] = {
            'flow': flow,
            'case': ws.case_name, 'route': ws.route,
            'write_screen': 'on', 'write_log': 'on',       # Change to on to see console output
            'log_folder': route_test_dir + '/output/',
            'log_file': ws.case_name + '.log',
            'route': route_test_dir + '/cases/'}

        ws.config['BeamLoader'] = {
            'unsteady': 'on',
            'orientation': ws.quat}

        ws.config['AerogridLoader'] = {
            'unsteady': 'on',
            'aligned_grid': 'off',                          # Changed from original case
            'mstar': ws.Mstar_fact * ws.M,
            'freestream_dir': ws.u_inf_direction,
            'wake_shape_generator': 'StraightWake',
            'wake_shape_generator_input': {'u_inf': ws.u_inf,
                                            'u_inf_direction': ws.u_inf_direction,
                                            'dt': ws.dt}}

        ws.config['StaticCoupled'] = {
            'print_info': 'on',
            'max_iter': 200,
            'n_load_steps': 10,                 #1
            'tolerance': 1e-10,
            'relaxation_factor': 0.35,            #0.
            'aero_solver': 'StaticUvlm',
            'aero_solver_settings': {
                'symmetry_condition': bool(n_surfaces % 2),
                'symmetry_plane': 1,
                'rho': ws.rho,
                'print_info': 'off',
                'horseshoe': 'off',
                'num_cores': 4,
                'n_rollup': 0,                #0
                'rollup_dt': ws.dt,
                'rollup_aic_refresh': 1,
                'rollup_tolerance': 1e-4,
                'velocity_field_generator': 'SteadyVelocityField',
                'velocity_field_input': {
                    'u_inf': ws.u_inf,
                    'u_inf_direction': ws.u_inf_direction}},
            'structural_solver': 'NonLinearStatic',
            'structural_solver_settings': {'print_info': 'off',
                                        'max_iterations': 150,
                                        'num_load_steps': 4,
                                        'delta_curved': 1e-1,
                                        'min_delta': 1e-10,
                                        'gravity_on': 'on',
                                        'gravity': 9.81}}

        ws.config['DynamicCoupled'] = {'structural_solver': 'NonLinearDynamicPrescribedStep',
                                'structural_solver_settings': {'print_info': 'off',
                                            'max_iterations': 950,
                                            'delta_curved': 1e-1,
                                            'newmark_damp': 5e-3,
                                            'gravity_on': True,
                                            'gravity': 9.81,
                                            'num_steps': ws.n_tstep,
                                            # 'initial_velocity': ws.u_inf,
                                            'dt': ws.dt},
                                'aero_solver': 'StepUvlm',
                                'aero_solver_settings': {'print_info': 'off',
                                            'num_cores': 4,
                                            'symmetry_condition': bool(n_surfaces % 2),
                                            'symmetry_plane': 1,
                                            'convection_scheme': 2,
                                            'gamma_dot_filtering': 6,
                                            'velocity_field_generator': 'GustVelocityField',
                                            'velocity_field_input': {'u_inf': u_inf,
                                                                    'u_inf_direction': [1., 0, 0],
                                                                    'gust_shape': '1-cos',
                                                                    'gust_parameters': {'gust_length': gust_length,
                                                                                        'gust_intensity': gust_intensity * u_inf},
                                                                    'offset': gust_offset,
                                                                    'relative_motion': 'on'},
                            'n_time_steps': ws.n_tstep,
                            'dt': ws.dt,
                            'cfl1': True},
                                'fsi_substeps': 200,
                                'minimum_steps': 1,
                                'relaxation_steps': 150,
                                'final_relaxation_factor': 0.5,
                                'n_time_steps': ws.n_tstep,
                                'dt': ws.dt,
                                'include_unsteady_force_contribution': 'on',
                                'postprocessors': ['BeamLoads', 'BeamPlot', 'AerogridPlot'],
                                'postprocessors_settings': {'BeamLoads': {'csv_output': 'off'},
                                                            'BeamPlot': {'include_rbm': 'on',
                                                                        'include_applied_forces': 'on'},
                                                            'AerogridPlot': {
                                                                'include_rbm': 'on',
                                                                'include_applied_forces': 'on',
                                                                'minus_m_star': 0},
                                                            }}

        ws.config['AerogridPlot'] = {'include_rbm': 'off',
                                    'include_applied_forces': 'on',
                                    'minus_m_star': 0}

        ws.config['BeamPlot'] = {'include_rbm': 'off',
                                'include_applied_forces': 'on'}

        ws.config['Modal'] = {'NumLambda': 20,
                            'rigid_body_modes': 'off',
                            'print_matrices': 'on', # 'keep_linear_matrices': 'on',            'write_dat': 'off',
                            'rigid_modes_cg': 'off',
                            'continuous_eigenvalues': 'off',
                            'dt': 0,
                            'plot_eigenvalues': False,
                            'max_rotation_deg': 15.,
                            'max_displacement': 0.15,
                            'write_modes_vtk': True,
                            'use_undamped_modes': True}

        ws.config['LinearAssembler'] = {'linear_system': 'LinearAeroelastic',
                                        'linear_system_settings': {
                                            'beam_settings': {'modal_projection': 'on',
                                                            'inout_coords': 'modes',
                                                            'discrete_time': 'on',
                                                            'newmark_damp': 0.5e-4,
                                                            'discr_method': 'newmark',
                                                            'dt': ws.dt,
                                                            'proj_modes': 'undamped',
                                                            'use_euler': 'off',
                                                            'num_modes': num_modes,
                                                            'print_info': 'on',
                                                            'gravity': 'on',
                                                            'remove_sym_modes': 'on',
                                                            'remove_dofs': []},
                                            'aero_settings': {'dt': ws.dt,
                                                            'ScalingDict': {'length': 0.5 * ws.c_ref,
                                                                            'speed': u_inf,
                                                                            'density': rho},
                                                            'integr_order': 2,
                                                            'density': ws.rho,
                                                            'remove_predictor': 'on',
                                                            'use_sparse': 'on',
                                                            'remove_inputs': ['u_gust'],
                                                            'rom_method': ['Krylov'],
                                                            'rom_method_settings': {'Krylov': rom_settings}},
                                            }}

        ws.config['AsymptoticStability'] = {'print_info': True,
                                            'velocity_analysis': [asym_v_min, asym_v_max, asym_v_num],
                                        'modes_to_plot': []}
        
        ws.config['FrequencyResponse'] = {'quick_plot': True,       # Saves a plot as an image for every input/output combination
                                          'frequency_unit': 'w',
                                          'num_freqs': freq_num,
                                          'frequency_bounds': [freq_min, freq_max],
                                          'frequency_spacing': 'log',
                                          }
        
        ws.config['AeroForcesCalculator'] = {'write_text_file': True,
                                             'screen_output': False}

        ws.config.write()

        case_data = sharpy.sharpy_main.main(['', ws.route + ws.case_name + '.sharpy'])     # Run SHARPy

        # Forces and Moments
        force_data = pandas.read_csv('./output/%s/forces/forces_aeroforces.txt' % case_name, delimiter=', ', index_col=False).to_dict()
        moment_data = pandas.read_csv('./output/%s/forces/moments_aeroforces.txt' % case_name, delimiter=', ', index_col=False).to_dict()

        for k in range(n_force_steps):
            forces_a_s[i, j, k, :] = [force_data['fx_steady_a'][k], force_data['fy_steady_a'][k], force_data['fz_steady_a'][k]]
            forces_g_s[i, j, k, :] = [force_data['fx_steady_G'][k], force_data['fy_steady_G'][k], force_data['fz_steady_G'][k]]
            moments_a_s[i, j, k, :] = [moment_data['mx_steady_a'][k], moment_data['my_steady_a'][k], moment_data['mz_steady_a'][k]]
            moments_g_s[i, j, k, :] = [moment_data['mx_steady_G'][k], moment_data['my_steady_G'][k], moment_data['mz_steady_G'][k]]
            forces_a_u[i, j, k, :] = [force_data['fx_unsteady_a'][k], force_data['fy_unsteady_a'][k]]
            forces_g_u[i, j, k, :] = [force_data['fx_unsteady_G'][k], force_data['fy_unsteady_G'][k], force_data['fz_unsteady_G'][k]]
            moments_a_u[i, j, k, :] = [moment_data['mx_unsteady_a'][k], moment_data['my_unsteady_a'][k], moment_data['mz_unsteady_a'][k]]
            moments_g_u[i, j, k, :] = [moment_data['mx_unsteady_G'][k], moment_data['my_unsteady_G'][k], moment_data['mz_unsteady_G'][k]]
        
        # Stability
        if 'AsymptoticStability' in ws.config['SHARPy']['flow']:
            for file_name in os.listdir('./output/%s/stability/' % case_name):
                if(fnmatch.fnmatch(file_name, '*.dat')):
                    velocity_analysis = np.loadtxt('./output/%s/stability/%s' % (case_name, file_name))

            u_inf_out = velocity_analysis[:, 0]
            eigs_r = velocity_analysis[:, 1]
            eigs_i = velocity_analysis[:, 2]

            n_modes = int(len(u_inf_out)/asym_v_num)

            n_unst_init_f = 0
            n_unst_init_d = 0

            im_div_t = 1e-5
            skip_n_vel = 70

            for k in range(n_modes):
                if eigs_r[k] >= 0:
                    if abs(eigs_i[k]) <= im_div_t:
                        n_unst_init_d += 1
                    else:
                        n_unst_init_f += 1
            
            # Divergence
            # Filtered to remove instabilities occuring at low velocity due to being far from linearisation point
            for k in range(asym_v_num):
                if k < skip_n_vel:
                    continue
                n_unst = 0
                for l in range(n_modes):
                    if eigs_r[k*n_modes + l] >= 0 and abs(eigs_i[k*n_modes + l]) <= im_div_t:
                        n_unst += 1
                if n_unst > n_unst_init_d:
                    v_div[i, j] = u_inf_out[k*n_modes + l]
                    break
            print("Divergence: ", v_div[i, j], " m/s")

            # Flutter (raw)
            # Error prone when using velocities far beyond what was linearised for
            # This can lead to random flutter points or a very low divergence velocity
            for k in range(asym_v_num):                                                          
                n_unst = 0
                for l in range(n_modes):
                    if eigs_r[k*n_modes + l] >= 0 and abs(eigs_i[k*n_modes + l]) >= im_div_t:
                        n_unst += 1
                if n_unst > n_unst_init_f:
                    v_flutter[i, j] = u_inf_out[k*n_modes + l]
                    break
            
            # Flutter (filtered)
            # Ignore low velocity results. Same effect as changing the lower limit for velocity sweep
            for k in range(asym_v_num):
                if k < skip_n_vel:
                    continue
                n_unst = 0
                for l in range(n_modes):
                    if eigs_r[k*n_modes + l] >= 0 and abs(eigs_i[k*n_modes + l]) >= im_div_t:
                        n_unst += 1
                if n_unst > n_unst_init_f:
                    v_flutter_filt[i, j] = u_inf_out[k*n_modes + l]
                    break
            print("Flutter: ", v_flutter[i, j], " m/s")

            # Number of unstable modes
            # Will include instability errors as with the raw flutter/divergence velocity
            for k in range(asym_v_num):
                for l in range(n_modes):
                    if eigs_r[k*n_modes + l] >= 0:
                        unst_modes[i, j, k] += 1

        # fig = plt.figure()
        # plt.scatter(eigs_r, eigs_i, c=u_inf_out, cmap='Blues')
        # cbar = plt.colorbar()
        # cbar.set_label('Free Stream Velocity, $u_\infty$ [m/s]')

        # plt.grid()
        # plt.xlim(-50, 5)
        # plt.ylim(0, 300)
        # plt.xlabel('Real Part, $\lambda$ [rad/s]')
        # plt.ylabel('Imag Part, $\lambda$ [rad/s]')
        # plt.show()

spio.savemat('case_output.mat', {"moments_a_s": moments_a_s, "moments_g_s": moments_g_s, "forces_a_s": forces_a_s, "forces_g_s": forces_g_s,\
                                 "moments_a_u": moments_a_u, "moments_g_u": moments_g_u, "forces_a_u": forces_a_u, "forces_g_u": forces_g_u,\
                            "n_unst_modes": unst_modes, "v_flutter": v_flutter, "v_flutter_filt": v_flutter_filt, "v_div": v_div,\
                                  "ang_flutter": ang_flutter, "pos_flutter": pos_flutter, "dt": dt})
