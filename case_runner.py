import numpy as np
import matplotlib.pyplot as plt
import os
import sys
import meshio
import fnmatch
import scipy.io as spio
import pandas

import flying_wings as wings                    #Goland wing defintion
import sharpy.sharpy_main 

import warnings
warnings.filterwarnings("ignore")               #Disable warnings to make output more readable

### Setup parameters
ang_h = np.deg2rad(np.linspace(-50, 50, 11))    #Tip sweep angles in degrees. Base case of (-50, 50, 11)
pos_frac_h = np.linspace(0.5, 0.95, 10)         #Tip sweep position in fraction of chord. Base case of (0.5, 0.95, 10)

u_inf = 40.                 # Velocity for static analysis
alpha_deg = 2.              # Define angle of attack for static aeroelastic analsis
rho = 1.02                  # Air density
M = 16                      # Number of chordwise panels
N = 20                      # Number of spanwise panels
M_star_fact = 10            # Length of the wake in chords.
c_ref = 1.8288              # Goland wing reference chord
num_modes =  8              # Number of vibration modes retained in the structural model.
n_vel_out = 500             # Number of velocities to use for stability analysis. Limits are defined in AsymptoticStability
n_surfaces = 1              # Number of wings
physical_time = 4          # Simulation runtime for dynamic coupled

### Output array allocation
v_flutter = np.zeros([len(ang_h), len(pos_frac_h)])
v_flutter_filt = np.zeros([len(ang_h), len(pos_frac_h)])
v_div = np.zeros([len(ang_h), len(pos_frac_h)])
v_max = np.zeros([len(ang_h), len(pos_frac_h)])
ang_flutter = np.zeros([len(ang_h), len(pos_frac_h)])
pos_flutter = np.zeros([len(ang_h), len(pos_frac_h)])
unst_modes = np.zeros([len(ang_h), len(pos_frac_h), n_vel_out])
moments_g = np.zeros([len(ang_h), len(pos_frac_h), 3])
moments_a = np.zeros([len(ang_h), len(pos_frac_h), 3])
forces_g = np.zeros([len(ang_h), len(pos_frac_h), 3])
forces_a = np.zeros([len(ang_h), len(pos_frac_h), 3])

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

for i in range(len(ang_h)):                 # Loop through sweep angles
    for j in range(len(pos_frac_h)):        # Loop through sweep positions

        print("\nCase ", j+i*len(pos_frac_h)+1, " of ", len(ang_h)*len(pos_frac_h))
        print("Angle: ", np.rad2deg(ang_h[i]), "\nPosition: ", pos_frac_h[j])

        ang_flutter[i, j] = ang_h[i]
        pos_flutter[i, j] = pos_frac_h[j]

        case_name = 'goland_ang_{}_pos_{}'.format("{:.2f}".format(np.rad2deg(ang_h[i])), "{:.2f}".format(pos_frac_h[j]))
        case_name = case_name.replace('.', '_')
        case_name = case_name.replace('-', 'n')  
        
        ### Generate Goland wing object
        ws = wings.GolandControlSurface(M=M,
                                        N=N,
                                        Mstar_fact=M_star_fact,
                                        u_inf=u_inf,
                                        alpha=alpha_deg,
                                        cs_deflection=[0, 0],
                                        rho=rho,
                                        sweep=0,
                                        ang_h=ang_h[i],
                                        pos_frac_h=pos_frac_h[j],
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
            'flow':
                ['BeamLoader', 'AerogridLoader',
                'StaticCoupled',
                'Modal',
                'AerogridPlot',
                'BeamPlot',
                'LinearAssembler',
                'AsymptoticStability',
                'DynamicCoupled',
                'AeroForcesCalculator'
                ],
            'case': ws.case_name, 'route': ws.route,
            'write_screen': 'on', 'write_log': 'on',       # Change to on to see console output
            'log_folder': route_test_dir + '/output/',
            'log_file': ws.case_name + '.log',
            'route': route_test_dir + '/cases/'}

        ws.config['BeamLoader'] = {
            'unsteady': 'off',
            'orientation': ws.quat}

        ws.config['AerogridLoader'] = {
            'unsteady': 'off',
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
                'rho': ws.rho,
                'print_info': 'off',
                'horseshoe': 'off',
                'num_cores': 4,
                'n_rollup': 100,                #0
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
                                            'convection_scheme': 2,
                                            'gamma_dot_filtering': 6,
                                            'velocity_field_generator': 'GustVelocityField',
                                            'velocity_field_input': {'u_inf': u_inf,
                                                                    'u_inf_direction': [1., 0, 0],
                                                                    'gust_shape': '1-cos',
                                                                    'gust_parameters': {'gust_length': 2,
                                                                                        'gust_intensity': 1 * u_inf},
                                                                    'offset': -2,
                                                                    'relative_motion': 0},
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
                                            'velocity_analysis': [1, 500, n_vel_out],
                                        'modes_to_plot': []}
        
        ws.config['AeroForcesCalculator'] = {'write_text_file': True,
                                             'screen_output': False,
                                             'q_ref':1/2*rho*u_inf**2,
                                             'b_ref': 2. * 6.096,
                                             'c_ref': c_ref,
                                             'S_ref': c_ref*2. * 6.096}

        ws.config.write()

        [case_data, prof_data] = sharpy.sharpy_main.main(['', ws.route + ws.case_name + '.sharpy'])     # Run SHARPy

        # Forces and Moments
        force_data = pandas.read_csv('./output/%s/forces/forces_aeroforces.txt' % case_name, delimiter=', ', index_col=False).to_dict()
        moment_data = pandas.read_csv('./output/%s/forces/moments_aeroforces.txt' % case_name, delimiter=', ', index_col=False).to_dict()

        forces_a[i, j, :] = [force_data['fx_steady_a'][0], force_data['fy_steady_a'][0], force_data['fz_steady_a'][0]]
        forces_g[i, j, :] = [force_data['fx_steady_G'][0], force_data['fy_steady_G'][0], force_data['fz_steady_G'][0]]
        moments_a[i, j, :] = [moment_data['mx_steady_a'][0], moment_data['my_steady_a'][0], moment_data['mz_steady_a'][0]]
        moments_g[i, j, :] = [moment_data['mx_steady_G'][0], moment_data['my_steady_G'][0], moment_data['mz_steady_G'][0]]
        
        # Stability
        for file_name in os.listdir('./output/%s/stability/' % case_name):
            if(fnmatch.fnmatch(file_name, '*.dat')):
                velocity_analysis = np.loadtxt('./output/%s/stability/%s' % (case_name, file_name))

        u_inf_out = velocity_analysis[:, 0]
        eigs_r = velocity_analysis[:, 1]
        eigs_i = velocity_analysis[:, 2]

        n_modes = int(len(u_inf_out)/n_vel_out)

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
        for k in range(n_vel_out):
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
        for k in range(n_vel_out):                                                          
            n_unst = 0
            for l in range(n_modes):
                if eigs_r[k*n_modes + l] >= 0 and abs(eigs_i[k*n_modes + l]) >= im_div_t:
                    n_unst += 1
            if n_unst > n_unst_init_f:
                v_flutter[i, j] = u_inf_out[k*n_modes + l]
                break
        
        # Flutter (filtered)
        # Ignore low velocity results. Same effect as changing the lower limit for velocity sweep
        for k in range(n_vel_out):
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
        for k in range(n_vel_out):
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

spio.savemat('case_output.mat', {"moments_a": moments_a, "moments_g": moments_g, "forces_a": forces_a, "forces_g": forces_g,\
                            "n_unst_modes": unst_modes, "v_flutter": v_flutter, "v_flutter_filt": v_flutter_filt, "v_div": v_div,\
                                  "ang_flutter": ang_flutter, "pos_flutter": pos_flutter})
