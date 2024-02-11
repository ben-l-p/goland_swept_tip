import shutil
import os
import warnings
import numpy as np
import matplotlib.pyplot as plt
import scipy.io as spio

import sharpy.sharpy_main

from case_data_extract import case_data_extract
from wing_generator import swept_tip_goland

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
        'Modal',
        'AerogridPlot',
        'BeamPlot',
        'DynamicCoupled',
        'AeroForcesCalculator',
        'LinearAssembler',
        'AsymptoticStability',
        ]

### Loop this code, varying parameters each run
case_name = 'test_case'
wing = swept_tip_goland(case_name, flow, np.deg2rad(45), 0.5, 2, 1)
case_data = sharpy.sharpy_main.main(['', wing.route + wing.case_name + '.sharpy'])
case_output.update({wing.case_name: case_data_extract(wing, case_data)})
### End of loop

spio.savemat('case_output.mat', case_output)