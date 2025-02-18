

import numpy as np
from interferopy.cube import Cube


def findclump(c,out):
    cube = Cube(c)
    catP, catN, candP, candN = cube.findclumps_full(output_file=out,
                                                    kernels=np.arange(3,20,2),
                                                    run_search=True,
                                                    run_crop=True,
                                                    SNR_min=3.0,
                                                    delta_offset_arcsec=2,
                                                    delta_freq=0.1,
                                                    run_fidelity=True,
                                                    fidelity_bins=np.arange(0,15,0.2),
                                                    min_SN_fit=4.0,fidelity_threshold=0.8,
                                                    verbose=True,
                                                    sextractor_param_file='',
                                                    ncores=1)