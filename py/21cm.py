#
# Details for a generalized 21cm experiment, similar to BMX.
#

from __future__ import division, print_function
import sys
import os.path
from TracerPk import TracerPk
import numpy as np
import h5py

class BMX(TracerPk):

    def __init__(self):
        self.da=h5py.File('sn_lowz_expA_50K.h5', 'r') # Data cube from Cosmic Visions
        zvals=self.zvals()
        min_z=zvals[0]
        max_z=zvals[-1]
        SNR=self.SNR()
        TracerPk.__init__(self, zvals, min_z, max_z, '21cm', SNR)

    def zvals(self):
        """Redshift values in data cube.
        """
        zs=[]
        for i in range(9):
            czs=self.da['sn_band_'+str(i)].attrs['z']
            zs.append(czs)
        zs=np.array(zs)
        zs=zs[::-1]
        return zs

    def SNR(self):
        """Signal-to-noise ratio values in data cube.
        """
        snr=[]
        for i in range(9):
            csnr=self.da['sn_band_'+str(i)].value
            snr.append(csnr)
        snr=np.array(snr)
        return snr

# Create an instance of this class.
BMX()
