#!/usr/bin/env python
##
## print limits
##
from __future__ import division, print_function
import sys,os
import numpy as np
sys.path=["py"]+sys.path
from S4Fisher import S4Fisher
from UBAO import UBFisher
from UCVLim import UCVLim
F=S4Fisher()
F2=S4Fisher()
#U=UBFisher(6,fname="FishData/desilz_14000_0.65_1.9_0.1_1_ELG1440_3_bao_0.expt")
#U=UBFisher(6)
U=UCVLim()
F2.addF(U)
for n in U.plist:
    print("%s sigma=%g %g %g"%(n,F.Error(n),F2.Error(n),U.Error(n)))

#import matplotlib.pyplot as plt
#di=np.sqrt(U.F.diagonal())
#Co=U.F/np.outer(di,di)
#plt.imshow(Co,interpolation='nearest')
#plt.colorbar()
#plt.show()
          
