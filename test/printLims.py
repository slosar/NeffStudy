#!/usr/bin/env python
##
## print limits
##
from __future__ import division, print_function
import sys,os
import numpy as np
sys.path=["py"]+sys.path
from S4Fisher import S4Fisher
from TracerPk import TracerPk
from CVFisher import CVFisher

F=S4Fisher()
F2=S4Fisher()
#U=TracerPk()
U=CVFisher(snrboost=1.0,zmin=1.0)
F2.addF(U)
for n in U.plist.nameList():
    print("%s sigma=%g %g %g"%(n,F.error(n),F2.error(n),U.error(n)))

#import matplotlib.pyplot as plt
#di=np.sqrt(U.F.diagonal())
#Co=U.F/np.outer(di,di)
#plt.imshow(Co,interpolation='nearest')
#plt.colorbar()
#plt.show()
          
