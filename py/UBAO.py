#
# Ultimate BAO matrix.
#

from __future__ import division, print_function
import sys
sys.path=['/home/anze/work/April/py']+sys.path
from LCDMCosmology import LCDMCosmology
#LCDMCosmology.rd_approx=147.238723337
LCDMCosmology.rd_approx="CuestaNeff"
import numpy as np
import scipy.linalg as la
from FishMat import FishMat

class UBFisher(FishMat):

    def __init__ (self, zmax=6.,fname="FishData/ultimate_41253_0.1_10_0.1_1_bao_0.expt"):
        da=open(fname).readlines()
        Nl=int(da[0].split()[-1])
        da=np.loadtxt(da[2:2+Nl])
        ## name, use, values, delta
        nal=[('tau', False, 0.07),
            ('omegac',True,0.11987),
            ('As',False,0),
            ('theta',True,1.0407781e-2),
            ('Neff',True,3.04),
            ('mnu',True, 0.06),
            ('omegab',True,0.022252),
            ('ns',False,0.96)]
        N=len(nal)
        F=np.zeros((N,N))
        Ce=self.getCosmo(nal,-1,0)
        eps=0.002
        for i1,(n1,u1,v1) in enumerate(nal):
            if not u1:
                continue
            Cep1=self.getCosmo(nal,i1,+eps)
            Cem1=self.getCosmo(nal,i1,-eps)
            for i2,(n2,u2,v2) in enumerate(nal):
                if (i2<i1):
                    continue
                if not u2:
                    continue
                print(i1,n1,i2,n2)
                if (i1!=i2):
                    Cep2=self.getCosmo(nal,i2,+eps)
                    Cem2=self.getCosmo(nal,i2,-eps)
                else:
                    Cep2=Cep1
                    Cem2=Cem1
                for z,DaE,HE,r in  da[:,:4]:
                    if (z>zmax):
                        continue
                    DaE*=0.01
                    HE*=0.01
                    m=np.array([Ce.DaOverrd(z),Ce.HIOverrd(z)])
                    Em=np.outer(m,m)*np.array([[DaE**2,DaE*HE*r],[DaE*HE*r,HE*HE]]) ## this is element wise mult
                    EmI=la.inv(Em)
                    D1=(np.array([Cep1.DaOverrd(z),Cep1.HIOverrd(z)])
                        -np.array([Cem1.DaOverrd(z),Cem1.HIOverrd(z)]))/(2*v1*eps) 
                    D2=(np.array([Cep2.DaOverrd(z),Cep2.HIOverrd(z)])
                        -np.array([Cem2.DaOverrd(z),Cem2.HIOverrd(z)]))/(2*v2*eps)
                    F[i1,i2]+=np.dot(np.dot(EmI,D1),D2)
        for i in range(N):
            for j in range(i,N):
                F[j,i]=F[i,j]
            
        print(F)
        FishMat.__init__(self,[x[0] for x in nal],F)
        
                
    
                
    def getCosmo(self,nal,p,o):
        omegac=nal[1][2]
        daord_cmb=0.01/(nal[3][2]*0.010601/1.04077)
        Neff=nal[4][2]
        mnu=nal[5][2]
        omegab=nal[6][2]
        if (p==1):
            omegac*=(1.+o)
        elif (p==3):
            daord_cmb*=(1+o/10.)
        elif (p==4):
            Neff*=(1+o)
        elif (p==5):
            mnu*=(1+o)
        elif (p==6):
            omegab*=(1+o)
        elif (p<0):
            pass
        else:
            print ("WTF?")
            stop()
        ## now need to find h
        hlow=0.6
        tlow=LCDMCosmology(omegab,(omegab+omegac+mnu/94.)/hlow**2,hlow,mnu,Neff).DaOverrd(1090)
        hhigh=0.8
        thigh=LCDMCosmology(omegab,(omegab+omegac+mnu/94.)/hhigh**2,hhigh,mnu,Neff).DaOverrd(1090)
        #print (tlow,thigh,daord_cmb)
        print (p,o,'XX')
        while True:
            #print (hlow,hhigh)
            hcen=(hlow+hhigh)/2
            tcen=LCDMCosmology(omegab,(omegab+omegac+mnu/94.)/hcen**2,hcen,mnu,Neff).DaOverrd(1090)
            if (tcen>daord_cmb):
                tlow=tcen
                hlow=hcen
            else:
                thigh=tcen
                hhigh=hcen
            if abs(hhigh-hlow)<1e-6:
                if (hcen<0.65) or (hcen>0.75):
                    print("SOMETHING BAD")
                    stop()
                break
        return LCDMCosmology(omegab,(omegab+omegac+mnu/94.)/hcen**2,hcen,mnu,Neff)

    


