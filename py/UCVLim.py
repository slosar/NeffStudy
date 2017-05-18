
##
## Ultimate BAO matrix
##

from __future__ import division, print_function
import sys
sys.path=['/home/anze/work/April/py']+sys.path
from LCDMCosmology import LCDMCosmology
#LCDMCosmology.rd_approx=147.238723337
LCDMCosmology.rd_approx="CuestaNeff"
import numpy as np
import scipy.linalg as la
from scipy.interpolate import interp1d
from FishMat import FishMat
import classy
print (classy.__file__)
from classy import Class


class UCVLim(FishMat):
    def __init__ (self, zmin=2, zmax=6, dz=0.1, kmax=2.0, kzscal=0.1,dk=0.01,fsky=0.15):
        nal=[('tau', False, 0.07),
            ('omegac',True,0.11987),
            ('As',False,1.0),
            ('theta',True,1.0407781e-2),
            ('Neff',True,3.04),
            ('mnu',True, 0.06),
            ('omegab',True,0.022252),
            ('ns',True,0.96475),
            ('bd',True,1.), ## we treat them here as one, but will marg z by z later
            ('be',True,1.)
        ]
        N=len(nal)
        F=np.zeros((N-2,N-2))
        self.kval=np.arange(dk/2,kmax,dk)
        self.Nk=len(self.kval)
        self.dk=dk
        self.dz=dz
        self.fsky=fsky
        self.kpar=np.outer(self.kval,np.ones(self.Nk))
        self.kperp=self.kpar.T
        self.kt=np.sqrt(self.kpar**2+self.kperp**2)
        self.mu=self.kpar/self.kt
        self.zval=np.arange(zmin+dz/2,zmax,dz)
        PkC=self.getCube(nal,None,0)
        PkEI=[]
        for z,P,nm in zip(self.zval,PkC,self.nmodes):
            PkE=P*P/nm
            knl= 0.04+0.016*(1+z)**2.2
            PkE[np.where(self.kt>knl)]=1e30
            PkEI.append(1/PkE)
        
        eps=0.002

        Pderivs=[]
        for i1,(n1,u1,v1) in enumerate(nal):
            print ("DER:",n1)
            if not u1:
                Ders=None
            else:
                PkP=self.getCube(nal,n1,+eps)
                PkM=self.getCube(nal,n1,-eps)
                Ders=[(P-M)/(2*v1*eps) for P,M in zip(PkP,PkM)]
            Pderivs.append(Ders)
                
        for zi,z in enumerate(self.zval):
            Fz=np.zeros((N,N))
            for i1,(n1,u1,v1) in enumerate(nal):
                if not u1:
                    continue
                for i2,(n2,u2,v2) in enumerate(nal):
                    if (i2<i1):
                        continue
                    if not u2:
                        continue
                    Fz[i1,i2]=(Pderivs[i1][zi]*PkEI[zi]*Pderivs[i2][zi]).sum()
                    Fz[i2,i1]=Fz[i1,i2]
            Fz=np.array(Fz)
            Fzm=la.inv(la.inv(Fz+np.diag([1e-30]*N))[:N-2,:N-2])
            F+=Fzm#Fz[:N-2,:N-2]
                    
        print(F)
        FishMat.__init__(self,[x[0] for x in nal[:-2]],F)
        
                
    
                
    def getCube(self,nal,p,o):
        for n,u,v in nal:
            if n==p:
                v*=(1+o)
            if n=='omegac':
                omegac=v
            elif n=='As':
                As=v
            elif n=='theta':
                theta=v
            elif n=='Neff':
                Neff=v
            elif n=='mnu':
                mnu=v
            elif n=='omegab':
                omegab=v
            elif n=='ns':
                ns=v
            elif n=='bd':
                bd=v
            elif n=='be':
                be=v
            
        print ("o=",o)

        zstr=",".join(map(str,self.zval+[self.zval[-1]+2]))
        pars = {
                'output': 'mPk',
                'P_k_max_h/Mpc': self.kval[-1]+3.0,
                '100*theta_s' : 100*theta,
                'tau_reio': 0.07,       
                'omega_cdm': omegac,     
                'A_s': 2.204e-9,          
                'N_ur': Neff-1,  
                'N_ncdm': 1.0,           
                'm_ncdm': mnu,       
                'omega_b': omegab,     
                'n_s': ns,          
                'z_pk' : zstr,
        }
        cosmo = Class()
        cosmo.set(pars)
        cosmo.compute()
        print (cosmo.sigma8(),cosmo.h(),cosmo.Neff(),cosmo.Omega_m())
        bg=cosmo.get_background()
        zs=bg['z']
        #Da=interp1d(zs,cosmo.h*bg['comov. dist.'])## in flat universe, comoving angular is... in Mpc/h
        #Hi=interp1d(zs,cosmo.h()*1./(bg['H [1/Mpc]'])) # in Mpc/h
        Da=interp1d(zs,bg['comov. dist.'])## cosmo.pk is actually all Mpc units
        Hi=interp1d(zs,1./(bg['H [1/Mpc]'])) # 

        if (p<0):
            self.Da_fid=Da
            self.Hi_fid=Hi
            ## now also calculate nmodes
            self.nmodes=[]
            
            for z in self.zval:
                zlow=z-self.dz/2
                zhigh=z+self.dz/2
                V=self.fsky*4*np.pi/3*(Da(zhigh)**3-Da(zlow)**3)
                Vk=2*np.pi*self.kperp*self.dk*self.dk
                cnm=V*Vk/(2*(2*np.pi)**3)
                self.nmodes.append(cnm)
        pkl=[]
        for z in self.zval:
            kperp_t=self.kperp/Da(z)*self.Da_fid(z) ## we are observing radians, so..
            kpar_t=self.kpar/Hi(z)*self.Hi_fid(z)
            kt=np.sqrt(kperp_t**2+kpar_t**2)
            #print ("z=",z)
            cpk=[cosmo.pk(k,z) for k in kt.flatten()]
            cpk=np.array(cpk).reshape(kt.shape)
            cpk*=(bd+be*self.mu**2)**2

            pkl.append(cpk)
        return pkl

    
        
