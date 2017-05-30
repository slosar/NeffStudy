#
# Wrapper for class
#
from __future__ import division, print_function
from classy import Class
import copy
from scipy.interpolate import interp1d
from ParameterVec import DefaultParamList, ParamList, Parameter
import sys
import numpy as np

class PkDiffer:

    def __init__ (self,pl, zvals, kvals, kperp, kpar):
        """ returns a list of Pks, each list containins 3D PS"""
        self.zstr=",".join(map(str,zvals+[zvals[-1]+2]))
        self.kvals=kvals
        self.kperp=kperp
        self.kpar=kpar
        self.zvals=zvals
        self.plist=copy.deepcopy(pl)
        self.cosmo = Class()
        self.ComputeCosmo(pl)
        bg=self.cosmo.get_background()
        zs=bg['z']
        Da=interp1d(zs,bg['comov. dist.'])## cosmo.pk is actually all Mpc units
        Hi=interp1d(zs,1./(bg['H [1/Mpc]'])) # 
        self.Da_fid=Da
        self.Hi_fid=Hi
        self.cube0=self.getCube(pl,'store_fid')

    def getDerivative(self, pa, frac):
        de=[]
        # if pa.name=='theta':
        #     ufrac=frac/50.
        # else:
        #     ufrac=frac
        ufrac=frac
        npl=copy.deepcopy(self.plist)
        for fa in [+1,-1]:
            nval=pa.value*(1.+1.0*fa*ufrac)
            npl.setValue(pa.name,nval)
            if "b_" in pa.name:
                mode="use_fid"
            else:
                mode="normal"
                self.ComputeCosmo(npl)
            de.append(self.getCube(npl,mode))
        toret=[(dp-dm)/(2*pa.value*fa*ufrac) for dp,dm in zip(de[1],de[0])]
        return toret
    
    
        
    def getCube(self,pl,mode='normal'):
        """mode defines caching of power spectra for biases
        mode can be 'store_fid', 'use_fid' or normal"""
        bg=self.cosmo.get_background()
        zs=bg['z']
        Da=interp1d(zs,bg['comov. dist.'])## cosmo.pk is actually all Mpc units
        Hi=interp1d(zs,1./(bg['H [1/Mpc]'])) # 
        if (mode=='store_fid'):
            self.Da_fid=Da
            self.Hi_fid=Hi
            self.cpk_cached, self.mu_cached=[], []
        pkl = []
        for i,z in enumerate(self.zvals):
            if (mode=='use_fid'):
                cpk=self.cpk_cached[i]
                mu=self.mu_cached[i]
            else:
                kperp_t=self.kperp/Da(z)*self.Da_fid(z) ## we are observing radians, so..
                kpar_t=self.kpar/Hi(z)*self.Hi_fid(z)
                kt=np.sqrt(kperp_t**2+kpar_t**2)
                mu=kpar_t/kt
                #mu=self.kpar/np.sqrt(self.kperp**2+self.kpar**2)
                cpk=[self.cosmo.pk(k,z) for k in kt.flatten()]
                cpk=np.array(cpk).reshape(kt.shape)
                if mode=='store_fid':
                    self.cpk_cached.append(cpk)
                    self.mu_cached.append(mu)
            f=self.growth_f(z)
            bpk=cpk*(pl.value('b_delta_'+str(i))+pl.value('b_eta_'+str(i))*f*mu**2)**2
            pkl.append(bpk)

        return pkl

    def growth_f(self,z):
        da=0.01
        a=1./(1.+z)
        gp,g,gm=[self.cosmo.scale_independent_growth_factor(1./ia-1.) for ia in [a+da,a,a-da]]
        f=a*(gp-gm)/(2*g*da)
        return f

    def ComputeCosmo(self,pl):
        #del self.cosmo
        #self.cosmo = Class()
        pars = {
                'output': 'mPk',
                'P_k_max_h/Mpc': self.kvals[-1]+3.0,
                '100*theta_s' : 100*pl.value('theta'),
                'tau_reio': pl.value('tau'),
                'omega_cdm': pl.value('omegac'),     
                'A_s': pl.value('As'),          
                'N_ur': pl.value('Neff')-1,  
                'N_ncdm': 1.0,           
                'm_ncdm': pl.value('mnu'),       
                'omega_b': pl.value('omegab'),     
                'n_s': pl.value('ns'),          
                'z_pk' : self.zstr,
        }
        self.cosmo.set(pars)
        #print ("Calling class compute...",end='')
        self.cosmo.compute()
        #print ("done")

        
 
