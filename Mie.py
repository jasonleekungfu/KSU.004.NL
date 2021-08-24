'''
========================================================================

    File Name:      Mie.py
    Author:         Jason Li
    Description:    Mie theory calculation module.
                    (Requires mathLib.py)
    Usage:          
                    [Use as module]:
                    
                    Mie
                        Plane wave Mie scattering methods.
                    Mie.__init__(n, I0=1, omega=1.722*eV, R=50*nm, \
                        Ri=0*nm, mu=1, maxOrder=5)
                        Constructor.
                    Mie.update()
                        Update private attributes, after any changes.
                    Mie.calMieCoefficients()
                        Calculate Mie coefficients.
                    Mie.EFieldS(co)
                        Calculate electric field at given spherical 
                        coodinates.
                    Mie.EFieldC(co)
                        Calculate electric field at given Cartesian 
                        coodinates.
                    
                    MiePulse
                        Pulse version scattering methods.
                    MiePulse.__init__(n, I0=1, omega=1.722*eV, R=50*nm,\
                        Ri=0*nm, mu=1, maxOrder=5, tau=2.47*fsec, \
                        nPW = 10)
                        Constructor.
                    MiePulse.calMieCoefficients()
                        Calculate Mie coefficients.
                    MiePulse.EFieldS(co)
                        Calculate electric field at given spherical 
                        coodinates.
                    MiePulse.EFieldC(co)
                        Calculate electric field at given Cartesian 
                        coodinates.
                    
                    
                    [Run separately]:
                    
                    main()
                        Main program.
                    plotField()
                        Plot electric field strength distribution.
                    plotCurve()
                        Plot enhancement and phase shift curves.
                    plotPulse()
                        Plot electric field strength distribution as a \
                        function of time for a pulse.
                    plotPulsePt():
                        Plot electric field strength at one point as a \
                        function of time.
                    plotPulseSurface():
                        Plot electric field strength on surface as a \
                        function of time.
                        
========================================================================
'''

import debug as db
from const import *
from mathLib import *
from scipy import interpolate
from scipy.linalg import solve
import matplotlib
import matplotlib.pyplot as plt
#from numba import jit


#=======================================================================
#                           Plane wave version
#=======================================================================

class Mie:
    """Mie calculation class. (Plane wave version)"""
    
    #-------------------------------------------------------------------
    #                           Public methods
    #-------------------------------------------------------------------
    
    def __init__(self, n, ni, I0=1, omega=1.722*eV, R=50*nm, Ri=0*nm, \
                    mu=1, maxOrder=5):
        """
        Constructor.
        Input:  n           --Shell complex permittivity.
                ni          --Core complex permittivity.
                I0          --Incident wave intensity. (Default: 1 a.u.)
                omega       --Incident center wavelength. 
                                (Default: 1.722 eV / 720 nm)
                R           --Nanosphere radius. (Default: 50 nm)
                Ri          --Nanoshell inner radius, 0 if solid. \
                                (Default: 0 nm)
                mu          --Nanosphere permeability. (Default: 1)
                maxOrder    --Maximum order to keep. (Default: 5)
        Output: oMie        --Mie calculation object.
        """
        
        self.n = n                  # Shell complex permittivity.
        self.ni = ni                # Core complex permittivity.
        self.I0 = I0                # Incident wave intensity.
        self.omega = omega          # Incident center wavelength.
        self.R = R                  # Nanosphere radius.
        self.Ri = Ri                # Nanoshell inner radius.
        self.mu = mu                # Nanosphere permeability.
        self.maxOrder = maxOrder    # Maximum order to keep.
        
        self.__ar = np.zeros(maxOrder, dtype=complex)
        self.__at = np.ones(maxOrder, dtype=complex)
        self.__br = np.zeros(maxOrder, dtype=complex)
        self.__bt = np.ones(maxOrder, dtype=complex)
        
        if (self.Ri>0): #Initiate shell parameters
            self.__c1 = np.zeros(maxOrder, dtype=complex)
            self.__c2 = np.zeros(maxOrder, dtype=complex)
            self.__d1 = np.zeros(maxOrder, dtype=complex)
            self.__d2 = np.zeros(maxOrder, dtype=complex)
        
        self.update()
        
    def update(self, nlr=0, nlri=0):
        """
        Update private attributes, after any changes.
        Input: nlr  --Nonlinear response term of shell, if applicable
               nlri --Nonlinear response term of core, if applicable
        """
        self.__E0 = np.sqrt(self.I0)
        self.__k = self.omega / c * (self.getn(self.n)+nlr)
        self.__ki = self.omega / c * (self.getn(self.ni)+nlri)
        self.__k0 = self.omega / c
        self.__wl = 2*Pi*c / self.omega
        
    def calMieCoefficients(self):
        """Calculate Mie coefficients."""
            
        # Whether it is a shell or nanoparticle
        if (self.Ri > 0):
            
            rho = self.__k0 * self.R
            nrho = self.__k * self.R
            Nk = self.__k / self.__k0
            
            rho1 = self.__ki * self.Ri
            nrho1 = self.__k * self.Ri
            Nk1 = self.__k / self.__ki
            
            for i in range(self.maxOrder):
                n=i+1
                A = np.array([ \
                    [0,-spherical_hn1(n,rho),0,0,\
                     spherical_jn(n,nrho), spherical_yn(n,nrho),0,0], \
                    [0,0,0,-Nk*self.__rzPrime(spherical_hn1, n, rho), \
                     0,0,self.__rzPrime(spherical_jn, n, nrho), \
                     self.__rzPrime(spherical_yn, n, nrho)], \
                    [0,-self.mu*self.__rzPrime(spherical_hn1, n, rho),\
                     0,0,self.__rzPrime(spherical_jn, n, nrho), \
                     self.__rzPrime(spherical_yn, n, nrho),0,0], \
                    [0,0,0,-self.mu*spherical_hn1(n,rho),0,0,\
                     Nk*spherical_jn(n,nrho),Nk*spherical_yn(n,nrho)], \
                    [spherical_jn(n,rho1),0,0,0,-spherical_jn(n,nrho1),\
                     -spherical_yn(n,nrho1),0,0], \
                    [0,0,Nk1*self.__rzPrime(spherical_jn, n, rho1),0,0,\
                     0,-self.__rzPrime(spherical_jn, n, nrho1),\
                     -self.__rzPrime(spherical_yn, n, nrho1)], \
                    [self.mu*self.__rzPrime(spherical_jn, n, rho1),0,\
                     0,0,-self.__rzPrime(spherical_jn, n, nrho1),\
                     -self.__rzPrime(spherical_yn, n, nrho1),0,0], \
                    [0,0,self.mu*spherical_jn(n,rho1),0,0,0,\
                     -Nk1*spherical_jn(n,nrho1),\
                     -Nk1*spherical_yn(n,nrho1)]\
                    ])
                B = np.array([ \
                    spherical_jn(n,rho), \
                    Nk*self.__rzPrime(spherical_jn, n, rho),\
                    self.mu*self.__rzPrime(spherical_jn, n, rho),\
                    self.mu*spherical_jn(n,rho),\
                    0,\
                    0,\
                    0,\
                    0\
                    ])
                (self.__at[i], \
                 self.__ar[i], \
                 self.__bt[i], \
                 self.__br[i], \
                 self.__c1[i], \
                 self.__c2[i], \
                 self.__d1[i], \
                 self.__d2[i]) = solve(A, B)
            
        else:
            
            rho = self.__k0 * self.R
            nrho = self.__k * self.R
            Nk = self.__k / self.__k0
            
            for i in range(self.maxOrder):
                n=i+1
                A = np.array([ \
                    [spherical_jn(n,nrho), -spherical_hn1(n,rho)], \
                    [self.__rzPrime(spherical_jn, n, nrho), \
                    -self.mu*self.__rzPrime(spherical_hn1, n, rho)] \
                    ])
                B = np.array([spherical_jn(n,rho), \
                    self.mu*self.__rzPrime(spherical_jn, n, rho)])
                (self.__at[i], self.__ar[i]) = solve(A, B)
                A = np.array([ \
                    [Nk*spherical_jn(n,nrho), \
                    -self.mu*spherical_hn1(n,rho)], \
                    [self.__rzPrime(spherical_jn, n, nrho), \
                    -Nk*self.__rzPrime(spherical_hn1, n, rho)] \
                    ])
                B = np.array([self.mu*spherical_jn(n,rho), \
                    Nk*self.__rzPrime(spherical_jn, n, rho)])
                (self.__bt[i], self.__br[i]) = solve(A, B)
    
    def EFieldS(self, co):
        """
        Calculate electric field at given spherical coodinates.
        Input:  co      --Coordinates in spherical.
        Output: EFieldS --Electric field in spherical.
        """
        
        # Initialization
        res = [0, 0, 0]
    
        # Eliminate singularity
        if (abs(co[1]) < small):
            co[1] = small
        elif (abs(co[1]-Pi) < small):
            co[1] = Pi - small
        
        # Calculate Legendre functions of all orders
        Legendres = self.__calLegendres(1, self.maxOrder, co)
        
        # Determin shell or nanoparticle
        if (self.Ri > 0):
            
            # Determine inside / shell /outside
            if (co[0]>self.R):
            
                # Calculate Bessel functions of all orders
                BesselsI = self.__calBessels(spherical_jn, 1, \
                                self.maxOrder, co)
                BesselsR = self.__calBessels(spherical_hn1, 1, \
                                self.maxOrder, co)
            
                # Calculate field
                for i in range(self.maxOrder):
                    n=i+1
                    self.__mfunc(n, co, BesselsI, Legendres, \
                        1j**n*(2*n+1)/n/(n+1)*self.__E0, res)
                    self.__nfunc(n, co, BesselsI, Legendres, \
                        1j**i*(2*n+1)/n/(n+1)*self.__E0, res)
                    self.__mfunc(n, co, BesselsR, Legendres, \
                        1j**n*(2*n+1)/n/(n+1)*self.__ar[i]*self.__E0, \
                        res)
                    self.__nfunc(n, co, BesselsR, Legendres, \
                        1j**i*(2*n+1)/n/(n+1)*self.__br[i]*self.__E0, \
                        res)
                
            elif (co[0]<self.Ri):
            
                # Calculate Bessel and Legendre functions of all orders
                BesselsT = self.__calBessels(spherical_jn, 1, \
                                self.maxOrder, co)
            
                # Calculate field
                for i in range(self.maxOrder):
                    n=i+1
                    self.__mfunc(n, co, BesselsT, Legendres, \
                        1j**n*(2*n+1)/n/(n+1)*self.__at[i]*self.__E0, \
                        res)
                    self.__nfunc(n, co, BesselsT, Legendres, \
                        1j**i*(2*n+1)/n/(n+1)*self.__bt[i]*self.__E0, \
                        res)
                
            else:
            
                # Calculate Bessel functions of all orders
                Bessels1 = self.__calBessels(spherical_jn, 1, \
                                self.maxOrder, co)
                Bessels2 = self.__calBessels(spherical_yn, 1, \
                                self.maxOrder, co)
            
                # Calculate field
                for i in range(self.maxOrder):
                    n=i+1
                    self.__mfunc(n, co, Bessels1, Legendres, \
                        1j**n*(2*n+1)/n/(n+1)*self.__c1[i]*self.__E0, \
                        res)
                    self.__nfunc(n, co, Bessels1, Legendres, \
                        1j**i*(2*n+1)/n/(n+1)*self.__d1[i]*self.__E0, \
                        res)
                    self.__mfunc(n, co, Bessels2, Legendres, \
                        1j**n*(2*n+1)/n/(n+1)*self.__c2[i]*self.__E0, \
                        res)
                    self.__nfunc(n, co, Bessels2, Legendres, \
                        1j**i*(2*n+1)/n/(n+1)*self.__d2[i]*self.__E0, \
                        res)
            
        else:
            
            # Determine inside / outside
            if (co[0]>self.R):
            
                # Calculate Bessel functions of all orders
                BesselsI = self.__calBessels(spherical_jn, 1, \
                                self.maxOrder, co)
                BesselsR = self.__calBessels(spherical_hn1, 1, \
                                self.maxOrder, co)
            
                # Calculate field
                for i in range(self.maxOrder):
                    n=i+1
                    self.__mfunc(n, co, BesselsI, Legendres, \
                        1j**n*(2*n+1)/n/(n+1)*self.__E0, res)
                    self.__nfunc(n, co, BesselsI, Legendres, \
                        1j**i*(2*n+1)/n/(n+1)*self.__E0, res)
                    self.__mfunc(n, co, BesselsR, Legendres, \
                        1j**n*(2*n+1)/n/(n+1)*self.__ar[i]*self.__E0, \
                        res)
                    self.__nfunc(n, co, BesselsR, Legendres, \
                        1j**i*(2*n+1)/n/(n+1)*self.__br[i]*self.__E0, \
                        res)
            
            else:
            
                # Calculate Bessel and Legendre functions of all orders
                BesselsT = self.__calBessels(spherical_jn, 1, \
                                self.maxOrder, co)
            
                # Calculate field
                for i in range(self.maxOrder):
                    n=i+1
                    self.__mfunc(n, co, BesselsT, Legendres, \
                        1j**n*(2*n+1)/n/(n+1)*self.__at[i]*self.__E0, \
                        res)
                    self.__nfunc(n, co, BesselsT, Legendres, \
                        1j**i*(2*n+1)/n/(n+1)*self.__bt[i]*self.__E0, \
                        res)
        
        # Return
        return res
    
    def EFieldC(self, co):
        """
        Calculate electric field at given Cartesian coodinates.
        Input:  co      --Coordinates in Cartesian.
        Output: EFieldS --Electric field in Cartesian.
        """
        
        coS = CtS(co)
        return(vStC(coS, self.EFieldS(coS)))
        
        """
        # Transform return value into Cartesian coordinates
        tranA = [   [ np.sin(coS[1]) * np.cos(coS[2]), \
                      np.cos(coS[1]) * np.cos(coS[2]), \
                     -np.sin(coS[2])], \
                    [ np.sin(coS[1]) * np.sin(coS[2]), \
                      np.cos(coS[1]) * np.sin(coS[2]), \
                      np.cos(coS[2])], \
                    [ np.cos(coS[1]), -np.sin(coS[1]), 0]  ]
        
        # Return
        return blas.cgemv(1, tranA, self.EFieldS(coS))
        """
        
        
    #-------------------------------------------------------------------
    #                           Private methods
    #-------------------------------------------------------------------
    
    def getn(self, n):
        """Get complex indeces of refraction."""
        x = np.array(list(n.keys()))
        x.sort()
        y1 = np.array(list(map(lambda xx: n[xx], x))).real
        y2 = np.array(list(map(lambda xx: n[xx], x))).imag
        tck1 = interpolate.splrep(x, y1, s=0)
        tck2 = interpolate.splrep(x, y2, s=0)
        xnew = np.array([self.omega/eV])
        ynew1 = interpolate.splev(xnew, tck1, der=0)
        ynew2 = interpolate.splev(xnew, tck2, der=0)
        return ynew1[0]+ynew2[0]*1j
        
    def __rzPrime(self, BesselFunc, n, rho):
        """Calculate rzPrime."""
        return BesselFunc(n, rho) + \
                rho * BesselFunc(n, rho, derivative=True)
           
    def __calBessels(self, BesselFunc, m, n, co):
        """Calculate the results of Bessel function and rzPrime."""
        nList = np.arange(n+1)
        
        # Determin shell or nanoparticle
        if (self.Ri > 0):
            
            if (co[0]>self.R or co[0]<self.Ri):
                Bessel0 = BesselFunc(nList, self.__k0*co[0])
                Bessels = [Bessel0 , Bessel0 + self.__k0*co[0] * \
                            BesselFunc(nList, self.__k0*co[0], \
                            derivative=True)]
            else:
                Bessel0 = BesselFunc(nList, self.__k*co[0])
                Bessels = [Bessel0 , Bessel0 + self.__k*co[0] * \
                            BesselFunc(nList, self.__k*co[0], \
                            derivative=True)]
            
        else:
            
            if (co[0]>self.R):
                Bessel0 = BesselFunc(nList, self.__k0*co[0])
                Bessels = [Bessel0 , Bessel0 + self.__k0*co[0] * \
                            BesselFunc(nList, self.__k0*co[0], \
                            derivative=True)]
            else:
                Bessel0 = BesselFunc(nList, self.__k*co[0])
                Bessels = [Bessel0, Bessel0 + self.__k*co[0] *\
                            BesselFunc(nList, self.__k*co[0], \
                            derivative=True)]
        return Bessels
           
    def __calLegendres(self, m, n, co):
        """Calculate the results of Legendre function."""
        return lpmnt(m, n, co[1])
    
    def __mfunc(self, n, co, Bessels, Legendres, coe, res):
        """Special function m"""
        
        # Theta part
        res[1] = res[1] + coe/np.sin(co[1]) * Bessels[0][n] \
                    * Legendres[0][1][n] * np.cos(co[2])
            
        # Phi part
        res[2] = res[2] - coe * Bessels[0][n] * Legendres[1][1][n] \
                    * np.sin(co[2])
    
    def __nfunc(self, n, co, Bessels, Legendres, coe, res):
        """Special function n"""
        
        # Initialization
        cosPhi = np.cos(co[2])
        
        # Ignore if at center
        if (co[0]==0):
            return res
                
        # Modify coefficient
        # Determin shell or nanoparticle
        if (self.Ri > 0):
            
            if (co[0]>self.R):
                coe = coe/co[0]/self.__k0
            elif (co[0]<self.Ri):
                coe = coe/co[0]/self.__ki
            else:
                coe = coe/co[0]/self.__k
            
        else:
            
            if (co[0]>self.R):
                coe = coe/co[0]/self.__k0
            else:
                coe = coe/co[0]/self.__k
        
        
        # Radial part
        res[0] = res[0] + coe * n*(n+1) * Bessels[0][n] \
                     * Legendres[0][1][n] * cosPhi
            
        # Theta part
        res[1] = res[1] + coe * Bessels[1][n] * Legendres[1][1][n] \
                    * cosPhi
            
        # Phi part
        res[2] = res[2] - coe * Bessels[1][n]/np.sin(co[1]) \
                         * Legendres[0][1][n] * np.sin(co[2])



#=======================================================================
#                             Pulse version
#=======================================================================

class MiePulse:
    """Mie calculation class. (Pulse version, remain to be changed)"""
    
    #-------------------------------------------------------------------
    #                           Public methods
    #-------------------------------------------------------------------
    
    def __init__(self, n, I0=1, omega=1.722*eV, R=50*nm, Ri=0*nm, \
                    mu=1, maxOrder=5, tau=2.47*fsec, nPW = 10):
        """
        Constructor.
        Input:  n           --Material complex permittivity.
                I0          --Incident wave intensity. (Default: 1 a.u.)
                omega       --Incident center wavelength. 
                                (Default: 1.722 eV / 720 nm)
                R           --Nanosphere radius. (Default: 50 nm)
                Ri          --Nanoshell inner radius, 0 if solid. \
                                (Default: 0 nm)
                mu          --Nanosphere permeability. (Default: 1)
                maxOrder    --Maximum order to keep. (Default: 5)
                tau         --FWHIM. (Default: 2.47 fs)
                nPW         --Plane wave components. (Default: 10)
        Output: oMiePulse   --Mie calculation object (Pulse version).
        """
        
        # Attributes
        self.n = n                  # Material complex permittivity.
        self.I0 = I0                # Incident wave intensity.
        self.omega = omega          # Incident center wavelength. 
        self.R = R                  # Nanosphere radius.
        self.Ri = Ri                # Nanoshell inner radius.
        self.mu = mu                # Nanosphere permeability.
        self.maxOrder = maxOrder    # Maximum order to keep.
        self.tau = tau              # FWHIM.
        self.nPW = nPW              # Plane wave components.
        
        self.__E0 = np.sqrt(I0)
        self.__PW = []
        
        # Calculate plane wave components
        tauOmega = 4*np.log(2)/tau
        dOmega   = 1.5*tauOmega/nPW
        for i in range(-nPW, nPW+1):
            amp = np.sqrt(2*np.log(2)/Pi)/tauOmega * \
                    np.exp(-2*np.log(2)*(i*dOmega/tauOmega)**2) * \
                    self.__E0*dOmega
            self.__PW.append(Mie(n, I0=amp, omega=omega+i*dOmega, R=R, \
                    Ri=Ri, mu=mu, maxOrder=maxOrder))
        
        #-------------------------------------------------
        # Center wavelength approximation (Remove later)
        self.__PWC = Mie(n, I0=I0, omega=omega, R=R, Ri=Ri, mu=mu, \
                    maxOrder=maxOrder)
        #-------------------------------------------------
        
    def calMieCoefficients(self):
        """Calculate Mie coefficients."""
            
        #-------------------------------------------------
        # Center wavelength approximation (Remove later)
        self.__PWC.calMieCoefficients()
        #-------------------------------------------------
        
        for i in range(0, 2*self.nPW+1):
            self.__PW[i].calMieCoefficients()
    
    def EFieldS(self, co, t):
        """
        Calculate electric field at given spherical coodinates.
        Input:  co      --Coordinates in spherical.
                t       --Time
        Output: EFieldS --Electric field in spherical.
        """
        
        # Initialization
        res = [0, 0, 0]
        
        #-------------------------------------------------
        # Center wavelength approximation (Remove later)
        
        tz = co[0]*np.cos(co[1])/c
                 
        # Gaussian shape
        res = blas.cscal( Gaussian(t-tz, self.tau)*\
                    np.exp(-1j*self.omega*t),  self.__PWC.EFieldS(co))
        
        """   
        # Cosine square shape
        if abs(t-tz)>self.tau:
            res = [0, 0, 0]
        else:
            res = blas.cscal( np.cos((t-tz)/self.tau*Pi/2)**2 * \
                    np.exp(-1j*self.omega*t),  self.__PWC.EFieldS(co))
        """
        
        return res
        #-------------------------------------------------
        
        # Eliminate singularity
        for i in range(0, 2*self.nPW+1):
            res = blas.caxpy(self.__PW[i].EFieldS(co), \
                    res, \
                    a = np.exp(-1j*self.__PW[i].omega*t))
        
        # Return
        return res
        
    
    def EFieldC(self, co, t):
        """
        Calculate electric field at given Cartesian coodinates.
        Input:  co      --Coordinates in Cartesian.
                t       --Time
        Output: EFieldS --Electric field in Cartesian.
        """
        
        coS = CtS(co)
        return(vStC(coS, self.EFieldS(coS, t)))
        
        """
        # Transform return value into Cartesian coordinates
        tranA = [   [ np.sin(coS[1]) * np.cos(coS[2]), \
                      np.cos(coS[1]) * np.cos(coS[2]), \
                     -np.sin(coS[2])], \
                    [ np.sin(coS[1]) * np.sin(coS[2]), \
                      np.cos(coS[1]) * np.sin(coS[2]), \
                      np.cos(coS[2])], \
                    [ np.cos(coS[1]), -np.sin(coS[1]), 0]  ]
        
        # Return
        return blas.cgemv(1, tranA, self.EFieldS(coS, t))
        """
    

#=======================================================================
#                   Plot field (if run separately)
#=======================================================================

def findMaxEnh(oMie, minR, maxR, minTheta, maxTheta):
    from scipy.linalg import norm
    nsteps = 20
    maxEnh = 0
    maxCoS = []
    for r in np.linspace(minR, maxR, nsteps):
        for theta in np.linspace(minTheta, maxTheta, nsteps):
            intensity = norm(oMie.EFieldS([r,theta,0]))**2 / oMie.I0
            if (intensity > maxEnh):
                maxEnh = intensity
                maxCoS = [r,theta,0]
    return(maxEnh, maxCoS)
    
def plotCurve():
    """Plot enhancement vs intensity at given wavelength."""
    import setup as st
    from scipy.linalg import norm
    
    # Initialization
    ksi = (-9.1+0.35j)*1e-19  #2e-19
    #aryIntensity = np.array([0.08, 0.18, 0.4, 0.8, 1.65, 8, ])# 1400 nm
    aryIntensity = np.linspace(0.08, 8, 100)                # 780 nm
    
    #-------------------------------------------------------------------
    #                           Core-shell
    #-------------------------------------------------------------------
    oMie = Mie(st.nAu, st.nSiO2, omega=st.omegaIR, R=st.rs, \
                Ri=st.ri)
    with open("Output/enh_I-shell.dat", "w") as f:
        
        for i in aryIntensity:
            
            # Initial trial intensity
            IIR = i*1e12*Wpcm2
            enhEff = 1
            
            # Loop several times to get converged intensity and 
            # nonlinear response
            for j in range(1):
                
                # Calculate nonlinear response term
                n2 = ksi*283/oMie.getn(oMie.n)/oMie.getn(oMie.n).real
                nlr = n2*1e4/Wpcm2 * IIR * enhEff
                #print(n2)
            
                # Update parameters and calculate Mie coefficients
                oMie.I0 = IIR
                oMie.update(nlr=nlr)
                oMie.calMieCoefficients()
            
                # Calculate field enhancement
                maxEnh, maxCoS = findMaxEnh(oMie, oMie.R+small, \
                        oMie.R+nm, 80*deg, 100*deg)
                
                # Get next improved intensity
                #IIR = i*1e12*Wpcm2 * maxEnh
                #enhEff = (enhEff + (\
                #          norm(oMie.EFieldS([oMie.R-small,0,0]))**2/oMie.I0+\
                #          2*norm(oMie.EFieldS([oMie.R-small,Pi/2,0]))**2/oMie.I0)\
                #          /3)/2
                #print(norm(oMie.EFieldS([oMie.R-small,0,0]))**2/oMie.I0, \
                #      norm(oMie.EFieldS([oMie.R-small,Pi/2,0]))**2/oMie.I0, \
                #      i)
                enhEff = (enhEff + norm(oMie.EFieldS([oMie.R-small,0,0]))**2/oMie.I0)/2
                print(enhEff, i)
            
            # Write converged results
            line = "{:f}\t{:f}\t{:f}\t{:f}\n".format(i, maxEnh, \
                    maxCoS[0]/nm, maxCoS[1]/deg)
            f.write(line)
    
    #-------------------------------------------------------------------
    #                            Solid
    #-------------------------------------------------------------------
    oMie = Mie(st.nAu, st.nAu, omega=st.omegaIR, R=st.rs, \
                Ri=0, I0=st.IIR)
    with open("Output/enh_I-solid.dat", "w") as f:
        
        for i in aryIntensity:
            
            # Initial trial intensity
            IIR = i*1e12*Wpcm2
            enhEff = 1
            
            # Loop several times to get converged intensity and 
            # nonlinear response
            for j in range(1):
                
                # Calculate nonlinear response term
                n2 = ksi*283/oMie.getn(oMie.n)/oMie.getn(oMie.n).real
                nlr = n2*1e4/Wpcm2 * IIR * enhEff
            
                # Update parameters and calculate Mie coefficients
                oMie.I0 = IIR
                oMie.update(nlr=nlr)
                oMie.calMieCoefficients()
            
                # Calculate field enhancement
                maxEnh, maxCoS = findMaxEnh(oMie, oMie.R+small, \
                        oMie.R+nm, 80*deg, 100*deg)
                
                # Get next improved intensity
                #IIR = i*1e12*Wpcm2 * maxEnh
                enhEff = (enhEff + norm(oMie.EFieldS([oMie.R-small,0,0]))**2/oMie.I0)/2
                print(enhEff, i)
            
            # Write converged results
            line = "{:f}\t{:f}\t{:f}\t{:f}\n".format(i, maxEnh, \
                    maxCoS[0]/nm, maxCoS[1]/deg)
            f.write(line)

def plotCurve2():
    """Plot enhancement vs wavelength at given intensity."""
    import setup as st
    
    # Initialization
    ksi = (-9.1+0.35j)*1e-19  #2e-19
    aryIntensity = np.array([0.08, 0.3, 1.2, 2.7, 8])
    iI = 0
    intensity = aryIntensity[iI] * 1e12*Wpcm2
    
    #-------------------------------------------------------------------
    #                           Core-shell
    #-------------------------------------------------------------------
    oMie = Mie(st.nAu, st.nSiO2, omega=st.omegaIR, R=st.rs, \
                Ri=st.ri, maxOrder=15)
    with open("Output/enh_wl-shell{:d}_2.dat".format(iI), "w") as f:
        
        for i in range(400, 900, 10):
            
            # Initial trial intensity
            IIR = intensity
            enhEff = 1
            oMie.omega = 2*Pi*c/(i*nm)
            
            # Loop several times to get converged intensity and 
            # nonlinear response
            for j in range(1):
                
                # Calculate nonlinear response term
                n2 = ksi*283/oMie.getn(oMie.n)/oMie.getn(oMie.n).real
                nlr = n2*1e4/Wpcm2 * IIR * enhEff
            
                # Update parameters and calculate Mie coefficients
                oMie.I0 = IIR
                oMie.update(nlr=nlr)
                oMie.calMieCoefficients()
            
                # Calculate field enhancement
                maxEnh, maxCoS = findMaxEnh(oMie, oMie.R+small, \
                        oMie.R+nm, 80*deg, 100*deg)
                
                # Get next improved intensity
                IIR = intensity * maxEnh
            
            # Write converged results
            line = "{:f}\t{:f}\t{:f}\t{:f}\n".format(i, maxEnh, \
                    maxCoS[0]/nm, maxCoS[1]/deg)
            f.write(line)
    
    #-------------------------------------------------------------------
    #                            Solid
    #-------------------------------------------------------------------
    oMie = Mie(st.nAu, st.nAu, omega=st.omegaIR, R=st.rs, \
                Ri=0, maxOrder=5)
    with open("Output/enh_wl-solid{:d}_2.dat".format(iI), "w") as f:
        
        for i in range(400,900, 10):
            
            # Initial trial intensity
            IIR = intensity
            enhEff = 1
            oMie.omega = 2*Pi*c/(i*nm)
            
            # Loop several times to get converged intensity and 
            # nonlinear response
            for j in range(1):
                
                # Calculate nonlinear response term
                n2 = ksi*283/oMie.getn(oMie.n)/oMie.getn(oMie.n).real
                nlr = n2*1e4/Wpcm2 * IIR * enhEff
            
                # Update parameters and calculate Mie coefficients
                oMie.I0 = IIR
                oMie.update(nlr=nlr)
                oMie.calMieCoefficients()
            
                # Calculate field enhancement
                maxEnh, maxCoS = findMaxEnh(oMie, oMie.R+small, \
                        oMie.R+nm, 80*deg, 100*deg)
                
                # Get next improved intensity
                IIR = intensity * maxEnh
            
            # Write converged results
            line = "{:f}\t{:f}\t{:f}\t{:f}\n".format(i, maxEnh, \
                    maxCoS[0]/nm, maxCoS[1]/deg)
            f.write(line)
    
def plotCurve2e():
    """Plot enhancement vs wavelength at given intensity, with error"""
    import setup as st
    
    # Initialization
    ksi = 9.1e-19  #2e-19
    ndr = 20
    drs = 3.5*nm
    dri = 2*nm
    aryIntensity = np.array([0.08, 0.3, 1.2, 2.7, 8])
    iI = 4
    intensity = aryIntensity[iI] * 1e12*Wpcm2
    
    #-------------------------------------------------------------------
    #                           Core-shell
    #-------------------------------------------------------------------
    oMie = Mie(st.nAu, st.nSiO2, omega=st.omegaIR, R=st.rs, \
                Ri=st.ri, maxOrder=10)
    with open("Output/enh_wl-shell{:d}.dat".format(iI), "w") as f:
        
        for i in range(400,900, 10):
            
            oMie.omega = 2*Pi*c/(i*nm)
            enh = []
            
            # Loop several times to find the max and min enhancement
            for j in range(ndr+1):
                
                # Calculate nonlinear response term
                n2 = abs(ksi*283/oMie.getn(oMie.n)/\
                        oMie.getn(oMie.n).real)*1j
                nlr = n2 * aryIntensity[iI]*1e16 #* 3
            
                # Update parameters and calculate Mie coefficients
                oMie.I0 = intensity
                oMie.R = st.rs - drs + j * 2*drs/ndr
                oMie.Ri = st.ri + dri - j * 2*dri/ndr
                oMie.update(nlr=nlr)
                oMie.calMieCoefficients()
            
                # Calculate field enhancement
                maxEnh, maxCoS = findMaxEnh(oMie, oMie.R+small, \
                        oMie.R+nm, 80*deg, 100*deg)
                enh.append(maxEnh)
            
            # Write converged results
            line = "{:f}\t{:f}\t{:f}\n".format(i, max(enh), min(enh))
            f.write(line)

def plotCurve3():
    """Plot field enhancement at inner shell vs intensity at given \
    wavelength."""
    import setup as st
    
    # Initialization
    ksi = 7.6e-19
    aryIntensity = np.linspace(0.08, 8, 100)                # 780 nm
    
    #-------------------------------------------------------------------
    #                           Core-shell
    #-------------------------------------------------------------------
    oMie = Mie(st.nAu, st.nSiO2, omega=st.omegaIR, R=st.rs, \
                Ri=st.ri, maxOrder=15)
    with open("Output/trans_I-shell.dat", "w") as f:
        
        for i in aryIntensity:
            
            # Initial trial intensity
            IIR = i*1e12*Wpcm2
            
            # Loop several times to get converged intensity and 
            # nonlinear response
            for j in range(1):
                
                # Calculate nonlinear response term
                n2 = abs(ksi*283/oMie.getn(oMie.n)/\
                        oMie.getn(oMie.n).real)*1j
                nlr = n2 * i*1e16 #* 3
            
                # Update parameters and calculate Mie coefficients
                oMie.I0 = IIR
                oMie.update(nlr=nlr)
                oMie.calMieCoefficients()
            
                # Calculate field enhancement
                maxEnh, maxCoS = findMaxEnh(oMie, oMie.Ri+small, \
                        oMie.Ri+nm, small, 20*deg)
                        #oMie.Ri+small, 80*deg, 100*deg)
                
                # Get next improved intensity
                IIR = i*1e12*Wpcm2 * maxEnh
            
            # Write converged results
            line = "{:f}\t{:f}\t{:f}\t{:f}\n".format(i, maxEnh*i, \
                    maxCoS[0]/nm, maxCoS[1]/deg)
            f.write(line)
    
    #-------------------------------------------------------------------
    #                            Solid
    #-------------------------------------------------------------------
    oMie = Mie(st.nAu, st.nAu, omega=st.omegaIR, R=st.rs, \
                Ri=st.ri, I0=st.IIR, maxOrder=5)
    with open("Output/trans_I-solid.dat", "w") as f:
        
        for i in aryIntensity:
            
            # Initial trial intensity
            IIR = i*1e12*Wpcm2
            
            # Loop several times to get converged intensity and 
            # nonlinear response
            for j in range(1):
                
                # Calculate nonlinear response term
                n2 = abs(ksi*283/oMie.getn(oMie.n)/\
                        oMie.getn(oMie.n).real)*1j
                nlr = n2 * i*1e16 #* 3
            
                # Update parameters and calculate Mie coefficients
                oMie.I0 = IIR
                oMie.update(nlr=nlr)
                oMie.calMieCoefficients()
            
                # Calculate field enhancement
                maxEnh, maxCoS = findMaxEnh(oMie, oMie.Ri+small, \
                        oMie.Ri+nm, small, 20*deg)
                        #oMie.Ri+small, 80*deg, 100*deg)
                
                # Get next improved intensity
                IIR = i*1e12*Wpcm2 * maxEnh
            
            # Write converged results
            line = "{:f}\t{:f}\t{:f}\t{:f}\n".format(i, maxEnh*i, \
                    maxCoS[0]/nm, maxCoS[1]/deg)
            f.write(line)

def plotCurve4():
    """Plot enhancement vs diameter for linear response."""
    import setup as st
    
    ksi = 7.6e-19  #2e-19
    
    #-------------------------------------------------------------------
    #                            Solid
    #-------------------------------------------------------------------
    oMie = Mie(st.nAu, st.nAu, omega=st.omegaIR, R=st.rs, \
                Ri=st.ri, maxOrder=5)
    with open("Output/enh_D-solid.dat", "w") as f:
        
        for i in range(50,450, 5):
                
            # Calculate nonlinear response term
            IIR = 2e12*Wpcm2
            n2 = abs(ksi*283/oMie.getn(oMie.n)/\
                        oMie.getn(oMie.n).real)*1j
            nlr = n2 * IIR*1e16 #* 3
            
            # Update parameters and calculate Mie coefficients
            oMie.I0 = IIR
            oMie.R = i*nm/2
            oMie.update(nlr=nlr)
            oMie.calMieCoefficients()
            
            # Calculate field enhancement
            maxEnh, maxCoS = findMaxEnh(oMie, oMie.R+small, \
                        oMie.R+nm, 40*deg, 140*deg)
            
            # Write converged results
            line = "{:f}\t{:f}\t{:f}\t{:f}\n".format(i, maxEnh, \
                    maxCoS[0]/nm, maxCoS[1]/deg)
            f.write(line)
    """Plot enhancement vs wavelength at given intensity, with error"""
    import setup as st
    
    # Initialization
    ksi = 7.6e-19  #2e-19
    ndr = 20
    drs = 3.5*nm
    dri = 2*nm
    aryIntensity = np.array([0.08, 0.3, 1.2, 2.7, 8])
    iI = 1
    intensity = aryIntensity[iI] * 1e12*Wpcm2
    
    #-------------------------------------------------------------------
    #                           Core-shell
    #-------------------------------------------------------------------
    oMie = Mie(st.nAu, st.nSiO2, omega=st.omegaIR, R=st.rs, \
                Ri=st.ri, maxOrder=10)
    with open("Output/enh_wl-shell{:d}.dat".format(iI), "w") as f:
        
        for i in range(400,900, 10):
            
            oMie.omega = 2*Pi*c/(i*nm)
            enh = []
            
            # Loop several times to find the max and min enhancement
            for j in range(ndr+1):
                
                # Calculate nonlinear response term
                n2 = abs(ksi*283/oMie.getn(oMie.n)/\
                        oMie.getn(oMie.n).real)*1j
                nlr = n2 * aryIntensity[iI]*1e16 #* 3
            
                # Update parameters and calculate Mie coefficients
                oMie.I0 = intensity
                oMie.R = st.rs - drs + j * 2*drs/ndr
                oMie.Ri = st.ri + dri - j * 2*dri/ndr
                oMie.update(nlr=nlr)
                oMie.calMieCoefficients()
            
                # Calculate field enhancement
                maxEnh, maxCoS = findMaxEnh(oMie, oMie.R+small, \
                        oMie.R+nm, 80*deg, 100*deg)
                enh.append(maxEnh)
            
            # Write converged results
            line = "{:f}\t{:f}\t{:f}\n".format(i, max(enh), min(enh))
            f.write(line)

def plotDepth():
    """Plot skin depth vs intensity at given wavelength."""
    import setup as st
    
    # Initialization
    ksi = 9.1e-19
    aryIntensity = np.linspace(0.08, 8, 100)                # 780 nm
    
    #-------------------------------------------------------------------
    #                           Core-shell
    #-------------------------------------------------------------------
    oMie = Mie(st.nAu, st.nSiO2, omega=st.omegaIR, R=st.rs, \
                Ri=st.ri, maxOrder=15)
    print(oMie.getn(oMie.ni))
    with open("Output/depth.dat", "w") as f:
        
        for i in aryIntensity:
            
            # Initial trial intensity
            IIR = i*1e12*Wpcm2
            
            # Calculate nonlinear response term
            n2 = abs(ksi*283/oMie.getn(oMie.n)/\
                    oMie.getn(oMie.n).real)*1j
            nlr = n2 * i*1e16 #* 3
            depth = c/2/oMie.omega/(oMie.getn(oMie.n)+nlr).imag/nm
            
            # Write results
            line = "{:f}\t{:f}\n".format(i, depth)
            f.write(line)
    
def plotField():
    """Plot field for spherical and infinite cylindar."""
    import setup as st
    from scipy.linalg import norm
    from scipy.special import jv, jvp, hankel1, h1vp
    
    # Initialization
    ksi = 7.6e-19
    aryIntensity = np.array([0.08, 0.3, 8])
    iIntensity = 1
    nr = 40
    ntheta = 60
    r = np.linspace(0, 3*st.rs, nr)
    theta = np.linspace(-Pi/2, Pi/2, ntheta)
    enh = np.zeros([nr, ntheta*2])
    """
    nn = 50
    lim = 2*st.rs
    x = np.linspace(-lim, lim, nn)
    z = np.linspace(-lim, lim, nn)
    enh = np.zeros([nn, nn*2])
    """
    
    #-------------------------------------------------------------------
    #                            Sphere
    #-------------------------------------------------------------------
    
    with open("Output/EField{:d}.dat".format(iIntensity), "w") as f:
    
        #---------------------------------------------------------------
        #                         Core-Shell
        #---------------------------------------------------------------
        
        # Initialization
        oMie = Mie(st.nAu, st.nSiO2, omega=st.omegaIR, R=st.rs, \
                Ri=st.ri, maxOrder=15)
                
        # Calculate non-linear response and update Mie coefficients
        n2 = abs(ksi*283/oMie.getn(oMie.n)/oMie.getn(oMie.n).real)*1j
        nlr = n2 * aryIntensity[iIntensity]*1e16 #* 3
        oMie.update(nlr=nlr)
        oMie.calMieCoefficients()
        
        # Plot field enhancement
        for j in range(int(ntheta)):
            for i in range(nr):
                enh[i,j] = norm(oMie.EFieldS([r[i], theta[j], 0]))**2
                line = "{:f}\t{:f}\t{:f}\n".format(r[i]/nm, theta[j], \
                        enh[i,j])
                f.write(line)
            f.write("\n")
        """
        # Plot field enhancement
        for j in range(int(nn/2)):
            for i in range(nn):
                enh[i,j] = norm(oMie.EFieldC([x[i], 0, z[j]]))**2
                line = "{:f}\t{:f}\t{:f}\n".format(x[i]/nm, z[j]/nm, \
                        enh[i,j])
                f.write(line)
            f.write("\n")
        """
        #---------------------------------------------------------------
        #                           Solid
        #---------------------------------------------------------------
        
        # Initialization
        oMie = Mie(st.nAu, st.nAu, omega=st.omegaIR, R=st.rs)
                
        # Calculate non-linear response and update Mie coefficients
        n2 = abs(ksi*283/oMie.getn(oMie.n)/oMie.getn(oMie.n).real)*1j
        nlr = n2 * aryIntensity[iIntensity]*1e16 #* 3
        oMie.update(nlr=nlr)
        oMie.calMieCoefficients()
        
        # Plot field enhancement
        for j in range(int(ntheta)):
            for i in range(nr):
                enh[i,j+ntheta] = norm(oMie.EFieldS([r[i], theta[j], \
                                    0]))**2
                line = "{:f}\t{:f}\t{:f}\n".format(r[i]/nm, \
                        Pi+theta[j], enh[i,j+ntheta])
                f.write(line)
            f.write("\n")
        """
        # Plot field enhancement
        for j in range(int(nn/2),nn):
            for i in range(nn):
                enh[i,j] = norm(oMie.EFieldC([x[i], 0, z[j]]))**2
                line = "{:f}\t{:f}\t{:f}\n".format(x[i]/nm, z[j]/nm, \
                        enh[i,j])
                f.write(line)
            f.write("\n")
        """
            
    plt.pcolormesh(enh, cmap='nipy_spectral')
    plt.show()
    
    """
    #-------------------------------------------------------------------
    #                            Cylindar
    #-------------------------------------------------------------------
    
    oMie = Mie(st.nAu, st.nAu, omega=st.omegaIR, R=R, maxOrder=maxOrder)
    k0 = st.omegaIR / c
    k1 = st.omegaIR / c * (oMie.getn(st.nAu))
    ai = np.zeros(2*maxOrder+1, dtype=complex)
    ae = np.zeros(2*maxOrder+1, dtype=complex)
    for n in range(-maxOrder, maxOrder+1):
        A = np.array([  [jv(n,k1*R),      -hankel1(n,k0*R)],
                        [k1*jvp(n,k1*R),  -k0*h1vp(n,k0*R)]])
        B = np.array([jv(n,k0*R), k0*jvp(n,k0*R)])
        (ai[n], ae[n]) = solve(A, B)
        
    with open("Output/EField,dat", "w") as f:
        
        for i in range(-lim, lim, 10):
            for j in range(-lim, lim, 10):
                res = 0
                for n in range(-maxOrder, maxOrder+1):
                    r = norm(np.array([i*nm,j*nm]))
                    if (r>R):
                        res += jv(n, k0*r) + ae[n]*hankel1(n, k0*r) 
                    else:
                        res += ai[n]*jv(n, k0*r)
                line = "{:f}\t{:f}\t{:f}\n".format(i, j, \
                        abs(res)**2)
                f.write(line)
            f.write("\n")
    """

def exportEpsilon():
    """Export epsilon."""
    import setup as st
    
    # Initialization
    ksi = 2e-19
    aryIntensity = np.array([0, 0.08, 0.6, 1.2, 2.7, 8])      # 780 nm
    omegaAry = np.array(list(st.nAu.keys()))
    omegaAry.sort()
                
    # Loop
    for i in aryIntensity:
        
        # Create output file
        with open("Output/epsilon/epsilon-I{:f}.txt".format(i), "w") \
            as fepsilon, open("Output/epsilon/n-I{:f}.txt".format(i), \
            "w") as fn:
        
            # Loop through frequencies
            for omega in omegaAry:
                
                # Calculate nonlinear response term
                n2 = abs(ksi*283/(st.nAu[omega])/st.nAu[omega].real)*1j
                nlr = n2 * i*1e16 * 3
                
                # Calculate n
                n = st.nAu[omega] + nlr
                
                # Calculate epsilon
                epsilon = n**2
                
                # Export
                fepsilon.write("{:f}\t{:f}\t{:f}\n".format(\
                    omega*241798.93, epsilon.real, epsilon.imag))
                fn.write("{:f}\t{:f}\t{:f}\n".format(\
                    omega*241798.93, n.real, n.imag))
        
  
    from scipy.linalg import norm
    maxEnh = 0
    maxCoS = []
    for r in np.linspace(minR, maxR, 10):
        for theta in np.linspace(minTheta, maxTheta, 10):
            intensity = norm(oMie.EFieldS([r,theta,0]))**2 / oMie.I0
            if (intensity > maxEnh):
                maxEnh = intensity
                maxCoS = [r,theta,0]
    return(maxEnh, maxCoS)

def plotEnhVsD_SiO2():
    """Plot enhancement vs diameter for SiO2 linear response."""
    import setup as st
    from scipy.linalg import norm
    
    #-------------------------------------------------------------------
    #                            Solid
    #-------------------------------------------------------------------
    oMie = Mie(st.nSiO2, st.nSiO2, omega=st.omegaIR, R=st.rs)
    with open("Output/enh_D-SiO2.dat", "w") as f:
        
        for i in range(5,501, 5):
            
            # Update parameters and calculate Mie coefficients
            oMie.R = i*nm/2
            oMie.calMieCoefficients()
            
            # Calculate field enhancement
            #maxEnh, maxCoS = findMaxEnh(oMie, oMie.R+small, \
            #            oMie.R+nm, 30*deg, 140*deg)
            
            # Write converged results
            #line = "{:f}\t{:f}\n".format(i, \
            #            norm(oMie.EFieldS([oMie.R+small,Pi/2,0]))**2 / oMie.I0)
            line = "{:f}\t{:f}\t{:f}\t{:f}\n".format(i, maxEnh, \
                    maxCoS[0]/nm, maxCoS[1]/deg)
            f.write(line)
              
def main():
    """Plot field, if run separately."""
    db.timeStamp()
    plotCurve2()
    db.timeStamp()
	
if __name__ == '__main__':
    main()
