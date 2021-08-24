'''
========================================================================

    File Name:      mathLib.py
    Author:         Jason Li
    Description:    Customized mathematical library.
                    (This module by default includes numpy).
    Usage:          spherical_hn1(n, z, derivative=False)
                        Spherical Hankel function of 1st kind.
                    lpmnt(m, n, theta):
                        Associate Legendre polynomial with theta 
                        arguement.
                    Gaussian(t, tau, omega=0, phi=0, beta=0):
                        Evaluate Gaussian envelope.
                    StC(co):
                        Spherical -> Cartesian transformation.
                    CtS(co):
                        Cartesian -> Spherical transformation.
                    vStC(co, v, isComplex=True)
                        Vector Spherical -> Cartesian transformation.
                    vCtS(co, v, isComplex=True)
                        Vector Cartesian -> Spherical transformation.
                    aCoulomb(r, t=0, charge=-1):
                        Coulomb acceleration calculation.
                    
                    RK4
                        4th order Runge-Kutta method.
                    RK4.__init__(self, a=aCoulomb, nt=1000, dt=asec, \
                        isComplex=False):
                        Constructor.
                    RK4.getPath(self, i=None):
                        Obtain path by time index.
                    RK4.run(self, vecPhase, t0=0):
                        Propagate from initial position and time.
                    
========================================================================
'''

from const import *
from scipy.special import spherical_jn, spherical_yn, lpmn, sph_harm
import numpy as np
import scipy.linalg.blas as blas
    


#=======================================================================
#                           Misc Methods
#=======================================================================

def spherical_hn1(n, z, derivative=False):
    """
    Spherical Hankel function of 1st kind.
    Input:  n               --Degree.
            z               --Arguement (complex compatible).
            derivative      --Return derivative. (Defaul: False)
    Output: spherical_hn1   --Spherical Hankel function of the 1st kind.
    """
    return spherical_jn(n, z, derivative) + \
            1j * spherical_yn(n, z, derivative)

def lpmnt(m, n, theta):
    """
    Associate Legendre polynomial with theta arguement.
    Input:  m       --Rank.
            n       --Degree.
            theta   --Theta arguement.
    Output: lpmnt   --Associate Legendre polynomial.
    """
    
    # Calculate derivative of theta
    res = list(lpmn(m, n, np.cos(theta)))
    res[1] = -res[1]*np.sin(theta)
    
    return res

def Gaussian(t, tau, omega=0, phi=0, beta=0):
    """
    Evaluate Gaussian envelope.
    Input:  t           --Time.
            tau         --FWHIM.
            omega       --Center frequency. (Default: 0)
            phi         --Additional phase. (Default: 0)
            beta        --Chirp. (Default: 0)
    Output: Gaussian    --Gaussian envelope.
    """
    return np.exp(-2*np.log(2)*(t/tau)**2)*np.cos(omega*t+beta*t**2+phi)
    
def StC(co):
    """
    Spherical -> Cartesian transformation.
    Input:  co  --Spherical coordinates.
    Output:     --Cartesion coordinates.
    """
    (r, theta, phi) = co
    return [r*np.sin(theta)*np.cos(phi), \
            r*np.sin(theta)*np.sin(phi), \
            r*np.cos(theta)]

def CtS(co):
    """
    Cartesian -> Spherical transformation.
    Input:  co  --Cartesion coordinates.
    Output:      --Spherical coordinates.
    """
    (x, y, z) = co
    
    # r
    r = np.sqrt(x**2 + y**2 + z**2)
    
    # phi
    if (abs(x)<=small):
        if (y > small):
            phi = Pi/2
        elif (y < -small):
            phi = 3*Pi/2
        else:
            phi = 0
    else:
        phi = np.arctan2(y,x)
    
    # theta
    if (abs(x)<=small and abs(y)<=small and abs(z)<=small):
        theta = 0
    else:
        theta = np.arccos(z/r)
        
    return [r, theta, phi]
    
def vStC(co, v, isComplex=True):
    """
    Vector Spherical -> Cartesian transformation.
    Input:  co          --Spherical coordinates.
            v           --Vector in spherical.
            isComplex   --Whether the vector is complex. (Default: True)
    Output:             --Vector in Cartesian.
    """
    
    # Transformation matrix
    tranA = [   [ np.sin(co[1]) * np.cos(co[2]), \
                  np.cos(co[1]) * np.cos(co[2]), \
                 -np.sin(co[2])], \
                [ np.sin(co[1]) * np.sin(co[2]), \
                  np.cos(co[1]) * np.sin(co[2]), \
                  np.cos(co[2])], \
                [ np.cos(co[1]), -np.sin(co[1]), 0]  ]
    
    # Return
    if isComplex:
        return blas.cgemv(1, tranA, v)
    else:
        return blas.sgemv(1, tranA, v)
    
def vCtS(co, v, isComplex=True):
    """
    Vector Cartesian -> Spherical transformation.
    Input:  co          --Spherical coordinates.
            v           --Vector in Cartesian.
            isComplex   --Whether the vector is complex. (Default: True)
    Output:             --Vector in spherical.
    """
    
    # Transformation matrix
    tranA = [   [ np.sin(co[1]) * np.cos(co[2]), \
                  np.sin(co[1]) * np.sin(co[2]), \
                  np.cos(co[1])], \
                [ np.cos(co[1]) * np.cos(co[2]), \
                  np.cos(co[1]) * np.sin(co[2]), \
                 -np.sin(co[1])], \
                [-np.sin(co[2]),  np.cos(co[2]), 0]  ]
    
    # Return
    if isComplex:
        return blas.cgemv(1, tranA, v)
    else:
        return blas.dgemv(1, tranA, v)

def aCoulomb(r, t=0, charge=-1):
    """
    Coulomb acceleration calculation.
    Input:  r, t        --Spatiotemporal coordinates.
            charge      --Coulomb charge. (Default: -1)
    Output: aCoulomb    --Acceleration vector.
    """
    
    return(r*charge/blas.dznrm2(r)**3)


#=======================================================================
#                           RK4 Methods
#=======================================================================

class RK4:
    """4th order Runge-Kutta method."""
    
    #-------------------------------------------------------------------
    #                           Public methods
    #-------------------------------------------------------------------
    
    def __init__(self, a=aCoulomb, nt=1000, dt=asec, \
        isComplex=False):
        """
        Constructor.
        Input:  a           --Acceleration function a(vecR, t), \
                              returns a vector. (Default: aCoulomb)
                nt          --Time propagation steps. (Default: 1000)
                dt          --Time step length. (Default: 0.1 fs)
                isComplex   --Whether is complex path. (Default: False)
        Output: oRK4        --Runge-Kutta calculation object.
        """
        
        self.a = a
        self.nt = nt
        self.dt = dt
        self.isComplex = isComplex
        
        # Initialize data arrays.
        if (isComplex):
            self.__path = np.zeros([nt, 6], dtype=complex)
            self.__tmp = np.zeros(6, dtype=complex)
            self.__k1 = np.zeros(6, dtype=complex)
            self.__k2 = np.zeros(6, dtype=complex)
            self.__k3 = np.zeros(6, dtype=complex)
            self.__k4 = np.zeros(6, dtype=complex)
        else:
            self.__path = np.zeros([nt, 6], dtype=float)
            self.__tmp = np.zeros(6, dtype=float)
            self.__k1 = np.zeros(6, dtype=float)
            self.__k2 = np.zeros(6, dtype=float)
            self.__k3 = np.zeros(6, dtype=float)
            self.__k4 = np.zeros(6, dtype=float)
            
    def getPath(self, i=None):
        """
        Obtain path by time index.
        Input:  i           --Time index. Return full path if omitted.
        Output: path        --Path (or path array) in phase space.
        """
        
        if (i):
            return(self.__path[i])
        else:
            return(self.__path)
    
    def run(self, vecPhase, t0=0):
        """
        Propagate from initial position and time.
        Input:  vecPhase    --Initial vector in phase space.
                t0          --Initial time.
        """
        
        # Initialize
        self.__path[0] = vecPhase
        
        # Propagate
        if (self.isComplex):
            for i in range(self.nt-1):
                self.__runStepC(i, t0)
        else:
            for i in range(self.nt-1):
                self.__runStepR(i, t0)
        
    
    #-------------------------------------------------------------------
    #                           Private methods
    #-------------------------------------------------------------------
    
    def __runStepR(self, i=0, t0=0):
        """
        Propagate once from step i (Real version).
        Input:  i           --Time index. (Default: 0)
                t0          --Initial time. (Default: 0)
        """
        
        # Initialize
        t = t0 + i * self.dt
        ho2 = self.dt/2
        
        # 4th order Runge-Kutta
        np.copyto(self.__tmp, self.__path[i])
        self.__f(self.__tmp, t, self.__k1)
        np.copyto(self.__tmp, self.__path[i])
        self.__f(blas.daxpy(self.__k1,self.__tmp, a=ho2), t+ho2,\
                    self.__k2)
        np.copyto(self.__tmp, self.__path[i])
        self.__f(blas.daxpy(self.__k2,self.__tmp, a=ho2), t+ho2,\
                    self.__k3)
        np.copyto(self.__tmp, self.__path[i])
        self.__f(blas.daxpy(self.__k3,self.__tmp, a=self.dt), t+self.dt,\
                    self.__k4)
        np.copyto(self.__path[i+1], self.__path[i])
        blas.daxpy(self.__k2,self.__k1,a=2)
        blas.daxpy(self.__k3,self.__k4,a=2)
        blas.daxpy(self.__k1,self.__k4)
        blas.daxpy(self.__k4, self.__path[i+1], a=self.dt/6)
    
    def __runStepC(self, i=0, t0=0):
        """
        Propagate once from step i (Complex version).
        Input:  i           --Time index. (Default: 0)
                t0          --Initial time. (Default: 0)
        """
        
        # Initialize
        t = t0 + i * self.dt
        ho2 = self.dt/2
        
        # 4th order Runge-Kutta
        np.copyto(self.__tmp, self.__path[i])
        self.__f(self.__tmp, t, self.__k1)
        np.copyto(self.__tmp, self.__path[i])
        self.__f(blas.zaxpy(self.__k1,self.__tmp, a=ho2), t+ho2,\
                    self.__k2)
        np.copyto(self.__tmp, self.__path[i])
        self.__f(blas.zaxpy(self.__k2,self.__tmp, a=ho2), t+ho2,\
                    self.__k3)
        np.copyto(self.__tmp, self.__path[i])
        self.__f(blas.zaxpy(self.__k3,self.__tmp, a=self.dt), t+self.dt,\
                    self.__k4)
        np.copyto(self.__path[i+1], self.__path[i])
        blas.zaxpy(self.__k2,self.__k1,a=2)
        blas.zaxpy(self.__k3,self.__k4,a=2)
        blas.zaxpy(self.__k1,self.__k4)
        blas.zaxpy(self.__k4, self.__path[i+1], a=self.dt/6)
    
    
    def __f(self, vecPhase, t, res):
        """
        f function with vector in phase space as i/o.
        Input:  vecPhase    --Vector in phase space.
                t           --Time.
        Output: res         --Results.
        """
        
        res[[0,1,2]] = vecPhase[[3,4,5]]
        res[[3,4,5]] = self.a(vecPhase[[0,1,2]], t)


#=======================================================================
#                               Test
#=======================================================================

def main():
    """Test function, if run separately."""
    oRK = RK4(nt=10000, dt=0.1*asec)
    oRK.run([1,0.2,0,-1.4,0,0])
    with open("Output/test.dat", "w") as f:
        path = oRK.getPath()
        for vecPhase in path:
            f.write("{:f}\t{:f}\n".format(vecPhase[0], vecPhase[1]))
    
	
if __name__ == '__main__':
    main()
