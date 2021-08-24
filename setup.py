'''
========================================================================

    File Name:      setup.py
    Author:         Jason Li
    Description:    Read setup files, and convert units.
    Usage:          
    
========================================================================
'''

# Read material index of refraction
from Input.setup import *
from Input.nAu import *
from Input.nSiO2 import *
from Input.nVac import *
from const import *

#-----------------------------------------------------------------------
#                       Enviroment setup
#-----------------------------------------------------------------------

rs = rs*nm                  # Sphere radius  [nm]
ri = ri*nm                  # Shell inner radius (0 if solid)  [nm]
V0 = V0*eV                  # Depth of spherical square well [eV]
Fermi = Fermi*eV            # Fermi level (negative number) [eV]
tauXUV = tauXUV*asec        # XUV duration  [as]
omegaXUV = omegaXUV*eV      # XUV frequency  [eV]
IXUV = IXUV*Wpcm2           # XUV intensity [W/cm^2]
tauIR  = tauIR*fsec         # IR duration  [fs]
lambdaIR = lambdaIR*nm      # IR wavelength  [nm]
IIR = IIR*Wpcm2             # IR intensity  [W/cm^2]
mfp  = mfp*Ang              # Mean free path  [Angstrom]

bXUV = 1.386294362/tauXUV**2    # XUV b factor square
omegaIR = 2*Pi*c/lambdaIR       # IR frequency  [eV]

#-----------------------------------------------------------------------
#                       Numerical grid setup
#-----------------------------------------------------------------------
    
rStart = rStart*nm              # Radial integral start point [nm]
rEnd = rEnd*nm                  # Radial integral end point [nm]
dr = (rEnd - rStart) / (nr-1)   # Radial integral interval

tauMin = tauMin*fsec            # XUV-NIR time delay minimal [fs]
tauMax = tauMax*fsec            # XUV-NIR time delay maximal [fs]
dtau = (tauMax-tauMin)/(ntau-1) # XUV-NIR time delay interval [fs]

EPEMin = EPEMin*eV          # Photoelectron mimumum energy [eV]
EPEMax = EPEMax*eV          # Photoelectron maximum energy [eV]

nangle = ntheta*nphi        # Photoelectron angluar step number
    
#-----------------------------------------------------------------------
#                       Initial states setup
#-----------------------------------------------------------------------

nISTot = nIS*mIS            # Segments of initial states
    
#-----------------------------------------------------------------------
#                       Final states setup
#-----------------------------------------------------------------------

thetaMin = thetaMin * deg
thetaMax = thetaMax * deg
