'''
========================================================================

    File Name:      setup.py
    Author:         Jason Li
    Description:    Input setup file (50 nm Au sphere)
    Usage:          
    
========================================================================
'''

#-----------------------------------------------------------------------
#                           Enviroment setup
#-----------------------------------------------------------------------

rs = 25.0                       # Sphere radius  [nm]
V0 = -13.1                      # Depth of spherical square well [eV]
Fermi = -5.1                    # Fermi level [eV]
tauXUV = 200                    # XUV duration  [as]
omegaXUV = 105                  # XUV frequency  [eV]
IXUV = 1e10                     # XUV intensity [W/cm^2]
tauIR  = 2.47256                # IR duration  [fs]
omegaIR = 1.722                 # IR frequency  [eV]
IIR = 1e12                      # IR intensity  [W/cm^2]
mfp  = 4.4                      # Mean free path  [Angstrom]

#-----------------------------------------------------------------------
#                       Numerical grid setup
#-----------------------------------------------------------------------

calStreaking = False        # Whether to recalculate streaked spectra.
    
rStart = 22.00000001        # Radial integral start point [nm]
rEnd = 25.50000001          # Radial integral end point [nm]
nr = 141                    # Radial integral step number

tauMin = -4                 # XUV-NIR time delay minimal [fs]
tauMax = 4                  # XUV-NIR time delay maximal [fs]
ntau = 81                   # XUV-NIR time delay steps

EPEMin = 60                 # Photoelectron mimumum energy [eV]
EPEMax = 130                # Photoelectron maximum energy [eV]
nEPE = 32                   # Photoelectron energy step number

ntheta = 11                 # Photoelectron direction theta step number
nphi = 10                   # Photoelectron direction phi step number
    
#-----------------------------------------------------------------------
#                       Initial states setup
#-----------------------------------------------------------------------

calInitStatesSearch = False     # Whether to search initial states
calInitStatesSample = False     # Whether to sample initial states
nIS = 10                        # Segments of initial states
mIS = 10                        # Sampled initial states per segment
    
#-----------------------------------------------------------------------
#                       Final states setup
#-----------------------------------------------------------------------

calFinStates = True     # Whether to calculate Volkov phase
nrGV = 7                # Volkov calculation radial step number
nEPEGV = 5              # Volkov calculation final energy step number
ntGV = 21               # Volkov calculation time step number
ntIntGV = 3200          # Volkov calculation time integeral step number\
                        #   (of 4*FWHM)
    
#-----------------------------------------------------------------------
#                           Export setup
#-----------------------------------------------------------------------

ntauEx = 241            # Exported XUV-NIR time delay steps
nEPEEx = 96             # Exported XUV-NIR time delay time steps
    
#-----------------------------------------------------------------------
#                   Field reconstruction setup
#-----------------------------------------------------------------------

reconField = True       # Whether to reconstruct electric field
    
#-----------------------------------------------------------------------
#                           Recyle
#-----------------------------------------------------------------------

calTotalField = True        # Whether to turn on plasmonic field
thetaPE = 90                # Detector (/final moment) theta angle [deg]
phiPE = 0                   # Detector (/final moment) phi angle [deg]

