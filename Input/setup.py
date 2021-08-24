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

rs = 77                       # Sphere radius  [nm]
ri = 57                         # Shell inner radius (0 if solid)  [nm]
V0 = -13.1                      # Depth of spherical square well [eV]
Fermi = -5.1                    # Fermi level [eV]
tauXUV = 200                    # XUV duration  [as]
omegaXUV = 50                   # XUV frequency  [eV]
IXUV = 1e10                     # XUV intensity [W/cm^2]
tauIR  = 2.47256                # IR duration  [fs]
lambdaIR = 780                  # IR wavelength  [nm]
IIR = 1e10                      # IR intensity  [W/cm^2]
mfp  = 4.4                      # Mean free path  [Angstrom]

#-----------------------------------------------------------------------
#                       Numerical grid setup
#-----------------------------------------------------------------------

calStreaking = False        # Whether to recalculate streaked spectra.
    
rStart = 47.00000001        # Radial integral start point [nm]
rEnd = 50.50000001          # Radial integral end point [nm]
nr = 141                    # Radial integral step number

tauMin = -4                 # XUV-NIR time delay minimal [fs]
tauMax = 4                  # XUV-NIR time delay maximal [fs]
ntau = 81                   # XUV-NIR time delay steps

EPEMin = 25                 # Photoelectron mimumum energy [eV]
EPEMax = 60                 # Photoelectron maximum energy [eV]
nEPE = 32                   # Photoelectron energy step number
    
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

calFinStates = False    # Whether to calculate Volkov phase
nrGV = 7                # Volkov calculation radial step number
nEPEGV = 5              # Volkov calculation final energy step number
ntGV = 21               # Volkov calculation time step number
ntIntGV = 3200          # Volkov calculation time integeral step number\
                        #   (of 4*FWHM)

ntheta = 6              # PE direction theta step number (4)
nphi = 8                # PE direction phi step number (4)
thetaMin = 45           # Minimal theta [deg]
thetaMax = 135          # Maximal theta [deg]
    
#-----------------------------------------------------------------------
#                           Export setup
#-----------------------------------------------------------------------

ntauEx = 241            # Exported XUV-NIR time delay steps (241/81)
nEPEEx = 96             # Exported XUV-NIR time delay time steps (96/64)
nthetaEx = 151          # Reconstruction theta step number (101/151)
nphiEx = 181            # Reconstruction phi step number (121/181)
    
#-----------------------------------------------------------------------
#                   Field reconstruction setup
#-----------------------------------------------------------------------

reconField = True       # Whether to reconstruct electric field
anaField = False        # Whether to calculate analytical electric field

