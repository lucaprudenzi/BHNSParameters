import configparser

import plot_binary_par as pbp
import error_module as em
import plot_ringdown_par as prp

# Configuration parameter file
config = configparser.ConfigParser()
config.read('config.ini')

# See [PHENOM_MIX]
pbp.plotPhenomMix(config)

# See [THETA_F]
pbp.plotThetaF(config)

# See [CHI_F_VARTHETA]
pbp.plotChiF_VarTheta(config)

# See [CHI_F_VARCHI]
pbp.plotChiF_VarChi(config)

# See [M_OUT]
pbp.plotMOut(config)

# See [M_F]
pbp.plotMf(config)

# Require Ringdown module
# See [SIGMA_MF]
em.sigmaMF(config)

# See [PARAMETER_ZONES] 
pbp.ParameterZones(config)

#Require Ringdown module
# See [PLOT_QNM]
#prp.PlotQNM(config)

# See [PLOT_AMP_QNM] par
#prp.PlotAmplitude(config)
