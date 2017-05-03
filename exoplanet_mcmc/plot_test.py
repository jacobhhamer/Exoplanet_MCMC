import exoplanet_mcmc.exoplanet_mcmc as exmc
from exoplanet_mcmc.data import beta_pic_b as dat
import numpy as np

obs=[[dat.ts, dat.del_RA, dat.del_DEC]]
b=exmc.System(19.2, obs)
f=b.plot_sky_positions()
f.show()