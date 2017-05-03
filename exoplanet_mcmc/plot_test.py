import exoplanet_mcmc.exoplanet_mcmc as exmc
import beta_pic_b as data
import numpy as np

obs=[[data.ts, data.del_RA, data.del_DEC]]
b=exmc.System(19.2, obs)
f=b.plot_sky_positions()
f.show()