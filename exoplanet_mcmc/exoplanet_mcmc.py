import numpy as np 
import astropy.constants as cons 
import astropy.units as u 
import emcee as mc
import matplotlib.pyplot as plt

class Planet(object):

    def __init__(self, TIMES, DELTA_RAS, DELTA_DECS):
        self.times=TIMES
        self.delta_ras=DELTA_RAS
        self.delta_decs=DELTA_DECS

class System(object):

    def __init__(self, DISTANCE, PLANETS_OBS, RA=0, DEC=0):
        self.ra=RA
        self.dec=DEC
        self.dist=DISTANCE
        self.planets={}
        for i in range(len(PLANETS_OBS[0])):
            self.planets['planet {0}'.format(str(i))]=Planet(PLANETS_OBS[j] for j in range(3))

    def plot_sky_positions(self):
        fig, ax=plt.subplots(1,1)
        ax.scatter(self.ra, self.dec)
        for p in self.planets:
            for i, t in enumerate(p.times):
                RA=self.ra+p.delta_ras[i]
                DEC=self.dec+p.delta_decs[i]
                ax.scatter(RA, DEC, color='k')
                ax.annotate(p.key+' t={0}'.format(str(t)), (RA-p.delta_ras[i]/10., DEC-p.delta_decs[i]/10.))
        return fig