import numpy as np 
import astropy.constants as cons 
import astropy.units as u 
import emcee as mc
import matplotlib.pyplot as plt

class Planet(object):

    def __init__(self, PLANET_OBS=False):
        if PLANET_OBS!=False:
            self.times=PLANET_OBS[0]
            self.delta_ras=PLANET_OBS[1]
            self.delta_decs=PLANET_OBS[2]



class System(object):

    def __init__(self, DISTANCE, PLANETS_OBS=False, RA=0, DEC=0):
        self.ra=RA
        self.dec=DEC
        self.distance=DISTANCE
        self.planets={}
        if PLANETS_OBS!=False:
            for i in range(len(PLANETS_OBS)):
                self.planets['Planet {0}'.format(str(i))]=Planet(PLANETS_OBS[i])

    def plot_sky_positions(self):
        fig, ax=plt.subplots(1,1)
        ax.scatter(self.ra, self.dec, marker='*')
        #offset=np.max(ax.dataLim.bounds)/10.
        ax.annotate('Host Star',
                    xy=(self.ra, self.dec),
                    xytext=(-10, -20),
                    textcoords=('offset pixels'),
                    size=9)
        ax.set_xlabel('RA [mas]')
        ax.set_ylabel('DEC [mas]')
        for key in self.planets.keys():
            planet=self.planets[key]
            for i, t in enumerate(planet.times):
                RA=self.ra+planet.delta_ras[i]
                DEC=self.dec+planet.delta_decs[i]
                ax.scatter(RA, DEC, color='k', s=10)
                ax.annotate(key+' {0}'.format(str(t)),
                            xy=(RA, DEC),
                            xytext = (10, 0),
                            textcoords = ('offset pixels'),
                            size = 9)
        return fig

def generate_sample_system(DISTANCE=5, MSTAR=.95, MPLANET=5, a=10, e=.2, NUMBER_OBS=20):
    true_anomalies=2.*np.pi*np.random.rand(NUMBER_OBS)
    test_system=System(DISTANCE)
    test_system.distance=test_system.distance*u.pc.to('m')
    test_system.mstar=MSTAR*u.Msun.to('kg')
    test_system.planets['Planet 0']=Planet()
    planet0=test_system.planets['Planet 0']
    planet0.mass=MPLANET*u.Mjup.to('kg')
    planet0.a=a*u.au.to('m')
    planet0.e=e
    planet0.r=[]
    planet0.sky_sep=[]
    planet0.theta=true_anomalies
    for theta in true_anomalies:
        r_theta=planet0.a*(1.-planet0.e**2)/(1.+planet0.e*np.cos(theta))
        planet0.r.append(r_theta)
        sep=np.arctan(r_theta/test_system.distance)*u.rad.to('mas')
        planet0.sky_sep.append(sep)
    return test_system
