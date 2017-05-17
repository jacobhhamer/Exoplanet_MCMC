import numpy as np 
import astropy.constants as cons 
import astropy.units as u 
import emcee
import corner
import matplotlib.pyplot as plt
import scipy.optimize as op


def mas_unc_to_m_unc(distance, mas_unc):
    """ Converts mas uncertainties to meters at the plane of the system.
        Arguments:
            distance - distance to the exoplanet in meters
            mas_unc - the uncertainty on the observations in mas
    """
    radians_unc = mas_unc.to('rad')
    m_unc = np.tan(radians_unc) * distance
    return m_unc.to('m')

def lnlike(pars, data, info):
    """ Likelihood function that returns the chi-squared distance between the model, using a set of orbital parameters passed by MCMC, and the observed values.
        Arguments:
            pars - passed during MCMC; tuple of unitless parameters in order:
                total_mass [in kg]
                eccentricity
                semimajor_axis [in m]
                time_of_periastron [between -1 and 1, units of period]
                argument_of_periastron [-2pi to 2pi]
            data - data to run mcmc on; tuple of data with units in order:
                true_anomalies [rad]
                radii [m]
            info - tuple of info in order of:
                time corresponding to data [sec]
                sigma_theta [rad]
                sigma_r [m]
    """
    mpm, e, a, time_periastron, arg_periastron = pars #Unpack parameters
    mpm=mpm*u.kg #Total mass
    a=a*u.m #Semimajor axis
    T=np.sqrt(4.*np.pi**2.*a**3./(cons.G*mpm)) #Period
    times, sigma_theta, sigma_r = info #Unpacking info (time values, uncertainties)
    thetas, rs=data #Unpacking data
    mean_anomalies=(2.*np.pi)*(times-time_periastron*T)/T #Calculating mean anomalies
    chis=np.zeros(len(times)*2) #For each time, two distances to calculate
    for i in range(len(times)):
        k=0.85
        mean_anomaly=mean_anomalies[i]
        eccentric_anomaly_0=mean_anomaly+np.sign(np.sin(mean_anomaly*u.rad))*k*e
        eccentric_anomaly_i=eccentric_anomaly_0
        for j in range(10):
            fppp=e*np.cos(eccentric_anomaly_i*u.rad)
            fpp=e*np.sin(eccentric_anomaly_i*u.rad)
            fp=1.-e*np.cos(eccentric_anomaly_i*u.rad)
            f=eccentric_anomaly_i-e*np.sin(eccentric_anomaly_i*u.rad)-mean_anomaly
            d1=-f/fp
            d2=-f/(fp+.5*d1*fpp)
            d3=-f/(fp+.5*d2*fpp+(1/6.)*d2**2.*fppp)
            eccentric_anomaly_next=eccentric_anomaly_i+d3
            eccentric_anomaly_i=eccentric_anomaly_next
        cos_true_anomaly=(np.cos(eccentric_anomaly_next*u.rad)-e)/(1.-e*np.cos(eccentric_anomaly_next*u.rad))
        true_anomaly=2.*np.arctan(np.sqrt((1+e)/(1-e))*np.tan(eccentric_anomaly_next*u.rad/2.))
        true_r=a*(1.-e**2.)/(1.+e*np.cos(true_anomaly+arg_periastron*u.rad))
        chis[2*i]=(thetas[i]-true_anomaly)/sigma_theta
        chis[2*i+1]=(rs[i]-true_r)/sigma_r
    chis=np.array(chis)
    return -0.5*(np.sum(chis*chis))

class System(object):

    def __init__(self, mstar, mplanet, distance, semimajor_axis, eccentricity, time_periastron, argument_of_periastron, timesteps):
        self.mstar=mstar*cons.M_sun
        self.mplanet=mplanet*cons.M_jup
        self.distance=distance*cons.pc
        self.a=semimajor_axis*cons.au
        self.e=eccentricity
        self.tau=time_periastron
        self.arg=argument_of_periastron*u.rad
        mu=cons.G*(self.mstar+self.mplanet)
        self.T=np.sqrt(4.*np.pi**2.*a**3./mu)
        n=2.*np.pi/self.T
        self.times=np.linspace(0 ,T , timesteps)
        mean_anomalies-n*(self.times-self.tau)
        self.rs=[]
        self.true_anomalies=[]
        for i in range(timesteps):
            k = 0.85  # Suggested by Murray&Dermott
            mean_anomaly = mean_anomalies[i]
            eccentric_anomaly_0 = mean_anomaly + np.sign(np.sin(mean_anomaly * u.rad)) * k * self.e  # First guess
            eccentric_anomaly_i = eccentric_anomaly_0

            for i in range(5):  # For low eccentricities, takes very few (~1) iterations to converge
                fppp = self.e * np.cos(eccentric_anomaly_i * u.rad)
                fpp = self.e * np.sin(eccentric_anomaly_i * u.rad)
                fp = 1. - self.e * np.cos(eccentric_anomaly_i * u.rad)
                f = eccentric_anomaly_i - self.e * np.sin(eccentric_anomaly_i * u.rad) - mean_anomaly
                d1 = -f / fp
                d2 = -f / (fp + .5 * d1 * fpp)
                d3 = -f / (fp + .5 * d2 * fpp + (1 / 6.) * d2 ** 2. * fppp)

                eccentric_anomaly_next = eccentric_anomaly_i + d3
                eccentric_anomaly_i = eccentric_anomaly_next

            cos_true_anomaly = (np.cos(eccentric_anomaly_next * u.rad) - self.e) / (1. - self.e * np.cos(eccentric_anomaly_next * u.rad))
            r = self.a * (1. - self.e ** 2.) / (1. + self.e * cos_true_anomaly)
            true_anomaly = 2. * np.arctan(np.sqrt((1 + self.e) / (1 - self.e)) * np.tan(eccentric_anomaly_next * u.rad / 2.))

            self.rs.append(r.value)
            self.true_anomalies.append(true_anomaly.value)

            self.rs=np.array(self.rs)*u.m
            self.true_anomalies=np.array(self.true_anomalies)*u.rad

    def plot_orbit(self):
        """ Method to plot the orbit, with no noise added. Called on an instance of System.
        """
        fig, ax=plt.subplots(1,1)
        ax.set_aspect('equal')
        ax.scatter(0, 0, marker='*', color='k', s=30)
        for i, t in enumerate(self.times):
            ax.scatter(self.rs[i] * np.cos(self.true_anomalies[i] + self.arg), self.rs[i] * np.sin(self.true_anomalies[i] + self.arg), color='c', s=5)
        return fig

    def generate_mcmc_sample(self, mas_unc, sigma_true_anomaly, sigma_mass, indices):
        """ For an instance of system, generates a sample to run MCMC on.
            Arguments:
                mas_unc - uncertainty in observations in mas, no units attached
                sigma_true_anomaly - uncertainty in observations of position angle in degrees
                sigma_mass - uncertainty in the stellar mass as a fraction of true stellar mass
                indices - slice of the time, true_anomaly, and radius arrays to use
        """
        self.sigma_r=mas_unc_to_m_unc(self.distance, mas_unc*u.mas)
        self.sigma_true_anom=(sigma_true_anomaly*u.deg).to('rad')
        sigma_mass=sigma_mass*self.mstar
        self.noisy_rs=self.rs[indices[0], indices[1]]+sigma_r*np.random.normal(loc=0, size=(indices[1]-indices[0]))
        self.noisy_true_anomalies=self.true_anomalies[indices[0], indices[1]]+sigma_true_anom*np.random.normal(loc=0, size=(indices[1]-indices[0]))
        self.noisy_mass=self.mstar+sigma_mass*np.random.normal(loc=0, size=1)
        self.mcmc_times=self.times[indices[0], indices[1]]

    def plot_mcmc_sample(self):
        """ Method to plot the data which the MCMC will run on. Called on an instance of System for which generate_mcmc_sample has been run.
        """
        fig, ax=plt.subplots(1,1)
        ax.set_aspect('equal')
        ax.scatter(0, 0, marker='*', color='k', s=30)
        for i, t in enumerate(self.times):
            ax.scatter(self.noisy_rs[i] * np.cos(self.noisy_true_anomalies[i] + self.arg), self.noisy_rs[i] * np.sin(self.noisy_true_anomalies[i] + self.arg), color='c', s=5)
        return fig

    def generate_first_guess(self):
        """ Method to generate the first guess [initial position in parameter space] for the MCMC, by minimizing the likelihood function.
        """
        nll = lambda *args: -lnlike(*args)
        result = op.minimize(nll, [(self.noisy_mass + self.mplanet).value, self.e, self.a.value, self.tau, self.arg], (self.data, self.info))
        minv = result["x"]
        minv[0] = minv[0][0]
        minv = np.array(minv, dtype=float)
        self.pos0=minv

    def lnprior(self, pars):
        """ Flat priors in each parameter, inspired by those used for Beta Pic b in (Wang et al. 2016: 1607.05272)
            Arguments:
            pars - passed during MCMC; tuple of UNITLESS parameters in order:
                total_mass [in kg]
                eccentricity
                semimajor_axis [in m]
                time_of_periastron [between -1 and 1]
                argument_of_periastron [-2pi to 2pi]
        """
        mpm, e, a, time_periastron, arg_periastron = pars
        if ((self.noisy_mass*.5).value < mpm < (self.noisy_mass*2.).value) and (0.0001 < e < 0.99) and ((self.a/2.).value < a < (self.a*2.).value) and (-1 < time_periastron < 1) and ( -2 * np.pi < arg_periastron < 2 * np.pi):
            return 0.0
        return -np.inf

    def lnprob(self, pars, data, info):
        """ Likelihood function that returns the chi-squared distance between the model, using a set of orbital parameters passed by MCMC, and the observed values.
            Arguments:
                pars - passed during MCMC; tuple of unitless parameters in order:
                    total_mass [in kg]
                    eccentricity
                    semimajor_axis [in m]
                    time_of_periastron [between -1 and 1, units of period]
                    argument_of_periastron [-2pi to 2pi]
                data - data to run mcmc on; tuple of data with units in order:
                    true_anomalies [rad]
                    radii [m]
                info - tuple of info in order of:
                    time corresponding to data [sec]
                    sigma_theta [rad]
                    sigma_r [m]
        """
        lp = self.lnprior(pars)
        if ~(np.isfinite(lp)):
            return -np.inf
        return lp + lnlike(pars, data, info)

    def runmcmc(self, p0spread, nwalker, nsteps):
        """ Method to run emcee on an instance of System for which generate_sample_data has been run.
            Arguments:
                p0spread - fraction of range of parameter space which is used as variance of normal distribution to sample from
                nwalker - number of MCMC walkers
                nsteps - number of MCMC steps
        """
        self.data=[self.noisy_true_anomalies, self.noisy_rs]
        self.info=[self.mcmc_times, self.sigma_true_anom, self.sigma_r]
        self.generate_first_guess()
        ndim=len(self.pos0)
        pos = [self.pos0 + p0spread * np.array([1.5*self.noisy_mass.value, 1., 1.5*self.a.value, 2., 4. * np.pi]) * np.random.randn(ndim) for i in range(nwalkers)]
        sampler = emcee.EnsembleSampler(nwalkers=nwalker, dim=ndim, lnpostfn=lnprob, args=(self.data, self.info))
        sampler.run_mcmc(pos0=pos, N=nsteps)
        self.mcmc=sampler

    def walker_plot(self):
        """ Method to plot the walkers in each parameter space. Called on an instance of System for which runmcmc has been run.
        """
        f, axs = plt.subplots(self.pos0.shape[0], 1, figsize=(6, 12))
        for i, ax in enumerate(f.axes):
            for j in range(self.mcmc.chain.shape[0]):
                ax.plot(self.mcmc.chain[j, :, i], color='k', alpha=0.2)
                ax.set_ylabel(['M', 'e', 'a', 'epoch', 'arg'][i])
        return f

    def corner_plot(self):
        """ Method to plot the corner plot for the MCMC. Called on an instance of System for which runmcmc has been run.
        """
        samples = self.mcmc.chain[:, :, :].reshape((-1, self.pos0.shape[0]))
        fig = corner.corner(samples[:, :], labels=["M", "e", "a", "tau", "arg"],
                            truths=[(self.mstar + self.mplanet).value, self.e, self.a.value, self.tau, self.arg])
        return fig