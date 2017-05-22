import numpy as np
import astropy.constants as cons
import astropy.units as u
import emcee
import corner
import matplotlib.pyplot as plt
import scipy.optimize as op
from matplotlib.ticker import FormatStrFormatter


def mas_unc_to_m_unc(distance, mas_unc):
    """ Converts mas uncertainties to meters at the plane of the system.
        Arguments:
            distance - distance to the exoplanet in meters
            mas_unc - the uncertainty on the observations in mas"""
    radians_unc = mas_unc.to('rad') #convert milli-arcseconds to radians
    m_unc = np.tan(radians_unc) * distance #convert radians to parsec in the plane of the system
    return m_unc.to('m') #convert parsecs to meters


def lnlike(pars, data, info):
    """ Likelihood function that returns the chi-squared distance between the model, using a set of orbital parameters
     passed by MCMC, and the observed values.
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
                sigma_r [m]"""
    mpm, e, a, time_periastron, arg_periastron = pars  # Unpack parameters
    times, sigma_theta, sigma_r = info  # Unpacking info
    thetas, rs = data  # Unpacking data

    mpm = mpm * u.kg  # Total mass
    a = a * u.m  # Semimajor axis
    T = np.sqrt(4. * np.pi ** 2. * a ** 3. / (cons.G * mpm))  # Period
    mean_anomalies = (2. * np.pi) * (times - time_periastron * T) / T  # Calculating mean anomalies
    chis = np.zeros(len(times) * 2)  # For each time, two distances to calculate

    #Iterative scheme to solve Kepler's Eq. as a function of time
    #described on Solar System Dynamics pg 36(Murray, Dermott)
    for i in range(len(times)):
        k = 0.85 # Suggested by Murray&Dermott
        mean_anomaly = mean_anomalies[i]
        eccentric_anomaly_0 = mean_anomaly + np.sign(np.sin(mean_anomaly * u.rad)) * k * e #First guess
        eccentric_anomaly_i = eccentric_anomaly_0 

        for j in range(3): # For low eccentricities, takes very few (~1) iterations to converge
            fppp = e * np.cos(eccentric_anomaly_i * u.rad)
            fpp = e * np.sin(eccentric_anomaly_i * u.rad)
            fp = 1. - e * np.cos(eccentric_anomaly_i * u.rad)
            f = eccentric_anomaly_i - e * np.sin(eccentric_anomaly_i * u.rad) - mean_anomaly
            d1 = -f / fp
            d2 = -f / (fp + .5 * d1 * fpp)
            d3 = -f / (fp + .5 * d2 * fpp + (1 / 6.) * d2 ** 2. * fppp)
            eccentric_anomaly_next = eccentric_anomaly_i + d3
            eccentric_anomaly_i = eccentric_anomaly_next

        true_anomaly = 2. * np.arctan(np.sqrt((1 + e) / (1 - e)) * np.tan(eccentric_anomaly_next * u.rad / 2.)) #Angle around orbit
        true_r = a * (1. - e ** 2.) / (1. + e * np.cos(true_anomaly + arg_periastron * u.rad)) #Radial distance
        chis[2 * i] = (thetas[i] - true_anomaly) / sigma_theta #distance between model and observations
        chis[2 * i + 1] = (rs[i] - true_r) / sigma_r #distance between model and observations
    chis = np.array(chis)
    return -0.5 * (np.sum(chis * chis)) #chi-sruared distance

def random_orbit():
    """ Instantiates a system with random argument values."""
    sys=System(mstar=np.random.randint(1,5), mplanet=np.random.randint(1, 10), distance=5., semimajor_axis=np.random.randint(1,25), eccentricity=np.random.random(), argument_of_periastron=2.*np.pi*np.random.random(), time_periastron=0, timesteps=25)
    return sys

class System(object):
    def __init__(self, mstar, mplanet, distance, semimajor_axis, eccentricity, time_periastron, argument_of_periastron, timesteps):
        """ Instantiates a System
            Arguments [all without units attached]:
                mstar - stellar mass in units of Msun
                mplanet - planetary mass in units of Mjup
                distance - distance to system in parsecs
                semimajor_axis - in au
                eccentricity
                time_of_periastron [between -1 and 1, units of period]
                argument_of_periastron [-2pi to 2pi]
                timesteps - number of timesteps in an orbit"""
        self.mstar = mstar * cons.M_sun #Attach units
        self.mplanet = mplanet * cons.M_jup #Attach units
        self.distance = distance * cons.pc #Attach units
        self.a = semimajor_axis * cons.au #Attach units
        self.e = eccentricity
        self.tau = time_periastron 
        self.arg = argument_of_periastron * u.rad #Attach units

        mu = cons.G * (self.mstar + self.mplanet)
        self.T = np.sqrt(4. * np.pi ** 2. * self.a ** 3. / mu)
        n = 2. * np.pi / self.T #Mean motion
        self.times = np.linspace(0, self.T, timesteps)
        mean_anomalies = n * (self.times - self.tau)
        rs = []
        true_anomalies = []

        #Iterative scheme to solve Kepler's Eq. as a function of time
        #described on Solar System Dynamics pg 36(Murray, Dermott)
        for i in range(timesteps):
            k = 0.85  # Suggested by Murray&Dermott
            mean_anomaly = mean_anomalies[i]
            eccentric_anomaly_0 = mean_anomaly + np.sign(np.sin(mean_anomaly * u.rad)) * k * self.e  # First guess
            eccentric_anomaly_i = eccentric_anomaly_0

            for i in range(3):  # For low eccentricities, takes very few (~1) iterations to converge
                fppp = self.e * np.cos(eccentric_anomaly_i * u.rad)
                fpp = self.e * np.sin(eccentric_anomaly_i * u.rad)
                fp = 1. - self.e * np.cos(eccentric_anomaly_i * u.rad)
                f = eccentric_anomaly_i - self.e * np.sin(eccentric_anomaly_i * u.rad) - mean_anomaly
                d1 = -f / fp
                d2 = -f / (fp + .5 * d1 * fpp)
                d3 = -f / (fp + .5 * d2 * fpp + (1 / 6.) * d2 ** 2. * fppp)

                eccentric_anomaly_next = eccentric_anomaly_i + d3
                eccentric_anomaly_i = eccentric_anomaly_next

            cos_true_anomaly = (np.cos(eccentric_anomaly_next * u.rad) - self.e) / (
            1. - self.e * np.cos(eccentric_anomaly_next * u.rad))
            r = self.a * (1. - self.e ** 2.) / (1. + self.e * cos_true_anomaly)
            true_anomaly = 2. * np.arctan(
                np.sqrt((1 + self.e) / (1 - self.e)) * np.tan(eccentric_anomaly_next * u.rad / 2.))

            rs.append(r.value)
            true_anomalies.append(true_anomaly.value)

            self.rs = np.array(rs) * u.m
            self.true_anomalies = np.array(true_anomalies) * u.rad

    def plot_orbit(self):
        """ Method to plot the orbit, with no noise added. Called on an instance of System."""
        fig, ax = plt.subplots(1, 1)
        ax.set_aspect('equal')
        ax.scatter(0, 0, marker='*', color='k', s=100, label='Host Star') #Plot host star at 0,0
        ax.scatter(self.rs.to('au') * np.cos(self.true_anomalies + self.arg),
                   self.rs.to('au') * np.sin(self.true_anomalies + self.arg), color='#ff7119', s=25, label='Exoplanet') #Plot the exoplanet position at each timestep
        lims=ax.dataLim.get_points() #Get the extrema of the data
        ymin=lims[0][1]
        ymax=lims[1][1]
        xmin=lims[0][0]
        xmax=lims[1][0]
        ax.set_xticks([xmin, 0, xmax]) #Only plot the extrema of the orbit and 0
        ax.set_yticks([ymin, 0, ymax]) #Only plot the extrema of the orbit and 0
        ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'+' AU')) #Set 1 decimal precison and add AU to labels
        ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'+' AU')) #Set 1 decimal precison and add AU to labels
        handles, labels=ax.get_legend_handles_labels()
        for key in ax.spines:
            ax.spines[key].set_visible(False) #Remove borders
        circle=plt.Circle((0,0), radius=1., color='g', fill=False, linewidth=2) #Plot circle of 1AU for reference
        ax.add_artist(circle)
        ax.grid(True, alpha=0.4) #Show grid
        ax.legend(handles=handles+[circle], labels=['Host Star', 'Exoplanet', "Earth's Orbit"], loc=1)   
        return fig

    def generate_mcmc_sample(self, mas_unc, sigma_true_anomaly, sigma_mass, indices):
        """ For an instance of system, generates a sample to run MCMC on.
            Arguments:
                mas_unc - uncertainty in observations in mas, no units attached
                sigma_true_anomaly - uncertainty in observations of position angle in degrees
                sigma_mass - uncertainty in the stellar mass as a fraction of true stellar mass
                indices - slice of the time, true_anomaly, and radius arrays to use
        """
        self.sigma_r = mas_unc_to_m_unc(self.distance, mas_unc * u.mas)
        self.sigma_true_anom = (sigma_true_anomaly * u.deg).to('rad')
        sigma_mass = sigma_mass * self.mstar
        self.noisy_rs = self.rs[indices[0]: indices[1]] + self.sigma_r * np.random.normal(loc=0,
                                                                                     size=(indices[1] - indices[0])) #Add gaussian noise to radial distance
        self.noisy_true_anomalies = self.true_anomalies[indices[0]: indices[1]] + self.sigma_true_anom * np.random.normal(
            loc=0, size=(indices[1] - indices[0])) #Add gaussian noise to angles
        self.noisy_mass = self.mstar + sigma_mass * np.random.randn() #Add gaussian noise to stellar mass
        self.mcmc_times = self.times[indices[0]: indices[1]] #Select times of mcmc data

    def plot_mcmc_sample(self):
        """ Method to plot the data which the MCMC will run on. Called on an instance of System for which
        generate_mcmc_sample has been run. """
        fig, ax = plt.subplots(1, 1)
        ax.set_aspect('equal')
        ax.scatter(0, 0, marker='*', color='k', s=100, label='Host Star')
        keep=np.logical_and(self.times>=self.mcmc_times.min(),self.times<=self.mcmc_times.max())
        ax.scatter(self.noisy_rs.to('au') * np.cos(self.noisy_true_anomalies + self.arg),
                   self.noisy_rs.to('au') * np.sin(self.noisy_true_anomalies + self.arg), color='#ff7119', s=25, label='Noise Added (Observations)')
        ax.scatter(self.rs[keep].to('au') * np.cos(self.true_anomalies[keep] + self.arg),
                   self.rs[keep].to('au') * np.sin(self.true_anomalies[keep] + self.arg), color='k', s=25, alpha=.3, label='True Values')
        handles, labels=ax.get_legend_handles_labels()
        lims=ax.dataLim.get_points() #Get the extrema of the data
        ymin=lims[0][1]
        ymax=lims[1][1]
        xmin=lims[0][0]
        xmax=lims[1][0]
        ax.set_xticks([xmin, 0, xmax]) #Only plot the extrema of the orbit and 0
        ax.set_yticks([ymin, 0, ymax])
        ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'+' AU')) #Set 1 decimal precison and add AU to labels
        ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'+' AU'))
        for key in ax.spines:
            ax.spines[key].set_visible(False) #Remove borders
        circle=plt.Circle((0,0), radius=1., color='g', fill=False, linewidth=2) #Plot circle of 1AU for reference
        ax.add_artist(circle)
        ax.grid(True, alpha=0.4) #Show grid
        ax.legend(handles=handles+[circle], labels=['Host Star', 'Noise Added (Observations)', 'True Values', "Earth's Orbit"], loc=1)   
        return fig


    def generate_first_guess(self):
        """ Method to generate the first guess [initial position in parameter space] for the MCMC, by minimizing the
        likelihood function. """
        self.data = [self.noisy_true_anomalies, self.noisy_rs] #To be passed to probability functions
        self.info = [self.mcmc_times, self.sigma_true_anom, self.sigma_r] #To be passed to probability functions
        nll = lambda *args: -lnlike(*args) #Dummy function to pass to minimize
        result = op.minimize(fun=nll, x0=[(self.noisy_mass + self.mplanet).value, self.e, self.a.value, self.tau, self.arg.value],
                             args=(self.data, self.info), bounds=[(0, 100.*self.noisy_mass.value), (.001, .999), (0, 100*self.a.value), (-1,1),(-2.*np.pi, 2*np.pi)])
        #Minimize within reasonable bounds
        minv = result["x"] #Call the resulting parameters
        minv = np.array(minv, dtype=float)
        self.pos0 = minv

    def lnprior(self, pars):
        """ Flat priors in each parameter, inspired by those used for Beta Pic b in (Wang et al. 2016: 1607.05272)
            Arguments:
            pars - passed during MCMC; tuple of UNITLESS parameters in order:
                total_mass [in kg]
                eccentricity
                semimajor_axis [in m]
                time_of_periastron [between -1 and 1]
                argument_of_periastron [-2pi to 2pi]"""
        mpm, e, a, time_periastron, arg_periastron = pars #Unpack parameters
        #Flat priors in all 5 parameters
        if ((self.noisy_mass * .5).value < mpm < (self.noisy_mass * 2.).value) and (0.0001 < e < 0.99) and (
                (self.a / 2.).value < a < (self.a * 2.).value) and (-1 < time_periastron < 1) and (
                    -2 * np.pi < arg_periastron < 2 * np.pi):
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
                    sigma_r [m]"""
        lp = self.lnprior(pars)
        if ~(np.isfinite(lp)):
            return -np.inf
        return lp + lnlike(pars, data, info)

    def runmcmc(self, p0spread, nwalker, nburn, nsteps):
        """ Method to run emcee on an instance of System for which generate_sample_data has been run.
            Arguments:
                p0spread - fraction of range of parameter space which is used as variance of normal distribution to sample from
                nwalker - number of MCMC walkers
                nburn - number of burn-in steps
                nsteps - number of MCMC steps
            Note: For the EnsembleSampler default affine invariant algorithm, it is good to use many walkers and fewer steps."""
        ndim = len(self.pos0) #Dimensions are length of initial parameter space position vector
        pos = np.array([self.pos0 + p0spread * np.array(
            [1.5 * self.noisy_mass.value, 1., 1.5 * self.a.value, 2., 4. * np.pi]) * np.random.randn(ndim) for i in
               range(nwalker)]) #Define initial position for each walker using gaussian noise with variance of p0spread * (range of prior in each dimension)
        for i in range(ndim):
            print(pos[:,i].min(), pos[:,i].max()) #Print the min and max of each parameter
        sampler = emcee.EnsembleSampler(nwalkers=nwalker, dim=ndim, lnpostfn=self.lnprob, args=(self.data, self.info))
        for i, result in enumerate(sampler.sample(pos, iterations=nburn)):
            if i==nburn-1: #Save end result of burn-in steps
                newpos, prob, state = result
            if (i+1) % (nburn/10.) == 0: #Print progress every 10 percent
                print("Burn-in progress: {0:5.1%}".format(float(i) / nburn)) #With a total of 5 characters and 1 decimal precision, print percent progress
        sampler.reset() #Dump burn-in information
        for i, result in enumerate(sampler.sample(newpos, iterations=nsteps)): #Start from where the walkers ended up after burn-in
            if (i+1) % (nsteps/10.) == 0: #Print progress every 10 percent
                print("MCMC progress: {0:5.1%}".format(float(i) / nsteps)) #With a total of 5 characters and 1 decimal precision, print percent progress
        self.mcmc = sampler
        if (np.mean(sampler.acceptance_fraction) < .25) or (np.mean(sampler.acceptance_fraction) > .5):
            print("The mean acceptance fraction was {0:.3f}. This value suggests the MCMC was not tuned well.".format(np.mean(sampler.acceptance_fraction)))


    def walker_plot(self):
        """ Method to plot the walkers in each parameter space. Called on an instance of System for which runmcmc has been run."""
        f, axs = plt.subplots(self.pos0.shape[0], 1, figsize=(4, 12)) #a subplot for each dimension
        for i, ax in enumerate(f.axes):
            for j in range(self.mcmc.chain.shape[0]): #Plot each walker
                ax.plot(self.mcmc.chain[j, :, i], color='k', alpha=0.2) #As indices are chain[dimension, steps, walker]
                ax.set_ylabel(['M', 'e', 'a', 'tau', 'arg'][i])
        return f

    def corner_plot(self):
        """ Method to plot the corner plot for the MCMC. Called on an instance of System for which runmcmc has been run."""
        samples = self.mcmc.chain[:, :, :].reshape((-1, self.pos0.shape[0])) #Add all of the walkers together, so you have achain for each dimension
        fig = corner.corner(samples[:, :], labels=["M", "e", "a", "tau", "arg"],
                            truths=[(self.mstar + self.mplanet).value, self.e, self.a.value, self.tau, self.arg.value])
        return fig
