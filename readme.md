# Exoplanet_MCMC Demonstration

The goal of this notebook is to demonstrate the effect of different properties of astrometric data on the posterior distributions resulting from the MCMC package I have written. I will show the effect of having more or less densely sampled data, greater or fewer points included, and smaller or larger uncertainties on the observations.

## Using the module exoplanet_mcmc

To get to a corner plot of the orbital parameters of a system, only four functions must be used:

### Create an instance of System

```print(exmc.System.__init__.__doc__)
 Instantiates a System
            Arguments [all without units attached]:
                mstar - stellar mass in units of Msun
                mplanet - planetary mass in units of Mjup
                distance - distance to system in parsecs
                semimajor_axis - in au
                eccentricity
                time_of_periastron [between -1 and 1, units of period]
                argument_of_periastron [-2pi to 2pi]
                timesteps - number of timesteps in an orbit

sys=exmc.System(mstar=2., mplanet=2., semimajor_axis=2., eccentricity=.2, distance=5., argument_of_periastron=np.pi/4., time_periastron=0, timesteps=40)

sys.plot_orbit()```

![alt text](https://github.com/jacobhhamer/Exoplanet_MCMC/test_orbit.pdf "Orbit")

### Pick a sample of the orbit to add noise to and use as observations for the MCMC

```print(exmc.System.generate_mcmc_sample.__doc__)
 For an instance of system, generates a sample to run MCMC on.
            Arguments:
                mas_unc - uncertainty in observations in mas, no units attached
                sigma_true_anomaly - uncertainty in observations of position angle in degrees
                sigma_mass - uncertainty in the stellar mass as a fraction of true stellar mass
                indices - slice of the time, true_anomaly, and radius arrays to use
sys.generate_mcmc_sample(indices=[10,20], mas_unc=5., sigma_true_anomaly=3., sigma_mass=.03)

sys.plot_mcmc_sample()```

![alt text](https://github.com/jacobhhamer/Exoplanet_MCMC/test_mcmc_sample.pdf "MCMC Sample")

### Generate the first guess of the parameters for the sample observations

```print(exmc.System.generate_first_guess.__doc__)
 Method to generate the first guess [initial position in parameter space] for the MCMC, by minimizing the
        likelihood function.

sys.generate_first_guess()
```

### Run the MCMC on the sample of the data

```print(exmc.System.runmcmc.__doc__)
 Method to run emcee on an instance of System for which generate_sample_data has been run.
            Arguments:
                p0spread - fraction of range of parameter space which is used as variance of normal distribution to sample from
                nwalker - number of MCMC walkers
                nburn - number of burn-in steps
                nsteps - number of MCMC steps
            Note: For the EnsembleSampler default affine invariant algorithm, it is good to use many walkers and fewer steps.

sys.runmcmc(p0spread=.01, nwalker=500, nburn=200, nsteps=400)

3.8032890891e+30 4.17956644614e+30
0.172149097695 0.230489784341
287480747958.0 313270933742.0
-0.0728754135833 0.0707064628404
-0.344134986222 0.42200270733
Burn in progress:  9.5%
Burn in progress: 19.5%
Burn in progress: 29.5%
Burn in progress: 39.5%
Burn in progress: 49.5%
Burn in progress: 59.5%
Burn in progress: 69.5%
Burn in progress: 79.5%

Burn in progress: 89.5%
Burn in progress: 99.5%
MCMC progress:  9.8%
MCMC progress: 19.8%
MCMC progress: 29.8%
MCMC progress: 39.8%
MCMC progress: 49.8%
MCMC progress: 59.8%
MCMC progress: 69.8%
MCMC progress: 79.8%
MCMC progress: 89.8%
MCMC progress: 99.8%
```

###Plot the results

```sys.walker_plot()```

![alt text](https://github.com/jacobhhamer/Exoplanet_MCMC/walker_plot.pdf "Walker Plot: Should look like white noise if run for enough steps")

```sys.corner_plot()```

![alt text](https://github.com/jacobhhamer/Exoplanet_MCMC/corner_plot.pdf "Corner Plot")


