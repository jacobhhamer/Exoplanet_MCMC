from exoplanet_mcmc.exoplanet_mcmc import *
import numpy as np

def test_system1():
	sys=random_orbit()
	assert sys.mstar.value>0

def test_system2():
	sys=random_orbit()
	assert sys.mplanet.value>0

def test_system3():
	sys=random_orbit()
	assert (sys.e>0)&(sys.e<1)

def test_system4():
	sys=random_orbit()
	assert sys.a.value>0

def test_system5():
	sys=random_orbit()
	assert (sys.tau>-1)&(sys.tau<1)

def test_system6():
	sys=random_orbit()
	assert (sys.arg.value>-2*np.pi)&(sys.arg.value<2*np.pi)