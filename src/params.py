from numpy import *

npart = 729
maxt = 20.
ratio = .3
dt = .001
H = 10
MU = .01

print "npart=\t\t %e" % npart
visc = 1.5e-4

print "max_t=\t\t %e" % 20
aradius = 1.35e-4 # micro m
partradius = 3.0/2 # micro m

mperpartical = partradius**3/aradius**3
print "mperpartical=\t %e " % mperpartical

Na = 6.0221413e+23
dens = 2.5e-12 # g/microm^3
molar_m = 241 # g / mol
mass_molecule = molar_m/Na # g
print "mass_molecule=\t %e " % mass_molecule
mass_particle = mass_molecule*mperpartical # g
print "mass_particle=\t %e " % mass_particle
part_dens = dens/mass_particle #particles/microm
print "part_dens=\t %e " % part_dens
volume = 1/part_dens*npart
print "volume=\t\t %e " % volume
print "X=\t\t %e " % (volume)**(1/3.)
print "Y=\t\t %e " % (volume)**(1/3.)
print "Z=\t\t %e " % (volume)**(1/3.)
print "total_mass=\t %e " % (npart*mass_particle)
# print "dt=\t\t %e" % (1e-8*mass_particle*3*pi*visc*partradius*2)
print "checkpoint_interval=\t %e" % 50
print "dt=\t\t %e" % dt
print "H=\t\t %e" % H
print "MU=\t\t %e" % MU
print "ratio=\t\t %e" % ratio
print "maxt=\t\t %e" % maxt

