

aradius = 1.35e-4 # micro m
partradius = 3.0 # micro m

mperpartical = partradius**3/aradius**3
print "mperpartical:\t %e moleculs/partical" % mperpartical

Na = 6.0221413e+23
dens = 2.5e-15 # g/microm^3
molar_m = 241 # g / mol
mass_molecule = molar_m/Na # g
print "mass_molecule:\t %e g" % mass_molecule
mass_particle = mass_molecule*mperpartical # g
print "mass_particle:\t %e g" % mass_particle
part_dens = dens/mass_particle #particles/microm
print "part_dens:\t %e particles/microm^3" % part_dens
volume = 1/part_dens*512
print "volume:\t\t %e microm^3" % volume

