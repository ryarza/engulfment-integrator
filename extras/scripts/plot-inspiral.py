import numpy as np
import matplotlib.pyplot as plt
import argparse
import rytools

cgs = rytools.units.get_cgs()

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--inputfile", help = "Input file", required = True)
parser.add_argument("-s", "--scalarsfile", help = "Scalars file", required = True)
parser.add_argument("-p", "--profile", help = "Stellar profile used in integration, in MESA form (not HDF5!)", required = True)
parser.add_argument("-od", "--outputdir", help = "Directory to store plots", required = False, default = './')
args = parser.parse_args()


d = rytools.planets.load_orbint_output(args.inputfile, args.scalarsfile, profile_path = args.profile, compute_derived = {'deltae_from_forces': True, 'L_emergent': True})

print("Inspiral time:                    %.2e yr"    % ( d['time'][-1] / cgs['YEAR'] ) )
print("Number of orbits:                 %i" % ( d['theta'][-1] / 2 / np.pi ) )
print("Initial orbital energy:           %.2e erg"   % d['eorb'][0])
print("Final orbital energy:             %.2e erg"   % d['eorb'][-1])
print("Final separation:                 %.2e Rsun"  % ( d['r'][-1] / cgs['RSUN'] ))
print("Density at final separation:      %.2e g/cc"  % ( d['rho'][-1] ))
print("Maximum energy deposition rate:   %.2e erg/s" % ( np.amax(-d['edot_from_forces']['total'])))

fig, ax = plt.subplots()
ax.plot(d['x'] / cgs['RSUN'], d['y'] / cgs['RSUN'])
ax.set_xlabel(r'$x/R_\odot$')
ax.set_ylabel(r'$y/R_\odot$')
plt.savefig('x_vs_y.png')
plt.close()

fig, ax = plt.subplots()
ax.loglog()
ax.plot(d['time'] / cgs['YEAR'], d['r'] / cgs['RSUN'])
ax.set_xlabel(r'Time in years')
ax.set_ylabel(r'$r/R_\odot$')
plt.savefig('t_vs_r.png')
plt.close()

fig, ax = plt.subplots()
ax.loglog()
ax.plot(d['r'] / cgs['RSUN'], - d['edot_from_forces']['total'])
ax.set_xlabel(r'$r/R_\odot$')
ax.set_ylabel(r'Energy deposition rate in erg/s')
plt.savefig('a_vs_edot.png')
plt.close()

fig, ax = plt.subplots()
ax.loglog()
ax.plot(d['r'] / cgs['RSUN'], d['rsb'] / d['ra'])
ax.set_xlabel(r'$r/R_\odot$')
ax.set_ylabel(r'$R_{sb}/R_a$')
plt.savefig('a_vs_rsb_over_ra.png')
plt.close()

fig, ax = plt.subplots()
ax.loglog()
ax.plot(d['time'] / cgs['YEAR'], d['L_emergent'])
ax.axhline(np.amax(d['prof']['L']), color = 'black', ls = 'dashed', label = 'Stellar luminosity')
ax.set_xlabel('Time in years')
ax.set_ylabel('Luminosity in erg/s')
ax.legend()
plt.savefig('t_vs_L.png')
plt.close()
