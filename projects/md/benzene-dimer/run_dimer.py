from openmm.app import *
from openmm import *
from openmm.unit import *
from sys import stdout
import glob

## This is a useful tool that will allow us to 
#  save the energies from each term in the force field
#  Do not edit this
class ForceReporter(object):
    def __init__(self, file, reportInterval):
        self._out = open(file, 'w')
        self._reportInterval = reportInterval

        for i, f in enumerate(system.getForces()):
            if i in [0,1,2,3,6]:
                state = simulation.context.getState(getEnergy=True, groups={i})
                self._out.write("%25s\t" % f.getName())
        self._out.write("\n")

    def __del__(self):
        self._out.close()

    def describeNextReport(self, simulation):
        steps = self._reportInterval - simulation.currentStep%self._reportInterval
        return  {'steps': steps, 'periodic': None, 'include':['forces']}

    def report(self, simulation, state):
        #forces = state.getForces().value_in_unit(kilojoules/mole/nanometer)
        #for f in forces:
        #    self._out.write('%g %g %g\n' % (f[0], f[1], f[2]))

        for i, f in enumerate(system.getForces()):
            if i in [0,1,2,3,6]:
                state = simulation.context.getState(getEnergy=True, groups={i})
                self._out.write("%25s\t" % state.getPotentialEnergy())
                #print(f.getName(), state.getPotentialEnergy())
        self._out.write("\n")

psf = CharmmPsfFile('bz-dimer.psf')
pdb = PDBFile('bz-dimer.pdb')


pfiles = glob.glob('toppar/*')
params = CharmmParameterSet(*pfiles)

system = psf.createSystem(params, nonbondedMethod=NoCutoff,
        nonbondedCutoff=None, constraints=HBonds)


## This block of code adds a spherical potential
#  around our system, to prevent our molecules from
#  flying apart and never reuniting 
force = CustomExternalForce('100*max(0, r-2)^2; r=sqrt(x*x+y*y+z*z)')
system.addForce(force)
for i, f in enumerate(system.getForces()):
    f.setForceGroup(i)
for i in range(system.getNumParticles()):
    force.addParticle(i, [])


# Here we use a more accurate integrator, where we still
# specify the temperature, friction coefficient, and timestep size
integrator = NoseHooverIntegrator(10*kelvin, 1/picosecond, 0.001*picoseconds)
simulation = Simulation(psf.topology, system, integrator)
simulation.context.setPositions(pdb.positions)
simulation.minimizeEnergy(10000)

simulation.reporters.append(PDBReporter('output.pdb', 1000))
simulation.reporters.append(StateDataReporter(stdout, 1000, step=True,
        potentialEnergy=True, temperature=True))

simulation.reporters.append(ForceReporter('file.txt',1000))
for i, f in enumerate(system.getForces()):
    state = simulation.context.getState(getEnergy=True, groups={i})
    print(f.getName(), state.getPotentialEnergy())

simulation.step(1000000)


