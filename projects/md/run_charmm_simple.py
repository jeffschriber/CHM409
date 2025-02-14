from openmm.app import *
from openmm import *
from openmm.unit import *
from sys import stdout
import glob

## Load the psf and pdb files
psf = CharmmPsfFile('step3_pbcsetup.psf')
pdb = PDBFile('step3_pbcsetup.pdb')

# This is where the parameter files are stored
# You will need this copied in the run directory,
# with the tip216.crd file removed
pfiles = glob.glob('toppar/*')
params = CharmmParameterSet(*pfiles)

# If you want to run your simulation with periodic boundary conditions,
# you need to specify the box size
# CHARMMGUI makes the inputs suitable for PBC, you
# just need to find the dimensions when running this custom script
psf.setBox(60,60,60)

# This creates the system. See the OpenMM documentation for
# details about specific keywords
system = psf.createSystem(params, nonbondedMethod=NoCutoff,
        nonbondedCutoff=1*nanometer, constraints=HBonds)

#This specifies the integrator, which also sets the temperature,
# friction coefficient, and timestep of the simulation
# This will run an NVT simulation
integrator = LangevinMiddleIntegrator(300*kelvin, 1/picosecond, 0.004*picoseconds)

simulation = Simulation(psf.topology, system, integrator)
simulation.context.setPositions(pdb.positions)

# Here, we are running an energy minimization to prep our structure
# This will get rid of any large forces that could
# spoil the simulation 
simulation.minimizeEnergy(maxIterations=100)

# these are reporters, that will append the pdb structure
# to a file and print statistics at desired intervals
# Let's store the interval as a variable

interval = 1 
simulation.reporters.append(PDBReporter('output.pdb', interval))
simulation.reporters.append(StateDataReporter(stdout, interval, step=True,
        potentialEnergy=True, temperature=True, speed=True))

# Lastly, we specify the number of steps to take, and then take them
# This is how you control the total simulation time
nsteps = 2
simulation.step(nsteps)


