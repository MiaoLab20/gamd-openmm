from openmm.app import *
from openmm import *
from openmm.unit import *
from sys import stdout
import time

cuda_device_index = "0"
num_steps=10000
timestep = 0.002*picoseconds

pdb = PDBFile("rna_fixed.pdb")
forcefield = ForceField("amoeba2018.xml")
system = forcefield.createSystem(pdb.topology, nonbondedMethod=PME, nonbondedCutoff=1*nanometer, constraints=HBonds)
integrator = LangevinIntegrator(300*kelvin, 1/picosecond, timestep)
platform = Platform.getPlatformByName("CUDA")
properties = {"CudaDeviceIndex": cuda_device_index, "CudaPrecision": "mixed"}
simulation = Simulation(pdb.topology, system, integrator, platform, properties)
simulation.context.setPositions(pdb.positions)
simulation.minimizeEnergy()
simulation.reporters.append(PDBReporter("output.pdb", 1000))
simulation.reporters.append(StateDataReporter(stdout, 1000, step=True, potentialEnergy=True, temperature=True))

start_time = time.time()
simulation.step(num_steps)
total_time = time.time() - start_time
simulation_in_ns = num_steps * timestep.value_in_unit(picoseconds) * 1e-3
total_time_in_days = total_time / (86400.0)
ns_per_day = simulation_in_ns / total_time_in_days
print("Equilibration benchmark:", ns_per_day, "ns/day")
