#!/usr/bin/python3

from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from sys import stdout
from sys import exit
import os
import traceback

from gamd_langevin import *


def createPositionsFile(integrator, thefile):
    positions = integrator.get_coordinates()
    f = open(thefile, 'w')
    f.write("particle, x, y, z\n")
    for i in range(len(positions)):
       f.write(str(i) + ", " + str(positions[i][0]) + ", " + str(positions[i][1])+ ", " + str(positions[i][2]) +"\n" )
    f.close()

def createGamdLog(gamdLog, filename):
    with  open(filename, 'w') as f:
        keys = gamdLog[0].keys()
        for header in keys[:-1]:
            f.write(header + ", ")
        f.write(keys[-1] + "\n")
        for entry in gamdLog:
            for header in keys[:-1]:
                f.write(str(entry[header]) + ", ")
            f.write(str(entry[keys[-1]]) + "\n")


def createGamdSimulationFromAmberFiles(prmtopfile, inpcrdfile, lowerBound=True, mode='TotalBoost'):
    prmtop = AmberPrmtopFile(prmtopfile)
    inpcrd = AmberInpcrdFile(inpcrdfile)
    system = prmtop.createSystem(nonbondedMethod=PME, nonbondedCutoff=1*nanometer, constraints=HBonds)
    if mode == 'DualBoost':
        # Move torsion terms to force group 1
        group = 1
        for force in system.getForces():
            if force.__class__.__name__ != 'PeriodicTorsionForce':
                force.setForceGroup(group)
                break
        
        if lowerBound:
            print("group:", group)
            integrator = GamdDualBoostPotentialLangevinIntegratorLowerBound(2.0 * femtoseconds, number_of_steps_in_stage_1,
                                                                     number_of_steps_in_stage_2, ntave,
                                                                     6.0 * kilocalories_per_mole, 6.0 * kilocalories_per_mole,
                                                                     temperature=300*kelvin, gamma=1/picosecond, group=group)
        else:
            integrator = GamdDualBoostPotentialLangevinIntegratorUpperBound(2.0 * femtoseconds, number_of_steps_in_stage_1,
                                                                     number_of_steps_in_stage_2, ntave,
                                                                     6.0 * kilocalories_per_mole, 6.0 * kilocalories_per_mole,
                                                                     temperature=300*kelvin, gamma=1/picosecond, group=group)
                                                                     
    elif mode == 'DihedralBoost':
        # Move torsion terms to force group 1
        group = 1
        for force in system.getForces():
            if force.__class__.__name__ != 'PeriodicTorsionForce':
                force.setForceGroup(group)
                break
        
        if lowerBound:
            print("group:", group)
            integrator = GamdGroupBoostPotentialLangevinIntegratorLowerBound(2.0 * femtoseconds, number_of_steps_in_stage_1,
                                                                     number_of_steps_in_stage_2, ntave,
                                                                     6.0 * kilocalories_per_mole,
                                                                     temperature=300*kelvin, gamma=1/picosecond, group=group)
        else:
            integrator = GamdGroupBoostPotentialLangevinIntegratorUpperBound(2.0 * femtoseconds, number_of_steps_in_stage_1,
                                                                     number_of_steps_in_stage_2, ntave,
                                                                     6.0 * kilocalories_per_mole,
                                                                     temperature=300*kelvin, gamma=1/picosecond, group=group)
    elif mode == 'TotalBoost':
        if lowerBound:
            integrator = GamdTotalBoostPotentialLangevinIntegratorLowerBound(2.0 * femtoseconds, number_of_steps_in_stage_1,
                                                                     number_of_steps_in_stage_2, ntave,
                                                                     6.0 * kilocalories_per_mole,
                                                                     temperature=300*kelvin, gamma=1/picosecond)
        else:
            integrator = GamdTotalBoostPotentialLangevinIntegratorUpperBound(2.0 * femtoseconds, number_of_steps_in_stage_1,
                                                                     number_of_steps_in_stage_2, ntave,
                                                                     6.0 * kilocalories_per_mole,
                                                                     temperature=300*kelvin, gamma=1/picosecond)
    else:
        raise Exception, "'mode' option not allowed: %s" % mode
        
    print("Global Variables:")
    for i in range(integrator.getNumGlobalVariables()):
        print(i, integrator.getGlobalVariableName(i))
    
    print("DOF Variables:")
    for i in range(integrator.getNumPerDofVariables()):
        print(i, integrator.getPerDofVariableName(i))
        
    print("Computations:")
    for i in range(integrator.getNumComputations()):
        print(i, integrator.getComputationStep(i))
    
    #exit()
    
    simulation = Simulation(prmtop.topology, system, integrator)
    simulation.context.setPositions(inpcrd.positions)
    if inpcrd.boxVectors is not None:
        simulation.context.setPeriodicBoxVectors(*inpcrd.boxVectors)
    simulation.minimizeEnergy()

    return [simulation, integrator]


def createGamdSimulationFromPdbFile(pdbfile, prmtopfile, lowerBound=True, mode='TotalBoost'):
    pdb = PDBFile(pdbfile)
    prmtop = AmberPrmtopFile(prmtopfile)
    system = prmtop.createSystem(nonbondedMethod=PME, nonbondedCutoff=1 * nanometer, constraints=HBonds)
    if mode == 'DualBoost':
        # Move torsion terms to force group 1
        group = 1
        for force in system.getForces():
            if force.__class__.__name__ != 'PeriodicTorsionForce':
                force.setForceGroup(group)
                break
        
        if lowerBound:
            print("group:", group)
            integrator = GamdDualBoostPotentialLangevinIntegratorLowerBound(2.0 * femtoseconds, number_of_steps_in_stage_1,
                                                                     number_of_steps_in_stage_2, ntave,
                                                                     6.0 * kilocalories_per_mole, 6.0 * kilocalories_per_mole,
                                                                     temperature=300*kelvin, gamma=1/picosecond, group=group)
        else:
            integrator = GamdDualBoostPotentialLangevinIntegratorUpperBound(2.0 * femtoseconds, number_of_steps_in_stage_1,
                                                                     number_of_steps_in_stage_2, ntave,
                                                                     6.0 * kilocalories_per_mole, 6.0 * kilocalories_per_mole,
                                                                     temperature=300*kelvin, gamma=1/picosecond, group=group)
                                                                     
    elif mode == 'DihedralBoost':
        # Move torsion terms to force group 1
        group = 1
        for force in system.getForces():
            if force.__class__.__name__ != 'PeriodicTorsionForce':
                force.setForceGroup(group)
                break
        
        if lowerBound:
            print("group:", group)
            integrator = GamdGroupBoostPotentialLangevinIntegratorLowerBound(2.0 * femtoseconds, number_of_steps_in_stage_1,
                                                                     number_of_steps_in_stage_2, ntave,
                                                                     6.0 * kilocalories_per_mole,
                                                                     temperature=300*kelvin, gamma=1/picosecond, group=group)
        else:
            integrator = GamdGroupBoostPotentialLangevinIntegratorUpperBound(2.0 * femtoseconds, number_of_steps_in_stage_1,
                                                                     number_of_steps_in_stage_2, ntave,
                                                                     6.0 * kilocalories_per_mole,
                                                                     temperature=300*kelvin, gamma=1/picosecond, group=group)
    elif mode == 'TotalBoost':
        if lowerBound:
            integrator = GamdTotalBoostPotentialLangevinIntegratorLowerBound(2.0 * femtoseconds, number_of_steps_in_stage_1,
                                                                     number_of_steps_in_stage_2, ntave,
                                                                     6.0 * kilocalories_per_mole,
                                                                     temperature=300*kelvin, gamma=1/picosecond)
        else:
            integrator = GamdTotalBoostPotentialLangevinIntegratorUpperBound(2.0 * femtoseconds, number_of_steps_in_stage_1,
                                                                     number_of_steps_in_stage_2, ntave,
                                                                     6.0 * kilocalories_per_mole,
                                                                     temperature=300*kelvin, gamma=1/picosecond)
    else:
        raise Exception, "'mode' option not allowed: %s" % mode
        
    simulation = Simulation(prmtop.topology, system, integrator)
    simulation.context.setPositions(pdb.positions)

    simulation.minimizeEnergy()

    return [simulation, integrator]


def create_output_directories(directories):
    for dir in directories:
        os.makedirs(dir, 0o755)


def upperBoundFunctionWrapper(a_function):
    def return_function(prmtopfile, coordinatesfile):
        a_function(prmtopfile, coordinatesfile, False)
    return return_function
  
def lowerBoundGroupFunctionWrapper(a_function):
    def return_function(prmtopfile, coordinatesfile):
        a_function(prmtopfile, coordinatesfile, False, True)
    return return_function
  
def upperBoundGroupFunctionWrapper(a_function):
    def return_function(prmtopfile, coordinatesfile):
        a_function(prmtopfile, coordinatesfile, False, True)
    return return_function


def main():

    output_directory = "output"
    coordinates_file = './data/md-4ns.rst7'
    prmtop_file =  './data/dip.top'

    total_number_of_steps = number_of_steps_in_stage_1 + number_of_steps_in_stage_2 + number_of_steps_in_stage_3

    function_dictionary = {'total-lb': {'amber' : createGamdSimulationFromAmberFiles,
                                        'pdb': createGamdSimulationFromPdbFile},
                            'total-ub': {'amber' : upperBoundFunctionWrapper(createGamdSimulationFromAmberFiles),
                                         'pdb': upperBoundFunctionWrapper(createGamdSimulationFromPdbFile) },
                           'group-lb': {'amber' : lowerBoundGroupFunctionWrapper(createGamdSimulationFromAmberFiles),
                                        'pdb': lowerBoundGroupFunctionWrapper(createGamdSimulationFromPdbFile)},
                            'group-ub': {'amber' : upperBoundGroupFunctionWrapper(createGamdSimulationFromAmberFiles),
                                         'pdb': upperBoundGroupFunctionWrapper(createGamdSimulationFromPdbFile) }}


    create_output_directories([output_directory, output_directory + "/states/", output_directory + "/positions/",
                               output_directory + "/pdb/", output_directory + "/checkpoints"])

    (simulation, integrator) = createGamdSimulationFromAmberFiles(prmtop_file, coordinates_file,lowerBound=lowerBound, mode=mode)


    simulation.saveState(output_directory + "/states/initial-state.xml")
    simulation.reporters.append(PDBReporter(output_directory + '/output.pdb', 100))
    simulation.reporters.append(StateDataReporter(stdout, 100, step=True, temperature=True,
                                                  potentialEnergy=True, totalEnergy=True, volume=True))
    gamdLog = []
    simulation.context.setVelocitiesToTemperature(300.0 * kelvin)
    for step in range(total_number_of_steps):
        simulation.step(1)
        
        if step % ntave == 0:
            print("Vmax:", integrator.getGlobalVariableByName("Vmax"), "Vmin:", integrator.getGlobalVariableByName("Vmin"))
            print("Vavg:", integrator.getGlobalVariableByName("Vavg"), "boostPotential:", integrator.getGlobalVariableByName("boostPotential"))
            if mode == 'DihedralBoost':
                print("groupForceScalingFactor", integrator.getGlobalVariableByName("groupForceScalingFactor"))
                print("E:", integrator.getGlobalVariableByName("E"), "k0:", integrator.getGlobalVariableByName("k0"))
                print("groupEnergy:", integrator.getGlobalVariableByName("groupEnergy"), "totalEnergy:", integrator.getGlobalVariableByName("totalEnergy"))
            elif mode == 'TotalBoost':
                print("totalForceScalingFactor", integrator.getGlobalVariableByName("totalForceScalingFactor"))
                print("E:", integrator.getGlobalVariableByName("E"), "k0:", integrator.getGlobalVariableByName("k0"))
                print("currentEnergy:", integrator.getGlobalVariableByName("currentEnergy"))
            elif mode == 'DualBoost':
                print("totalForceScalingFactor", integrator.getGlobalVariableByName("totalForceScalingFactor"))
                print("E:", integrator.getGlobalVariableByName("E"), "k0:", integrator.getGlobalVariableByName("k0"))
                print("groupForceScalingFactor", integrator.getGlobalVariableByName("groupForceScalingFactor"))
                print("EGroup:", integrator.getGlobalVariableByName("EGroup"), "k0Group:", integrator.getGlobalVariableByName("k0Group"))
                print("groupEnergy:", integrator.getGlobalVariableByName("groupEnergy"), "totalEnergy:", integrator.getGlobalVariableByName("totalEnergy"))
                
            simulation.saveState(output_directory + "/states/" + str(step) + ".xml")
            simulation.saveCheckpoint(output_directory + "/checkpoints/" + str(step) + ".bin")
            if mode == 'DihedralBoost':
                gamdLog.append({'total_nstep': step,
                            'Unboosted-Potential-Energy': integrator.get_current_potential_energy(),
                            'Group-Force-Weight': integrator.get_group_force_scaling_factor(),
                            'Boost-Energy-Potential': integrator.get_boost_potential() })
            elif mode == 'TotalBoost':
                gamdLog.append({'total_nstep': step,
                            'Unboosted-Potential-Energy': integrator.get_current_potential_energy(),
                            'Total-Force-Weight': integrator.get_total_force_scaling_factor(),
                            'Boost-Energy-Potential': integrator.get_boost_potential() })
            elif mode == 'DualBoost':
                gamdLog.append({'total_nstep': step,
                            'Unboosted-Potential-Energy': integrator.get_current_potential_energy(),
                            'Total-Force-Weight': integrator.get_total_force_scaling_factor(),
                            'Group-Force-Weight': integrator.get_group_force_scaling_factor(),
                            'Boost-Energy-Potential': integrator.get_boost_potential() })
                
            
            createPositionsFile(integrator, output_directory + '/positions/coordinates-' + str(step) + '.csv')

    createGamdLog(gamdLog, output_directory + "/gamd.log")

mode = "TotalBoost"
#mode = "DihedralBoost"
#mode = "DualBoost"
lowerBound = True
ntave = 1000 # 1000
number_of_steps_in_stage_1 = 10000 #10000
number_of_steps_in_stage_2 = 100000 # 100000
number_of_steps_in_stage_3 = 100000 # 100000

if __name__ == "__main__":
    main()
