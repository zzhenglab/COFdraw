"""Independent ideal-gas fixture, run in a child process to isolate native code."""
import json
from pathlib import Path
import sys

import raspalib as r

reports=[]
for fractions in ([1.], [.2,.8]):
    pseudo=[r.PseudoAtom(name=f"Ideal{i}",framework_type=False,mass=40.,charge=0.,atomic_number=18) for i in range(len(fractions))]
    force=r.ForceField(pseudo_atoms=pseudo,parameters=[r.VDWParameters(0.,3.) for _ in fractions],use_charge=False)
    components=[r.Component(force_field=force,component_name=f"Ideal{i}",critical_temperature=150.,critical_pressure=4.8e6,acentric_factor=0.,
                            defined_atoms=[r.Atom((0.,0.,0.),charge=0.,type=i)],
                            particle_probabilities=r.MCMoveProbabilities(translation_probability=.5,swap_cbmc_probability=1.),
                            fugacity_coefficient=len(fractions)*fraction) for i,fraction in enumerate(fractions)]
    system=r.System(force_field=force,simulation_box=r.SimulationBox(30.,30.,30.),components=components,
                    external_temperature=300.,external_pressure=1e6,initial_number_of_molecules=[0]*len(fractions))
    mc=r.MonteCarlo(number_of_cycles=30000,number_of_initialization_cycles=2000,systems=[system],random_seed=20261117,output_to_files=False)
    mc.run()
    rows=[]
    for i,fraction in enumerate(fractions):
        mean,ci=mc.systems[0].average_loadings.average_loading_number_of_molecules(i)
        expected=fraction*1e6*(30e-10)**3/(1.380649e-23*300)
        rows.append({"physical_fraction":fraction,"mean_molecules":mean,"confidence95":ci,"ideal_gas_PV_kT":expected})
    reports.append(rows)
Path(sys.argv[1]).write_text(json.dumps(reports,indent=2))
