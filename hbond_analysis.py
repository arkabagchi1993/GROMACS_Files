import numpy as np
import MDAnalysis as mda

xtcfile = 'md100_center.xtc'
topfile = 'md100_rescale.tpr'
u = mda.Universe(topfile, xtcfile)

from MDAnalysis.analysis.hydrogenbonds.hbond_analysis import HydrogenBondAnalysis as HBA
print("Hydrogen bonds between MMP2 and Rolipram at each nanosecond")

# See https://manual.gromacs.org/current/reference-manual/analysis/hydrogen-bonds.html
gromacs_cutoff = 3.5
gromacs_cutoff_angle = 150

start_frame = 9999
end_frame = -1  # Process frames from start_frame down to 0
# Instruct script to change start frame and write
for frame in range(start_frame, end_frame, -100):
    print("Current Frame:", frame)
    
    ligand_name = 'resname LIG'
    protein_name = 'protein'

    hbonds = HBA(universe=u, \
                 d_a_cutoff=gromacs_cutoff, \
                 d_h_a_angle_cutoff=gromacs_cutoff_angle, \
                 between=[ligand_name, protein_name], \
                 update_selections=False)

    hbonds.run(start=frame, stop=frame + 1, step=1)  # Process one frame at a time
    for hbond in hbonds.results.hbonds:
        frame, donor_ix, hydrogen_ix, acceptor_ix = hbond[:4].astype(int)
        u.trajectory[frame]
        atoms = u.atoms[[donor_ix, hydrogen_ix, acceptor_ix]]
        print(atoms)

print("Processing complete")
