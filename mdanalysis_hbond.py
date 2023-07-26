%%time
import numpy as np
import MDAnalysis as mda

xtcfile='md100_center.xtc'
topfile='md100_rescale.tpr'
u = mda.Universe(topfile, xtcfile)

%%time
from MDAnalysis.analysis.hydrogenbonds.hbond_analysis import HydrogenBondAnalysis as HBA

#See https://manual.gromacs.org/current/reference-manual/analysis/hydrogen-bonds.html
gromacs_cutoff = 3.5
gromacs_cutoff_angle = 150

start_frame = 9199
end_frame = start_frame + 1

ligand_name = 'resname LIG'
protein_name = 'protein'

hbonds = HBA(universe=u,\
             d_a_cutoff=gromacs_cutoff,\
             d_h_a_angle_cutoff=gromacs_cutoff_angle,\
             between=[ligand_name, protein_name],
            update_selections=False)

hbonds.run(start=start_frame, stop=end_frame,step=1)

first_hbond = hbonds.results.hbonds[0]
frame, donor_ix, hydrogen_ix, acceptor_ix = first_hbond[:4].astype(int)
u.trajectory[frame]
atoms = u.atoms[[donor_ix, hydrogen_ix, acceptor_ix]]
print(atoms)
