**###GROMACS UNIVERSAL TUTORIAL###**

# Basics of KUHPC #
1. `ssh <id>@182.16.16.10` (from your Terminal (Linux) or Powershell (Windows)
2. create a directory (which will be your working directory) by `mkdir <directory name>`
3. enter the directory by `cd <directory name>`
4. Upload essential files like the .pbs files (pbs scripts), .mdp files (provided in the gromacs tutorial) and .pdb files (of the receptor and ligand) by WinSCP (for windows) or SFTP (Linux).
5. You can list the files in your directory by typinf `ls`.
6. To remove any file type `rm <filename>`.
7. To copy files type `cp /path/to/source/file /path/to/destination/file`
8. To rename or move files type `mv /path/tp/source/file /path/to/destination/file`



# REQUIRED-SOFTWARES-FOR-GROMACS #

#########################################################################
-	For adding hydrogens to the ligand and export .mol2 file 	-
-	Install Avogadro						-

######################################################################### 
-	For visualization of .pdb and .gro files			-
-	Install PyMol							-

#########################################################################
-	For final MD trajectory visualization				-	
-	Install VMD (Requires high spec computer for proper visualization

#########################################################################
-	For xvg file visualization					-
-	Install grace (for Linux)					-

#########################################################################
-	Use Matplot (from jupyter notebook) (for Windows/Linux)		-
-	(Detailed instructions given below)				-

#########################################################################
-	If you are simulating a protein with side chains, sometimes fixing of side-chains is required, otherwise `gmx pdb2gmx` gives an error
-	For that install spdbv software (Swiss PDB Viewer Software)
-	There you can fix the side chains from 'tools'
-	This software is available for windows and mac and not for linux.

#########################################################################


#########################################################################
# STEPS FOR MD AS FOLLOWS #
#########################################################################

# 1. PREPARE LIGAND IN AVOGADRO #
1.1. Open the ligand pdb file and go to "Build" --> "Hydrogens" --> "Add Hydrogens" (This will add hydrogen molecules to the ligand, which might appear as new white dots)
1.2. Then go to "Files" --> "Export" --> "Molecule". Select the "file type" as "All Files" and save the molecule as "LIG.mol2". (The name "LIG" can be changed as per your wish, but then subsequent changes in the next commands are necesary. Whenever there is something written as "LIG", that must be changed with the new name given by you. Renaming the files will not have any impact on the MD run.


# 2. Correction to be made in LIG.mol2 #
###Open LIG.mol2 by using gedit command or simply opening file in any text editor###
2.1.	"@<TRIPOS>MOLECULE" make sure this is the first line in file
	delete the header and empty space if you have to	
2.2.      "@<TRIPOS>MOLECULE” there will be name after this line maybe xxx.pdb or ****** or anything else
	change it to LIG (Or the name given by you instead of "LIG") 
	
2.3.	bond orders "@<TRIPOS>BOND" will be arranged differently in each file 
arrange them in specific order to avoid errors use "sort_mol2_bonds.pl" script.
	
	perl sort_mol2_bonds.pl LIG.mol2 LIG.mol2


3. Go to SwissParam "http://www.swissparam.ch/" and upload the 'LIG.mol2 file'

4. Download the .zip folder
(Sometimes you will not be able to download the .zip file, then follow these instructions:
	
i. In your working directory, type `wget <URL of the Swissparam result page>` (This will download a file named "index.html".

ii. Open this index.html file with "nano" or any other text editor of your choice. To open the file with nano, type `nano index.html`.

iii. In this file, you will find a line, `Results can be found in the following zip file: <a class="sib_link"href="http://www.swissparam.ch/results/514710792/<filename>.zip"><filename>.zip</a>. <br><br>`.

(Where 'filename' will be the name of your .mol2 file). From this line copy the `http://www.swissparam.ch/results/514710792/<filename>.zip` part.

iv. Now exit the text editor and in your working directory, type `wget <the link you just copied>` (Eg: wget http://www.swissparam.ch/results/514710792/LIG.zip)

v. This will now download a .zip file

vi. To unzip the .zip file, type 

	unzip <filename>.zip

vii. Make sure now you have the LIG.pdb, LIG.itp in your working directory.



# THINGS TO BE NOTED #


#########################################################################
-	In KUHPC, Gromacs runs as MPI, so for each command to work, 	-
-	you have to initiate the `gmx` command as `gmx_mpi`.		-
#########################################################################


# Using Sifimage of GROMACS from Nvidia NGC catalog #

1. Download the SIFIMAGE and place it in 
	~/.config/sifdir/gromacs_2022.3.sif

2. Now you can use the following command--


		LD_LIBRARY_PATH="" singularity run --nv -B ${PWD}:/workspace ~/.config/sifdir/gromacs_2022.3.sif bash -c "cd /workspace && gmx hbond -s md100_rescale.tpr -f md100_center.xtc -num hb.xvg -tu ns"


**(N.B.: The inverted commas are part of the command in this case, do not remove it)**

3. Here you have to replace the `gmx hbond -s md100_rescale.tpr -f md100_center.xtc -num hb.xvg -tu ns` with your command of interest, such as `gmx pdb2gmx -f REC.pdb -o REC.gro`.



###########################################################################

# Begin the MD simulation #

###########################################################################

	gmx_mpi pdb2gmx -f REC.pdb -ignh

Choose the options:

	8 (CHARMM27)


	1 (TIP3P)

*Point to note*: If you are simulating a protein with multiple side chains, sometimes fixing the side-chains of the PDB acquired from "www.rcsb.org" is required. For that use the "spdbv" software as described below.

**This will generate a file named "conf.gro".** You can also specify the name of the output gro file by `-o name.gro` option.

	gmx_mpi editconf -f LIG.pdb -o LIG.gro

# If you are simulating Protein-Protein or Protein-DNA/RNA then #

Generate the LIgand topology with
	
 	gmx_mpi pdb2gmx -f LIG.pdb -o LIG.gro -ignh

This will generate three files, LIG.gro, posre.itp and topol.top. Now rename the posre.itp as "posre_LIG.itp" and topol.top as "LIG.top"
# Edit the LIG.top file as described below: #
Remove from the top of the file,
	
 	; Include forcefield parameters
 	#include amberGS.f/forcefiled.itp"
And remove from the bottom of the file the protein mentioned below, this will ensure that no itp file is called twice, as they are already called at the topol.top file.
	
 	; Include Position restraint file
	#ifdef POSRES
	#include "posre.itp"
	#endif

	; Include water topology
	#include "amber96.ff/tip3p.itp"

	#ifdef POSRES_WATER
	; Position restraint for each water oxygen
	[ position_restraints ]
	;  i funct	 fcx        fcy        fcz
   	1    1	1000	   1000       1000
	#endif

	; Include topology for ions
	#include "amber96.ff/ions.itp"

After that you have to also remove the part from the `LIG.top` file, described below: (Which will ensure that the RNA/DNA/Protein molecule is called only once in the `topol.top` file)
For Protein/RNA simulation:

	[molecules]
 	RNA_chain_A    1

For Protein-DNA simulation:

	[molecules]
 	DNA_chain_A    1

For Protein-Protein simulation

	[molecules]
 	Protein_chain_A	   1

# Build the .gro file for the complex #
**Now you have to copy the "conf.gro" and "LIG.gro" file into a "complex.gro" file**
To do that, follow the steps-----
	1. Type `cp ./conf.gro ./complex.gro`  (This will make a copy the "conf.gro" and name it as "complex.gro")
	2. Now open LIG.gro file in text editor by typing `nano LIG.gro` and copy from the third line of the file to the second last line of the file.
	3. To open the "complex.gro" file in text editor, type `nano complex.gro`
	4. Now go to the last line of the file, you will find some numbers (cordinates) written there. Place the cursor there and paste the copied part from the LIG.gro file. OR you can just place the cursor at the described place and type `Ctrl+R` and the `LIG.gro` and press enter. This will insert the whole `LIG.gro` file there. Then you just have to remove the unnecessary lines.
	5. Now, you will find a `molecule number` at the top of each conf.gro and LIG.gro files (mentioned as just numbers, such as 4265 in case of conf.gro and 60 in case of LIG.gro). Add these to numbers (which will be for eg. 4285) and simply replace the molecule number of "complex.gro" file with the added value.
	
	
*(You can check the receptor and ligand by downloading the "complex.gro" and opening it in PyMol)*

# EDIT THE FOLLOWING in topol.top #
`nano topol.top`
add 

	; Include ligand topology 
	#include "LIG.itp"

	below- Include forcefield parameters
	#include "amberGS.ff/forcefield.itp")

# If you are performing the simulation with Protein-Protein or Protein-DNA/RNA #

Then add,
	
 	; Include ligand topology
 	#include "LIG.top"

below
	
 	; Include forcefield parameters
 	#include "amberGS.ff/forcefield.itp"


AT THE BOTTOM OF THE SAME FILE PERFORM FOLLOWING CHANGES
(add `LIG 1`
align exactly below-
	
 	Protein_chain_A     1)

SO, IT WILL LOOK LIKE--
	
	Protein_chain_A		1
	LIG			1

For Protein-DNA/RNA simulation `LIG` will be replaced by `DNA_chain_A` and `RNA_chain_A` respectively.

# EDIT THE FOLLOWING in lig.itp  or LIG.top #

	[ moleculetype ]
	; Name nrexcl
	lig_gmx2 3
TO

	[ moleculetype ]
	; Name nrexcl
	LIG 3
*(in certain cases this will already be `LIG 3` so for such case no change is needed)*

----------
Then to place the system or complex within a box, use the following `gmx editconf` command:

	gmx_mpi editconf -f complex.gro -d 1.0 -bt triclinic -o box.gro 

*(You can also change `triclinic` to `dodecahedron` as per your requirement)*
Next step is the solvation of the complex..

Now add solvent into the box created above by the `gmx solvate` command:

	gmx_mpi solvate -cp box.gro -cs spc216.gro -p topol.top -o solv.gro

# Build the .tpr file for the ions, that will be used for total charge neutralization #

Use the [ions.mdp](MDP_Files/ions.mdp) file for the command.
	
	gmx_mpi grompp -f ions.mdp -c solv.gro -p topol.top -o ions.tpr

(OR)
	
	gmx_mpi grompp -f ions.mdp -c solv.gro -maxwarn 2 -p topol.top -o ions.tpr

*(The `-maxwarn 2` option is sometimes required to ignore the warnings)*

Then, use the following `gmx genion` command to add required number of NA and CL to neutralize the charge of the whole system. You can also specify the number of atoms you want to add by seeing the `qtot` number described in the `[ATOMS]` section of the topol.top file. Otherwise `gmx genion` will automatically calculate the total charge of the system and add required NA and CL into the system.

	gmx_mpi genion -s ions.tpr -p topol.top -conc 0.1 -neutral -o solv_ions.gro

(Select Option) This will ensure that `gmx genion` replaces "SOL" to add NA and CL

	15 (For "SOL")

# Energy Minimization #
Then, go for energy minimization. To build the energy minimization tpr file (em.tpr) by `gmx grompp`, use the [em.mdp](MDP_Files/em.mdp). And use the following command:
	
	gmx_mpi grompp -f em.mdp -c solv_ions.gro -p topol.top -o em.tpr

(OR)
	
	gmx_mpi grompp -f em.mdp -c solv_ions.gro -maxwarn 2 -p topol.top -o em.tpr

**(The `-maxwarn 2` option is sometimes required to ignore the warnings)**
Then for final energy minimization

	gmx_mpi mdrun -v -deffnm em

**(For this you can use the [gromacs_em.pbs](PBS_Files/gromacs_em.pbs) script)**
In that script you can change the walltime, output name, job name as per your requirement.

# Making index file for ligand #
Now make index files.
 
	gmx_mpi make_ndx -f LIG.gro -o index_LIG.ndx

(Select options)
	
	> 0 & ! a H*
 	> q

# Making the position restraint file for the ligand #

Now, make the position restraint file for the ligand. It is not always required for Protein-Protein and Protein-DNA/RNA simulation as the `gmx pdb2gmx` already builds a position restraint file for the ligand (which was previously described to rename as "posre_LIG.itp".
	
	gmx_mpi genrestr -f LIG.gro -n index_LIG.ndx -o posre_LIG.itp -fc 1000 1000 1000

(Select group)

	
	> select group "3"

	
Now, open topol.top file
at the end of the document 


	; Include Position restraint file
	#ifdef POSRES
	#include "posre.itp"
	#endif
	
Modify it as
	
	; Include Position restraint file
	#ifdef POSRES
	#include "posre.itp"
	#include "posre_LIG.itp"
	#endif

**For Protein-Protein and Protein-DNA/RNA simulation, the posre_LIG.itp file needs to be included at a different place, which is at the top of the file where,**

	; Include ligand topology
 	#include "LIG.top"

is mentioned. Modify it as,

	; Include ligand topology
 	#include "LIG.top
  	#ifdef POSRES
   	#include "posre_LIG.itp"
    	#endif


# Making other Index file for the Complex system from the whole System #
This index file is specifically useful for Protein-ligand or Protein-DNA/RNA/Protein complex simulation, as the new index option created as `Protein_LIG` or `Protein_DNA` or `Protein_RNA` or `Protein_Protein` can be used in post-processing steps to specify the `Protein_LIG` or `Protein_DNA` or `Protein_RNA` or `Protein_Protein` complexes rather that specifying the whole system involving SOL and other molecules.

	gmx_mpi make_ndx -f em.gro -o index.ndx

(Select Options) Which will signify Protein and LIG

	
	> 1 | 13
	> q

**The option might change to the following for Protein-DNA/RNA simulation** Which will signify Protein and DNA/RNA

	> 1 | 12
 	> q
# NVT MINIMIZATION #
**Remember to edit the [nvt.mdp](MDP_Files/nvt.mdp) file to insert proper tc coupling groups.**
For Protein-Ligand simulation, choose tc groups as 

	Protein_LIG Water_and_ions
For Protein-DNA simulation, choose tc groups as

	Protein_DNA Water_and_ions
For Protein-RNA simulation, choose rc groups as 

	Protein_RNA Water_and_ions
For Protein in water (Protein only) simulation, choose tc groups as

	Protein Non-protein

Then use the following `gmx grompp` command to generate `nvt.tpr` file


	gmx_mpi grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -n index.ndx -maxwarn 2 -o nvt.tpr
	
AND THEN
	
	gmx_mpi mdrun -deffnm nvt
**(For this you can use [gromacs_nvt.pbs](PBS_Files/gromacs_nvt.pbs) script)**


# NPT MINIMIZATION #
Remember to edit the [npt.mdp](MDP_Files/npt.mdp) file to insert proper tc coupling groups.
For Protein-Ligand simulation, choose tc groups as 

	Protein_LIG Water_and_ions
For Protein-DNA simulation, choose tc groups as

	Protein_DNA Water_and_ions
For Protein-RNA simulation, choose rc groups as 

	Protein_RNA Water_and_ions
For Protein in water (Protein only) simulation, choose tc groups as

	Protein Non-protein

Then use the following `gmx grompp` command to generate `npt.tpr` file

	gmx_mpi grompp -f npt.mdp -c nvt.gro -r nvt.gro -p topol.top -n index.ndx -maxwarn 2 -o npt.tpr
	
AND THEN
	
	gmx_mpi mdrun -deffnm npt
**(For this you have to use [gromacs_npt.pbs](PBS_Files/gromacs_npt.pbs) script)**


# FEW THINGS TO KEEP IN MIND WHILE USING THE `npt.mdp` and `md.mdp` FILE #
1. Pressure coupling can be changed to `anisotropic` if required (See references).
2. `Refcoord scaling` can be removed if required.
3. There are several other pressure coupling groups that can be applied with gromacs. You can utilize them as per your requirement.


# FINAL MD RUN/PRODUCTION #
`nano md.mdp` (Change MD RUN TIME as per your need).
-	Check for all the parameters in the [md.mdp](MDP_Files/md.mdp) file to match the previously used [nvt.mdp](MDP_Files/nvt.mdp) and [npt.mdp](MDP_Files/npt.mdp) files.		-



		gmx_mpi grompp -f md.mdp -c npt.gro -t npt.cpt -p topol.top -n index.ndx -maxwarn 2 -o md.tpr

AND THEN
	
	gmx_mpi mdrun -deffnm md
**(For this you can use [md_complex.pbs](PBS_Files/md_complex.pbs) script or the [md_cpu.pbs](PBS_Files/md_cpu.pbs) script. Remember to change the MD_NAME, Job name, walltime, output filename accordingly)**



# RESUME MD RUN FROM CHECKPOINT #

**If the run somehow stops, it will generate a checkpoint file (`MD_NAME.cpt`). You can use this file to run the [resume_md.pbs](PBS_Files/resume_md.pbs) script, which will resume the run from the point where it stopped**

# EXTEND MD RUN FROM LAST CHECKPOINT #

If you want to extend the run any further, follow the steps:

Build the tpr file for extended time
	
	
	gmx_mpi convert-tpr -s md_10.tpr -extend 10000 -o md_20.tpr

Change the value after `-extend` accordingly (the value indicates the amount of time increased in picoseconds)

Then run the md-run from the last checkpoint file created:
	
	gmx_mpi mdrun -v -deffnm md_20 -cpi md_10.cpt -noappend

You can also find a checkpoint file named `md_prev.cpt`.


**You can use the [extend_md.pbs](PBS_Files/extend_md.pbs) for this extension.**

#################################################################################
#										#
# THIS IS THE END OF FINAL MD RUN #						
#										#
#################################################################################

----------------------------------------------------------------------------------------

----------------------------------------------------------------------------------------
########################################################################################



#################################################################################
#										#
# POST-PROCESSING #								
#										#
#################################################################################



# Recentering and Rewrapping Coordinates #

	gmx_mpi trjconv -s md.tpr -f md.xtc -o md_center.xtc -center -pbc mol -ur compact
#Choose `Protein` for centering and `System` for output.

# Point to Note if you are extending a simulation #


As described earlier, extending a MD simulation involves extending the md.tpr file with `gmx convert-tpr` option. After extending the run with `-noappend` (using `extend_md.pbs` script), you need to concentrate the new trajectory (which in this case will be named something like "md_20.part0002.xtc" with the existing one, using `gmx trjcat`. Use the following command for that:

	gmx_mpi trjcat -f md_10.xtc md_20.part0002.xtc -o md_20.xtc

After that you can proceed for centering with the previously mentioned command.


# Dumping pdb at different time frames #
To extract the first frame (t = 0 ns) of the trajectory, use `trjconv` and  `-dump` with the recentered trajectory:

	gmx_mpi trjconv -s md.tpr -f md_center.xtc -o start.pdb -dump 0
(Here "0" refers to 0 picoseconds)
#To extract any time point frame (such as t = 10 ns) of the trajectory, use `trjconv` and `-dump` with the recentered trajectory:

	gmx_mpi trjconv -s md.tpr -f md_center.xtc -o start.pdb -dump 10000
(Here "10000" means 10000 picoseconds = 10 ns)



# RMSD Calculations #

	gmx_mpi rms -s md.tpr -f md_center.xtc -o rmsd.xvg

(OR)

	gmx_mpi rms -s md.tpr -f md_center.xtc -o rmsd.xvg -tu ns 

(Select Options respectively) (Backbone and Backbone)

	4

	4
**(Here `-tu ns` option ensures time unit to be ns)**

# RMSD Calculation for Protein-Protein Simulation #

For Protein-Protein simulation, it is important to split the two chains before calculating the RMSDs of the two proteins. So, for that you have to make new index for the splitted chains, with the following command:

	gmx_mpi make_ndx -f md.gro -o chain_CA.ndx

Type the following option

 	> splitch 3

and then

	>q

Then utilize the following command to calculate RMSD

	gmx_mpi rms -s md.tpr -f md_center.xtc -o rmsd.xvg -tu ns

Choose the following option

	> 19 (C-alpha_chain1)

AND then,

	> 19 (C-alpha_chain1)

For the next chain, utilize the same `gmx rms` command and choose

	> 20 (C-alpha_chain2)

AND then,

	>20 (C-alpha_chain2)




# RMSF Calculations #

	gmx_mpi rmsf -s md.tpr -f md_center.xtc -o rmsf.xvg

(Select Option) (Backbone)

	4

**(OR)**
For residue specific RMSF calculation, use:

	gmx_mpi rmsf -s md.tpr -f md_center.xtc -o rmsf.xvg -res

# RMSF Calculation for Protein-Protein Simulation #

For calculation of RMSF of the two protein chains, first make index for the full chains

	gmx_mpi make_ndx -f md.gro -o chain.ndx

Type the following option

	> splitch 1

AND then,

	> q

Then utilize the following command to calculate the RMSF. It is important to note that for this step, do not use `-res` option, as the residue number of two chains might overlap. So, it is better to calculate RMSF atom wise.

	gmx_mpi rmsf -s md.tpr -f md_center.xtc -o rmsf.xvg

Choose option

	> 19 (Protein_chain1)

AND for the next chain,

	> 20 (Protein_chain2)


# Calculating No.of hydrogen bonds #

	gmx_mpi hbond -s md.tpr -f md_center.xtc -num hb.xvg

(OR)

	gmx_mpi hbond -s md.tpr -f md_center.xtc -num hb.xvg -tu ns

(Select Option) (Protein and LIG) **(Options may change depending on the run, such as 4 & 12 for Protein-DNA/RNA simulation)**

	1

	13

**Note: These `gmx_mpi hbond` command does not work in kuhpc, so follow th next steps to use sifimage of latest version of GROMACS and run it on GPU.
# For running the hbond with Gromacs sifimage #
	
	ssh kuhpcgn1

(OR)

	ssh kuhpcgn2

cd to the working directory and then type:

	LD_LIBRARY_PATH="" singularity run --nv -B ${PWD}:/workspace ~/.config/sifdir/gromacs_2022.3.sif bash -c "cd /workspace && gmx hbond -s md.tpr -f md_center.xtc -num hb.xvg -tu ns"

Then select options

# Hydrogen Bonds for Protein-Protein Complex #

Use the following command:

	LD_LIBRARY_PATH="" singularity run --nv -B ${PWD}:/workspace ~/.config/sifdir/gromacs_2022.3.sif bash -c "cd /workspace && gmx hbond -s md.tpr -f md_center.xtc -n chain.ndx -num hb.xvg -tu ns"

Then choose options

	> 19 (Protein_chain1)
 	> 20 (Protein_chain2)



# Calculation of Gyration Radius #

	gmx_mpi gyrate -s md.tpr -f md_center.xtc -o gyrate1.xvg

#Choose the group of your choice

	4

# ENERGY Calculations #

	gmx_mpi energy -f md.edr -o energy1.xvg

#Choose the option of your choice


#########################################################################
# Analysis and identification of H-bond residues #
#########################################################################

-------------- With Matplot (Jupyter notebook)-----------------
1. Install jupyter notebook in a conda environment by typing, `conda install jupyter notebook` or `pip install jupyter notebook`.
2. Install `MDAnalysis` in the same conda environment.
3. Open jupyter notebook and open Python3 kernel (from 'New').
4. Paste the lines from the [mdanalysis_hbond.py](mdanalysis_hbond.py).
5. Give the proper input files in the `/path/to/input/files` section (the final trajectory file named `md_center.xtc` and the final tpr file named `md.tpr`) and change the `start=` frame according to your query.
6. Click `run`.

**OR**

If you want to analyse all the hydrogen bonds at a time, you can use the [gromacs_hbond_analysis.pbs](PBS_Files/gromacs_hbond_analysis.pbs) script. This script utilizes the [hbond_analysis.py](hbond_analysis.py) python script to write all the hydrogen bonds at each nanosecond in a text file. You can change the name of the name of the text file by changing the `hb_2R_all.txt` in the `python hbond_analysis.py >> hb_2R_all.txt` section of the [gromacs_hbond_analysis.pbs](PBS_Files/gromacs_hbond_analysis.pbs) script. You will also have to keep in mind to change the name of the python environment in the `python3` section of the `source opt/anaconda3/bin/activate ${PBS_O_HOME}/.conda/envs/python3` line in the [gromacs_hbond_analysis.pbs](PBS_Files/gromacs_hbond_analysis.pbs) script. You will also have to enter the path or name of the xtc file and tpr file in the `/path/to/xtc/file` and `/path/to/tpr/file` section of the [hbond_analysis.py](hbond_analysis.py) python script.


#########################################################################
# Visualization of xvg files #
#########################################################################



-------------- With Grace in Linux ----------------------------
1. Download the xvg files in your local computer.
2. Open Terminal and `cd` to the folder where you downloaded the files.
3. Type

		xmgrace <filename>.xvg


-------------- With Matplot (Jupyter notebook)-----------------
1. Install jupyter notebook in a conda environment by typing, `conda install jupyter notebook` or `pip install jupyter notebook`.
2. Open jupyter notebook and open Python3 kernel (from 'New').
3. Paste the lines in the [matplot.py](matplot.py) or [matplo_multidata.py](matplot_multidata.py) and click `run`.
(Change the filename, x label, y label and title of plot according to your requirement)





