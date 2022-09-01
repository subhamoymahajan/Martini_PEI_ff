#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=4
#SBATCH --time=03:30:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=subhamoy@ualberta.ca
module load gromacs/5.1.4
module load scipy-stack
export  GMX_MAXCONSTRWARN=-1
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"
set -e

run_CG_sim () { #Argument 1 is $IDX, and Argument 2 if given is _th or _K
	#You should be in jacobian_1 folder right now 
        #In dir main
        #Copy the initial configuration
        cp cg_ini.gro CG$1$2/
        cd CG$1$2 #In dir CG*

	#Create a new box
        gmx editconf -f cg_ini.gro -o cgPEI_newbox.gro -c -d 1.5 -bt cubic
	#Change the molecule name in the topology file (as the itp file has different molecule name) and copy 
        sed "s/bPEI2K/bPEI2K_$1$2/g" ../common/topol.top > topol.top
	#Add pokarizable water to the initial structure
        gmx solvate -cp cgPEI_newbox.gro -cs ../common/polarize-water.gro -o cgPEI_solv.gro -p topol.top
        
        n=`grep 'WP' cgPEI_solv.gro| wc -l`  #Count the number of water molecules added
        echo "PW $n" >> topol.top #Append the number of water molecules in the topology file 

	#Add neutralizing ions to the system
        gmx grompp -f ../common/ions.mdp -c cgPEI_solv.gro -p topol.top -o ions.tpr -maxwarn 1
        echo "3"| gmx genion -s ions.tpr -o cgPEI_solv_ions.gro -p topol.top -pname K -nname CL -neutral
        
	#Energy minimization
        gmx grompp -f ../common/em.mdp -c cgPEI_solv_ions.gro -p topol.top -o em.tpr -maxwarn 1
        srun gmx_mpi mdrun -v -deffnm em
        
	#NPT equilibriation with constraints
        gmx grompp -f ../common/npt.mdp -c em.gro -p topol.top -o npt.tpr
        srun gmx_mpi mdrun -v -deffnm npt -rdd 2
        
	#NPT without constraints
        gmx grompp -f ../common/md.mdp -c npt.gro -p topol.top -o md_1.tpr 
        srun gmx_mpi mdrun -v -deffnm md_1 -rdd 2
        
	#Generate CG bonded distributipns
        bash ../common/calc_dist.sh 
	python ../common/gen_png.py
	cd comparing_pngs  #In dir CG*/comparing_pngs
	python ../../common/create_latex.py
	pdflatex all_images.tex
	cd ../../ #In dir main

}

run_CG_sim 4
        
        
        

