#!/bin/bash

gen_cg_dist (){
        #In the CG DIR.
        echo "Generating CG distributions"
        set -e
        CG_T0=$1   #First argument specifies how many last nanoseconds ti use to calculated pei.xtc
	CG_TN=$2
	AA_dir=$3

        source ${AA_dir}/init.sh &>>output.log
	nbead=`grep 'Bead' ${AA_dir}/aa2cg.ndx | wc -l`
	n=$(( nbead - 1 ))
	echo "2" | gmx traj -f md_1.xtc -s md_1.tpr -oxt pei.xtc -b $CG_T0 -e $CG_TN &>>output.log 
	##Check what is "2"
	mkdir -p  bonded_distribution

	gmx polystat -s md_1.tpr -f pei.xtc -n ${AA_dir}/e2e.ndx -o polystat.xvg -mw &>>output.log

        if [ ! -f ${AA_dir}/polystat.xvg ] #AA polystat
        then
	    gmx polystat -s md_1.tpr -f ${AA_dir}/mapped_cg.* -n ${AA_dir}/e2e.ndx -o ${AA_dir}/polystat.xvg -mw &>> output.log 
	fi
	
	bonds=$(grep '\[' ${AA_dir}/cg_bonds.ndx | awk '{print $2}')
	i=0
	echo "Generating bond length distribution"
	for item in $bonds
	do
	    if [ ! -f bonded_distribution/bond_${item}.xvg ]
	    then
	        echo "$i" | gmx distance -f pei.xtc -n ${AA_dir}/cg_bonds.ndx -len 0.35 -tol 0.8 -oh bonded_distribution/bond_${item}.xvg &>>output.log
	    fi
	    i=$(( i + 1))
	done
	echo "DONE"
	
	
	angles=$(grep '\[' ${AA_dir}/cg_angs.ndx | awk '{print $2}')
	i=0
	echo "Generating bond angle distribution"
	for item in $angles
	do
	    if [ ! -f bonded_distribution/angle_${item}.xvg ]
	    then
 	    	echo "$i" | gmx angle -f pei.xtc -n ${AA_dir}/cg_angs.ndx -od bonded_distribution/angle_${item}.xvg &>>output.log
	    fi
	    i=$(( i + 1))
	done
	echo "DONE"

	dih=$(grep '\[' ${AA_dir}/cg_dihs.ndx | awk '{print $2}')
	i=0
	echo "Generating dihedral angle distribution"
	for item in $dih
	do
	    if [ ! -f bonded_distribution/dih_${item}.xvg ]
	    then
		    echo "$i" | gmx angle -type dihedral -f pei.xtc -n ${AA_dir}/cg_dihs.ndx -od bonded_distribution/dih_${item}.xvg &>>output.log
	    fi
	    i=$(( i + 1))
	done
	echo "DONE"
}

gen_aa_dist(){
        set -e
	TRR=$1
	TPR=$2
	T1=$3
	T2=$4
	AADIR=$5
        source $AADIR/init.sh &>>output.log
	#python aatop_2_cg.py #Generate CG itp file, CG mapping scheme, and reference CG distribution index file
	N=`grep 'Bead' $AADIR/aa2cg.ndx | wc -l`
	#Mapped CG trajectories
	typ=${TRR: -3}
	echo -e "\nMapping AA trajectories to CG: "
	if [ ! -f $AADIR/mapped_cg.$typ ]
	then
		if [ ${TRR:0:1} == '.' ]
		then
    		     seq 0 $((N-1)) | gmx traj -f $AADIR/$TRR -s $AADIR/$TPR -oxt $AADIR/mapped_cg.$typ -n $AADIR/aa2cg.ndx -com -ng $N -b $T1 -e $T2 &>> output.log
		else 
    		     seq 0 $((N-1)) | gmx traj -f $TRR -s $TPR -oxt $AADIR/mapped_cg.$typ -n $AADIR/aa2cg.ndx -com -ng $N -b $T1 -e $T2 &>> output.log
		fi

	fi
	echo "DONE"

        #last timestep gro
	echo "Configuring last trajectory of AA to CG: "
	if [ ! -f $AADIR/cg_ini.gro ]
	then
		if [ ${TRR:0:1} == '.' ]
		then
		      seq 0 $((N-1)) | gmx traj -f $AADIR/$TRR -s $AADIR/$TPR -oxt $AADIR/cg_ini.gro -n $AADIR/aa2cg.ndx -com -ng $N -b $T2 &>> output.log
		else
		      seq 0 $((N-1)) | gmx traj -f $TRR -s $TPR -oxt $AADIR/cg_ini.gro -n $AADIR/aa2cg.ndx -com -ng $N -b $T2 &>> output.log
	        fi
	fi
	echo "DONE"

	mkdir -p bonded_distribution
	bonds=$(grep '\[' $AADIR/cg_bonds.ndx | awk '{print $2}')
	echo "Generating reference bond length distribution: "
	i=0
	for item in $bonds
	do
	    if [ ! -f $AADIR/bonded_distribution/bond_${item}.xvg ]
	    then
	        echo "$i" | gmx distance -f $AADIR/mapped_cg.${typ} -n $AADIR/cg_bonds.ndx -len 0.35 -tol 0.8 -oh $AADIR/bonded_distribution/bond_${item}.xvg &>> output.log
	    fi
	    i=$(( i + 1))
	done
	echo "DONE"
	
	angles=$(grep '\[' $AADIR/cg_angs.ndx | awk '{print $2}')
	i=0
	echo "Generating reference bond angle distribution: "
	for item in $angles
	do
	    if [ ! -f $AADIR/bonded_distribution/angle_${item}.xvg ]
	    then
	        echo "$i" | gmx angle -f $AADIR/mapped_cg.${typ} -n $AADIR/cg_angs.ndx -od $AADIR/bonded_distribution/angle_${item}.xvg &>> output.log
	    fi
	    i=$(( i + 1))
	done
	echo "DONE"
	
	dih=$(grep '\[' $AADIR/cg_dihs.ndx | awk '{print $2}')
	i=0
	echo "Generating reference dihedral angle distribution: "
	for item in $dih
	do
	    if [ ! -f $AADIR/bonded_distribution/dih_${item}.xvg ]
	    then
	        echo "$i" | gmx angle -type dihedral -f $AADIR/mapped_cg.${typ} -n $AADIR/cg_dihs.ndx -od $AADIR/bonded_distribution/dih_${item}.xvg &>> output.log
	    fi
	       i=$(( i + 1))
	done
	echo "DONE"
	
}

run_CG_sim () {
        CGDIR=$1 #CG dir relative to current dir
	AADIR=$2 #AA relative to CGDIR or absolute path
	CG_T0=$3
	CG_TN=$4
	GPU=$5
	ION=$6
	EM=$7
        NPT=$8
	MD=${9}
        WAT=${10}

        echo " "
        source init.sh &>> output.log
        #In main directory
	if [ -f $CGDIR/md_1.gro ] #Constrained NPT simulation has finished.
	then
		echo "Skiping $CGDIR Simulations: Finished previously"
		cd $CGDIR
		#Generate CG bonded distributipns
		if [ ! -f comparing_pngs/all_images.pdf ] #If comapring pdf is not generated
		then
	                gen_cg_dist $CG_T0 $CG_TN $AADIR
		fi
		cd $AADIR #In dir main
	else
		echo "Starting to run $CGDIR Simulations"
        	cd $CGDIR #In dir CG*
        	cp $AADIR/cg_ini.gro .

		#Create a new box
		echo "Creating a new box"
                cp cg_ini.gro cgPEI_newbox.gro
        	#gmx editconf -f cg_ini.gro -o cgPEI_newbox.gro -c -d 2.0 -bt cubic &>> output.log
		#Change the molecule name in the topology file (as the itp file has different molecule name) and copy 
        
		#Add pokarizable water to the initial structure
		echo "Adding polarizable water"
		gmx solvate -cp cgPEI_newbox.gro -cs ${AADIR}/${WAT} -o cgPEI_solv.gro -p cgtopol.top -radius 0.12 >> output.log 2>&1
                n=`grep 'WP' cgPEI_solv.gro | wc -l`
		if [ -f cgtopol.top ]
		then
		    sed -i '/^PW.*/d' cgtopol.top
		    sed -i '/^CL.*/d' cgtopol.top
       	 	    echo "PW $n" >> cgtopol.top #Append the number of water molecules in the topology file 
		else
		    exit
		fi
		#Add neutralizing ions to the system
		echo "Adding ions"
        	gmx grompp -f $AADIR/$ION -c cgPEI_solv.gro -p cgtopol.top -o ions.tpr -maxwarn 1 &>> output.log
        	echo "3"| gmx genion -s ions.tpr -o cgPEI_solv_ions.gro -p cgtopol.top -pname K+ -nname CL- -neutral  &>> output.log
		
                sed 's/martini_v2\.2refP\.itp/martini_v2\.2refP_constr\.itp/g' cgtopol.top > cgtopol_constr.top	
        
		#Energy minimization
		echo "Running energy minimization"
		if [ ! -f em.tpr ]
		then
        	    gmx grompp -f $AADIR/$EM -c cgPEI_solv_ions.gro -p cgtopol_constr.top -o em.tpr -maxwarn 1 &>> output.log
		fi

                for i in {1..5}
		do
                    if [ ! -f em.gro ]
		    then
		        if [ $GPU != "1" ]
	      	        then
     	      	            echo "runnning on CPU"
	                else 
	      	            echo "running on CPU+GPU"
			fi
                        export OMP_NUM_THREADS=32
	      	        gmx mdrun -v -deffnm em -ntomp 32  &>> output.log
		    else 
		        break
	      	    fi
		done

		#NPT equilibriation with constraints
		echo "Running NPT equilibration"
		for i in {1..5}
		do
		    if [ ! -f npt.gro ]
		    then
        	          gmx grompp -f $AADIR/$NPT -c em.gro -p cgtopol.top -o npt.tpr &>> output.log
                          export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"
		          srun gmx_mpi mdrun -deffnm npt -v  &>> output.log
                          if [ $? -ne 0 ]
                          then
                               export OMP_NUM_THREADS=8
		               gmx mdrun -deffnm npt -v -ntomp 8 &>> output.log
                          fi
		    else
		        break
		    fi
		done

		#NPT without constraints
		echo "Running NPT Production Run"
		for i in {1..5}
		do 
		    if [ ! -f md_1.gro ]
		    then
        	         gmx grompp -f $AADIR/$MD -c npt.gro -p cgtopol.top -o md_1.tpr  &>> output.log
                         export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"
                         srun gmx_mpi mdrun -v -deffnm md_1 &>> output.log
                         if [ $? -ne 0 ]
                         then
                             export OMP_NUM_THREADS=8
		             gmx mdrun -ntomp 8 -deffnm md_1 -v  &>> output.log
                         fi
		    else
		        break
		    fi
		done
       
                gen_cg_dist $CG_T0 $CG_TN ..
		cd $AADIR #In dir main
	fi
}

conv_itp (){
     fname=$1

     N=`grep -n '\[.*constraints.*\]' martini_v${fname}.itp | sed 's/:.*//g'`
     N=$((N-1))

     head -n $N martini_v${fname}.itp > martini_v${fname}_constr.itp     
     echo "[ bonds ]" >> martini_v${fname}_constr.itp
     echo ";  i     j   funct   length   force const.;" >> martini_v${fname}_constr.itp
     echo "   1     2    1       0.14    100000" >> martini_v${fname}_constr.itp
     echo "   1     3    1       0.14    100000" >> martini_v${fname}_constr.itp
     echo " "  >> martini_v${fname}_constr.itp
     echo "[ angles ]" >> martini_v${fname}_constr.itp
     echo ";   i    j   k   funct  angle    fc" >> martini_v${fname}_constr.itp
     echo "    2    1   3    2     0.0     4.2" >> martini_v${fname}_constr.itp
     echo " "  >> martini_v${fname}_constr.itp
     echo "[ exclusions ]" >> martini_v${fname}_constr.itp
     echo "1   2 3" >> martini_v${fname}_constr.itp
     echo "2   3" >> martini_v${fname}_constr.itp

}       

gen_edr_prop (){
        CG_T0=$1   #First argument specifies how many last nanoseconds ti use to calculated pei.xtc
	CG_TN=$2
        Prop=$3
        source init.sh &>>output.log

	echo "$Prop" | gmx energy -f md_1.edr -o $Prop.xvg -b $CG_T0 -e $CG_TN &>>output.log 
     
}

gen_msd (){
        CG_T0=$1   #First argument specifies how many last nanoseconds ti use to calculated pei.xtc
	CG_TN=$2
        MOL=$3
        source init.sh &>>output.log
        echo "$MOL" | gmx msd -f md_1.trr -s md_1.tpr -o msd_${MOL}.xvg -b $CG_T0 -e $CG_TN -mol -mw &>> output.log 

}

gen_msd_N (){
        CG_T0=$1   #First argument specifies how many last nanoseconds ti use to calculated pei.xtc
	CG_TN=$2
        MOL=$3
        N=$4
        echo "$CG_T0 $CG_TN $MOL $N"
        source init.sh &>>output.log
        for (( i=1; i<=$N; i++))
        do
            CG_T1=$((CG_T0+(i-1)*(CG_TN-CG_T0)/N))
            CG_T2=$((CG_T0+i*(CG_TN-CG_T0)/N))
            a=$((CG_T1/1000))
            b=$((CG_T2/1000))
            echo "$MOL" | gmx msd -f md_1.trr -s md_1.tpr -o msd_${MOL}_${a}-${b}ns.xvg -b $CG_T1 -e $CG_T2 -mol -mw &>> output.log
        done 

} 
 
"$@"
