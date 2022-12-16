# Martini\_PEI\_ff

## References

If you use newly generated forcefield parameters or the automated coarse-graining algorithm, please cite,

[1] Subhamoy Mahajan and Tian Tang, "Automated Parameterization of Coarse-Grained Polyethylenmine Under Martini Framework", *ChemRxiv* **2022** DOI:[10.26434/chemrxiv-2022-1fp1v](https://doi.org/10.26434/chemrxiv-2022-1fp1v)

If you are using old forcefield parameters, please cite,

[2] Subhamoy Mahajan and Tian Tang, "Martini Coarse-Grained Model for Polyethylenimine", *J. Comput. Chem.* **2019**, 40 (3), 607-618. 

All-atom simulation references.

[3] Chongbo Sun, Tian Tang, and Hasan Uludag, Javier E. Cuervo, "Molecular Dynamics Simulations of DNA/PEI Complexes: Effect of PEI Branching and Protonation State", *Biophys. J.* **2011**,100(11), 2754-2763.  

[4] Chongbo Sun, Tian Tang, and Hasan Uludaǧ, "Molecular Dynamics Simulations for Complexation of DNA with 2 KDa PEI Reveal Profound Effect of PEI Architecture on Complexation." *J. Phys. Chem. B* **2012**, 116 (8), 2405–2413.


## Installation

```bash
python setup.py install
```

For copyright issues, the martini topology files are not provided with the python module. Please do the following in your local copy of coarsen. This needs to be done only once.

1) Copy all .itp files from Maritni website and place them the coarsen folder.
2) For itp files `martini_v2.2refP.itp`, `martini_v2.2P.itp`, etc. make a copy as `martini_v2.2refP_constr.itp`, `martini_v2.2P_constr.itp`, etc. These files will be used in energy minimization.
3) For energy minimization .itp files you can connect the constraints portion of polarizable water models, and replace them with harmonic bonds. For example, in `martini_v2.2refP_constr.itp`

convert the portion below,
  
```bash
[constraints]
;  i     j   funct   length
   1     2    1       0.14
   1     3    1       0.14

; for minimization purposes constraints might be replaced by stiff bonds:
;
;[bonds]
;  i     j   funct   length   force const.
;   1     2    1       0.14    50000
;   1     3    1       0.14    50000

``` 

to, 
 
```bash
;[constraints]
;  i     j   funct   length
;   1     2    1       0.14
;   1     3    1       0.14

; for minimization purposes constraints might be replaced by stiff bonds:
;
[bonds]
;  i     j   funct   length   force const.
   1     2    1       0.14    100000
   1     3    1       0.14    100000
```
You can alter the force constant to reach a stable initial configuration.


## Quickstart
### 1. Generate CG reference trajectories

Convert AA trajectories to reference CG trajectories.

```bash 
coarsen aa2cg -f AA.trr -s AA.tpr -p AA.top -b T1 -e T2 -n aa2cg.ndx -o cgtopol.top -x parameters.dat 
``` 

> **_Note1:__** T1 and T2 are time in ps (integer).
>
> **_Note2:__** `aa2cg.ndx` and `cgtopol.top` are output files.
> 
> **_Note3:__** parameter `peiname`, `init`, `bond_small`, `bond_large`, `bond_ymax`, `ang_ymax`, `dih_ymax` are required in `parameters.dat`

Reads:
- AA files (`AA.trr` or `AA.xtc`, `AA.tpr`, `AA.top` or `AA.itp` (of only PEI)).
- Parameters (`parameters.dat`)


Writes:
- CG mapping scheme (`aa2cg.ndx`)
- AA and CG molecular structure (`aa_struct.pickle`, `cg_struct.pickle`)
- groupings for CG bonded distributions (`cg_bonds.ndx`, `cg_angles.ndx`, `cg_dihs.ndx`)
- pickled list of parmeterized (`param.pickle`) and unparameterized (`unparam.pickle`) bonded distributions
- AA bonded distributions (`Bonded_distributions/`). Contains individual .png and combined .pdf.
- Curve-fits of AA distributions and predictions of unparameterized paramters (`ref_param_fit/`)
- CG forcefield  (`CG1/cgff1.pickle`) and topology (`CG1/cgtopol.top`, `CG1/peiname1.itp`) for first iteration.   
## 2. Iteratively determine bonded CG parameters

```bash
coarsen parameterize -x parameters.dat -b T1 -e T2
```
> **_Note1:__** T1 and T2 are CG time in ps (integer).
>
> **_Note2:__** parameter `start_iter`, `max_iter`, `fa`, `wa`, `thc`, `kc`, `fd`, `peiname`, `init`, `bond_small`, `bond_large`, `bond_ymax`, `ang_ymax`, `dih_ymax`, `gpu` (default: 0), `em_mdp` (default: em.mdp), `ions_mdp` (default: ions.mdp), `npt_mdp` (default: npt.mdp), `md_mdp` (default: md.mdp ) are required in `parameters.dat`


Reads:
- CG structure (`cg_struct.pickle`)
- `unparam.pickle`
- .mdp files specified in `parameters.dat`.

Writes:
- GROMACS output and log files for each iteration `CG*` and subiterations `CG*_th` and `CG*_K`.
- Bonded distribtuions `CG*/Bonded_distribution/`
- Comparison of reference CG and CG distributions `CG*/comparing_pngs/*.png`, `CG*/comparing_pngs/all_images.pdf`
- Updates `Cost.pickle` for the current iterations.


### 3. Generate CG distributions

```bash
coarsen gen_cg_dist -b T1 -e T2 -d directory
```

Writes:
- Bonded probability distributions `bonded_distribution/*.xvg`
- GROMACS .xtc file containing only pei trajectories.
- Radius of gyration and end-to-end distance of the CG-polymer (`polystat.xvg`) and AA-polymer (`../polystat.xvg`) 


### 4. Calculate cost 

```bash
coarsen add_cost -x parameters.dat -k string -f directory
```
            
### 5. SMILE to CG simulation files

```bash
coarsen smile2cg SMILE_STRING 
``` 
String should contain `t`, `sq`, `s`, `pq`, `p` representing tertiary, protonated-secondary, secondary. protonated-primary, and primary beads respectively. 

Branches can be specified between `(` and `)` blocks.

Additionally repeating blocks can be specified between `[` and `]`.

### 6. Generate Reports

To generate PNGs of all bonded distribution and a pdf containing all bonded distribution use,

```bash
coarsen gen_report -x parameters.dat
```
The command reads bonded distributions (saved as .xvg) in `dir/bonded_distribution/`. All the directories `dir` should be specified 
in the `dirs` parameter. 
The legends can be specified using `labels` and number of legend columns using `ncols`. Other parameters read being used are: 
`bond_small`, `bond_large`, `bond_ymax`, `ang_ymax`, `dih_ymax`.
 

### 7. Calculate PEI properties

All properties are automatically stored in `prop.pickle` in the local directory.

```bash
coarsen get_prop pei -f cg_struc.pickle 
```
Reports molecular properties of PEI. 

```bash
coarsen get_prop polystat
```
Reports average and standard deviation of end-to-end distance and radius of gyration.

```bash
coarsen get_prop edr [prop] -b T1 -e T2 -x parameters.dat
```
Calculates a property from `md_1.edr` file. This command runs the `gmx energy` and stores the average and std of the property in `prop.pickle`. It has currently been tested for calculating volume (`[prop] = Volume`).

```bash
coarsen get_prop conc -i nmol
``` 
Calculates concentration of PEI in the simulation. Note that all PEIs should be identical for this calculation. Property of PEI (`get_prop pei`) and volume of simulation box (`get_prop edr Volume`) must have been previously executed. 

`-i nmol` is the number of PEI molecules.

```bash
coarsen get_prop diff_scaling -x parameters.dat
```

Determines the diffusion scaling factor. Currently works for Martini 2.1P and Martini Ref2.2P. Calculation of concenteration must have been previously calculated. 

The parameter `martini ` = `2.2refP`, `2.1P`, or `2.2P` is required. The default parameter is `2.1P`.

```bash
coarsen get_prop diff -b T1 -e T2 -k mol_name -n tfit
```

Calculates diffusion coefficient of PEI. Calculation of diffusion scaling factor must have been previously calculated.

- `-k mol_name`: Molecule name of PEI.
- `-n tfit`: The fitting time in nanoseconds.
- `-s show`: Shows the output rather than saving it as PNG.

```bash
coarsen get_prop diff_N -b T1 -e T2 -k PEI -n tfit -i partitions
```

Calculates diffusion coefficient of PEI by creating `partitions` partitions of the time `T1`-`T2`. This is used to calculate 
the standard deviation of diffusion coefficient.

- `-s show`: Shows the output rather than saving it as PNG.
- `-k [name]`: Specify the name of the molecule in GROMACS index files.
 

## GROMACS type options

- `-f` : Input file name, typically for .xtc .trr
- `-s` : tpr file name  (.tpr)
- `-p` : topol file name (.topol or .itp)
- `-n` : index file (.ndx)
- `-b` : Being time of the simulation.
- `-e` : End time of the simulation.

## Module specific options

- `--paramfile` : Parameter file name for coarsen module
- `--output` : Outpul file name
- `--dir` : Name of a directory
- `--key` : A string input
- `--id` : Integer input 

## Parameters for paramfile

- `init` : Text separated by ';' (end-of-line) for loading modules associated with Gromacs.
- `cgff_curr` : Path to current CG forcefield parameters as pickled file. Default value is the 2019 forcefield.  
- `bond_small` : Smallest bond length for plotting. Default is 0.3 nm.
- `bond_large` : Largest bond length for plotting. Default is 0.6 nm.
- `bond_ymax` : Largest probability density value of bond length distribution for plotting.
- `ang_ymax` : Largest probability density value of bond angle distribution for plotting.
- `dih_ymax` : Largest probability density value of dihedral angle distribution for plotting.
- `dih_initial` : text file containing dihedral angle name followed by the number of periodic potential to use. 
   If the file is not specified default value of 8 is used for each dihedral angle. 
- `dirs`: Name of directories to generate report of bonded distribution.
- `labels`: Labels to use in legend of bonded distribution report.
- `ncols`: Number of columns to use in legends for bond length, bond angless, and dihedral angles respectively. 
- `gpu`: 1 if simulations are performed using GPU, 0 (default) if it is not.
- `cost_tol`: Default is 1E-8.
- `ions_mdp`: GROMACS .mdp file for adding ions
- `em_mdp`: GROMACS .mdp file for energy minimization 
- `npt_mdp`: GROMACS .mdp file for constrainted NPT simulation.
- `md_mdp`: GROMACS .mdp file for unconstrainted NPT simulation.
- `fb` : Gradient descent step for bond length parameter optimization. Not used.
- `wb` : Weight assigning relative importance to mean or standard deviation of bond length distributions. 1.0 implies standard deviation is not important, and 0.5 implies both mean and standard deviation are equally important. Typically, 0.75 is prefered.
- `wa` : Weight assigning relative importance to mean or standard deviation of bond angle distributions.
- `fa` : Gradient descent step for bond angle parameter optimization.
- `thc` : Change in equilibrium theta parameter for calculating Jacobian
- `kc` : Change in angle force constant parameter for calculating Jacobian
- `fd` : Gradient descent step for dihedral angle parameter optimisation. 
- `pos_prec` : Default is 3, based on GROMACS .gro files. 
- `start_iter` : Starting iteration number. Should be greater than 0.
- `max_iter` : Last iteration number. Should be greater than 1.
- `peiname` : Name of PEI.

## Pickle formats

### cgff.pickle
- `info`: String containing the storage format.
- bonded distribution name: For bond lengths (b0,Kb), for bond angles (a0,Ka), Dihedrals [number of functions, list of [phi, Kd, n]] 
              
### Cost.pickle
- Keys are iteration numbers (natural numbers) or custom strings (Ex. 'Best').
- Cost[key] has several keys: 
    - `all`, `param`, `unparam`: Cost associated with all, parameterized or unparameterized bonded distribution.
    - `all_bond`, `all_ang`, `all_dih`: Cost associated with either bond length, bond angle or dihedral angle distriutions.
    - `param_bond`, `param_ang`, `param_dih`, `unparam_bond`, `unparam_ang`, `unparam_dih`: Similar as above.
    - bonded distributions: Cost associated with a specific bonded distribution name.
    - `Rg` or `Re`: Cost associated with radius of gyration or end-to-end distance.

### param.pickle
- `bonds`: list of parameterized bond lenght names
- `angs`: list of parameterized bond angle names
- `dihs`: list of parameterized dihedral angle names

### unparam.pickle
- `bonds`: A dictionary of bond length names : associated gradient descent step modifier.
- `angs`: A dictionary of bond angle names : associated gradient descent step modifier.
- `dihs`: A dictionary of dihedral angle names: [ number of functions, associated gradient descent step modifier]

### aa\_struct.pickle
- A networkx Graph(). 
- Nodes are atom indeces. Node attributes include `atom_name`, `atom_type`, `res_name`, `charge`, and `mass`.
- Edges are covalent bonds.

### cg\_struct.pickle
- A networkx DiGraph(). 
- Nodes are atom indeces. Node attributes include `bead_name`, `charge`, and `mass`.
- Edges are covalent bonds. The bonds start from N terminal of one bead to C terminal of another.
 
### prop.pickle
- `MolWt`: Molecular weight of PEI in Da.
- `Charge`: Charge of PEI.
- `pr`: Protonation ratio.
- `N_tsp`: Number of t, s, and p beads. sq and pq are treated as s and p respectively.
- `Re`: end-to-end distance avg, std 
- `Rg`: radius of gyration avg, std 
- `Rg_eig`: index 0 is list of average eigen values, index 1 is list of std of eigen values.
- `shape`: `asphericity`,`acylindricity`, relative shape anisotropy `shape_anisotropy`.
- `Volume`: Volume of simulation box averaged over user-specified simulation time range (nm3), std.
- `conc`: Concentration of PEI (g/L).
- `D`: Diffusion coefficent (1E-5 cm2/s). 
- `D_err`: Error in diffusion coefficent (1E-5 cm2/s). 
- `D_scaled`: Scaled diffusion coefficent (1E-5 cm2/s).  
- `D_scaled_err`: Error in scaled diffusion coefficent (1E-5 cm2/s). 
