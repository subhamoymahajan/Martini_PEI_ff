# Martini\_PEI\_ff

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

Additionally repeating blocks can be specified between `{` and `}`.

### 6. Calculate PEI properties

```bash
coarsen get_prop pei 
```
Reports molecular properties of PEI. Also updates the properties in `prop.pickle`.

```bash
coarsen get_prop polystat
```
Reports average and standard deviation of end-to-end distance and radius of gyration.
Also updates the properties in `prop.pickle`.

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
- `dih_initial` : 
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
- `Re`: end-to-end distance avg, std 
- `Rg`: radius of gyration avg, std 
- `Rg_eig`: index 0 is list of average eigen values, index 1 is list of std of eigen values.
- `shape`: `asphericity`,`acylindricity`, relative shape anisotropy `shape_anisotropy`. 
