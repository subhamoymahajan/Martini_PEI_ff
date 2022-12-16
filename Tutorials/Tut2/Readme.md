# Tutorial 2: Iteratively Coarse-Grain PEI

## a. Convert AA topology to CG topology

```bash
coarsen aa2cg -f md_1_reduced.xtc -s md_1_pei.tpr -p aa_pei.itp -b 50000 -e 100000 -n aa2cg.ndx -o cgtopol.top -x parameters.dat 
```

Parameters used: `martini`, and `pei_name`. 

Optional parameters: `cgff_cur`, `dih_initial`, and `init`.


Ouput generated:
- AA to CG mapping scheme: `aa2cg.ndx`
- Index files for bonded distributions: `cg_bonds.ndx`, `cg_angs.ndx`, `cg_dihs.ndx`.
- Reference distributions `bonded_distributions` and figures `ref_param_fit`.
- Index fiels for end-to-end distance of : `e2e.ndx`.
- ITP files: `martimi_v2.2refP.itp`, `martini_v2.2refP_constr.itp`, `martini_v2.0_ions.ipt`.
- networkx graphs for AA and CG structures: `aa_struct.pickle`, `cg_struct.pickle`
- Mapped CG trajectories: `mapped_cg.xtc`, `cg_ini.gro`.
- Parameterized and unparameterized bonded distributions: `param.pickle`, `unparam.pickle`.
- Gromacs loading script : `init.sh`.
- log file for debugging `output.log` 


## b. Parameterize

Name the mdp files for energy minimization, NPT equilibration, and MD production as `em.mdp`, `npt.mdp`, and `md.mdp`. Otherwise,
specify the mdp file names as parameters `em_mdp`, `npt_mdp`, and `md_mdp`. Three iterations take with this small system provided 
takes ~30 min with 32 CPUs.


```bash
coarsen parameterize -x parameters.dat -b 2500 -e 5000
``` 

Output generated: 
- `CG[1-3]` iteration folders.
- `CG[1-3]_th`, `CG[1-3]_K` sub-iteration folders.
- Cost funciton values associated with each iteration : `Cost.pickle`
- End-to-end distance and radius gyration of simulations: `polystat.xvg` , `CG[1-3]/polystat.xvg` , etc.
- Each `CG[1-3]` folde contains:
    - Reference distributions `bonded_distributions` and comparison figures `comparing_pngs`.
    - CG forcefield: `cgff1.pickle`, `tutPEI1.itp`.
    - intial strucutre `cg_ini.gro`, boxed structure `cgPEI_newbox.gro` (currently identical to initial structure), solvated PEI `cgPEI_solv`, and solvated-neutralized PEI `cgPEI_solv_ions.gro`
    - Energy minimization files `em.*`, constrained NPT simulation files `npt.*`, MD production files `md_1.*`.
    - trajectories of only PEI `pei.xtc`. 
    - log file for debugging `output.log` 


Several iterations and long simulation time would typically be required to coarse grain a new bonded parameters. For generating new bonded parameters for publication use the lastest forcefield file.
