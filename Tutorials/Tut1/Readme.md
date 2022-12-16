# Tutorial 1: CG simulations using existing forcefield parameters

### SMILE strings of PEIs in [2,3]
**600 Da**

PL23: `ssqssssssqsssspq`

SL23: `t(pq)sst(pq)ssst(pq)sp`

MB23: `sqst(ssspq)sst(p)spq`

HB23: `st(t(p)spq)t(spq)t(p)spq`

PL46: `ssqssqssqssqssqsspq`

SL46: `t(pq)ssqt(pq)ssqst(pq)spq`

MB46: `sqst(sqsspq)sqst(p)spq`

HB46: `sqt(t(pq)spq)t(spq)t(pq)spq`

### SMILE strings of PEI in [1,4]

**2 kDa**

t(spq)sqt(sspq)t(pq)t(sqt(pq)spq)sqt(spq)t(sspq)t(pq)t(sqt(pq)spq)sqt(sqt(pq)spq)t(spq)pq

### SMILE string of PEI in [1]

**2 kDa - MilliporeSigma**

SL2k: [[t(spq)]4t(sspq)t(spq)]2t(spq)p

MB2k: t(spq)[t(t(sspq)t(spq)pq)]2t(t(sspq)t(spq)p)t(t(sspq)t(spq)pq)sspq

HB2k: t(spq)t(t(spq)t(sspq)t(spq)pq)t(t(spq)t(sspq)t(spq)p)t(t(spq)t(sspq)t(spq)pq)sspq

**25 kDa - MilliporeSigma**

SL25k: t(sspq)[t(spq)t(sspq)[t(spq)]2]17t(sspq)sspq

MB25k: [t([t(sspq)]2t(spq)pq) t(spq)t([t(sspq)]2t(spq)pq)t(sspq)]7pq

HB25k: [t([[t(spq)]3t(sspq)]2sspq)t(spq)]7pq

## Generating topology file

```bash
coarsen smile2cg 't(pq)ssqt(pq)ssqst(pq)spq' -x parameters.dat -o cg_ini.gro -p cg_topol.top
```

- Creates topology file `cg_topol.top`. 
- The ITP file name of PEI, `tutPEI.itp`, is based on the `peiname` parameter in `parameters.dat`.

- The parameter `init` stores a command that is needed to load Gromacs or any other necessary softwares in your local system. The tutorial uses commands necessary for compute canada servers. 
The commands are stored in a local file `init.sh` that is run every time Gromacs commands are necesary.
If multiple commands are necessary to load Gromacs, use `init = cmd1 ; cmd2 ; cmd3`.


- `cgff_curr = cgff_2022` or `cgff_curr = cgff_2019` uses the pickled CG forfield file in the `coarsen` folder. To use a custom forcefield parameters supply the path of a pickled file 'cgff = path/myff.pickle' 

- Martini version and water model can be defined using `martini = 2.2refP`, `martini = 2.2P` or `martini = 2.1`, etc.

- Index files are generated for evaluting bonded CG distributions: `cg_bonds.ndx`, `cg_angs.ndx`, and `cg_dihs.ndx`.

- Bonded distributions with known parameters are stored in `param.pickle` and unknown parameters are stored in `unparam.pickle`.

- To determine end-to-end distance of the longest linear chain, `e2e.ndx` is generated.

- ITP files from the `coarsen` folder are copied to local directory.

- The directed graph of the CG-PEI is stored in `cg_struct.pickle`.


Now, one can run desired simulations with .mdp files.
