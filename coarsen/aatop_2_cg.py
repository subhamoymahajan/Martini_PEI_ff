#    Copyright 2022 SUBHAMOY MAHAJAN 
#    
#    This file is part of Martini_PEI_ff software.
# 
#    InSilicoMicroscopy is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.)

import networkx as nx
import numpy as np
from scipy.interpolate import interp1d
import os 
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
tol=1E-6 # A very small number

def find_bonds(CGStruc):
    """ Find bonds from a CG networkx structure

    Parameters
    ----------
    CGStruc: networkx DiGraph()
        Nodes are numbered from 1 to n. Each node has attribute 'bead_name',
        which can take values 't', 'tq', 's', 'sq', 'p', or 'pq'/

    Returns
    -------
    bonds: 2D list
        each item is a list containing index of two nodes and the corresponding
        bond name. [i, j, name]
    """
    bonds=[]
    for ni in CGStruc:
        next_beads=list(CGStruc[ni])
        for nj in next_beads:
            name=CGStruc.nodes[ni]['bead_name']+CGStruc.nodes[nj]['bead_name']
            bonds.append([ni, nj, name])
    return bonds

def find_angles(CGStruc):
    """ Find angles from CG networkx structure

    Parameters
    ----------
    CGStruc: networkx DiGraph()
        Nodes are numbered from 1 to n. Each node has attribute 'bead_name',
        which can take values 't', 'tq', 's', 'sq', 'p', or 'pq'

    Returns
    -------
    bonds: 2D list
        each item is a list containing index of two nodes and the corresponding
        angle name. [i, j, k, name]

    """
    angles=[]
    #types='pst'
    types=['pq', 'p', 'sq', 's', 'tq', 't' ]
    for nj in CGStruc: # iterate over beads that would be the central bead in
                       # the bond angle.
        prev_n=list(CGStruc.predecessors(nj)) # This would be bead 'i' in angle
                                              # ijk
        succ_n=list(CGStruc.successors(nj)) # This would form bead 'k' in angle
                                            # ijk
        namej=CGStruc.nodes[nj]['bead_name'] 

        #normal angles i->j->k
        if len(prev_n)==1:
            namei=CGStruc.nodes[prev_n[0]]['bead_name'] #each bead has at most
                                                        #one predecessor
            for nk in succ_n: 
                namek=CGStruc.nodes[nk]['bead_name']
                angles.append([prev_n[0],nj,nk,namei+namej+namek])

        #N-type angles
        if len(succ_n)==2: # N-type angles only possible for tertiary central 
                           # bead, which would have two successors.
            namei=CGStruc.nodes[succ_n[0]]['bead_name']
            namek=CGStruc.nodes[succ_n[1]]['bead_name']
            # i<-j->k can be named both Nijk or Nkji. 
            # To make the name unique, degree of i is less than or equal to k.
#            if types.find(namei[0]) < types.find(namej[0]):
            
            if types.index(namei) < types.index(namek):
                angles.append([succ_n[0],nj,succ_n[1],'N'+namei+namej+namek])
            else:
                angles.append([succ_n[1],nj,succ_n[0],'N'+namek+namej+namei])

    return angles

def find_dihs(CGStruc):
    """ Find dihedral angles from CG networkx structure

    Parameters
    ----------
    CGStruc: networkx DiGraph()
        Nodes are numbered from 1 to n. Each node has attribute 'bead_name',
        which can take values 't', 'tq', 's', 'sq', 'p', or 'pq'

    Returns
    -------
    bonds: 2D list
        each item is a list containing index of two nodes and the corresponding
        angle name. [i, j, k, l, name]

    """
    dih=[]
    for nj in CGStruc: #iterate over bead 'j' in dihedral ijkl
        prev_n=list(CGStruc.predecessors(nj)) #bead 'i' in dihedral ijkl
        succ_n=list(CGStruc.successors(nj)) #bead 'k' in dihedral ijkl
        namej=CGStruc.nodes[nj]['bead_name']

        #normal angles i->j->k->l
        if len(prev_n)==1:
            namei=CGStruc.nodes[prev_n[0]]['bead_name'] # beads have at most 
                                                        # one predecessor.
            for nk in succ_n:
                namek=CGStruc.nodes[nk]['bead_name']
                succ_n2=list(CGStruc.successors(nk)) #bead 'l' in ijkl
                for nl in succ_n2:
                    namel=CGStruc.nodes[nl]['bead_name']
                    dih.append([prev_n[0],nj,nk,nl,namei+namej+namek+namel])

        #N-type angles i<-j->k->l or k<-j->i->l
        if len(succ_n)==2:
            namei=CGStruc.nodes[succ_n[0]]['bead_name']
            namek=CGStruc.nodes[succ_n[1]]['bead_name']
            #N(id0) i<-j->k->l
            succ_n2=list(CGStruc.successors(succ_n[1])) #bead 'l' in ijkl
            for nl in succ_n2:
                namel=CGStruc.nodes[nl]['bead_name']
                dih.append([succ_n[0],nj,succ_n[1],nl,'N'+namei+namej+namek+namel])

            #N(id2) k<-j->i->l
            succ_n2=list(CGStruc.successors(succ_n[0])) #bead 'l' in kjil
            for nl in succ_n2:
                namel=CGStruc.nodes[nl]['bead_name']
                dih.append([succ_n[1],nj,succ_n[0],nl,'N'+namek+namej+namei+namel])
    return dih

def topol2graph(topol): 
    """ Read AA Gromacs topology and convert to a network structure

    Parameters
    ----------
    topol: str
        Gromacs topology filename.

    Returns
    -------
    Structure: networkx Graph()
        Graph representation of AA-PEI topology.
        Nodes have attributes:
            'atom_name': name of the atom
            'atom_type': type of the atom (Currently not used)
            'charge': charge of the atom (Currently not used)
            'mass': mass of the atom 
            'res_name': residue name of the atom (Currently not used)
    """
    f=open(topol,"r")
    Structure=nx.Graph()
    read_cnt=0
    readid='None'
    for lines in f:
        foo=lines.split()
        if len(foo)==0:
            continue
        elif foo[0]==";" or foo[0]=="#" or foo[0][0]=="#": #Ignore comments
            continue
        elif foo[0]=="[":
            if readid in ["atoms", "bonds"]:
                #finished reading atoms or bonds
                read_cnt+=1
            readid=foo[1] #Currently reading "readid"

        elif readid=="atoms": #reading atoms information
            #foo[0] is ID of the atom
            #foo[1] is the atom type
            #foo[3] is the residue name
            #foo[4] is the atom name
            #foo[6] is the charge
            #foo[7] is the mass
            Structure.add_node(int(foo[0]), atom_name=foo[4], atom_type=foo[1], 
                 charge=float(foo[6]), mass=float(foo[7]), res_name=foo[3])

        elif readid=="bonds": #read covalent bond information.
            Structure.add_edge(int(foo[0]),int(foo[1]))

        if read_cnt==2: #read the atoms and bond information.
            break

    find_begin(Structure) #Find the begining Carbon.
    return Structure

def bound_atoms(Structure,node,atom):
    """ Cound the number of bound atoms of type 'atom'.

    Parameters
    ----------
    Structure: networkx Graph()
        AA Structure as graph.
    node: node of Graph()
        node for which the bound carbon atoms are counted. Each node has an 
        attribute 'atom_name' with data type str.
    atom: str
        The 'atom' type that is counted.

    Returns
    -------
    cnt: int
        Number of bound carbon atoms.
    """
    con_atoms=list(Structure[node]) # connected atoms.
    cnt=0
    for ni in con_atoms:
        if Structure.nodes[ni]['atom_name'][0]==atom:
            cnt+=1
    return cnt

def find_begin(Structure):
    """ Determine the begining carbon atom of a PEI structure.
    PEI only contains one carbon that is bound to 3 hydrogen and 1 carbon.

    Parameters
    ----------
    Structure: networkx Graph()
        AA Structure as graph.

    Returns
    -------
    Structure: networkx Graph()
        A new graph attribute 'start' is added.
    """
    for node in Structure:
        if Structure.nodes[node]['atom_name'][0]=='C':
            num_H=bound_atoms(Structure,node,'H')
            num_C=bound_atoms(Structure,node,'C')
            if num_H==3 and num_C==1:
                 Structure.graph['start']=node
                 break
    return Structure

def get_cgnode(Structure,node):
    """ Get the corresponding CG bead type for a given carbon atom C1, in 
    C1-C2-N1

    Parameters
    ----------
    Structure: networkx Graph()
        AA Structure as graph.
    node: node of Graph()
        node for which the CG bead type is determined.

    Returns
    -------
    cg_name: str
        Name of the bead 't', tq', 's', 'sq', 'p', 'pq'.
    next_C2: list of int
        IDs of C1 atoms in the next beads.
    mass: float
        Total mass of atoms in the bead
    charge: int
        Charge of bead
    bead_atoms: list of int
        IDs of atoms in the bead.
    """
    bead_atoms=[node] # AA atoms in the CG bead.
    next_nodes=list(Structure[node]) # Atoms connected to C1.
    for n in next_nodes: #Add all atoms connected to C1 except N
        if Structure.nodes[n]['atom_name'][0]!='N':
            bead_atoms.append(n)

    for n in next_nodes:
        if Structure.nodes[n]['atom_name'][0]=='C':
            next_C=n #atom C2.
            break

    next_nodes=list(Structure[next_C]) #Atoms connected to C2.
    for n in next_nodes:
        if n==node: #do not repeat node C1.
            continue
        bead_atoms.append(n) 

    for n in next_nodes:
        if Structure.nodes[n]['atom_name'][0]=='N':
            next_N=n #atom N1.
            break
    
    num_H=bound_atoms(Structure,next_N,'H') # get number of next hydrogens
    num_C=bound_atoms(Structure,next_N,'C') # get number of next carbons
    name='pst'
    charge=0
    cg_name=name[num_C-1] # num_C is degree of nitrogen.
    if num_C+num_H==4: # 4 covalent bonds implies protonation.
         cg_name+='q'
         charge=1

    next_C2=[] #C1 in next beads
    next_nodes=list(Structure[next_N])
    for n in next_nodes:
        if Structure.nodes[n]['atom_name'][0]!='C':
            bead_atoms.append(n) # add all hydrogen attached to N1.
    mass=0
    for ni in bead_atoms:
        mass+=Structure.nodes[ni]['mass']


    for n in next_nodes:
        if Structure.nodes[n]['atom_name'][0]=='C':
            if n != next_C: #not equal to C2 in current bead.
                next_C2.append(n)
    
    return cg_name, next_C2, mass, charge, bead_atoms

def AA2CG(Structure,ndx='aa2cg.ndx'):
    """ Converts AA structure to CG structure and writes the mapping scheme.

    Parameters
    ----------
    Structure: networkx Graph()
        AA Structure as graph.
    ndx: str
        File name to write the mapping scheme in gromacs .ndx format.

    Returns
    -------
    CG Structure: networkx DiGraph.
        Nodes have attributes:
            'bead_name': 't', 'tq', 's', 'sq', 'p', or 'pq'
            'mass': total mass of the bead.
            'charge': charge of the bead, 0 or 1.

    Writes
    ------
    ndx: aa2cg.ndx
        CG mapping scheme.
    """
    CGStruc=nx.DiGraph() #Since CG beads are assymetric, a DiGraph is used.
    mapping={} #map C1 atom IDs in C1-C2-N1 to CG atom IDs.
    cnt=1 #CG atom ID.
    queue=[Structure.graph['start']] #Starting carbon atom.
    w=open(ndx,'w') 
    beads={}

    while len(queue)>0:
        cur_node=queue.pop()
        cg_name,Cid,M,Q,bead_atoms=get_cgnode(Structure,cur_node)
        if cur_node not in CGStruc: #This is only run for the first bead.
            #w.write('[ Bead'+str(cnt)+' ] \n')
            CGStruc.add_node(cur_node,bead_name=cg_name,mass=M,charge=Q)
            mapping[cur_node]=cnt
            cnt+=1
        else:
            #w.write('[ Bead'+str(mapping[cur_node])+' ] \n')
            CGStruc.nodes[cur_node]['bead_name']=cg_name
            CGStruc.nodes[cur_node]['mass']=M
            CGStruc.nodes[cur_node]['charge']=Q
        for ni in Cid:
            CGStruc.add_edge(cur_node,ni)
            mapping[ni]=cnt
            cnt+=1
        beads[mapping[cur_node]]=sorted(bead_atoms)
        #for ns in bead_atoms:
        #    beads[mapping[cur_node]].append(ns
            #w.write(str(ns)+' ')
        #w.write('\n')

        queue+=list(Cid)#add new beads to queue
    for key in sorted(beads.keys()):
        w.write('[ Bead'+str(key)+' ]\n')
        for x in beads[key]:
            w.write(str(x)+' ')
        w.write('\n')
    w.close()

    return nx.relabel_nodes(CGStruc,mapping) #relabel the C1 IDs to CG bead IDs 1 to n.

def write_e2e(CGStruc):
    """ Write index file for calculating end-to-end distance of the polymer

    Parameters
    ----------
    CGStruc: networkx Graph
        CG chemical structure in network form.
    
    Writes
    ------
    e2e.ndx: A GROMACS .ndx file
 
    """
    large_len=0
    end_node=0
    for ni in CGStruc.nodes():
        if ni == 1:
            continue
        try:
            plen=nx.shortest_path_length(CGStruc,source=1,target=ni)
            if plen>large_len:
                large_len=plen
                end_node=ni
        except:
            continue
    w=open('e2e.ndx','w')
    w.write('[ PEI ]\n')
    w.write('1 ')
    cnt=1
    for ni in CGStruc.nodes():
        if ni==1 or ni==end_node:
            continue
        w.write(str(ni)+' ')
        cnt+=1
        if cnt%10==0:
            w.write('\n')
    w.write(str(end_node)+' ')
    w.write('\n')
    w.close()


def write_CGtopol(cg_idx,top_file,CGStruc,cgff_pickle,PEINAME='PEI'):
    """ Writes the Gromoacs topology file for CG-PEI.

    Parameters
    ----------
    cg_dir: str
        directory name of the CG iteration
    top_file: str
        filename of the CG topology file.
    CGStruc: networkx DiGraph()
        CG Structure 
    cgff_pickle: pickled dictionary
        Contains CG bonded parameters
    PEINAME: str
        Name to be given to the PEI.

    Writes
    ------
    [filename]: CG topology file. .itp
    cg_bonds.ndx: bonded pairs for generating probability distribution
    cg_angles.ndx: angle tuples for generating prob. dist.
    cg_dihs.ndx: dihedral tuples for generating prob. dist.
    unparam.pickle:  Unparameterized bonds/angles/dihedrals
    param.pickle: Parameterized bonds/angles/dihedrals 
    """
    #Writing .top file
    if cg_idx==-1:
        top=open(top_file,'w')
        itp=open(PEINAME+'.itp','w')
        top.write('#include "martini_v2.2refP.itp"\n' + 
                   '#include "martini_v2.0_ions.itp"\n' +
                   '#include "'+PEINAME+'.itp"\n\n[ system ]\n'+ 
                   ';name\n'+PEINAME[:4]+'\n\n[ molecules ]\n'+ 
                   ';name\tnumber\n'+PEINAME+'\t1\n')
        if type(cgff_pickle)==str:
            cgff=nx.read_gpickle(cgff_pickle)
        else:
            cgff=cgff_pickle
        itp.write('[ moleculetype ]\n; molname\tnrexcl\n '+str(PEINAME)+
            ' \t1\n\n[ atoms ]\n;id\t type\t resnr\t residue\t atom\t cgnr\t ' + 
            'charge\t Mass\n')
    else:
        top=open('CG'+str(cg_idx)+'/'+top_file,'w')
        itp=open('CG'+str(cg_idx)+'/'+PEINAME+str(cg_idx)+'.itp','w')
        top.write('#include "../martini_v2.2refP.itp"\n' + 
                  '#include "../martini_v2.0_ions.itp"\n' +
                  '#include "'+PEINAME+str(cg_idx)+'.itp"\n\n[ system ]\n'+ 
                  ';name\n'+PEINAME[:4]+'\n\n[ molecules ]\n'+ 
                  ';name\tnumber\n'+PEINAME+str(cg_idx)+'\t1\n')
        cgff=nx.read_gpickle('CG'+str(cg_idx)+'/cgff'+str(cg_idx)+'.pickle')
        itp.write('[ moleculetype ]\n; molname\tnrexcl\n '+str(PEINAME)+str(cg_idx)+
            ' \t1\n\n[ atoms ]\n;id\t type\t resnr\t residue\t atom\t cgnr\t ' + 
            'charge\t Mass\n')
    top.close() #Number of water molecules and salts will be written by other script.

    bead_type={'t': 'P1', 's': 'P1', 'p': 'P1','tq':'Qd', 'sq':'Qd', 'pq':'Qd'}
    bead_name={'t': 'Ter', 's': 'Sec', 'sq': 'QSec', 'p': 'Pri', 'pq': 'QPri'}

    cnt=1
    for n in CGStruc:
        CGStruc.nodes[n]['id']=cnt
        bname=CGStruc.nodes[n]['bead_name']

        itp.write(str(cnt)+'\t'+ #ID
            bead_type[bname]+'\t'+ #Type
            str(cnt)+'\tPEI\t'+ #Resnr Residue
            bead_name[bname]+'\t'+ #atom
            str(cnt)+'\t'+ #cgnr
            str(CGStruc.nodes[n]['charge'])+'\t'+ #charge
            str(round(CGStruc.nodes[n]['mass'],3))+'\n' #mass
            )
        cnt+=1

    bonds_byname={}
    angs_byname={}
    dihs_byname={}
    #Bonds
    bonds=find_bonds(CGStruc)
    itp.write('\n\n[ bonds ]\n;i\tj\tfunct\tlength\tforce.c\n')
    for i in range(len(bonds)):
        if bonds[i][2] not in bonds_byname:
            bonds_byname[bonds[i][2]]=[]
        bonds_byname[bonds[i][2]].append(bonds[i][:2])

        if bonds[i][2] not in cgff:
            raise Exception('Missing parameters'+bonds[i][2]+', generate new forcefield')
        else:
            itp.write(str(bonds[i][0])+'\t' + str(bonds[i][1])+'\t1\t' + \
                str(cgff[bonds[i][2]][0])+'\t'+str(cgff[bonds[i][2]][1]) + \
                '\t#'+str(bonds[i][2])+'\n')
    
    #Angles
    itp.write('\n\n[ angles ]\n;i\tj\tk\tfunct\tangle\tforce.c\n')
    angles=find_angles(CGStruc)
    for i in range(len(angles)):
        if angles[i][3] not in angs_byname:
            angs_byname[angles[i][3]]=[]
        angs_byname[angles[i][3]].append(angles[i][:3])

        if angles[i][3] not in cgff:
            raise Exception('Missing parameters for '+angles[i][3]+', generate new forcefield')
        else:
            itp.write(str(angles[i][0])+'\t' + str(angles[i][1])+'\t' + \
                str(angles[i][2])+'\t2\t' + str(cgff[angles[i][3]][0])+'\t' + \
                str(cgff[angles[i][3]][1])+'\t#'+str(angles[i][3])+'\n')
    
    #Dihedrals
    dihs=find_dihs(CGStruc)
    itp.write('\n\n[ dihedrals ]\n;i\tj\tk\tl\tfunc\tangle\tforce.c\tmul\n')
    for i in range(len(dihs)):
        if dihs[i][4] not in dihs_byname:
            dihs_byname[dihs[i][4]]=[]
        dihs_byname[dihs[i][4]].append(dihs[i][:4])

        itp.write(';'+str(dihs[i][4])+'\n')
        if dihs[i][4] not in cgff:
            raise Exception('Missing parameters for '+dihs[i][4]+', generate new forcefield')
        else:
            iters=cgff[dihs[i][4]][0]
            for j in range(iters):
                itp.write(str(dihs[i][0])+'\t' + str(dihs[i][1])+'\t' + \
                    str(dihs[i][2])+'\t' + str(dihs[i][3])+'\t1\t' + \
                    str(cgff[dihs[i][4]][1][j][0])+'\t' + 
                    str(cgff[dihs[i][4]][1][j][1])+'\t' + 
                    str(cgff[dihs[i][4]][1][j][2])+'\n')
            itp.write('\n')
    itp.close() 
        
def get_dist(filename):
    """ Read probability distribution generated by Gromacs in .xvg format

    Parameters
    ----------
    filename: str
        filename of .xvg file containting probability distribution.
    
    Returns
    -------
    data: index 0 is x, 1 is y
    """
    f=open(filename,'r')
    data=[]
    for lines in f:
        foo=lines.split()
        if len(foo)==0:
            continue
        if foo[0][0]=='#' or foo[0][0]=='@':
            continue
        foo=[float(x) for x in foo]
        data.append(foo)
    f.close()
    data=np.array(data)
    return data

def gen_unparam(CGStruc,cgff_cur_pickle,dih_initial=None):
    """ Categorizes bonded distributions to parameterized and unparameterized 
    
    Parameters
    ----------
    CGStruc: networkx DiGraph()
        CG Structure 
    cgff_cur_pickle: pickled dictionary
        Contains CG bonded parameters

    Writes
    ------
    "unparam.pickle" : Pickled dictionary containing names of unparameterized
                       bonded parameters.
    "param.pickle" : Pickled dictionary containing names of parameterized
                     bonded parameters.
    "cg_bonds.ndx" : GROMACS .ndx file containing grouped CG bonds.
    "cg_angs.ndx" : GROMACS .ndx file containing grouped CG angles.
    "cg_dihs.ndx" : GROMACS .ndx file containing grouped CG dihedrals.
    """
    if cgff_cur_pickle is None:
        cgff_cur={}
    elif type(cgff_cur_pickle)==str:
        cgff_cur=nx.read_gpickle(cgff_cur_pickle)
    else:
        cgff_cur=cgff_cur_pickle
    
    bonds=find_bonds(CGStruc)
    angs=find_angles(CGStruc)
    dihs=find_dihs(CGStruc)
    unparam={'bonds':{},'angs':{},'dihs':{}}
    param={'bonds':[],'angs':[],'dihs':[]}
    bonds_byname={}
    angs_byname={}
    dihs_byname={}
    for data in bonds:
        if data[2] not in bonds_byname:
            bonds_byname[data[2]]=[]
        bonds_byname[data[2]].append(data[:2])

        if data[2] not in cgff_cur:
            if data[2] not in unparam['bonds']:
                unparam['bonds'][data[2]]=1.0
        else:
            if data[2] not in param['bonds']:
                param['bonds'].append(data[2])
    for data in angs:
        if data[3] not in angs_byname:
            angs_byname[data[3]]=[]
        angs_byname[data[3]].append(data[:3])

        if data[3] not in cgff_cur:
            if data[3] not in unparam['angs']:
                unparam['angs'][data[3]]=1.0
        else:
            if data[3] not in param['angs']:
                param['angs'].append(data[3])
    for data in dihs:
        if data[4] not in dihs_byname:
            dihs_byname[data[4]]=[]
        dihs_byname[data[4]].append(data[:4])

        if data[4] not in cgff_cur:
            if data[4] not in unparam['dihs']:
                unparam['dihs'][data[4]]=[8,1.0]
        else:
            if data[4] not in param['dihs']:
                param['dihs'].append(data[4])

    if dih_initial is not None:
        f=open(dih_initial,'r')
        for lines in f:
            foo=lines.split()
            unparam['dihs'][foo[0]]=[int(foo[2]),1.0]

    nx.write_gpickle(unparam,'unparam.pickle')
    nx.write_gpickle(param,'param.pickle')

    ##Write mapping for bonds ###
    w=open('cg_bonds.ndx','w')
    for name in bonds_byname:
        w.write('[ '+name+' ]\n')
        bonds_byname[name].reverse()
        for data in bonds_byname[name]:
            for i in range(2):
                w.write(str(data[i])+' ')
            w.write('\n')
        w.write('\n')
    w.close()
    ##Write mapping for angs ###
    w=open('cg_angs.ndx','w')
    for name in angs_byname:
        w.write('[ '+name+' ]\n')
        angs_byname[name].reverse()
        for data in angs_byname[name]:
            for i in range(3):
                w.write(str(data[i])+' ')
            w.write('\n')
        w.write('\n')
    w.close()
    ##Write mapping for dihs ###
    w=open('cg_dihs.ndx','w')
    for name in dihs_byname:
        w.write('[ '+name+' ]\n')
        dihs_byname[name].reverse()
        for data in dihs_byname[name]:
            for i in range(4):
                w.write(str(data[i])+' ')
            w.write('\n')
        w.write('\n')
    w.close()

def gen_new_cgff(CGStruc,cgff_cur_pickle,cgff_out,Temp=300):
    """ Generates initial guess for missing bonded parameters.
    It is only used for the first iteration.
    
    Parameters
    ----------
    CGStruc: networkx DiGraph()
        CG Structure 
    cgff_cur_pickle: pickled dictionary/str
        Contains CG bonded parameters. If datatype is str, it refers to the name
        of the pickled file.
    cgff_out: pickled dictionary
        Contains all CG bonded parameters needed for simulation.

    Writes
    ------
    [cgff_out]: a pickled dictionary containing bonded parameters.
    """
    #Write code to gen new ff.
    if cgff_cur_pickle is None:
        cgff_cur={}
    elif type(cgff_cur_pickle)==str:
        cgff_cur=nx.read_gpickle(cgff_cur_pickle)
    else:
        cgff_cur=cgff_cur_pickle

    cgff_new=gen_bond_params(cgff_cur,cgff_pickle_out=None,Temp=Temp)
    cgff_new=gen_ang_params(cgff_cur,cgff_pickle_out=None,Temp=Temp)
    cgff_new=gen_dih_params(cgff_new,cgff_pickle_out=None,Temp=Temp)
    nx.write_gpickle(cgff_new,cgff_out)

#To generate initial dihedral parameters from reference AA distributions.
# Procedure:
# 1) Reference potential energy of a dihedral angle potential (V_{dih}^AA) is calulated from the dihedral angle distribution P_{dih}^{AA} using the equation
#    V_{dih}^AA=-k_B*T*ln(P_{dih}^{AA})
# 2) V_{dih}^AA is fit using fourier series of 'n' terms  a_0 + \sum_{i=1}^{i=n} ( a_i \cos(i*x) + b_i \sin(i*x) )
# 4) parameters a_i and b_i (for i= 1 to n) are used to obtain initial parameters for CG dihedral angle potential 
#    V_{dih}^{CG} (\phi_{ijkl})=  \sum_{w=1}^{w=2n} K_{d,w} \[ 1 + \cos(n_w*\phi_{ijkl} - \phi_{eq,w}) \]
#    where, n_w = greatest integer function of (w+1)/2 and \phi_{eq,w}= 0 or 180 when w is odd and -90 or 90 when w is even.
#    using the relations 
#    K_{d,2*i-1} = abs(a_i), \phi_{eq,2*i-1}= 90 - sign(a_i)*90
#    K_{d,2*i} = abs(b_i), \phi_{eq,2*i}= sign(b_i)*90
#    Note: If Fourier series with 'n' parameters are used, we would get parameters for dihedral angle potential for sum of '2*n' periodic potentials.

def mean_pdf(data):
    """ Calculates mean of a probability distribution function

    Parameters
    ----------
    data: 2D numpy array
        data[:,0] is bond lenght/bond angle/dihedral angle, data[:,1] is probability.

    Returns
    -------
    mean of the distribution
    """
    return np.sum(np.multiply(data[:,0],data[:,1]))/np.sum(data[:,1])

def var_pdf(data):
    """ Calculates variance of a probability distribution function

    Parameters
    ----------
    data: 2D numpy array
        data[:,0] is bond lenght/bond angle/dihedral angle, data[:,1] is probability.

    Returns
    -------
    variance of the distribution
    """
    return np.sum(np.multiply(np.square(data[:,0]),data[:,1]))/np.sum(data[:,1])-mean_pdf(data)**2

def std_pdf(data):
    """ Calculates standard deviation of a probability distribution function

    Parameters
    ----------
    data: 2D numpy array
        data[:,0] is bond lenght/bond angle/dihedral angle, data[:,1] is probability.

    Returns
    -------
    standard deviation of the distribution
    """
    return np.sqrt(var_pdf(data))

def dih_name_to_param(dih_name,cgff):
    """ Returns Kd, phi and n of the dihedral angle for a given name

    Parameters
    ----------
    dih_name: str
        name of the dihedral angle
    cgff: dictionary
        CG forcefield dictionary

    Returns
    -------
    Kd: force constants
    phi: phase angles
    1n_dih: multiplicities
    """
    foo=cgff[dih_name][1]
    Kd=np.zeros(cgff[dih_name][0])
    phi=np.zeros(cgff[dih_name][0])
    n_dih=np.zeros(cgff[dih_name][0])
    for i in range(len(Kd)):
        phi[i]=foo[i][0]
        Kd[i]=foo[i][1]
        n_dih[i]=foo[i][2]
    return Kd,phi,n_dih

def dih_param_to_fourier(Kd,phi):
    """ Convert dihedral parameters Kd, phi (assuming n is 1, 1, 2, 2, ... ) to
        Fourier series coefficients.

    Parameters
    ----------
    Kd: array
        Force constants for the dihedral angle potential
    phi: array
        Phase constants for the dihedral angle potential

    Returns
    ------
    popt: Coefficients of the Fourier series. The first coefficient is the 
          constant term. 
    """
    popt=np.zeros(len(Kd)+1)
    for i in range(len(Kd)):
        if i%2==0:
            if phi[i]==0:
                popt[i+1]=Kd[i]
            else:
                popt[i+1]=-Kd[i]
        if i%2==1:
            if phi[i]==90:
                popt[i+1]=Kd[i]
            else:
                popt[i+1]=-Kd[i]
    return popt

def fourier_series(x,*a):
    """Generate fourier series.
    Typically used for dihedral angles.

    Parameters
    ---------
    x: array
        Typically this value is dihedral angle.
    *a: array
        coefficients of fourier series
    Returns
    -------
    foo: an array containing the fourier series value
    """
    foo=0
    i=0
    for args in a:
        if i==0:
            foo+=args
        elif i%2==1:
            foo+=args*np.cos(np.pi*(int((i+1)/2)*x)/180.0)
        else:
            foo+=args*np.sin(np.pi*(int(i/2)*x)/180.0)
        i+=1
    return foo

def bond_fit(x,a0,a1,a2):
    """ bond angle distribution used to fit data.
    Angle potential is = k/2(r-r0)^2 
                       = k/2(r^2  - 2*r*r0 + r0**2)
    th0: equilibrium bond angle parameter
    k: force constant for bond angle distribution

    Parameters
    ----------
    x: array
        bond angle
    a0: float
        0.5*k*r0**2 +C
    a1: float
        -k*r0
    a2: float
        k/2
    Returns
    -------
        Returns anlge potential for given a0,a1,a2
    """
    return a0+a1*x+a2*np.square(x)

def ang_fit(x,a0,a1,a2):
    """ bond angle distribution used to fit data.
    Angle potential is = k/2*(cos(th)-cos(th0))^2 
                       = k/2(1/2+cos^2(th0)) - kcos(th0)cos(th) + k/4cos(2*th)
    th0: equilibrium bond angle parameter
    k: force constant for bond angle distribution

    Parameters
    ----------
    x: array
        bond angle
    a0: float
        0.5*k*(0.5+cos^2(th0))+C
    a1: float
        -k*cos(th0)
    a2: float
        k/4
    Returns
    -------
        Returns anlge potential for given a0,a1,a2
    """
    return a0+a1*np.cos(np.pi*x/180.0)+a2*np.cos(2*np.pi*x/180.0)

def conv_dih_par(a):
    """ Convert dihedral angle parameters in Fourier series format to CG 
        forcefield format

    Parameters
    ----------
    a: array
        Fourier series coefficients
    Returns
    -------
    data: Dihedral angle in forcefield format.
          data[0] is number of periodic functions
          data[1] contains "data[0]" arrays [Phi_eq, Kd, n].
    """
    n=len(a)
    params=[]
    for i in range(1,n):
        foo=[]
        if i%2==1:
            foo.append(int(90-np.sign(a[i])*90)) #Phi eq
        else:
            foo.append(int(np.sign(a[i])*90)) #Phi eq
        foo.append(round(abs(a[i]),5)) #Kd
        foo.append(int((i-1)/2.0+0.005)+1)
        params.append(foo)
    return [n-1,params]

def conv_bond_par(a):
    """ Convert bond potential fit parameters a0, a1, a2 to force constant k
        and equilibrium bond length r0

    k=a2*2
    r0=-a1/k

    Parameters
    ----------
    a: array
        array containing a0, a1, a2
    Returns
    -------
    Kb: Force constant
    r0: Equilibrium angle 
    """
    Kb=round(a[2]*2,0) #K/2=a[2]
    r0=round(-a[1]/Kb,3) #ro
    return Kb, r0

def conv_ang_par(a):
    """ Convert angle potential fit parameters a0, a1, a2 to force constant k
        and equilibrium angle th0

    k=a2*4
    th0=acos(0.25*a1/a2)

    Parameters
    ----------
    a: array
        array containing a0, a1, a2
    Returns
    -------
    Ka: Force constant
    th0: Equilibrium angle 
    """
    Ka=round(a[2]*4,0) #K/4=a[2]
    foo=-a[1]/Ka #cos(th0)
    if foo>1:#Cos cant be more than 1
        th0=0
    elif foo<-1: #Cos cant be less than -1
        th0=180
    else:#Get angle tho0
        th0=round(np.arccos(foo)*180.0/np.pi,0)
    return Ka,th0

def ucg_dih_ref_dist(x,params):
    """ Generate CG dihedral potential energy distribution from parameters

    Parameters
    ----------
    x: array of floats
        Dihedral angle
    params: 2D array 
        Dihedral angle parameters in forcefield format.

    Returns
    -------
        Dihedral angle potential energy distribution
    """
    u=np.zeros(len(x))
    for i in range(params[0]):
        u=np.add(u,params[1][i][1]*(1+np.cos(np.pi*(params[1][i][2]*x-params[1][i][0])/180.0)))
    return u

def ucg_bond_ref_dist(x,ka,r0):
    """ Generate CG angle potential energy distribution from parameters

    Parameters
    ----------
    x: array of floats
        bond angle 
    ka: float
        Force constant
    r0: float
        Equilibrium bond angle
    Returns
    -------
        Bond angle potential energy distribution
    """
    u=np.zeros(len(x))
    u=0.5*ka*np.square(x-r0)
    return u

def ucg_ang_ref_dist(x,ka,th0):
    """ Generate CG angle potential energy distribution from parameters

    Parameters
    ----------
    x: array of floats
        bond angle 
    ka: float
        Force constant
    th0: float
        Equilibrium bond angle
    Returns
    -------
        Bond angle potential energy distribution
    """
    u=np.zeros(len(x))
    u=0.5*ka*np.square(np.cos(np.pi*x/180.0)-np.cos(np.pi*th0/180.0))
    return u

def fill_data(data,typ='dih'): 
    """ Fill the dihedral angle data to ensure dihedral angle is between 
        -180 to 180 both included

    Parameters
    ----------
    data: 2D array
        data[:,0] is dihedral angle, data[:,1] is the probability distribution
    Returns
    -------
    foo: 2D array with complete dihedral angle range.
    """
    if typ=='dih':
        min_val=-180
    elif typ=='ang':
        min_val=0
    max_val=180

    if data[0,0]==-min_val and data[-1,0]==max_val:
        return data
    else:
        y=data[:,1].tolist()
        y=[0]*int(data[0,0]-min_val)+y+[0]*int(max_val-data[-1,0])
        foo=np.zeros((max_val-min_val+1,2))
        foo[:,0]=np.linspace(min_val,max_val,max_val-min_val+1)
        foo[:,1]=np.array(y)
        return foo

def gen_dih_params(cgff_pickle,cgff_pickle_out=None,Temp=300):
    """ Generate initial guess for dihedral angle parameters

    Parameters
    ----------
    unparam_pickle: dictionary
        Dictionary containing unparameterized bonded distribution names.
    cgff_pickle: dictionary/string
        Initial forcefield pickle. Can be a dictionary or a string. The string 
        is the name of the pickled file that stores the dictionary.
    cgff_pickle_out: dictionary/string (Optional)
        Returns dictionary or if the value is a string writes a pickled file. 
        (Default None)
    Returns
    -------
    cgff: dictionary
        if cgff_pickle_out is None (default).
    """
    unparam=nx.read_gpickle('unparam.pickle')
    kbt=2.479*Temp/298.0 #kJ/mol  

    if type(cgff_pickle)==str:
        cgff=nx.read_gpickle(cgff_pickle)
    else:
        cgff=cgff_pickle

    for name in unparam['dihs']:
        data_aa=get_dist('bonded_distribution/dih_'+name+'.xvg') #Get AA dihedral angle distribution
        data_aa=fill_data(data_aa) #make the data complete that is x belongs to (-180,180)
        Uaa=-kbt*np.log(data_aa[:,1]+tol) #Get potential energy of bonded distribution
        p0=[max(Uaa)]*(unparam['dihs'][name][0]+1)
        popt,pcov=curve_fit(fourier_series,data_aa[:,0],Uaa,p0) #Fit dihedral angle data with Fourier series 
        params=conv_dih_par(popt) #Convert Fourier series parameters to dihedral angle force field parameters
        cgff[name]=params
        #Check the CG and AA bonded distribution
        Ucg=ucg_dih_ref_dist(data_aa[:,0],params)
        Ucg=np.array(Ucg)
        Pcg=np.exp(-Ucg/kbt)
        Pcg=Pcg/np.sum(Pcg)
    
        plt.figure
        plt.plot(data_aa[:,0],data_aa[:,1],'k')
        plt.plot(data_aa[:,0],Pcg,'r')
        plt.savefig('ref_param_fit/dih_'+name+'.png',dpi=300)
        plt.clf()
        #print(name+' done')
    if cgff_pickle_out is None:
        return cgff
    elif type(cgff_pickle_out)==str:
        nx.write_gpickle(cgff,cgff_pickle_out)

def gen_bond_params(cgff_pickle,cgff_pickle_out=None,Temp=300):
    """ Generate initial guess for bond angle parameters

    Parameters
    ----------
    unparam_pickle: dictionary
        Dictionary containing unparameterized bonded distribution names.
    cgff_pickle: dictionary/string
        Initial forcefield pickle. Can be a dictionary or a string. The string 
        is the name of the pickled file that stores the dictionary.
    cgff_pickle_out: dictionary/string (Optional)
        Returns dictionary or if the value is a string writes a pickled file. 
        (Default None)
    Returns
    -------
    cgff: dictionary
        if cgff_pickle_out is None (default).
    """
    kbt=2.479*Temp/298.0 #kJ/mol  
    unparam=nx.read_gpickle('unparam.pickle')

    if type(cgff_pickle)==str:
        cgff=nx.read_gpickle(cgff_pickle)
    else:
        cgff=cgff_pickle

    for name in unparam['bonds']:
        data_aa=get_dist('bonded_distribution/bond_'+name+'.xvg') #Read AA angle distributions
        for i in range(len(data_aa)):
            if data_aa[i,1]>0:
                break
        for j in range(len(data_aa)-1,-1,-1):
            if data_aa[j,1]>0:
                break
        data_aa=data_aa[max(0,i-10):min(len(data_aa),j+11),:]
                
        Uaa=-kbt*np.log(data_aa[:,1]+tol) #Get the potential energy of angle distribution from Boltzmann inversion
        p0=[225,750,2500]
        popt,pcov=curve_fit(bond_fit,data_aa[:,0],Uaa,p0)  #Generate angle parameters from curve fitting
        Kb,r0=conv_bond_par(popt) #Generate angle parameters for force field
        cgff[name]=[r0,Kb]
        #Plot the fit to a png file
        Ucg=ucg_bond_ref_dist(data_aa[:,0],Kb,r0)
        dx=0.5*(data_aa[2,0]-data_aa[0,0])
        Ucg=np.array(Ucg)
        Pcg=np.exp(-Ucg/kbt)
        Pcg=Pcg/(np.sum(Pcg)*dx)

    
        plt.figure
        plt.plot(data_aa[:,0],data_aa[:,1],'k')

def gen_ang_params(cgff_pickle,cgff_pickle_out=None,Temp=300):
    """ Generate initial guess for bond angle parameters

    Parameters
    ----------
    unparam_pickle: dictionary
        Dictionary containing unparameterized bonded distribution names.
    cgff_pickle: dictionary/string
        Initial forcefield pickle. Can be a dictionary or a string. The string 
        is the name of the pickled file that stores the dictionary.
    cgff_pickle_out: dictionary/string (Optional)
        Returns dictionary or if the value is a string writes a pickled file. 
        (Default None)
    Returns
    -------
    cgff: dictionary
        if cgff_pickle_out is None (default).
    """
    kbt=2.479*Temp/298.0 #kJ/mol  
    unparam=nx.read_gpickle('unparam.pickle')

    if type(cgff_pickle)==str:
        cgff=nx.read_gpickle(cgff_pickle)
    else:
        cgff=cgff_pickle

    for name in unparam['angs']:
        data_aa=get_dist('bonded_distribution/angle_'+name+'.xvg') #Read AA angle distributions
        Uaa=-kbt*np.log(data_aa[:,1]+tol) #Get the potential energy of angle distribution from Boltzmann inversion
        p0=[900,200,250]
        popt,pcov=curve_fit(ang_fit,data_aa[:,0],Uaa,p0)  #Generate angle parameters from curve fitting
        Ka,th=conv_ang_par(popt) #Generate angle parameters for force field
        if th<50:
            th=50
        cgff[name]=[th,Ka]
        #Plot the fit to a png file
        Ucg=ucg_ang_ref_dist(data_aa[:,0],Ka,th)
        Ucg=np.array(Ucg)
        Pcg=np.exp(-Ucg/kbt)
        Pcg=Pcg/np.sum(Pcg)
    
        plt.figure
        plt.plot(data_aa[:,0],data_aa[:,1],'k')
        plt.plot(data_aa[:,0],Pcg,'r')
        plt.savefig('ref_param_fit/angle_'+name+'.png',dpi=300)
        plt.clf()
        #print(name+' done')
    if cgff_pickle_out is None:
        return cgff
    elif type(cgff_pickle_out)==str:
        nx.write_gpickle(cgff,cgff_pickle_out)

def mean_pdf(data):
    """ Mean of a probability density function.

    Parameters
    ----------
    data: 2D numpy array
         probability density function
    
    Returns
    -------
    mean of the probability density function (float)
    """
    return np.sum(np.multiply(data[:,0],data[:,1]))/np.sum(data[:,1])

def var_pdf(data):
    """ Variance of a probability density function.

    Parameters
    ----------
    data: 2D numpy array
         probability density function
    
    Returns
    -------
    Variance of the probability density function (float)
    """
    return np.sum(np.multiply(np.square(data[:,0]),data[:,1]))/np.sum(data[:,1])-mean_pdf(data)**2

def std_pdf(data):
    """ Standard deviation of a probability density function.

    Parameters
    ----------
    data: 2D numpy array
         probability density function
    
    Returns
    -------
    Standard deviation of the probability density function (float)
    """
    return np.sqrt(var_pdf(data))

def dih_name_to_param(dih_name,cgff):
    """ Finds the parameters associated with the dihedral name

    Parameters
    ----------
    dih_name: str
        Name of the dihedral angle
    cgff: list
        CG forcefield parameters

    Returns
    -------
    Dihedral angle parameters K, phi, n
    """

    foo=cgff[dih_name][1]
    Kd=np.zeros(cgff[dih_name][0])
    phi=np.zeros(cgff[dih_name][0])
    n_dih=np.zeros(cgff[dih_name][0])
    for i in range(len(Kd)):
        phi[i]=foo[i][0]
        Kd[i]=foo[i][1]
        n_dih[i]=foo[i][2]
    return Kd,phi,n_dih

def dih_param_to_fourier(Kd,phi):
    """ Convert dihedral angle parameters to Fourier series form.

    Parameters
    ----------
    Kd: array of float
        Force constants of dihedral angle
    phi: array of integers
        Angle phis (0, -90, 90, 180) of dihedral angle

    Returns
    -------
    Array of fourier series parameters.
    """

    popt=np.zeros(len(Kd)+1)
    for i in range(len(Kd)):
        if i%2==0:
            if phi[i]==0:
                popt[i+1]=Kd[i]
            else:
                popt[i+1]=-Kd[i]
        if i%2==1:
            if phi[i]==90:
                popt[i+1]=Kd[i]
            else:
                popt[i+1]=-Kd[i]
    return popt

def update_bonded_params(fb, fa, fd, wb, wa, aa_dir, idx, Temp=300, cost_tol=1E-8):
    """ Update bonded parameters for next iteration

    Since all bond length parameters have been determined, they are not 
    updated.

    Parameters
    ----------
    fa: float
        Gradient descent step parameter for bond angle distributions.
    fd: float
        Gradient descent step parameter for dihedral angle distributions.
    wb: float
        Weight [0,1] that assigns relative importance to match mean and 
        standard deviation of CG bond length distributions to that of 
        reference CG distribution.
    wa: float
        Weight [0,1] that assigns relative importance to match mean and 
        standard deviation of CG bond angle distributions to that of 
        reference CG distribution.
    aa_dir: str
        Path to the directory containing directors CG*
    idx: int
        Index of the new iteration.
    cost_tol: float
        Cost tolerance below which parameters are not changed
        
    Returns
    -------
    change: (int) 
        0 if some parameters have changed
       -1 if no parameters have changed.
    """
    kbt=2.479*Temp/298.0 # kJ/mol 
    unparam=nx.read_gpickle(aa_dir + '/unparam.pickle')
    err_stat=nx.read_gpickle(aa_dir + '/Cost.pickle')
    ff={}
    for x in ["", "_th", "_K"]:
        try:
            ff['CG'+x]=nx.read_gpickle(aa_dir+'/CG'+str(idx-1)+x+'/cgff'+str(idx-1)+x+'.pickle')
        except:
            pass

    #cost=err_stat[idx-1]['cost']
    #if idx>2:
    #    old_cost=err_stat[idx-2]['cost']

    change=False
    for bond in unparam['bonds']:
        if err_stat[idx-1][bond]<cost_tol: #cost[bond]<cost_tol:
            continue

        if idx>2:
            if err_stat[idx-1][bond]>err_stat[idx-2][bond]: #cost[bond]>old_cost[bond]:
                unparam['bonds'][bond]*=0.9
        data={}
        mean={}
        std={}

        data['AA']=get_dist(aa_dir+'/bonded_distribution/bond_'+bond+'.xvg')
        mean['AA']=mean_pdf(data['AA'])
        std['AA']=std_pdf(data['AA'])
        foo_idx=idx-1

        for x in ["", "_th", "_K"]:
            data['CG'+x]=get_dist(aa_dir+'/CG'+str(foo_idx)+x+'/bonded_distribution/bond_'+bond+'.xvg')
            mean['CG'+x]=mean_pdf(data['CG'+x])
            std['CG'+x]=std_pdf(data['CG'+x])
                
        dr=ff['CG_th'][bond][0]-ff['CG'][bond][0]
        dK=ff['CG_K'][bond][1]-ff['CG'][bond][1]

        if dr<0.001:
            dmudr=0
            dsigdr=0
        else:
            dmudr=(mean['CG_th']-mean['CG'])/dr
            dsigdr=(std['CG_th']-std['CG'])/dr

        if dK<0.1:
            dmudK=0
            dsigdK=0
        else:
            dmudK=(mean['CG_K']-mean['CG'])/dK
            dsigdK=(std['CG_K']-std['CG'])/dK

        d_param=[0,0]
        cnt=0
        if dr>0:
            d_param[0]-= unparam['bonds'][bond]*fb*2*(wb*(mean['CG'] - 
                mean['AA'])*dmudr + (1-wb)*(std['CG'] - std['AA'])*dsigdr)
            cnt+=1
        if dK>0:
            d_param[1]-= unparam['bonds'][bond]*fb*2*(wb*(mean['CG'] - 
                mean['AA'])*dmudK + (1-wb)*(std['CG'] - std['AA'])*dsigdK)*100.0
            cnt+=1

        if abs(d_param[0])>0.05:
            d_param[0]=0.05*np.sign(d_param[0])
        if abs(d_param[1])>500:
            d_param[1]=500.*np.sign(d_param[1])

        if cnt>0:
            ff['CG'][bond][0]=round(d_param[0]+ff['CG'][bond][0],3)
            ff['CG'][bond][1]=round(min(max(d_param[1]+ff['CG'][bond][1],10),40000),1)
            change=True
        
    for ang in unparam['angs']:
        if err_stat[idx-1][ang]<cost_tol: #cost[ang]<cost_tol:
            continue

        if idx>2:
            if err_stat[idx-1][ang]>err_stat[idx-2][ang]: #cost[ang]>old_cost[ang]:
                unparam['angs'][ang]*=0.9
        data={}
        mean={}
        std={}

        data['AA']=get_dist(aa_dir+'/bonded_distribution/angle_'+ang+'.xvg')
        mean['AA']=mean_pdf(data['AA'])
        std['AA']=std_pdf(data['AA'])
        foo_idx=idx-1

        for x in ["", "_th", "_K"]:
            data['CG'+x]=get_dist(aa_dir+'/CG'+str(foo_idx)+x+'/bonded_distribution/angle_'+ang+'.xvg')
            mean['CG'+x]=mean_pdf(data['CG'+x])
            std['CG'+x]=std_pdf(data['CG'+x])
                
        dth=ff['CG_th'][ang][0]-ff['CG'][ang][0]
        dK=ff['CG_K'][ang][1]-ff['CG'][ang][1]

        if dth<0.1:
            dmudth=0
            dsigdth=0
        else:
            dmudth=(mean['CG_th']-mean['CG'])/dth
            dsigdth=(std['CG_th']-std['CG'])/dth

        if dK<0.1:
            dmudK=0
            dsigdK=0
        else:
            dmudK=(mean['CG_K']-mean['CG'])/dK
            dsigdK=(std['CG_K']-std['CG'])/dK

        d_param=[0,0]
        cnt = 0
        if dth>0:
            d_param[0]-= unparam['angs'][ang]*fa*2*(wa*(mean['CG'] - 
                mean['AA'])*dmudth + (1-wa)*(std['CG'] - std['AA'])*dsigdth)
            cnt+=1
        if dK>0:
            d_param[1]-= unparam['angs'][ang]*fa*2*(wa*(mean['CG'] - 
                mean['AA'])*dmudK + (1-wa)*(std['CG'] - std['AA'])*dsigdK)*100.0
            cnt+=1

        if abs(d_param[0])>10:
            d_param[0]=10.*np.sign(d_param[0])
        if abs(d_param[1])>200:
            d_param[1]=200.*np.sign(d_param[1])

        if cnt>0:
            ff['CG'][ang][0]=round(min(max(d_param[0]+ff['CG'][ang][0],30),180),1)
            ff['CG'][ang][1]=round(min(max(d_param[1]+ff['CG'][ang][1],10),10000),1)
            change=True
        
    for dih in unparam['dihs']:
        if err_stat[idx-1][dih]<cost_tol: #cost[dih]<cost_tol:
            continue

        if idx>2:
            if err_stat[idx-1][dih]>err_stat[idx-2][dih]: #cost[dih]>old_cost[dih]:
                unparam['dihs'][dih][1]*=0.9

        #Get bonded distritions
        data={}
        U={} #Get potential energy associated with bonded distributions 
        data['AA']=get_dist(aa_dir+'/bonded_distribution/dih_'+dih+'.xvg')
        data['AA']=fill_data(data['AA'])
        U['AA']=-kbt*np.log(np.array(data['AA'][:,1])+tol) 
        data['CG']=get_dist(aa_dir+'/CG'+str(idx-1)+'/bonded_distribution/dih_'+dih+'.xvg')
        data['CG']=fill_data(data['CG'])
        U['CG']=-kbt*np.log(np.array(data['CG'][:,1])+tol) 
        dU=np.subtract(U['AA'],U['CG']) # {e0}
        p0=[1]*(ff['CG'][dih][0]+1)
        popt,pcov=curve_fit(fourier_series,data['AA'][:,0],dU,p0)
            
        Kd_old,phi_old,n_old=dih_name_to_param(dih,ff['CG'])
        param_old=dih_param_to_fourier(Kd_old,phi_old)
        max_popt=max(np.abs(popt[1:]))
        if max_popt*fd*unparam['dihs'][dih][1]>0.5:
            popt=popt*0.5/max_popt
        if max_popt*fd*unparam['dihs'][dih][1]>1E-5:
            param_new=np.add(param_old,fd*popt*unparam['dihs'][dih][1])
            ff['CG'][dih]=conv_dih_par(param_new)
            change=True

    if change:
        nx.write_gpickle(ff['CG'],aa_dir+'/CG'+str(idx)+'/cgff'+str(idx)+'.pickle')
        nx.write_gpickle(unparam,'unparam.pickle')
        return 0
    else:
        return -1

def update_th_params(rc,thc,idx,aa_dir,cost_tol=1E-6):
    """ Update only theta parameters of bond angle distribution for determining 
        Jacobian

    Parameters
    ----------
    thc: float
        Theta parameters are increased by thc. If current theta parameter is 
        >= 179, new theta parameters is 180-thc. New theta parameter is >= 30.
    idx: int
        Index of current iteration.
    aa_dir: str
        Path to the directory containting CG* directories.
    cost_tol: float
        Cost tolerance below which theta parameters are not changed.

    Returns
    -------
    changed: int
        0 if some parameters are changed
       -1 if no parameters are changed
   
    """
    unparam=nx.read_gpickle('unparam.pickle')
    cgff=nx.read_gpickle(aa_dir+'/CG'+str(idx)+'/cgff'+str(idx)+'.pickle')
    err_stat=nx.read_gpickle(aa_dir + '/Cost.pickle')
    changed=False
    for bond in unparam['bonds']:
        if err_stat[idx][bond]<cost_tol:
            continue
        else:
            cgff[bond][0]=cgff[bond][0]+rc
        changed=True

    for ang in unparam['angs']:
        if err_stat[idx][ang]<cost_tol:
            continue
        if cgff[ang][0]>=179:
            cgff[ang][0]=180-thc
        else:
            cgff[ang][0]=min(180,max(cgff[ang][0]+thc,30))
        changed=True

    if changed:
        nx.write_gpickle(cgff,aa_dir+'/CG'+str(idx)+'_th/cgff'+str(idx)+'_th.pickle')
        return 0
    else:
        return -1

def update_K_params(Kbc,Kac,idx,aa_dir,cost_tol=1E-6):
    """ Update only force constant parameters of bond angle distribution for 
        determining Jacobian

    Parameters
    ----------
    kc: float
        Force constant parameter are increased by kc. New force constant
        is >= 10 kJ/mol. 
    idx: int
        Index of current iteration.
    aa_dir: str
        Path to the directory containting CG* directories.
    cost_tol: float
        Cost tolerance below which theta parameters are not changed.

    Returns
    -------
    changed: int
        0 if some parameters are changed
       -1 if no parameters are changed
   
    """
    unparam=nx.read_gpickle(aa_dir + '/unparam.pickle')
    cgff=nx.read_gpickle(aa_dir + '/CG' + str(idx) + '/cgff' + str(idx) + '.pickle')
    err_stat=nx.read_gpickle(aa_dir + '/Cost.pickle')

    changed=False
    for ang in unparam['angs']:
        if err_stat[idx][ang]<cost_tol:
            continue
        cgff[ang][1]=max(Kac+cgff[ang][1],10)
        changed=True

    for bond in unparam['bonds']:
        if err_stat[idx][bond]<cost_tol:
            continue
        cgff[bond][1]=Kbc+cgff[bond][1]
        changed=True

    if changed:
        nx.write_gpickle(cgff,aa_dir+'/CG'+str(idx)+'_K/cgff'+str(idx)+'_K.pickle')
        return 0
    else:
        return -1

