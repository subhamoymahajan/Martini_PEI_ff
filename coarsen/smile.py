import numpy as np
import networkx as nx
from scipy.optimize import fsolve
from scipy.signal import find_peaks
from . import aatop_2_cg
import math

def rep2smile(rep):
    if '[' not in rep:
        return rep
        
    i1=-1
    i2=-1
    num=0
    for i in range(len(rep)):
        if rep[i]=='[':
            i1=i
        elif rep[i]==']':
            i2=i
            j=1
            while i+j<len(rep):
                if not rep[i+j].isdigit():
                    break
                j+=1
            num=int(rep[i+1:i+j])
            break
    return rep2smile(rep[:i1]+rep[i1+1:i2]*num+rep[i2+j:])

def smile2top(smile):
    # CGStruc: networkx DiGraph.
    #     Nodes have attributes:
    #         'bead_name': 't', 'tq', 's', 'sq', 'p', or 'pq'
    #         'mass': total mass of the bead.
    #         'charge': charge of the bead, 0 or 1.
    if '[' in smile:
        smile=rep2smile(smile)
    i=0
    cnt=1
    CGStruc=nx.DiGraph()
    tsp='tsp'
    prevs=[]
    Mtot=0
    Qtot=0
    Ntsp=[0,0,0]
    for i in range(len(smile)):
        if smile[i] in ['(',')','q',"'"]:
            continue
        name=smile[i]
        Q=0
        M=1.00784*4+12.0107*2+14.0067 #C2H4N
        M+=tsp.index(smile[i])*1.00784
        if i==0:
            M+=1.00784
       
        if i+1<len(smile):
            if smile[i+1]=='q':
                name+='q'
                Q=1
                M+=1.00784
        CGStruc.add_node(cnt,bead_name=name,mass=M,charge=Q)
        Mtot+=M
        Qtot+=Q
      #  print('cnt:',cnt)
      #  print('bead name:',name)
      #  print('Mass:',M)
      #  print('Charge:',Q)
      #  print('prevs:',prevs)
        if cnt>1:
            if smile[i-1]!=')':
                CGStruc.add_edge(cnt-1,cnt)
                #print('Added edge '+str(cnt-1)+' -> '+str(cnt))

            else:
                idx=prevs.pop()
                CGStruc.add_edge(idx,cnt)
      #          print('Added edge from prevs:'+str(idx)+' -> '+str(cnt))
      #          print('Updated prevs:',prevs)
        if name=='t' or name=='tq':
            Ntsp[0]+=1
            prevs.append(cnt)
      #      print('added to prevs:',prevs)
        elif name=='sq' or name=='s':
            Ntsp[1]+=1
        elif name=='pq' or name=='p':
            Ntsp[2]+=1

        cnt+=1
      #  print(' ')
    print('PROPERTIES OF GENERATED MOLECULE')
    print('--------------------------------')
    print('Molecular Weight : '+str(Mtot))
    print('Charge : '+str(Qtot))
    max_=max(Ntsp)
    ratio_=[x/max_ for x in Ntsp ]
    print('Pri : Sec: Ter = ' + str(round(ratio_[2],2)) +
        ': ' + str(round(ratio_[1],2)) + ': ' + 
        str(round(ratio_[0],2)))
    print('Protonation ratio: '+str(int(0.5+100*Qtot/np.sum(Ntsp[:])))+'%')
    print(' \n\n')
    return CGStruc

##tests 1
#sqt(spq)ssqspq
#1->2->5->6->7->8
#   |->3->4
#
#print('Smile: sqt(spq)ssqspq')
#CGS=smile2top('sqt(spq)ssqspq')

##test 2
#t(st(spq)pq)ssqspq
#1->7->8->9->10
#|->2->3->6
#      |->4->5
#
#print('Smile: t(st(spq)pq)ssqspq')
#CGS=smile2top('t(st(spq)pq)ssqspq')
#
#nodes=CGS.nodes
#for n in nodes:
#    print(n,CGS.nodes[n]['bead_name'],CGS.nodes[n]['mass'],CGS.nodes[n]['charge'])
#for e1,e2 in CGS.edges():
#    print(e1,'->',e2)

def get_name(CGStruc,idxs):

    if len(idxs)==2:#Bonds
        if CGStruc.has_edge(idxs[0],idxs[1]):
            return [ CGStruc.nodes[idxs[0]]['bead_name'] + 
                CGStruc.nodes[idxs[1]]['bead_name']]

        elif CGStruc.has_edge(idxs[1],idxs[0]):
            return [CGStruc.nodes[idxs[1]]['bead_name'] + 
                CGStruc.nodes[idxs[0]]['bead_name']]

        else:
            raise Exception('No bond for ',idxs)

    if len(idxs)==3: #Angles
        if CGStruc.has_edge(idxs[0],idxs[1]) and \
            CGStruc.has_edge(idxs[1],idxs[2]): #Normal
            return [CGStruc.nodes[idxs[0]]['bead_name'] + 
                CGStruc.nodes[idxs[1]]['bead_name'] +
                CGStruc.nodes[idxs[2]]['bead_name']]

        elif CGStruc.has_edge(idxs[1],idxs[0]) and \
            CGStruc.has_edge(idxs[2],idxs[1]): #Normal
            return [CGStruc.nodes[idxs[2]]['bead_name'] + 
                CGStruc.nodes[idxs[1]]['bead_name'] +
                CGStruc.nodes[idxs[0]]['bead_name']]

        elif CGStruc.has_edge(idxs[1],idxs[0]) and \
            CGStruc.has_edge(idxs[1],idxs[2]): #Ntype
            return [ 'N' + CGStruc.nodes[idxs[0]]['bead_name'] + 
                CGStruc.nodes[idxs[1]]['bead_name'] +
                CGStruc.nodes[idxs[2]]['bead_name'] , 
                'N' + CGStruc.nodes[idxs[2]]['bead_name'] + 
                CGStruc.nodes[idxs[1]]['bead_name'] +
                CGStruc.nodes[idxs[0]]['bead_name']] 

    if len(idxs)==4: #Dih
        if CGStruc.has_edge(idxs[0],idxs[1]) and \
            CGStruc.has_edge(idxs[1],idxs[2]) and \
            CGStruc.has_edge(idxs[2],idxs[3]): #Normal
            return [CGStruc.nodes[idxs[0]]['bead_name'] + 
                CGStruc.nodes[idxs[1]]['bead_name'] +
                CGStruc.nodes[idxs[2]]['bead_name'] +
                CGStruc.nodes[idxs[3]]['bead_name']]

        elif CGStruc.has_edge(idxs[1],idxs[0]) and \
             CGStruc.has_edge(idxs[2],idxs[1]) and \
            CGStruc.has_edge(idxs[3],idxs[2]): #Normal
            return [CGStruc.nodes[idxs[3]]['bead_name'] + 
                CGStruc.nodes[idxs[2]]['bead_name'] + 
                CGStruc.nodes[idxs[1]]['bead_name'] +
                CGStruc.nodes[idxs[0]]['bead_name']]

        elif CGStruc.has_edge(idxs[1],idxs[0]) and \
            CGStruc.has_edge(idxs[1],idxs[2]) and \
            CGStruc.has_edge(idxs[2],idxs[3]): #Ntype
            return [ 'N' + CGStruc.nodes[idxs[0]]['bead_name'] + 
                CGStruc.nodes[idxs[1]]['bead_name'] +
                CGStruc.nodes[idxs[2]]['bead_name'] +
                CGStruc.nodes[idxs[3]]['bead_name']]

        elif CGStruc.has_edge(idxs[1],idxs[0]) and \
            CGStruc.has_edge(idxs[2],idxs[1]) and \
            CGStruc.has_edge(idxs[2],idxs[3]): #Ntype
            return [ 'N' + CGStruc.nodes[idxs[3]]['bead_name'] + 
                CGStruc.nodes[idxs[2]]['bead_name'] +
                CGStruc.nodes[idxs[1]]['bead_name'] +
                CGStruc.nodes[idxs[0]]['bead_name']]

def get_phi(params):
    phis=np.arange(-180,180)
    U=aatop_2_cg.ucg_dih_ref_dist(phis,params)
    peaks,_ =find_peaks(-U)
    foo=[[phis[x],U[x]] for x in peaks]
    foo=sorted(foo,key=lambda item: item[1])
    return [x[0] for x in foo]
    #phis={}
    #for i in range(params[0]):
    #    for m in range(0,params[1][i][2]+2):
    #        phi=int((2*m*180+params[1][i][0])/params[1][i][2]+0.5) #maxima
    #        if phi>360 or phi<0:
    #            continue
    #        if phi not in phis:
    #            phis[phi]=params[1][i][1]
    #        else:
    #            phis[phi]+=params[1][i][1]

    #        phi=int((2*m*180+180+params[1][i][0])/params[1][i][2]+0.5) #minima
    #        if phi>=360 or phi<0:
    #            continue
    #        if phi not in phis:
    #            phis[phi]=-params[1][i][1]
    #        else:
    #            phis[phi]-=params[1][i][1]
    #foo=sorted(phis.items(),key=lambda item: item[1])
    #print(foo)
    #return [x[0] for x in foo]
    #sorted lowest value to largest. lowest potential energy is more likely.

def get_dih_pos(p1,p2,p3,b12,b23,b34,th123,th234,phi,pos_prec):
#    print('p1',p1)
#    print('p2',p2)
#    print('p3',p3)
#    print('b12',b12)
#    print('b23',b23)
#    print('b34',b34)
#    print('th123',th123)
#    print('th234',th234)
#    print('phi',phi)
    p1=np.array(p1)
    p2=np.array(p2)
    p3=np.array(p3)
    
    def eqns(p):
       p=np.array(p)
       e1=np.dot(p-p3,p2-p3) - b23*b34*np.cos(th234*np.pi/180.0)
       foo1=np.cross(p2-p1,p3-p2)
       foo2=np.cross(p3-p2,p-p3)
       e2=math.atan2(b23*np.dot(p2-p1,foo2), np.dot(foo1,foo2))*180/np.pi-phi
       e3=np.dot(p-p3,p-p3) - b34**2
       return (e1,e2,e3)
    x,y,z=fsolve(eqns,(p3[0]+0.35,p3[1]+0.35,p3[2]+0.35))

#    print(x,y,z)
    return round(x,pos_prec), round(y,pos_prec), round(z,pos_prec)

def gen_ini_cord(CGStruc,cgff_pickle,pos_prec,smile,outname='cg_ini.gro'):
    print('Generating Initial Structure: \nWarning on initial structure can be '
        'ignored which will be corrected during energy minimization') 
    if type(cgff_pickle)==str:
        cgff=nx.read_gpickle(cgff_pickle)
    else:
        cgff=cgff_pickle
    # cgff
    # Bonds= [b0,k], Angles=[a0,k], Dihedrals=[number of functions, #[phi, K, n]
    unparam=nx.read_gpickle('unparam.pickle')
    param=nx.read_gpickle('param.pickle')
    nodes=list(CGStruc.nodes())
    N=len(nodes)
    names=[CGStruc.nodes[x]['bead_name'] for x in range(1,N+1)]

    pos=np.zeros((N,3))
    # Head node 1 is in origin.

    pos[1,0]=round(cgff[get_name(CGStruc,[1,2])[0]][0],pos_prec) #x-coord is the bond length param
    # Node 2 has to be connected to head node. Node 2 is on x-axis
    # 1->2 

    if CGStruc.has_edge(2,3): #1->2->3
        name=get_name(CGStruc,[1,2,3])
        b=cgff[get_name(CGStruc,[2,3])[0]][0]
    elif CGStruc.has_edge(1,3):#   3<-1->2
        name=get_name(CGStruc,[2,1,3])
    else:
        raise Exception('First three beads are not connected')

    for x in name:
        if x in cgff:
            th=cgff[x][0]
            break

    # Bead 3 is on xy plane. 
    if CGStruc.has_edge(1,3):
        b=cgff[get_name(CGStruc,[1,3])[0]][0]
        pos[2,0]=round(b*np.cos(th*np.pi/180.0),3)
        pos[2,1]=round(b*np.sin(th*np.pi/180.0),3)
    elif CGStruc.has_edge(2,3):
        b=cgff[get_name(CGStruc,[2,3])[0]][0]
        pos[2,0]=round(pos[1,0]-b*np.cos(th*np.pi/180.0),3)
        pos[2,1]=round(b*np.sin(th*np.pi/180.0),3)
    
    queue=[1]
    while len(queue)>0:
        bead=queue.pop()
        succs=list(CGStruc.successors(bead)) 
        queue+=succs
        if bead in [1,2,3]:
            continue
        if np.sum(pos[bead-1,:])>1E-3: #Beads are assigned positions based on their number
            continue
        
        # preds3[0] -> preds2[0] -> preds1[0] -> bead
        #               \                 \ 
        #                -> preds2_s[0]    -> preds1_s[0] -> preds1_ss[0]
        #                                        \
        #                                         -> preds1_ss[1]
        #             
        preds1=list(CGStruc.predecessors(bead))
        preds2=[]
        preds1_s=[]
        preds1_ss=[]
        if len(preds1)>0: #Preds1 will always have one bead. 
            preds2=list(CGStruc.predecessors(preds1[0])) #Exactly one bead
            preds1_s=list(CGStruc.successors(preds1[0])) #[e,Bead] or [Bead]
            preds1_s.remove(bead) 
            if len(preds1_s)>0:
                preds1_ss=list(CGStruc.successors(preds1_s[0]))
        preds3=[]
        preds2_s=[]
        if len(preds2)>0: #Preds2 will have exactly one bead
            preds3=list(CGStruc.predecessors(preds2[0])) #[a] or None
            preds2_s=list(CGStruc.successors(preds2[0])) #[c,d] or [c]
            preds2_s.remove(preds1[0]) #[d] or None

        #First Check normal dihedral angle: preds3[0]->preds2[0]->preds1[0]->bead
        if len(preds3)>0:
            b12=cgff[get_name(CGStruc,[preds3[0],preds2[0]])[0]][0]
            b23=cgff[get_name(CGStruc,[preds2[0],preds1[0]])[0]][0]
            b34=cgff[get_name(CGStruc,[preds1[0],bead])[0]][0]
            foo=get_name(CGStruc,[preds3[0],preds2[0],preds1[0]])
            for x in foo:
                if x in cgff:
                    th123=cgff[x][0]
                    break
            foo=get_name(CGStruc,[preds2[0],preds1[0],bead])
            for x in foo:
                if x in cgff:
                    th234=cgff[x][0]
                    break

            params=cgff[get_name(CGStruc,[preds3[0],preds2[0],preds1[0],bead])[0]]
            phis=get_phi(params)
            for i in range(len(phis)):
                pos[bead-1,:]=get_dih_pos(pos[preds3[0]-1,:], pos[preds2[0]-1,:], 
                    pos[preds1[0]-1,:], b12, b23, b34, th123, th234, phis[i],
                    pos_prec)
                if len(preds1_s)>0:
                    d2=np.sum(np.square(pos[bead-1,:]-pos[preds1_s[0]-1,:]))
                    if d2 > 0.05:
                        break
                else:
                    break
        #Try N-type angle 1: preds2[0]->preds1[0]->bead
        #                     \
        #                      -> preds2_s[0]
        elif len(preds2_s)>0:
            b12=cgff[get_name(CGStruc,[preds2_s[0],preds2[0]])[0]][0]
            b23=cgff[get_name(CGStruc,[preds2[0],preds1[0]])[0]][0]
            b34=cgff[get_name(CGStruc,[preds1[0],bead])[0]][0]

            foo=get_name(CGStruc,[preds2_s[0],preds2[0],preds1[0]])
            for x in foo:
                if x in cgff:
                    th123=cgff[x][0]
                    break
            foo=get_name(CGStruc,[preds2[0],preds1[0],bead])
            for x in foo:
                if x in cgff:
                    th234=cgff[x][0]
                    break
            params=cgff[get_name(CGStruc,[preds2_s[0],preds2[0],preds1[0],bead])[0]]
            foo=get_name(CGStruc,[preds2_s[0],preds2[0],preds1[0]])
            for x in foo:
                if x in cgff:
                    th123=cgff[x][0]
                    break
            foo=get_name(CGStruc,[preds2[0],preds1[0],bead])
            for x in foo:
                if x in cgff:
                    th234=cgff[x][0]
                    break
            phis=get_phi(params)
            for i in range(len(phis)):
                pos[bead-1,:]=get_dih_pos(pos[preds2_s[0]-1,:], pos[preds2[0]-1,:], 
                    pos[preds1[0]-1,:], b12, b23, b34, th123, th234, phis[i],
                    pos_prec)

                if len(preds1_s)>0:
                    d2=np.sum(np.square(pos[bead-1,:]-pos[preds1_s[0]-1,:]))
                    if d2 > 0.05:
                        break
                else:
                    break
        #Try N-type angle 2:  preds1[0]->bead
        #                      \
        #                       ->preds1_s[0]->preds1_ss[0]
        elif len(preds1_ss)>0:
            b12=cgff[get_name(CGStruc,[preds1_ss[0],preds1_s[0]])[0]][0]
            b23=cgff[get_name(CGStruc,[preds1_s[0],preds1[0]])[0]][0]
            b34=cgff[get_name(CGStruc,[preds1[0],bead])[0]][0]
            foo=get_name(CGStruc,[preds1_ss[0],preds1_s[0],preds1[0]])
            for x in foo:
                if x in cgff:
                    th123=cgff[x][0]
                    break
            foo=get_name(CGStruc,[preds1_s[0],preds1[0],bead])
            for x in foo:
                if x in cgff:
                    th234=cgff[x][0]
                    break
            params=cgff[get_name(CGStruc,[preds1_ss[0],preds1_s[0],preds1[0],bead])[0]]
            phis=get_phi(params)
            pos[bead-1,:]=get_dih_pos(pos[preds1_ss[0]-1,:], pos[preds1_s[0]-1,:], 
                pos[preds1[0]-1,:], b12, b23, b34, th123, th234, phis[0], pos_prec)
    foo=np.amin(pos,axis=0)
    pos=pos-foo
    write_gro(CGStruc,pos,smile,outname)


def write_gro(CGStruc,pos,smile,outname='cg_ini.gro'):
    N=len(list(CGStruc.nodes()))
    w=open(outname,'w')
    w.write('Initialized by Coarsen SMILE = '+smile+' Needs energy minimization\n   '+str(N)+'\n')
    name={'t':'Ter', 'tq':'QTer', 's':'Sec', 'sq':'QSec', 'p':'Pri', 'pq':'QPri'}
    for i in range(1,N+1):
       w.write('%5d'%i+'  PEI'+'%5s'%name[CGStruc.nodes[i]['bead_name']]+'%5d'%i)
       w.write('%8.3f'%pos[i-1,0]+'%8.3f'%pos[i-1,1]+'%8.3f'%pos[i-1,2]+'\n')
    box=np.amax(pos,axis=0)
    box=[round(x+0.001,3) for x in box]
    w.write(' '+str(box[0])+' '+str(box[1])+' '+str(box[2])+'\n')
    w.close()

def suggest_nodes(Mn,sec_ter,PR):
    m=[1.00785,44.07514,43.06713,42.0595]
    #Ax=B
    n3=int((Mn-m[1]-m[0]*(1.+PR/100.))/(m[1]+sec_ter*m[2]+m[3]+(2+sec_ter)*m[0]*PR/100.) + 0.5)
    n2=int(n3/sec_ter + 0.5)
    n1=n3+1
    nH=int((2.+sec_ter)*n3*PR/100.+PR/100. +0.5 )
    print('Number of protonated: '+str(nH))
    print('Number of primary: '+str(n1))
    print('Number of secondary: '+str(n2))
    print('Number of tertiary: '+str(n3))
