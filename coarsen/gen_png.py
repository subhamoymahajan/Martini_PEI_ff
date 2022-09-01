import glob
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

def import_data(filename):
    #IMPORT BOND ANGLE AND DIHEDRAL DISTRIBUTIONS
    f=open(filename,'r')
    x=[]
    y=[]
    for lines in f:
        foo=lines.split()
        if foo[0][0]=='#' or foo[0][0]=='@':
            continue
        x.append(float(foo[0]))
        y.append(float(foo[1]))
    return [x,y]

if not os.path.exists('comparing_pngs'):
    os.mkdir('comparing_pngs')

filenames=glob.glob('bonded_distribution/*.xvg')
print(filenames[0])
for f in filenames:
    cur_folder=os.path.abspath('.').split('/')[-1]
    if '_' in cur_folder:
        foo=cur_folder.split('_')[-2]
        print(foo)
    else:
        foo=cur_folder
    n=int(foo.split('G')[-1])
    #n=int(os.path.abspath('.').split('G')[-1])
    #print(n)
    c=np.linspace(0.5,1,3)
    cm=[(a,1-a,a) for a in c]
    plt.figure()
    print('here')
    dist_name=f.split('/')[1].split('.')[0]
    cg_data=import_data(f)
    old={}
    for i in range(1,n):
        if i < n-3:
            continue
        old['cgold_data'+str(i)]=import_data('../CG'+str(i)+'/'+f)

    aa_data=import_data('../'+f)
    if dist_name[0]=='b':
        plt.xlabel('Bond length (nm)')
    elif dist_name[0]=='a':
        plt.xlabel('Bond angle (deg)')
        plt.xlim(0,180)
    elif dist_name[0]=='d':
        plt.xlabel('Dihedral angle (deg)')
        plt.xlim(-180,180)
    plt.ylabel('Probability')
    plt.title(dist_name)

    plt.plot(aa_data[0],aa_data[1],'k',label='AA')
    for i in range(1,n):
        if i < n-3:
            continue
        plt.plot(old['cgold_data'+str(i)][0],old['cgold_data'+str(i)][1],color=cm[(i-1)%3],linestyle='--',label='CG'+str(i))
    plt.plot(cg_data[0],cg_data[1],'r',label='CG')
    plt.legend()
    plt.savefig('comparing_pngs/'+dist_name+'.png',dpi=300)
    plt.clf()
    print(dist_name+' done')

