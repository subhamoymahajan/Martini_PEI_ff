import os
import matplotlib.pyplot as plt
import numpy as np
import networkx as nx
from . import aatop_2_cg
from scipy.optimize import curve_fit
kbt=aatop_2_cg.kbt
tol=aatop_2_cg.tol

def gen_mol_prop(cgstruc_pickle,props):
    if type(cgstruc_pickle)==str:
        cgstruc=nx.read_gpickle(cgstruc_pickle)
    else:
        cgstruc=cgstruc_pickle
    M=0
    Q=0
    Ntsp=np.zeros(3)
    for node in cgstruc.nodes():
        if cgstruc.nodes[node]['bead_name'] in ['tq', 't']:
            Ntsp[0]+=1
        elif cgstruc.nodes[node]['bead_name'] in ['sq', 's']:
            Ntsp[1]+=1
        elif cgstruc.nodes[node]['bead_name'] in ['pq', 'p']:
            Ntsp[2]+=1
        M+=cgstruc.nodes[node]['mass']
        Q+=cgstruc.nodes[node]['charge']
    print(Ntsp)

    props['MolWt']=round(M)
    props['Charge']=Q
    props['pr']=round(Q/np.sum(Ntsp),3)
    props['N_tsp']=Ntsp

    print('Molecular Weight : '+str(round(M)))
    print('Charge           : '+str(Q))
    print('Protonation ratio: '+str(round(Q*100/np.sum(Ntsp),1)))
    Ntsp=Ntsp/Ntsp[0]
    print('Pri/Sec/Ter      : '+str(round(Ntsp[2],2))+'/'+str(round(Ntsp[1],2))+'/1')
    nx.write_gpickle(props,'prop.pickle')
       
def gen_png2(dirs, labels, bond_small=None, bond_large=None, bond_ymax=None, 
    ang_ymax=None, dih_ymax=None, ncols_val=[1,1,1]):
    import matplotlib
    cmap=matplotlib.cm.get_cmap('hsv')
    os.system('mkdir -p comparing_pngs')
    unparam=nx.read_gpickle(dirs[0]+'unparam.pickle')
    param=nx.read_gpickle(dirs[0]+'param.pickle')
    cols=['k']
    for i in range(2,len(dirs)):
        cols.append(cmap((i-2.)/float(len(dirs)-2)))

    for typ in ["bonds", "angs"]:
        typ1=typ[:-1]
        if typ1=="ang":
            typ1="angle"

        for name in list(unparam[typ].keys())+param[typ]:
            fig,ax=plt.subplots(1,1,figsize=(3,3))

            for i in range(len(labels)):
                data=aatop_2_cg.get_dist(dirs[i+1]+'bonded_distribution/' + typ1 + 
                    '_' + name + '.xvg')
                ax.plot(data[:,0], data[:,1], color=cols[i], label=labels[i])

            if typ=='bonds':
                ax.set_xlabel('Bond length '+name+' (nm)')
                ax.set_ylim(bottom=0)
                if bond_ymax is not None:
                    ax.set_ylim(top=bond_ymax)
                if bond_small is not None and bond_large is not None:
                    ax.set_xlim(bond_small,bond_large)

            else: 
                ax.set_xlabel('Bond angle '+name+' (deg)')
                ax.set_xlim(0,180)
                ax.set_xticks(np.arange(0,181,30))
                ax.set_ylim(bottom=0)
                if ang_ymax is not None:
                    ax.set_ylim(top=ang_ymax)

            if name in unparam[typ]:
                ax.xaxis.label.set_color('red')
            ax.set_ylabel('Probability')
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            plt.tight_layout()
            if typ == "bonds": 
                plt.legend(loc='upper right', fontsize='small', ncol=ncols_val[0])
            if typ == "angs": 
                plt.legend(loc='upper right', fontsize='small', ncol=ncols_val[1])
            plt.savefig('comparing_pngs/' + typ1 + '_' + name + '.png', dpi=300)
            print('Writing: ' + typ1 + '_' + name + '.png') 
            plt.clf()
            plt.close()

    for name in list(unparam['dihs'].keys()) + param['dihs']:
        fig,ax=plt.subplots(1, 1, figsize=(3,3))
        for i in range(len(labels)):
            data=aatop_2_cg.get_dist(dirs[i+1]+'bonded_distribution/dih_' + name + '.xvg')
            ax.plot(data[:,0], data[:,1], color=cols[i], label=labels[i])
        ax.set_xlim(-180,180)
        ax.set_xticks(np.arange(-180,181,60))
        ax.set_xticklabels(np.arange(-180,181,60),rotation=45)
        ax.set_ylim(bottom=0)
        if dih_ymax is not None:
            ax.set_ylim(top=dih_ymax)
        ax.set_xlabel('Dihedral angle '+name+' (deg)')
        ax.set_ylabel('Probability')
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        if name in unparam['dihs']:
            ax.xaxis.label.set_color('red')
        plt.tight_layout()
        plt.legend(loc='upper right', fontsize='small', ncol=ncols_val[2])
        plt.savefig('comparing_pngs/dih_' + name + '.png',dpi=300)
        print('Writing: dih_' + name + '.png') 
        plt.clf()
        plt.close()

def gen_png(cg_idx, bond_small=None, bond_large=None, bond_ymax=None, 
    ang_ymax=None, dih_ymax=None):
    os.system('mkdir -p CG'+str(cg_idx)+'/comparing_pngs')
    unparam=nx.read_gpickle('unparam.pickle')
    param=nx.read_gpickle('param.pickle')
    c=np.linspace(0.5,1,3)
    cm=[(a,1-a,a) for a in c]
    for typ in ["bonds", "angs"]:
        typ1=typ[:-1]
        if typ1=="ang":
            typ1="angle"

        for name in list(unparam[typ].keys())+param[typ]:
            fig,ax=plt.subplots(1,1,figsize=(3,3))
            cg_data=aatop_2_cg.get_dist('CG' + str(cg_idx) + \
                '/bonded_distribution/' + typ1 + '_' + name + '.xvg')
            aa_data=aatop_2_cg.get_dist('bonded_distribution/' + typ1 + \
                '_' + name + '.xvg')
            ax.plot(aa_data[:,0], aa_data[:,1], 'k', label='AA')
            ax.plot(cg_data[:,0], cg_data[:,1], 'r', label='CG' + str(cg_idx))
            if cg_idx>1:
                for i in range(cg_idx-1,max(cg_idx-3,0),-1):
                    cg_data_old=aatop_2_cg.get_dist('CG' + str(i) + \
                        '/bonded_distribution/' + typ1 + '_' + name + '.xvg')
                    ax.plot(cg_data_old[:,0], cg_data_old[:,1], color=cm[i%3], \
                        label='CG' + str(i))
            if typ=='bonds':
                ax.set_xlabel('Bond length '+name+' (nm)')
                ax.set_ylim(bottom=0)
                if bond_ymax is not None:
                    ax.set_ylim(top=bond_ymax)
                if bond_small is not None and bond_large is not None:
                    ax.set_xlim(bond_small,bond_large)
            else: 
                ax.set_xlabel('Bond angle '+name+' (deg)')
                ax.set_xlim(0,180)
                ax.set_xticks(np.arange(0,181,30))
                ax.set_ylim(bottom=0)
                if ang_ymax is not None:
                    ax.set_ylim(top=ang_ymax)
            if name in unparam[typ]:
                ax.xaxis.label.set_color('red')
            ax.set_ylabel('Probability')
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            plt.tight_layout()
            plt.legend(loc='upper right', fontsize='small')
            plt.savefig('CG' + str(cg_idx) + '/comparing_pngs/' + typ1 + \
                '_' + name + '.png', dpi=300)
            print('Writing: ' + typ1 + '_' + name + '.png') 
            plt.clf()
            plt.close()

    for name in list(unparam['dihs'].keys()) + param['dihs']:
        fig,ax=plt.subplots(1, 1, figsize=(3,3))
        cg_data=aatop_2_cg.get_dist('CG' + str(cg_idx) + \
            '/bonded_distribution/dih_' + name + '.xvg')
        aa_data=aatop_2_cg.get_dist('bonded_distribution/dih_' + name + '.xvg')
        ax.plot(aa_data[:,0], aa_data[:,1], 'k', label='AA')
        ax.plot(cg_data[:,0], cg_data[:,1], 'r', label='CG' + str(cg_idx))
        if cg_idx > 1:
            for i in range(cg_idx-1,max(cg_idx-3,0),-1):
                cg_data_old=aatop_2_cg.get_dist('CG' + str(i) + \
                    '/bonded_distribution/dih_' + name + '.xvg')
                ax.plot(cg_data_old[:,0], cg_data_old[:,1], color=cm[i%3], 
                    label='CG' + str(i))
        ax.set_xlim(-180,180)
        ax.set_xticks(np.arange(-180,181,60))
        ax.set_ylim(bottom=0)
        if dih_ymax is not None:
            ax.set_ylim(top=dih_ymax)
        ax.set_xlabel('Dihedral angle '+name+' (deg)')
        ax.set_ylabel('Probability')
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        if name in unparam['dihs']:
            ax.xaxis.label.set_color('red')
        plt.tight_layout()
        plt.legend(loc='upper right', fontsize='small')
        plt.savefig('CG' + str(cg_idx) + '/comparing_pngs/dih_' + name + \
            '.png',dpi=300)
        print('Writing: dih_' + name + '.png') 
        plt.clf()
        plt.close()

def create_latex(cg_dir,aa_dir='',fig_per_row=3):
    unparam=nx.read_gpickle(aa_dir+'unparam.pickle')
    param=nx.read_gpickle(aa_dir+'param.pickle')

    w=open(cg_dir + 'comparing_pngs/all_images.tex', 'w')
    
    w.write('\\documentclass{article}\n\n\\usepackage{graphicx}\n' + \
        '\\usepackage[left=2cm, top=2cm, bottom=2cm, right=2cm]{geometry}\n' + \
        '\\usepackage{caption}\n\\usepackage{xcolor}\n\\usepackage{subcapti' + \
        'on}\n\\usepackage{float}\n\n\\begin{document}\n\n')
    
    w.write('\\section{Bonds}\n\n')
    cnt=0
    for name in list(unparam['bonds'].keys()) + param['bonds']:
        if cnt%fig_per_row==0:
            w.write('\t\\begin{figure}[H]\n\t\t\\centering\n')
        w.write('\t\t\\begin{subfigure}{0.3\\textwidth}\n')
        if name in unparam['bonds']:
            w.write('\t\t\t\\captionsetup{font={color=red,bf}')
        w.write('\t\t\t\\includegraphics[width=\\textwidth]{bond_' + name + \
            '.png}\n\t\t\\end{subfigure}\n')

        if cnt%fig_per_row==fig_per_row-1 or cnt==len(unparam['bonds'].keys()) + \
            len(param['bonds'])-1:
            w.write('\t\\end{figure}\n\n')
        cnt+=1
    
    w.write('\\section{Angles}\n\n')
    cnt=0
    for name in list(unparam['angs'].keys())+param['angs']:
        if cnt%fig_per_row==0:
            w.write('\t\\begin{figure}[H]\n\t\t\\centering\n')
        w.write('\t\t\\begin{subfigure}{0.3\\textwidth}\n')
        if name in unparam['angs']:
            w.write('\t\t\t\\captionsetup{font={color=red,bf}}\n')
        w.write('\t\t\t\\includegraphics[width=\\textwidth]{angle_' + name + \
            '.png}\n\t\t\\end{subfigure}\n')
        if cnt%fig_per_row==fig_per_row - 1 or cnt==len(unparam['angs'].keys()) + \
            len(param['angs']) - 1:
            w.write('\t\\end{figure}\n\n')
        cnt+=1
    
    w.write('\\section{Dihedrals}\n\n')
    cnt=0
    for name in list(unparam['dihs'].keys()) + param['dihs']:
        if cnt%fig_per_row==0:
            w.write('\t\\begin{figure}[H]\n\t\t\\centering\n')
        w.write('\t\t\\begin{subfigure}{0.3\\textwidth}\n')
        if name in unparam['dihs'].keys():
            w.write('\t\t\t\\captionsetup{font={color=red,bf}}\n')
        w.write('\t\t\t\\includegraphics[width=\\textwidth]{dih_' + name + \
            '.png}\n\t\t\\end{subfigure}\n')
        if cnt%fig_per_row==fig_per_row-1 or cnt==len(unparam['dihs'].keys())+ \
            len(param['dihs']) - 1:
            w.write('\t\\end{figure}\n\n')
        cnt+=1
    w.write('\\end{document}')
    w.close()

def add_cost(wb,wa,aa_dir,cg_dir,idx):
    aa_dir_b=aa_dir+'/bonded_distribution/'
    cg_dir_b=cg_dir+'/bonded_distribution/'
    unparam=nx.read_gpickle(aa_dir+'/unparam.pickle')
    param=nx.read_gpickle(aa_dir+'/param.pickle')
    ff=nx.read_gpickle(cg_dir+'/cgff'+str(idx)+'.pickle')

    if not os.path.exists('Cost.pickle'):
        Stats={}
    else:
        Stats=nx.read_gpickle('Cost.pickle')

    Stats[idx]={}
    for a in ['all', 'param', 'unparam']:
        for b in ['', '_bond', '_ang', '_dih']:
            Stats[idx][a+b]=0

    data={}
    mean={}
    std={}
    U={}
    for typ in ["bonds", "angs"]:
        if typ=="bonds":
            w=wb
        else:
            w=wa
        for name in list(unparam[typ].keys())+param[typ]:
            typ2=typ[:-1]
            if typ=='angs':
                typ2='angle'
            data['AA']=aatop_2_cg.get_dist(aa_dir_b+typ2+'_'+name+'.xvg')
            mean['AA']=aatop_2_cg.mean_pdf(data['AA'])
            std['AA']=aatop_2_cg.std_pdf(data['AA'])
            data['CG']=aatop_2_cg.get_dist(cg_dir_b+typ2+'_'+name+'.xvg')
            mean['CG']=aatop_2_cg.mean_pdf(data['CG'])
            std['CG']=aatop_2_cg.std_pdf(data['CG'])
            
            foo=w*(mean['CG']-mean['AA'])**2 +(1-w)*(std['CG']-std['AA'])**2
            Stats[idx][name]=foo
            Stats[idx]['all']+=foo
            Stats[idx]['all_'+typ[:-1]]+=foo
            
            if name in unparam[typ]:
                Stats[idx]['unparam']+=foo
                Stats[idx]['unparam_'+typ[:-1]]+=foo
            else:
                Stats[idx]['param']+=foo
                Stats[idx]['param_'+typ[:-1]]+=foo

    for dih in list(unparam['dihs'].keys())+param['dihs']:
        data['AA']=aatop_2_cg.get_dist(aa_dir_b+'dih_'+dih+'.xvg')
        data['AA']=aatop_2_cg.fill_data(data['AA'])
        U['AA']=-kbt*np.log(np.array(data['AA'][:,1])+tol) 
        p0=[1]*(ff[dih][0]+1)
        popt1,pcov1=curve_fit(aatop_2_cg.fourier_series,data['AA'][:,0],U['AA'],p0)

        data['CG']=aatop_2_cg.get_dist(cg_dir_b+'dih_'+dih+'.xvg')
        data['CG']=aatop_2_cg.fill_data(data['CG'])
        U['CG']=-kbt*np.log(np.array(data['CG'][:,1])+tol) 
        p0=[1]*(ff[dih][0]+1)
        popt2,pcov2=curve_fit(aatop_2_cg.fourier_series,data['AA'][:,0],U['CG'],p0)

        foo=np.sum(np.square(np.subtract(popt1[1:],popt2[1:])))
        Stats[idx][dih]=foo
        Stats[idx]['all']+=foo
        Stats[idx]['all_dih']+=foo
        if dih in unparam['dihs']:
            Stats[idx]['unparam']+=foo
            Stats[idx]['unparam_dih']+=foo
        else:
            Stats[idx]['param']+=foo
            Stats[idx]['param_dih']+=foo
    #Cost Rg
    Rg_aa=aatop_2_cg.get_dist(aa_dir+'/polystat.xvg')
    Rg_cg=aatop_2_cg.get_dist(cg_dir+'/polystat.xvg')
    aa_avg=np.average(Rg_aa[:,2])
    aa_avg2=np.average(np.square(Rg_aa[:,2]))
    aa_std=np.sqrt(aa_avg2-aa_avg**2)

    cg_avg=np.average(Rg_cg[:,2])
    cg_avg2=np.average(np.square(Rg_cg[:,2]))
    cg_std=np.sqrt(cg_avg2-cg_avg**2)
    Stats[idx]['Rg']=(aa_avg-cg_avg)**2+(aa_std-cg_std)**2

    #Cost Re
    Re_aa=aatop_2_cg.get_dist(aa_dir+'/polystat.xvg')
    Re_cg=aatop_2_cg.get_dist(cg_dir+'/polystat.xvg')
    aa_avg=np.average(Re_aa[:,1])
    aa_avg2=np.average(np.square(Re_aa[:,1]))
    aa_std=np.sqrt(aa_avg2-aa_avg**2)

    cg_avg=np.average(Re_cg[:,1])
    cg_avg2=np.average(np.square(Re_cg[:,1]))
    cg_std=np.sqrt(cg_avg2-cg_avg**2)
    Stats[idx]['Re']=(aa_avg-cg_avg)**2+(aa_std-cg_std)**2
    nx.write_gpickle(Stats,'Cost.pickle')

def calc_polystat(props):
    Rg=aatop_2_cg.get_dist('polystat.xvg')
    props['Re']=[round(np.average(Rg[:,1]),4),
                 round(np.std(Rg[:,1]),4)]
    print('Re : '+str(props['Re'][0])+' +- '+str(props['Re'][1]))

    props['Rg']=[round(np.average(Rg[:,2]),4),
                 round(np.std(Rg[:,2]),4)]
    print('Rg : '+str(props['Rg'][0])+' +- '+str(props['Rg'][1]))

    props['Rg_eig']=[[round(np.average(Rg[:,3]),4),
                      round(np.average(Rg[:,4]),4),
                      round(np.average(Rg[:,5]),4)],
                     [round(np.std(Rg[:,3]),4),
                      round(np.std(Rg[:,4]),4),
                      round(np.std(Rg[:,5]),4)]]
    props['shape']={}
    foo=sorted(props['Rg_eig'][0])
    print('Eigen value of Rg: ',foo)

    props['shape']['asphericity']=round((1.5*foo[-1]**2 - 0.5*props['Rg'][0]**2),4)
    props['shape']['acylindricity']=round((foo[1]-foo[0])*(foo[1]+foo[0]),4)
    props['shape']['shape_anisotropy']=round((props['shape']['asphericity']**2 + \
        0.75*props['shape']['acylindricity']**2)/(props['Rg'][0]**4),4)
    nx.write_gpickle(props,'prop.pickle')

def calc_diff_scaling(data,martini):
    if martini not in ['2.1P-dna', '2.2refP','2.P','2.1P']:
        print('Scaling has not been determined.')
        return
    if martini=='2.1P-dna' or martini=='2.P':
        martini='2.1P'
    props=nx.read_gpickle('prop.pickle')
    pr=props['pr']
    conc=props['conc']
    print('Using martini'+martini)
    print('PR: ',pr)
    print('Conc.: ',conc)
    num=(data['P1_AA']['intercept'] + conc*data['P1_AA']['coef'])*(1-pr) + \
        (data['Qd_AA']['intercept'] + conc*data['Qd_AA']['coef'])*pr
    den=(data['P1_'+martini]['intercept'] + conc*data['P1_'+martini]['coef'])*(1-pr) + \
        (data['Qd_'+martini]['intercept'] + conc*data['Qd_'+martini]['coef'])*pr

    props['scaling']=round(num/den,2)
    print('Diffusion scaling: '+str(props['scaling']))
    nx.write_gpickle(props,'prop.pickle')

def calc_Diff(fname,tmax,tmul=0.001,tunit='ns',show=False):

    data=aatop_2_cg.get_dist('msd_'+fname+'.xvg')
    try:
        props=nx.read_gpickle('prop.pickle')
    except:
        props={}
    data[:,0]=data[:,0]*tmul
    dt=data[1,0]-data[0,0]
    N1=int(tmax/dt+0.5)
    data=data[:N1,:]

    from scipy import stats
    
    a,b,r,p,std_err=stats.linregress(data[:,0],data[:,1])

    print("Diffusion coefficient : "+str(round(a/6.,4))+' +- '+
        str(round(std_err/6.,4))+' 1E-5 cm2/s')
    props['D']=round(a/6.,4)
    props['D_err']=round(std_err/6.,4)
    try:
        print("Diffusion coefficient : "+str(round(a*props['scaling']/6.,4))+' +- '+
            str(round(std_err*props['scaling']/6.,4))+' 1E-5 cm2/s')
        props['D_scaled']=round(a*props['scaling']/6.,4)
        props['D_scaled_err']=round(std_err*props['scaling']/6.,4)
    except:
        pass
    nx.write_gpickle(props,'prop.pickle')
    y_pred=data[:,0]*a+b

    import matplotlib
    matplotlib.use('TkAgg')
    
    fig,ax=plt.subplots(1,1,figsize=(4,4))
    ax.plot(data[:,0],data[:,1],color='k')
    ax.plot(data[:,0],y_pred,color='r',linestyle='dashed')
    
    foo_dx=[0.01,0.02,0.04,0.05,0.1,0.2,0.4,0.5,1,2,4,5,10,20,40,50]
    nticks_foo=np.array([abs(int(tmax/x+0.5)-5) for x in foo_dx])
    idx=np.argmin(nticks_foo)
    dx=foo_dx[idx]
    xmax=(int(tmax/dx+0.5)+1)*dx

    nticks_foo=np.array([abs(int(max(data[:,1])/x+0.5)-5) for x in foo_dx])
    idx=np.argmin(nticks_foo)
    dy=foo_dx[idx]
    ymax=(int(max(data[:,1])/dy+0.5)+2)*dy
    
    
    ax.set_ylabel(r'MSD (nm$^2$)')
    ax.set_xlabel('Time ('+tunit+')')
    ax.set_xticks(np.arange(0,xmax,dx))
    ax.set_xlim(0,tmax)
    ax.set_yticks(np.arange(0,ymax,dy))
    ax.set_ylim(0,ymax-dy)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    lines=ax.get_lines()
    plt.tight_layout()
    ax.legend([lines[1]],['fit'], fontsize=10, bbox_to_anchor=(1.0,1.0), framealpha=1.0,
               edgecolor='k')
    if show:
        plt.show()
    else:
        plt.savefig('msd_'+fname+'.png',dpi=1200)

def calc_Diff_N(fname,tmax,T0,TN,N,tmul=0.001,tunit='ns',show=False):

    datas=[]
    As=[]
    Bs=[]
    from scipy import stats
    for i in range(N):
        t1=int(T0+i*int((TN-T0)/N))
        t2=int(T0+(i+1)*int((TN-T0)/N))
        data=aatop_2_cg.get_dist('msd_'+fname+'_'+str(int(t1/1000))+
            '-'+str(int(t2/1000))+'ns.xvg')
        data[:,0]=data[:,0]*tmul
        dt=data[1,0]-data[0,0]
        N1=int(tmax/dt+0.5)
        data=data[:N1,:]
        datas.append(data)
        a,b,r,p,std_err=stats.linregress(data[:,0],data[:,1])
        As.append(a)
        Bs.append(b)
    As=np.array(As)
    a_avg=np.average(As)
    a_std=np.std(As)

    Bs=np.array(Bs)
    b_avg=np.average(Bs)
    b_std=np.std(Bs)

    Davg=a_avg/6.
    Derr=a_std/6.

    datas=np.array(datas)
    foo=datas.shape
    data_comp=np.zeros((foo[1],3))
    data_comp[:,0]=datas[0,:,0]
    data_comp[:,1]=np.average(datas[:,:,1],axis=0)
    data_comp[:,2]=np.std(datas[:,:,1],axis=0)

    try:
        props=nx.read_gpickle('prop.pickle')
    except:
        props={}

    

    props['D']=round(Davg,4)
    props['D_err']=round(Derr/6.,4)
    print("Diffusion coefficient : "+str(props['D'])+' +- '+
        str(props['D_err'])+' 1E-5 cm2/s')

    try:
        props['D_scaled']=round(props['D']*props['scaling'],4)
        props['D_scaled_err']=round(props['D_err']*props['scaling'],4)
        print("Diffusion coefficient : "+str(props['D_scaled'])+' +- '+
            str(props['D_scaled_err'])+' 1E-5 cm2/s')
    except:
        pass

    nx.write_gpickle(props,'prop.pickle')
    y_pred=data[:,0]*a_avg+b_avg
    y_p1=y_pred - (data[:,0]*a_std+b_std)
    y_p2=y_pred + (data[:,0]*a_std+b_std)

    import matplotlib
    matplotlib.use('TkAgg')
    
    fig,ax=plt.subplots(1,1,figsize=(4,4))
    ax.plot(data_comp[:,0],data_comp[:,1],color='k')
    ax.fill_between(data_comp[:,0],data_comp[:,1]-data_comp[:,2], 
        y2=data_comp[:,1]+data_comp[:,2], alpha=0.5, facecolor='k', edgecolor=None)

    ax.plot(data_comp[:,0],y_pred,color='r',linestyle='dashed')
    ax.fill_between(data_comp[:,0],y_p1, y2=y_p2,  
        alpha=0.5, facecolor='r', edgecolor=None)
    
    foo_dx=[0.01,0.02,0.04,0.05,0.1,0.2,0.4,0.5,1,2,4,5,10,20,40,50]
    nticks_foo=np.array([abs(int(tmax/x+0.5)-5) for x in foo_dx])
    idx=np.argmin(nticks_foo)
    dx=foo_dx[idx]
    xmax=(int(tmax/dx+0.5)+1)*dx

    nticks_foo=np.array([abs(int(max(data[:,1])/x+0.5)-5) for x in foo_dx])
    idx=np.argmin(nticks_foo)
    dy=foo_dx[idx]
    ymax=(int(max(data[:,1])/dy+0.5)+2)*dy
    
    
    ax.set_ylabel(r'MSD (nm$^2$)')
    ax.set_xlabel('Time ('+tunit+')')
    ax.set_xticks(np.arange(0,xmax,dx))
    ax.set_xlim(0,tmax)
    ax.set_yticks(np.arange(0,ymax,dy))
    ax.set_ylim(0,ymax-dy)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    lines=ax.get_lines()
    plt.tight_layout()
    ax.legend([lines[1]],['fit'], fontsize=10, bbox_to_anchor=(1.0,1.0), framealpha=1.0,
               edgecolor='k')
    if show:
        plt.show()
    else:
        plt.savefig('msd_'+fname+'.png',dpi=1200)
         
