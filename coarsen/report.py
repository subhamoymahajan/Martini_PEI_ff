import os
import matplotlib.pyplot as plt
import numpy as np
import networkx as nx
from . import aatop_2_cg
from scipy.optimize import curve_fit
kbt=aatop_2_cg.kbt
tol=aatop_2_cg.tol

def gen_mol_prop(cgstruc_pickle):
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

    print('Molecular Weight : '+str(round(M)))
    print('Charge           : '+str(Q))
    print('Protonation ratio: '+str(round(Q*100/np.sum(Ntsp),1)))
    Ntsp=Ntsp/Ntsp[0]
    print('Pri/Sec/Ter      : '+str(round(Ntsp[2],2))+'/'+str(round(Ntsp[1],2))+'/1')
        

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

def create_latex(cg_idx,fig_per_row=3):
    unparam=nx.read_gpickle('unparam.pickle')
    param=nx.read_gpickle('param.pickle')

    w=open('CG' + str(cg_idx) + '/comparing_pngs/all_images.tex', 'w')
    
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
    aa_avg=np.average(Rg_aa[:,1])
    aa_avg2=np.average(np.square(Rg_aa[:,1]))
    aa_std=np.sqrt(aa_avg2-aa_avg**2)
    cg_avg=np.average(Rg_cg[:,1])
    cg_avg2=np.average(np.square(Rg_cg[:,1]))
    cg_std=np.sqrt(cg_avg2-cg_avg**2)
    Stats[idx]['Rg']=(aa_avg-cg_avg)**2+(aa_std-cg_std)**2

    nx.write_gpickle(Stats,'Cost.pickle')
    
