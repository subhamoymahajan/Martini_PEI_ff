#    Copyright 2022 SUBHAMOY MAHAJAN 
#    
#    This file is part of InSilicoMicroscopy software.
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

from optparse import OptionParser
from coarsen import aatop_2_cg
from coarsen import report
from coarsen import smile
import networkx as nx
import os
import subprocess
import numpy as np

def main():
    parser = OptionParser()
    parser.add_option('-f', '--file', dest="file", metavar="FILE", 
                      type="str",
                      help="Provide filename")
    parser.add_option('-s', '--tpr', dest="tpr_file", metavar="FILE", 
                      type="str",
                      help="TPR filename")
    parser.add_option('-x', '--paramfile', dest="paramfile", metavar="FILE", 
                      type="str",
                      help="Prameter filename")
    parser.add_option('-p', '--topol', dest="top", metavar="FILE", 
                      type="str", 
                      help="topology file")
    parser.add_option('-n', '--ndx', dest="ndx", metavar="FILE", 
                      type="str", help="Provide output ndx")
    parser.add_option('-o', '--output', dest="outname", metavar="FILE", 
                      type="str", help="Provide output filename")
    parser.add_option('-b', '--begint', dest="begin", metavar="CONST", 
                      type="int", help="Start time in ps")
    parser.add_option('-e', '--endt', dest="end", metavar="CONST", 
                      type="int", help="End time in ps")
    parser.add_option('-d', '--dir', dest="dir", metavar="DIRECTORY", 
                      type="str", help="Path of a directory")
    parser.add_option('-k', '--key', dest="key", metavar="STRING", 
                      type="str", help="A string of text")
    parser.add_option('-i', '--id', dest="id", metavar="INT", 
                      type="int", help="An integer ID")

#    parser.add_option('-s', '--multiprocess', dest="mprocess",
#                      action="store_true", default=False, 
#                      help="Use for multiprocessing")
    
    options, remainder = parser.parse_args()
    params={}
    params['cgff_2019']=os.path.dirname(__file__)+'/cgff_2019.pickle'
    params['cgff_curr']=os.path.dirname(__file__)+'/cgff_2019.pickle'
    params['dih_initial']=None
    params['bond_small']=0.3
    params['bond_large']=0.6
    params['dih_ymax']=0.001
    params['ang_ymax']=0.01
    params['bond_ymax']=0.05
    params['gpu']=0
    params['cost_tol']=1E-8
    params['npt_mdp']='npt.mdp'
    params['em_mdp']='em.mdp'
    params['ions_mdp']='ions.mdp'
    params['md_mdp']='md.mdp'
    params['wa']=0.5
    params['wb']=0.5
    params['pos_prec']=3
    params['start_iter']=1
    if options.paramfile != None:
        f=open(options.paramfile,'r')
        for lines in f:
            if lines[0]=='#':
                continue
            foo=lines.split('=')
            var=foo[0].strip()
            if var in ["peiname", "trr", "init", "npt_mdp", "ions_mdp", 
                "em_mdp", "md_mdp", "dih_initial", "cgff_curr"]:
                params[var]=foo[1].strip()
            elif var in ["max_iter", "gpu", "start_iter"]:
                params[var]=int(foo[1])
            elif var in ["thc", "kc", "kdc", "fa", "fd", 'bond_small', 'wa',
                'bond_large', 'dih_ymax', 'ang_ymax', 'bond_ymax', 'wb',
                'cost_tol']:
                params[var]=float(foo[1])
    if "xtc" not in params:
        if "trr" in params:
            params["xtc"]=params["trr"]
                
    if "init" in params:
        f=open('init.sh','w')
        f.write('#!/bin/bash\n')
        foo=params['init'].split(';')
        for x in foo:
            f.write(x+'\n')
        f.close()
    if remainder[0]=='aa2cg':
        subprocess.run('cp ' + os.path.dirname(__file__) + '/martini_v2.0_ions.itp .', shell=True, check=True)
        subprocess.run('cp ' + os.path.dirname(__file__) + '/martini_v2.2refP.itp .', shell=True, check=True)
        subprocess.run('cp ' + os.path.dirname(__file__) + '/martini_v2.2refP_constr.itp .', shell=True, check=True)
        #Read aatopol.top amd return all-atom structure as networkx Graph()
        aa_struct=aatop_2_cg.topol2graph(options.top) 
        # Returns aa_topol as networkx Graph()
        nx.write_gpickle(aa_struct,'aa_struct.pickle')
        # Read aa_struct , a networkx Graph(), write mapping scheme index file 
        # aa2cg.ndx and return CG structure as networkx DiGraph()
        cg_struct=aatop_2_cg.AA2CG(aa_struct,ndx=options.ndx) 
        nx.write_gpickle(cg_struct,'cg_struct.pickle')

        # Writes: unparam.pickle, param.pickle, dictionaries containing 
        # unparameterized and parameterized bonded distributions.
        aatop_2_cg.gen_unparam(cg_struct,params['cgff_curr'],
            dih_initial=params['dih_initial'])
        aatop_2_cg.write_e2e(cg_struct)
        # Generates AA bonded distriutions
        if not os.path.exists('CG1/'+options.outname):
            subprocess.run('bash ' + os.path.dirname(__file__) + 
                '/run_cg_sim.sh gen_aa_dist ' + options.file + ' ' +
                options.tpr_file + ' ' + str(options.begin) + ' ' + 
                str(options.end) + ' ' + os.getcwd(), shell=True, check=True)

            subprocess.run('rm -f bonded_distribution/#*', shell=True, 
            check=True)

            # Generate initial guess for unparameterized CG distributions
            subprocess.run('mkdir -p CG1', shell=True)
            subprocess.run('mkdir -p ref_param_fit', shell=True)

            aatop_2_cg.gen_new_cgff(cg_struct,params['cgff_curr'],'CG1/cgff1.pickle')
            aatop_2_cg.write_CGtopol(1,options.outname,cg_struct, None, params['peiname'])

    if remainder[0]=='parameterize':
        #run CG1
        cg_struct=nx.read_gpickle('cg_struct.pickle')
        if params['start_iter']==1:
            subprocess.run('bash ' + os.path.dirname(__file__) + '/run_cg_sim.sh ' +
                'run_CG_sim CG1 .. ' + str(options.begin) + ' ' + str(options.end) + 
                ' ' + os.path.dirname(__file__) + ' ' + str(params['gpu']) + ' ' +
                params['ions_mdp'] + ' ' +  params['em_mdp'] + ' ' + 
                params['npt_mdp'] + ' ' + params['md_mdp'], shell=True, check=True)
                
            if not os.path.exists('CG1/comparing_pngs/all_images.pdf'):
                report.gen_png(1, bond_small=params['bond_small'], 
                    bond_large=params['bond_large'], bond_ymax=params['bond_ymax'],
                    ang_ymax=params['ang_ymax'], dih_ymax=params['dih_ymax'])
                report.create_latex(1, fig_per_row=3)
                os.chdir('CG1/comparing_pngs')
                print('Writing result PDF')
                subprocess.run('pdflatex all_images.tex &>> output.log',shell=True, 
                check=True)
                os.chdir('../../')
            report.add_cost(params['wb'],params['wa'],'.','CG1',1)

        for i in range(max(2,params['start_iter']),params['max_iter']+1):
            subprocess.run('mkdir -p CG'+str(i), shell=True)
            subprocess.run('mkdir -p CG'+str(i-1)+'_th', shell=True)
            subprocess.run('mkdir -p CG'+str(i-1)+'_K', shell=True)
            # Run th step 
            run=0
            if not os.path.exists('CG' + str(i-1) + '_th/cgff' + str(i-1) + \
                '_th.pickle'): 
                run=aatop_2_cg.update_th_params(params['thc'], i-1, '.', 
                    params['cost_tol'])
            if run==0 and not os.path.exists('CG' +str(i-1)+ '_th/cgtopol.top'):
                aatop_2_cg.write_CGtopol(str(i-1) + '_th', 'cgtopol.top', 
                    cg_struct, None, PEINAME=params['peiname'])

            if run==0:
                subprocess.run('bash ' + os.path.dirname(__file__) +
                    '/run_cg_sim.sh run_CG_sim CG'+str(i-1)+'_th .. ' +
                    str(options.begin) + ' ' + str(options.end) + ' ' + 
                    os.path.dirname(__file__) + ' ' + str(params['gpu']) + ' ' +
                    params['ions_mdp'] + ' ' +  params['em_mdp'] + ' ' + 
                    params['npt_mdp'] + ' ' + params['md_mdp'], shell=True,
                    check=True)
            # Run K step
            if not os.path.exists('CG' + str(i-1) + '_K/cgff' + str(i-1) + \
                '_K.pickle'):
                run=aatop_2_cg.update_K_params(params['kc'], i-1, '.',
                    params['cost_tol'])
                #write topol, .itp
            if run==0 and not os.path.exists('CG' +str(i-1)+ '_K/cgtopol.top'):
                aatop_2_cg.write_CGtopol(str(i-1) + '_K', 'cgtopol.top',
                    cg_struct, None, PEINAME=params['peiname'])

            if run==0:
                subprocess.run('bash ' + os.path.dirname(__file__) + 
                    '/run_cg_sim.sh run_CG_sim CG'+str(i-1)+'_K .. ' +
                    str(options.begin) + ' ' + str(options.end)+ ' ' +
                    os.path.dirname(__file__) + ' ' + str(params['gpu']) + ' ' +
                    params['ions_mdp'] + ' ' +  params['em_mdp'] + ' ' + 
                    params['npt_mdp'] + ' ' + params['md_mdp'], shell=True, 
                    check=True)

            # Run new step
            if not os.path.exists('CG' + str(i) + '/cgff' + str(i) + '.pickle'):
                #update cgff
                run=aatop_2_cg.update_bonded_params(params['fa'], params['fd'], 
                        params['wb'], params['wa'], '.', i, params['cost_tol'])
            if run==-1:
                print('Tolerance Reached: Parameterization Complete')
                return
            elif not os.path.exists('CG' +str(i)+ '/cgtopol.top'):
                aatop_2_cg.write_CGtopol(i, 'cgtopol.top', cg_struct, None,
                        PEINAME=params['peiname'])
                #write topol, .itp
 
            subprocess.run('bash ' + os.path.dirname(__file__) + \
                '/run_cg_sim.sh run_CG_sim CG'+str(i)+' .. ' + \
                str(options.begin) + ' ' + str(options.end) + ' ' +
                os.path.dirname(__file__) + ' ' + str(params['gpu']) + ' ' +
                params['ions_mdp'] + ' ' +  params['em_mdp'] + ' ' + 
                params['npt_mdp'] + ' ' + params['md_mdp'], shell=True, 
                check=True)

            if not os.path.exists('CG'+str(i)+'/comparing_pngs/all_images.pdf'):
                os.system('mkdir -p CG'+str(i)+'/comparing_pngs')
                report.gen_png(i, bond_small=params['bond_small'], 
                    bond_large=params['bond_large'], bond_ymax=params['bond_ymax'],
                ang_ymax=params['ang_ymax'], dih_ymax=params['dih_ymax'])
                report.create_latex(i, fig_per_row=3)
                os.chdir('CG'+str(i)+'/comparing_pngs')
                print('Writing result PDF')
                subprocess.run('pdflatex all_images.tex &>> output.log',shell=True, 
                check=True)
                os.chdir('../../')
            report.add_cost(params['wa'],params['wb'],'.','CG'+str(i),i)

    if remainder[0]=='gen_cg_dist':
        subprocess.run('bash ' + os.path.dirname(__file__) + '/run_cg_sim.sh ' +
            'gen_cg_dist ' + str(options.begin) + ' ' + str(options.end) + ' ' + 
            options.dir, shell=True, check=True) 
       #Example: from the current CG directory.
       # coarsen gen_cg_dist -b [t0] -e [tn] -d [AA dir] 

    if remainder[0]=='add_cost':
        try:
            val=int(options.key)
        except:
            val=options.key
        report.add_cost(params['wb'], params['wa'], '.', options.file, val)
        #Example: coarsen add_cost -x parameters.dat -k 10 (or string) -f [dir]
             
    if remainder[0]=='smile2cg':
        subprocess.run('cp ' + os.path.dirname(__file__) + '/martini_v2.0_ions.itp .', shell=True, check=True)
        subprocess.run('cp ' + os.path.dirname(__file__) + '/martini_v2.2refP.itp .', shell=True, check=True)
        subprocess.run('cp ' + os.path.dirname(__file__) + '/martini_v2.2refP_constr.itp .', shell=True, check=True)

        cg_struct=smile.smile2top(remainder[1])
        nx.write_gpickle(cg_struct,'cg_struct.pickle')
        aatop_2_cg.gen_unparam(cg_struct,params['cgff_curr'],
            dih_initial=params['dih_initial'])
        aatop_2_cg.write_e2e(cg_struct)
        smile.gen_ini_cord(cg_struct,params['cgff_curr'], params['pos_prec'], remainder[1], options.outname)
        aatop_2_cg.write_CGtopol(-1,options.top,cg_struct,params['cgff_curr'],params['peiname'])

    if remainder[0]=='get_prop':
        if remainder[1]=='pei':
            if options.file is None:
                report.gen_mol_prop('cg_struct.pickle')
            else:
                report.gen_mol_prop(options.file)
    if remainder[0]=='cg_fac':
        smile.smile2cgfac(remainder[1])

if __name__=='__main__':
    main()