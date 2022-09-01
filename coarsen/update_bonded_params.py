import sys
import numpy as np
import os
import importlib
import copy
from gen_bonded_params import * 


def gen_list2(filename,typ):
    f=open(filename,'r')
    arr=[]
    if typ=='dih':
        vals=[]
    for lines in f:
        foo=lines.split()
        if len(foo)<1:
            continue
        arr.append(typ+'_'+foo[0])
        if typ=='dih':
            vals.append([float(foo[1]),int(foo[2])])
    if typ=='dih':
        return arr,vals
    else:
        return arr
        
def mean_pdf(data):
    return np.sum(np.multiply(data[:,0],data[:,1]))/np.sum(data[:,1])

def var_pdf(data):
    return np.sum(np.multiply(np.square(data[:,0]),data[:,1]))/np.sum(data[:,1])-mean_pdf(data)**2

def std_pdf(data):
    return np.sqrt(var_pdf(data))

def dih_name_to_param(dih_name,cgff):
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

def update_bonded_params(thc,kc,kdc,fa,fd): #REQUIRES THE FILES add_cgtopol_*_K.py and dd_cgtopol_*_th.py files to be copied to the current directory
#K for bond angles changes at fa*0.1
    ang_list=gen_list2('../unparam_ang.dat','angle')
    dih_list,vals=gen_list2('../unparam_dih.dat','dih')
    if not os.path.exists('../CG'+str(N+1)):
        os.mkdir('../CG'+str(N+1))
    outname_ang='../CG'+str(N+1)+'/add_cgtopol_ang'+str(N+1)+'.py'
    w=open(outname_ang,'w')
    w.write("import sys\nimport os\npath_loc=os.path.abspath('../')\nsys.path.insert(1,path_loc)\nfrom cg_topol_vals import topol\n")
    w.close()

    outname_dih='../CG'+str(N+1)+'/add_cgtopol_dih'+str(N+1)+'.py'
    w=open(outname_dih,'w')
    w.write("from add_cgtopol_ang"+str(N+1)+" import topol\n")
    w.close()

    #Find Jacobian for Angle parameters 
    #w=1 #Difference in standard deviation is scaled by a factor of w to make it less important for distrbution matching 
    for i in range(len(ang_list)):
        #  [  d(e_meu)    d(e_meu)  ]
        #  [   d(th)        d(K)    ]
        #J=[                        ]
        #  [ d(e_sig)    d(e_sig)   ]
        #  [  d(th)        d(K)     ]    a 2x2 matrix

        # e_meu=meu_CG-meu_AA
        # e_sig=var_CG/var_AA-1

        J=np.zeros((2,2))
        ang_name=ang_list[i].split('_')[1]
        #Read bonded distributions
        data_aa=import_data('../bonded_distribution/'+ang_list[i]+'.xvg')
        data_cg=import_data('bonded_distribution/'+ang_list[i]+'.xvg')
        data_cg_th=import_data('../CG'+str(N)+'_th/bonded_distribution/'+ang_list[i]+'.xvg')
        data_cg_K=import_data('../CG'+str(N)+'_K/bonded_distribution/'+ang_list[i]+'.xvg')
       
        #print(ang_name)
        #Calculate mean of bonded distributions
        mean_aa=mean_pdf(data_aa[0],data_aa[1])
        mean_cg=mean_pdf(data_cg[0],data_cg[1])
        mean_cg_th=mean_pdf(data_cg_th[0],data_cg_th[1])
        mean_cg_K=mean_pdf(data_cg_K[0],data_cg_K[1])
        #print('mean',mean_aa,mean_cg,mean_cg_th,mean_cg_K)
        #Calculate standard deviation of bonded distributions
        var_aa=var_pdf(data_aa[0],data_aa[1])
        var_cg=var_pdf(data_cg[0],data_cg[1])
        var_cg_th=var_pdf(data_cg_th[0],data_cg_th[1])
        var_cg_K=var_pdf(data_cg_K[0],data_cg_K[1])
        #print('var',var_aa,var_cg,var_cg_th,var_cg_K)
        #Determine the jacobian
        dth=topol_th[ang_name][0]-topol_cur[ang_name][0]
        dK=topol_K[ang_name][1]-topol_cur[ang_name][1]
        #print('dth, dK',dth,dK)
        J[0,0]=np.divide(mean_cg_th-mean_cg,dth)
        J[0,1]=np.divide(mean_cg_K-mean_cg,dK)
        J[1,0]=np.divide(var_cg_th-var_cg,dth)
        J[1,1]=np.divide(var_cg_K-var_cg,dK)
        #d_param=[-np.divide(mean_cg-mean_aa,J[0,0]),-np.divide(mean_cg-mean_aa,J[0,1])]
        d_param=[-(mean_cg-mean_aa)/J[0,0],-(var_cg-var_aa)/J[1,1]]

        if np.abs(mean_cg-mean_aa)>2:
            if (topol_cur[ang_name][0]-mean_cg)*(mean_aa-mean_cg)>0:#K must increase
                if d_param[1]<0:
                     d_param[1]=1000

        if abs(mean_cg-mean_aa)>2:
            th=fa*d_param[0]+topol_cur[ang_name][0]
            if th<mean_aa-30:
                th=mean_aa-30
            if th>180:
                th=180
            Ka=fa*0.1*d_param[1]+topol_cur[ang_name][1]
            if Ka<5:
                Ka=5
        else: 
            th=topol_cur[ang_name][0]
            Ka=topol_cur[ang_name][1]
            #print(ang_name,int(th),int(Ka))
        write_ang_topol(outname_ang,ang_list[i],int(Ka),int(th))  #CHECK if the parameters actually change with this approximation
     
     # Find updated dihedral angle parameters 
     # {a} is parameters used for the i^th CG simulations.i.e. {a}={a_1, ... , a_j, ... a_n} 
     # {b_cg} is the parameters obtaiend from fitting the CG bonded distributions 
     # {b_cg+Δ} is the parameters obtaiend from fitting the CG_K bonded distributions 
     # {b_cg-Δ} is the parameters obtaiend from fitting the CG_th bonded distributions 
     # {b_aa} is the parameters obtaiend from fitting the AA bonded distributions 
     # {e} is {b_aa}-{b_cg}.
     # {e}'={e} + [J]{Δa} ; [J] is approximated to be a diagonal function
     # For {e}'=0 we get, [J]{Δa}=-{e}
     # J_(j,j) Δa_j = (b_cg)_j - (b_aa)_j
     # Δa_j = ( (b_cg)_j - (b_aa)_j )/de_j/da_j = -2Δ( (b_cg)_j - (b_aa)_j )/( (b_cg+Δ)_j - (b_cg-Δ)_j )
     #
    for i in range(len(dih_list)):
        
        dih_name=dih_list[i].split('_')[1]
        if dih_name[0]=='N':
            Kd,phi,n=dih_name_to_param(dih_name,topol_cur)
            write_dih_topol(outname_dih,dih_list[i],Kd,phi)

        else:
            #Get bonded distritions
            data_aa=import_data('../bonded_distribution/'+dih_list[i]+'.xvg')
            data_cg=import_data('bonded_distribution/'+dih_list[i]+'.xvg')
            data_cg_th=import_data('../CG'+str(N)+'_th/bonded_distribution/'+dih_list[i]+'.xvg') #Th is not changed K is changed by a negative amount
            data_cg_K=import_data('../CG'+str(N)+'_K/bonded_distribution/'+dih_list[i]+'.xvg')

            #Fill the bonded distributions so that all the distributions are from -180 to 180 degree
            xaa,yaa=fill_data(np.array(data_aa[0]),np.array(data_aa[1]))
            xcg,ycg=fill_data(np.array(data_cg[0]),np.array(data_cg[1]))
            xcg_th,ycg_th=fill_data(np.array(data_cg_th[0]),np.array(data_cg_th[1]))
            xcg_K,ycg_K=fill_data(np.array(data_cg_K[0]),np.array(data_cg_K[1]))

            #Get potential energy associated with bonded distributions 
            Uaa=-kbt*np.log(np.array(yaa)+tol) 
            Ucg=-kbt*np.log(np.array(ycg)+tol) 
            Ucg_th=-kbt*np.log(np.array(ycg_th)+tol) 
            Ucg_K=-kbt*np.log(np.array(ycg_K)+tol) 
            dU1=np.subtract(Ucg,Uaa) 
            dU2=np.subtract(Ucg_K,Ucg_th)
            
            p0=[1]*(topol_cur[dih_name][0]+1)
            popt1,pcov1=scipy.optimize.curve_fit(fourier_series,xaa,dU1,p0)
            popt2,pcov2=scipy.optimize.curve_fit(fourier_series,xaa,dU2,p0)

            Kd_old,phi_old,n_old=dih_name_to_param(dih_name,topol_cur)
            popt_old=dih_param_to_fourier(Kd_old,phi_old) 
            # CHECK if CG2 using jacobian and CG2 using just difference (test_new4) simulations is in the same direction

            d_popt=-2*kdc*np.divide(popt1,popt2)
            popt_new=np.zeros(len(d_popt))
            for j in range(len(d_popt)):
                if abs(popt1[j])>1E-4:
                    popt_new[j]=popt_old[j]+fd*d_popt[j]
                else:
                    popt_new[j]=popt_old[j]
            Kd,phi=conv_dih_par(popt_new)
#            print(dih_list[i]+' written')
            write_dih_topol(outname_dih,dih_list[i],Kd,phi)

def update_jack_params(thc,kc,kdc):
    ang_list=gen_list2('../unparam_ang.dat','angle')
    dih_list,vals=gen_list2('../unparam_dih.dat','dih')
    
    #Increase equilibrium bond angles by thc, and decrease dihedral parameters angles (in fourier series form) by kdc
    ############## Create files to write the force field #######################
    if not os.path.exists('../CG'+str(N)+'_th'):
        os.mkdir('../CG'+str(N)+'_th')
    outname_ang='../CG'+str(N)+'_th/add_cgtopol_ang'+str(N)+'_th.py'
    w=open(outname_ang,'w')
    w.write("import sys\nimport os\npath_loc=os.path.abspath('../')\nsys.path.insert(1,path_loc)\nfrom cg_topol_vals import topol\n")
    w.close()

    outname_dih='../CG'+str(N)+'_th/add_cgtopol_dih'+str(N)+'_th.py'
    w=open(outname_dih,'w')
    w.write("from add_cgtopol_ang"+str(N)+"_th import topol\n")
    w.close()
    ##############################################################################

    for i in range(len(ang_list)):
        ang_name=ang_list[i].split('_')[1]
        #Change ang parameters only for th.
        update_ang=topol_cur[ang_name][0]+thc
        if update_ang > 180:
            update_ang=180
        write_ang_topol(outname_ang,ang_list[i],topol_cur[ang_name][1],update_ang)
    for i in range(len(dih_list)):
        dih_name=dih_list[i].split('_')[1]
        Kd,phi,n_dih=dih_name_to_param(dih_name,topol_cur)
        if dih_name[0]!='N':
            #Convert to fourier series for small change in the forcefield parameters
            fourier_param=dih_param_to_fourier(Kd,phi)
            fourier_param=np.array(fourier_param)-kdc
            #Convert back to Gromacs diherdral potential
            Kd,phi=conv_dih_par(fourier_param)
        write_dih_topol(outname_dih,dih_list[i],Kd,phi)
    
    #################Increase bond angle force constants by kc and Increase dihedral parameters (in fourier form) by kdc #################
    ########create files for forcefield definitions###########
    if not os.path.exists('../CG'+str(N)+'_K'):
        os.mkdir('../CG'+str(N)+'_K')
    outname_ang='../CG'+str(N)+'_K/add_cgtopol_ang'+str(N)+'_K.py'
    w=open(outname_ang,'w')
    w.write("import sys\nimport os\npath_loc=os.path.abspath('../')\nsys.path.insert(1,path_loc)\nfrom cg_topol_vals import topol\n")
    w.close()

    outname_dih='../CG'+str(N)+'_K/add_cgtopol_dih'+str(N)+'_K.py'
    w=open(outname_dih,'w')
    w.write("from add_cgtopol_ang"+str(N)+"_K import topol\n")
    w.close()
    ########################################################33
    for i in range(len(ang_list)):
        ang_name=ang_list[i].split('_')[1]
        #Change ang parameters only for K.
        write_ang_topol(outname_ang,ang_list[i],topol_cur[ang_name][1]+kc,topol_cur[ang_name][0])
    for i in range(len(dih_list)):
        dih_name=dih_list[i].split('_')[1]
        Kd,phi,n_dih=dih_name_to_param(dih_name,topol_cur)
        if dih_name[0]!='N':
            #######Increase parameters in fourier series form by kdc
            fourier_param=dih_param_to_fourier(Kd,phi)
            fourier_param=np.array(fourier_param)+kdc
            #######Convert back to GROMACS dihedral potential
            Kd,phi=conv_dih_par(fourier_param)
        write_dih_topol(outname_dih,dih_list[i],Kd,phi)


if __name__ == '__main__':
    state=sys.argv[1]
    cur_dir=os.getcwd()
    cur_folder=cur_dir.split('/')
    cur_folder=cur_folder[-1]
    if '_' in cur_folder:#Either in CG*_th or CG*_K
        foo=cur_folder.split('_')
        N=int(foo[0][2:]) #Iteration number 
        jac_type=foo[1] #Theta is changed or K is changed
    else: #In CG iteration not for jacobian calculation
        N=int(cur_folder[2:]) #Iteration number
        
    foo=importlib.import_module('add_cgtopol_dih'+str(N))
    topol_cur=copy.deepcopy(foo.topol)
    if state == 'update':
        sys.path.insert(1,'../CG'+str(N)+'_th/')
        foo1=importlib.import_module('add_cgtopol_dih'+str(N)+'_th')
        topol_th=copy.deepcopy(foo1.topol)
    
        sys.path.insert(1,'../CG'+str(N)+'_K/')
        foo2=importlib.import_module('add_cgtopol_dih'+str(N)+'_K')
        topol_K=copy.deepcopy(foo2.topol)

        update_bonded_params(5,50,0.0001,0.5,0.5)
    elif state=='jacobian': #In CG iteration not for jacobian calculation
        update_jack_params(5,50,0.0001)
    else:
        print('To generate jacobian type : python '+sys.argv[0]+' jacobian\n To run the next iteration type: python '+sys.argv[0]+' update')
    
