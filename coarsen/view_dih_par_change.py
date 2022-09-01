import sys
sys.path.insert(1,'./common/')
from update_bonded_params import dih_name_to_param, dih_param_to_fourier
import importlib
import copy
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('TkAgg')
from matplotlib import cm
import numpy as np

def view_dih_par_change(dih_name,I,J):
    params=[]
    for i in range(I,J+1):
        sys.path.insert(1,'./CG'+str(i)+'/')
        foo=importlib.import_module('add_cgtopol_dih'+str(i))
        top=copy.deepcopy(foo.topol)
        Kd,phi,n_dih=dih_name_to_param(dih_name,top)
        dih_fourier=dih_param_to_fourier(Kd,phi)
        params.append(dih_fourier)
    params=np.array(params)
    x=np.arange(len(params))
    #print(params)
    #print(params[:,1])
    for i in range(len(params[0])):
        plt.figure()
        plt.plot(x,params[:,i],color=cm.jet(i))
        plt.show()

view_dih_par_change('tttt',1,8)
