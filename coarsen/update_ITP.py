import os 
import sys
from add_cgtopol_dih import topol

def add_ffparam(filename,dist,typ,params):
    f=open(filename,'a') #Append to the new itp file
    if typ=='ang': #If force field parameter is for an angle
        for i in range(3):
            f.write(params[i]+'\t') #Write the angle i, j, k values
        #Funtion type = 2 \t Angle \t Force constant \t #Distribution name new 
        f.write('2\t'+str(topol[dist][0])+'\t'+str(topol[dist][1])+'\t#'+dist+' new\n')
    elif typ=='dih':#If force field parameter is for an dihedral
        for i in range(topol[dist][0]):
            for j in range(4):
                f.write(params[j]+'\t') #Write the dihedral i, j, k, l values
            #Funtion type = 1 \t Angle psi \t Force constant \t periodicity n \t #new 
            f.write('1\t'+str(topol[dist][1][i][0])+'\t'+str(topol[dist][1][i][1])+'\t'+str(topol[dist][1][i][2])+' ; new\n')
    f.close()




def update_itp(filename,idx,outpath):
    f=open(filename,'r') #Read old itp file
    molname=filename.split('.')[0]  #Get the molecule name for itp filename
    newitp=molname+'_'+str(idx) #Generate the name of new itp file by appending the iteration number
    w=open(outpath+newitp+'.itp','w') #Open new itp file for writing
    cnt=0
    for lines in f:
        foo=lines.split()
        if len(foo)<2: #If empty lines, [] directivels  
            w.write(lines) 
            continue
        elif foo[0]==molname:
            w.write(newitp+'\t'+foo[1]+'\n') #Change molecule name based based on the iteration
        elif foo[1]=='angles' or foo[1]=='dihedrals': #Change cnt values to which is used to insert bonded parameters
            cnt+=1
            w.write(lines)
        elif foo[cnt+2]==';Determine': #If bonded parameters is not determined
            params=foo[0:cnt+2] # Idx for bonded parameters. for bonds i,j ; for angles i,j,k ; for dihedrals i,j,k,l 
            dist_name=foo[-1] #save the name of bonded distribution
            w.close()
            if cnt==1:#Add the parameters for angles
                add_ffparam(outpath+newitp+'.itp',dist_name,'ang',params)
            elif cnt==2:#Add the parameters for dihedrals
                add_ffparam(outpath+newitp+'.itp',dist_name,'dih',params)
            w=open(outpath+newitp+'.itp','a')
        else:
            w.write(lines)
    w.close()
if not os.path.exists(sys.argv[3]):
    os.mkdir(sys.argv[3])
update_itp(sys.argv[1],sys.argv[2],sys.argv[3]) #Argv[1]=old itp file Argv[2]=iteration number Argv[3]=Iteration folder location
#Example: python update_ITP.py bPEI2K.itp 1 './CG1/'


