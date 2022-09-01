#Create a latex file containing all the images:
import glob
num_per_fig=3

def gen_list2(filename,typ):
    f=open(filename,'r')
    arr=[]
    if typ=='dih':
        vals=[]
    for lines in f:
        foo=lines.split()
        arr.append(foo[0])
        if typ=='dih':
            vals.append([float(foo[1]),int(foo[2])])
    if typ=='dih':
        return arr,vals
    else:
        return arr

def get_list(filename):
    f=open(filename,'r')
    a=[]
    for lines in f:
        foo=lines[0:-1]
        a.append(foo)
    return a

w=open('all_images.tex','w')

w.write('\\documentclass{article}\n\n\\usepackage{graphicx}\n\\usepackage[left=2cm, top=2cm, bottom=2cm, right=2cm]{geometry}\n\\usepackage{caption}\n\\usepackage{xcolor}\n\\usepackage{subcaption}\n\\usepackage{float}\n\n\\begin{document}\n\n')

w.write('\\section{Bonds}\n\n')
files=glob.glob('bond_*.png')
cnt=0
for f in files:
    if cnt%num_per_fig==0:
        w.write('\t\\begin{figure}[H]\n\t\t\\centering\n')
    dist_name=f.split('.')[0].split('_')[1]
    w.write('\t\t\\begin{subfigure}{0.3\\textwidth}\n\t\t\t\\includegraphics[width=\\textwidth]{'+f+'}\n\t\t\t\\caption{'+dist_name+'}\n\t\t\t\\label{fig:'+dist_name+'}\n\t\t\\end{subfigure}\n')
    if cnt%num_per_fig==num_per_fig-1 or cnt==len(files)-1:
        w.write('\t\\end{figure}\n\n')
    cnt+=1

w.write('\\section{Angles}\n\n')
files=glob.glob('angle_*.png')
ang_list=get_list('../../unparam_ang.dat')
print(ang_list)
cnt=0
for f in files:
    if cnt%num_per_fig==0:
        w.write('\t\\begin{figure}[H]\n\t\t\\centering\n')
    dist_name=f.split('.')[0].split('_')[1]

    w.write('\t\t\\begin{subfigure}{0.3\\textwidth}\n')
    if dist_name in ang_list:
        w.write('\t\t\t\\captionsetup{font={color=red,bf}}\n')
    w.write('\t\t\t\\includegraphics[width=\\textwidth]{'+f+'}\n\t\t\t\\caption{'+dist_name+'}\n\t\t\t\\label{fig:'+dist_name+'}\n\t\t\\end{subfigure}\n')
    if cnt%num_per_fig==num_per_fig-1 or cnt==len(files)-1:
        w.write('\t\\end{figure}\n\n')
    cnt+=1

w.write('\\section{Dihedrals}\n\n')
files=glob.glob('dih_*.png')
dih_list,vals=gen_list2('../../unparam_dih.dat','dih')
print(dih_list)
cnt=0
for f in files:
    if cnt%num_per_fig==0:
        w.write('\t\\begin{figure}[H]\n\t\t\\centering\n')
    dist_name=f.split('.')[0].split('_')[1]
    w.write('\t\t\\begin{subfigure}{0.3\\textwidth}\n')
    if dist_name in dih_list:
        w.write('\t\t\t\\captionsetup{font={color=red,bf}}\n')
    w.write('\t\t\t\\includegraphics[width=\\textwidth]{'+f+'}\n\t\t\t\\caption{'+dist_name+'}\n\t\t\t\\label{fig:'+dist_name+'}\n\t\t\\end{subfigure}\n')
    if cnt%num_per_fig==num_per_fig-1 or cnt==len(files)-1:
        w.write('\t\\end{figure}\n\n')
    cnt+=1
w.write('\\end{document}')

