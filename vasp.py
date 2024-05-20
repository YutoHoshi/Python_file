import os,sys

if __name__=="__main__":
    name='Pb'
    vca=[]
    core_in_f=1
    ispin=1
    voskwn=0 #0 considers magnetizm, 1 doesnot
    slab=0
    mag=[]
    soc=0
    ladu=0
    gga="PS"
    umat=[]#[[-1,0,0],[2,3,0.3],[-1,0,0],[-1,0,0]]
#    exc='paw_PBE_52'
    knum_scf=[12,12,12]
    knum_band=[20,0,0]
    encut=600  #550eV corresponds about 40Ry
    ediff=0.00001
    nbands=120
    nelm=200
    nelmin=2
    nedos=20
    lorbit=0 #0,1=none,2=use rwigs(no phase), 10 or 11=not use rwigs(+phase in 11)
    ismear=1 #-5=line. tetrahedron
    sigma=0.01
    nedos=20
    par=[8,4] #kpar,npa
    point=[['$\Gamma$',0,0,0],['X',0.0,0.0,0.5],['N',0.5,0,0],['$\Gamma$',0,0,0],['Z',0.25,0.25,-0.25]]
    with open("kpoint_forme.py",mode="w") as kt:
        kt.write(f'point = {point}\n')
    isym=1
#   optimization
    pstress=0
    isif=3
    nsw = 200
    ibrion=2
    ediffg = -0.05
    opt=[pstress,isif,nsw,ibrion,ediffg]
#wannier90 inputs
    wpara=[8,0.2,-0.2,-8] #winmax, frozmax, frozmin,winmin
    knum_wan=[12,12,12]
    dis_num_iter=1000
    num_iter=300
    dis_mix_ratio=0.4
    write_FS=0 # 0=do not plot Fermi Surface, 1=do
    FS_num=50 #the default value is 50
    wan_orb=0
    iprint=1  #0 to 3, minuteness of information in output file. default is 1

    pdpos=[0.25,0.341997,0.143442,0.75,0.658003,0.856558, 0.25,0.1434416,0.341997,0.75,0.856558,0.658003]
    prj=[
        ['Pd1:','f='+str(pdpos[0])+','+str(pdpos[1])+','+str(pdpos[2])+',','dx2-y2'],
        ['Pd1:','f='+str(pdpos[0])+','+str(pdpos[1])+','+str(pdpos[2])+',','dz2'],
        ['Pd2:','f='+str(pdpos[3])+','+str(pdpos[4])+','+str(pdpos[5])+',','dx2-y2'],
        ['Pd2:','f='+str(pdpos[3])+','+str(pdpos[4])+','+str(pdpos[5])+',','dz2'],
        ['Pd3:','f='+str(pdpos[6])+','+str(pdpos[7])+','+str(pdpos[8])+',','dx2-y2'],
        ['Pd3:','f='+str(pdpos[6])+','+str(pdpos[7])+','+str(pdpos[8])+',','dz2'],
        ['Pd4:','f='+str(pdpos[9])+','+str(pdpos[10])+','+str(pdpos[11])+',','dx2-y2'],
        ['Pd4:','f='+str(pdpos[9])+','+str(pdpos[10])+','+str(pdpos[11])+',','dz2']
    ]

    nwan=len(prj) 
    ef=0.0
    if os.path.isfile("OUTCAR"):
        with open("OUTCAR") as outcar:
            otlines=outcar.readlines()
            for ot in otlines:
                if "Fermi energy:" in ot:
                    dd=ot.split()
                    ef=float(dd[2])
                    print("Fermi Energy is "+str(ef)+" eV")
                    break
    else:
        print("OUTCAR does not exists: we set the Fermi energy as ZERO")

        
#-----------------------------


import os,sys

def read_poscar(flag):
    tmp=0; latvec=[]
    fposcar=open('POSCAR','r')
    fposcar.readline() #fake
    fposcar.readline() #fake
    for i in range(3):  #read lattice vector as 3*3 matrix
        tmp=fposcar.readline().split()
        latvec.append([float(tmp[0]),float(tmp[1]),float(tmp[2])])
    tmp=fposcar.readline().split() #store the kinds of atom
    atom=[] 
    for i in range(len(tmp)):
        atom.append(tmp[i])
    if flag=='POSIT': #read also atomic position
        tmp=fposcar.readline().split() #read the number of each atom
        atom_num=[]
        anum=0
        for i in range(len(tmp)): # store the number of different atoms
            atom_num.append(int(tmp[i]))  
        fposcar.readline() #fake
        atom_list=[] 
        for i in range(len(atom)): # loop of atomic kind 
            for j in range(atom_num[i]): # loop of same atom
                tmp=fposcar.readline().split()
                atom_list.append([atom[i],float(tmp[0]),float(tmp[1]),float(tmp[2])]) # store the atomic positions
        fposcar.close()
        atom_name_num=[]
        for i in range(len(atom)):
            atom_name_num.append([atom[i],atom_num[i]])
        return(latvec,atom_list,atom_name_num)
    if flag=='SYMBOL':
        fposcar.close()
        return(atom)
    else:
        print('ERROR IN read_poscar!')
        sys.exit(1)

def get_rwings():
    rings=[]
    for f in open('POTCAR'):
        if(f.find('RWIGS')!=-1):
            tmp=f.split()
            rwings.append(float(tmp[5])) #please check by typing "grep RWIGS POTCAR"
    return(rwigs)

#def make_potcar(exc,core_in_f):
def make_potcar(core_in_f):
    def coref(atom):
        f3_list=('Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Lu')
        f2_list=('Eu','Er','Yb')
        for i in range(len(atom)):
            for j in range(len(f3_list)):
                if(atom[i]==f3_list[j]):
                    atom[i]=atom[i]+'_3'
        return(0)
    def svpv(atom,pv):
        sv_list=['Ba','Ca','Cs','K','Nb','Ra','Rb','Sr','Y','Zr']
        pv_list=['Ca','K','Nb','Rb']
        for i in range(len(atom)):
            if pv==1:
                for j in range(len(pv_list)):
                    if atom[i]==pv_list[j]:
                        atom[i]=atom[i]+'_pv'
            for j in range(len(sv_list)):
                    if atom[i]==sv_list[j]:
                        atom[i]=atom[i]+'_sv'
        return(0)
    mssg='begin makeing POTCAR'
    atom=read_poscar('SYMBOL')  #use only atomic symbols
    mssg=mssg+'\n this compound consists of \n'
    for i in range(len(atom)):
        mssg=mssg+atom[i]+''
        info=svpv(atom,1)
        if core_in_f==1:
            mssg=mssg+'\n f-orbitals are kept inside core'
            info=coref(atom)

    hm=os.environ['HOME']
    pot_link=hm+"/vasppot"
    dum=''
    for i in range(len(atom)):
        dum=dum+'%s/%s/POTCAR '%(pot_link,atom[i])
    dum=dum+'> POTCAR'
    os.system('cat %s'%dum)
    mssg=mssg+'\n POSCAR has been made'
    print(mssg)

def make_incar(projectname,encut,istart,icharg,ispin,mag,nelm,nelmin,ediff,nbands,lorbit,nedos,ismear,sigma,
               voskwn,soc,ldau,umat,slab,vca,flag_wan,par,opt,isym,ef,
               prj,wpara,sympoint,dis_num_iter,num_iter,dis_mix_ratio,iprint,wan_orb):
    f=open('INCAR','w')
    f.write('''SYSTEM   =%s
PREC     =Normal
ENCUT    =%d
ISTART   =%d
ICHARG   =%d
ISPIN    =%d
NELM     =%d
NELMIN   =%d
EDIFF    =%f
NBANDS   =%d
LORBIT   =%d
NEDOS    =%d
ISMEAR   =%d
SIGMA    =%5.3f
WOSKWN   =%d
KPAR     =%d
NPAR     =%d
GGA      =%s
    \n'''%(projectname,encut,istart,icharg,ispin,nelm,nelmin,ediff,nbands,lorbit,nedos,ismear,sigma,voskwn,par[0],par[1],gga))

    if(lorbit==2):
        rwigs=get_rwigs()
        f.write('RWIGS    =')
        for i in range(len(rwigs)):
            f.write(' %f'%rwigs[i])
        f.write('\n')

    if(len(opt)!=0 and isif!=0 ):
        f.write('''PSTRESS  =%d
ISIF     =%d
NSW      =%d
IBRION   =%d
EDIFFG   =%f
ISYM     =%d\n'''%(pstress,isif,nsw,ibrion,ediffg,isym))
        
    if(len(mag)!=0 and soc==1):
        f.write("MAGMON   = ")
        for i in range(len(mag)):
            f.write(" %d %d %d "%(mag[i][0],mag[i][1],mag[i][2]))
        f.write("\n")
        
    if(len(mag)!=0 and soc==0):
        f.write("MAGMON   = ")
        for i in range(len(mag)):
            f.write(" %s "%(mag[i]))
        f.write("\n")

    if(soc==1):
        f.write("ISYM     =-1\n")
        
    if(ldau==1):
        f.write("""LMAXMIX  = 4\nLDAU     =.TRUE.  \nLDAUL    =""")
        if(len(umat)==0):
            print('ERROR IN LDA+U INPUTS!')
            sys.exit(1)
        for i in range(len(umat)):
            f.write(" %d "%(umat[i][0]))
        f.write("\nLDAUU    =")
        for i in range(len(umat)):
            f.write(" %f "%(umat[i][1]))
        f.write("\nLDAUJ    =")
        for i in range(len(umat)):
            f.write(" %f "%(umat[i][2]))
        f.write("\n")

    if(len(vca)!=0):
        f.write('VCA    =')
        for i in range(len(vca)):
            f.write(' %f '%vca[i])
        f.write('\n')
    if(flag_wan==1):
#        f.write("""LWANNIER90      = .TRUE.\nLWRITE_MMN_AMN  = .TRUE.\nLWANNIER90_RUN  = .TRUE.\n""")
        f.write("""LWANNIER90      = .TRUE.
        LWRITE_MMN_AMN  = .TRUE. \n""") #WANNIER90_RUN  = .TRUE.\n""")
        f.write("""NUM_WANN = %d \n"""%len(prj))
########################begin wannier90.win replacement##########################
        f.write('''WANNIER90_WIN ="
iprint = %d
!write_proj   = .true.
write_hr   =.true.
write_rmn  =.true.
write_tb   =.true.
write_xyz    = .true.
translate_home_cell = .true.
dis_win_max  = %f
dis_froz_max = %f
dis_froz_min = %f
dis_win_min  = %f
dis_num_iter = %d
num_iter     = %d
dis_mix_ratio= %f\n
        '''%(iprint,(ef+wpara[0]),(ef+wpara[1]),(ef+wpara[2]),(ef+wpara[3]),dis_num_iter,num_iter,dis_mix_ratio))

        if(wan_orb==1):
            f.write('wannier_plot=.True.\nwannier_plot_supercell=4\n\n')

        f.write('''
use_bloch_phases = %s
fermi_energy     = %15.10f
bands_plot       = %s\n'''%('.false.',ef,'.true.'))
        if(write_FS==1):  #plot Fermi surface
            f.write('fermi_surface_plot = .true. \n')
        if(FS_num!=0):
            f.write('fermi_surface_num_points = %d\n'%(FS_num))

        f.write('Begin projections\n') #specify the wannier projection funciton
        for i in range(len(prj)):
            f.write('%s: %s\n'%(prj[i][1],prj[i][2]))
        f.write('End projections\n\n')                

        f.write('''Begin Kpoint_Path\n''')
        print(sympoint)
        for i in range(len(sympoint)-1):
            if sympoint[i][0]=='GAMMA':
                sympoint[i][0]=='G'
            if sympoint[i+1][0]=='GAMMA':
                sympoint[i+1][0]=='G'
            f.write('%7s %7.4f %7.4f %7.4f %7s  %7.4f %7.4f %7.4f \n'
                    %(sympoint[i][0],sympoint[i][1],sympoint[i][2],sympoint[i][3],
                      sympoint[i+1][0],sympoint[i+1][1],sympoint[i+1][2],sympoint[i+1][3]))
        f.write('End Kpoint_Path\n\n')
########################end wannier90.win replacement##########################
        f.write('''"\n''')
        if(wan_orb==1):
            f.write('LWRITE_UNK      = .TRUE.')
    if(flag_wan==-1):
        f.write("""LWANNIER90      = .TRUE.\nLWRITE_MMN_AMN  = .FALSE.\nLWANNIER90_RUN  = .FALSE.\n""")

        
    f.close()
def make_kpoints(kpoint,sympoint,line):
    f=open('KPOINTS','w')
    f.write('KPOINTS\n')
    if line==1:
        f.write('%d\nLine-mode\nRec\n'%kpoint[0])
        for i in range(len(sympoint)-1):
            f.write('  %f %f %f ! %s\n'%(sympoint[i][1],sympoint[i][2],sympoint[i][3],sympoint[i][0]))
            f.write('  %f %f %f ! %s\n'%(sympoint[i+1][1],sympoint[i+1][2],sympoint[i+1][3],sympoint[i+1][0]))
            if i!=len(sympoint)-2:
                f.write('\n')
    else:
        f.write('0\nGamma\n')
        f.write(' %d %d %d\n 0 0 0\n'%(kpoint[0],kpoint[1],kpoint[2]))
    f.close()        

def main(nelm,nelmin,nedos,knum_scf):
    flag_scf=0; flag_dos=0; flag_band=0; flag_wan=0; flag_potcar=0; flag_win=0
    if(len(sys.argv)>1):
        for i in range(1,len(sys.argv)):
            if(sys.argv[i].find('-scf')!=-1):
                flag_scf=1
            elif(sys.argv[i].find('-band')!=-1):
                flag_band=1
            elif(sys.argv[i].find('-wannier')!=-1 or sys.argv[i].find('-wan')!=-1):
                flag_wan=1
            elif(sys.argv[i].find('-dos')!=-1):
                flag_dos=1
            elif(sys.argv[i].find('-pot')!=-1):
                flag_potcar=1
            elif(sys.argv[i].find('-make_win')!=-1):
                flag_win=1
                
    if(flag_potcar==1):
#        make_potcar(exc,core_in_f)
        make_potcar(core_in_f)
        dum1=1; dum2=11 ; wan=0; line=0
        
    if(flag_win==1):
        dum1=0; dum2=2 ; wan=0; line=0
        wan=-1
        nelm=1
        nelmin=0
        nedos=20
        make_kpoints(knum_wan,point,line)
        for i in range(len(knum_scf)):
            knum_scf[i]=knum_wan[i]

    dumsmear=ismear
    if(flag_scf==1):
        dum1=1; dum2=1 ; wan=0; line=0
        make_kpoints(knum_scf,point,line)
    if(flag_dos==1):
        dum1=1; dum2=11 ; wan=0; line=0
    if(flag_band==1):
        dum1=1; dum2=11 ; wan=0; line=1; dumsmear=1
        make_kpoints(knum_band,point,line)
    if(flag_wan==1):
        dum1=1; dum2=11 ; wan=1; line=0
        print(point)
        make_kpoints(knum_wan,point,line)

    make_incar(name,encut,dum1,dum2,ispin,mag,nelm,nelmin,ediff,nbands,lorbit,nedos,dumsmear,sigma,
               voskwn,soc,ladu,umat,slab,vca,wan,par,opt,isym,ef,
               prj,wpara,point,dis_num_iter,num_iter,dis_mix_ratio,iprint,wan_orb)
   
if __name__=="__main__":
    main(nelm,nelmin,nedos,knum_scf)
        

