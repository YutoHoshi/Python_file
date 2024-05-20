import numpy as np
import matplotlib.pyplot as plt
import math
import numpy.linalg as LA
from kpoint_forme import point
import os, sys, copy

N = 1 #原子数
num_points = 80 #メッシュ(k点座標間)
rprm = 4 #ホッピングは小数第3位まで
rcut = [50, 50, 50] #カットオフ
figname="vasp-original"
os.makedirs(figname, exist_ok=True)
os.makedirs("flexin", exist_ok=True)
flexhop="flexin/ham2.flex"
flexvec="flexin/rvec.flex"


#E_f設定
if len(sys.argv) == 2:#コマンドライン引数でE_fを指定
    ef=float(sys.argv[1])
    print("The Fermi energy is given via command line.")
elif os.path.isfile("scf/OUTCAR"):
    with open("scf/OUTCAR",mode="r") as outcar:
        for line in outcar:
            if "Fermi energy" in line:
                dd = line.split()
                ef = float(dd[2])
                print("The Fermi energy is read from OUTCAR file.")
                break
else:
    ef=0
    print("no input for Fermi energy: ef=0")

print("Fermi energy is ",ef)


#基本並進ベクトルa,逆基本並進ベクトルbを取得
a=np.zeros((3,3))#, dtype = np.complex128)
b=np.zeros((3,3))#, dtype = np.complex128)
with open("scf/OUTCAR") as out:
    outcar=out.readlines()

for i,line in enumerate(outcar):
    if "  direct lattice vectors                    reciprocal lattice vectors"  in line:
        dd1=outcar[i+1].split()
        dd2=outcar[i+2].split()
        dd3=outcar[i+3].split()
for i in range(3):
    a[0,i]=float(dd1[i])
    a[1,i]=float(dd2[i])
    a[2,i]=float(dd3[i])
    b[0,i]=float(dd1[i + 3])
    b[1,i]=float(dd2[i + 3])
    b[2,i]=float(dd3[i + 3])
print("a=",a)
print("b=",b, "\n")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#mainham書き込み
mainham=open("mainham",mode="w")
wn=[] #1,2 etc.
with open("wannier90_hr.dat") as f:
    for i,line in enumerate(f):
        dam = line.split()
        if i == 0:
            writtentime = dam[2] + dam[3] + dam[4]
        elif i == 1:
            nw = int(dam[0])
        elif i == 2:
            num_sites = int(dam[0])
        elif i < math.ceil(num_sites / 15) + 3: #15→横並び数,+3→開始位置分
            for value in dam:
                wn.append(value)
        else:
            mainham.write(line) 
mainham.close()

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#k_vを計算

#vasp.pyから座標取得
points = point
print("points=",points)
kname = []
for i in range(len(points)):
	kname.append(points[i][0])
	
#波数ベクトルkのリストを作る
k_v = [] 
L_list = [] #各lのノルムを考慮した横軸座標
k_list = [0] #波数座標の値　[0,X,N,G,Z]
L_point = 0 #初期値
k_point = 0
for i in range(len(points)-1):
    lx = points[i+1][1]-points[i][1]
    ly = points[i+1][2]-points[i][2]
    lz = points[i+1][3]-points[i][3]
    K = lx*b[0] + ly*b[1] + lz*b[2]
    L = np.sqrt(K[0]**2 + K[1]**2 + K[2]**2)
    k_point += L
    k_list.append(k_point)
    
    l_x=np.linspace(points[i][1],points[i+1][1],num_points)
    l_y=np.linspace(points[i][2],points[i+1][2],num_points)
    l_z=np.linspace(points[i][3],points[i+1][3],num_points)
    for j in range(num_points):
        k = l_x[j]/N*b[0] + l_y[j]/N*b[1] + l_z[j]/N*b[2]
        k_v.append(k)
        
        L_list.append(L_point)
        L_point += L / num_points
  
#print(L_list)
#print(len(L_list))
#print("k_v=",k_v)
#print("len(k_v)=", len(k_v))


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#mainham_wnを作り、各値から固有値計算

#ホッピングの値を整形してmainham_wnを生成
mainham_wn = open("mainham_wn",mode="w") 
with  open("mainham") as f:
    for i,line in enumerate(f):
        dam = line.split()
        dam[5] = round(float(dam[5]) / float(wn[int(i / (nw**2))]), rprm)#丸め込み
        dam[6] = round(float(dam[6]) / float(wn[int(i / (nw**2))]), rprm)
        mainham_wn.write("   "+dam[0]+"    "+dam[1]+"    "+dam[2]+"    "+dam[3]+"    "+dam[4]+"    "+str(dam[5])+"    "+str(dam[6])+"\n")        
mainham_wn.close()

#ホッピングのリストを作る
T_list = [] #格子ベクトルのリスト
H = np.zeros((nw, nw), dtype = np.complex128)
H_list = []
count = 0 #ハミルトニアンの区切りをつけるためのカウント
hflex=open(flexhop,mode="w")
rflex=open(flexvec,mode="w")
with open ("mainham_wn") as f:
        for i, line in enumerate(f):
            dam = line.split()
            n = np.array([float(dam[0]), float(dam[1]),float(dam[2])])
            T = n[0] * a[0] + n[1] * a[1] + n[2] * a[2]   
            H[int(dam[3])-1, int(dam[4])-1] = complex(float(dam[5]), float(dam[6]))
            hflex.write("   "+dam[0]+"    "+dam[1]+"    "+dam[2]+"    "+dam[3]+"    "+dam[4]+"    "+dam[\
5]+"    "+dam[6]+"\n")
            count += 1
            if count == nw ** 2:
                T_list.append(T)
                H_list.append(H)#(copy.copy(H))のほうが良い？
                H = np.zeros((nw, nw), dtype = np.complex128)
                rflex.write(dam[0]+"    "+dam[1]+"    "+dam[2]+"\n")
                count = 0
mainham_wn.close()

#固有値計算
ene = []
eve = [] #単に固有ベクトルではなくて規格化して成分ごとの割合
for i in range (len(k_v)):
    H_sum = np.zeros((nw,nw), dtype = np.complex128)#ハミルトニアン
    for p in range(len(T_list)):
        FT = np.exp(2.0j * np.pi * np.dot((k_v[i]), T_list[p]))
        H_sum += H_list[p] * FT
    w,v = LA.eigh(H_sum)
    ene.append(w - ef)
    v_norm = np.abs(v)
    V = np.round(v_norm, 3)
    eve.append(V)

#print(k_list)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#以下描画関係
#wan単色
plt.ylim(-3, 4)
plt.xlim(0, max(L_list))
plt.yticks(np.arange(-3, 3, 1))
plt.xticks(k_list, kname)
for i in  k_list:
    plt.axvline(i, color='k', linestyle='--')
plt.axhline(y=0, linestyle="solid", color="k")
plt.plot(L_list, ene, linestyle='solid',c="deeppink",linewidth=1, alpha=1)
plt.savefig(figname+"/wan2.png",format="png",dpi=300)
plt.close()
#all色分け
plt.ylim(-3, 4)
plt.xlim(0, max(L_list))
plt.yticks(np.arange(-3, 4, 1), fontsize=14)
plt.xticks(k_list, kname, fontsize=14)
plt.ylabel('ENERGY [eV]', fontsize=14, labelpad=10)
for i in  k_list:
    plt.axvline(i, color='k', linestyle='--')
plt.axhline(y=0, linestyle="solid", color="k")
for i in range(len(L_list)):
    for j in range(nw):
        color = [(eve[i][0][j]+eve[i][2][j])/2,(eve[i][1][j]+eve[i][3][j])/2,0]
        plt.scatter(L_list[i], ene[i][j],c=[color],marker='o',alpha=0.4)
        red_legend = plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='red', markersize=10, label='dx2-y2')
        green_legend = plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='green', markersize=10, label='dz2')
        plt.legend(handles=[red_legend, green_legend])
plt.savefig(figname+"/all.png",format="png",dpi=300)
