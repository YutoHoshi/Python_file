#!/bin/bash
#============ SBATCH Directives =======
#SBATCH -p sa
#SBATCH -t 24:00:00
#SBATCH --rsc p=1:t=32:c=64:m=64G
#SBATCH -o output.o
#SBATCH -J 0Pd1
vasp=${HOME}/vasp.6.3.2/bin/vasp_std
wannier90=${HOME}/wannier90-3.1.0/wannier90.x

#module load intel
#module load mkl
#module load intel-mpi

#export OMP_PLACES=cores
#export OMP_PROC_BIND=close
export OMP_STACKSIZE=512m
export MPI_GROUP_MAX=20000
export MPI_COMM_MAX=1000

#pot計算
echo "begin make_pot"
python3 vasp.py -pot
echo "end make_pot"

#scf計算
if [ -d "opt" ]; then
  echo "opt フォルダが存在します。opt/CONTCAR を POSCAR にコピーします。"
  cp opt/CONTCAR POSCAR
else
  echo "opt フォルダが存在しません。"
fi

echo "begin scf_loop"
python3 vasp.py -scf
srun $vasp
mkdir -p scf
cp * scf
echo "end scf_loop"

#band計算
echo "begin band_cal"
cp scf/*CAR .
python3 vasp.py -band
srun $vasp
mkdir -p bands
cp * bands
echo "end band_cal"
date
cp scf/OUTCAR .

#lwannierワニエ軌道の準備計算
echo "begin lwan"
cp scf/KPOINTS KPOINTS
cp scf/*CAR .
cp scf/CHG CHG 
python3 vasp.py -wannier
srun $vasp
cp wannier90.win wannier90.win_VASPmade
echo "end lwan"

#wannier計算
echo "begin wannier90"
srun ${wannier90} wannier90
echo "end wannier90"

# スクリプトの終了
echo "自動化が完了しました。"
