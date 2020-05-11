#!/bin/csh

set N=48_XT8_4
set Sz=16.0

set folder = /scratch/sciteam/hiteshjc/chED
set inpath = input_files/Kagome/${N}/Sz_${Sz}
set outpath = output_files/Kagome/${N}/Sz_${Sz}
set ex_file = /u/sciteam/hiteshjc/projrepo/ED_Krishna/analysis/pspincorr

mkdir ${folder}/${inpath}/spin_spin
mkdir ${folder}/${outpath}/spin_spin

set Ch = 0.0
set J2 = 0.0
set J1 = 1.0

set kx = 0
set ky = 0
set a1=0.00000000000
set a2=0.00000000000


foreach Jz  (0.00)
foreach ev1 (1)
foreach ev2 (1)
cat>${folder}/${inpath}/spin_spin/Spin_spin_Jz_${Jz}_J2_${J2}_Ch_${Ch}_kx_${kx}_ky_${ky}_a1_${a1}_a2_${a2}_ev1_${ev1}_ev2_${ev2}.txt  << EOFm

N=48

Sz=${Sz}

T1=[2, 3, 4, 5, 6, 7, 0, 1, 9, 10, 11, 8, 14, 15, 16, 17, 18, 19, 12, 13, 21, 22, 23, 20, 26, 27, 28, 29, 30, 31, 24, 25, 33, 34, 35, 32, 38, 39, 40, 41, 42, 43, 36, 37, 45, 46, 47, 44]

T2=[12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]


[Order of triangles is important]
Triangles = [ [0,1,8], [2,3,9], [4,5,10], [6,7,11], [9,14,13], [10,16,15], [11,18,17], [8,12,19], [12,13,20], [14,15,21], [16,17,22], [18,19,23], [21,26,25],[22,28,27], [23,30,29], [20,24,31], [24,25,32], [26,27,33], [28,29,34], [30,31,35], [33,38,37], [34,40,39], [35,42,41], [32,36,43], [36,37,44], [38,39,45], [40,41,46], [42,43,47], [45,2,1], [46,4,3], [47,6,5], [44,0,7] ]

pEV1_file=${folder}/${outpath}/EV_back/Jz_${Jz}_J2_${J2}_Ch_${Ch}_kx_${kx}_ky_${ky}_a1_${a1}_a2_${a2}_pEV1.txt
pEV2_file=${folder}/${outpath}/EV_back/Jz_${Jz}_J2_${J2}_Ch_${Ch}_kx_${kx}_ky_${ky}_a1_${a1}_a2_${a2}_pEV2.txt

ev1= ${ev1}
ev2= ${ev2}

EOFm

cat >${folder}/${inpath}/spin_spin/job_spin_spin_Jz_${Jz}_J2_${J2}_Ch_${Ch}_kx_${kx}_ky_${ky}_a1_${a1}_a2_${a2}_ev1_${ev1}_ev2_${ev2} << EOFm

#!/bin/bash
# Batch Queue Script
#PBS -l walltime=00:30:00
#PBS -l nodes=1:ppn=32
#PBS -A jpz
#PBS -q debug
#PBS -j oe
# Option: -N, names the job
#PBS -N ${N}_spin_spin_Jz_${Jz}_J2_${J2}_Ch_${Ch}_kx_${kx}_ky_${ky}

cd ${folder}
export OMP_NUM_THREADS=32
aprun -n 1 -d 32 ${ex_file} ${folder}/${inpath}/spin_spin/Spin_spin_Jz_${Jz}_J2_${J2}_Ch_${Ch}_kx_${kx}_ky_${ky}_a1_${a1}_a2_${a2}_ev1_${ev1}_ev2_${ev2}.txt ${folder}/${outpath}/spin_spin/Spin_spin_Jz_${Jz}_J2_${J2}_Ch_${Ch}_kx_${kx}_ky_${ky}_a1_${a1}_a2_${a2}_ev1_${ev1}_ev2_${ev2}.txt 

EOFm

chmod 755 ${folder}/${inpath}/spin_spin/job_spin_spin_Jz_${Jz}_J2_${J2}_Ch_${Ch}_kx_${kx}_ky_${ky}_a1_${a1}_a2_${a2}_ev1_${ev1}_ev2_${ev2}
qsub ${folder}/${inpath}/spin_spin/job_spin_spin_Jz_${Jz}_J2_${J2}_Ch_${Ch}_kx_${kx}_ky_${ky}_a1_${a1}_a2_${a2}_ev1_${ev1}_ev2_${ev2}

end
end
end #Jz
