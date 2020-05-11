#!/bin/csh

set N=30
set Sz=0.0

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

#foreach Jz  (-0.10 -0.20 -0.30 -0.40 -0.49)
foreach Jz  (-0.60 -0.58 -0.55 -0.52 -0.50 -0.495 -0.49 -0.45 -0.40 -0.35 -0.30 -0.20 -0.10 0.00 0.10 0.20 0.30 0.40 0.50 0.60 0.70 0.80 0.90 1.00)
#foreach Jz  (-0.495 -0.45 -0.35 -0.52 -0.55 -0.58)
foreach ev1 (1)
foreach ev2 (1)
cat>${folder}/${inpath}/spin_spin/Spin_spin_Jz_${Jz}_J2_${J2}_Ch_${Ch}_kx_${kx}_ky_${ky}_ev1_${ev1}_ev2_${ev2}.txt  << EOFm

N=30

Sz=${Sz}

[Translations for computing momentum eigenvalues]
T1=[26, 3, 28, 29, 5, 0, 8, 9, 10, 1, 2, 12, 13, 4, 16, 17, 18, 6, 7, 20, 11, 23, 24, 25, 14, 15, 27, 19, 21, 22]

T2=[4, 7, 8, 9, 12, 13, 15, 16, 17, 18, 6, 19, 20, 11, 22, 23, 24, 25, 14, 26, 27, 2, 3, 28, 29, 21, 5, 0, 10, 1]

[Order of triangles is important]
Triangles = [[0,1,10], [26,3,2], [1,2,4], [3,28,5], [4,7,6], [5,9,8], [18,6,11], [7,8,12], [9,10,13], [11,14,25], [12,16,15], [13,18,17], [14,15,19], [16,17,20], [19,22,21], [20,24,23], [22,23,26], [24,25,27], [27,29,28], [29,21,0] ]

pEV1_file=${folder}/${outpath}/EV/Jz_${Jz}_J2_0.0_Ch_0.00_ChCh_0.00_kx_${kx}_ky_${ky}_pEV${ev1}.txt
pEV2_file=${folder}/${outpath}/EV/Jz_${Jz}_J2_0.0_Ch_0.00_ChCh_0.00_kx_${kx}_ky_${ky}_pEV${ev2}.txt

ev1= ${ev1}
ev2= ${ev2}

EOFm

cat >${folder}/${inpath}/spin_spin/job_spin_spin_Jz_${Jz}_J2_${J2}_Ch_${Ch}_kx_${kx}_ky_${ky}_ev1_${ev1}_ev2_${ev2} << EOFm

#!/bin/bash
# Batch Queue Script
#PBS -l walltime=00:30:00
#PBS -l nodes=1:ppn=32
#PBS -A jpz
#PBS -q low
#PBS -j oe
# Option: -N, names the job
#PBS -N ${N}_spin_spin_Jz_${Jz}_J2_${J2}_Ch_${Ch}_kx_${kx}_ky_${ky}

cd ${folder}
export OMP_NUM_THREADS=32
aprun -n 1 -d 32 ${ex_file} ${folder}/${inpath}/spin_spin/Spin_spin_Jz_${Jz}_J2_${J2}_Ch_${Ch}_kx_${kx}_ky_${ky}_ev1_${ev1}_ev2_${ev2}.txt ${folder}/${outpath}/spin_spin/Spin_spin_Jz_${Jz}_J2_${J2}_Ch_${Ch}_kx_${kx}_ky_${ky}_ev1_${ev1}_ev2_${ev2}.txt 

EOFm

chmod 755 ${folder}/${inpath}/spin_spin/job_spin_spin_Jz_${Jz}_J2_${J2}_Ch_${Ch}_kx_${kx}_ky_${ky}_ev1_${ev1}_ev2_${ev2}
qsub ${folder}/${inpath}/spin_spin/job_spin_spin_Jz_${Jz}_J2_${J2}_Ch_${Ch}_kx_${kx}_ky_${ky}_ev1_${ev1}_ev2_${ev2}

end
end
end #Jz
