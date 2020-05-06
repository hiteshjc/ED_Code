import os
import copy
import cmath
import numpy as N
import sys

def print_cmatrix(matr):
	rows=matr.shape[0]
	cols=matr.shape[1]
	#print rows,cols
	for i in range(rows):
		for j in range(cols):
			if (j!=cols-1): print "( %+8f %+8f j)" %( matr[i,j].real, matr[i,j].imag) ,
			if (j==cols-1): print "( %+8f %+8f j)" %( matr[i,j].real, matr[i,j].imag)


def make_input_and_job(Jz,Sz,executable,inpathname,outpathname,wf_1_name,wf_2_name):
	input_name= inpathname+"/Dot/dotprod_wf1_"+wf_1_name+"_wf2_"+wf_2_name+".txt"
	output_name= outpathname+"/Dot/dotprod_wf1_"+wf_1_name+"_wf2_"+wf_2_name+".txt"
	ifile=open(input_name,"w")
	ifile.write("N=42 \n\n")
	ifile.write("Sz="+str(Sz)+" \n\n")
	ifile.write("Jz="+Jz+" \n\n")

	ifile.write("[Translations for computing momentum eigenvalues] \n")
	ifile.write("T1=[5, 3, 4, 8, 9, 6, 7, 15, 10, 11, 12, 13, 14, 19, 20, 16, 17, 18, 26, 21, 22, 23, 24, 25, 29, 30, 27, 28, 36, 31, 32, 33, 34, 35, 39, 40, 37, 38, 0, 41, 1, 2] \n\n")

	ifile.write("T2=[15, 12, 13, 14, 19, 16, 17, 18, 20, 21, 22, 23, 24, 25, 29, 26, 27, 28, 36, 30, 31, 32, 33, 34, 35, 39, 37, 38, 0, 40, 41, 1, 2, 3, 4, 8, 5, 6, 7, 9, 10, 11] \n\n")

	str1="pEV1_file="+outpathname+"/EV/"+wf_1_name+".txt \n"
	str2="pEV2_file="+outpathname+"/EV/"+wf_2_name+".txt \n"
	ifile.write(str1)
	ifile.write(str2)
	ifile.close()

	job_name= inpathname+"/Dot/job_dotprod_wf1_"+wf_1_name+"_wf2_"+wf_2_name
	jfile=open(job_name,"w")
	jfile.write("#!/bin/bash \n")
	jfile.write("# Batch Queue Script \n")
	jfile.write("#PBS -l walltime=00:10:00 \n")
	jfile.write("#PBS -l nodes=1:ppn=32 \n")
	jfile.write("#PBS -A jpz \n")
	jfile.write("#PBS -q high \n")
	jfile.write("#PBS -j oe \n")
	jfile.write("# Option: -N, names the job \n")
	#jfile.write("cd "+folder+" \n")
	jfile.write("export OMP_NUM_THREADS=32 \n")
	str1="aprun -n 1 -d 32 "+executable+"  "+input_name+"   "+output_name+" \n"
	jfile.write(str1) 
	jfile.close()
	return input_name, job_name 
	
def get_dot_string(fname):
	dotval="(0,0)"
	fopen=True
	try:
		f=open(fname,'r')
	except IOError:
		fopen=False

	#print fopen	
	if (fopen==True):
		for line in f:
			line=line.strip()
			if (line[0:28]=="Dot product of psi1 and psi2"):dotval=line.split(" ")[-1]
		f.close()
	return dotval

def get_dot_val(wf_1_name,wf_2_name):
	fname="/scratch/sciteam/hiteshjc/chED/output_files/Kagome/42b/Sz_14.0/Dot/dotprod_wf1_"+wf_1_name+"_wf2_"+wf_2_name+".txt"
	#/scratch/sciteam/hiteshjc/chED/output_files/Kagome/42b/Sz_14.0/Dot/dotprod_wf1_Jz_0.00_J2_0.0_Ch_0.0_kx_0_ky_0_theta1_0.000_theta2_0.000_pEV1_wf2_Jz_0.00_J2_0.0_Ch_0.0_kx_0_ky_0_theta1_0.698_theta2_0.000_pEV1.txt
	#print fname
	r=fname.split("_")
	#print r
	#a11=r[17]
	#a21=r[19]
	#a12=r[33]
	#a22=r[35]
	dotval=get_dot_string(fname)
	#print dotval
	dotval=dotval.split(",")
	dotval[0]=dotval[0].strip()
	dotval[1]=dotval[1].strip()
	a=float((dotval[0].split("("))[1])
	b=float((dotval[1].split(")"))[0])
	cdot=complex(a,b)
	return cdot

def get_overlap(coeffs_left,coeffs_right,overlap_matrix):
	overlap=complex(0,0)
	for i,cl in enumerate(coeffs_left):
		for j,cr in enumerate(coeffs_right):
			temp=cl
			temp=temp.conjugate()	
			overlap=overlap+(overlap_matrix[i,j]*cr*temp)
	abs_overlap=abs(overlap)
	phase_overlap=cmath.phase(overlap)
	#print abs_overlap,phase_overlap
	return abs_overlap,phase_overlap


def maximize_and_get_overlap(prev_coeffs,overlap_matrix):
	present_coeffs=[]
	new_overlaps=[]
	num_wfs=len(prev_coeffs)
	for i in range(num_wfs):
		overlap=complex(0,0)
		for j in range(num_wfs):
			temp=prev_coeffs[j]
			temp=temp.conjugate()
			overlap=overlap+overlap_matrix[j,i]*(temp)	
		new_overlaps.append(overlap)

	amatrix=N.zeros((num_wfs,num_wfs),dtype=complex)
	
	for i in range(num_wfs):
		for j in range(num_wfs):
			temp=new_overlaps[i]
			temp=temp.conjugate()
			amatrix[i,j]=temp*new_overlaps[j]	

	eigs,vecs=N.linalg.eigh(amatrix)
	#print "Eigenvalues of amatrix = ",eigs
	#print vecs
	overlap=complex(0,0)
	for i in range(num_wfs):
		#overlap=overlap+(new_overlaps[i]*vecs[num_wfs-1][i])
		overlap=overlap+(new_overlaps[i]*vecs[i][num_wfs-1])
	#print "Overlap =",overlap
	abs_overlap=abs(overlap)
	phase_overlap=cmath.phase(overlap)
	for i in range(num_wfs):
		present_coeffs.append(vecs[i][num_wfs-1])
	#print present_coeffs
	#print abs_overlap,phase_overlap
	return present_coeffs,abs_overlap,phase_overlap

xy=sys.argv[1]
Ch=str(sys.argv[2])
if (xy=="x"):
	xs=["0.00","0.698","1.396","2.094","2.793","3.491","4.189","4.887","5.585","6.283185307179586"]
	ylow="0.00"
	trajectory_x1=[]
	for x in xs: trajectory_x1.append([x,ylow])

if (xy=="y"):
	ys=["0.00","0.698","1.396","2.094","2.793","3.491","4.189","4.887","5.585","6.283185307179586"]
	xlow="0.00"
	trajectory_x1=[]
	for y in ys: trajectory_x1.append([xlow,y])

trajectory=[]
for t in trajectory_x1: trajectory.append(t)
wf_trajectory=[]
wf_pairs=[]
for t in trajectory:
	a1=t[0]
	a2=t[1]
	wf_trajectory.append([a1,a2])
num_wfs=len(wf_trajectory)

for i,wf in enumerate(wf_trajectory):
	if (i!=num_wfs-1):
		wf_1_name=wf_trajectory[i]
		wf_2_name=wf_trajectory[i+1]
	else:
		wf_1_name=wf_trajectory[i]
		wf_2_name=wf_trajectory[0]
	if (wf_1_name!=wf_2_name):
		wf_pairs.append([wf_1_name,wf_2_name])

print wf_pairs
#stop

pi=3.14159
cumu_phase=0.0
coeffs=[]
if (xy=="x"):
	c1=0.74
	p1=5.215044
	sqc=N.sqrt(1.0-(c1*c1))
	phase=N.exp(complex(0,1)*p1)
	coeffs.append([c1,sqc*phase])
if (xy=="y"):
	c1=1.000000000000000000000
	p1=0.00
	sqc=N.sqrt(abs(1.0-(c1*c1)))
	phase=N.exp(complex(0,1.000000)*p1)
	coeffs.append([c1,sqc*phase])

evpairs=[[1,1],[1,2],[2,1],[2,2]]
num_wfs=len(coeffs[0])
inpathname="/scratch/sciteam/hiteshjc/chED/input_files/Kagome/42b/Sz_14.0"
outpathname="/scratch/sciteam/hiteshjc/chED/output_files/Kagome/42b/Sz_14.0"
executable="/u/sciteam/hiteshjc/projrepo/ED_Krishna/analysis/pdotprod"

for n,wf_pair in enumerate(wf_pairs):
	a11=wf_pair[0][0]
	a21=wf_pair[0][1]
	a12=wf_pair[1][0]
	a22=wf_pair[1][1]
	overlap_matrix=N.zeros((num_wfs,num_wfs),dtype=complex)
	for evpair in evpairs:
		i=evpair[0]
		j=evpair[1]
		if(i==1): wf_1_name="Jz_0.00_J2_0.0_Ch_"+Ch+"_kx_0_ky_0_theta1_"+a11+"_theta2_"+a21+"_pEV1"
		#/scratch/sciteam/hiteshjc/chED/output_files/Kagome/42b/Sz_14.0/Dot/dotprod_wf1_Jz_0.00_J2_0.0_Ch_0.0_kx_0_ky_0_theta1_0.000_theta2_0.000_pEV1_wf2_Jz_0.00_J2_0.0_Ch_0.0_kx_0_ky_0_theta1_0.698_theta2_0.000_pEV1.txt
		if(i==2): wf_1_name="Jz_0.00_J2_0.0_Ch_"+Ch+"_kx_7_ky_0_theta1_"+a11+"_theta2_"+a21+"_pEV1"
		if(j==1): wf_2_name="Jz_0.00_J2_0.0_Ch_"+Ch+"_kx_0_ky_0_theta1_"+a12+"_theta2_"+a22+"_pEV1"
		if(j==2): wf_2_name="Jz_0.00_J2_0.0_Ch_"+Ch+"_kx_7_ky_0_theta1_"+a12+"_theta2_"+a22+"_pEV1"
		cdot=get_dot_val(wf_1_name,wf_2_name)
		if (cdot==complex(0,0)): 
			cdot=get_dot_val(wf_2_name,wf_1_name)
			if (cdot==complex(0,0)): 
				print 
				print "Phase unavailable for",wf_1_name,wf_2_name
				print "Phase unavailable...launching job"
				Jz="0.00"
				#input_name,job_name=make_input_and_job(Jz,14.0,executable,inpathname,outpathname,wf_1_name,wf_2_name)
				#print input_name,job_name
				#os.system("chmod 755 "+job_name)
				#os.system("qsub "+job_name)
			else: cdot=cdot.conjugate()
		overlap_matrix[j-1,i-1]=cdot
	print "====================================================================================================================="
	print
	print "Overlap matrix"
	print_cmatrix(overlap_matrix)

	if (n!=len(wf_pairs)-1):
		present_coeffs=copy.deepcopy(coeffs[n]) 
		new_coeffs,abs_overlap,phase_overlap=maximize_and_get_overlap(present_coeffs,overlap_matrix)
		coeffs.append(copy.deepcopy(new_coeffs))
		print 
		print "(",a11,a21,") --> ","(",a12,a22,")","Mag overlap =",abs_overlap#,"Phase overlap =",phase_overlap/(2*pi),"*2pi"
		#abs_overlap,phase_overlap=get_overlap(coeffs[n],coeffs[n+1],overlap_matrix)
		print "New (optimal) coeffs = ",new_coeffs
	else:
		abs_overlap,phase_overlap=get_overlap(coeffs[n],coeffs[0],overlap_matrix)
		print 
		print "(",a11,a21,") --> ","(",a12,a22,")","Mag overlap =",abs_overlap#,"Phase overlap =",phase_overlap/(2*pi),"*2pi"
	cumu_phase=cumu_phase+phase_overlap
