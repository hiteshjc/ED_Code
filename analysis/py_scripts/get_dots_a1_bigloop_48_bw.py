import os
import copy
import cmath

def get_dot(fname):
	dotval="(0,0)"
	fopen=True
	try:
		f=open(fname,'r')
	except IOError:
		fopen=False
	
	if (fopen==True):
		for line in f:
			line=line.strip()
			if (line[0:28]=="Dot product of psi1 and psi2"):dotval=line.split(" ")[-1]
		f.close()
	return dotval


#xs=["0.001","0.25","0.50","0.75","1.00","1.25","1.5708"]
xs=["0.001","0.25","0.50","0.75","1.00","1.25","1.5708","1.75","2.00","2.25","2.50","2.75","3.00","3.25","3.50","3.75","4.00","4.25","4.50","4.75","5.00","5.25","5.50","5.75","6.00"]
ylow="0.00"
trajectory_x1=[]
for x in xs: trajectory_x1.append([x,ylow])
trajectory=[]
for t in trajectory_x1: trajectory.append(t)
wf_trajectory=[]
wf_pairs=[]
for t in trajectory:
	a1=t[0]
	a2=t[1]
	wf_trajectory.append("Jz_0.00_J2_0.0_Ch_0.0_kx_0_ky_0_a1_"+a1+"_a2_"+a2+"_pEV1")
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
		#print wf_1_name,wf_2_name

pi=3.14159
totaldot=complex(1,0.0)
for n,wf_pair in enumerate(wf_pairs):
	wf_1_name=wf_pair[0]
	wf_2_name=wf_pair[1]
	fname="/scratch/sciteam/hiteshjc/chED/output_files/Kagome/48_XT8_4/Sz_16.0/Dot/dotprod_wf1_"+wf_1_name+"_wf2_"+wf_2_name+".txt"
	r=fname.split("_")
	#print r
	a11=r[17]
	a21=r[19]
	a12=r[33]
	a22=r[35]
	dotval=get_dot(fname)
	dotval=dotval.split(",")
	dotval[0]=dotval[0].strip()
	dotval[1]=dotval[1].strip()
	#print dotval[0]
	#print dotval[1]
	a=float((dotval[0].split("("))[1])
	b=float((dotval[1].split(")"))[0])
	cdot=complex(a,b)
	print 
	print "(",a11,a21,") --> ","(",a12,a22,")",dotval,cdot
	#print wf_1_name,wf_2_name,dotval,cdot
	if (cdot!=complex(0,0)): totaldot=totaldot*cdot
	phase1=cmath.phase(cdot)
	if (phase1<0.0): phase1=(2*pi)+phase1
	phase2=cmath.phase(totaldot)
	if (phase2<0.0): phase2=(2*pi)+phase2
	print "(",a11,a21,") --> ","(",a12,a22,")",abs(cdot),phase1/(2*pi)
	print
	print "Cumulative Chern number/phase (div by 2pi included))", phase2/(2*pi)
#	inpathname="/home/hiteshjc/project-cse/big_in/Kagome/48/Sz_16.0"
#	outpathname="/home/hiteshjc/project-cse/big_out/Kagome/48/Sz_16.0"
#	input_name,job_name=make_input_and_job(1.0,16.0,executable,inpathname,outpathname,wf_1_name,wf_2_name)
#	os.system("chmod 755 "+job_name)
#	os.system("qsub "+job_name)
