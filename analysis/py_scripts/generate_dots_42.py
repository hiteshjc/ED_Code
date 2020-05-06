import os
import copy
import cmath
import numpy as N
import sys

def make_input_and_job(Jz,Sz,executable,inpathname,outpathname,wf_1_name,wf_2_name):
	input_name= inpathname+"/Dot/dotprod_wf1_"+wf_1_name+"_wf2_"+wf_2_name+".txt"
	output_name= outpathname+"/Dot/dotprod_wf1_"+wf_1_name+"_wf2_"+wf_2_name+".txt"
	ifile=open(input_name,"w")
	ifile.write("N=42_uniformrotate \n\n")
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
	
	if (fopen==True):
		for line in f:
			line=line.strip()
			if (line[0:28]=="Dot product of psi1 and psi2"):dotval=line.split(" ")[-1]
		f.close()
	return dotval

def get_dot_val(wf_1_name,wf_2_name):
	fname="/scratch/sciteam/hiteshjc/chED/output_files/Kagome/42_uniformrotate/Sz_14.0/Dot/dotprod_wf1_"+wf_1_name+"_wf2_"+wf_2_name+".txt"
	r=fname.split("_")
	print r
	a11=r[17-1]
	a21=r[19-1]
	a12=r[33-1]
	a22=r[35-1]
	print a11,a21,a12,a22
	dotval=get_dot_string(fname)
	dotval=dotval.split(",")
	dotval[0]=dotval[0].strip()
	dotval[1]=dotval[1].strip()
	a=float((dotval[0].split("("))[1])
	b=float((dotval[1].split(")"))[0])
	cdot=complex(a,b)
	return cdot

def find_directed_point(point,direction,points):
	b=0
	mind=10000
	sendpoint=["N","N"]
	if (direction=="right"):
		for p in points:
			if (p!=point and point[1]==p[1] and abs(float(point[0])-float(p[0]))<mind and float(point[0])<float(p[0])):
				sendpoint=copy.deepcopy(p)
				mind=abs(float(point[0])-float(p[0]))
				b=1
	if (direction=="up"):
		for p in points:
			if (p!=point and point[0]==p[0] and abs(float(point[1])-float(p[1]))<mind and float(point[1])<float(p[1])):
				sendpoint=copy.deepcopy(p)
				mind=abs(float(point[1])-float(p[1]))
				b=1
	return sendpoint,b

def make_plaqs_and_unique_bonds_given_points(points):
	plaqs=[]
	bonds=[]
	for point in points:
		rpoint,br=find_directed_point(point,"right",points)
		#print "Rpoint",point,rpoint
		upoint,bu=find_directed_point(point,"up",points)
		#print "Upoint",point,upoint
		if (br==1 and bu==1):
			dpoint,bd=find_directed_point(upoint,"right",points)
			if (bd==1):
				#print "Dpoint",point,dpoint
				area=(float(rpoint[0])-float(point[0]))
				area=area*(float(upoint[1])-float(point[1]))
				plaqs.append([[point,rpoint,dpoint,upoint],area])
				bonds.append([point,rpoint])	
				bonds.append([rpoint,dpoint])	
				bonds.append([dpoint,upoint])	
				bonds.append([upoint,point])	
	unibonds=[]
	for bond in bonds:
		b=0
		for unibond in unibonds:
			if (bond==unibond or bond[0]==unibond[1] and bond[1]==unibond[0]):
				b=1
		if (b==0):unibonds.append(bond)
	return plaqs,unibonds

def find_dot(bond,unique_bonds_and_dots):
	for u in unique_bonds_and_dots:
		if (bond[0]==u[0] and bond[1]==u[1]): return u[2]
		if (bond[0]==u[1] and bond[1]==u[0]): return u[2].conjugate()
	return complex(0,0)

kx=str(sys.argv[1])
Jz=str(sys.argv[2])
inpathname="/scratch/sciteam/hiteshjc/chED/input_files/Kagome/42_uniformrotate/Sz_14.0"
outpathname="/scratch/sciteam/hiteshjc/chED/output_files/Kagome/42_uniformrotate/Sz_14.0"
executable="/u/sciteam/hiteshjc/projrepo/ED_Krishna/analysis/pdotprod"

unique_bonds=[]
all_points=[]
pi=3.14159265359
cumuphase=0.0
xs=["0.000","0.698","1.396","2.094","2.793","3.491","4.189","4.887","5.585","6.283185307179586"]
for x in xs:
	unique_bonds.append([[x,"0.000"],["0.000","0.000"]])
unique_bonds_and_dots=[]
print "Getting dot products of bonds..."
for num,bond in enumerate(unique_bonds):
	#print num,bond
	a1=bond[0][0]
	a2=bond[0][1]
	wf_1_name="Jz_"+Jz+"_J2_0.0_Ch_0.0_kx_"+kx+"_ky_0_theta1_"+a1+"_theta2_"+a2+"_pEV1"
	a1=bond[1][0]
	a2=bond[1][1]
	wf_2_name="Jz_"+Jz+"_J2_0.0_Ch_0.0_kx_"+kx+"_ky_0_theta1_"+a1+"_theta2_"+a2+"_pEV1"
	dot_val=get_dot_val(wf_1_name,wf_2_name)
	if (dot_val==complex(0,0)):
		conjdot_val=get_dot_val(wf_2_name,wf_1_name)
		dot_val=conjdot_val.conjugate()
		if (dot_val==complex(0,0)):
			print "Phase unavailable...launching job"
			input_name,job_name=make_input_and_job(Jz,14.0,executable,inpathname,outpathname,wf_1_name,wf_2_name)
			print input_name,job_name
			os.system("chmod 755 "+job_name)
			os.system("qsub "+job_name)
	unique_bonds_and_dots.append([bond[0],bond[1],dot_val])
