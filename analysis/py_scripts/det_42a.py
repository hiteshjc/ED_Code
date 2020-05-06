import os
import copy
import cmath
import numpy as N

def make_input_and_job(Jz,Sz,executable,inpathname,outpathname,wf_1_name,wf_2_name):
	input_name= inpathname+"/Dot/dotprod_wf1_"+wf_1_name+"_wf2_"+wf_2_name+".txt"
	output_name= outpathname+"/Dot/dotprod_wf1_"+wf_1_name+"_wf2_"+wf_2_name+".txt"
	ifile=open(input_name,"w")
	ifile.write("N=42 \n\n")
	ifile.write("Sz="+str(Sz)+" \n\n")
	ifile.write("Jz="+str(Jz)+" \n\n")

	ifile.write("[Translations for computing momentum eigenvalues] \n")
	ifile.write("T1=[6, 3, 4, 5, 10, 11, 7, 8, 9, 19, 12, 13, 14, 15, 16, 17, 18, 24, 25, 20, 21, 22, 23, 34, 26, 27, 28, 29, 30, 31, 32, 33, 38, 39, 35, 36, 37, 0, 40, 41, 1, 2] \n\n")

	ifile.write("T2=[19, 15, 16, 17, 18, 24, 20, 21, 22, 23, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 0, 38, 39, 40, 41, 1, 2, 3, 4, 5, 10, 6, 7, 8, 9, 11, 12, 13, 14] \n\n")

	str1="pEV1_file="+outpathname+"/EV/"+wf_1_name+".txt \n"
	str2="pEV2_file="+outpathname+"/EV/"+wf_2_name+".txt \n"
	ifile.write(str1)
	ifile.write(str2)
	ifile.close()

	job_name= inpathname+"/Dot/job_dotprod_wf1_"+wf_1_name+"_wf2_"+wf_2_name
	jfile=open(job_name,"w")
	jfile.write("#!/bin/bash \n")
	jfile.write("# Batch Queue Script \n")
	jfile.write("#PBS -l walltime=00:20:00 \n")
	jfile.write("#PBS -l nodes=1:ppn=32 \n")
	jfile.write("#PBS -A jpz \n")
	jfile.write("#PBS -q normal \n")
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
	fname="/scratch/sciteam/hiteshjc/chED/output_files/Kagome/42a/Sz_14.0/Dot/dotprod_wf1_"+wf_1_name+"_wf2_"+wf_2_name+".txt"
	#print fname	
	r=fname.split("_")
	#print r
	a11=r[17-1]
	a21=r[19-1]
	a12=r[33-1]
	a22=r[35-1]
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

def find_dot_matrix(bond,unique_bonds_and_dots):
	default=N.zeros((6,6),dtype=complex)
	default=N.mat(default)
	for u in unique_bonds_and_dots:
		if (bond[0]==u[0] and bond[1]==u[1]): return u[2]
		if (bond[0]==u[1] and bond[1]==u[0]): return u[2].conjugate()
	return default


def find_dot(bond,unique_bonds_and_dots):
	for u in unique_bonds_and_dots:
		if (bond[0]==u[0] and bond[1]==u[1]): return u[2]
		if (bond[0]==u[1] and bond[1]==u[0]): return u[2].conjugate()
	return complex(0,0)

inpathname="/scratch/sciteam/hiteshjc/chED/input_files/Kagome/42a/Sz_16.0"
outpathname="/scratch/sciteam/hiteshjc/chED/output_files/Kagome/48_XT8_4uniformgridfinemult/Sz_16.0"
executable="/u/sciteam/hiteshjc/projrepo/ED_Krishna/analysis/pdotprod"

all_points=[]
pi=3.1415926535897932384626433
cumuphase=0.0
xs=["-0.3927","-0.35","-0.30","-0.25","-0.20","-0.15","-0.10","-0.05","-0.0001","0.0001","0.05","0.10","0.15","0.20","0.25","0.30","0.35","0.3927"]
ys=["-0.3927","-0.35","-0.30","-0.25","-0.20","-0.15","-0.10","-0.05","-0.0001","0.0001","0.05","0.10","0.15","0.20","0.25","0.30","0.35","0.3927"]
wfn=[1,2,3,4,5,6]
evpairs=[]
for wf1 in wfn:
	for wf2 in wfn:
		evpairs.append([wf1,wf2])
print evpairs
nums=[]
for evpair in evpairs:
	if (evpair[0] not in nums): 
		nums.append(evpair[0])	
num_wfs=len(nums)
print "num_wfs=",num_wfs
for x in xs:
	for y in ys:
		all_points.append([x,y])
#print all_points
unique_bonds_and_dots=[]
plaqs,unique_bonds=make_plaqs_and_unique_bonds_given_points(all_points)
print 
print 
print "Enumerating plaquettes"
for num,plaq in enumerate(plaqs):
	print num,plaq[0]
print 
print
job=0 
print "Getting dot products of bonds"
for num,bond in enumerate(unique_bonds):
	#print num,bond
	a11=bond[0][0]
	a21=bond[0][1]
	a12=bond[1][0]
	a22=bond[1][1]
	overlap_matrix=N.zeros((num_wfs,num_wfs),dtype=complex)
	overlap_matrix=N.mat(overlap_matrix)
	for evpair in evpairs:
		i=evpair[0]
		j=evpair[1]
		wf_1_name="Jz_0.00_J2_0.0_Ch_0.0_kx_0_ky_0_a1_"+a11+"_a2_"+a21+"_pEV"+str(i)
		wf_2_name="Jz_0.00_J2_0.0_Ch_0.0_kx_0_ky_0_a1_"+a12+"_a2_"+a22+"_pEV"+str(j)
		#print wf_1_name,wf_2_name
		dot_val=get_dot_val(wf_1_name,wf_2_name)
		if (dot_val==complex(0,0)):
			conjdot_val=get_dot_val(wf_2_name,wf_1_name)
			dot_val=conjdot_val.conjugate()
			if (dot_val==complex(0,0)):
				job=job+1
				print "Phase unavailable...launching job"
				if (job>=8000 and job<9000):
					input_name,job_name=make_input_and_job(0.0,16.0,executable,inpathname,outpathname,wf_1_name,wf_2_name)
					print input_name,job_name
					os.system("chmod 755 "+job_name)
					os.system("qsub "+job_name)
		overlap_matrix[j-1,i-1]=dot_val
	print "Bond",bond
	unique_bonds_and_dots.append([bond[0],bond[1],overlap_matrix])

print
print 
print "Unique bond and dot products"
#for u in unique_bonds_and_dots: print "Overlap ",float(u[0][0]),float(u[0][1]),float(u[1][0]),float(u[1][1]),abs(u[2])
areatot=0
cumuphase1=0.0
cumuphase2=0.0
for num,plaq in enumerate(plaqs):
	B=complex(1,0)
	area=plaq[1]
	for i,point in enumerate(plaq[0]):
		#print (i+1)%4
		dot_matrix=find_dot_matrix([point,plaq[0][(i+1)%4]],unique_bonds_and_dots)
		if (i==0):Bmatrix=copy.deepcopy(dot_matrix)
		if (i==1 or i==2 or i==3):Bmatrix=Bmatrix*dot_matrix
		#if (i==2 or i==3):Bmatrix=Bmatrix*N.linalg.inv(dot_matrix)
	B=N.linalg.det(Bmatrix)
	#B=(B/abs(B)) 
	#B=B*dot_val
	Bdgamma=(cmath.phase(B)/(2*pi)) # factor of 16 to go from a1,a2 to theta1, theta2 -> full flux
	if (float(plaq[0][0][0])<=float(plaq[0][0][1])): cumuphase1=cumuphase1+Bdgamma
	if (float(plaq[0][0][0])>float(plaq[0][0][1])): cumuphase2=cumuphase2+Bdgamma
	areatot=areatot+(area)
	cumuphase=cumuphase+Bdgamma
	print "Total Theta for plaquette ",'%03i' %num,"[",'%+05f' % (0.5*(float(plaq[0][0][0])+float(plaq[0][1][0]))),",",'%+05f' %(0.5*(float(plaq[0][0][1])+float(plaq[0][3][1]))),"]"," = ",'%+08f' %(Bdgamma)
print 
print "===================================================================================================="
print "Total Bdgamma (lower)/(2*pi)(a1,a2)    for all plaquettes = ",'%+10f' %cumuphase1
print "Total Bdgamma (upper)/(2*pi)(a1,a2)    for all plaquettes = ",'%+10f' %cumuphase2
print "Total Bdgamma (both )/(2*pi)(a1,a2)    for all plaquettes = ",'%+10f' %(cumuphase)
print "Total Area(a1,a2)                      for all plaquettes = ",'%+10f' %areatot
print "====================================================================================================="
