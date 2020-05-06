import os

def get_energies(fname):
	energies=[]
	f2=fname.split("_")
	kx=f2[9+2]
	ky=f2[11+2]
	a1=f2[13+2]
	a2=(f2[15+2].split(".txt"))[0]
	fopen=True
	try:
		f=open(fname,'r')
	except IOError:
		fopen=False
	
	energy_found=False
	special_line=""
	if (fopen==True):
		for line in f:
			line=line.strip()
			if (line[0:19]=="Eigenvalue number 0"):
				energy_found=True
				print kx,ky,a1,a2,((((line.split("="))[-1].strip()).split(","))[0]).split("(")[1]
		f.close()
		if (energy_found==False):
			f=open(fname,'r')
			for line in f:
				line=line.strip()
				if (line[0:7]=="Eigs[0]"):
					special_line=line
					energy_found=True
			f.close()
			line=special_line
			print kx,ky,a1,a2,((((line.split("="))[-1].strip()).split(","))[0]).split("(")[1],"NC"

	return energies

a1s=["0.01","0.10","0.20","0.25","0.30","0.40","0.50","0.60","0.70","0.75","0.80","0.90","1.00","1.10","1.20","1.25","1.30","1.40","1.50","1.5708"]
a2s=["0.00","0.10","0.20","0.25","0.30","0.40","0.50","0.60","0.70","0.75","0.80","0.90","1.00","1.10","1.20","1.25","1.30","1.40","1.50","1.5708"]

for a1 in a1s:
	for a2 in a2s:
		fname="/scratch/sciteam/hiteshjc/chED/output_files/Kagome/48_XT8_4/Sz_16.0/Jz_0.00_J2_0.0_Ch_0.0_kx_0_ky_0_a1_"+str(a1)+"_a2_"+str(a2)+".txt"
		energies=get_energies(fname)
