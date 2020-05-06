import os
os.system("grep scratch tt >t3")
f=open("t3","r")
stri="["
for line in f:
	line=line.strip()
	line=line.split("_")
	print line[15]
	if (line[15][0]!="-"):stri=stri+"\""+line[15]+"\","
print stri+"]"
