import sys
x =  sys.argv[1]
file = open(x,"r")
lines = file.readlines()
for line in lines:
	if ">" in line:
		strings = line.split(" ")[0]
		print(strings)
	else:
		print(line.rstrip())

