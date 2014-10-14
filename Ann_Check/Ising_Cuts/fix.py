import os

new_file = open("Arith_Exact_Miles",'w')
cmd = os.popen("grep 2.0000000000 *")

data_all = cmd.readlines()
data_all.pop()
for data in data_all:

	line = data.replace(':', ' ')
	line = line.split()
	print line
	line = '\t'.join([line[0],line[2]]) + '\n'
	new_file.write(line)

new_file.close()

