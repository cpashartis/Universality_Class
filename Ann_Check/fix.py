import os

cmd = os.popen("grep 2.0000000000 *")
for data in cmd.readlines():

	line = data.lstrip(':')
	print line

