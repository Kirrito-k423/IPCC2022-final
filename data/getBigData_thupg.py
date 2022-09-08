from icecream import ic
import re
import time
from tqdm import tqdm
ic.disable()
filename="thupg1.spice"
# filename="test.spice"
write_filename="thupg.mtx"
matchLayer="2"
sorted_dot=dict()
write_list=[]

dot_number=1
fp = open(filename, "r")
num_file = sum([1 for i in open(filename, "r")])
for line in tqdm(fp.readlines(),total=num_file): #读取每一行
	matchObj = re.match( r'^R([0-9]*) (.*) (.*) (.*)$', line, re.M|re.I)
	if matchObj:
		ic(matchObj[2][1])
		if matchObj[2][1]==matchLayer:
			ic(matchObj[2],matchObj[3],matchObj[4])
			raw_point1 = matchObj[2]
			raw_point2 = matchObj[3]
			weight = float(matchObj[4])
			ic(raw_point1, raw_point2, weight)

			if raw_point1 not in sorted_dot:
				sorted_dot[raw_point1] = dot_number
				dot_number += 1
			
			if raw_point2 not in sorted_dot:
				sorted_dot[raw_point2] = dot_number
				dot_number += 1

			point1 = sorted_dot[raw_point1]
			point2 = sorted_dot[raw_point2]
			ic(point1, point2)
			write_list.append([point1 ,point2, weight])
			write_list.append([point2 ,point1, weight])

write_list.sort(key = lambda x: x[1])
ic(len(write_list))

fw = open(write_filename+"_layer_"+matchLayer, "w")
fw.write("{} {} {}\n".format(dot_number-1,dot_number-1,len(write_list)))
for line in tqdm(write_list):
	# write each item on a new line
	fw.write("%d %d %s\n"%(line[0], line[1], line[2]))


	