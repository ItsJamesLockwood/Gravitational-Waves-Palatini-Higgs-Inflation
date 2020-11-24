import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

#wd = os.getcwd()
#file_name = "data/Dataset-Friday_screen.log"
file_name = "data/trouble2_screen.log"

'''
def repair_header(file_path,rep_str=False,out_path=False):
	with open(file_path,'r') as file:
		data = file.readlines()
	if rep_str==False:
		rep_str = "rho/3H^2"
	else:
		pass
	first_index = data[0].index('rho') + len(rep_str)
	data[0] = data[0][:first_index] + data[0][first_index+1:]

	if out_path == False:
		if '.' in file_path:
			point_ind = file_path.index('.')
			out_path = file_path[:point_ind] + '_out' + file_path[point_ind:]
		else:
			out_path = file_path + '_out'
	else:
		pass
	with open(out_path,'w') as outfile:
		for line in data:
			outfile.write(line)

	return out_path, 
out_path = repair_header(file_name)
'''


old_col_names = [chr(i) for i in range(ord('a'),ord('a')+10)]
col_names = ['a',
			'h',
			'omega',
			'pratio',
			'kratio',
			'gratio',
			'mean1',
			'rms1']

data = pd.read_csv(file_name,delim_whitespace=True,skiprows=1,names=col_names,index_col=False)
try:
	data = data[~data.iloc[:,0].str.contains("[a-zA-Z]").fillna(False)]
except AttributeError:
	pass

data = data.apply(pd.to_numeric)

print(data.columns)
print(data.shape)
print(data.dtypes)
print(data.iloc[:5])
print(data.dtypes)

plt.plot(data['h'])
#plt.yscale('log')
#plt.semilogy(data['b'])
plt.show()