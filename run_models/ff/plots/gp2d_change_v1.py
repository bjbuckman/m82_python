import numpy as np
import math

#CHANGE GALDEF PARAMETERS
init_galdef = 'm82_05a'
galdefname = 'galdef_54_'+init_galdef
gpfilename = 'gp.'+init_galdef+'.job'

suffix = 'ff'
# free-free
dir_in = "'../gp_out/'"

num_params = 2
num_galdef = np.zeros(num_params).astype(int)

num_galdef_tot = int(np.prod(num_galdef))
num_start = np.zeros(num_params).astype(int)

#STUFF TO GRAPH
num_graph_galdef = np.zeros(num_params).astype(int)
num_graph_galdef[0] = 1 #
num_graph_galdef[1] = 5 #

num_graph_start = np.zeros(num_params).astype(int)
num_graph_start[0] = 10
num_graph_start[1] = 16

new_val = []

##############################################
##############################################
##############################################

MODEL_NUMBER = 1

if MODEL_NUMBER == 1:
	num_params = 2
	num_galdef = np.zeros(num_params).astype(int)
	num_galdef[0] = 1 #
	num_galdef[1] = 20 #
	
	num_galdef_tot = int(np.prod(num_galdef))
	num_start = np.zeros(num_params).astype(int)
	num_start[0] = 10
	num_start[1] = 10
	
	## argument 1
	# z_min = 1.e27
	# z_max = 1.e29
	# pow_diff = math.log10(z_max/z_min)/(num_galdef[0]-1)
	
	arr_1 = [1.]
	
	## argument 2
	ff_min = 1.e0
	ff_max = 1.e-4
	pow_diff = math.log10(ff_max/ff_min)/(num_galdef[1]-1)
	
	arr_2 = []
	for ii in range(0,num_galdef[1]):
		arr_2.append(ff_min*10.**(pow_diff*ii))	


		
new_val.append(arr_1)
new_val.append(arr_2)

VAR_PLOT = []
VAR_PLOT.append(r"")
VAR_PLOT.append(r"$C_{ff}=$ ")

PLOT_TEXT = "'''dfree-free'''"


NORM_PLOT = 0 # #0 -> USE CHI_NORM, #1 -> Straight GP
######################################

E2FLAG = 1 #LEAVE! NOT WRITTEN YET

new_val_arr = np.zeros([num_params, num_graph_galdef.max()])
for n in range(0,num_params):
	for m in range(0,num_graph_galdef[n]):
		new_val_arr[n][m] = new_val[n][m+num_graph_start[n]-num_start[n]]

GALDEF = init_galdef
CHANGED = suffix

N=np.zeros(num_params).astype(int)
galdef_name = []
for n in range(0,num_galdef_tot):
	#Getting suffix for files
	N_index=N+num_start
	N_suffix = ''
	for s in range(0,num_params):
		N_suffix += str(N_index[s])
		
	galdef_name.append(GALDEF+CHANGED+N_suffix)
	
	#Update N
	N[num_params-1] += 1
	for i in (num_params-1 -np.array(range(0,num_params))):
		if N[i] >= num_galdef[i]:
			N[i] = 0
			if i >= 1:
				N[i-1] += 1
				
linechange = []
newline = []

linechange.append('num_params ')
newline.append('num_params = '+str(num_params))

linechange.append('num_galdef ')
newline.append('num_galdef = np.'+str(np.array_repr(num_galdef)))

linechange.append('BEGIN ')
newline.append('BEGIN = np.'+str(np.array_repr(num_start)))

linechange.append('GALDEF ')
newline.append("GALDEF = '"+str(GALDEF)+"'")

linechange.append('CHANGED ')
newline.append("CHANGED = '"+str(CHANGED)+"'")

linechange.append('new_val ')
newline.append('new_val = np.'+str(np.array_repr(new_val_arr)))

linechange.append('VAR_PLOT ')
newline.append('VAR_PLOT = '+str(VAR_PLOT))

linechange.append('num_graph_start ')
newline.append('num_graph_start = np.'+str(np.array_repr(num_graph_start)))

linechange.append('num_graph_galdef ')
newline.append('num_graph_galdef = np.'+str(np.array_repr(num_graph_galdef)))

linechange.append('E2FLAG ')
newline.append('E2FLAG = '+str(E2FLAG))

linechange.append('NORM_PLOT ')
newline.append('NORM_PLOT = '+str(NORM_PLOT))

linechange.append('plot_text ')
newline.append('plot_text = '+str(PLOT_TEXT))

linechange.append('dir_in ')
newline.append('dir_in = '+dir_in)

files = []
files.append('gp2d_plot_v2.py')

print( files)
print( linechange)
print( newline)

for i in range(0,len(files)):
	source = open(files[i],'r')
	out = open('ext_'+files[i],'w')
	source.seek(0)
	for line in source:
		for j in range(0,len(linechange)):
			if line.startswith(linechange[j]):
				line = newline[j]+'\n'
		out.write(line)
	source.close()
	out.close()
	


