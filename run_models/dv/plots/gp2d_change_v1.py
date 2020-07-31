import numpy as np
import math

#CHANGE GALDEF PARAMETERS
init_galdef = 'm82_02b'
galdefname = 'galdef_54_'+init_galdef
gpfilename = 'gp.'+init_galdef+'.job'

suffix = 'dv'
dir_in = "'../gp_out/'"

num_params = 2
num_galdef = np.zeros(num_params).astype(int)
# num_galdef[0] = 6 #
# num_galdef[1] = 6 #
# num_galdef[2] = 4 #
# num_galdef[3] = 3 #
# num_galdef[4] = 2 #
# num_galdef[5] = 2 #
# num_galdef[6] = 2 #

num_galdef_tot = int(np.prod(num_galdef))
num_start = np.zeros(num_params).astype(int)
# num_start[0] = 30
# num_start[1] = 30
# num_start[2] = 30
# num_start[3] = 10
# num_start[4] = 10
# num_start[5] = 10
# num_start[6] = 10

#STUFF TO GRAPH
num_graph_galdef = np.zeros(num_params).astype(int)
num_graph_galdef[0] = 9 #
num_graph_galdef[1] = 1 #
# num_graph_galdef[2] = 1 #
# num_graph_galdef[3] = 1 #
# num_graph_galdef[4] = 1 #
# num_graph_galdef[5] = 2 #
# num_graph_galdef[6] = 1 #

num_graph_start = np.zeros(num_params).astype(int)
num_graph_start[0] = 10
num_graph_start[1] = 20
# num_graph_start[2] = 10
# num_graph_start[3] = 12
# num_graph_start[4] = 10
# num_graph_start[5] = 10
# num_graph_start[6] = 10

new_val = []

##############################################
##############################################
##############################################

MODEL_NUMBER = 1

if MODEL_NUMBER == 1:
	num_params = 2
	num_galdef = np.zeros(num_params).astype(int)
	num_galdef[0] = 9 #
	num_galdef[1] = 11 #
	
	num_galdef_tot = int(np.prod(num_galdef))
	num_start = np.zeros(num_params).astype(int)
	num_start[0] = 10
	num_start[1] = 10
	
	## arr 1
	arr_1_min = 1.e26
	arr_1_max = 1.e30
	arr_1_diff = np.log10(arr_1_max/arr_1_min)/(num_galdef[0]-1)
	
	arr_1 = []
	for ii in range(0,num_galdef[0]):
		arr_1.append(arr_1_min*10.**(ii*arr_1_diff))
	
	## arr 2
	arr_2_min = 0.
	arr_2_max = 2000.
	arr_2_diff = (arr_2_max-arr_2_min)/(num_galdef[1]-1)
	
	arr_2 = []
	for ii in range(0,num_galdef[1]):
		arr_2.append(arr_2_min+ii*arr_2_diff)
		
new_val.append(arr_1)
new_val.append(arr_2)

# new_val.append([1.0e-1, 4.e-1, 1.e0, 3.e0, 7.e0, 1.e1, 5.e1]) #FF
# new_val.append([50., 100., 200., 400., 600., 800., 1000.]) #WIND
# new_val.append([2.15, 2.20, 2.25]) #z_scale
# new_val.append([-1., -2., -3.]) #B_pow
# new_val.append([0.05, 0.2]) #r_wind
# new_val.append([2000., 2500.]) #v_0
# new_val.append([2.1, 2.2]) #alpha

VAR_PLOT = []
VAR_PLOT.append(r"$D=$ ")
VAR_PLOT.append(r"$V_0=$ ")
# VAR_PLOT.append(r"$ff=$ ")
# VAR_PLOT.append("$p_{B}=$")
# VAR_PLOT.append("$r_{w}=$")
# VAR_PLOT.append("$V_{0}=$")
# VAR_PLOT.append(r"$\alpha=$")

PLOT_TEXT = "'''diffusion coefficient, wind velocity'''"


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
# files.append('gp2d_plot_core_v1.py')
# files.append('gp2d_plot_core_v1.py')

print files
print linechange
print newline

for i in range(0,len(files)):
	source = open(files[i],'r')
	out = open('ext'+files[i],'w')
	source.seek(0)
	for line in source:
		for j in range(0,len(linechange)):
			if line.startswith(linechange[j]):
				line = newline[j]+'\n'
		out.write(line)
	source.close()
	out.close()
	


