import numpy as np
import math
import os

#CHANGE GALDEF PARAMETERS
init_galdef = 'm82_05aw'
galdefname = 'galdef_54_'+init_galdef
gpfilename_f = 'gp.'+init_galdef
gpfilename = 'gp.'+init_galdef+'.job'

suffix = 'ff'

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
	ff_min = 1.e3
	ff_max = 1.e-1
	pow_diff = math.log10(ff_max/ff_min)/(num_galdef[1]-1)
	
	arr_2 = []
	for ii in range(0,num_galdef[1]):
		arr_2.append(ff_min*10.**(pow_diff*ii))


new_val.append(arr_1)
new_val.append(arr_2)

print '# '+str(arr_1)
print '# '+str(arr_2)

linechange = [] 
newline = [] 
midline = []
oidline = []
endline = []

#LINES PARAM1
cline = []
nline = []
mline = []
oline = []
eline = []

cline.append('AAAAAAAAAAAA')
nline.append('AAAAAAAAAAAA')
mline.append('')
oline.append('') 
eline.append('')

linechange.append(cline)
newline.append(nline)
midline.append(mline)
oidline.append(oline)
endline.append(eline)

#LINES PARAM2
cline = []
nline = []
mline = []
oline = []
eline = []

cline.append('HII_clumping_factor ')
nline.append('HII_clumping_factor  = ')
mline.append('')
oline.append('') 
eline.append('')

linechange.append(cline)
newline.append(nline)
midline.append(mline)
oidline.append(oline)
endline.append(eline)

###################

source_galdef = open(galdefname, 'r')
source_gp = open(gpfilename,'r')

if not os.path.isdir('./'+init_galdef+suffix):
	os.mkdir('./'+init_galdef+suffix)

index=np.zeros([num_galdef_tot,num_params])
index_n=index
N=np.zeros(num_params).astype(int)
for gt in range(0,num_galdef_tot):
	
	#Getting suffix for files
	N_index=N+num_start
	N_suffix = ''
	for s in range(0,num_params):
		N_suffix += str(N_index[s])
	
	#New galdef file
	new_galdefname = galdefname+suffix+N_suffix
	new_galdef = open('./'+init_galdef+suffix+'/'+new_galdefname, 'w')
	
	#Writing new galdef file
	source_galdef.seek(0)
	for line in source_galdef:
		for ii in range(0,len(linechange)):
			for ij in range(0,len(linechange[ii])):
				if line.startswith(linechange[ii][ij]):
					
					if linechange[ii][ij] == 'B_field_parameters ':
						line = newline[ii][ij]+str(new_val[ii][N[ii]])+midline[ii][ij]+oidline[ii][ij]+endline[ii][ij]+'\n'
					elif linechange[ii][ij] == 'gas_model_parameters ':
						line = newline[ii][ij]+str(new_val[ii][N[ii]])+midline[ii][ij]+oidline[ii][ij]+endline[ii][ij]+'\n'
					else:
						line = newline[ii][ij]+str(new_val[ii][N[ii]])+endline[ii][ij]+'\n'
		new_galdef.write(line)
	new_galdef.close()
	
	#New job file
	new_gpname = gpfilename_f+suffix+N_suffix+'.job'
	new_gp = open('./'+init_galdef+suffix+'/'+new_gpname, 'w')
	
	#Writing new job file
	source_gp.seek(0)
	for line in source_gp:
		if line.startswith('GALPROPNAME'):
			line = 'GALPROPNAME='+init_galdef+suffix+N_suffix+'\n'
		if line.startswith('#PBS -N'):
			line = '#PBS -N '+new_gpname+'\n'
		new_gp.write(line)	
	new_gp.close()
	print 'cp '+'./'+new_galdefname+' /users/PCON0003/cond0064/galprop/galdef/'+'\n'
	print 'qsub '+new_gpname+'\n'

	#Update N
	N[num_params-1] += 1
	for i in (num_params-1 -np.array(range(0,num_params))):
		if N[i] >= num_galdef[i]:
			N[i] = 0
			if i >= 1:
				N[i-1] += 1
	
source_galdef.close()

