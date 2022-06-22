#! /usr/bin/python 
#mario xerri 6/11/22

from sys import exit
import os
import pandas as pd
import numpy as np
from math import cos
from math import sin 
from math import sqrt
import diffpy.structure as dp
from niggli import niggli_reduce
import subprocess


df = pd.read_csv("spreadsheet_formatted_13_33_04.csv",low_memory =False).fillna("x")
df.columns = ["id_code","id_code2","a","b","c","alpha","beta","gamma", "spg_type","spg1","spg2","spg3","spg4","?"]
df_trimmed = df.iloc[3:].reset_index()

#fixing formatting and getting space group in accessible format
df_trimmed["a"] = df_trimmed["a"].astype(float)
df_trimmed["b"] = df_trimmed["b"].astype(float)
df_trimmed["c"] = df_trimmed["c"].astype(float)
df_trimmed["alpha"] = df_trimmed["alpha"].astype(float)
df_trimmed["beta"] = df_trimmed["beta"].astype(float)
df_trimmed["gamma"] = df_trimmed["gamma"].astype(float)
df_trimmed["spg_concat"] = (df_trimmed["spg1"] + " " + df_trimmed["spg2"] + " " + df_trimmed["spg3"] + " "  + df_trimmed["spg4"]).str.rstrip(" x x")
df_trimmed["z"] = df_trimmed["spg_concat"].str.split(" ").str[-1]
df_trimmed["spg_concat"] = df_trimmed.apply(lambda row:row["spg_concat"].rstrip(row["z"]),axis = 1)

df_cell_param =  df_trimmed.iloc[:,[3,4,5,6,7,8,9,15,16]]

#global dataframes
df_reduced = pd.DataFrame()
df_reduced_vectors = pd.DataFrame()
df_g6 = pd.DataFrame()
df_tau = pd.DataFrame()
df_output = pd.DataFrame()
df_g6redprim = pd.DataFrame(columns = ["r","s","t","u","v","w","spg_type","spg_concat","z"])

def toPrimitive(cell):
	""" This converts the non-primitive cells and primitive cells to niggli reduced primitive G6 representation 
	    The code called was devleoped by Herbert Bernstein and Lawrence Andrews and can be found in https://github.com/yayahjb/ncdist

	    Parameters: pandas DataFrame of a b c  alpha beta gamma spg_type spg z 
	"""

	cell_string = " ".join(str(item) for item in cell)
	bash_script = subprocess.Popen(["/home/mxerri/SULI/ncdist21J/ncdist/dc7unsrtdist_mat --info"],shell = True,executable = "/bin/bash",stdin = subprocess.PIPE,stdout = subprocess.PIPE)
	bash_script.stdin.write(cell_string)
	out = bash_script.communicate()[0]
	bash_script.stdin.close()
	primitive =(((out.split("\n")[1]).split(":")[1]).lstrip()).split()
	return([float(x) for x in primitive])

#Calls the toPrimitive Function for every row and determines the G6 reduced preimitive representations
df_g6redprim =  df_cell_param.apply(lambda row: pd.Series(toPrimitive([row["spg_type"],row["a"],row["b"],row["c"],row["alpha"],row["beta"],row["gamma"]])\
+["P*",row["spg_concat"],row["z"]],index = df_g6redprim.columns) if row["spg_type"]!= "P" else pd.Series(toPrimitive([row["spg_type"],row["a"],row["b"],row["c"],row["alpha"],\
row["beta"],row["gamma"]])+["P",row["spg_concat"],row["z"]],index = df_g6redprim.columns),axis = 1)


def square(i):
	"""helper function  for finding magnitude of vector
	"""
	return np.dot(i,i)

#Not being used anymore becasue im implementing code from https://github.com/yayahjb/ncdist 
def niggli_reduction():
	""" 
	Thiis function niggli takes the 6 cell parameters and outputs the niggli reduced g6 vectors 

	I  utilize code by atztogo in his niggli repository on github to convert abc basis vectors to the niggli reduced form:
	https://github.com/atztogo/niggli/tree/develop/python/niggli

	I utilize code by diffpy in the diffpy.structure repository to convert cell parameters to abc basis vectors:
	https://github.com/diffpy/diffpy.structure
	"""
	basis_matrix_array  = (df_cell_param.apply(lambda row:np.transpose((dp.Lattice(row[0],row[1],row[2],row[3],row[4],row[5]).base)) ,axis = 1)).to_numpy(dtype =np.ndarray)

	reduced_basis_matrix_array = map(niggli_reduce,basis_matrix_array)
	red_a = []
	red_b = []
	red_c = []

	#this sorts ||a|| , ||b||,||c||:  a being smallest and c being largest  
	for matrix in reduced_basis_matrix_array:
		unsrted_abc_array = [matrix[:,0],matrix[:,1],matrix[:,2]]
		srted_abc_array = sorted(unsrted_abc_array,key = square)

		red_a.append(srted_abc_array[0])
		red_b.append(srted_abc_array[1])
		red_c.append(srted_abc_array[2])

	red_a_s = pd.Series(red_a)
	red_b_s = pd.Series(red_b)
	red_c_s = pd.Series(red_c)

	df_reduced_vectors = pd.concat([red_a_s,red_b_s,red_c_s], axis = 1)
	df_reduced_vectors.columns = ["red_a","red_b","red_c"]

	#calculates g6 parameters 
	df_reduced_vectors["r"] = df_reduced_vectors.apply(lambda row: np.dot(row["red_a"],row["red_a"]), axis = 1)
	df_reduced_vectors["s"] = df_reduced_vectors.apply(lambda row: np.dot(row["red_b"],row["red_b"]), axis = 1)
	df_reduced_vectors["t"] = df_reduced_vectors.apply(lambda row: np.dot(row["red_c"],row["red_c"]), axis = 1)
	df_reduced_vectors["u"] = df_reduced_vectors.apply(lambda row: 2*np.dot(row["red_b"],row["red_c"]), axis = 1)
	df_reduced_vectors["v"] = df_reduced_vectors.apply(lambda row: 2*np.dot(row["red_a"],row["red_c"]), axis = 1)
	df_reduced_vectors["w"] = df_reduced_vectors.apply(lambda row: 2*np.dot(row["red_a"],row["red_b"]), axis = 1)

	global df_g6
	df_g6 = df_reduced_vectors.iloc[:,3:]

	g6_to_dc7()

df_dc7 = pd.DataFrame()

def g6_to_dc7():
	"""converts g6 parameters to dc7unsrt parameters for +++ case
	Parameters: dataframe
	Predoncition: dataframe is a type pd.DataFrame object
	"""
	global df_dc7
	global df_g6redprim
	df_dc7["dc7_1"] =df_g6redprim["r"]**(.5)
	df_dc7["dc7_2"] =df_g6redprim["s"] **(.5)
	df_dc7["dc7_3"] =df_g6redprim["t"]**(.5)
	df_dc7["dc7_4"] =(df_g6redprim.apply(lambda row: row["s"] + row["t"] - abs(row["u"]),axis = 1))**(.5) 
	df_dc7["dc7_5"] =(df_g6redprim.apply(lambda row: row["r"] + row["t"] - abs(row["v"]),axis = 1))**(.5)
	df_dc7["dc7_6"] =df_g6redprim.apply(lambda row: row["r"] + row["s"] - abs(row["w"]),axis = 1) **(.5)

	df_dc7["dc7_7"] = df_g6redprim.apply(lambda row: min(row["r"] + row["s"]+row["t"] + row["u"] + row["v"] + row["w"], \
	row["r"] + row["s"]+row["t"] + row["u"] - row["v"] - row["w"],row["r"] + row["s"]+row["t"] - row["u"] + row["v"] - row["w"], \
	row["r"] + row["s"]+row["t"] - row["u"] - row["v"] + row["w"]),axis = 1)**(.5)

	recover_uvw(df_dc7)


df_g6_reco = pd.DataFrame()
def recover_uvw(df_dc7):
	""" Recovers g6 uvw parameters from dc7unsrt
	Preconditions dataframe that contains dc7unsrt data
	"""
	global df_g6_reco
	global df_tau
	df_g6_reco["u_reco"] = (df_dc7["dc7_2"]**2+df_dc7["dc7_3"]**2-df_dc7["dc7_4"]**2).round(6) 
	df_g6_reco["v_reco"] = (df_dc7["dc7_1"]**2+df_dc7["dc7_3"]**2-df_dc7["dc7_5"]**2).round(6)
	df_g6_reco["w_reco"] = (df_dc7["dc7_1"]**2+df_dc7["dc7_2"]**2-df_dc7["dc7_6"]**2).round(6)
	df_g6_reco["u_reco"] = df_g6_reco["u_reco"].abs()
	df_g6_reco["v_reco"] = df_g6_reco["v_reco"].abs()
	df_g6_reco["w_reco"] = df_g6_reco["w_reco"].abs()

	df_g6_reco["r"] = df_g6redprim["r"]
	df_g6_reco["s"] = df_g6redprim["s"]
	df_g6_reco["t"] = df_g6redprim["t"] 
	df_g6_reco["u_og"] = df_g6redprim["u"]
	df_g6_reco["v_og"] = df_g6redprim["v"]
	df_g6_reco["w_og"] = df_g6redprim["w"] 

	tau = df_g6_reco["r"] + df_g6_reco["s"] + df_g6_reco["t"] - (df_g6_reco["u_reco"].abs()+df_g6_reco["v_reco"].abs()+df_g6_reco["w_reco"].abs())
	tau= tau.round(6)

	df_tau = tau
	df_tau.columns = ["tau"]
	print(pd.concat([df_dc7["dc7_7"]**2,df_tau],axis = 1))

	#checks to see if tau equals  7th element of dcunsrt and determines sign of u,v,w  for g6 
	df_g6_reco.loc[tau == (df_dc7["dc7_7"]**2).round(6), "case"] = "---"
	df_g6_reco.loc[tau != (df_dc7["dc7_7"]**2).round(6), "case"] = "+++"
	# adjusts the sign of  u,v,w after determining if dc7unsrt_7 == tau 
	df_g6_reco.loc[df_g6_reco["case"]== "+++","u"] = df_g6_reco["u_reco"].round(6)
	df_g6_reco.loc[df_g6_reco["case"]== "---","u"] =-1 *  df_g6_reco["u_reco"] .round(6)
	df_g6_reco.loc[df_g6_reco["case"]== "+++","v"] = df_g6_reco["v_reco"].round(6)
	df_g6_reco.loc[df_g6_reco["case"]== "---","v"] =-1 *  df_g6_reco["v_reco"] .round(6)
	df_g6_reco.loc[df_g6_reco["case"]== "+++","w"] = df_g6_reco["w_reco"].round(6)
	df_g6_reco.loc[df_g6_reco["case"]== "---","w"] =-1 *  df_g6_reco["w_reco"] .round(6)

#	print(df_g6_reco)

df_outlier = pd.DataFrame()

def find_outliers():
	""" Looks thorugh dataset to find cases where  u,v,w were recoverred incorrectly
	"""
	df_outlier =  df_g6_reco[(df_g6redprim["u"].round(6)!=df_g6_reco["u"].round(6))|(df_g6redprim["v"].round(6)!=df_g6_reco["v"].round(6))| \
	(df_g6redprim["w"].round(6)!=df_g6_reco["w"].round(6))]


	df_outlier_final = (pd.concat([df_outlier,df_cell_param.iloc[:,0:7]],axis = 1,join = "inner")).iloc[:,3:]
#	print(df_outlier_final)

def outputting_df():
	""" Outputs a dataframe showing a b c alpha beta gamma, original niggli reduced g6 vectors and 
	the recovered g6 vectors 
	"""
	global df_output 
	df_output = pd.concat([df_cell_param,df_g6redprim.iloc[:,0:6],df_dc7,df_tau,df_g6_reco["case"],df_g6_reco.iloc[:,10:13]],axis = 1) 

	df_output.to_csv("/home/mxerri/SULI/cellParamsG6RecoveryCompleted.csv",index = False)


#function/method calls

g6_to_dc7()
find_outliers()
outputting_df()
