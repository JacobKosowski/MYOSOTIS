import numpy as np
from datetime import datetime, timedelta
import params_clean as params
import os

init_time = 0
loop1_time = 0
loop2_time = 0
sedload_time = 0
data_save_time = 0
total_time = 0

#2D Runtime models
#Feather/Parallel SEDs
Ptotal_ims = [7.31211070e-05, 5.81019865e-03, 2.92794206e+01]
Ploop1_ims = [1.11221497e-04, 6.86773320e+00]
Psed_ims   = [9.29250349e-05, 8.86136853e+00]
Ploop2_ims = [7.32512615e-05, 5.11335343e-03, 1.33812489e+01]

Ptotal_nst = [7.47156672e-03, 2.47731275e+01]
Ploop1_nst = [4.01540759e-05, 6.70256429e+00]
Psed_nst   = [4.13868111e-04, 4.73288399e+00]
Ploop2_nst = [6.99843293e-03, 1.32249851e+01]



def time_estimate(p_ims,p_nst,nstars):
	imsize = np.sqrt(params.xpix*params.ypix)
	return (np.polyval(p_ims,imsize)+np.polyval(p_nst,nstars))/2

def prediction(nstars):
	print("-----------------------------------------------------------------")
	print("Runtime estimates for "+str(nstars)+" stars and "+str(params.xpix)+"x"+str(params.ypix)+" pixels")

	print("Loop 1:", "{:.5}".format(time_estimate(Ploop1_ims,Ploop1_nst,nstars)/60), "[min]")

	print("Loading SEDs:", "{:.5}".format(time_estimate(Psed_ims,Psed_nst,nstars)/60), "[min]")

	print("Loop 2:", "{:.5}".format(time_estimate(Ploop2_ims,Ploop2_nst,nstars)/60), "[min]")

	print("Total:", "{:.5}".format(time_estimate(Ptotal_ims,Ptotal_nst,nstars)/60), "[min]")

	print("Estimated time of code completion:", (datetime.now() + timedelta(seconds=time_estimate(Ptotal_ims,Ptotal_nst,nstars))).strftime("%H:%M:%S %d/%m/%Y"))
	print("-----------------------------------------------------------------")

def profiling():
	print("-----------------------------------------------------------------")
	print("Initiation:", "\t\t","{:.5}".format(init_time/60), "[min]")
	print("Loop 1:", "\t\t","{:.5}".format(loop1_time/60), "[min]")
	print("Loading SEDs ("+str(params.Ncpus)+" cpus):","\t","{:.5}".format(sedload_time/60), "[min]")
	print("Loop 2:", "\t\t","{:.5}".format(loop2_time/60), "[min]")
	print("Saving data to file:","\t","{:.5}".format(data_save_time/60), "[min]")
	print("Total:", "\t\t\t","{:.5}".format(total_time/60), "[min]")
	print("Time of code completion:", datetime.now().strftime("%H:%M:%S"))
	print("-----------------------------------------------------------------")

def basic_profile():
	print("-----------------------------------------------------------------")
	print('Simulation time: ', "{:.5}".format(total_time/60.) ,'[min]')
	print("Time of code completion:", datetime.now().strftime("%H:%M:%S"))
	print("-----------------------------------------------------------------")

def file_write(fn,data):
	    data = np.array(data)
	    f = open(fn, "a")

	    for r in data:
	        f.write(str(r))
	        f.write("\t")
	    f.write("\n")

def profile_output_vs_nstars(filename,nfovstars):

	if not os.path.exists("Runtime Profiles"):
		os.makedirs("Runtime Profiles")

	run = [nfovstars,total_time,init_time,loop1_time,sedload_time,loop2_time,data_save_time]

	file_write("Runtime Profiles/"+filename,run)

def profile_output_vs_imsize(filename):

	if not os.path.exists("Runtime Profiles"):
		os.makedirs("Runtime Profiles")

	imsizesqrt = np.sqrt(params.xpix*params.ypix)
	run = [imsizesqrt,total_time,init_time,loop1_time,sedload_time,loop2_time,data_save_time]
	
	file_write("Runtime Profiles/"+filename,run)

