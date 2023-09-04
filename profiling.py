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


def runtime_estimate(size,Nstars,stage=3):
	"""Provides an estimate for the runtime of MYOSOTIS with an expected 16 cores, parallelized.

	Size should be given in number of pixels of either the length or width (assumes square final result)

	Nstars should be given in thousands of stars (e.g. for 50,000 stars --> Nstars=50)

	Each stage corresponds to a different section of the code
	Stage=0: 1st main loop
	Stage=1: Loading SEDs
	Stage=2: Building Image
	Stage=3: Total

	Time returned is and APPROXIMATION based on a model of previous runs. It is meant to simply provide a ball-park estimate for the total runtime."""

	X = size
	Y = Nstars
    
	if stage==0:
	    
	    C = [5.62184725e+00, 1.53117324e-05, 1.49842181e-01]
	    Z = C[0] + C[1]*X + C[2]*Y
	    
	elif stage==1:
	    
	    C = [1.54714513e+00, -5.58590566e-06,  1.08888029e-01]
	    Z = C[0] + C[1]*X + C[2]*Y
	    
	elif stage==2:
	    
	    C = [ 5.43866311e+01, 2.83294317e-02, 5.76937366e+00, -6.83704301e-06,
	         -2.89454956e-03, -1.08497326e-01, 7.57686955e-10, 1.56529921e-06,
	          1.36483786e-05, 6.29317268e-04]
	    Z = C[0] + C[1]*X + C[2]*Y + C[3]*X**2 + C[4]*X*Y + C[5]*Y**2 + C[6]*X**3 + C[7]*X**2*Y + C[8]*X*Y**2 + C[9]*Y**3
	    
	elif stage==3:

	    C = [ 6.19789079e+01, 2.81152808e-02, 5.99643648e+00, -6.69729698e-06,
	         -2.89315400e-03, -1.07575622e-01, 7.39577146e-10, 1.56617635e-06,
	          1.35789600e-05, 6.24962899e-04]
	    Z = C[0] + C[1]*X + C[2]*Y + C[3]*X**2 + C[4]*X*Y + C[5]*Y**2 + C[6]*X**3 + C[7]*X**2*Y + C[8]*X*Y**2 + C[9]*Y**3


	return Z

def prediction(nstars):
	print("-----------------------------------------------------------------")
	print("Runtime estimates for "+str(nstars)+" stars and "+str(params.xpix)+"x"+str(params.ypix)+" pixels, using "+str(params.Ncpus)+" cpus")

	print("Loop 1:\t", "{:.5}".format(runtime_estimate(params.xpix,nstars/1000,stage=0)/60), "[min]")

	print("Loading SEDs:\t", "{:.5}".format(runtime_estimate(params.xpix,nstars/1000,stage=1)/60), "[min]")

	print("Loop 2:\t", "{:.5}".format(runtime_estimate(params.xpix,nstars/1000,stage=2)/60), "[min]")

	print("Total:\t", "{:.5}".format(runtime_estimate(params.xpix,nstars/1000,stage=3)/60), "[min]")

	print("Estimated time of code completion:", (datetime.now() + timedelta(seconds=runtime_estimate(params.xpix,nstars/1000,stage=3))).strftime("%H:%M:%S %d/%m/%Y"))
	print("-----------------------------------------------------------------")

def profiling():
	print("-----------------------------------------------------------------")
	print("Runtimes using "+str(params.Ncpus)+" cpus")
	# print("Initiation:", "\t\t","{:.5}".format(init_time/60), "[min]")
	print("Loop 1:", "\t\t","{:.5}".format(loop1_time/60), "[min]")
	print("Loading SEDs","\t\t","{:.5}".format(sedload_time/60), "[min]")
	print("Loop 2:", "\t\t","{:.5}".format(loop2_time/60), "[min]")
	# print("Saving data to file:","\t","{:.5}".format(data_save_time/60), "[min]")
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

