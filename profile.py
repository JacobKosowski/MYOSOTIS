import numpy as np
from datetime import datetime, timedelta

init_time = 0
loop1_time = 0
loop2_time = 0
sedload_time = 0
data_save_time = 0
total_time = 0

Ptotal = [1.17546443e-02, 1.98501547e+01]
Ploop1 = [4.36024844e-05, 6.82053610e+00]
Ploop2 = [7.05361571e-03, 1.35331203e+01]
Psed = [0.00463701, -0.61997875]

def prediction(nstars):
	print("-----------------------------------------------------------------")
	print("Runtime estimates for "+str(nstars)+" stars")
	print("Loop 1:", "{:.5}".format(abs(np.polyval(Ploop1,nstars)/60)), "[min]")
	print("Loading SEDs:", "{:.5}".format(abs(np.polyval(Psed,nstars)/60)), "[min]")
	print("Loop 2:", "{:.5}".format(abs(np.polyval(Ploop2,nstars)/60)), "[min]")
	print("Total:", "{:.5}".format(abs(np.polyval(Ptotal,nstars)/60)), "[min]")
	print("Estimated time of code completion:", (datetime.now() + timedelta(seconds=np.polyval(Ptotal,nstars))).strftime("%H:%M:%S %d/%m/%Y"))
	print("-----------------------------------------------------------------")

def profiling():
	print("-----------------------------------------------------------------")
	print("Initiation:", "{:.5}".format(init_time/60), "[min]")
	print("Loop 1:","{:.5}".format(loop1_time/60), "[min]")
	print("Loading SEDs:","{:.5}".format(sedload_time/60), "[min]")
	print("Loop 2:","{:.5}".format(loop2_time/60), "[min]")
	print("Saving data to file:","{:.5}".format(data_save_time/60), "[min]")
	print("Total:","{:.5}".format(total_time/60), "[min]")
	print("Time of code completion:", datetime.now().strftime("%H:%M:%S"))
	print("-----------------------------------------------------------------")

def basic_profile():
	print("-----------------------------------------------------------------")
	print('Simulation time: ', "{:.5}".format(total_time/60.) ,'[min]')
	print("Time of code completion:", datetime.now().strftime("%H:%M:%S"))
	print("-----------------------------------------------------------------")

def profile_output(filename,nfovstars):
	run = [nfovstars,total_time,init_time,loop1_time,sedload_time,loop2_time,data_save_time]

	def file_write(fn,data):
	    data = np.array(data)
	    f = open(fn, "a")

	    for r in data:
	        f.write(str(r))
	        f.write("\t")
	    f.write("\n")

	file_write(filename,run)

