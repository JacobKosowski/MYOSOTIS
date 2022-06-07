import numpy as np
from datetime import datetime, timedelta

Ptotal = [1.17546443e-02, 1.98501547e+01]
Ploop1 = [4.36024844e-05, 6.82053610e+00]
Ploop2 = [7.05361571e-03, 1.35331203e+01]
Psed = [0.00463701, -0.61997875]

def prediction(nstars):
	print("----------------------------------------------------------")
	print("Runtime estimates for "+str(nstars)+" stars")
	print("Total:", "{:.3}".format(np.polyval(Ptotal,nstars)/60), "[min]")
	print("Loop 1:", "{:.3}".format(np.polyval(Ploop1,nstars)/60), "[min]")
	print("Loading SEDs:", "{:.3}".format(np.polyval(Psed,nstars)/60), "[min]")
	print("Loop 2:", "{:.3}".format(np.polyval(Ploop2,nstars)/60), "[min]")
	print("Estimated time of code completion:", (datetime.now() + timedelta(seconds=np.polyval(Ptotal,nstars))).strftime("%H:%M:%S"))
	print("----------------------------------------------------------")