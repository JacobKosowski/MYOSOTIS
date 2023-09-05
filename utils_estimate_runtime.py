import numpy as np
from datetime import datetime, timedelta
import sys
import profiling
import params_clean as params



# Provides an estimate for the runtime of MYOSOTIS with an expected 16 cores, parallelized.

# Leave arguments blank to use current size and number of stars from current params.

# Add arguments of size and number of stars for custom calculations.

# Size should be given in number of pixels of either the length or width (assumes square final result)
# Nstars should be given in thousands of stars (e.g. for 50,000 stars --> Nstars=50)

# Time returned is and APPROXIMATION based on a model of previous runs. It is meant to simply provide a ball-park estimate for the total runtime.





print("")
print("==============================================================")


if len(sys.argv)==1:

	xstar,_,_ = np.loadtxt(params.filestar,usecols=(0,1,2),unpack=True)
	nstars = len(xstar)

	est = profiling.runtime_estimate(params.xpix,nstars)

	print("Runtime estimate for "+str(nstars)+" stars and "+str(params.xpix)+"x"+str(params.ypix)+" pixels, using "+str(16)+" cpus")

else:

	est = (profiling.runtime_estimate(float(sys.argv[1]),float(sys.argv[2])))

	print("Runtime estimate for "+sys.argv[2]+" stars and "+sys.argv[1]+"x"+sys.argv[1]+" pixels, using "+str(16)+" cpus")

print("--------------------------------------------------------------")
print("Estimated Runtime:",est/60,"[min]")
print("Estimated time of code completion:", (datetime.now() + timedelta(seconds=est)).strftime("%H:%M:%S %d/%m/%Y"))
print("==============================================================")