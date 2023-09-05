#!/usr/bin/env python


try:
	from os.path import basename
	import os, shutil
except:
	print("os or shutil not installed")

try:
	import numpy as np
except:
	print("numpy not installed")

try:
	import pandas as pd
except:
	print("pandas not installed")

try:
	from math import sqrt, log10, log
except:
	print("math not installed")

try:
	from scipy.linalg import expm
except:
	print("scipy not installed")

try:
	import timeit
except:
	print("timeit not installed")

try:
	from datetime import datetime
except:
	print("datetime not installed")

try:
	import h5py
except:
	print("h5py not installed")


try:
	from astropy.io import fits
except:
	print("astropy not installed")

try:
	from numba import jit, prange
except:
	print("numba not installed")

try:
	from multiprocessing import Pool, Array
except:
	print("multiprocessing not installed")

try:
	import pyarrow
except:
	print("pyarrow not installed (for Feather")


try:
	import params_clean as params
	import functions
	import directories
	import constants
	import parallel_functions
	import profile
except:
	print("Myosotis core codes not found")
else:
	print("All imports good")
