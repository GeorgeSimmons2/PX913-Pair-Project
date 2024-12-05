from matplotlib import pyplot as plt
import netCDF4
import numpy as np


dataset = netCDF4.Dataset('traj.nc')

# line taken from the lecture 7

print('File contains the following variables: ')
for n, el in enumerate(dataset.variables.keys()):
    print(n, ': ', el)
