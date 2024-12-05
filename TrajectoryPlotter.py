from matplotlib import pyplot as plt
import netCDF4
import numpy as np


dataset = netCDF4.Dataset('traj.nc')

# line taken from the lecture 7

print('File contains the following variables: ')
for n, el in enumerate(dataset.variables.keys()):
    print(n, ': ', el)
    #print(dataset[el][:])



E_x  = dataset['x_Electric_field'][:]


x = dataset['x_trajectory'][:]
y = dataset['y_trajectory'][:]

plt.imshow(E_x,cmap = 'jet')
plt.savefig('E_x.jpg',dpi = 600)

plt.close()

plt.scatter(x,y)
plt.xlabel('X position (AU)')
plt.ylabel('Y position (AU)')

plt.savefig('Trajectory.jpg',dpi = 600)

plt.close()