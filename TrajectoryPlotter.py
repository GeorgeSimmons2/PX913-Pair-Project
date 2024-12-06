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




nx = len(E_x[0])
ny = len(E_x[:,0])

plt.imshow(E_x,cmap = 'jet')
cb = plt.colorbar(label = '$E_{x}$ (AU)')




x_ticks = np.linspace(0,nx-1,int(nx/10))
x_ticks_labels = np.linspace(-1,1,int(nx/10),)

x_ticks = np.around(x_ticks, decimals=2)
x_ticks_labels = np.around(x_ticks_labels, decimals=2)

yticks = np.linspace(0,ny-1,int(ny/10))
y_ticks_labels = np.linspace(-1,1,int(ny/10),)



plt.xticks(x_ticks,x_ticks_labels)
plt.yticks(yticks,y_ticks_labels)



plt.xlabel('X (AU)')
plt.ylabel('Y (AU)')

plt.title('Electric Field in X direction')


plt.savefig('E_x.jpg',dpi = 600)
plt.close()





plt.plot(x,y,'x',color = 'm',markersize = 1)
plt.xlabel('X (AU)')
plt.ylabel('Y (AU)')
plt.grid()
plt.tight_layout()
plt.title('Trajectory of a Particle for problem: ')


plt.savefig('Trajectory.jpg',dpi = 600)

plt.close()