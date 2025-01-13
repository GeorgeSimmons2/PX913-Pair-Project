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


pot  = dataset['Electric_potential'][:]

plt.imshow(pot)

plt.savefig('Potential.jpg',dpi = 600)
plt.close()




nx = len(E_x[0])
ny = len(E_x[:,0])

plt.imshow(E_x,origin = 'lower',cmap = 'RdBu',extent = [-1,1,-1,1])
cb = plt.colorbar()




# x_ticks = np.linspace(0,nx-1,int(nx/10))
# x_ticks_labels = np.linspace(-1,1,int(nx/10),)

# x_ticks = np.around(x_ticks, decimals=2)
# x_ticks_labels = np.around(x_ticks_labels, decimals=2)

# yticks = np.linspace(0,ny-1,int(ny/10))
# y_ticks_labels = np.linspace(-1,1,int(ny/10),)



# plt.xticks(x_ticks,x_ticks_labels)
# plt.yticks(yticks,y_ticks_labels)

cb.set_label(label = '$E_{x}$ (AU)',size= 12)
cb.ax.tick_params(labelsize=10)


plt.xlabel('X (arb. units)',fontsize= 12)
plt.ylabel('Y (arb. units)',fontsize= 12)

plt.title('Electric Field in X direction',fontsize = 15)

plt.savefig('E_x.jpg',dpi = 600)
plt.close()





plt.plot(x,y,'x',color = 'k',markersize = 1,linestyle = 'dotted',label = '$r_{init}$ = (0,0), $v_{init} = (0,0)$')
plt.xlabel('X (arb. units)',fontsize= 12)
plt.ylabel('Y (arb. units)',fontsize= 12)
plt.grid()
plt.legend(fancybox=True, framealpha=1, shadow=True, borderpad=1,fontsize = 10)
plt.gca().set_aspect('equal')
plt.tick_params(axis='both', which='major', labelsize=12, width=2)
plt.tick_params(axis='both', which='minor', length=4, color='k')
plt.title('Trajectory of a Particle for problem: ',fontsize = 15)
plt.minorticks_on()

plt.savefig('Trajectory.jpg',dpi = 600)

plt.close()