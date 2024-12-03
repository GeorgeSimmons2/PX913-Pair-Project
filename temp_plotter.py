from matplotlib import pyplot as plt
import csv
import numpy as np


fname = 'potential.csv'



pot = []
with open(fname,'r') as f:
    reader = csv.reader(f)
    
    for row in reader:
        sub = []
        for element in row:
            sub.append(float(element))
        pot.append(sub)

plt.imshow(pot,cmap = 'jet')
cb = plt.colorbar(label = 'Potential (AU)')


plt.savefig('potential.png')
plt.close()