import numpy as np
import matplotlib.pyplot as plt


plt.rc('text', usetex=True)
plt.rc('font', family='serif', size=22)

fileNumber = 6500

Lx = 2*np.pi
Ly = 2*np.pi

filename = './vortex_xy.%.6d' % fileNumber

print(filename)
data = np.loadtxt(filename)
fig, axs = plt.subplots(1,1)

l1, = plt.plot(data[(data[:,2] < 0),0]-Lx/2, data[(data[:,2] < 0),1]-Ly/2, 'o', markersize=5, color='red')
l2, = plt.plot(data[(data[:,2] > 0),0]-Lx/2, data[(data[:,2] > 0),1]-Ly/2, 'o', markersize=5,color='blue')

plt.xticks([], [])
plt.yticks([], [])

#axs.set_xticks(np.arange(-Lx/2, (Lx/2)+0.01, Lx/4))
#labels = ['$-\pi$', r'$-\pi/2$', r'$0$', r'$\pi/2$', r'$\pi$']
#axs.set_xticklabels(labels)
#axs.set_yticks(np.arange(-Ly/2, (Ly/2)+0.01,Ly/4))
#labels = ['$-\pi$', r'$-\pi/2$', r'$0$', r'$\pi/2$', r'$\pi$']
#axs.set_yticklabels(labels)
plt.xlim(-0.5*Lx, 0.5*Lx)
plt.ylim(-0.5*Ly, 0.5*Ly)
axs.set_box_aspect(1)
plt.tight_layout()
plt.savefig('vortex_xy.pdf')
plt.show()