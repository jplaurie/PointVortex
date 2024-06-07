import numpy as np
import matplotlib.pyplot as plt


plt.rc('text', usetex=True)
plt.rc('font', family='serif', size=18)

startFile = 0
endFile = 100
Lx = 2*np.pi
Ly = 4*np.pi
yaspect = 10
outTime = 1.0
boxNumber = 40
boxWidth = Ly/(boxNumber)


box = np.zeros((boxNumber,3))

for i in range(boxNumber):
    box[i,0] = (i-boxNumber/2)*boxWidth
   

for i in range(startFile,endFile,1):
    filename = '../../src_c++/data/vortex_xy.%.5d' % i;
    print(filename)
    data = np.loadtxt(filename)
   
    for j in range(len(data)):
        
        boxIndex = int((data[j,1]+(Ly/2)-np.pi)/boxWidth) 
        print(boxIndex)
        if data[j,2] > 0:   
            box[boxIndex,1] += 1
        elif data[j,2] < 0:
            box[boxIndex,2] -= 1


fig, axs = plt.subplots(2,1,figsize=(10,6) )

axs[0].fill_between(box[:,0], box[:,1], 0)
axs[0].fill_between(box[:,0], box[:,2], 0)
axs[1].fill_between(box[:,0], box[:,1]+box[:,2], 0)

#axs.set_xticks(np.arange(-Lx/2, (Lx/2)+0.01, Lx/4))
#labels = ['$-\pi$', r'$-\pi/2$', r'$0$', r'$\pi/2$', r'$\pi$']
#axs.set_xticklabels(labels)
#axs.set_yticks(np.arange(-Ly/2, (Ly/2)+0.01,Ly/4))
#labels = ['$-\pi$', r'$-\pi/2$', r'$0$', r'$\pi/2$', r'$\pi$']
#axs.set_yticklabels(labels)
axs[0].set_ylabel(r'$PDF(\Gamma_1), PDF(\Gamma_2)$')
axs[0].set_xlim(-0.5*Ly, 0.5*Ly)
axs[0].set_xlabel(r'$y$')
axs[0].set_xticks(np.arange(-Ly/2, (Ly/2)+0.01,Ly/4))
#labels = ['$-\pi$', r'$-\pi/2$', r'$0$', r'$\pi/2$', r'$\pi$']
#axs[0].set_xticklabels(labels)

#plt.ylim(-yaspect*0.5*Ly, yaspect*0.5*Ly)

axs[1].set_xlim(-0.5*Ly, 0.5*Ly)
axs[1].set_xlabel(r'$y$')
axs[1].set_ylabel(r'$PDF(\Gamma_1+\Gamma_2)$')
axs[1].set_xticks(np.arange(-Ly/2, (Ly/2)+0.01,Ly/4))
#labels = ['$-\pi$', r'$-\pi/2$', r'$0$', r'$\pi/2$', r'$\pi$']
#axs[1].set_xticklabels(labels)



plt.tight_layout()
#plt.legend(fontsize=16, loc=1)

plt.savefig('distribution.pdf')
