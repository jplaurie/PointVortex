import numpy as np
import matplotlib
#matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.animation import FFMpegWriter

plt.rc('text', usetex=True)
plt.rc('font', family='serif', size=22)

startFile = 500
endFile = 2000
Lx = 2*np.pi
Ly = 2*np.pi
yaspect = 1
outTime = 1.0

metadata = dict(title='Movie', artist='Jason Laurie')
writer = FFMpegWriter(fps=25,metadata=metadata)

fig, axs = plt.subplots(1,1)

l1, = plt.plot([], [], 'o', markersize=5, color='red')
l2, = plt.plot([], [], 'o', markersize=5,color='blue')
l3  = plt.figtext(0.65,0.9, [])

axs.set_xticks(np.arange(-Lx/2, (Lx/2)+0.01, Lx/4))
labels = ['$-\pi$', r'$-\pi/2$', r'$0$', r'$\pi/2$', r'$\pi$']
axs.set_xticklabels(labels)
axs.set_yticks(np.arange(-Ly/2, (Ly/2)+0.01,Ly/4))
labels = ['$-\pi$', r'$-\pi/2$', r'$0$', r'$\pi/2$', r'$\pi$']
axs.set_yticklabels(labels)
plt.xlim(-0.5*Lx, 0.5*Lx)
plt.ylim(-yaspect*0.5*Ly, yaspect*0.5*Ly)



with writer.saving(fig, "movie.mp4", 300):
    for i in range(startFile,endFile,1):
        filename = '../src_new/data/vortex_xy.%.5d' % i;
        print(filename)
        data = np.loadtxt(filename)
        time = 'time = %.1f' % i
        l1.set_data(data[(data[:,2] < 0),0], data[(data[:,2] < 0),1])
        l2.set_data(data[(data[:,2] > 0),0], data[(data[:,2] > 0),1])
        l3.set_text(time)
      
        writer.grab_frame()
       
