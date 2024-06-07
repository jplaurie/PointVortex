import numpy as np
import matplotlib
#matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.animation import FFMpegWriter

plt.rc('text', usetex=True)
plt.rc('font', family='serif', size=22)

startFile = 0
endFile = 1000
R=1

outTime = 1.0

metadata = dict(title='Movie', artist='Jason Laurie')
writer = FFMpegWriter(fps=25,metadata=metadata)

fig, axs = plt.subplots(1,1)

x=np.linspace(-R,R,1000)

l1, = plt.plot([], [], 'o', markersize=5, color='red')
l2, = plt.plot([], [], 'o', markersize=5,color='blue')
l3  = plt.figtext(0.65,0.9, [])
l4 = plt.plot(x,np.sqrt(R**2-x**2),color='black')
l5 = plt.plot(x,-np.sqrt(R**2-x**2),color='black')
axs.set_xticks(np.arange(-R, R+0.01, R/2))
labels = ['$-R$', r'$-R/2$', r'$0$', r'$R/2$', r'$R$']
axs.set_xticklabels(labels)
axs.set_yticks(np.arange(-R, R+0.01,R/2))
labels = ['$-R$', r'$-R/2$', r'$0$', r'$R/2$', r'$R$']
axs.set_yticklabels(labels)
axs.set_aspect('equal')
plt.xlim(-R, R)
plt.ylim(-R, R)




with writer.saving(fig, "movie_disc.mp4", 300):
    for i in range(startFile,endFile,1):
        filename = '../src_c++/data/vortex_xy.%.5d' % i;
        print(filename)
        data = np.loadtxt(filename)
        time = 'time = %.1f' % i
        l1.set_data(data[(data[:,2] > 0),0], data[(data[:,2] > 0),1])
        l2.set_data(data[(data[:,2] < 0),0], data[(data[:,2] < 0),1])
        l3.set_text(time)
      
        writer.grab_frame()
       
