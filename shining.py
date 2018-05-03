from scipy.integrate import odeint
%matplotlib inline
#import sdf
import matplotlib
import matplotlib as mpl
#matplotlib.use('agg')
import matplotlib.pyplot as plt
import numpy as np
from numpy import ma
from matplotlib import colors, ticker, cm
from matplotlib.mlab import bivariate_normal
from optparse import OptionParser
import os
from mpl_toolkits.mplot3d import Axes3D
import random
from mpl_toolkits import mplot3d
from matplotlib import rc
import matplotlib.transforms as mtransforms
import matplotlib.colors as mcolors 
import scipy.ndimage as ndimage

font = {'family'   : 'monospace',  
          'color'  : 'black',  
          'weight' : 'normal',  
          'size'   : 20,  
          }  

x=np.linspace(0.0,12.0,1001)
p=np.linspace(0,5,6)
x,p=np.meshgrid(x,p)
y=np.zeros_like(x)
d=np.zeros_like(x)
for j in range(6):
        y[j,:]=np.sin(x[j,:])*(j/4.0+1)
        d[j,:]=abs(np.cos(x[j,:]))
for i in np.arange(400,800,2):
        ax=plt.subplot()
        #### manifesting colorbar, changing label and axis properties ####
        #cbar=plt.colorbar(ticks=[-eee, -eee/2, 0, eee/2, eee])
        #cbar.set_label('Normalized electric field',fontdict=font)
        index=np.arange(i-50,i,2)
        #plt.scatter(x[index], y[index], c=d[index], s=abs(y[index])*20, cmap='rainbow', norm=mcolors.Normalize(vmin=0, vmax=1), edgecolors='None')
        #plt.scatter(x[:,index], y[:,index], c=np.arange(index.size)*1, s=np.arange(index.size)*1, cmap='rainbow', edgecolors='None')
        plt.scatter(x[:,index], y[:,index], c=np.tile(np.arange(index.size)*1, (x[:,0].size, 1)), s=np.tile(np.arange(index.size)*1, (x[:,0].size, 1)), cmap='rainbow', edgecolors='None')
        #plt.plot((t[index,:])/2/np.pi,np.sqrt(px[index,:]**2+py[index,:]**2+1),'--k',linewidth=2.5,label='No RR')
        #plt.legend(loc='upper right')
        #cbar=plt.colorbar(ticks=np.linspace(np.min(gamma), np.max(gamma), 5))
        #cbar=plt.colorbar(ticks=np.linspace(np.min(lgR), np.max(lgR), 5))
        #cbar.set_label(r'$R$', fontdict=font)#plt.xlim(47,53)
        #ax.xaxis.set_major_locator(ticker.MultipleLocator(0.5))
        plt.xlabel('x',fontdict=font)
        plt.ylabel('y',fontdict=font)
        plt.xticks(fontsize=20.0); plt.yticks(fontsize=20);
        plt.ylim(-3.025,3.025)
        plt.xlim(0,12)
        #plt.legend(loc='best')


        #plt.show()
        #plt.subplots_adjust(top=0.92, bottom=0.08, left=0.15, right=0.95, hspace=0.05, wspace=0.30)
        fig = plt.gcf()
        fig.set_size_inches(8, 6.5)
        #fig.set_size_inches(5, 4.5)
        fig.savefig('./shining/series'+str(i).zfill(4)+'.png',format='png',dpi=80)
        plt.close("all")
