#!/public/home/users/bio001/tools/python-2.7.11/bin/python
import sdf
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import numpy as np
import os
from numpy import ma
from matplotlib import colors, ticker, cm
from matplotlib.mlab import bivariate_normal
import matplotlib.colors as mcolors 
import scipy.ndimage as ndimage
 
if __name__ == "__main__":
  print ('This is main of module "test2d.py"')
  ######## Constant defined here ########
  pi        =     3.1415926535897932384626
  q0        =     1.602176565e-19 # C
  m0        =     9.10938291e-31  # kg
  v0        =     2.99792458e8    # m/s^2
  kb        =     1.3806488e-23   # J/K
  mu0       =     4.0e-7*np.pi       # N/A^2
  epsilon0  =     8.8541878176203899e-12 # F/m
  h_planck  =     6.62606957e-34  # J s
  wavelength=     1.0e-6
  frequency =     v0*2*pi/wavelength
  
  exunit    =     m0*v0*frequency/q0
  bxunit    =     m0*frequency/q0
  denunit    =     frequency**2*epsilon0*m0/q0**2
  jalf      =     4*np.pi*epsilon0*m0*v0**3/q0/wavelength**2
  print('electric field unit: '+str(exunit))
  print('magnetic field unit: '+str(bxunit))
  print('density unit nc: '+str(denunit))
  
  font = {'family' : 'monospace',  
          'color'  : 'black',  
          'weight' : 'normal',  
          'size'   : 20,  
          }  
# this is for constructing a transparent colorbar
  colors = [(1,0,0,c**1) for c in np.linspace(0,1,500)]
  cmapred = mcolors.LinearSegmentedColormap.from_list('mycmap', colors, N=256)
  colors = [(0,0,1,c**1) for c in np.linspace(0,1,500)]
  cmapblue = mcolors.LinearSegmentedColormap.from_list('mycmap', colors, N=256)
  series_c =  np.linspace(0,1,500)
  colors = [(0,0,0,c) for c in series_c**0.4]
  cmapblack = mcolors.LinearSegmentedColormap.from_list('mycmap', colors, N=256)
#from colour import Colo

  
##below is for norm colorbar
  class MidpointNormalize(mcolors.Normalize):
    def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
        self.midpoint = midpoint
        mcolors.Normalize.__init__(self, vmin, vmax, clip)

    def __call__(self, value, clip=None):
        # I'm ignoring masked values and all kinds of edge cases to make a
        # simple example...
        x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
        return np.ma.masked_array(np.interp(value, x, y)) 
##end for norm colorbar####

##below is for generating mid transparent colorbar
  c_orange = matplotlib.colors.colorConverter.to_rgba('orange')
  c_white_trans = matplotlib.colors.colorConverter.to_rgba('white',alpha = 0.0)
  cmap_orange = matplotlib.colors.LinearSegmentedColormap.from_list('ow_cmap',[c_white_trans,c_orange],128) 
  c_green = matplotlib.colors.colorConverter.to_rgba('green')
  c_white_trans = matplotlib.colors.colorConverter.to_rgba('white',alpha = 0.0)
  cmap_green = matplotlib.colors.LinearSegmentedColormap.from_list('gw_cmap',[c_white_trans,c_green],128) 
##end for transparent colorbar##
 
 

 
  ######### Parameter you should set ###########
  start   =  1  # start time
  stop    =  20  # end time
  step    =  1  # the interval or step
  
  #youwant field  ex,ey,ez,bx,by,bz,ex_averaged,bx_averaged...
  youwant =  ['jx','jy']#,'ey','ex','ey_averaged','bz','bz_averaged','electron_density','carbon_density']
  #youwant Derived electron_density,electron_ekbar...
  #youwant dist_fn electron_x_px, electron_py_pz, electron_theta_en...
  ######### Script code drawing figure ################
  from_path = './Data_m01/'
  to_path = './jpg_m01/'
  for n in range(start,stop+step,step):
    #### header data ####
    data = sdf.read(from_path+str(n).zfill(4)+".sdf",dict=True)
    header=data['Header']
    time=header['time']
    x  = data['Grid/Grid_mid'].data[0]/1.0e-6
    print('ok')
    y  = data['Grid/Grid_mid'].data[1]/1.0e-6
    y  = y[200:400]
    X, Y = np.meshgrid(x, y)
    
    for name in youwant:
        jj = data['Current/'+str.capitalize(name)].data/jalf
        jj = jj[:,200:400]
        positive = np.zeros_like(jj)
        negetive = np.zeros_like(jj)
        positive[jj>0] = jj[jj>0]
        negetive[jj<-0] = -jj[jj<-0]
        
        den=negetive        
        levels = np.logspace(-1,2,41)
        den[den > 0.999*np.max(levels)] =  0.999*np.max(levels)
        den[den < 1.001*np.min(levels)] =  1.001*np.min(levels)
        plt.contourf(X, Y, den.T, levels=levels, norm=mcolors.LogNorm(vmin=levels.min(), vmax=levels.max()), cmap=cmapred)
        #### manifesting colorbar, changing label and axis properties ####
        #cbar=plt.colorbar(ticks=np.linspace(-20.0, 20.0, 5))
        cbar=plt.colorbar(ticks=np.logspace(-1,2,4))
        cbar.set_label(r'negetive [$j_\alpha$]', fontdict=font)

        den=positive
        levels = np.logspace(-1,2,41)
        den[den > 0.999*np.max(levels)] =  0.999*np.max(levels)
        den[den < 1.001*np.min(levels)] =  1.001*np.min(levels)
        plt.contourf(X, Y, den.T, levels=levels, norm=mcolors.LogNorm(vmin=levels.min(), vmax=levels.max()), cmap=cmapblue)
        #### manifesting colorbar, changing label and axis properties ####
        #cbar=plt.colorbar(ticks=np.linspace(-20.0, 20.0, 5))
        cbar=plt.colorbar(ticks=np.logspace(-1,2,4))
        cbar.set_label(r'positive [$j_\alpha$]', fontdict=font)

        plt.xlabel('X [$\mu m$]',fontdict=font)
        plt.ylabel('Y [$\mu m$]',fontdict=font)
        plt.xticks(fontsize=20); plt.yticks(fontsize=20);
        plt.title(name+' at '+str(round(time/1.0e-15,6))+' fs',fontdict=font)
        fig = plt.gcf()
        fig.set_size_inches(12, 7)
        fig.savefig(to_path+name+str(n).zfill(4)+'.png',format='png',dpi=100)
        plt.close("all")
    print('finised '+str(round(100.0*(n-start+step)/(stop-start+step),4))+'%')
  
