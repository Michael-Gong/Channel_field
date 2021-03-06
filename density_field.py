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
##below is for generating mid transparent colorbar
  c_red = matplotlib.colors.colorConverter.to_rgba('red')
  c_blue= matplotlib.colors.colorConverter.to_rgba('blue')
  c_white_trans = matplotlib.colors.colorConverter.to_rgba('white',alpha = 0.0)
  cmap_rb = matplotlib.colors.LinearSegmentedColormap.from_list('rb_cmap',[c_red,c_white_trans,c_blue],128) 
  cmap_br = matplotlib.colors.LinearSegmentedColormap.from_list('rb_cmap',[c_blue,c_white_trans,c_red],128) 
##end for transparent colorbar##
 
##below is for norm colorbar
  class MidpointNormalize(colors.Normalize):
    def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
        self.midpoint = midpoint
        colors.Normalize.__init__(self, vmin, vmax, clip)

    def __call__(self, value, clip=None):
        # I'm ignoring masked values and all kinds of edge cases to make a
        # simple example...
        x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
        return np.ma.masked_array(np.interp(value, x, y)) 
##end for norm colorbar####

  ######### Parameter you should set ###########
  start   =  1  # start time
  stop    =  43  # end time
  step    =  1  # the interval or step
  from_path='./Data_hosing_pi/'
  from_path1='./jpg_hosing_pi/'
  to_path  ='./jpg_hosing_pi/'

#  xx_2d_outer = np.loadtxt(from_path1+'xx_2d_outer.txt') 
#  yy_2d_outer = np.loadtxt(from_path1+'yy_2d_outer.txt') 
#  gg_2d_outer = np.loadtxt(from_path1+'gg_2d_outer.txt') 

#  xx_2d_iner = np.loadtxt(from_path1+'xx_2d_iner.txt') 
#  yy_2d_iner = np.loadtxt(from_path1+'yy_2d_iner.txt') 
#  gg_2d_iner = np.loadtxt(from_path1+'gg_2d_iner.txt') 
#  youwant = ['electron_x_px','electron_density','electron_en','electron_theta_en','ey'] #,'electron_ekbar']
  #youwant field  ex,ey,ez,bx,by,bz,ex_averaged,bx_averaged...
  #youwant Derived electron_density,electron_ekbar...
  #youwant dist_fn electron_x_px, electron_py_pz, electron_theta_en...
  #if (os.path.isdir('jpg') == False):
  #  os.mkdir('jpg')
  ######### Script code drawing figure ################
  for n in range(start,stop+step,step):
    #### header data ####
    data = sdf.read(from_path+str(n).zfill(4)+".sdf",dict=True)
    header=data['Header']
    time=header['time']
    x  = data['Grid/Grid_mid'].data[0]/1.0e-6
    print('ok')
    y  = data['Grid/Grid_mid'].data[1]/1.0e-6
    X, Y = np.meshgrid(x, y) 
    
    name = 'Subset_high_e_density'
  #  den = data['Derived/Number_Density/electron_outside'].data/denunit
    den = data['Derived/Number_Density/electron'].data/denunit
  #  den = den+data['Derived/Number_Density/electron_outside_no'].data/denunit
  #  den = den+data['Derived/Number_Density/electronin'].data/denunit
    den = den[:,:]
    
    if np.min(den.T) == np.max(den.T):
            continue
    levels = np.logspace(-0.3, 2, 51)/2
    den.T[den.T > 49.999]=49.999 
    plt.contourf(X, Y, den.T, levels=levels, norm=mcolors.LogNorm(vmin=levels.min(), vmax=levels.max()), cmap='viridis')
    #### manifesting colorbar, changing label and axis properties ####
    cbar=plt.colorbar(ticks=np.logspace(0.0, 2.0, 3))
    cbar.set_label(r'$n_e\ [n_c]$', fontdict=font)

    name = 'ey'
    ex = data['Electric Field/'+str.capitalize(name)].data/exunit
    ex = ex[:,:]
    if np.min(ex.T) == np.max(ex.T):
            continue
    levels = np.linspace(-35.5, 35.5, 41)
    ex.T[ex.T < -35.499]=-35.499
    ex.T[ex.T >  35.499]= 35.499
    plt.contourf(X, Y, ex.T, levels=levels, cmap=cmap_br, norm=mcolors.Normalize(vmin=levels.min(), vmax=levels.max()))
    #### manifesting colorbar, changing label and axis properties ####
    cbar=plt.colorbar(ticks=np.linspace(-35.5, 35.5, 5))
    cbar.set_label(r'$E_y\ [m_ec\omega_0/e]$',fontdict=font)        
    plt.ylim(-12,12)
    if n >= 250:
        plt.xlim(0+(n-249.0)*0.1,30+(n-249.0)*0.1)
    else:
        plt.xlim(0,200)
    plt.xlabel('X [$\lambda_0$]',fontdict=font)
    plt.ylabel('Y [$\lambda_0$]',fontdict=font)
    plt.xticks(fontsize=20); plt.yticks(fontsize=20);

    fig = plt.gcf()
    fig.set_size_inches(32, 8)
    fig.savefig(to_path+'field_charge_'+str(n).zfill(4)+'.png',format='png',dpi=80)
    plt.close("all")


    print('finised '+str(round(100.0*(n-start+step)/(stop-start+step),4))+'%')

