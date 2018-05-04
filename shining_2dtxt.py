from scipy.integrate import odeint
#%matplotlib inline
import sdf
import matplotlib
import matplotlib as mpl
matplotlib.use('agg')
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
  from_path='./Data_new/'
  to_path='./jpg_shining/'
  data = sdf.read(from_path+"0299.sdf",dict=True)
  #grid_x = data['Grid/Particles/subset_high_e/electron'].data[0]/wavelength
  #work_x = data['Particles/Time_Integrated_Work_x/subset_high_e/electron'].data
  #work_y = data['Particles/Time_Integrated_Work_y/subset_high_e/electron'].data
  px = data['Particles/Px/subset_high_e/electron'].data/(m0*v0)
  py = data['Particles/Py/subset_high_e/electron'].data/(m0*v0)
  gg = (px**2+py**2+1)**0.5
  part13_id = data['Particles/ID/subset_high_e/electron'].data
  part13_id = part13_id[ (gg>2.0) ]
  print('part13_id size is ',part13_id.size,' max ',np.max(part13_id),' min ',np.min(part13_id))


  data = sdf.read(from_path+"0050.sdf",dict=True)
  grid_y = data['Grid/Particles/subset_high_e/electron'].data[1]/wavelength
  part00_id = data['Particles/ID/subset_high_e/electron'].data
  part00_id = part00_id[ ( abs(grid_y)<5 ) ]
  part13_id = np.intersect1d(part00_id,part13_id)
  print('after intersect 0050.sdf part_id size is ',part13_id.size,' max ',np.max(part13_id),' min ',np.min(part13_id))


  ######### Parameter you should set ###########
  start   =  50  # start time
  stop    =  299  # end time
  step    =  1  # the interval or step

  #  if (os.path.isdir('jpg') == False):
  #    os.mkdir('jpg')
  ######### Script code drawing figure ################
  for n in range(start,stop+step,step):
      #### header data ####
      data = sdf.read(from_path+str(n).zfill(4)+".sdf",dict=True)
      header=data['Header']
      time=header['time']
      if ( n==start ):
          part_id = data['Particles/ID/subset_high_e/electron'].data
      else:
          part_id = np.intersect1d(data['Particles/ID/subset_high_e/electron'].data, part_id)
      print('Particle_ID size is ',part_id.size,' max ',np.max(part_id),' min ',np.min(part_id))

  part_id = np.intersect1d(part_id,part13_id)
  print('After intersecting with 0013.sdf')
  print('Particle_ID size is ',part_id.size,' max ',np.max(part_id),' min ',np.min(part_id))

  #px_2d = np.zeros([part_id.size,stop-start+1])
  #py_2d = np.zeros([part_id.size,stop-start+1])
  xx_2d = np.zeros([part_id.size,stop-start+1])
  yy_2d = np.zeros([part_id.size,stop-start+1])
  gg_2d = np.zeros([part_id.size,stop-start+1])
  #work_x_2d = np.zeros([part_id.size,stop-start+1])
  #work_y_2d = np.zeros([part_id.size,stop-start+1])
  for n in range(start,stop+step,step):
      #### header data ####
      data = sdf.read(from_path+str(n).zfill(4)+".sdf",dict=True)
      px = data['Particles/Px/subset_high_e/electron'].data/(m0*v0)
      py = data['Particles/Py/subset_high_e/electron'].data/(m0*v0)
      gg = (px**2+py**2+1.0)**0.5
      grid_x = data['Grid/Particles/subset_high_e/electron'].data[0]/wavelength
      grid_y = data['Grid/Particles/subset_high_e/electron'].data[1]/wavelength
      #work_x = data['Particles/Time_Integrated_Work_x/subset_high_e/electron'].data
      #work_y = data['Particles/Time_Integrated_Work_y/subset_high_e/electron'].data
      temp_id = data['Particles/ID/subset_high_e/electron'].data

      #px = px[np.in1d(temp_id,part_id)]
      #py = py[np.in1d(temp_id,part_id)]
      gg     = gg[np.in1d(temp_id,part_id)]
      grid_x = grid_x[np.in1d(temp_id,part_id)]
      grid_y = grid_y[np.in1d(temp_id,part_id)]
      #work_x = work_x[np.in1d(temp_id,part_id)]
      #work_y = work_y[np.in1d(temp_id,part_id)]
      temp_id = temp_id[np.in1d(temp_id,part_id)]

      for ie in range(part_id.size):
          #px_2d[ie,n-start] = px[temp_id==part_id[ie]]
          #py_2d[ie,n-start] = py[temp_id==part_id[ie]]
          xx_2d[ie,n-start] = grid_x[temp_id==part_id[ie]]
          yy_2d[ie,n-start] = grid_y[temp_id==part_id[ie]]
          gg_2d[ie,n-start] = gg[temp_id==part_id[ie]]
          #work_x_2d[ie,n-start] = work_x[temp_id==part_id[ie]]
          #work_y_2d[ie,n-start] = work_y[temp_id==part_id[ie]]
      print('finised constructing 2darray'+str(round(100.0*(n-start+step)/(stop-start+step),4))+'%')
  np.savetxt('./jpg_shining/xx_2d.txt',xx_2d)
  np.savetxt('./jpg_shining/yy_2d.txt',yy_2d)
  np.savetxt('./jpg_shining/gg_2d.txt',gg_2d)



