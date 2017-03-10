#!/usr/bin/env python
# -------------------------------------------------------------------
# Python code to read and plot the following quantities
#  vertical wind speed
#  effective buoyancy
#  geostrophic imbalance
#  hydrostatic imbalance
#  Source of vertical momentum flux
#  Zonal wind
#  Medidional wind
#  density pert (r_prime)
#  b_prime
#  tracer
#
# This is repeated for a specified range of data directories
#
# Ross Bannister, Dec 2016
# -------------------------------------------------------------------

# ===================================================================
def max_field (field):
  # Subroutine to return maxmimum absolute value of field
  maxlev  = []
  for lev in field:
    maxlev.append(max(lev))
  globalmax = max(maxlev)
  return globalmax

def min_field (field):
  # Subroutine to return maxmimum absolute value of field
  minlev  = []
  for lev in field:
    minlev.append(min(lev))
  globalmin = min(minlev)
  return globalmin

def remove_extremes (field, low, high):
  # Set values less than low to low
  # Set values larger than high to high
  new_field = []
  for line in field:
    new_line = []
    for val in line:
      if val < low: 
        new_line.append(low)
      else:
        if val > high:
          new_line.append(high)
        else:
          new_line.append(val)
    new_field.append(new_line)
  return new_field

import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from matplotlib import colors, cm

# Set directories (each is a different experiment)
data_dir = []
data_dir.append('/export/diamet/raid1/ross/RuthsModel/ReferenceRun')
data_dir.append('/export/diamet/raid1/ross/RuthsModel/ReducedA')
data_dir.append('/export/diamet/raid1/ross/RuthsModel/ReducedB')
data_dir.append('/export/diamet/raid1/ross/RuthsModel/ReducedC')
data_dir.append('/export/diamet/raid1/ross/RuthsModel/IncreasedA')
data_dir.append('/export/diamet/raid1/ross/RuthsModel/IncreasedB')
data_dir.append('/export/diamet/raid1/ross/RuthsModel/IncreasedC')

datafile = 'Fields.nc'

# The time step of interest
tstep    = 6

# Set the domain dimensions
nlongs   = 360
nlevs    = 60


for datadir in data_dir:

  plotdir = datadir
  print 'INPUT DIRECTORY  : ', datadir
  print 'PLOT DIRECTORY   : ', plotdir

  # Deal with the vertical wind
  # -------------------------------------------------------------------
  # For this, we plot -0.4 to 0.4, but mask out anything above 1, below -1

  nc_file = Dataset(datadir + '/' + datafile)
  print '-------------------'
  #print nc_file
  field      = nc_file.variables['w'][tstep,:,:]
  longs_v    = nc_file.variables['longs_v'][:] / 1000.0
  full_level = nc_file.variables['full_level'][:] / 1000.0
  nc_file.close
  print 'w : ', field.shape
  minlev = -0.4
  maxlev =  0.4
  con_levs = np.linspace(minlev, maxlev, 11)
  minfield = min_field(field)
  maxfield = max_field(field)
  print 'Min value of this field ', minfield
  print 'Max value of this field ', maxfield
  mn       = np.mean(field)
  err      = field - mn
  rms      = np.sqrt(np.mean(err * err))
  print 'RMS of this field       ', rms
  print 'Contour levels: ', con_levs
  # These two lines are tricks for making sure that the scales remain the same for all plots
  field[0][0] = minlev ; field[0][1] = maxlev
  field = remove_extremes (field, minlev, maxlev)

  fig, ax = plt.subplots()
  cmap    = cm.get_cmap('seismic', len(con_levs))
  cax     = ax.contourf(longs_v, full_level, field, clevs=con_levs, cmap=cmap, extend='both')
  cbar    = fig.colorbar(cax, orientation='vertical', ticks=con_levs, cmap=cmap)
  cax     = ax.contour (longs_v, full_level, field, clevs=con_levs, colors='k')
  # Labels etc
  ax.set_title('w for t = ' + str(tstep*1800) + 's\nmin=' + str.format('%.2f' % minfield) + ', max=' + str.format('%.2f' % maxfield) + ', mean=' + str.format('%.5f' % mn) + ', rms=' + str.format('%.5f' % rms))
  ax.set_xlabel('Longitudinal distance (km)')
  ax.set_ylabel('Vertical distance (km)')
  #plt.show()
  plt.savefig(plotdir + '/w' + str(tstep) + '.eps', bbox_inches='tight')
  plt.close


  # Deal with the effective buoyancy
  # -------------------------------------------------------------------
  # For this we use the autoscale

  nc_file = Dataset(datadir + '/' + datafile)
  print '-------------------'
  #print nc_file
  field      = nc_file.variables['b_effective'][tstep,:,:]
  longs_v    = nc_file.variables['longs_v'][:] / 1000.0
  full_level = nc_file.variables['full_level'][:] / 1000.0
  nc_file.close
  print ' b_effective: ', field.shape
  #minlev = -0.04
  #maxlev =  0.04
  #con_levs = np.linspace(minlev, maxlev, 11)
  minfield = min_field(field)
  maxfield = max_field(field)
  print 'Min value of this field ', minfield
  print 'Max value of this field ', maxfield
  mn       = np.mean(field)
  err      = field - mn
  rms      = np.sqrt(np.mean(err * err))
  print 'mean of this field      ', mn
  print 'RMS of this field       ', rms
  #print 'Contour levels: ', con_levs

  fig, ax = plt.subplots()
  cmap    = cm.get_cmap('Greys', 11)
  cax     = ax.contourf(longs_v, full_level, field, cmap=cmap)
  cbar    = fig.colorbar(cax, orientation='vertical', cmap=cmap)
  cax     = ax.contour (longs_v, full_level, field, colors='k')
  # Labels etc
  ax.set_title(r'$\partial b^{\prime}/\partial z + A^2$ for t = ' + str(tstep*1800) + 's\nmin=' + str.format('%.2f' % minfield) + ', max=' + str.format('%.2f' % maxfield) + ', mean=' + str.format('%.5f' % mn) + ', rms=' + str.format('%.5f' % rms))
  ax.set_xlabel('Longitudinal distance (km)')
  ax.set_ylabel('Vertical distance (km)')
  #plt.show()
  plt.savefig(plotdir + '/b_effective' + str(tstep) + '.eps', bbox_inches='tight')
  plt.close



  # Deal with the geostrophic imbalance
  # -------------------------------------------------------------------
  # For this we plot -1 to 1, but mask out anything above 1, below -1

  nc_file = Dataset(datadir + '/' + datafile)
  print '-------------------'
  #print nc_file
  field      = nc_file.variables['geo_imbal'][tstep,:,:]
  longs_v    = nc_file.variables['longs_v'][:] / 1000.0
  half_level = nc_file.variables['half_level'][:] / 1000.0
  nc_file.close
  print ' geo_imbal: ', field.shape
  minlev = -1.0
  maxlev =  1.0
  con_levs = np.linspace(minlev, maxlev, 11)
  minfield = min_field(field)
  maxfield = max_field(field)
  print 'Min value of this field ', minfield
  print 'Max value of this field ', maxfield
  mn       = np.mean(field)
  err      = field - mn
  rms      = np.sqrt(np.mean(err * err))
  print 'mean of this field      ', mn
  print 'RMS of this field       ', rms
  print 'Contour levels: ', con_levs
  # These two lines are tricks for making sure that the scales remain the same for all plots
  field[0][0] = minlev ; field[0][1] = maxlev
  field = remove_extremes (field, minlev, maxlev)

  fig, ax = plt.subplots()
  cmap    = cm.get_cmap('seismic', len(con_levs))
  cax     = ax.contourf(longs_v, half_level, field, clevs=con_levs, cmap=cmap, extend='both')
  cbar    = fig.colorbar(cax, orientation='vertical', ticks=con_levs, cmap=cmap)
  cax     = ax.contour (longs_v, half_level, field, clevs=con_levs, colors='k')
  # Labels etc
  ax.set_title('Geostrophic imbalance for t = ' + str(tstep*1800) + 's\nmin=' + str.format('%.2f' % minfield) + ', max=' + str.format('%.2f' % maxfield) + ', mean=' + str.format('%.5f' % mn) + ', rms=' + str.format('%.5f' % rms))
  ax.set_xlabel('Longitudinal distance (km)')
  ax.set_ylabel('Vertical distance (km)')
  #plt.show()
  plt.savefig(plotdir + '/geo_imbal' + str(tstep) + '.eps', bbox_inches='tight')
  plt.close



  # Deal with the hydrostatic imbalance
  # -------------------------------------------------------------------
  # For this we plot -1 to 1, but mask out anything above 1, below -1

  nc_file = Dataset(datadir + '/' + datafile)
  print '-------------------'
  #print nc_file
  field      = nc_file.variables['hydro_imbal'][tstep,:,:]
  longs_v    = nc_file.variables['longs_v'][:] / 1000.0
  full_level = nc_file.variables['full_level'][:] / 1000.0
  nc_file.close
  print ' hydro_imbal: ', field.shape
  minlev = -1.0
  maxlev =  1.0
  con_levs = np.linspace(minlev, maxlev, 11)
  minfield = min_field(field)
  maxfield = max_field(field)
  print 'Min value of this field ', minfield
  print 'Max value of this field ', maxfield
  mn       = np.mean(field)
  err      = field - mn
  rms      = np.sqrt(np.mean(err * err))
  print 'mean of this field      ', mn
  print 'RMS of this field       ', rms
  print 'Contour levels: ', con_levs
  # These two lines are tricks for making sure that the scales remain the same for all plots
  field[0][0] = minlev ; field[0][1] = maxlev
  field = remove_extremes (field, minlev, maxlev)

  fig, ax = plt.subplots()
  cmap    = cm.get_cmap('seismic', 9) #len(con_levs))
  cax     = ax.contourf(longs_v, full_level, field, cmap=cmap) #clevs=con_levs, cmap=cmap, extend='both')
  cbar    = fig.colorbar(cax, orientation='vertical', cmap=cmap) #ticks=con_levs, cmap=cmap)
  cax     = ax.contour (longs_v, full_level, field, colors='k') #clevs=con_levs, colors='k')
  # Labels etc
  ax.set_title('Hydrostatic imbalance for t = ' + str(tstep*1800) + 's\nmin=' + str.format('%.2f' % minfield) + ', max=' + str.format('%.2f' % maxfield) + ', mean=' + str.format('%.5f' % mn) + ', rms=' + str.format('%.5f' % rms))
  ax.set_xlabel('Longitudinal distance (km)')
  ax.set_ylabel('Vertical distance (km)')
  #plt.show()
  plt.savefig(plotdir + '/hydro_imbal' + str(tstep) + '.eps', bbox_inches='tight')
  plt.close

  
  # Deal with the source of vertical momentum flux
  # -------------------------------------------------------------------
  # For this we use the autoscale, but with same size of +/- limits

  nc_file = Dataset(datadir + '/' + datafile)
  print '-------------------'
  #print nc_file
  field      = nc_file.variables['wmom_source'][tstep,:,:]
  longs_v    = nc_file.variables['longs_v'][:] / 1000.0
  full_level = nc_file.variables['full_level'][:] / 1000.0
  nc_file.close
  print ' wmom_source: ', field.shape
  #minlev = -1.0
  #maxlev =  1.0
  #con_levs = np.linspace(minlev, maxlev, 11)
  minfield = min_field(field)
  maxfield = max_field(field)
  print 'Min value of this field ', minfield
  print 'Max value of this field ', maxfield
  minlev   = -abs(max(abs(minfield), abs(maxfield)))
  maxlev   = abs(max(abs(minfield), abs(maxfield)))
  con_levs = np.linspace(minlev, maxlev, 11)
  mn       = np.mean(field)
  err      = field - mn
  rms      = np.sqrt(np.mean(err * err))
  print 'mean of this field      ', mn
  print 'RMS of this field       ', rms
  #print 'Contour levels: ', con_levs
  # These two lines are tricks for making sure that the scales remain the same for all plots
  field[0][0] = minlev ; field[0][1] = maxlev
  field = remove_extremes (field, minlev, maxlev)
  fig, ax = plt.subplots()
  cmap    = cm.get_cmap('seismic', len(con_levs))
  cax     = ax.contourf(longs_v, full_level, field, clevs=con_levs, cmap=cmap, extend='both')
  cbar    = fig.colorbar(cax, orientation='vertical', ticks=con_levs, cmap=cmap)
  cax     = ax.contour (longs_v, full_level, field, clevs=con_levs, colors='k')
  # Labels etc
  ax.set_title('Vertical momentum source for t = ' + str(tstep*1800) + 's\nmin=' + str.format('%.2f' % minfield) + ', max=' + str.format('%.2f' % maxfield) + ', mean=' + str.format('%.5f' % mn) + ', rms=' + str.format('%.5f' % rms))
  ax.set_xlabel('Longitudinal distance (km)')
  ax.set_ylabel('Vertical distance (km)')
  #plt.show()
  plt.savefig(plotdir + '/wmom_source' + str(tstep) + '.eps', bbox_inches='tight')
  plt.close


  # Deal with the zonal wind
  # -------------------------------------------------------------------
  # For this we use the autoscale

  nc_file = Dataset(datadir + '/' + datafile)
  print '-------------------'
  #print nc_file
  field      = nc_file.variables['u'][tstep,:,:]
  longs_u    = nc_file.variables['longs_u'][:] / 1000.0
  half_level = nc_file.variables['half_level'][:] / 1000.0
  nc_file.close
  print ' u: ', field.shape
  #minlev = -0.04
  #maxlev =  0.04
  #con_levs = np.linspace(minlev, maxlev, 11)
  minfield = min_field(field)
  maxfield = max_field(field)
  print 'Min value of this field ', minfield
  print 'Max value of this field ', maxfield
  mn       = np.mean(field)
  err      = field - mn
  rms      = np.sqrt(np.mean(err * err))
  print 'mean of this field      ', mn
  print 'RMS of this field       ', rms
  #print 'Contour levels: ', con_levs

  fig, ax = plt.subplots()
  cmap    = cm.get_cmap('Greys', 11)
  cax     = ax.contourf(longs_u, half_level, field, cmap=cmap)
  cbar    = fig.colorbar(cax, orientation='vertical', cmap=cmap)
  cax     = ax.contour (longs_u, half_level, field, colors='k')
  # Labels etc
  ax.set_title('u for t = ' + str(tstep*1800) + 's\nmin=' + str.format('%.2f' % minfield) + ', max=' + str.format('%.2f' % maxfield) + ', mean=' + str.format('%.5f' % mn) + ', rms=' + str.format('%.5f' % rms))
  ax.set_xlabel('Longitudinal distance (km)')
  ax.set_ylabel('Vertical distance (km)')
  #plt.show()
  plt.savefig(plotdir + '/u' + str(tstep) + '.eps', bbox_inches='tight')
  plt.close



  # Deal with the meridional wind
  # -------------------------------------------------------------------
  # For this we use the autoscale

  nc_file = Dataset(datadir + '/' + datafile)
  print '-------------------'
  #print nc_file
  field      = nc_file.variables['v'][tstep,:,:]
  longs_v    = nc_file.variables['longs_v'][:] / 1000.0
  half_level = nc_file.variables['half_level'][:] / 1000.0
  nc_file.close
  print ' v: ', field.shape
  #minlev = -0.04
  #maxlev =  0.04
  #con_levs = np.linspace(minlev, maxlev, 11)
  minfield = min_field(field)
  maxfield = max_field(field)
  print 'Min value of this field ', minfield
  print 'Max value of this field ', maxfield
  mn       = np.mean(field)
  err      = field - mn
  rms      = np.sqrt(np.mean(err * err))
  print 'mean of this field      ', mn
  print 'RMS of this field       ', rms
  #print 'Contour levels: ', con_levs

  fig, ax = plt.subplots()
  cmap    = cm.get_cmap('Greys', 11)
  cax     = ax.contourf(longs_v, half_level, field, cmap=cmap)
  cbar    = fig.colorbar(cax, orientation='vertical', cmap=cmap)
  cax     = ax.contour (longs_v, half_level, field, colors='k')
  # Labels etc
  ax.set_title('v for t = ' + str(tstep*1800) + 's\nmin=' + str.format('%.2f' % minfield) + ', max=' + str.format('%.2f' % maxfield) + ', mean=' + str.format('%.5f' % mn) + ', rms=' + str.format('%.5f' % rms))
  ax.set_xlabel('Longitudinal distance (km)')
  ax.set_ylabel('Vertical distance (km)')
  #plt.show()
  plt.savefig(plotdir + '/v' + str(tstep) + '.eps', bbox_inches='tight')
  plt.close



  # Deal with r_prime
  # -------------------------------------------------------------------
  # For this we use the autoscale

  nc_file = Dataset(datadir + '/' + datafile)
  print '-------------------'
  #print nc_file
  field      = nc_file.variables['r_prime'][tstep,:,:]
  longs_v    = nc_file.variables['longs_v'][:] / 1000.0
  half_level = nc_file.variables['half_level'][:] / 1000.0
  nc_file.close
  print ' r_prime: ', field.shape
  #minlev = -0.04
  #maxlev =  0.04
  #con_levs = np.linspace(minlev, maxlev, 11)
  minfield = min_field(field)
  maxfield = max_field(field)
  print 'Min value of this field ', minfield
  print 'Max value of this field ', maxfield
  mn       = np.mean(field)
  err      = field - mn
  rms      = np.sqrt(np.mean(err * err))
  print 'mean of this field      ', mn
  print 'RMS of this field       ', rms
  #print 'Contour levels: ', con_levs

  fig, ax = plt.subplots()
  cmap    = cm.get_cmap('Greys', 11)
  cax     = ax.contourf(longs_v, half_level, field, cmap=cmap)
  cbar    = fig.colorbar(cax, orientation='vertical', cmap=cmap)
  cax     = ax.contour (longs_v, half_level, field, colors='k')
  # Labels etc
  ax.set_title(r'$r^{\prime}$ for t = ' + str(tstep*1800) + 's\nmin=' + str.format('%.2f' % minfield) + ', max=' + str.format('%.2f' % maxfield) + ', mean=' + str.format('%.5f' % mn) + ', rms=' + str.format('%.5f' % rms))
  ax.set_xlabel('Longitudinal distance (km)')
  ax.set_ylabel('Vertical distance (km)')
  #plt.show()
  plt.savefig(plotdir + '/r_prime' + str(tstep) + '.eps', bbox_inches='tight')
  plt.close


  # Deal with b_prime
  # -------------------------------------------------------------------
  # For this we use the autoscale

  nc_file = Dataset(datadir + '/' + datafile)
  print '-------------------'
  #print nc_file
  field      = nc_file.variables['b_prime'][tstep,:,:]
  longs_v    = nc_file.variables['longs_v'][:] / 1000.0
  full_level = nc_file.variables['full_level'][:] / 1000.0
  nc_file.close
  print ' b_prime: ', field.shape
  #minlev = -0.04
  #maxlev =  0.04
  #con_levs = np.linspace(minlev, maxlev, 11)
  minfield = min_field(field)
  maxfield = max_field(field)
  print 'Min value of this field ', minfield
  print 'Max value of this field ', maxfield
  mn       = np.mean(field)
  err      = field - mn
  rms      = np.sqrt(np.mean(err * err))
  print 'mean of this field      ', mn
  print 'RMS of this field       ', rms
  #print 'Contour levels: ', con_levs

  fig, ax = plt.subplots()
  cmap    = cm.get_cmap('Greys', 11)
  cax     = ax.contourf(longs_v, full_level, field, cmap=cmap)
  cbar    = fig.colorbar(cax, orientation='vertical', cmap=cmap)
  cax     = ax.contour (longs_v, full_level, field, colors='k')
  # Labels etc
  ax.set_title(r'$b^{\prime}$ for t = ' + str(tstep*1800) + 's\nmin=' + str.format('%.2f' % minfield) + ', max=' + str.format('%.2f' % maxfield) + ', mean=' + str.format('%.5f' % mn) + ', rms=' + str.format('%.5f' % rms))
  ax.set_xlabel('Longitudinal distance (km)')
  ax.set_ylabel('Vertical distance (km)')
  #plt.show()
  plt.savefig(plotdir + '/b_prime' + str(tstep) + '.eps', bbox_inches='tight')
  plt.close


  # Deal with the tracer
  # -------------------------------------------------------------------
  # For this we use the autoscale

  nc_file = Dataset(datadir + '/' + datafile)
  print '-------------------'
  #print nc_file
  field      = nc_file.variables['tracer'][tstep,:,:]
  field0     = nc_file.variables['tracer'][0,:,:]          # Initial tracer
  longs_v    = nc_file.variables['longs_v'][:] / 1000.0
  half_level = nc_file.variables['half_level'][:] / 1000.0
  nc_file.close
  print ' tracer: ', field.shape
  #minlev = -0.04
  #maxlev =  0.04
  #con_levs = np.linspace(minlev, maxlev, 11)
  minfield = min_field(field)
  maxfield = max_field(field)
  print 'Min value of this field ', minfield
  print 'Max value of this field ', maxfield
  mn       = np.mean(field)
  err      = field - mn
  rms      = np.sqrt(np.mean(err * err))
  print 'mean of this field      ', mn
  print 'RMS of this field       ', rms
  #print 'Contour levels: ', con_levs

  fig, ax = plt.subplots()
  cmap    = cm.get_cmap('Greys', 11)
  cax     = ax.contourf(longs_v, half_level, field, cmap=cmap)
  cbar    = fig.colorbar(cax, orientation='vertical', cmap=cmap)
  cax     = ax.contour (longs_v, half_level, field, colors='k')
  cax     = ax.contour (longs_v, half_level, field0, colors='k')
  # Labels etc
  ax.set_title('tracer for t = ' + str(tstep*1800) + 's\nmin=' + str.format('%.2f' % minfield) + ', max=' + str.format('%.2f' % maxfield) + ', mean=' + str.format('%.5f' % mn) + ', rms=' + str.format('%.5f' % rms))
  ax.set_xlabel('Longitudinal distance (km)')
  ax.set_ylabel('Vertical distance (km)')
  #plt.show()
  plt.savefig(plotdir + '/tracer' + str(tstep) + '.eps', bbox_inches='tight')
  plt.close
