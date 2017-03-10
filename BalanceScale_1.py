#!/usr/bin/env python
# -------------------------------------------------------------------
# Python code to read the following quantities
#   b_prime
#   rho_prime
#   v
#   longs_v
#   half_level
#   full_level
# Do selective filtering of the fields
# Compute
#   Geostrophic imbalance
#   Hydrostatic imbalance
#
# This is repeated for a specified range of data directories
#
# Ross Bannister, Dec 2016
# -------------------------------------------------------------------


def filterfield (field_in, lengths, scale):
  # Subroutine to return filtered field (bandpass down to the given scale)
  import numpy as np

  if scale == 1.0:
    # No filtering required
    print 'No filtering required'
    field_out = field_in
  else:
    # FFT the field
    field_fft = np.fft.fft(field_in)

    # Modify the field
    print 'Filtering the field'
    for k in range(0, len(lengths)):
      modfactor       = np.arctan(lengths[k] - scale) / np.pi + 0.5
      field_fft[:,k] *= modfactor

    # Inverse FFT the field
    field_out = np.real(np.fft.ifft(field_fft))
  return field_out


def calc_geoimbal (rp, v, dist, f, C, nlongs, nlevs):
  # Calculate the geostrophic imbalance diagnostic
  import numpy as np

  # Set-up
  recipdx = 1.0 / (dist[1] - dist[0])
  gi_1     = np.zeros((nlevs,nlongs))
  gi_2     = np.zeros((nlevs,nlongs))

  # Calculate unnormalised geostrophic imbalance, terms 1 and 2
  for z in range(0,nlevs):
    for x in range(1,nlongs):
      gi_1[z,x] = C * recipdx * (rp[z,x] - rp[z,x-1])
      gi_2[z,x] = -1.0 * f * (v[z,x] + v[z,x-1]) / 2.0

  # Normalise
  # Calculate the rms of term 1
  mean    = np.mean(gi_1[:,1:nlongs])
  rms_1   = np.sqrt(np.mean( (gi_1[:,1:nlongs]-mean)*(gi_1[:,1:nlongs]-mean) ))
  # Calculate the rms of term 2
  mean    = np.mean(gi_2[:,1:nlongs])
  rms_2   = np.sqrt(np.mean( (gi_2[:,1:nlongs]-mean)*(gi_2[:,1:nlongs]-mean) ))
  gi      = (gi_1 + gi_2) / (rms_1 + rms_2)

  return gi


def calc_hydimbal (rp, bp, full_lev, half_lev, C, nlongs, nlevs):
  # Calculate the hydrostatic imbalance diagnostic
  import numpy as np

  # Set-up
  hi_1 = np.zeros((nlevs,nlongs))
  hi_2 = np.zeros((nlevs,nlongs))

  # Calculate unnormalised hydrostatic imbalance, terms 1 and 2
  for x in range(0,nlongs):
    for z in range(0,nlevs-1):
      hi_1[z,x] = C * (rp[z+1,x] - rp[z,x]) / (half_lev[z+1] - half_lev[z])
      hi_2[z,x] = -1.0 * bp[z,x]

  # Normalise
  # Calculate the rms of term 1
  mean    = np.mean(hi_1[0:nlevs,:])
  rms_1   = np.sqrt(np.mean( (hi_1[0:nlevs,:]-mean)*(hi_1[0:nlevs,:]-mean) ))
  # Calculate the rms of term 2
  mean    = np.mean(hi_2[0:nlevs,:])
  rms_2   = np.sqrt(np.mean( (hi_2[0:nlevs,:]-mean)*(hi_2[0:nlevs,:]-mean) ))
  hi      = (hi_1 + hi_2) / (rms_1 + rms_2)

  return hi




import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from matplotlib import colors, cm

# Set directories
base_dir = '/export/diamet/raid1/ross/RuthsModel'
data_dir = []                             ; param_C = []             ; linetwo = []              ; colour = []            ; label = []
data_dir.append(base_dir+'/ReferenceRun') ; param_C.append(10000.0)  ; linetwo.append('solid')   ; colour.append('k')	  ; label.append('Ref')
data_dir.append(base_dir+'/ReducedA')	 ; param_C.append(10000.0)  ; linetwo.append('dashed')  ; colour.append('0.5')   ; label.append('A-')
data_dir.append(base_dir+'/IncreasedA')   ; param_C.append(10000.0)  ; linetwo.append('dashed')  ; colour.append('k')	  ; label.append('A+')
data_dir.append(base_dir+'/ReducedB')	 ; param_C.append(10000.0)  ; linetwo.append('dotted')  ; colour.append('0.5')   ; label.append('B-')
data_dir.append(base_dir+'/IncreasedB')   ; param_C.append(10000.0)  ; linetwo.append('dotted')  ; colour.append('k')	  ; label.append('B+')
data_dir.append(base_dir+'/ReducedC')	 ; param_C.append(1000.0)   ; linetwo.append('dashdot') ; colour.append('0.5')   ; label.append('C-')
data_dir.append(base_dir+'/IncreasedC')   ; param_C.append(100000.0) ; linetwo.append('dashdot') ; colour.append('k')	  ; label.append('C+')

datafile      = 'Fields.nc'
f             = 0.0001

# The time step of interest
tstep_max     = 6
times         = np.linspace(0,float(tstep_max)/2.0,tstep_max+1)  # Times in hours

# Set the domain dimensions
nlongs        = 360
nlevs         = 60
Lx            = 1.5 * float(nlongs)

# Wavenumbers and scales
hori_wns      = np.linspace(0,nlongs-1,nlongs)
hori_lens     = np.zeros(nlongs)
hori_lens[1:nlongs/2+1] = Lx / hori_wns[1:nlongs/2+1]
hori_lens[0]  = Lx
for l in range(1, nlongs/2):
  hori_lens[nlongs-l] = hori_lens[l]

#for l in range(0,nlongs):
#  print l, hori_wns[l], hori_lens[l]


# Define the scales of interest (km)
scalelim = []                  ; lineone = []
scalelim.append(100.0)         ; lineone.append('solid')
scalelim.append(10.0)          ; lineone.append('dashed')
scalelim.append(1.0)           ; lineone.append('dotted')


# Set-up the arrays to store the results
print 'Number of directories ', len(data_dir)
print 'Number of time steps  ', tstep_max+1
print 'Number of scales      ', len(scalelim)
gimbalresults = np.zeros((len(data_dir), tstep_max+1, len(scalelim)))
himbalresults = np.zeros((len(data_dir), tstep_max+1, len(scalelim)))



for dirno in range(0,len(data_dir)):

  datadir = data_dir[dirno]
  print 'INPUT DIRECTORY  : ', datadir
  
  for time in range(0,tstep_max+1):
    print '  Dealing with timestep ', time

    # Read-in the relevant fields
    nc_file    = Dataset(datadir + '/' + datafile)
    print '  -------------------'
    bp         = nc_file.variables['b_prime'][time,:,:]
    rp         = nc_file.variables['r_prime'][time,:,:]
    v          = nc_file.variables['v'][time,:,:]
    longs_v    = nc_file.variables['longs_v'][:] / 1000.0
    half_level = nc_file.variables['half_level'][:] / 1000.0
    full_level = nc_file.variables['full_level'][:] / 1000.0
    nc_file.close

    for scale in range (0,len(scalelim)):
      print '    Dealing with lengthscale band ', scalelim[scale]

      # Filter the fields to remove specific scales
      print '    Filtering fields'
      bp_filt = filterfield(bp, hori_lens, scalelim[scale])
      rp_filt = filterfield(rp, hori_lens, scalelim[scale])
      v_filt  = filterfield(v,  hori_lens, scalelim[scale])
      
      #plt.figure(1)
      #plt.subplot(2,1,1)
      #plt.contourf(v)
      #plt.subplot(2,1,2)
      #plt.contourf(v_filt)
      #plt.show()

      # Calculate the diagnostics
      gimbal  = calc_geoimbal(rp_filt, v_filt, 1000.0*longs_v, f, param_C[dirno], nlongs, nlevs)
      himbal  = calc_hydimbal(rp_filt, bp_filt, 1000.0*full_level, 1000.0*half_level, param_C[dirno], nlongs, nlevs)

      mean       = np.mean(gimbal[:,1:nlongs])
      gimbal_rms = np.sqrt(np.mean( (gimbal[:,1:nlongs]-mean)*(gimbal[:,1:nlongs]-mean) ))
      
      mean       = np.mean(himbal[0:nlevs,:])
      himbal_rms = np.sqrt(np.mean( (himbal[0:nlevs,:]-mean)*(himbal[0:nlevs,:]-mean) ))


      # Store the diagnostics
      gimbalresults[dirno,time,scale] = gimbal_rms
      himbalresults[dirno,time,scale] = himbal_rms



# Plot the results - separate plots for each experiment
for dirno in range(0,len(data_dir)):

  plotdir = data_dir[dirno]



  # Plot the geostrophic imbalance
  plotfile = plotdir + '/Imbal.eps'

  fig, ax1 = plt.subplots()
  ax2      = ax1.twinx()

  print 'Plotting geostrophic imbalance'
  for scale in range (0,len(scalelim)):
    ax1.plot(times, gimbalresults[dirno,:,scale], color='black', linewidth='2', ls=lineone[scale], label='>' + str(scalelim[scale]) + 'km')
  ax1.legend(loc='upper left')
  ax1.set_xlabel('time (hours)', color='black')
  ax1.set_ylabel('Geostrophic imbalance', color='black')
  for tl in ax1.get_yticklabels():
    tl.set_color('black')


  print 'Plotting hydrostatic imbalance'
  for scale in range (0,len(scalelim)):
    ax2.plot(times, himbalresults[dirno,:,scale], color='0.6', linewidth='2', ls=lineone[scale], label='>' + str(scalelim[scale]) + 'km')
  ax2.legend(loc='lower right')
  ax2.set_ylabel('Hydrostatic imbalance', color='0.5')
  for tl in ax2.get_yticklabels():
    tl.set_color('0.5')

  plt.title('Geostrophic and hydrostatic imbalance', color='black')
  #plt.show()
  plt.savefig(plotfile, bbox_inches='tight')
  plt.close()

