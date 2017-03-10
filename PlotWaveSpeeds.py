#!/usr/bin/env python
# -------------------------------------------------------------------
# Python code to read data concerning horizontal wave freqiemcies
# and speeds and to plot them
# Ross Bannister, Dec 2016
# -------------------------------------------------------------------


# ===================================================================
def read_freq (nz=0, nlongs=0, filename=''):
  import numpy as np
  # Subroutine to read-in wave speeds
  print 'Reading: ', filename
  infile = open (filename, 'r')
  freq   = np.zeros(nlongs)

  # Read-in frequencies
  for k in range(0, nlongs, 1):
    line    = infile.readline()
    data    = line.split()
    freq[k] = data[nz]

  infile.close()
  return freq

# ===================================================================
def read_horiz_speed (nz=0, filename=''):
  import numpy as np
  # Subroutine to read-in wave speeds
  print 'Reading: ', filename
  infile = open (filename, 'r')

  # Read-in up to the relevant vertical wavenumber
  for l in range(0, nz+1, 1):
    line = infile.readline()

  speed = line.split()

  infile.close()
  return speed

# ===================================================================
def read_vert_speed (nx=0, filename=''):
  import numpy as np
  # Subroutine to read-in wave speeds
  print 'Reading: ', filename
  infile = open (filename, 'r')

  # Read-in up to the relevant horizontal wavenumber
  for l in range(0, nx+1, 1):
    line = infile.readline()

  speed = line.split()

  infile.close()
  return speed


# ===================================================================
# ===== START THE MAIN PROGRAM ======================================
# ===================================================================
import numpy as np
import matplotlib.pyplot as plt


# ===================================================================
# ===== SET-UP WORKING DIRECTORIES ==================================
# ===================================================================

# Data directory
data_dir = '/export/diamet/raid1/ross/RuthsModel/LinearAnalysis'
exp_name = 'Waves'

# Output directory for plots
plot_dir = data_dir

# Set the vertical and horizontal wavenumbers to plot
nz       = 10
nx       = 100

# Set the domain dimensions
nlongs   = 360
nlevs    = 60

print 'INPUT DIRECTORY  : ', data_dir
print 'PLOT DIRECTORY   : ', plot_dir



hori_wns = np.linspace(0,nlongs-1,nlongs)
vert_wns = np.linspace(0,nlevs-1,nlevs)

# ===================================================================
# ===== READ-IN WAVE FREQUENCIES (FROM 5 EXPERIMENTS)================
# ===================================================================

print 'Reading in gravity wave frequencies'
grav_freq01 = read_freq (nz, nlongs, data_dir + '/Run01/grav_frequency_' + exp_name + '.dat')
grav_freq02 = read_freq (nz, nlongs, data_dir + '/Run02/grav_frequency_' + exp_name + '.dat')
grav_freq03 = read_freq (nz, nlongs, data_dir + '/Run03/grav_frequency_' + exp_name + '.dat')
grav_freq04 = read_freq (nz, nlongs, data_dir + '/Run04/grav_frequency_' + exp_name + '.dat')
grav_freq05 = read_freq (nz, nlongs, data_dir + '/Run05/grav_frequency_' + exp_name + '.dat')

print 'Reading in acoustic wave frequencies'
acou_freq01 = read_freq (nz, nlongs, data_dir + '/Run01/acou_frequency_' + exp_name + '.dat')
acou_freq02 = read_freq (nz, nlongs, data_dir + '/Run02/acou_frequency_' + exp_name + '.dat')
acou_freq03 = read_freq (nz, nlongs, data_dir + '/Run03/acou_frequency_' + exp_name + '.dat')
acou_freq04 = read_freq (nz, nlongs, data_dir + '/Run04/acou_frequency_' + exp_name + '.dat')
acou_freq05 = read_freq (nz, nlongs, data_dir + '/Run05/acou_frequency_' + exp_name + '.dat')


# ===================================================================
# ===== PLOT WAVE FREQUENCIES =======================================
# ===================================================================

# Sensitivities of gravity wave frequencies to A
print 'Plotting gravity wave frequency sensitivities to A'
plot_file = plot_dir + '/GravityFreqA.eps'
plt.plot(hori_wns, grav_freq02, ls='solid',  color='black', linewidth='2', label=r'$A = 0.002$')
plt.plot(hori_wns, grav_freq01, ls='dashed', color='black', linewidth='2', label=r'$A = 0.02$')
plt.plot(hori_wns, grav_freq03, ls='dotted', color='black', linewidth='2', label=r'$A = 0.2$')
plt.xlim(xmin=0.0, xmax=360.0)
plt.legend(loc='upper left')
plt.xlabel('horizontal wavenumber index')
plt.ylabel('frequency (/s)')
plt.title('Gravity wave frequencies')
plt.savefig(plot_file)
plt.close()

# Sensitivities of acoustic wave frequencies to A
print 'Plotting acoustic wave frequency sensitivities to A'
plot_file = plot_dir + '/AcousticFreqA.eps'
plt.plot(hori_wns, acou_freq02, ls='solid',  color='black', linewidth='2', label=r'$A = 0.002$')
plt.plot(hori_wns, acou_freq01, ls='dashed', color='black', linewidth='2', label=r'$A = 0.02$')
plt.plot(hori_wns, acou_freq03, ls='dotted', color='black', linewidth='2', label=r'$A = 0.2$')
plt.xlim(xmin=0.0, xmax=360.0)
plt.legend(loc='upper left')
plt.xlabel('horizontal wavenumber index')
plt.ylabel('frequency (/s)')
plt.title('Acoustic wave frequencies')
plt.savefig(plot_file)
plt.close()

# ===================================================================

# Sensitivities of gravity wave frequencies to B
print 'Plotting gravity wave frequency sensitivities to BC'
plot_file = plot_dir + '/GravityFreqBC.eps'
plt.plot(hori_wns, grav_freq04, ls='solid',  color='black', linewidth='2', label=r'$BC = 10^1$')
plt.plot(hori_wns, grav_freq01, ls='dashed', color='black', linewidth='2', label=r'$BC = 10^2$')
plt.plot(hori_wns, grav_freq05, ls='dotted', color='black', linewidth='2', label=r'$BC = 10^3$')
plt.xlim(xmin=0.0, xmax=360.0)
plt.legend(loc='upper left')
plt.xlabel('horizontal wavenumber index')
plt.ylabel('frequency (/s)')
plt.title('Gravity wave frequencies')
plt.savefig(plot_file)
plt.close()

# Sensitivities of acoustic wave frequencies to BC
print 'Plotting acoustic wave frequency sensitivities to BC'
plot_file = plot_dir + '/AcousticFreqBC.eps'
plt.plot(hori_wns, acou_freq04, ls='solid',  color='black', linewidth='2', label=r'$BC = 10^1$')
plt.plot(hori_wns, acou_freq01, ls='dashed', color='black', linewidth='2', label=r'$BC = 10^2$')
plt.plot(hori_wns, acou_freq05, ls='dotted', color='black', linewidth='2', label=r'$BC = 10^3$')
plt.xlim(xmin=0.0, xmax=360.0)
plt.legend(loc='upper left')
plt.xlabel('horizontal wavenumber index')
plt.ylabel('frequency (/s)')
plt.title('Acoustic wave frequencies')
plt.savefig(plot_file)
plt.close()




# ===================================================================
# ===== READ-IN HORIZONTAL WAVE SPEEDS ==============================
# ===================================================================

print 'Reading in horizontal gravity wave speeds'
horiz_grav_speed01 = read_horiz_speed (nz, data_dir + '/Run01/hori_grav_speed_' + exp_name + '.dat')
horiz_grav_speed02 = read_horiz_speed (nz, data_dir + '/Run02/hori_grav_speed_' + exp_name + '.dat')
horiz_grav_speed03 = read_horiz_speed (nz, data_dir + '/Run03/hori_grav_speed_' + exp_name + '.dat')
horiz_grav_speed04 = read_horiz_speed (nz, data_dir + '/Run04/hori_grav_speed_' + exp_name + '.dat')
horiz_grav_speed05 = read_horiz_speed (nz, data_dir + '/Run05/hori_grav_speed_' + exp_name + '.dat')

print 'Reading in horizontal acoustic wave speeds'
horiz_acou_speed01 = read_horiz_speed (nz, data_dir + '/Run01/hori_acou_speed_' + exp_name + '.dat')
horiz_acou_speed02 = read_horiz_speed (nz, data_dir + '/Run02/hori_acou_speed_' + exp_name + '.dat')
horiz_acou_speed03 = read_horiz_speed (nz, data_dir + '/Run03/hori_acou_speed_' + exp_name + '.dat')
horiz_acou_speed04 = read_horiz_speed (nz, data_dir + '/Run04/hori_acou_speed_' + exp_name + '.dat')
horiz_acou_speed05 = read_horiz_speed (nz, data_dir + '/Run05/hori_acou_speed_' + exp_name + '.dat')



# ===================================================================
# ===== PLOT HORIZONTAL WAVE SPEEDS =================================
# ===================================================================


# Sensitivities of horizontal gravity wave speeds to A
print 'Plotting horizontal gravity wave frequency speeds to A'
plot_file = plot_dir + '/HorizGravitySpeedA.eps'
plt.plot(hori_wns[1:], horiz_grav_speed02[1:], ls='solid',  color='black', linewidth='2', label=r'$A = 0.002$')
plt.plot(hori_wns[1:], horiz_grav_speed01[1:], ls='dashed', color='black', linewidth='2', label=r'$A = 0.02$')
plt.plot(hori_wns[1:], horiz_grav_speed03[1:], ls='dotted', color='black', linewidth='2', label=r'$A = 0.2$')
plt.xlim(xmin=0.0, xmax=360.0)
plt.legend(loc='upper right')
plt.xlabel('horizontal wavenumber index')
plt.ylabel('speed (m/s)')
plt.title('Gravity wave speeds (horizontal)')
plt.savefig(plot_file)
plt.close()

# Sensitivities of horizontal acoustic wave speeds to A
print 'Plotting horizontal acoustic wave speed sensitivities to A'
plot_file = plot_dir + '/HorizAcousticSpeedA.eps'
plt.plot(hori_wns[1:], horiz_acou_speed02[1:], ls='solid',  color='black', linewidth='2', label=r'$A = 0.002$')
plt.plot(hori_wns[1:], horiz_acou_speed01[1:], ls='dashed', color='black', linewidth='2', label=r'$A = 0.02$')
plt.plot(hori_wns[1:], horiz_acou_speed03[1:], ls='dotted', color='black', linewidth='2', label=r'$A = 0.2$')
plt.xlim(xmin=0.0, xmax=360.0)
plt.legend(loc='upper left')
plt.xlabel('horizontal wavenumber index')
plt.ylabel('speed (m/s)')
plt.title('Acoustic wave speeds (horizontal)')
plt.savefig(plot_file)
plt.close()

# ===================================================================

# Sensitivities of horizontal gravity wave speeds to BC
print 'Plotting horizontal gravity wave speed sensitivities to BC'
plot_file = plot_dir + '/HorizGravitySpeedBC.eps'
plt.plot(hori_wns[1:], horiz_grav_speed04[1:], ls='solid',  color='black', linewidth='2', label=r'$BC = 10^1$')
plt.plot(hori_wns[1:], horiz_grav_speed01[1:], ls='dashed', color='black', linewidth='2', label=r'$BC = 10^2$')
plt.plot(hori_wns[1:], horiz_grav_speed05[1:], ls='dotted', color='black', linewidth='2', label=r'$BC = 10^3$')
plt.xlim(xmin=0.0, xmax=360.0)
plt.legend(loc='upper right')
plt.xlabel('horizontal wavenumber index')
plt.ylabel('speed (m/s)')
plt.title('Gravity wave speeds (horizontal)')
plt.savefig(plot_file)
plt.close()

# Sensitivities of horizontal acoustic wave speeds to BC
print 'Plotting horizontal acoustic wave speed sensitivities to BC'
plot_file = plot_dir + '/HorizAcousticSpeedBC.eps'
plt.plot(hori_wns[1:], horiz_acou_speed04[1:], ls='solid',  color='black', linewidth='2', label=r'$BC = 10^1$')
plt.plot(hori_wns[1:], horiz_acou_speed01[1:], ls='dashed', color='black', linewidth='2', label=r'$BC = 10^2$')
plt.plot(hori_wns[1:], horiz_acou_speed05[1:], ls='dotted', color='black', linewidth='2', label=r'$BC = 10^3$')
plt.xlim(xmin=0.0, xmax=360.0)
plt.legend(loc='upper left')
plt.xlabel('horizontal wavenumber index')
plt.ylabel('speed (m/s)')
plt.title('Acoustic wave speeds (horizontal)')
plt.savefig(plot_file)
plt.close()







# ===================================================================
# ===== READ-IN VERTICAL WAVE SPEEDS ================================
# ===================================================================

print 'Reading in vertical gravity wave speeds'
vert_grav_speed01 = read_vert_speed (nx, data_dir + '/Run01/vert_grav_speed_' + exp_name + '.dat')
vert_grav_speed02 = read_vert_speed (nx, data_dir + '/Run02/vert_grav_speed_' + exp_name + '.dat')
vert_grav_speed03 = read_vert_speed (nx, data_dir + '/Run03/vert_grav_speed_' + exp_name + '.dat')
vert_grav_speed04 = read_vert_speed (nx, data_dir + '/Run04/vert_grav_speed_' + exp_name + '.dat')
vert_grav_speed05 = read_vert_speed (nx, data_dir + '/Run05/vert_grav_speed_' + exp_name + '.dat')

print 'Reading in vertical acoustic wave speeds'
vert_acou_speed01 = read_vert_speed (nx, data_dir + '/Run01/vert_acou_speed_' + exp_name + '.dat')
vert_acou_speed02 = read_vert_speed (nx, data_dir + '/Run02/vert_acou_speed_' + exp_name + '.dat')
vert_acou_speed03 = read_vert_speed (nx, data_dir + '/Run03/vert_acou_speed_' + exp_name + '.dat')
vert_acou_speed04 = read_vert_speed (nx, data_dir + '/Run04/vert_acou_speed_' + exp_name + '.dat')
vert_acou_speed05 = read_vert_speed (nx, data_dir + '/Run05/vert_acou_speed_' + exp_name + '.dat')



# ===================================================================
# ===== PLOT VERTICAL WAVE SPEEDS ===================================
# ===================================================================


# Sensitivities of vertical gravity wave speeds to A
print 'Plotting vertical gravity wave frequency speeds to A'
plot_file = plot_dir + '/VertGravitySpeedA.eps'
plt.plot(vert_wns[1:], vert_grav_speed02[1:], ls='solid',  color='black', linewidth='2', label=r'$A = 0.002$')
plt.plot(vert_wns[1:], vert_grav_speed01[1:], ls='dashed', color='black', linewidth='2', label=r'$A = 0.02$')
plt.plot(vert_wns[1:], vert_grav_speed03[1:], ls='dotted', color='black', linewidth='2', label=r'$A = 0.2$')
plt.xlim(xmin=0.0, xmax=60.0)
plt.legend(loc='upper right')
plt.xlabel('vertical wavenumber index')
plt.ylabel('speed (m/s)')
plt.title('Gravity wave speeds (vertical)')
plt.savefig(plot_file)
plt.close()

# Sensitivities of vertical acoustic wave speeds to A
print 'Plotting vertical acoustic wave speed sensitivities to A'
plot_file = plot_dir + '/VertAcousticSpeedA.eps'
plt.plot(vert_wns[1:], vert_acou_speed02[1:], ls='solid',  color='black', linewidth='2', label=r'$A = 0.002$')
plt.plot(vert_wns[1:], vert_acou_speed01[1:], ls='dashed', color='black', linewidth='2', label=r'$A = 0.02$')
plt.plot(vert_wns[1:], vert_acou_speed03[1:], ls='dotted', color='black', linewidth='2', label=r'$A = 0.2$')
plt.xlim(xmin=0.0, xmax=60.0)
plt.legend(loc='lower right')
plt.xlabel('vertical wavenumber index')
plt.ylabel('speed (m/s)')
plt.title('Acoustic wave speeds (vertical)')
plt.savefig(plot_file)
plt.close()

# ===================================================================

# Sensitivities of vertical gravity wave speeds to BC
print 'Plotting vertical gravity wave speed sensitivities to BC'
plot_file = plot_dir + '/VertGravitySpeedBC.eps'
plt.plot(vert_wns[1:], vert_grav_speed04[1:], ls='solid',  color='black', linewidth='2', label=r'$BC = 10^1$')
plt.plot(vert_wns[1:], vert_grav_speed01[1:], ls='dashed', color='black', linewidth='2', label=r'$BC = 10^2$')
plt.plot(vert_wns[1:], vert_grav_speed05[1:], ls='dotted', color='black', linewidth='2', label=r'$BC = 10^3$')
plt.xlim(xmin=0.0, xmax=60.0)
plt.legend(loc='lower right')
plt.xlabel('vertical wavenumber index')
plt.ylabel('speed (m/s)')
plt.title('Gravity wave speeds (vertical)')
plt.savefig(plot_file)
plt.close()

# Sensitivities of vertical acoustic wave speeds to BC
print 'Plotting vertical acoustic wave speed sensitivities to BC'
plot_file = plot_dir + '/VertAcousticSpeedBC.eps'
plt.plot(vert_wns[1:], vert_acou_speed04[1:], ls='solid',  color='black', linewidth='2', label=r'$BC = 10^1$')
plt.plot(vert_wns[1:], vert_acou_speed01[1:], ls='dashed', color='black', linewidth='2', label=r'$BC = 10^2$')
plt.plot(vert_wns[1:], vert_acou_speed05[1:], ls='dotted', color='black', linewidth='2', label=r'$BC = 10^3$')
plt.xlim(xmin=0.0, xmax=60.0)
plt.legend(loc='center right')
plt.xlabel('vertical wavenumber index')
plt.ylabel('speed (m/s)')
plt.title('Acoustic wave speeds (vertical)')
plt.savefig(plot_file)
plt.close()

