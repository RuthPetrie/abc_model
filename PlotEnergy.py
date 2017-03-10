#!/usr/bin/env python
# -------------------------------------------------------------------
# Python code to read and plot total energy to show loss of conservation
#
# This is repeated for a specified range of data directories
#
# Ross Bannister, Dec 2016
# -------------------------------------------------------------------

# ===================================================================
import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset

# Set directories (each is a different experiment)
data_dir = []
data_dir.append('/export/diamet/raid1/ross/RuthsModel/ReferenceRun_dt_div_10')
data_dir.append('/export/diamet/raid1/ross/RuthsModel/ReducedA_dt_div_10')
data_dir.append('/export/diamet/raid1/ross/RuthsModel/ReducedB_dt_div_10')
data_dir.append('/export/diamet/raid1/ross/RuthsModel/ReducedC_dt_div_10')
data_dir.append('/export/diamet/raid1/ross/RuthsModel/IncreasedA_dt_div_10')
data_dir.append('/export/diamet/raid1/ross/RuthsModel/IncreasedB_dt_div_10')
data_dir.append('/export/diamet/raid1/ross/RuthsModel/IncreasedC_dt_div_10')

# Set legends for each experiment
legend = []           ;  thisls = []              ;  thiscol = []
legend.append('Ref')  ;  thisls.append('solid')   ;  thiscol.append('k')
legend.append('A-')   ;  thisls.append('dashed')  ;  thiscol.append('0.75')
legend.append('B-')   ;  thisls.append('dotted')  ;  thiscol.append('0.75')
legend.append('C-')   ;  thisls.append('dashdot') ;  thiscol.append('0.75')
legend.append('A+')   ;  thisls.append('dashed')  ;  thiscol.append('k')
legend.append('B+')   ;  thisls.append('dotted')  ;  thiscol.append('k')
legend.append('C+')   ;  thisls.append('dashdot') ;  thiscol.append('k')




datafile = 'Diagnostics.dat'

plotdir  = '/export/diamet/raid1/ross/RuthsModel'

# The last time considered (seconds)
tstep    = 10800.0

# Set-up the plot
fig, ax = plt.subplots()

counter = -1

for datadir in data_dir:

  counter        += 1
  input_file_name = datadir + '/' + datafile
  input_file      = open (input_file_name, 'r')

  print 'Dealing with ', input_file_name

  energy = []
  t      = 0.0
  time   = []

  while t < tstep:
    line     = input_file.readline()
    data     = line.split()
    e        = float(data[5])
    if (t == 0.0):
      e0 = e
    t        = float(data[1])
    energy.append(e/e0)
    time.append(t)

  input_file.close()

  # Add this line to the plot
  cax = ax.plot(time, energy, label=legend[counter], ls=thisls[counter], color=thiscol[counter])


# Labels etc
ax.set_title('Relative total energy')
ax.set_xlabel('time (s)')
ax.set_ylabel('Total energy (E/E0)')
ax.legend(loc='lower left')
#plt.show()
plt.savefig(plotdir + '/Energy_dt_div_10.eps')
plt.close()                                             
