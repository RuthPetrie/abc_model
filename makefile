# On Meteorology departmental linux pre commands:
# setup studio12
# export NAG_KUSARI_FILE=/opt/compilers/NAG/license.dat


# These compile and link options are for carrot (comment out for ubuntu and dept linux)
# The -O3 is for faster code
#CompileOpts=-O3 -dalign
#CompileOpts_NC=-dalign
#CompileOpts_NAG=-dalign -M/opt/tools/lib/nag_mod_dir
#LinkOpts=-lnetcdf -lhdf5_hl -lhdf5 -lnagfl90

# These compile and link options are for ubuntu (comment out for carrot and dept linux)
# The -O3 is for faster code
#CompileOpts=-O3
#CompileOpts_NC=
#LinkOpts=-I/usr/include/ -L/usr/lib/ -lnetcdff -lnetcdf

# These compile and link options are for departmental linux (comment out for carrot and ubuntu)
# The -O3 is for faster code
CompileOpts=-O3
CompileOpts_NC=-O3
CompileOpts_NAG=-M/opt/compilers/NAG/fnl6a04ddl/nagfl90_modules
LinkOpts=-I/opt/graphics/64/include -M/opt/graphics/64/include -L/opt/graphics/64/lib -lnetcdf-4.0 -lhdf5_hl -lhdf5 -L/opt/compilers/NAG/fnl6a04ddl/lib -R/opt/compilers/NAG/fnl6a04ddl/lib -lnagfl90_spl -xlic_lib=sunperf



# These two lines are for carrot and departmental linux (comment out for ubuntu)
Main.out: Main.o DefConsTypes.o UM_data_proc.o Functions.o Diagnostics.o Boundaries.o BoundaryMod.o Forcings.o Initialise.o Model.o ReadWrite_data.o Linear_Analysis.o
	f95 -o Main.out Main.o DefConsTypes.o UM_data_proc.o Functions.o Diagnostics.o Boundaries.o BoundaryMod.o Forcings.o Initialise.o Model.o ReadWrite_data.o Linear_Analysis.o $(LinkOpts)

# These two lines are for ubuntu (comment out for carrot and departmental linux)
#Main.out: Main.o DefConsTypes.o UM_data_proc.o Functions.o Diagnostics.o Boundaries.o BoundaryMod.o Forcings.o Initialise.o Model.o ReadWrite_data.o
#	f95 -o Main.out Main.o DefConsTypes.o UM_data_proc.o Functions.o Diagnostics.o Boundaries.o BoundaryMod.o Forcings.o Initialise.o Model.o ReadWrite_data.o $(LinkOpts)

Main.o: Main.f90 DefConsTypes.o
	f95 -c Main.f90 $(CompileOpts)

DefConsTypes.o: DefConsTypes.f90
	f95 -c DefConsTypes.f90 $(CompileOpts)

UM_data_proc.o: UM_data_proc.f90 DefConsTypes.o
	f95 -c UM_data_proc.f90 $(CompileOpts_NC)

Functions.o: Functions.f90 DefConsTypes.o
	f95 -c Functions.f90 $(CompileOpts)

Diagnostics.o: Diagnostics.f90 DefConsTypes.o
	f95 -c Diagnostics.f90 $(CompileOpts_NAG)

Boundaries.o: Boundaries.f90 DefConsTypes.o
	f95 -c Boundaries.f90 $(CompileOpts)

BoundaryMod.o: BoundaryMod.f90 DefConsTypes.o
	f95 -c BoundaryMod.f90 $(CompileOpts)

Forcings.o: Forcings.f90 DefConsTypes.o
	f95 -c Forcings.f90 $(CompileOpts)

Initialise.o: Initialise.f90 DefConsTypes.o
	f95 -c Initialise.f90 $(CompileOpts)

Model.o: Model.f90 DefConsTypes.o
	f95 -c Model.f90 $(CompileOpts)

ReadWrite_data.o: ReadWrite_data.f90 DefConsTypes.o
	f95 -c ReadWrite_data.f90 $(CompileOpts_NC)

# These two lines are for carrot and departmental linux (comment out for ubuntu)
Linear_Analysis.o: Linear_Analysis.f90 DefConsTypes.o
	f95 -c Linear_Analysis.f90 $(CompileOpts_NAG)

clean:
	rm -rf *.o *.out *.mod
