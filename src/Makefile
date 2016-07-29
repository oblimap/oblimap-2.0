# File name: Makefile
#
# Copyright (C) 2016 Thomas Reerink.
#
# This file is distributed under the terms of the
# GNU General Public License.
#
# This file is part of OBLIMAP 2.0
#
# OBLIMAP's scientific documentation and its first open source
# release (see the supplement) is published at:
# http://www.geosci-model-dev.net/3/13/2010/gmd-3-13-2010.html
#
# OBLIMAP is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# OBLIMAP is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with OBLIMAP. If not, see <http://www.gnu.org/licenses/>.
#
#
# OBLIMAP is maintained by:
#
# Thomas Reerink
# Institute for Marine and Atmospheric Research Utrecht (IMAU)
# Utrecht University
# Princetonplein 5
# 3584 CC Utrecht
# The Netherlands
#
# email: tjreerink@gmail.com
#


#    Makefile for compiling OBLIMAP
#
#    Execute options of this Makefile (in the src/ directory) e.g.:
#     make clean
#     make all
#     make oblimap_gcm_to_im_program
#     make oblimap_im_to_gcm_program
#     make oblimap_convert_program

#    For "make" syntax search on "make tutorial", or see:
#     http://www.gnu.org/software/make/manual/make.html

# Include below the Makefile.include which matches with your fortran compiler:
include Makefile.gfortran
#include Makefile.gfortran-mp-4.9-mac
#include Makefile.gfortran-lisa
#include Makefile.ifort-lisa
#include Makefile.ifort-staff.science
#include Makefile.ftn-cca-ecmwf
#include Makefile.ifort-stampede

# The OBJ_PATH and the MOD_PATH are defined:
OBJ_PATH = object-files
MOD_PATH = module-files

# Pattern rules for creating object files from *.f90 (fortran 90) source files
$(OBJ_PATH)/%.o : %.f90
	$(F90) $(F90FLAGS) -c -o $@ $<

# This rule is used to build an executable from the object files:
% :
	$(F90) $(F90FLAGS) -o $@ $^ $(LDFLAGS)



# =============================================================
#       Building the oblimap_gcm_to_im_program
# =============================================================

# Adding this program to the list of all programs:
all_programs += oblimap_gcm_to_im_program

# List of object files for oblimap_gcm_to_im_program:
oblimap_gcm_to_im_program_objects = \
        oblimap_configuration_module.o \
        oblimap_read_and_write_module.o \
        oblimap_projection_module.o \
        oblimap_scan_contributions_module.o \
        oblimap_mapping_module.o \
        oblimap_gcm_to_im_mapping_module.o \
        oblimap_gcm_to_im_program.o

# The executable (the target) depends on the object files:
oblimap_gcm_to_im_program : $(oblimap_gcm_to_im_program_objects:%=$(OBJ_PATH)/%)



# =============================================================
#       Building the oblimap_im_to_gcm_program
# =============================================================

# Adding this program to the list of all programs:
all_programs += oblimap_im_to_gcm_program

# List of object files for oblimap_im_to_gcm_program:
oblimap_im_to_gcm_program_objects = \
        oblimap_configuration_module.o \
        oblimap_read_and_write_module.o \
        oblimap_projection_module.o \
        oblimap_scan_contributions_module.o \
        oblimap_mapping_module.o \
        oblimap_im_to_gcm_mapping_module.o \
        oblimap_im_to_gcm_program.o

# The executable (the target) depends on the object files:
oblimap_im_to_gcm_program : $(oblimap_im_to_gcm_program_objects:%=$(OBJ_PATH)/%)



# =============================================================
#       Building the oblimap_convert_program
# =============================================================

# Adding this program to the list of all programs:
all_programs += oblimap_convert_program

# List of object files for oblimap_convert_program:
oblimap_convert_program_objects = \
        oblimap_configuration_module.o \
        oblimap_read_and_write_module.o \
        oblimap_projection_module.o \
        oblimap_scan_contributions_module.o \
        oblimap_convert_module.o \
        oblimap_convert_program.o

# The executable (the target) depends on the object files:
oblimap_convert_program : $(oblimap_convert_program_objects:%=$(OBJ_PATH)/%)



# =============================================================
#       Cleaning up everything
# =============================================================

.PHONY : clean

clean :
	rm -fv $(all_programs) $(OBJ_PATH)/*.o $(MOD_PATH)/*.mod *.mod



# =============================================================
#       Building everything
# =============================================================

.PHONY : all

all : $(all_programs)
