# Copyright (c) 2017, Institut National de Recherche en Informatique et en Automatique (Inria)
#               2017, EDF (https://www.edf.fr/)
#               2020, Institut National de Recherche en Informatique et en Automatique (Inria)
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# * Redistributions of source code must retain the above copyright notice, this
#   list of conditions and the following disclaimer.
#
# * Redistributions in binary form must reproduce the above copyright notice,
#   this list of conditions and the following disclaimer in the documentation
#   and/or other materials provided with the distribution.
#
# * Neither the name of the copyright holder nor the names of its
#   contributors may be used to endorse or promote products derived from
#   this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
# FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

cmake_minimum_required(VERSION 3.7.2)
project(melissa-heat-pde-example LANGUAGES C Fortran)

include(GNUInstallDirs)

find_package(Melissa CONFIG REQUIRED)
find_package(MPI REQUIRED)
find_package(PythonInterp 3.5.3 REQUIRED)
# show if paths are accurate
message(STATUS "Melissa_include include dir: ${Melissa_INCLUDE_DIR}")
message(STATUS "Melissa_library include dir: ${Melissa_LIBRARY}")

# fixes compilation on CentOS 8 with CMake 3.11.4, OpenMPI 4.0.2
separate_arguments(MPI_C_COMPILE_FLAGS)
separate_arguments(MPI_Fortran_COMPILE_FLAGS)

# remove leading, trailing whitespace when
# * CMake 3.7 and MPICH are used, or
# * on Devuan Ascii with Spack and OpenMPI.
if(MPI_C_LINK_FLAGS)
    string(STRIP ${MPI_C_LINK_FLAGS} MPI_C_LINK_FLAGS)
endif()
if(MPI_Fortran_LINK_FLAGS)
    string(STRIP ${MPI_Fortran_LINK_FLAGS} MPI_Fortran_LINK_FLAGS)
endif()


# this is a convenience operation for users following the README instructions
#configure_file(options.py options.py COPYONLY)
#configure_file(plot-results.py plot-results.py COPYONLY)


add_library(heat_utils STATIC heat_utils.f90)
target_compile_options(heat_utils PRIVATE ${MPI_Fortran_COMPILE_FLAGS})
target_include_directories(heat_utils PRIVATE ${MPI_Fortran_INCLUDE_PATH} ${Melissa_INCLUDE_DIR})
target_link_libraries(
    heat_utils PRIVATE ${MPI_Fortran_LINK_FLAGS} ${MPI_Fortran_LIBRARIES}
)

set(CMAKE_C_EXTENSIONS OFF)
set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED ON)
add_executable(heatc heat.c)
target_compile_options(heatc PRIVATE ${MPI_C_COMPILE_FLAGS})
target_include_directories(heatc PRIVATE ${MPI_C_INCLUDE_PATH} ${Melissa_INCLUDE_DIR})
target_link_libraries(
	heatc heat_utils ${Melissa_LIBRARY} ${MPI_C_LINK_FLAGS} ${MPI_C_LIBRARIES}
)

add_executable(heat_no_melissac heat_no_melissa.c)
target_compile_options(heat_no_melissac PRIVATE ${MPI_C_COMPILE_FLAGS})
target_include_directories(heat_no_melissac PRIVATE ${MPI_C_INCLUDE_PATH} ${Melissa_INCLUDE_DIR})
target_link_libraries(
	heat_no_melissac heat_utils ${Melissa_LIBRARY} ${MPI_C_LINK_FLAGS} ${MPI_C_LIBRARIES}
)

add_executable(heatf heat.f90)
target_compile_options(heatf PRIVATE ${MPI_Fortran_COMPILE_FLAGS})
target_include_directories(heatf PRIVATE ${MPI_Fortran_INCLUDE_PATH} ${Melissa_INCLUDE_DIR})
target_link_libraries(
	heatf heat_utils ${Melissa_LIBRARY} ${MPI_Fortran_LINK_FLAGS} ${MPI_Fortran_LIBRARIES}
)

install(
    TARGETS heatc heatf
    DESTINATION ${CMAKE_INSTALL_BINDIR}
)

