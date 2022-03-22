project (orbit-integrator)

cmake_minimum_required(VERSION 3.19)
enable_testing()

include_directories(src/)

# Flags needed for GSL and HDF5
set (LINK_FLAGS "-lgsl -lgslcblas -lm -lhdf5")
set (CMAKE_EXE_LINKER_FLAGS ${LINK_FLAGS})

configure_file(test/input_keplerian_rotation.txt input_keplerian_rotation.txt COPYONLY)
configure_file(test/input_keplerian_rotation_2.txt input_keplerian_rotation_2.txt COPYONLY)
file(COPY extras/stellar_profiles DESTINATION .)
file(COPY test/interp_data.h5 DESTINATION .)


add_executable(kep test/keplerian_rotation.c src/io.c src/model.c src/units.c src/star.c src/utils.c)
add_test(NAME keplerian_rotation COMMAND kep input_keplerian_rotation.txt)

add_executable(kep_2 test/keplerian_rotation_2.c src/io.c src/model.c src/units.c src/star.c src/utils.c)
add_test(NAME keplerian_rotation_2 COMMAND kep_2 input_keplerian_rotation_2.txt)

add_executable(interp test/interpolation.c src/io.c src/utils.c src/star.c src/model.c)
add_test(NAME interpolation COMMAND interp)
