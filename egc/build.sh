#!/bin/bash
set -eou pipefail
#
# Build / recompile opal
#
# Need new pybind `pip install --upgrade pybind11[global]`
# pybind11-config --cmakedir

declare -A codes_dir=()

codes_dir_setup() {
    declare todo=()
    declare d n p
    for n in prefix \
        etc \
        bashrc_d \
        bin \
        lib \
        include \
        share
    do
        case $n in
            bashrc_d)
                p=/etc/bashrc.d
                ;;
            prefix)
                p=
                ;;
            *)
                p=/$n
                ;;
        esac
        d=$HOME/.local$p
        todo+=( $d )
        codes_dir[$n]=$d
    done
    # if codes_is_common; then
    #     install_msg 'creating directories'
    #     mkdir -p "${todo[@]}"
    # else

    #     # common doesn't use pyenv_prefix, it creates it
        codes_dir[pyenv_prefix]=$(realpath "$(pyenv prefix)")
    # fi
}

codes_cmake() {
    mkdir build
    cd build
    declare t=Release
    if [[ ${CODES_DEBUG_FLAG:-} ]]; then
        t=Debug
    fi
    CLICOLOR=0 cmake -D CMAKE_RULE_MESSAGES:BOOL=OFF -D CMAKE_BUILD_TYPE:STRING="$t" "$@" ..
}

codes_make() {
    declare cmd=( make -j1 )
    if [[ $@ ]]; then
        cmd+=( "$@" )
    fi
    "${cmd[@]}"
}
codes_dir_setup

cd ~/src/radiasoft/mlopal
if [[ ! -d build ]]; then
    perl -pi -e 's{add_compile_options \(-Werror\)}{}' CMakeLists.txt
    perl -pi -e '
        # https://stackoverflow.com/a/20991533
        # boost is compiled multithreaded, because it does not mean "pthreads",
        # but just that the code takes a bit more care on values in static sections.
        # If we do not turn this ON, it will not find the variant compiled.
        s{(?<=Boost_USE_MULTITHREADED )OFF}{ON};
        # otherwise fails with -lmpi_mpifh not found, because
        # that is part of openmpi, not mpich
        s{.*mpi_mpifh.*}{};
        s{-fPIE}{};
        s{add_link_options.*-pie.*}{};
    ' CMakeLists.txt
    # need to specify CC and CXX otherwise build uses wrong
    # compiler.
    H5HUT_PREFIX="${codes_dir[prefix]}" \
        BOOST_DIR="${codes_dir[prefix]}" \
        HDF5_INCLUDE_DIR=/usr/include \
        HDF5_LIBRARY_DIR="$BIVIO_MPI_LIB" \
        MITHRA_INCLUDE_DIR="${codes_dir[include]}" \
        MITHRA_LIBRARY_DIR="${codes_dir[lib]}" \
        CC=mpicc CXX=mpicxx \
        codes_cmake \
        -D CMAKE_INSTALL_PREFIX="${codes_dir[prefix]}" \
        -D CMAKE_POSITION_INDEPENDENT_CODE=FALSE \
        -D ENABLE_OPAL_FEL=yes \
        -D ENABLE_SAAMG_SOLVER=TRUE \
        -D USE_STATIC_LIBRARIES=FALSE
fi
cd ~/src/radiasoft/mlopal/build
codes_make all


    # H5HUT_PREFIX="${codes_dir[prefix]}" \
    #     BOOST_DIR="${codes_dir[prefix]}" \
    #     HDF5_INCLUDE_DIR=/usr/include \
    #     HDF5_LIBRARY_DIR="$BIVIO_MPI_LIB" \
    #     MITHRA_INCLUDE_DIR="${codes_dir[include]}" \
    #     MITHRA_LIBRARY_DIR="${codes_dir[lib]}" \
    #     CC=mpicc CXX=mpicxx \
    # CLICOLOR=0 cmake -D CMAKE_RULE_MESSAGES:BOOL=OFF -D CMAKE_BUILD_TYPE:STRING="Debug" --debug-output .. \
    #     -D CMAKE_INSTALL_PREFIX="${codes_dir[prefix]}" \
    #     -D CMAKE_POSITION_INDEPENDENT_CODE=FALSE \
    #     -D ENABLE_OPAL_FEL=yes \
    #     -D ENABLE_SAAMG_SOLVER=TRUE \
    #     -D USE_STATIC_LIBRARIES=FALSE