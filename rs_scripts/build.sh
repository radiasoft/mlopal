#!/bin/bash
set -eou pipefail
#
# Build / recompile opal
#


# TODO(e-carlin): copied from codes.sh
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
    codes_dir[pyenv_prefix]=$(realpath "$(pyenv prefix)")
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

declare repo_root=$(git rev-parse --show-toplevel)
cd "$repo_root"
if [[ ! -d build ]]; then
    cd /
    # TODO(e-carlin): add back in
    # rpm2cpio /home/vagrant/jupyter/StaffScratch/e-carlin/src/radiasoft/chris_mlopal/rs_scripts/rscode-trilinos-20230829.182134-1.x86_64.rpm | cpio -ivd &> /dev/null
    cd "$repo_root"
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
cd "$repo_root"/build
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
