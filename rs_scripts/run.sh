#!/bin/bash
set -eou pipefail

# To build ~/src/radiasoft/mlopal/egc/build.sh

cd $(git rev-parse --show-toplevel)
# start program with r
# get backtrace with bt
gdb --args $PWD/build/src/opal $PWD/rs_scripts/opal.in
