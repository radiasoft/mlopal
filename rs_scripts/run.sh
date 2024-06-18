#!/bin/bash
set -eou pipefail

# To build ~/src/radiasoft/mlopal/egc/build.sh

cd "$(dirname "$0")"../
# start program with r
# get backtrace with bt
gdb --args $PWD/build/src/opal $PWD/rs_scripts/opal.in
