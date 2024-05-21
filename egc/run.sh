#!/bin/bash
set -eou pipefail

# To build ~/src/radiasoft/mlopal/egc/build.sh

cd "$(dirname "$0")"
# start program with r
# get backtrace with bt
gdb --args ~/src/radiasoft/mlopal/build/src/opal ~/src/radiasoft/mlopal/egc/opal.in