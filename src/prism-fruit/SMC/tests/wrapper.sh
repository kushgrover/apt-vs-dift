#!/bin/sh

PRISM_DIR="$(dirname $(dirname $(realpath $0)))"
timeout -s KILL -k 1s 5m "$PRISM_DIR/prism/bin/prism" $@
