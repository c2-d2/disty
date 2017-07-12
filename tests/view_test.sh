#! /usr/bin/env bash

set -u
set -e
set -o pipefail

readonly PROGNAME=$(basename $0)
readonly PROGDIR=$(dirname $0)
readonly ARGS="$@"
readonly NARGS="$#"

if [ $NARGS -ne 1 ]; then
	echo "usage: $PROGNAME options"
	exit 1
fi

head {test$1*,exp.test$1*}


