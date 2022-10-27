#!/usr/bin/env bash
function helptext() {
  echo "
  Usage: ./helper.sh args 

  Args Required -
  -----------------------------------------------------------------------------
  build_venv    : Generate python venv with all the dependencies.
  -----------------------------------------------------------------------------
  "
}

REQUIREMENTS_FILE="requirements.txt"
CODE_DIR="src"

function build_venv() {
    python3 -m venv venv
    source venv/bin/activate
    python -m pip install -r $REQUIREMENTS_FILE
}

if [ -z "$1" ]; then
    helptext
elif [[ $(declare -F) == *"$1"* ]]; then
  set -x
  $1 ${@:2}
fi