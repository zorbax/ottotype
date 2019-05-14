#!/bin/bash

tm() {

  local T=$(date +%s)
  local command="$@"
  $@
  rc=$?
  local R=$[$(date +%s)-$T]
  local D=$((R/60/60/24))
  local H=$((R/60/60%24))
  local M=$((R/60%60))
  local S=$((R%60))

  printf 'CMD-> %s\n' "$command"
  printf 'RUNTIME-> '
  (( $D > 0 )) && printf '%d d ' $D
  (( $H > 0 )) && printf '%d h ' $H
  (( $M > 0 )) && printf '%d m ' $M
  (( $D > 0 || $H > 0 || $M > 0 )) && printf 'and '
  printf '%d s\n' $S
  return $rc
}

error(){
  local parent_lineno="$1"
  local script="$2"
  local message="$3"
  local code="${4:-1}"

  if [[ -n "$message" ]] ; then
    echo -e "\n---------------------------------------\n"
    echo -e "ERROR in $script on or near line ${parent_lineno}; exiting with status ${code}"
    echo -e "MESSAGE:\n"
    echo -e "$message"
    echo -e "\n---------------------------------------\n"
  else
    echo -e "\n---------------------------------------\n"
    echo -e "ERROR in $script on or near line ${parent_lineno}; exiting with status ${code}"
    echo -e "\n---------------------------------------\n"
  fi

  exit "${code}"
}
