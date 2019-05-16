#!/bin/bash

end(){

  echo -e "\n\n========"
  echo "  DONE"
  echo -e "========\n\n"
  R=$(($(date +%s)-$T))
  D=$((R/60/60/24))
  H=$((R/60/60%24))
  M=$((R/60%60))
  S=$((R%60))

  printf 'CMD-> %s\n' "$0"
  printf 'RUNTIME-> '
  (( $D > 0 )) && printf '%d d ' $D
  (( $H > 0 )) && printf '%d h ' $H
  (( $M > 0 )) && printf '%d m ' $M
  (( $D > 0 || $H > 0 || $M > 0 )) && printf 'and '
  printf '%d s\n' $S
}

tm() {

  local T=$(date +%s)
  local command="$@"
  $@
  rc=$?
  local R=$(($(date +%s)-$T))
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

mkdir -p log
log_file=log/ottotype.log

if [ -f $log_file ];then
  rm -f $log_file
fi

cat << EOF > $log_file
====================================
Serotyping pipeline from SSB-CNRDOGM
====================================
EOF

echo -e "\nLOG FILE OTTOTYPE" >> $log_file
echo $(date +"%Y-%m-%d %H:%M") >> $log_file
