#!/bin/bash


start_tm(){
    # Get current time
    T=$(date +%s)
}

end_tm(){
    # Get ending time ddhhmm
    echo -e "\n\n========"
    echo "  DONE"
    echo -e "========\n\n"
    local R D H M S
    R=$(($(date +%s) - T))
    D=$((R/60/60/24))
    H=$((R/60/60%24))
    M=$((R/60%60))
    S=$((R%60))

    printf 'CMD-> %s\n' "$0"
    printf 'RUNTIME-> '
    (( D > 0 )) && printf '%d d ' ${D}
    (( H > 0 )) && printf '%d h ' ${H}
    (( M > 0 )) && printf '%d m ' ${M}
    (( D > 0 || H > 0 || M > 0 )) && printf 'and '
    printf '%d s\n' ${S}
}

tm() {
    # Get starting time ddhhmm
    local T command rc R D H M S
    T=$(date +%s)
    command=("$@")
    rc=$?
    R=$(($(date +%s) - T))
    D=$((R/60/60/24))
    H=$((R/60/60%24))
    M=$((R/60%60))
    S=$((R%60))

    printf 'CMD-> %s\n' "${command[@]}"
    printf 'RUNTIME-> '
    (( D > 0 )) && printf '%d d ' ${D}
    (( H > 0 )) && printf '%d h ' ${H}
    (( M > 0 )) && printf '%d m ' ${M}
    (( D > 0 || H > 0 || M > 0 )) && printf 'and '
    printf '%d s\n' ${S}
    return ${rc}
}

error(){
    # create error log file
    local parent_lineno script message code
    parent_lineno="$1"
    script="$2"
    message="$3"
    code="${4:-1}"

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