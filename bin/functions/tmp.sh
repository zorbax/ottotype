#!/bin/bash

tmpmk(){
    # Create temporary directory
    if [[ -O $PWD/tmp && -d $PWD/tmp ]]; then
        TMPDIR=$PWD/tmp
    else
        rm -rf $PWD/tmp 2> /dev/null
        mkdir -p $PWD/tmp
        TMPDIR=$(mktemp -d $PWD/tmp/XXXX)
    fi

    local TMP=$TMPDIR
    local TEMP=$TMPDIR
    export TMPDIR TMP TEMP
}

tmprm(){
    # Remove temporary directory
    TMPDIR=$PWD/tmp

    if [[ -O $TMPDIR && -d $TMPDIR ]]; then
        rm -rf $TMPDIR
    fi
}
