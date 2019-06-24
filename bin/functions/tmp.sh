#!/bin/bash

tmpmk(){
  if [[ -O $PWD/tmp && -d $PWD/tmp ]]; then
    TMPDIR=$PWD/tmp
  else
    rm -rf $PWD/tmp 2> /dev/null
    mkdir -p $PWD/tmp
    TMPDIR=$(mktemp -d $PWD/tmp/XXXX)
  fi

  TMP=$TMPDIR
  TEMP=$TMPDIR
  export TMPDIR TMP TEMP
}

tmprm(){
  if [[ -O $TMPDIR && -d $TMPDIR ]]; then
    rm -rf $TMPDIR/*
  fi
}
