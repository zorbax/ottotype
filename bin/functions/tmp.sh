#!/bin/bash

tmpmk(){
  if [[ -O $PWD/TMP && -d $PWD/TMP ]]; then
    TMPDIR=$PWD/TMP
  else
    rm -rf $PWD/TMP 2> /dev/null
    mkdir -p $PWD/TMP
    TMPDIR=$(mktemp -d $PWD/TMP/XXXX)
  fi

  TMP=$TMPDIR
  TEMP=$TMPDIR
  export TMPDIR TMP TEMP
}

tmprm(){
  TMPDIR=$PWD/TMP

  if [[ -O $TMPDIR && -d $TMPDIR ]]; then
    rm -rf $TMPDIR/*
  fi
}
