PREFIX  := $(PWD)
PROFILE := $(HOME)/.bashrc
PROJECT := "setTestProject"
NUMCPUS := $(nproc)
SHELL   := /bin/bash

# Derived variables
TMPDIR := build
TMPTARFILE=$(TMPDIR)/$(TARFILE)

# Style variables
T= "	"
T2=$(T)$(T)

default: install

help:
  @echo "Please see README.md for additional help"

all: install env

install: install-prerequisites
  @echo "Don't forget to set up update PATH to $(PWD)/scripts";
  @echo "'make env' performs this step for you"
  @echo "DONE: installation of Ottotype complete."

install-prerequisites: install-mkdir install-dockers
  @echo DONE installing prerequisites

clean: clean-tmp clean-symlinks clean-SOMETHING
  @echo "Remember to remove the line with PATH and Lyve-SET from $(PROFILE)"

install-mkdir:
  -mkdir build lib scripts

clean-tmp:
  rm -rfv $(TMPDIR)/*
