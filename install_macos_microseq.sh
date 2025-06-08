#!/usr/bin/env bash 

set -euo pipefail # exit on error, unset vars, pipe failures 

detect_arch() { # this here is a funciton that detects the architechture of your machine 
  uname -m 
} 
