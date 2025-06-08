#!/usr/bin/env bash 

set -euo pipefail # exit on error, unset vars, pipe failures 

detect_arch() { # this here is a funciton that detects the architechture of your machine 
  uname -m 
}

have_conda() { command -v conda >/dev/null 2>&1; } # return true if conda exists 

get_conda_subdir() { # asking existing conda what platform its solving for 
  conda config --show subdir 2>/dev/null | awk '/^subdir/ {print $2}' 
}

patch_needed=false # default: assume no fix required 
if have_conda; then # inspect if conda exists 
   curr=$(get_conda_subdir || true) # may be empty or unset here 
   [[ -z $curr ]] && curr=osx-64 # empty means conda uses its default 
   if [[ $curr != osx-64 ]]; then # when key  is missing treat as intel osx-64 
     echo "[installer] conda subdir is '$curr' - will patch to osx-64" 
     patch_needed=true 
   fi
fi  # fi closes both if statements here 

