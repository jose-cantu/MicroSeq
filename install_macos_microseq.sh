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

# --- short-circuit reuse existing Intel conda -----------
if { have_conda && [[ $(get_conda_subdir || echo osx-64) = osx-64 ]]; } \
   || [[ -d "$HOME/miniconda3" || -d "$HOME/anaconda3" ]]; then
  echo "[installer] Existing osx-64 conda detected so skipping bootstrp." 
  return 0 # script source'd so return is safer vs exit .....
fi # fall thorugh when conda missing or wrong-arch


# Detect CI / non-interacitve shells and pre-seed choice here 
if [[ ! -t 0 ]]; then # stdin is not a TTY -> scipt source by CI 
  choice=${MICROSEQ_CHOICE:-1} # allow override via env; default to 1 (miniconda)
else 
  choice="" # interactive session -> will use prompts I made below 
fi 

# ask the user where they wants miniconda or anaconda 
if [[ -z $choice ]]; then	
	echo 
	echo "Select Conda distribution:"
	echo " [1] Miniconda (roughly 80 MB) ~ beginners to conda"
	echo " [2] Anaconda (rouchly 4GB, bundles science stack) ~ peeps familiar with Conda" 
	read -rp "Choice 1/2 -> " choice 
fi

[[ $choice == 1 || $choice == 2 ]] || { echo "Abort - enter 1 or 2."; exit 1; }

# build download URL and local filename 
base_url=https://repo.anaconda.com 
file_miniconda=Miniconda3-latest-MacOSX-x86_64.sh 
file_anaconda=Anaconda3-latest-MacOSX-x86_64.sh 

if [[ $choice == 1 ]]; then # user picked Miniconda 
  inst_file=$file_miniconda
  url="$base_url/miniconda/$inst_file" 
  prefix="$HOME/miniconda3" # install here 
else
  inst_file=$file_anaconda
  url="$base_url/archive/$inst_file"
  prefix="$HOME/anaconda3" # install here ..... 
fi 


# download if not cached - idempotent which skips redownload on re-run 
[[ -f $inst_file ]] || curl -L "$url" -o "$inst_file" 

# run installer under ROsetta when needed 
run_cmd="bash"
[[ $(detect_arch) == arm64 ]] && run_cmd="arch -x86_64 bash" 
$run_cmd "$inst_file" -b -p "$prefix"
export PATH="$prefix/bin:$PATH"
echo "[installer] ${inst_file%%-*} installed to $prefix"

# patch ~/.condarc only when requried 
if $patch_needed; then 
  echo -e "\n# Pinned by MicroSeq installer\ nsubdir: osx-64" >> ~/.condarc 
  echo "[installer] Wrote osx-64 subdir to ~/.condarc" 
fi

# --- clone repo and build env and run wizard ------------
[[ -d MicroSeq ]] || git clone https://github.com/jose-cantu/MicroSeq.git # clone once 
cd MIcroSeq # enter repo 

# make 'conda acitvate' available inside a sourced script 
source "$(conda info --base)/etc/profile.d/conda.sh" 

# create env only if it doesn't exist 
conda env list | grep -q '^MicroSeq ' || conda env create -f config/environment.yml -n MicroSeq 

conda activate MicroSeq # put env's Python on PATH 
pip install -e . # editable install keeps repo and CLI in sync 

echo "[installer] running microseq-setup (may take 3-5 min) please wait until prompt which will take some time to show up ~ 1min" 
if [[ -n ${CI-} ]]; then # CI -> non-interactive wizard for github actions test 
  microseq-setup --quiet 
else 
  microseq-setup # local user that I want full prompts for 
fi 

