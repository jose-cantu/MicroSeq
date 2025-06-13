#!/usr/bin/env bash 

set -euo pipefail # exit on error, unset vars, pipe failures

detect_arch() { # this here is a funciton that detects the architechture of your machine 
  uname -m 
}

have_conda() { command -v conda >/dev/null 2>&1; } # return true if conda exists 

get_conda_subdir() { # asking existing conda what platform its solving for 
  conda config --show subdir 2>/dev/null | awk '/^subdir/ {print $2}' 
}

# --- prerequisite: ensure Git is available ----------------------------------
if ! command -v git >/dev/null 2>&1; then
  if [[ $(uname) == Darwin ]]; then
    echo "[installer] Git not found – installing Xcode Command-Line Tools…"
    xcode-select --install || true                # GUI pops up; no effect in CI
    echo "[installer] Re-run this script after the installer finishes Xcode Command Line Tools will help in installing git for you."
  else
    echo "[installer] Git not found."
    echo "           Run:  sudo apt install git   (or use your distro’s package manager)"
  fi
  exit 1
fi

patch_needed=false # default: assume no fix required 
if have_conda; then # inspect if conda exists 
   curr=$(get_conda_subdir || true) # may be empty or unset here 
   [[ -z $curr ]] && curr=osx-64 # empty means conda uses its default 
   if [[ $curr != osx-64 ]]; then # when key  is missing treat as intel osx-64 
	   echo "[installer] Found a non-Intel conda (subdir=$curr). Switching to Miniconda-x86_64." 
     patch_needed=true 
   fi
fi  # fi closes both if statements here

# --- helper: skip bootstrap if Intel conda already present -----------
bootstrap_done=false
_reuse_existing_conda() {                    # fn scope -> return legit
  if have_conda && [[ $(get_conda_subdir || echo osx-64) = osx-64 ]]; then
    echo "[installer] Existing osx-64 conda detected; skipping bootstrap."
    bootstrap_done=true                      # flip flag
  fi
}

# call helper early in main script flow
_reuse_existing_conda

run_cmd="bash" # default interpreter; ARM-Anaconda 
# block below can overired this later based on choice selected for conda selected 

# ---------------------------------------------------------------------
# run Miniconda / Anaconda bootstrap only when not already done
# ---------------------------------------------------------------------
if ! $bootstrap_done; then
  # Detect CI / non-interactive shell...
  if [[ ! -t 0 ]]; then                      # stdin is not a TTY -> CI
    choice=${MICROSEQ_CHOICE:-1}             # allow override via env; default 1
  else
    choice=""                                # interactive session -> will prompt
  fi

  # ask the user whether they want Miniconda or Anaconda
  if [[ -z $choice ]]; then
    echo
    echo "Select Conda distribution:"
    echo " [1] Miniconda (roughly 80 MB) ~ beginners to conda"
    echo " [2] Anaconda (roughly 4GB, bundles science stack) ~ peeps familiar with Conda"
    read -rp "Choice 1/2 -> " choice
  fi
  [[ $choice == 1 || $choice == 2 ]] || { echo "Abort - enter 1 or 2."; exit 1; }

  # helper: return latest filename for Mac x86_64 
  latest_anaconda() { # fn: echo newest installer 
      curl -s https://repo.anaconda.com/archive/ | # fetch HTML listing 
	      grep -oE 'Anaconda3-[0-9]{4}\.[0-9]+(-[0-9]+)?-MacOSX-x86_64\.sh' | # keep matching links using grep to find latest version  
	sort -V | tail -n 1 # version-aware sorting and newest version = last line 
} 

  # build download URL and local filename
  base_url=https://repo.anaconda.com
  file_miniconda=Miniconda3-latest-MacOSX-x86_64.sh
  file_anaconda=$(latest_anaconda) # dynamic filename 
  url="https://repo.anaconda.com/archive/$file_anaconda" 
  if [[ $choice == 1 ]]; then                # user picked Miniconda
    inst_file=$file_miniconda
    url="$base_url/miniconda/$inst_file"
    prefix="$HOME/miniconda3"                # install here
  else
    inst_file=$file_anaconda
    url="$base_url/archive/$inst_file"
    prefix="$HOME/anaconda3"                 # install here.....
  fi

  # adding a safeguard for Rosetta Anaconda on Apple-silicon only -------------  
  if [[ $choice == 2 ]]; then # Anaconda only here 
     if [[ $(uname -m) == arm64 ]] && ! pgrep -qx oahd;then 
	echo "[installer] Rosetta 2 required for the Intel-only Anaconda installer and I recognize you don't have that" 
	echo "            Will use... sudo softwareupdate --install-rosetta --agree-to-license" 
	echo "[installer] will now install rosetta for you so Anaconda installer will run..." 
	sudo /usr/sbin/softwareupdate --install-rosetta --agree-to-license 
     fi
     run_cmd=(arch -x86_64 /bin/bash) # Call Apple's universal bash making it expicit path (simple and deterministic avoiding crash with homebrew arm-only bash) 
   else 
     run_cmd=(bash) # Native interpreter for Miniconda or Intel host 
  fi 

  # download if not cached - idempotent which skips redownload on re-run
  # lookup official SHA-256 once per version (Anaconda published .sha256) 
  sha_url="$url.sha256" # append .sha256 to file URL 
  
  # ─── resume download ────────────────────────────────────────────────
  rc=0
  curl --retry 3 -C - -fL "$url" -o "$inst_file" || rc=$? # store curl code 
 
  if [[ $rc -ne 0 ]]; then 
     if [[ $rc -eq 22 || $rc -eq 33 ]]; then # 416 range error 
	echo "[installer] curl exit $rc (HTTP 416) verifying local file...." 
     else 
	exit $rc 
     fi 
  fi 

  # ── tier-1: per-file checksum ──────────────────────────────────────────────
  if curl -sfL "$sha_url" -o "$inst_file.sha"; then
    echo "[installer] downloaded per-file checksum"
  else
    echo "[installer] .sha256 missing – trying master SHA256SUMS"
    sums_url=https://repo.anaconda.com/archive/sha256sum.txt                # <- actual file name
    if curl -sfL "$sums_url" -o SHA256SUMS; then
        grep "  $inst_file\$" SHA256SUMS >"$inst_file.sha" || true
    fi
  fi

  # ── tier-3: warn if checksum still absent ─────────────────────────────────
  if [[ ! -f $inst_file.sha || ! -s $inst_file.sha ]]; then
    echo "[installer] WARNING: no checksum available – proceeding without hash verification"
  else
    if ! sha256sum -c "$inst_file.sha"; then
        echo "[installer] checksum failed – redownloading whole file"
        rm -f "$inst_file"
        curl --retry 3 -fL "$url" -o "$inst_file"
        sha256sum -c "$inst_file.sha"      # must pass now or abort
    fi
  fi


  "${run_cmd[@]}" "$inst_file" -b -p "$prefix" # exec installer directly... 
  source "$prefix/etc/profile.d/conda.sh"    # refresh functions in current shell
  export PATH="$prefix/bin:$PATH"
  echo "[installer] ${inst_file%%-*} installed to $prefix"

else
  echo "[installer] Using existing conda at $(command -v conda)"
fi

# patch ~/.condarc only when requried 
if $patch_needed; then 
  echo -e "\n# Pinned by MicroSeq installer\nsubdir: osx-64" >> ~/.condarc 
  echo "[installer] Wrote osx-64 subdir to ~/.condarc" 
fi

# --- clone repo (skip if already cloned/existing) and build env and run wizard ----  
if [[ -d .git ]]; then
  echo "[installer] Repo already present in current directory so no need to clone."
  repo_root="$PWD"
elif [[ -d MicroSeq/.git ]]; then
  echo "[installer] Repo ./MicroSeq already exists so reusing it."
  repo_root="MicroSeq"
else 
  echo "[installer] Clone MicroSeq..."
  git clone https://github.com/jose-cantu/MicroSeq.git 
  repo_root="MicroSeq" 
fi 
cd "$repo_root" # enter repo 

# make 'conda acitvate' available inside a sourced script 
source "$(conda info --base)/etc/profile.d/conda.sh"

conda_cmd="conda run --no-capture-stdin -n MicroSeq" # alias: env + stdin now fine...  

# create env only if it doesn't exist 
conda env list | grep -q '^MicroSeq ' || conda env create -f config/environment.yml -n MicroSeq 

conda activate MicroSeq # put env's Python on PATH

# ---- ensure Homebrew + expect is installed ----------
if [[ $(uname) == Darwin ]]; then 
  if ! command -v brew >/dev/null 2>&1; then # no homebrew 
    echo "[installer] Homebrew not found - soo installing that now for you."
    /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)" 
    eval "$(/opt/homebrew/bin/brew shellenv)" # add brew to PATH for arm Macs 
  fi

  if ! brew list --formula | grep -q '^expect$'; then  # if expect missing 
    echo "[installer] Installing 'expect' via Homebrew." 
    brew install expect 
  fi
fi

$conda_cmd pip install -e . # editable install keeps repo and CLI in sync (Pip install after so unbuffer works....)  

# if running in GitHub actions, expose this env's bin/ to later steps 
if [[ -n ${CI-} ]]; then # variable always set in Actions 
  echo "[installer] exporting $CONDA_PREFIX/bin to GITHUB_PATH"
  echo "$CONDA_PREFIX/bin" >> "$GITHUB_PATH" # makes microseq visible in next steps 
fi 


echo "[installer] running microseq-setup (may take 3-5 min) please wait until prompt which will take some time to show up ~ 1min" 
if [[ -n ${CI-} ]]; then # CI -> non-interactive wizard for github actions test 
  $conda_cmd microseq-setup --quiet 
else 
  $conda_cmd microseq-setup  # local user that I want full prompts for 
fi 

