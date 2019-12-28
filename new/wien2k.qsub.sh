#!/bin/bash
#$ -S /bin/bash
# start on current working dir
#$ -cwd
#$ -V
# job name
#$ -N wien2k
#$ -j Y
#$ -q all.q
# send to slack
#$ -M o1i3r8g8z8q9y5a9@kurokigroup.slack.com
# mail at beginning/end/on suspension
#$ -m bes

test $# -eq 1 || exit 1
case=$1

# RUN FLAGS
# =========
run_init=true
run_scf=true
run_band=true
run_so_scf=true
run_so_band=true

# WIEN2k
# ======
# for init_lapw
rmt_red=0
vxc=13 # DO NOT CHANGE without a special reason!
ecut=-8.0
rkmax=7.0
mix=0.2
n_k=1000
# 0: no shift, 1: shift k-mesh
shift_k=0 # no shift

# for run_lapw
# energy convergence
ec=0.0001 # 0.1 mRy
# charge convergence
cc=0.000050 # 50 μe
n_iter=50
# regenerate in1 file after N iter
in1new= # disabled

# for gen_klist_band.py
pearson=tP # Pearson symbol
bz_syms=XGXMGZ # k path
bz_div_n=100 # BZの分割数

# ENVIRONMENT
# ===========
W2K=/usr/local/WIEN2k_18.2
W2T=~/w2ktools/bin/w2t

set -eu
trap 'echo ERROR: $0:$LINENO exit $?' ERR INT

ef () {
  cat $case.scf | grep -E '^:FER' | tail -n1 | perl -pe 's/^[^=]+=\s*//'
}

scf_args=
[[ -n "$ec" ]] && scf_args="$scf_args -ec $ec"
[[ -n "$cc" ]] && scf_args="$scf_args -cc $cc"
[[ -n "$n_iter" ]] && scf_args="$scf_args -i $n_iter"
[[ -n "$in1new" ]] && scf_args="$scf_args -in1new $in1new"

echo "==> WIEN2k for case: $case"

# ディレクトリ名とCASEを合わせる
mkdir $case &> /dev/null || true
cp $case.cif $case
cd $case

if $run_init; then
  echo "--> generate master input $case.struct"
  $W2K/cif2struct $case.cif

  echo "--> run init"
  $W2K/init_lapw -b -red $rmt_red -vxc $vxc -ecut $ecut -rkmax $rkmax -mix $mix -numk $n_k

  echo -e "$n_k\n$shift_k" | $W2K/x kgen
fi

if $run_scf; then
  echo "--> run SCF"
  $W2K/run_lapw $scf_args

  echo -n "cp "; cp -v $case.dayfile $case.dayfile_scf
fi

if $run_band; then
  echo "--> run band structure calc"
  # output: case.klist_band
  $W2T gen_klist_band $case $pearson $bz_syms $bz_div_n

  # input: case.klist_band
  # output: case.scf1, case.energy
  $W2K/x lapw1 -band
  # backup for restarting spaghetti
  echo -n "cp "; cp -v $case.energy $case.energy_band

  # update charges
  $W2K/x lapw2 -qtl -band

  rm $case.insp 2> /dev/null || true
  $W2K/x spaghetti || true # generate case.insp
  sed -i $case.insp -e 's/0\.xxxx/'$(ef)'/'
  $W2K/x spaghetti # generate case.spaghetti_ene

  # input: case.struct, case.klist_band, case.spaghetti_ene, case.qtl
  # output: spaghetti.npz
  $W2T spaghetti $case
fi

if $run_so_scf; then
  if [[ ! -f $case.inso ]]; then
    echo "ERROR $case.inso not found"
    echo "run initso to generate $case.inso"
    exit 1
  fi

  echo -e "$n_k\n$shift_k" | $W2K/x kgen -so

  echo "--> run SO SCF"
  $W2K/run_lapw -so $scf_args
fi

if $run_so_band; then
  echo "--> run SO band"
  # output: case.klist_band
  $W2T gen_klist_band $case $pearson $bz_syms $bz_div_n

  # input: case.klist_band
  # output: case.scf1, case.energy
  $W2K/x lapw1 -band

  # output: case.energyso
  $W2K/x lapwso
  echo -n "cp "; cp -v $case.energyso $case.energyso_band

  # update charges
  $W2K/x lapw2 -c -so -qtl -band

  rm $case.insp 2> /dev/null || true
  $W2K/x spaghetti || true # generate case.insp
  sed -i $case.insp -e 's/0\.xxxx/'$(ef)'/'
  $W2K/x spaghetti -so # generate case.spaghetti_ene

  # input: case.struct, case.klist_band, case.spaghetti_ene, case.qtl
  # output: spaghetti.npz
  $W2T spaghetti $case
fi

echo "--> DONE (WIEN2k for case: $case)"
