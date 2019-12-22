#!/bin/bash
#$ -S /bin/bash
# start on current working dir
#$ -cwd
#$ -V
# job name
#$ -N boltztrap
#$ -j Y
#$ -q all.q
# send to slack
#$ -M o1i3r8g8z8q9y5a9@kurokigroup.slack.com
# Mail at beginning/end/on suspension
#$ -m bes

test $# -eq 1 || exit 1
case=$1

# BoltzTraP
# =========
n_k=100000
# 0: no shift, 1: shift k-mesh
shift_k=1 # shift

# ENVIRONMENT
# ===========
W2K=/usr/local/WIEN2k_18.2
W2T=~/w2ktools/bin/w2t

set -eu
trap 'echo ERROR: $0:$LINENO exit $?' ERR INT

get_ef () {
  grep '^:FER' "${case}.scf2" | sed -r 's/^.+=\s+(.+)$/\1/'
}

get_num_elec () {
  grep '^:NOE' "${case}.scf2" | sed -r 's/^.+=\s+(.+)$/\1/'
}

make_intrans () {
  cat << EOF > "${case}.intrans"
WIEN       # Format of DOS
0 0 0 0.0  # iskip (not presently used), idebug, setgap, shiftgap
$(get_ef) 0.0005 0.3 $(get_num_elec) # Fermilevel (Ry), energygrid, energy span around Fermilevel, #of electrons
CALC       # CALC (calculate expansion coeff) | NOCALC (read from file)
5          # lpfac (number of latt-points per k-point)
BOLTZ      # run mode (only BOLTZ is supported)
0.15       # efcut (energy range of chemical potential)
800. 50.   # Tmax, temperature grid
-1         # energyrange of bands given individual DOS output sig_xxx and dos_xxx (xxx is band number)
HISTO
0 0 0 0 0
2
1E20 -1E20
EOF
}

echo "==> BoltzTraP for case: $case"
cd $case
echo -e "$n_k\n$shift_k" | $W2K/x kgen

echo "--> x lapw1"
$W2K/x lapw1

make_intrans
echo "--> x_trans BoltzTraP"
$W2K/x_trans BoltzTraP

echo "--> DONE (BoltzTraP for case: $case)"
