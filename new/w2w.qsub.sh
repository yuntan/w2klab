#!/bin/bash
#$ -S /bin/bash
# start on current working dir
#$ -cwd
#$ -V
# job name
#$ -N w2w
#$ -j Y
#$ -q all.q
# send to slack
#$ -M o1i3r8g8z8q9y5a9@kurokigroup.slack.com
# mail at beginning/end/on suspension
#$ -m bes

[[ $# -eq 2 ]] || exit 1

w2k_case=$1
case=$2

# Wien2Wannier
# ============
n_k=1000 # num. of k points
shift=0 # shift k points, should be 0

# ENVIRONMENT
# ===========
W2K=/usr/local/WIEN2k_18.2
W2T=~/w2ktools/bin/w2t

set -eu
trap 'echo ERROR: $0:$LINENO exit $?' ERR INT

echo "==> Wien2Wannier for case: $case"
source "$case.w2w"
cd $w2k_case
$W2K/prepare_w2wdir $case
cd $case

echo "--> generate k-mesh (kgen)"
echo -e "$n_k\n$shift" | $W2K/x kgen -fbz

echo "--> get band indices"
$W2K/x findbands -all $emin $emax
cat $case.outputfind | tail -n3

line=$(cat $case.outputfind | tail -n1 | sed 's/^at any k:\s*//')
i_band_min=$(echo $line | awk '{print $1}')
i_band_max=$(echo $line | awk '{print $2}')

echo "--> generate $case.inwf"
$W2K/write_inwf -bands $i_band_min $i_band_max $projs

echo "--> generate $case.win"
$W2K/write_win
cat - $case.win <<EOF >tmp
! $case.w2w
dis_froz_min = $dis_froz_min
dis_froz_max = $dis_froz_max

EOF
mv tmp $case.win

echo "--> update k-points"
$W2K/x wannier90 -pp

echo "--> update energy and vector (lapw1)"
$W2K/x lapw1

echo "--> compute initial projection"
$W2K/x w2w

echo "--> compute MLWF (wannier90)"
$W2K/x wannier90

echo "--> DONE (Wien2Wannier for case: $case)"
# echo "$z $(~/anaconda3/bin/python3 ../get_hopping.py $case)" >> ../hoppings.dat
