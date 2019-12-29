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

[[ $# -ne 2 ]] && exit 1
w2k_case=$1
case=$2

# WIEN2WANNIER
# ============
n_k=1000 # num. of k points
shift_k=0 # shift k points, should be 0

# ENVIRONMENT
# ===========
W2K=/usr/local/WIEN2k_18.2
W2T=~/w2ktools/bin/w2t

set -eu
trap 'echo ERROR: $0:$LINENO exit $?' ERR INT

echo "==> WIEN2WANNIER for case: $case"
source $case.w2w
$with_so && so="-so" || so=

cd $w2k_case
$W2K/prepare_w2wdir $case
cd $case

echo "--> generate k-mesh (kgen)"
echo -e "$n_k\n$shift_k" | $W2K/x kgen $so -fbz

echo "--> get band indices"
$W2K/x findbands $so -all $e_min $e_max
cat $case.outputfind | tail -n3

line=$(cat $case.outputfind | tail -n1 | sed 's/^at any k:\s*//')
i_band_min=$(echo $line | awk '{print $1}')
i_band_max=$(echo $line | awk '{print $2}')

echo "--> write $case.inwf"
if $with_so; then
  $W2K/write_inwf -bands $i_band_min $i_band_max $projs $projs
else
  $W2K/write_inwf -bands $i_band_min $i_band_max $projs
fi

echo "--> write $case.win"
$W2K/write_win -fresh
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

if $with_so; then
  echo "--> update energy and vector (lapwso)"
  $W2K/x lapwso
fi

echo "--> compute initial projection"
if $with_so; then
  $W2K/x w2w -so -up; $W2K/x w2w -so -dn
else
  $W2K/x w2w
fi

echo "--> compute MLWF (wannier90)"
$W2K/x wannier90 $so

echo "--> DONE (Wien2Wannier for case: $case)"
# echo "$z $(~/anaconda3/bin/python3 ../get_hopping.py $case)" >> ../hoppings.dat
