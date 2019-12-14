#!/bin/bash
#$ -S /bin/bash
# start on current working dir
#$ -cwd
#$ -V
# job name
#$ -N wien2k-w2w
#$ -j Y
#$ -q all.q
# stdout/stderr
#$ -o stdout
#$ -e stderr
# send to slack
#$ -M o1i3r8g8z8q9y5a9@kurokigroup.slack.com
# mail at beginning/end/on suspension
#$ -m bes
