#!/bin/sh
# compare lengths of the F0 founder from the simulated data (GFA)
# to the reconstruction (walk tsv)
if [ $# -ne 2 ]; then
	echo usage: "$0 sim.gfa founders.txt" 1>&2
	exit 1
fi
x=$(grep founder_seq0 $2 |
	tr '<>' \\n |
	sed 1d |
	wc -l
)
y=$(grep F0 $1 |
	tr , \\n |
	wc -l
)
echo -n $2 $x $1 $y:' '
(test $x -le $y && echo shorter or equal || echo longer)
