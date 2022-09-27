#!/bin/bash -e
if [ $# -ne 2 ]; then
	echo usage: "$0 RAWGFA XTRACTREGEX" 1>&2
	exit 1
fi
graph=$1
regex=$2

# check if any haplotypes can be extracted and used at all
if [[ `sed -n "/^P\s$2/p" $1 | wc -l` == 0 ]]; then
	echo "no haplotypes in input gfa, or none matching regex!" >1&2
	exit 2
fi

# check for multiple sinks/sources
awk '
/^P/{
	gsub("[\\-\\+,]"," ",$3)
	n = split($3, a, " ")
	s[a[1]] = 1
	S[a[n]] = 1
} END {
	if(length(s) != 1 || length(S) != 1){
		print "unhandled multiple sources and/or sinks!"
		for(i in s)
			print "s", i
		for(i in S)
			print "S", i
		exit 3
	}
}' $1
