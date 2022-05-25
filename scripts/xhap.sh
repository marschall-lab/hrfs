#!/bin/bash
# extract haplotypes from gfa with +- path format and convert to <> format
# usage: $0 [-p pattern] [gfa..]
# -p: pattern to filter path records by
PAT=H
if [ $1 = "-p" ]; then
	PAT=$2
		shift 2
fi
if [ -z $PAT ]; then
	echo usage: "$0 [-p pattern] [gfa..]" 1>&2
	exit 1
fi
awk -v "p=^$PAT" -v 'OFS=\t' '
$1 == "P" && $2 ~ p{
	n = split($3, a, ",")
	r = ""
	for(i=1; i<=n; i++){
		m = split(a[i], b, "+")
		s = ">"
		if(m == 1){
			m = split(a[i], b, "-")
			s = "<"
		}
		r = r s b[1]
	}
	print $2, r
}
' "$@"
