#!/bin/bash
# usage: $0 [WALKFILE..]
# >1>2 walk 1-liner notation to path entries in fully specified gfa file
sed 's/\([<>]\)/\n\1/g' "$@" \
	| awk \
'
$1 ~ /^[<>]/{
	dir = $1 ~ /^>/ ? "+" : "-";
	n = substr($1, 2)
	S[n] = 1
	if(i > 0){
		L[last "\t" n "\t" dir] = 1
		P[pn] = P[pn] "," n dir
	}else
		P[pn] = n dir
	last = n "\t" dir
	i++
}
$1 !~ /^[<>]/{
	pn = $1
	i = 0
}
END{
	OFS = "\t"
	print "H", "VN:Z:1.0"
	for(i in S)
		print "S", i, "*"
	for(i in L)
		print "L", i, "0M"
	for(i in P)
		print "P", i, P[i], "*"
}
'
