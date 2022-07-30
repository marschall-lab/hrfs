#!/bin/awk -f
# usage: $0 [WIDEWALKFILE..]
# convert wide ("horizontal", default) format founder output to full gfa1
NF == 1{
	fn = $1
	fi = 0
}
NF != 1{
	subp = fn "_" fi++ "_" $1
	split($2, a, ">")
	last = ""
	for(i=1; i<=length(a); i++){
		split(a[i], b, "<")
		for(j=1; j<=length(b); j++){
			S[b[j]] = 1
			# first element is always +, all others are always -
			dir = j == 1 ? "+" : "-"
			new = b[j] "\t" dir
			if(last != ""){
				L[last "\t" new] = 1
				P[subp] = P[subp] "," b[j] dir
			}else
				P[subp] = b[j] dir
			last = new
		}
	}
}
END{
	OFS = "\t"
	sort = "sort -V"
	print "H", "VN:Z:1.0"
	for(i in S)
		print "S", i, "*" | sort
	close(sort)
	for(i in L)
		print "L", i, "0M" | sort
	close(sort)
	for(i in P)
		print "P", i, P[i], "*" | sort
}
