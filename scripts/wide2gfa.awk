#!/usr/bin/env -S awk -f
# usage: $0 [WIDEWALKFILE..]
# convert wide ("horizontal", default) format founder output to full gfa1
NF == 1{
	fn = $1
	fi = 0
}
NF != 1{
	subp = sprintf("%s_%03d_%s", fn, fi++, $1)
	gsub(">", " >")
	gsub("<", " <")
	last = ""
	for(i=2; i<=NF; i++){
		s = substr($i, 2)
		S[s] = 1
		d = substr($i, 1, 1) == ">" ? "+" : "-"
		l = s d
		if(last != ""){
			L[last "\t" l] = 1
			P[subp] = P[subp] "," l
		}else
			P[subp] = l
		last = l
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
