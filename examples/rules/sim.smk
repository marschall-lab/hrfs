def fprange(a, b, step):
	while a <= b:
		yield a
		a = round(a + step, 5)

def sfill(x, max):
	return str(x).zfill(len(str(abs(max))))

def prange(c):
	return fprange(config[c][0], config[c][1], config[c][2])

def sampnames(c):
	return [sfill(s, config[c]) for s in range(1, config[c]+1)]

rule check_length:
	input:
		ref = f"{{sample}}.gfa",
		f1 = f"{{sample}}.founders.noopt.txt",
		f2 = f"{{sample}}.founders.final.txt",
	output:
		f"{{sample}}.founders.stats",
	shell:
		f"{SHDIR}/fndcmp.sh {{input.ref}} {{input.f1}} >{{output}}; "
		f"{SHDIR}/fndcmp.sh {{input.ref}} {{input.f2}} >>{{output}}"

rule generate_graph:
	output:
		"sim.{nnodes}_{rdup}_{rinv}_{nhap}_{i}.gfa",
	log:
		"sim.{nnodes}_{rdup}_{rinv}_{nhap}_{i}.gfa.log",
	benchmark:
		"sim.{nnodes}_{rdup}_{rinv}_{nhap}_{i}.gfa.perf",
	run:
		shell("export RUST_LOG={RUST_LOG}; "
			"{RUST_BIN}/hapsim -l {wildcards.nnodes} "
			"-d {wildcards.rdup} -r {wildcards.rinv} "
			"{wildcards.nhap} >{output} 2>{log}"
		)
