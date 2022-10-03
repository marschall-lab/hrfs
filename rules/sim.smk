# define an interval instead of a list of values
def fseq(a, b, step):
	while a <= b:
		yield a
		a = round(a + step, 5)

def frange(c):
	return fseq(config[c][0], config[c][1], config[c][2])

# generate all simulation sample names
def sfill(x, max):
	return str(x).zfill(len(str(abs(max))))

def sampnames(n):
	return [sfill(s, n) for s in range(1, n+1)]

rule generate_graph:
	output:
		f"{DATADIR}/sim.{{nnodes}}_{{rdup}}_{{rinv}}_{{nhap}}_{{i}}.gfa"
	log:
		f"{OUTDIR}/sim/sim.{{nnodes}}_{{rdup}}_{{rinv}}_{{nhap}}_{{i}}.log"
	benchmark:
		f"{OUTDIR}/sim/sim.{{nnodes}}_{{rdup}}_{{rinv}}_{{nhap}}_{{i}}.prof"
	shell:
		f"{RUSTBIN}/hapsim"
		f"	-l {{wildcards.nnodes}}"
		f"	-d {{wildcards.rdup}}"
		f"	-r {{wildcards.rinv}}"
		f"	{{wildcards.nhap}}"
		f"	>{{output}} 2>{{log}}"

rule generate:
	input:
		expand(f"{DATADIR}/sim.{{nnodes}}_{{rdup}}_{{rinv}}_{{nhap}}_{{i}}.gfa",
			nnodes = config["nnodes"],
			rdup = config["dup_ratio"],
			rinv = config["inv_ratio"],
			nhap = config["nhaplotypes"],
			i = sampnames(config["nsamples"]),
		),

rule expandyourface:
	input:
		expand(f"{DATADIR}/sim.{{nnodes}}_{{rdup}}_{{rinv}}_{{nhap}}_{{i}}.gfa",
			nnodes = config["nnodes"],
			rdup = config["dup_ratio"],
			rinv = config["inv_ratio"],
			nhap = config["nhaplotypes"],
			i = sampnames(config["nsamples"]),
		),
		expand(f"{OUTDIR}/{{nf}}/{{sample}}.min.full.gfa",
			nf = config.get("nforced", 0),
			sample = expand(f"sim.{{nnodes}}_{{rdup}}_{{rinv}}_{{nhap}}_{{i}}",
				nnodes = config["nnodes"],
				rdup = config["dup_ratio"],
				rinv = config["inv_ratio"],
				nhap = config["nhaplotypes"],
				i = sampnames(config["nsamples"]),
			),
		),
