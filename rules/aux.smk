rule clean:
	shell:
		f"rm -rf gurobi.log runs .err .out"

rule nuke:
	shell:
		f"rm -rf gurobi.log runs .err .out .snakemake __pycache__"

rule compress_flow:
	input:
		f"{OUTDIR}/{{pset}}/{{sample}}.flow.lp",
	output:
		f"{OUTDIR}/{{pset}}/{{sample}}.flow.lp.gz",
	shell:
		f"gzip -9f {{output}} {{input}}"

rule compress_min:
	input:
		f"{OUTDIR}/{{pset}}/{{sample}}.min.lp",
	output:
		f"{OUTDIR}/{{pset}}/{{sample}}.min.lp.gz",
	shell:
		f"gzip -9f {{output}} {{input}}"

rule compress_flow_zq:
	input:
		f"{OUTDIR}/{{pset}}/{{sample}}.flow.lp",
	output:
		f"{OUTDIR}/{{pset}}/{{sample}}.flow.lp.zpaq",
	threads:
		8
	shell:
		f"zpaq a {{output}} {{input}}"
		f"	-m2 -t8"
		"	>/dev/null 2>&1"
		f" && rm {{input}}"

rule compress_min_zq:
	input:
		f"{OUTDIR}/{{pset}}/{{sample}}.min.lp",
	output:
		f"{OUTDIR}/{{pset}}/{{sample}}.min.lp.zpaq",
	threads:
		8
	shell:
		f"zpaq a {{output}} {{input}}"
		f"	-m2 -t8"
		"	>/dev/null 2>&1"
		f" && rm {{input}}"
