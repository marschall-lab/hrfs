rule check_length:
	input:
		ref = f"{{sample}}.gfa",
		f1 = f"{{sample}}.flow.founders.txt",
		f2 = f"{{sample}}.min.founders.txt",
	output:
		f"{{sample}}.founders.stats",
	shell:
		f"{SHDIR}/fndcmp.sh {{input.ref}} {{input.f1}} >{{output}}"
		f"; {SHDIR}/fndcmp.sh {{input.ref}} {{input.f2}} >>{{output}}"
