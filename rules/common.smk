DEBUG = config.get("debug", False)
if DEBUG:
	os.environ["RUST_LOG"] = "debug"
	os.environ["RUST_BACKTRACE"] = "1"
else:
	os.environ["RUST_LOG"] = "info"
	os.environ["RUST_BACKTRACE"] = "0"
XHAPEXP = config.get("xhap_regex", ".")

BASEDIR = workflow.basedir
OUTDIR = os.path.join(BASEDIR, config.get("outdir", "."))
DATADIR = os.path.join(BASEDIR, "data")
REPODIR = config.get("repodir")
SHDIR = os.path.join(REPODIR, "scripts")
BINDIR = os.path.join(REPODIR, "bin")
RUSTBIN = os.path.join(REPODIR, "target", "debug" if DEBUG else "release")
ENVDIR = os.path.join(REPODIR, "env")

envvars:
	"RUST_LOG",
	"RUST_BACKTRACE",

wildcard_constraints:
	sample="[^/]+",
	pset="[^/]+",

include: "gurobi.smk"
include: "aux.smk"

rule all:
	input:
		expand(f"{OUTDIR}/{{nf}}/{{sample}}.min.full.gfa",
			nf = config.get("nforced", 0),
			sample = glob_wildcards(f"{DATADIR}/{{sample}}.gfa").sample,
		),
	default_target: True

rule check_input:
        input:
                f"{DATADIR}/{{sample}}.gfa"
        output:
                f"{OUTDIR}/.{{sample}}.ok"
        shell:
                f"{SHDIR}/check_input.sh"
                f"      {{input}}"
                f"      {XHAPEXP}"
                f"      >{{output}}"

rule extract_haplotypes:
	input:
		g = f"{DATADIR}/{{sample}}.gfa",
                o = f"{OUTDIR}/.{{sample}}.ok"
	output:
		f"{OUTDIR}/{{sample}}.haplotypes.txt"
	shell:
		f"{SHDIR}/xhap.sh "
		f"	-p {XHAPEXP}"
		f"	{{input.g}}"
		f"	>{{output}}"

rule write_founder_flow_lp:
	input:
		f"{DATADIR}/{{sample}}.gfa"
	output:
		f"{OUTDIR}/{{nf}}/{{sample}}.flow.lp"
	log:
		f"{OUTDIR}/{{nf}}/log/{{sample}}.flow.lp.log"
	benchmark:
		f"{OUTDIR}/{{nf}}/log/{{sample}}.flow.lp.prof"
	params:
		nf = lambda wc: f"-f {wc.nf}" if wc.nf != "" else ""
	shell:
		f"{RUSTBIN}/mkflow {{params}} {{input}} >{{output}} 2>{{log}}"

rule solve_founder_flow:
	input:
		f"{OUTDIR}/{{pset}}/{{sample}}.flow.lp"
	output:
		f"{OUTDIR}/{{pset}}/{{sample}}.flow.sol"
	log:
		f"{OUTDIR}/{{pset}}/log/{{sample}}.flow.sol.log"
	benchmark:
		f"{OUTDIR}/{{pset}}/log/{{sample}}.flow.sol.prof"
	envmodules:
		"gurobi",
	benchmark:
		f"{{sample}}.flow.sol.perf"
	threads:
		GUROBI_THREADS
	params:
		t = GUROBI_THREADS
	shell:
		f"gurobi_cl"
		f"	ResultFile={{output}}"
		f"	LogFile={{log}}"
		f"	Threads={{params.t}}"
		f"	{{input}} >/dev/null"

rule construct_founder_seqs:
	input:
		g = f"{DATADIR}/{{sample}}.gfa",
		s = f"{OUTDIR}/{{pset}}/{{sample}}.flow.sol",
	output:
		f"{OUTDIR}/{{pset}}/{{sample}}.flow.founders.txt"
	log:
		f"{OUTDIR}/{{pset}}/log/{{sample}}.flow.founders.log"
	benchmark:
		f"{OUTDIR}/{{pset}}/log/{{sample}}.flow.founders.prof"
	shell:
		f"{RUSTBIN}/flow2seq"
		f"	{{input.s}}"
		f"	>{{output}} 2>{{log}}"
		f"; {SHDIR}/check_founder_solution.py"
		f"	{{input.g}}"
		f"	{{output}} >>{{log}} 2>&1"

rule flow_to_gfa:
	input:
		f"{OUTDIR}/{{pset}}/{{sample}}.flow.founders.txt"
	output:
		f"{OUTDIR}/{{pset}}/{{sample}}.flow.founders.gfa"
	shell:
		f"{SHDIR}/walk2path.sh {{input}} >{{output}}"

rule get_recombination_number:
	input:
		f = f"{OUTDIR}/{{pset}}/{{sample}}.flow.founders.txt",
		h = f"{OUTDIR}/{{sample}}.haplotypes.txt",
	output:
		f"{OUTDIR}/{{pset}}/{{sample}}.flow.nrecomb.txt"
	log:
		f"{OUTDIR}/{{pset}}/log/{{sample}}.flow.nrecomb.log"
	benchmark:
		f"{OUTDIR}/{{pset}}/log/{{sample}}.flow.nrecomb.prof"
	params:
		ntrials = 100000
	shell:
		f"{RUSTBIN}/min_random"
		f"	-n {{params.ntrials}}"
		f"	{{input.f}}"
		f"	{{input.h}}"
		f"	>{{output}} 2>{{log}}"

rule write_minimization_lp:
	input:
		f = f"{OUTDIR}/{{pset}}/{{sample}}.flow.founders.txt",
		h = f"{OUTDIR}/{{sample}}.haplotypes.txt",
		r = f"{OUTDIR}/{{pset}}/{{sample}}.flow.nrecomb.txt",
	output:
		f"{OUTDIR}/{{pset}}/{{sample}}.min.lp"
	log:
		f"{OUTDIR}/{{pset}}/log/{{sample}}.min.lp.log"
	benchmark:
		f"{OUTDIR}/{{pset}}/log/{{sample}}.min.lp.prof"
	shell:
		f"{RUSTBIN}/mkmin"
		f"	{{input.f}}"
		f"	{{input.h}}"
		f"	>{{output}} 2>{{log}}"

rule solve_minimization:
	input:
		f"{OUTDIR}/{{pset}}/{{sample}}.min.lp"
	output:
		f"{OUTDIR}/{{pset}}/{{sample}}.min.sol"
	log:
		f"{OUTDIR}/{{pset}}/log/{{sample}}.min.sol.log"
	benchmark:
		f"{OUTDIR}/{{pset}}/log/{{sample}}.min.sol.prof"
	envmodules:
		"gurobi"
	resources:
		runtime_hrs = 12
	threads:
		GUROBI_THREADS
	params:
		t = GUROBI_THREADS
	shell:
		f"gurobi_cl"
		f"	Threads={{params.t}}"
		f"	ResultFile={{output}}"
		f"	TimeLimit={GUROBI_SOLTIME}"
		f"	LogFile={{log}}"
		f"	{{input}} >/dev/null"

rule construct_minimal_founders_output_compact:
	input:
		s = f"{OUTDIR}/{{pset}}/{{sample}}.min.sol",
		g = f"{DATADIR}/{{sample}}.gfa",
	output:
		f"{OUTDIR}/{{pset}}/{{sample}}.min.founders.compact.txt"
	log:
		f"{OUTDIR}/{{pset}}/log/{{sample}}.min.founders.compact.log"
	shell:
		f"{RUSTBIN}/min2seq -c"
		f"	{{input.s}}"
		f"	>{{output}} 2>{{log}}"

rule construct_minimal_founders_output_long:
	input:
		s = f"{OUTDIR}/{{pset}}/{{sample}}.min.sol",
		g = f"{DATADIR}/{{sample}}.gfa",
	output:
		f"{OUTDIR}/{{pset}}/{{sample}}.min.founders.long.txt"
	log:
		f"{OUTDIR}/{{pset}}/log/{{sample}}.min.founders.long.log"
	shell:
		f"{RUSTBIN}/min2seq -l"
		f"	{{input.s}}"
		f"	>{{output}} 2>{{log}}"

rule construct_minimal_founders_output_wide:
	input:
		s = f"{OUTDIR}/{{pset}}/{{sample}}.min.sol",
		g = f"{DATADIR}/{{sample}}.gfa",
	output:
		f"{OUTDIR}/{{pset}}/{{sample}}.min.founders.txt"
	log:
		f"{OUTDIR}/{{pset}}/log/{{sample}}.min.founders.log"
	shell:
		f"{RUSTBIN}/min2seq"
		f"	{{input.s}}"
		f"	>{{output}} 2>{{log}}"

rule construct_minimal_founders_output_gfa:
	input:
		f = f"{OUTDIR}/{{pset}}/{{sample}}.min.founders.txt",
		c = f"{OUTDIR}/{{pset}}/{{sample}}.min.founders.compact.txt",
		l = f"{OUTDIR}/{{pset}}/{{sample}}.min.founders.long.txt",
	output:
		f"{OUTDIR}/{{pset}}/{{sample}}.min.founders.gfa"
	shell:
		f"{SHDIR}/wide2gfa.awk {{input.f}} >{{output}}"

rule construct_minimal_founders_output_full_gfa:
	input:
		h = f"{DATADIR}/{{sample}}.gfa",
		m = f"{OUTDIR}/{{pset}}/{{sample}}.min.founders.gfa",
		c = f"{OUTDIR}/{{pset}}/{{sample}}.min.founders.compact.txt",
		f = f"{OUTDIR}/{{pset}}/{{sample}}.flow.founders.txt",
	output:
		f"{OUTDIR}/{{pset}}/{{sample}}.min.full.gfa",
	shell:
		f"cp {{input.h}} {{output}}"
		f"; {SHDIR}/walk2path.sh {{input.f}}"
		f"	| sed -n 's/^P\t/&flow_/p' >>{{output}}"
		f"; sed -n 's/^P\t/#&min_/p' {{input.m}} >>{{output}}"
		f"; {SHDIR}/walk2path.sh {{input.c}}"
		f"	| sed -n 's/^P\t/&min_/p' >>{{output}}"
