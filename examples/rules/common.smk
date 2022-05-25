def globgfa(dir):
	return [(f.split(".gfa"))[0] for f in os.listdir(dir) if f.endswith(".gfa")]

def globgfa_hap(dir):
	fn = [(f.path,f.name) for f in os.scandir(dir) if f.name.endswith(".gfa")]
	return [(n.split(".gfa"))[0] for (p,n) in fn if subprocess.getoutput(DATADIR + "/xhap.sh -p H " + p) != ""]

GUROBI="gurobi_cl"
if "gurobi_root" in config:
	GUROBI = f"GRB_LICENSE_FILE={config['gurobi_lic']} {GUROBI}"
	GUROBI = f"PATH=$PATH:{config['gurobi_root']}/bin {GUROBI}"
	GUROBI = f"LD_LIBRARY_PATH={config['gurobi_root']}/lib {GUROBI}"

DATADIR = workflow.basedir + "/" + config.get('datadir')
SHDIR = workflow.basedir + "/" + config.get('shelldir')
BINDIR = workflow.basedir + "/" + config.get('bindir')
RUST_LOG = "debug" if config['debug'] else "info"
RUST_BIN = f"{BINDIR}/{'debug' if config['debug'] else 'release'}"
if not "solve_time_limit" in config:
	MAXSOLTIME = "Infinity"
else:
	MAXSOLTIME = 60 * config['solve_time_limit']

# ouch
if not "STRIP" in globals():
	STRIP = ""


rule write_founder_flow_lp:
	input:
		f"{{sample}}.gfa"
	output:
		f"{{sample}}.flow.lp"
	log:
		f"{{sample}}.flow.lp.log"
	benchmark:
		f"{{sample}}.flow.lp.perf"
	shell:
		f"export RUST_LOG={RUST_LOG}; "
		f"{RUST_BIN}/mkflow {{input}} >{{output}} 2>{{log}}"

rule solve_founder_flow:
	input:
		f"{{sample}}.flow.lp"
	output:
		f"{{sample}}.flow.sol"
	log:
		f"{{sample}}.flow.sol.log"
	benchmark:
		f"{{sample}}.flow.sol.perf"
	shell:
		f"{GUROBI} ResultFile={{output}} LogFile={{log}} Threads=1 {{input}} >/dev/null"

rule construct_founder_seqs:
	input:
		sol = f"{{sample}}.flow.sol",
		graph = f"{{sample}}.gfa",
	output:
		f"{{sample}}.founders.noopt.txt"
	log:
		f"{{sample}}.founders.noopt.log"
	benchmark:
		f"{{sample}}.founders.noopt.perf"
	shell:
		f"export RUST_LOG={RUST_LOG}; "
		f"{RUST_BIN}/flow2seq {{input.sol}} >{{output}} 2>{{log}}; "
		f"{SHDIR}/check_founder_solution.py "
			f"{{input.graph}} {{output}} >>{{log}}"

rule get_recombination_number:
	input:
		fnd = f"{{sample}}.{STRIP}founders.noopt.txt",
		haps = f"{{sample}}.haps.tsv",
	output:
		f"{{sample}}.nrecomb.txt"
	log:
		f"{{sample}}.nrecomb.log"
	benchmark:
		f"{{sample}}.nrecomb.perf"
	shell:
		f"export RUST_LOG={RUST_LOG}; "
		f"{RUST_BIN}/min_random -n 100000 {{input.fnd}} {{input.haps}} "
			f">{{output}} 2>{{log}}"

rule extract_haplotypes:
	input:
		f"{{sample}}.gfa"
	output:
		f"{{sample}}.haps.tsv"
	shell:
		f"{SHDIR}/xhap.sh -p {config['xhap_regex']} "
			f"{{input}} >{{output}}; "
		f"test 0 -ne $(du -b {{output}} | cut -f1)"

rule write_minimization_lp:
	input:
		fnd = f"{{sample}}.{STRIP}founders.noopt.txt",
		haps = f"{{sample}}.haps.tsv",
		rc = rules.get_recombination_number.output,
	output:
		f"{{sample}}.min.lp"
	log:
		f"{{sample}}.min.lp.log"
	benchmark:
		f"{{sample}}.min.lp.perf"
	shell:
		f"export RUST_LOG={RUST_LOG}; "
		f"{RUST_BIN}/mkmin {{input.fnd}} {{input.haps}} "
			f">{{output}} 2>{{log}}"

rule solve_minimization:
	input:
		f"{{sample}}.min.lp"
	output:
		f"{{sample}}.min.sol"
	log:
		f"{{sample}}.min.sol.log"
	benchmark:
		f"{{sample}}.min.sol.perf"
	shell:
		f"{GUROBI} ResultFile={{output}} LogFile={{log}} TimeLimit={MAXSOLTIME} Threads=1 "
		f"{{input}} >/dev/null"

rule construct_minimal_founders:
	input:
		f"{{sample}}.min.sol"
	output:
		f"{{sample}}.founders.final.txt"
	log:
		f"{{sample}}.founders.final.log"
	benchmark:
		f"{{sample}}.founders.final.perf"
	shell:
		f"export RUST_LOG={RUST_LOG}; "
		f"{RUST_BIN}/min2seq {{input}} >{{output}} 2>{{log}}"

rule compress_flow:
	input:
		f = f"{{sample}}.flow.lp",
		c = f"{{sample}}.founders.noopt.txt"
	output:
		f"{{sample}}.flow.lp.gz"
	shell:
		f"gzip -9 {{input.f}}"

rule compress_min:
	input:
		f = f"{{sample}}.min.lp",
		c = f"{{sample}}.founders.final.txt"
	output:
		f"{{sample}}.min.lp.gz"
	shell:
		f"gzip -9 {{input.f}}"

rule compress_flow_zq:
	input:
		f = f"{{sample}}.flow.lp",
		c = f"{{sample}}.founders.noopt.txt"
	output:
		f"{{sample}}.flow.lp.zpaq"
	shell:
		f"zpaq a {{output}} {{input.f}} -m1 -t1 >/dev/null 2>&1 "
		f"&& rm {{input.f}}"

rule compress_min_zq:
	input:
		f = f"{{sample}}.min.lp",
		c = f"{{sample}}.founders.final.txt"
	output:
		f"{{sample}}.min.lp.zpaq"
	shell:
		f"zpaq a {{output}} {{input.f}} -m1 -t1 >/dev/null 2>&1 "
		f"&& rm {{input.f}}"
