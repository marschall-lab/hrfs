# Template Switch Founder Set

Analysis of the set of founder sequences under the template switch model

## Requirements

* `rust` version >= 1.60
* `python` version >= 3.7
* `snakemake`

## Build software

```
cargo build --manifest-path Cargo.toml --release --offline
```

## Run examples

Example with 4 available CPU cores.

### Simulation experiments

```
cd examples/experiments/sim
snakemake -j 4
```

### Example from paper

```
cd examples/experiments/examples
snakemake -k -j 4
```

### Clean up

Use the `clean` snakemake target:

```
snakemake -k -j 1 clean
```


## Usage

The following uses the `examples' experiment as reference.  It demonstrates the
software's typical usage with the provided `snakemake` workflows.

### Configuration

Experiments reside in their own respective directories and are configured via a `config.yaml` file, used to configure simulation and analysis parameters.
Paths should be left as-is unless changing directory structure.
Remaining recognized parameters:

- `debug` (boolean): toggles verbose debugging output
- `xhap_regex` (string): regular expression used to select haplotype paths in the input GFA files
- `solve_time_limit` (integer, minutes): time limit for the `gurobi` optimization steps
- `nnodes` (integer list): number of nodes in the graph
- `dup_ratio` (list of floats ∈ [0;1]): duplications ratio
- `inv_ratio` (list of floats ∈ [0;1]): inversions among duplications ratio
- `nhaplotypes` (integer): number of haplotypes to generate
- `nsamples` (integer): number of replicates per parameter set

Data used by the experiment should reside in a subdirectory under `examples/data`.


### Programs

- `hapsim`: generate simulated founder set, haplotypes, and their variation graph
- `subgr`: select subset of haplotypes and resulting subgraph from a GFA file
- `mkflow`: write to file flow linear program to solve
- `flow2seq`: reconstruct founder set sequences from flow solution
- `min_random`: estimate number of recombinations in flow solution by random assignment trials
- `mkmin`: write to file minimization program to solve
- `min2seq`: reconstruct founder set sequences from minimization solution


### Output

Most relevant output, by file extension:

- `.gfa`: user-provided GFA, or one generated by the simulator
- `.lp`, `.sol`: linear program and solution of flow program and recombination minimization
- `.nrecomb.txt`: number of recombinations in flow solution after random assignment trials
- `.founders.noopt.txt`: minimal founder sequences set reconstructed from flow solution
- `.founders.final.txt`: minimal founder sequences set after minimizing their number of recombinations


## License

This software is distributed under the MIT license.  For more details, see the
LICENSE file.
