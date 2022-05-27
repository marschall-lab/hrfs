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
