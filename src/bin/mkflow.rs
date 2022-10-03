/* standard use */
use std::io;
use std::io::prelude::*;

/* crate use */
use clap::Parser;
use gfa::{gfa::GFA, parser::GFAParser};
use handlegraph::{
    handle::{Direction, Edge, Handle},
    handlegraph::*,
    hashgraph::HashGraph,
};
use rustc_hash::FxHashSet;

/* private use */
use founderset as ff;

#[derive(clap::Parser, Debug)]
#[clap(
    version = "0.1",
    author = "Daniel Doerr <daniel.doerr@hhu.de>",
    about = "Generate linear program to compute size of founder set"
)]
pub struct Command {
    #[clap(
        short = 'f',
        long = "nfounder",
        help = "Force fixed constant number of founder sequences"
    )]
    nfounder: Option<usize>,

    #[clap(index = 1, help = "graph in GFA1 format", required = true)]
    pub graph: String,
}

fn write_lp<W: io::Write>(
    graph: &HashGraph,
    nfounder: Option<usize>,
    out: &mut io::BufWriter<W>,
) -> Result<(), io::Error> {
    // find sources & sinks
    //
    // - sources are all nodes that in default orientation have no edges to their left
    // - sinks are all nodes that in default orientation have no edges to their right
    let mut sources: FxHashSet<Handle> = FxHashSet::default();
    let mut sinks: Vec<Handle> = Vec::new();
    for v in graph.handles() {
        if graph.degree(v, Direction::Right) == 0 {
            sinks.push(v);
        } else if graph.degree(v, Direction::Left) == 0 {
            sources.insert(v);
        }
    }

    log::info!(
        "identified {} sources: {}\nidentified {} sinks: {}",
        sources.len(),
        ff::hv2seq(&sources, ","),
        sinks.len(),
        ff::v2seq(&sinks, ",")
    );

    // write objective
    writeln!(out, "Minimize")?;
    writeln!(
        out,
        "{}",
        graph
            .edges()
            .map(|Edge(u, v)| {
                if u == v.flip() {
                    format!("o{}_i{}", ff::v2extstr(&u), ff::v2extstr(&v.flip()))
                } else {
                    format!(
                        "o{}_i{} + o{}_i{}",
                        ff::v2extstr(&u),
                        ff::v2extstr(&v.flip()),
                        ff::v2extstr(&v.flip()),
                        ff::v2extstr(&u)
                    )
                }
            })
            .collect::<Vec<String>>()
            .join(" + ")
    )?;

    writeln!(out, "Subject To")?;

    // write flow equations
    for v in graph.handles() {
        //
        // match flow with capacity
        //

        if graph.degree(v, Direction::Left) > 0 {
            writeln!(
                out,
                "i{} - {} = 0",
                ff::v2extstr(&v.flip()),
                graph
                    .neighbors(v, Direction::Left)
                    .map(|u| format!("o{}_i{}", ff::v2extstr(&u), ff::v2extstr(&v.flip())))
                    .collect::<Vec<String>>()
                    .join(" - ")
            )?;
            writeln!(
                out,
                "o{} - {} = 0",
                ff::v2extstr(&v.flip()),
                graph
                    .neighbors(v, Direction::Left)
                    .map(|u| format!("o{}_i{}", ff::v2extstr(&v.flip()), ff::v2extstr(&u)))
                    .collect::<Vec<String>>()
                    .join(" - ")
            )?;
        }

        if graph.degree(v, Direction::Right) > 0 {
            writeln!(
                out,
                "i{} - {} = 0",
                ff::v2extstr(&v),
                graph
                    .neighbors(v, Direction::Right)
                    .map(|u| format!("o{}_i{}", ff::v2extstr(&u.flip()), ff::v2extstr(&v)))
                    .collect::<Vec<String>>()
                    .join(" - ")
            )?;
            writeln!(
                out,
                "o{} - {} = 0",
                ff::v2extstr(&v),
                graph
                    .neighbors(v, Direction::Right)
                    .map(|u| format!("o{}_i{}", ff::v2extstr(&v), ff::v2extstr(&u.flip())))
                    .collect::<Vec<String>>()
                    .join(" - ")
            )?;
        }

        //
        // ensure flow conservation
        //

        writeln!(
            out,
            "i{} - o{} = 0",
            ff::v2extstr(&v),
            ff::v2extstr(&v.flip())
        )?;
        writeln!(
            out,
            "i{} - o{} = 0",
            ff::v2extstr(&v.flip()),
            ff::v2extstr(&v)
        )?;
    }

    //
    // write flow constraints
    //
    for Edge(u, v) in graph.edges() {
        if u == v.flip() {
            log::debug!(
                "edge {}-{} is self-edge",
                ff::v2extstr(&u),
                ff::v2extstr(&v.flip())
            );
            writeln!(
                out,
                "o{}_i{} >= 1",
                ff::v2extstr(&u),
                ff::v2extstr(&v.flip())
            )?;
        } else {
            writeln!(
                out,
                "o{}_i{} + o{}_i{} >= 1",
                ff::v2extstr(&u),
                ff::v2extstr(&v.flip()),
                ff::v2extstr(&v.flip()),
                ff::v2extstr(&u)
            )?;
        }
    }

    //
    // write sink constrains
    //
    for v in sinks.iter() {
        writeln!(out, "i{} = 0", ff::v2extstr(&v))?;
    }
    if let Some(nf) = nfounder {
        if nf > 0 {
            writeln!(
                out,
                "{}",
                format!(
                    "{} = {}",
                    sinks
                        .iter()
                        .map(|v| "i".to_owned() + &ff::v2extstr(&v.flip()))
                        .collect::<Vec<_>>()
                        .join(" + "),
                    nf
                )
            )?;
        }
    }

    //
    // write source constrains
    //
    for v in sources.iter() {
        writeln!(out, "o{} = 0", ff::v2extstr(&v.flip()))?;
    }

    for v in graph.handles() {
        writeln!(out, "o{} >= 0", ff::v2extstr(&v))?;
        writeln!(out, "i{} >= 0", ff::v2extstr(&v))?;
        writeln!(out, "o{} >= 0", ff::v2extstr(&v.flip()))?;
        writeln!(out, "i{} >= 0", ff::v2extstr(&v.flip()))?;
    }

    for Edge(u, v) in graph.edges() {
        writeln!(
            out,
            "o{}_i{} >= 0",
            ff::v2extstr(&u),
            ff::v2extstr(&v.flip())
        )?;
        if u != v.flip() {
            writeln!(
                out,
                "o{}_i{} >= 0",
                ff::v2extstr(&v.flip()),
                ff::v2extstr(&u)
            )?;
        }
    }

    //
    // write integer constraints
    //
    writeln!(out, "\nGeneral")?;
    for v in graph.handles() {
        writeln!(out, "o{}", ff::v2extstr(&v))?;
        writeln!(out, "i{}", ff::v2extstr(&v))?;
        writeln!(out, "o{}", ff::v2extstr(&v.flip()))?;
        writeln!(out, "i{}", ff::v2extstr(&v.flip()))?;
    }

    for Edge(u, v) in graph.edges() {
        writeln!(out, "o{}_i{}", ff::v2extstr(&u), ff::v2extstr(&v.flip()))?;
        if u != v.flip() {
            writeln!(out, "o{}_i{}", ff::v2extstr(&v.flip()), ff::v2extstr(&u))?;
        }
    }
    writeln!(out, "End")?;

    Ok(())
}

fn main() -> Result<(), io::Error> {
    env_logger::init();

    // print output to stdout
    let mut out = io::BufWriter::new(std::io::stdout());

    // initialize command line parser & parse command line arguments
    let params = Command::parse();

    log::info!("loading graph {}", &params.graph);
    let parser = GFAParser::new();
    let gfa: GFA<usize, ()> = parser.parse_file(&params.graph).unwrap();

    log::info!("constructing handle graph");
    let graph = HashGraph::from_gfa(&gfa);

    log::info!("writing linear program");
    write_lp(&graph, params.nfounder, &mut out)?;
    out.flush()?;

    log::info!("done");
    Ok(())
}
