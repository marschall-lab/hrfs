/* standard use */
use std::io;
use std::io::prelude::*;
use std::str;

/* crate use */
use clap::Parser;
use gfa::{gfa::GFA, parser::GFAParser};
use handlegraph::{
    handle::{Edge, Handle},
    hashgraph::HashGraph,
    mutablehandlegraph::AdditiveHandleGraph,
    pathhandlegraph::{GraphPathNames, PathId},
};
use regex::Regex;
use rustc_hash::FxHashSet;

/* private use */
use founderset as ff;

#[derive(clap::Parser, Debug)]
#[clap(
    version = "0.1",
    author = "Daniel Doerr <daniel.doerr@hhu.de>",
    about = "Extract graph and path annotations surrounding a given path segment"
)]
pub struct Command {
    #[clap(index = 1, help = "graph in GFA1 format", required = true)]
    pub graph: String,
    #[clap(
        short = 'o',
        long = "only",
        help = "Only report subgraph that is induced by set of paths (\"P\" lines) whose names match given regular expression",
        default_value = ".*"
    )]
    pub paths_only: String,

    #[clap(
        short = 'l',
        long = "length",
        help = "Minimum length of any considered path",
        default_value = "1"
    )]
    pub min_length: usize,
}

fn identify_traversable_subgraph(
    graph: &HashGraph,
    re_path_names: &Regex,
    min_length: &usize,
) -> (
    FxHashSet<Handle>,
    FxHashSet<Handle>,
    FxHashSet<Edge>,
    Vec<PathId>,
) {
    let mut nodes: FxHashSet<Handle> = FxHashSet::default();
    let mut source_sinks: FxHashSet<Handle> = FxHashSet::default();
    let mut edges: FxHashSet<Edge> = FxHashSet::default();

    let mut paths: Vec<PathId> = Vec::new();
    for path_id in graph.paths.keys() {
        let path_name_vec = graph.get_path_name_vec(*path_id).unwrap();
        let path_name = str::from_utf8(&path_name_vec[..]).unwrap();
        if &graph.get_path(path_id).unwrap().len() >= min_length
            && re_path_names.is_match(&path_name)
        {
            let path = graph.get_path(path_id).unwrap();
            nodes.extend(path.nodes.iter().map(|x| x.forward()));
            for i in 0..path.nodes.len() - 1 {
                edges.insert(Edge(path.nodes[i], path.nodes[i + 1]));
                edges.insert(Edge(path.nodes[i + 1].flip(), path.nodes[i].flip()));
            }
            //
            // XXX detect the direction of the path; all paths must go in the same direction,
            // otherwise the LP will be at best infeasible and at worst return meaningless results
            //
            // direction is determined in a *very* naive way by counting the number of reversed
            // nodes along the path
            let no_rev_nodes = path
                .nodes
                .iter()
                .map(|x| if x.is_reverse() { 1 } else { 0 })
                .sum::<usize>();
            log::debug!(
                "path {} of length {} has {} reversed nodes",
                &path_name,
                &path.nodes.len(),
                &no_rev_nodes
            );
            let so = *path.nodes.first().unwrap();
            let si = *path.nodes.last().unwrap();
            log::debug!(
                "source and sink of path {} are {}{} and {}{}, respectively",
                &path_name,
                ff::v2str(&so),
                &so.unpack_number(),
                ff::v2str(&si),
                &si.unpack_number()
            );
            source_sinks.insert(so);
            source_sinks.insert(si.flip());
            paths.push(path_id.clone());
        }
    }
    (nodes, source_sinks, edges, paths)
}

fn add_source_sink(
    graph: &mut HashGraph,
    source_sinks: &FxHashSet<Handle>,
    subgraph_nodes: &mut FxHashSet<Handle>,
    subgraph_edges: &mut FxHashSet<Edge>,
) {
    let integrator = graph.append_handle(b"*");
    let source = graph.append_handle(b"*");
    let sink = graph.append_handle(b"*");

    for v in source_sinks.iter() {
        let e = Edge(integrator, *v);
        graph.create_edge(e);
        subgraph_edges.insert(e);
        subgraph_edges.insert(Edge(v.flip(), integrator.flip()));
    }

    graph.create_edge(Edge(source, integrator));
    graph.create_edge(Edge(sink.flip(), integrator));
    subgraph_edges.insert(Edge(source, integrator));
    subgraph_edges.insert(Edge(integrator.flip(), source.flip()));
    subgraph_edges.insert(Edge(sink.flip(), integrator));
    subgraph_edges.insert(Edge(integrator.flip(), sink));

    log::info!(
        "added source node >{}, sink node >{} and integrator node >{}",
        source.unpack_number(),
        sink.unpack_number(),
        integrator.unpack_number()
    );

    subgraph_nodes.insert(integrator);
    subgraph_nodes.insert(source);
    subgraph_nodes.insert(sink);
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
    let mut graph = HashGraph::from_gfa(&gfa);

    log::info!(
        "searching for paths matching regular expression \"{}\"",
        &params.paths_only
    );
    let re = Regex::new(&params.paths_only).unwrap();
    let (mut nodes, source_sinks, mut edges, paths) =
        identify_traversable_subgraph(&graph, &re, &params.min_length);

    log::info!(
        "identified {} source/sinks: {}",
        source_sinks.len(),
        ff::hv2seq(&source_sinks, ", ")
    );

    add_source_sink(&mut graph, &source_sinks, &mut nodes, &mut edges);

    log::info!("printing subgraph induced by path selection");
    ff::write_subset_gfa(&graph, &nodes, &edges, &paths, &mut out)?;
    out.flush()?;

    log::info!("done");
    Ok(())
}
