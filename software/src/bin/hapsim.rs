/* standard use */
use std::cmp::Reverse;
use std::collections::{BinaryHeap, HashSet};
use std::convert::TryInto;
use std::io;
use std::io::prelude::*;

/* crate use */
use clap::Parser;
use gfa::{gfa::GFA, parser::GFAParser};
use handlegraph::{
    handle::{Direction, Edge, Handle, NodeId},
    handlegraph::IntoNeighbors,
    hashgraph::HashGraph,
    mutablehandlegraph::AdditiveHandleGraph,
    pathhandlegraph::{embedded_paths::MutableGraphPaths, GraphPathsRef, IntoPathIds, PathId},
};
use rand::{
    distributions::{Bernoulli, Distribution, Uniform},
    rngs::ThreadRng,
    Rng,
};
use rustc_hash::FxHashSet;
use simd_adler32::Adler32;

/* private use */
use founderset as ff;

#[derive(clap::Parser, Debug)]
#[clap(
    author = "Konstantinn Bonnet <bonnetk@hhu.de>",
    about = "Simulate haplotype generation from a single founder sequence and output an annotated GFA1 file"
)]
struct Args {
    #[clap(
        short = 'f',
        long = "nfounder",
        default_value = "1",
        help = "Number of founder sequences"
    )]
    nfounder: usize,

    #[clap(
        short = 'l',
        long = "length",
        default_value = "1000",
        help = "Number of nodes in founder sequence (excl. duplicates, sources, sinks)"
    )]
    nnodes: usize,

    #[clap(
        short = 'n',
        long = "ndups",
        default_value = "50",
        help = "Number of duplications",
        conflicts_with = "dupratio"
    )]
    ndups: usize,

    #[clap(
        short = 'd',
        long = "dupratio",
        help = "Ratio of duplications",
        conflicts_with = "ndups"
    )]
    dupratio: Option<f64>,

    #[clap(
        short = 'r',
        long = "revdupratio",
        default_value = "0.5",
        help = "Ratio of reverse duplications in founders, [0;1]"
    )]
    revdupratio: f64,

    #[clap(
        short = 'G',
        long = "graphonly",
        conflicts_with = "infile",
        help = "Only generate graphs"
    )]
    nohap: bool,

    #[clap(
        short = 'i',
        long = "infile",
        default_value = "",
        help = "Read founder sequences from a GFA1 file instead of generating them"
    )]
    infile: String,

    nhaplotypes: usize,
}

fn hash_founders(g: &mut HashGraph, hashes: &mut HashSet<u32>) {
    let mut hash = Adler32::new();
    for p in &g.paths {
        hash.reset();
        let p = g.get_path(p.0).unwrap();
        hash.write(
            p.nodes
                .iter()
                .map(|v| {
                    format!(
                        "{}{}",
                        v.unpack_number(),
                        if v.is_reverse() { '-' } else { '+' }
                    )
                })
                .collect::<Vec<String>>()
                .join("")
                .as_bytes(),
        );
        hashes.insert(hash.finish());
    }
}

/* generate haplotype walks */
fn mk_haplotypes(
    g: &mut HashGraph,
    args: &Args,
    rng: &mut ThreadRng,
) -> (FxHashSet<Handle>, FxHashSet<Edge>) {
    log::info!("mk_haplotypes: nhaplotypes = {}", args.nhaplotypes);

    let mut ns: FxHashSet<Handle> = FxHashSet::default();
    let mut es: FxHashSet<Edge> = FxHashSet::default();
    let mut subns: FxHashSet<Handle> = FxHashSet::default();
    let mut subes: FxHashSet<Edge> = FxHashSet::default();

    let mut hashes = HashSet::new();
    let source = Handle::from(NodeId(0u64));
    let sink = Handle::from(NodeId((g.graph.len() - 1).try_into().unwrap()));

    /* consider founders as preexisting haplotypes to reject duplicates */
    hash_founders(g, &mut hashes);

    let mut hash = Adler32::new();
    let mut toterr = 0;
    let mut nerr = 0;
    let mut nh = 0;
    while nh < args.nhaplotypes {
        subns.clear();
        subes.clear();
        /* give up if we keep failing to generate new unique haplotypes */
        if nerr > 100 {
            log::debug!("giving up after {} errors in a row, {} total", nerr, toterr);
            break;
        }
        let s = "H".to_owned() + &nh.to_string();
        let p = g.create_path(s.as_bytes(), false).unwrap();
        hash.reset();

        let mut reject = false;
        let mut n = source;
        loop {
            g.path_append_step(p, n);
            let s = n.unpack_number().to_string() + (if n.is_reverse() { "-" } else { "+" });
            hash.write(s.as_bytes());
            if n.unpack_number() == sink.unpack_number() {
                break;
            }
            if n.unpack_number() == source.unpack_number() && n.is_reverse() {
                /* just reject the haplotype if we reached the source again
                 * this is the simplest solution, instead of trying to fix
                 * a random walk manually */
                reject = true;
                break;
            }

            let dir = if n.is_reverse() {
                Direction::Left
            } else {
                Direction::Right
            };
            let x = rng.gen_range(0..g.degree(n, dir));
            let m = g.neighbors(n, Direction::Right).nth(x).unwrap();
            subns.insert(n);
            subns.insert(m);
            subes.insert(Edge::edge_handle(n, m));
            n = m;
        }

        let h = hash.finish();
        if reject || hashes.contains(&h) {
            /* pathid not removed from known id's */
            let q = &g.get_path_ref(p).unwrap().name.clone();
            g.path_id.remove(q);
            g.destroy_path(p);
            nerr += 1;
            toterr += 1;
            if reject {
                log::debug!(
                    "rejecting haplotype returning to source; {} errors in a row, {} total",
                    nerr,
                    toterr
                );
            } else {
                log::debug!(
                    "rejecting duplicate haplotype; {} errors in a row, {} total",
                    nerr,
                    toterr
                );
            }
            continue;
        }
        assert!(!n.is_reverse());

        hashes.insert(h);
        ns.extend(&subns);
        es.extend(&subes);
        nerr = 0;
        nh += 1;
    }
    (ns, es)
}

fn mk_node(g: &mut HashGraph, id: usize) -> Handle {
    let h: Handle;
    let id = id.try_into().unwrap();
    match g.get_node(&NodeId(id)) {
        Some(_) => {
            h = Handle::from(NodeId(id));
        }
        None => {
            h = g.create_handle(b"*", id);
        }
    }
    h
}

/* create graph from generated founder sequences;
 * does not make assumptions on the number of sinks or sources,
 * the number of founders, or whether relabeling is necessary
 */
fn mk_graph(fnd: Vec<Vec<(usize, bool)>>) -> HashGraph {
    log::info!("mk_graph: using {} founder sequences", fnd.len());
    let mut g = HashGraph::new();
    for (i, f) in fnd.iter().enumerate() {
        let name = "F".to_owned() + &i.to_string();
        let p = g.create_path(name.as_bytes(), false).unwrap();
        let mut h = mk_node(&mut g, f[0].0);
        g.path_append_step(p, h);
        let mut hp = h;
        for (id, dir) in f.iter().skip(1) {
            h = mk_node(&mut g, *id);
            if !dir {
                h = h.flip();
            }
            g.create_edge(Edge(hp, h));
            g.path_append_step(p, h);
            hp = h;
            if id % 1000 == 0 {
                log::info!("processing node {} from founder {}", id, i);
            }
        }
    }
    g
}

/*
 * there's no need to handle multiple sources and sinks, one can just add
 * more artificial nodes to make sure that assumptions are respected:
 * source only has a single outgoing edge, sink only has a single
 * incoming edge.  the first nodeid must always be 0, the last nodeid
 * must always be length+1.
 * neither is not checked for!
 */
fn read_founders(args: &Args) -> HashGraph {
    let p = GFAParser::new();
    let gfa: GFA<usize, ()> = p.parse_file(&args.infile).unwrap();
    HashGraph::from_gfa(&gfa)
}

/* generate founder in O(nnodes + log ndups);  nodes are (index, direction),
 * where true = forward */
fn mk_one_founder(f: &mut Vec<(usize, bool)>, args: &Args, ndups: usize, rng: &mut ThreadRng) {
    log::info!(
        "mk_one_founder: generate founder with {} nodes, {} duplications with {} inversion ratio",
        args.nnodes,
        ndups,
        args.revdupratio
    );

    let mut dups = BinaryHeap::new();
    let mut pos = Vec::with_capacity(ndups);
    let ur = Uniform::from(1..args.nnodes);
    let br = Bernoulli::new(args.revdupratio).unwrap();

    /* sample duplicate positions in founder sequence as min-heap */
    for _ in 0..ndups {
        dups.push(Reverse(ur.sample(rng)));
    }
    /* sample duplicate nodes */
    for _ in 0..ndups {
        pos.push(ur.sample(rng));
    }
    /* source */
    f.push((0, true));
    /* founder sequence */
    let mut p = 1;
    let mut dp = 1;
    let mut draw = false;
    if let Some(Reverse(x)) = dups.pop() {
        dp = x;
        draw = true;
    }
    while p <= args.nnodes {
        f.push((p, true));
        while draw && dp == p {
            /* actual duplicate; avoid it linking source and make sure it
             * references existing node, otherwise we can forward and
             * induce more founders */
            let r = pos.pop().unwrap();
            let dir = br.sample(rng);
            f.push((r, !dir));
            if !dups.is_empty() {
                if let Some(Reverse(x)) = dups.pop() {
                    dp = x
                }
            } else {
                draw = false;
            }
        }
        p += 1;
        if p % 1000 == 0 {
            log::info!("processing node {} out of {}", p, args.nnodes);
        }
    }
    assert!(
        p == args.nnodes + 1 && dups.is_empty() && f.len() == args.nnodes + 1 + ndups,
        "bug: p {} dups {} dp {} f {} args.nnodes {}",
        p,
        dups.len(),
        dp,
        f.len(),
        args.nnodes
    );
    f.push((p, true));
}

fn mk_founders(args: &Args, ndups: usize, rng: &mut ThreadRng) -> Vec<Vec<(usize, bool)>> {
    log::info!(
        "mk_founders: generating {} founder sequences",
        args.nfounder
    );

    let mut f = vec![Vec::with_capacity(args.nnodes + 2 + ndups); args.nfounder];
    /* currently assumes only one founder is generated,
     * all others would be independent */
    for i in &mut f {
        mk_one_founder(i, args, ndups, rng);
    }
    f
}

fn main() -> Result<(), std::io::Error> {
    env_logger::init();
    let mut args = Args::parse();

    /* mk_graph also disregards multiple founders unlike the rest */
    if args.nfounder != 1 {
        log::warn!("only single founder generation is supported");
        args.nfounder = 1;
    }
    if args.revdupratio < 0.0 || args.revdupratio > 1.0 {
        panic!("invalid duplication inversion ratio");
    }
    let ndups = match args.dupratio {
        Some(d) => {
            if d < 0.0 {
                panic!("invalid duplication ratio");
            }
            (args.nnodes as f64 * d) as usize
        }
        None => args.ndups,
    };

    let mut rng = rand::thread_rng();
    let mut g;
    if args.infile.is_empty() {
        let fnd = mk_founders(&args, ndups, &mut rng);
        g = mk_graph(fnd);
    } else {
        g = read_founders(&args);
    }

    let mut out = io::BufWriter::new(std::io::stdout());
    if args.nhaplotypes < 1 || args.nohap {
        ff::write_gfa(&g, &mut out)?;
    } else {
        let (ns, es) = mk_haplotypes(&mut g, &args, &mut rng);
        ff::write_subset_gfa(
            &g,
            &ns,
            &es,
            &g.path_ids().collect::<Vec<PathId>>(),
            &mut out,
        )?;
    }
    out.flush()?;

    log::info!("done");
    Ok(())
}
