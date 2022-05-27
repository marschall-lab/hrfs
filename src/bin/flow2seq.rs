/* standard use */
use std::io::{self, Write};

/* crate use */
use clap::Parser;
use handlegraph::handle::Handle;
use itertools::Itertools;
use petgraph::{dot::Dot, graphmap::UnGraphMap};
use rustc_hash::{FxHashMap, FxHashSet};

/* private use */
use founderset as ff;

#[derive(clap::Parser, Debug)]
#[clap(
    version = "0.1",
    author = "Daniel Doerr <daniel.doerr@hhu.de>",
    about = "Transform network flow solution of founder sequence problem into a founder set"
)]
pub struct Command {
    #[clap(index = 1, help = "Solution file from Gurobi run", required = true)]
    pub flow_solution: String,
}

pub fn extract_linear_components(
    edges: &mut FxHashMap<ff::Extremity, FxHashSet<(ff::Extremity, usize)>>,
    flow: &ff::Flow,
) -> Vec<Vec<Handle>> {
    let mut res: Vec<Vec<Handle>> = Vec::new();

    log::info!("extracting linear components");
    for v in flow.sources.iter() {
        let w = flow.nodes.get(&(v.clone(), ff::Direction::In)).unwrap();
        log::info!("extracting {} linear components start at source {}", w, v);
        for _ in 0..*w {
            let s = ff::extract_random_walk_from_flow(
                edges,
                &ff::Extremity {
                    id: v.id,
                    etype: match v.etype {
                        ff::ExtremityType::Head => ff::ExtremityType::Tail,
                        ff::ExtremityType::Tail => ff::ExtremityType::Head,
                    },
                },
            );
            log::debug!("linear flow sequence: {:?}", ff::v2seq(&s, ""));
            res.push(s);
        }
    }
    log::info!("collected a total of {} linear compoonents", res.len());

    res
}

pub fn extract_circular_components(
    edges: &mut FxHashMap<ff::Extremity, FxHashSet<(ff::Extremity, usize)>>,
) -> Vec<Vec<Handle>> {
    let mut res: Vec<Vec<Handle>> = Vec::new();

    while !edges.is_empty() {
        let v = edges.keys().next().unwrap().clone();
        let w: usize = edges.get(&v).unwrap().iter().map(|x| x.1).sum();
        if w > 0 {
            log::info!("extracting circular component starting at {}", v);
            let c = ff::extract_random_walk_from_flow(edges, &v);
            log::debug!("circular flow sequence: {:?}", ff::v2seq(&c, ""));
            res.push(c);
        } else {
            edges.remove(&v);
        }
    }
    log::info!("collected a total of {} circular compoonents", res.len());

    res
}

pub fn construct_component_graph(components: &Vec<Vec<Handle>>) -> UnGraphMap<usize, Vec<u64>> {
    log::info!("constructing component graph");

    let mut graph: UnGraphMap<usize, Vec<u64>> = UnGraphMap::new();

    //                      node id ->component id
    let mut hits: FxHashMap<u64, Vec<usize>> = FxHashMap::default();

    for (i, c) in components.iter().enumerate() {
        for v in c.iter() {
            // record only first occurence of a node in a component
            hits.entry(v.unpack_number())
                .and_modify(|e| {
                    if e.last().unwrap() != &i {
                        e.push(i);
                    }
                })
                .or_insert(vec![i]);
        }
    }

    for (&vid, co_occurrences) in hits.iter() {
        for (&cid1, &cid2) in co_occurrences.iter().tuple_combinations() {
            if !graph.contains_node(cid1) {
                graph.add_node(cid1);
            }
            if !graph.contains_node(cid2) {
                graph.add_node(cid2);
            }
            match graph.edge_weight_mut(cid1, cid2) {
                None => {
                    graph.add_edge(cid1, cid2, vec![vid]);
                }
                Some(w) => w.push(vid),
            }
        }
    }

    log::debug!(
        "components: {}",
        components
            .iter()
            .map(|v| ff::v2seq(v, ""))
            .collect::<Vec<String>>()
            .join(", ")
    );
    log::debug!("component graph:\n\n{:?}", Dot::with_config(&graph, &[]));

    graph
}

fn build_founder_sequences(
    graph: &UnGraphMap<usize, Vec<u64>>,
    components: &Vec<Vec<Handle>>,
    n_linear: usize,
) -> Vec<Vec<Handle>> {
    log::info!("building founder sequences");

    let mut res: Vec<Vec<Handle>> = components[..n_linear].to_vec();

    let neighbors: Vec<Vec<usize>> = (0..components.len())
        .map(|v| graph.neighbors(v).collect())
        .collect();
    let mut queue: Vec<(usize, usize)> = (0..n_linear).zip(vec![0; n_linear]).collect();

    let mut merged_with: Vec<usize> = (0..n_linear).collect();
    merged_with.extend(vec![std::usize::MAX; components.len() - n_linear]);

    while !queue.is_empty() {
        let (u, e) = queue.pop().unwrap();

        if e < neighbors[u].len() {
            if merged_with[neighbors[u][e]] == std::usize::MAX {
                let v = neighbors[u][e];

                let mut donor = components[v].clone();
                let recipient = &mut res[merged_with[u]];

                // just pick the first one
                let recombination_node = graph.edge_weight(u, v).unwrap()[0];

                // prepare donor
                let mut i = 0;
                while donor[i].unpack_number() != recombination_node {
                    i += 1;
                }
                // remove last node -- which must be identical to the first node
                donor.pop();
                log::debug!(
                    "original sequence (recombination node has index {}): {}",
                    i,
                    ff::v2seq(&donor, "")
                );
                // rotate so that the recombination node will be the last element of the vector
                donor.rotate_left(i + 1);
                log::debug!(
                    "rotated sequence (recombination node is last): {}",
                    ff::v2seq(&donor, "")
                );

                // move cursor to recombination node of donor, where the recepient will be inserted
                let mut j = 0;
                while recipient[j].unpack_number() != recombination_node {
                    j += 1;
                }
                log::debug!(
                    "recombinations nodes: {} -- {}",
                    ff::v2str(&recipient[j]),
                    ff::v2str(donor.last().unwrap())
                );

                // if orientations of the donor and recipient handle do not agree, then reverse donor
                if recipient[j].is_reverse() != donor.last().unwrap().is_reverse() {
                    donor.reverse();
                    for x in 0..donor.len() {
                        donor[x] = donor[x].flip();
                    }
                    // move recombination node to the back again
                    donor.rotate_left(1);
                    log::debug!(
                        "reversed donor sequence (again recombination node is last): {}",
                        ff::v2seq(&donor, "")
                    );
                }
                for x in (0..donor.len()).rev() {
                    recipient.insert(j + 1, donor[x]);
                }
                log::debug!("merged sequence: {}", ff::v2seq(&recipient, ""));

                merged_with[v] = merged_with[u];
                queue.push((u, e + 1));
                queue.push((v, 0));
            } else {
                queue.push((u, e + 1));
            }
        }
    }

    res
}

fn main() -> Result<(), io::Error> {
    env_logger::init();

    // print output to stdout
    let mut out = io::BufWriter::new(io::stdout());

    // initialize command line parser & parse command line arguments
    let params = Command::parse();

    let f: ff::Flow = ff::read_flow(&params.flow_solution)?;
    f.log_sources_sinks();

    let mut edges = f.edges.clone();
    let mut components = extract_linear_components(&mut edges, &f);
    let linear_n = components.len();
    components.extend(extract_circular_components(&mut edges));
    let graph = construct_component_graph(&components);

    let founders = build_founder_sequences(&graph, &components, linear_n);
    ff::write_founders(&founders, &mut out)?;
    out.flush()?;

    log::info!("done");
    Ok(())
}
