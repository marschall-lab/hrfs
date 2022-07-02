/* standard use */
use std::cmp::min;
use std::fs;
use std::io::{self, BufRead, Write};
use std::iter::FromIterator;

/* crate use */
use clap::Parser;
use handlegraph::handle::Handle;
use itertools::Itertools;
use petgraph::{graphmap::DiGraphMap, Incoming, Outgoing};
use rustc_hash::{FxHashMap, FxHashSet};

/* private use */
use founderset as ff;

/* variables table type indices */
const VINT: usize = 0;
const VBIN: usize = VINT + 1;
const VNIL: usize = VBIN + 1;

#[derive(clap::Parser, Debug)]
#[clap(
    version = "0.1",
    author = "Konstantinn Bonnet <bonnetk@hhu.de>, Daniel Doerr <daniel.doerr@hhu.de>",
    about = "Generate ILP that minimizes recombinations between an (unoptimized) founder set and a set of haplotype sequences"
)]
pub struct Command {
    #[clap(index = 1, required = true, help = "founder sequences")]
    pub founder_set: String,

    #[clap(index = 2, required = true, help = "haplotype sequences")]
    pub haplotype_set: String,
}

/* https://www.gurobi.com/documentation/9.5/refman/lp_format.html
 * - constants must be on right hand side
 * - spaces are significant
 * - no constant variables
 * - no strict LT/GT: <, <= and resp. >, >= are strictly equivalent
 */

fn write_obj<W: io::Write>(
    color_conservation_vars: &FxHashSet<String>,
    out: &mut io::BufWriter<W>,
) -> Result<(), io::Error> {
    writeln!(
        out,
        "T - {}",
        color_conservation_vars
            .into_iter()
            .cloned()
            .collect::<Vec<String>>()
            .join(" - ")
    )
}

fn write_con_match<W: io::Write>(
    g: &DiGraphMap<ff::Node, ff::EdgeType>,
    vars: &mut [FxHashSet<String>],
    out: &mut io::BufWriter<W>,
) -> Result<(), io::Error> {
    writeln!(out, "\\\n\\ matching constraints for \"dashed\" edges\n\\")?;
    g.nodes().try_for_each(|u| {
        let mut y = Vec::new();
        let s = g
            .neighbors_directed(
                u,
                match u.direction {
                    ff::Direction::In => Outgoing,
                    ff::Direction::Out => Incoming,
                },
            )
            .map(|v| {
                let (xu, xv) = if u.direction == ff::Direction::In {
                    (u, v)
                } else {
                    (v, u)
                };
                let s = format!("x_{}_{}", xu, xv);
                let yu = format!("y_{}", xu);
                let yv = format!("y_{}", xv);
                vars[VBIN].insert(s.clone());
                vars[VBIN].insert(yu.clone());
                vars[VBIN].insert(yv.clone());
                y.push(format!("2 {} - {} - {} <= 0", s, yu, yv));
                s
            })
            .collect::<Vec<String>>()
            .join(" + ");
        for c in &y {
            writeln!(out, "{}", c)?;
        }
        writeln!(out, "{} - y_{} = 0", s, u)
    })
}

fn write_con_flow_solid<W: io::Write>(
    g: &DiGraphMap<ff::Node, ff::EdgeType>,
    vars: &mut [FxHashSet<String>],
    out: &mut io::BufWriter<W>,
) -> Result<(), io::Error> {
    writeln!(out, "\\\n\\ flow constraints for \"solid\" edges\n\\")?;
    g.all_edges()
        .filter(|(_, _, ref t)| **t == ff::EdgeType::Solid)
        .try_for_each(|(u, v, _)| {
            vars[VINT].insert(format!("f_{}", u));
            vars[VINT].insert(format!("f_{}", v));
            writeln!(out, "f_{} - f_{}  = 0", v, u)
        })
}

fn write_solid_edges<W: io::Write>(
    g: &DiGraphMap<ff::Node, ff::EdgeType>,
    vars: &mut [FxHashSet<String>],
    out: &mut io::BufWriter<W>,
) -> Result<(), io::Error> {
    g.all_edges()
        .filter(|(_, _, ref t)| **t == ff::EdgeType::Solid)
        .try_for_each(|(u, v, _)| {
            vars[VBIN].insert(format!("x_{}_{}", u, v));
            writeln!(out, "y_{} - x_{}_{} = 0", u, u, v)?;
            writeln!(out, "y_{} - x_{}_{} = 0", v, u, v)
        })
}

fn write_con_flow_dashed<W: io::Write>(
    g: &DiGraphMap<ff::Node, ff::EdgeType>,
    totflow: usize,
    vars: &mut [FxHashSet<String>],
    out: &mut io::BufWriter<W>,
) -> Result<(), io::Error> {
    writeln!(out, "\\\n\\ flow constraints for \"dashed\" edges\n\\")?;
    g.all_edges()
        .filter(|(_, _, ref t)| **t == ff::EdgeType::Dashed)
        .try_for_each(|(u, v, _)| {
            vars[VINT].insert(format!("f_{}", u));
            vars[VINT].insert(format!("f_{}", v));
            vars[VBIN].insert(format!("x_{}_{}", u, v));
            writeln!(
                out,
                "f_{} - f_{} + {} x_{}_{} <= {}",
                v,
                u,
                totflow,
                u,
                v,
                totflow + 1
            )?;
            writeln!(
                out,
                "f_{} - f_{} + {} x_{}_{} <= {}",
                u,
                v,
                totflow,
                u,
                v,
                totflow - 1
            )
        })
}

fn write_con_flow_source<W: io::Write>(
    sources: &FxHashSet<ff::Node>,
    vars: &mut [FxHashSet<String>],
    out: &mut io::BufWriter<W>,
) -> Result<(), io::Error> {
    writeln!(out, "\\\n\\ flow constraints for \"source\" nodes\n\\")?;
    sources.iter().try_for_each(|u| {
        vars[VINT].insert(format!("f_{}", u));
        writeln!(out, "f_{} = 0", u)
    })
}

fn write_con_color_contiguity<W: io::Write>(
    graph: &DiGraphMap<ff::Node, ff::EdgeType>,
    sources: &FxHashSet<ff::Node>,
    haplotype: &Vec<Handle>,
    haplotype_id: usize,
    color_vars: &mut FxHashMap<ff::Node, FxHashSet<(usize, usize)>>,
    color_conservation_vars: &mut FxHashSet<String>,
    vars: &mut [FxHashSet<String>],
    pick_up_nodes: bool,
    out: &mut io::BufWriter<W>,
) -> Result<(), io::Error> {
    // by construction, cur_nodes always maintain a set of incoming nodes
    let mut cur_nodes: FxHashSet<ff::Node> = FxHashSet::default();
    cur_nodes.extend(sources.iter().cloned());

    let mut handle2in: FxHashMap<Handle, Vec<ff::Node>> = FxHashMap::default();
    for v in graph.nodes() {
        if v.direction == ff::Direction::In {
            let vh = Handle::pack(v.node, v.etype == ff::ExtremityType::Head);
            handle2in
                .entry(vh)
                .and_modify(|x| x.push(v))
                .or_insert(vec![v]);
        }
    }

    writeln!(
        out,
        "\\\n\\ color contiguity constraints for haplotype {} \n\\",
        haplotype_id
    )?;

    for ((i, uu), (j, vv)) in haplotype.iter().enumerate().tuple_windows() {
        if pick_up_nodes && cur_nodes.is_empty() {
            cur_nodes.extend(handle2in.get(uu).or(Some(&Vec::new())).unwrap());
            log::info!("picking up {} flow nodes at {}:{} of haplotype {} for writing color contiguity constraints", cur_nodes.len(), ff::v2str(uu), i, haplotype_id);
        }
        let mut next_nodes = FxHashSet::default();
        for u in cur_nodes.iter() {
            assert!(
                u.direction == ff::Direction::In
                    && u.node == uu.unpack_number()
                    && u.etype
                        == if uu.is_reverse() {
                            ff::ExtremityType::Head
                        } else {
                            ff::ExtremityType::Tail
                        },
                "variation graph node {} does not correspond to flow graph node {}",
                ff::v2str(uu),
                u
            );
            for w in graph.neighbors(*u) {
                assert!(graph.neighbors(w).count() == 1,
                    "flow graph outgoing inner node {} does not have exactly one adjacency but {}: {}",
                    w, graph.neighbors(w).count(), graph.neighbors(w).join(", "));
                let v = graph.neighbors(w).next().unwrap();
                let v_etype = if vv.is_reverse() {
                    ff::ExtremityType::Head
                } else {
                    ff::ExtremityType::Tail
                };
                if v.node == vv.unpack_number() && v.etype == v_etype {
                    let t = format!("t_{}_{}_{}_{}", u, w, haplotype_id, i);
                    let x = format!("x_{}_{}", u, w);
                    let c1 = format!("c_{}_{}_{}", u, haplotype_id, i);
                    let c2 = format!("c_{}_{}_{}", w, haplotype_id, i);
                    let c3 = format!("c_{}_{}_{}", v, haplotype_id, j);
                    writeln!(out, "3 {} - {} - {} - {} <= 0", t, x, c1, c2)?;
                    writeln!(out, "{} - {} = 0", c2, c3)?;
                    vars[VBIN].insert(t.clone());
                    vars[VBIN].insert(c1);
                    vars[VBIN].insert(c2);
                    color_vars
                        .entry(*u)
                        .or_insert(FxHashSet::default())
                        .insert((haplotype_id, i));
                    color_vars
                        .entry(w)
                        .or_insert(FxHashSet::default())
                        .insert((haplotype_id, i));
                    color_vars
                        .entry(v)
                        .or_insert(FxHashSet::default())
                        .insert((haplotype_id, j));
                    color_conservation_vars.insert(t);
                    next_nodes.insert(v);
                }
            }
        }
        cur_nodes = next_nodes;
    }

    let i = haplotype.len() - 1;
    for v in cur_nodes.iter() {
        for w in graph.neighbors(*v) {
            let t = format!("t_{}_{}_{}_{}", v, w, haplotype_id, i);
            let x = format!("x_{}_{}", v, w);
            let c1 = format!("c_{}_{}_{}", v, haplotype_id, i);
            let c2 = format!("c_{}_{}_{}", w, haplotype_id, i);
            writeln!(out, "3 {} - {} - {} - {} <= 0", t, x, c1, c2)?;
            vars[VBIN].insert(t.clone());
            vars[VBIN].insert(c1);
            vars[VBIN].insert(c2);
            color_vars
                .entry(*v)
                .or_insert(FxHashSet::default())
                .insert((haplotype_id, i));
            color_vars
                .entry(w)
                .or_insert(FxHashSet::default())
                .insert((haplotype_id, i));
            color_conservation_vars.insert(t);
        }
    }

    Ok(())
}

fn write_con_color_singularity<W: io::Write>(
    color_vars: &FxHashMap<ff::Node, FxHashSet<(usize, usize)>>,
    out: &mut io::BufWriter<W>,
) -> Result<(), io::Error> {
    writeln!(out, "\\\n\\ color singularity constraints \n\\")?;

    for (v, cs) in color_vars.iter() {
        writeln!(
            out,
            "{} = 1",
            cs.into_iter()
                .map(|(h, i)| format!("c_{}_{}_{}", v, h, i))
                .join(" + ")
        )?;
    }

    Ok(())
}

fn write_con_flow_matching<W: io::Write>(
    g: &DiGraphMap<ff::Node, ff::EdgeType>,
    flowmap: &FxHashMap<Handle, usize>,
    out: &mut io::BufWriter<W>,
) -> Result<(), io::Error> {
    writeln!(out, "\\\n\\ flow node matching constraints \n\\")?;

    // intermediary data structure to access nodes of extremities in a convenient way
    let mut nodes: FxHashMap<(u64, ff::ExtremityType), Vec<ff::Node>> = FxHashMap::default();
    g.nodes()
        .for_each(|v| nodes.entry((v.node, v.etype)).or_insert(Vec::new()).push(v));

    // intermediary data structure that reports the total flow for each handle, no matter the flow
    // orientation
    let mut totflow: FxHashMap<Handle, usize> = FxHashMap::default();
    flowmap.iter().for_each(|(v, f)| {
        totflow
            .entry(v.forward())
            .and_modify(|x| *x += f)
            .or_insert(*f);
    });

    for (v, f) in totflow.iter() {
        writeln!(
            out,
            "{} = {}",
            nodes
                .get(&(v.unpack_number(), ff::ExtremityType::Tail))
                .unwrap()
                .iter()
                .map(|u| format!("y_{}", u))
                .join(" + "),
            f
        )?;
        writeln!(
            out,
            "{} = {}",
            nodes
                .get(&(v.unpack_number(), ff::ExtremityType::Head))
                .unwrap()
                .iter()
                .map(|u| format!("y_{}", u))
                .join(" + "),
            f
        )?;
    }

    Ok(())
}

fn write_con_flow_adj_matching<W: io::Write>(
    g: &DiGraphMap<ff::Node, ff::EdgeType>,
    flowmap: &FxHashMap<(Handle, Handle), usize>,
    out: &mut io::BufWriter<W>,
) -> Result<(), io::Error> {
    writeln!(out, "\\\n\\ flow adjacency matching constraints \n\\")?;

    // intermediary data structure that reports the total flow for adjacency, no matter the flow
    // orientation
    let mut total_flow: FxHashMap<(Handle, Handle), usize> = FxHashMap::default();
    flowmap.into_iter().for_each(|(e, c)| {
        total_flow
            .entry(normalize(e))
            .and_modify(|x| *x += c)
            .or_insert(*c);
    });

    // intermediary data structure to access adjacencies in a convenient way
    let mut adjs: FxHashMap<(Handle, Handle), Vec<(ff::Node, ff::Node)>> = FxHashMap::default();
    g.all_edges()
        .filter(|(_, _, ref t)| **t == ff::EdgeType::Solid)
        .for_each(|(u, v, _)| {
            adjs.entry(normalize(&(
                Handle::pack(u.node, u.etype == ff::ExtremityType::Tail),
                Handle::pack(v.node, v.etype == ff::ExtremityType::Head),
            )))
            .or_insert(Vec::new())
            .push((u, v));
        });

    for (e, f) in total_flow.iter() {
        writeln!(
            out,
            "{} = {}",
            adjs.get(e)
                .unwrap()
                .iter()
                .map(|(u, v)| format!("x_{}_{}", u, v))
                .join(" + "),
            f
        )?;
    }

    Ok(())
}

fn write_bounds<W: io::Write>(
    g: &DiGraphMap<ff::Node, ff::EdgeType>,
    totflow: usize,
    out: &mut io::BufWriter<W>,
) -> Result<(), io::Error> {
    writeln!(out, "Bounds")?;
    g.nodes()
        .try_for_each(|u| writeln!(out, "0 <= f_{} <= {}", u, totflow))
}

fn write_vars<W: io::Write>(
    vars: &[FxHashSet<String>],
    out: &mut io::BufWriter<W>,
) -> Result<(), io::Error> {
    if vars[VINT].len() > 0 {
        writeln!(out, "Generals")?;
        vars[VINT].iter().try_for_each(|v| writeln!(out, "{}", v))?;
    }
    if vars[VBIN].len() > 0 {
        writeln!(out, "Binary")?;
        vars[VBIN].iter().try_for_each(|v| writeln!(out, "{}", v))?;
    }
    Ok(())
}

fn write_lp<W: io::Write>(
    g: &DiGraphMap<ff::Node, ff::EdgeType>,
    haplotypes: &Vec<(String, Vec<Handle>)>,
    flowmap: &FxHashMap<(Handle, Handle), usize>,
    out: &mut io::BufWriter<W>,
) -> Result<(), io::Error> {
    let mut vars: [FxHashSet<String>; VNIL] = [FxHashSet::default(), FxHashSet::default()];
    let mut tmp = io::BufWriter::new(Vec::new());
    let mut color_vars: FxHashMap<ff::Node, FxHashSet<(usize, usize)>> = FxHashMap::default();
    let mut color_conservation_vars: FxHashSet<String> = FxHashSet::default();
    let mut node_multimap = map_node_multiplicity(flowmap);

    let (src, snk) = ff::term_nodes(&g);
    log::info!(
        "identified {} source nodes: {}",
        src.len(),
        src.iter()
            .map(|x| format!("{}", x))
            .collect::<Vec<String>>()
            .join(", ")
    );
    log::info!(
        "identified {} sink nodes: {}",
        snk.len(),
        snk.iter()
            .map(|x| format!("{}", x))
            .collect::<Vec<String>>()
            .join(", ")
    );

    snk.iter().for_each(|x| {
        node_multimap
            .entry(Handle::pack(x.node, x.etype == ff::ExtremityType::Tail))
            .and_modify(|y| *y += 1)
            .or_insert(1);
    });
    log::debug!(
        "node multiplicity map: {}",
        node_multimap
            .iter()
            .map(|(u, c)| format!("{}:{}", ff::v2str(u), c))
            .collect::<Vec<String>>()
            .join(", ")
    );

    let totflow: usize = node_multimap.values().sum::<usize>();
    for (i, (_, hap)) in haplotypes.iter().enumerate() {
        write_con_color_contiguity(
            g,
            &src,
            hap,
            i + 1,
            &mut color_vars,
            &mut color_conservation_vars,
            &mut vars,
            false,
            &mut tmp,
        )?;

        // sources should be the outgoing nodes of the sinks; but since the outgoing flow of the
        // sinks is by definition 0, such nodes cannot exist, so the list of sources is empty
        write_con_color_contiguity(
            g,
            &FxHashSet::default(),
            &ff::reverse_seq(hap),
            i + 1,
            &mut color_vars,
            &mut color_conservation_vars,
            &mut vars,
            true,
            &mut tmp,
        )?;
    }
    vars[VINT].insert("T".to_string());

    // write objective function
    writeln!(out, "Minimize")?;
    write_obj(&color_conservation_vars, out)?;
    writeln!(out, "Subject To")?;
    writeln!(out, "T = {}", totflow)?;
    // write constraints
    write_con_match(g, &mut vars, out)?;
    write_con_flow_matching(g, &node_multimap, out)?;
    write_con_flow_adj_matching(g, flowmap, out)?;
    write_solid_edges(g, &mut vars, out)?; // fixed, used in min2seq
    write_con_flow_solid(g, &mut vars, out)?;
    write_con_flow_dashed(g, totflow, &mut vars, out)?;
    write_con_flow_source(&src, &mut vars, out)?;

    out.write(&tmp.into_inner()?)?;
    write_con_color_singularity(&color_vars, out)?;

    // write bounds and collected variables
    write_bounds(g, totflow, out)?;
    write_vars(&vars, out)?;
    writeln!(out, "End")?;
    Ok(())
}

fn construct_graph_from_adj_multiplicities(
    multiplicities: &FxHashMap<(Handle, Handle), usize>,
) -> DiGraphMap<ff::Node, ff::EdgeType> {
    let mut count: FxHashMap<(Handle, ff::Direction), usize> = FxHashMap::default();

    let mut graph: DiGraphMap<ff::Node, ff::EdgeType> = DiGraphMap::new();

    for ((u, v), &c) in multiplicities.into_iter() {
        for _ in 0..c {
            let uo = ff::Node {
                node: u.unpack_number(),
                direction: ff::Direction::Out,
                etype: {
                    if u.is_reverse() {
                        ff::ExtremityType::Tail
                    } else {
                        ff::ExtremityType::Head
                    }
                },
                id: *count
                    .entry((*u, ff::Direction::Out))
                    .and_modify(|x| *x += 1)
                    .or_insert(1)
                    - 1,
            };
            let vi = ff::Node {
                node: v.unpack_number(),
                direction: ff::Direction::In,
                etype: {
                    if v.is_reverse() {
                        ff::ExtremityType::Head
                    } else {
                        ff::ExtremityType::Tail
                    }
                },
                id: *count
                    .entry((*v, ff::Direction::In))
                    .and_modify(|x| *x += 1)
                    .or_insert(1)
                    - 1,
            };
            graph.add_node(uo);
            graph.add_node(vi);
            graph.add_edge(uo, vi, ff::EdgeType::Solid);
        }
    }

    let sources: FxHashSet<(Handle, usize)> =
        FxHashSet::from_iter(count.iter().filter_map(|((v, d), c)| {
            if *d == ff::Direction::Out && !count.contains_key(&(*v, ff::Direction::In)) {
                Some((*v, *c))
            } else {
                None
            }
        }));
    log::info!(
        "identified sources: {}",
        sources.iter().map(|(v, _)| ff::v2str(v)).join(", ")
    );

    sources.iter().for_each(|(v, c)| {
        (0..*c).for_each(|i| {
            graph.add_node(ff::Node {
                node: v.unpack_number(),
                direction: ff::Direction::In,
                etype: if v.is_reverse() {
                    ff::ExtremityType::Head
                } else {
                    ff::ExtremityType::Tail
                },
                id: i,
            });
        });
        count.insert((v.clone(), ff::Direction::In), *c);
    });

    let sinks: FxHashSet<(Handle, usize)> =
        FxHashSet::from_iter(count.iter().filter_map(|((v, d), c)| {
            if *d == ff::Direction::In && !count.contains_key(&(*v, ff::Direction::Out)) {
                Some((*v, *c))
            } else {
                None
            }
        }));
    log::info!(
        "identified sinks: {}",
        sinks.iter().map(|(v, _)| ff::v2str(v)).join(", ")
    );

    sinks.iter().for_each(|(v, c)| {
        (0..*c).for_each(|i| {
            graph.add_node(ff::Node {
                node: v.unpack_number(),
                direction: ff::Direction::Out,
                etype: if v.is_reverse() {
                    ff::ExtremityType::Tail
                } else {
                    ff::ExtremityType::Head
                },
                id: i,
            });
        });
        count.insert((v.clone(), ff::Direction::Out), *c);
    });

    // construct dashed edges
    for ((v, d), ci) in count.iter() {
        if *d == ff::Direction::In {
            let co = count.get(&(*v, ff::Direction::Out)).expect(&format!(
                "count map does not contain entry {}::Out",
                ff::v2str(v)
            ));
            // construct temporary data structure to prohibit switching between nodes whose other ends
            // point the exact same handle
            let ext_l = if v.is_reverse() {
                ff::ExtremityType::Head
            } else {
                ff::ExtremityType::Tail
            };
            let ext_r = if v.is_reverse() {
                ff::ExtremityType::Tail
            } else {
                ff::ExtremityType::Head
            };
            let vid = v.unpack_number();

            let orig: Vec<(u64, ff::ExtremityType)> = (0..*ci)
                .map(|i| {
                    let vr = ff::Node {
                        node: vid,
                        direction: ff::Direction::In,
                        etype: ext_l,
                        id: i,
                    };
                    match graph.neighbors_directed(vr, Incoming).next() {
                        None => (u64::MAX, ff::ExtremityType::Tail),
                        Some(ul) => (ul.node, ul.etype),
                    }
                })
                .collect();
            let mut orig_count: FxHashMap<(u64, ff::ExtremityType), usize> = FxHashMap::default();
            orig.iter().for_each(|z| {
                orig_count.entry(*z).and_modify(|x| *x += 1).or_insert(1);
            });

            let dest: Vec<(u64, ff::ExtremityType)> = (0..*co)
                .map(|i| {
                    let vr = ff::Node {
                        node: vid,
                        direction: ff::Direction::Out,
                        etype: ext_r,
                        id: i,
                    };
                    match graph.neighbors_directed(vr, Outgoing).next() {
                        None => (u64::MAX, ff::ExtremityType::Tail),
                        Some(ul) => (ul.node, ul.etype),
                    }
                })
                .collect();
            let mut dest_count: FxHashMap<(u64, ff::ExtremityType), usize> = FxHashMap::default();
            dest.iter().for_each(|z| {
                dest_count.entry(*z).and_modify(|x| *x += 1).or_insert(1);
            });

            let mut dest_idx: FxHashMap<
                ((u64, ff::ExtremityType), (u64, ff::ExtremityType)),
                (usize, usize),
            > = FxHashMap::default();

            for (i, j) in (0..*ci).cartesian_product(0..*co) {
                if orig_count[&orig[i]] != dest_count[&dest[j]]
                    || !dest_idx.contains_key(&(orig[i], dest[j]))
                    || {
                        let (ii, jj) = dest_idx.get(&(orig[i], dest[j])).unwrap();
                        *ii < i && *jj < j
                    }
                {
                    let vl = ff::Node {
                        node: vid,
                        direction: ff::Direction::In,
                        etype: ext_l,
                        id: i,
                    };
                    let vr = ff::Node {
                        node: vid,
                        direction: ff::Direction::Out,
                        etype: ext_r,
                        id: j,
                    };
                    dest_idx.insert((orig[i], dest[j]), (i, j));
                    graph.add_edge(vl, vr, ff::EdgeType::Dashed);
                }
            }
        }
    }

    //    log::debug!("graph: \n{}", Dot::new(&graph));
    graph
}

fn read_founderseq_adj_multiplicity<R: io::Read>(
    data: io::BufReader<R>,
) -> Result<FxHashMap<(Handle, Handle), usize>, String> {
    let mut multiplicity: FxHashMap<(Handle, Handle), usize> = FxHashMap::default();

    for line_op in data.lines() {
        if let Ok(line) = line_op {
            let walk = ff::parse_walk(&line)?;
            // insert first element of walk
            walk.into_iter().tuple_windows().for_each(|a| {
                multiplicity.entry(a).and_modify(|x| *x += 1).or_insert(1);
            });
        }
    }
    Ok(multiplicity)
}

fn map_node_multiplicity(fmap: &FxHashMap<(Handle, Handle), usize>) -> FxHashMap<Handle, usize> {
    let mut mmap: FxHashMap<Handle, usize> = FxHashMap::default();
    fmap.iter().for_each(|((u, _), k)| {
        mmap.entry(*u).and_modify(|x| *x += k).or_insert(*k);
    });
    mmap
}

fn map_haplotype_adj_multiplicity(
    haps: &Vec<(String, Vec<Handle>)>,
) -> FxHashMap<(Handle, Handle), usize> {
    let mut res: FxHashMap<(Handle, Handle), usize> = FxHashMap::default();
    for (_, s) in haps.into_iter() {
        s.into_iter().tuple_windows().for_each(|(u, v)| {
            res.entry((*u, *v)).and_modify(|x| *x += 1).or_insert(1);
        })
    }
    res
}

fn normalize(e: &(Handle, Handle)) -> (Handle, Handle) {
    if (e.0.unpack_number() == e.1.unpack_number() && e.0 == e.1 && !e.0.is_reverse())
        || e.0.unpack_number() < e.1.unpack_number()
    {
        (e.0.clone(), e.1.clone())
    } else {
        (e.1.flip(), e.0.flip())
    }
}

fn merge_adj_multiplicities(
    flowmap: &FxHashMap<(Handle, Handle), usize>,
    hapmap: &FxHashMap<(Handle, Handle), usize>,
) -> FxHashMap<(Handle, Handle), usize> {
    log::info!("merging adjacency multiplicities");

    let mut total_flow: FxHashMap<(Handle, Handle), usize> = FxHashMap::default();
    flowmap.into_iter().for_each(|(e, c)| {
        total_flow
            .entry(normalize(e))
            .and_modify(|x| *x += c)
            .or_insert(*c);
    });
    log::debug!(
        "total flow multiplicity map: {}",
        total_flow
            .iter()
            .map(|((u, v), c)| format!("{}{}:{}", ff::v2str(u), ff::v2str(v), c))
            .collect::<Vec<String>>()
            .join(", ")
    );

    let mut merged: FxHashMap<(Handle, Handle), usize> = FxHashMap::default();
    FxHashSet::from_iter(flowmap.keys().chain(hapmap.keys()))
        .into_iter()
        .for_each(|&e| {
            merged.insert(
                e,
                std::cmp::max(
                    *flowmap.get(&e).or(Some(&0)).unwrap(),
                    min(
                        *total_flow.get(&normalize(&e)).expect(&format!(
                            "edge {}{} not contained in total flow map",
                            ff::v2str(&e.0),
                            ff::v2str(&e.1)
                        )),
                        *hapmap.get(&e).or(Some(&0)).unwrap(),
                    ),
                ),
            );
        });
    merged
}

fn read_haplotypes<R: io::Read>(
    data: io::BufReader<R>,
) -> Result<Vec<(String, Vec<Handle>)>, String> {
    let mut res: Vec<(String, Vec<Handle>)> = Vec::new();

    for line_op in data.lines() {
        if let Ok(line) = line_op {
            res.push((
                line[..line.find('\t').unwrap()].to_string(),
                ff::parse_walk(&line)?,
            ));
        }
    }

    Ok(res)
}

fn main() -> Result<(), std::io::Error> {
    env_logger::init();
    // initialize command line parser & parse command line arguments
    let params = Command::parse();

    log::info!(
        "reading adjacency multiplicities from founder set {}",
        params.founder_set
    );
    let founder_data = io::BufReader::new(fs::File::open(&params.founder_set)?);
    let flow_multi = read_founderseq_adj_multiplicity(founder_data).unwrap();
    log::debug!(
        "flow multiplicity map: {}",
        flow_multi
            .iter()
            .map(|((u, v), c)| format!("{}{}:{}", ff::v2str(u), ff::v2str(v), c))
            .collect::<Vec<String>>()
            .join(", ")
    );
    //    let totflow: usize = g.node_count() / 2;
    //    log::info!("flow graph total flow: {}", totflow);

    log::info!(
        "reading adjacency multiplicities from haplotype set {}",
        params.haplotype_set
    );
    let hap_data = io::BufReader::new(fs::File::open(&params.haplotype_set)?);
    let haplotypes = read_haplotypes(hap_data).unwrap();
    let hap_multimap = map_haplotype_adj_multiplicity(&haplotypes);
    log::debug!(
        "haplotype multiplicity map: {}",
        hap_multimap
            .iter()
            .map(|((u, v), c)| format!("{}{}:{}", ff::v2str(u), ff::v2str(v), c))
            .collect::<Vec<String>>()
            .join(", ")
    );

    log::info!("merging multiplicities between the two sets");
    let merged_multimap = merge_adj_multiplicities(&flow_multi, &hap_multimap);
    log::debug!(
        "merged multiplicity map: {}",
        merged_multimap
            .iter()
            .map(|((u, v), c)| format!("{}{}:{}", ff::v2str(u), ff::v2str(v), c))
            .collect::<Vec<String>>()
            .join(", ")
    );

    log::info!("constructing graph from multiplicity map");
    let g = construct_graph_from_adj_multiplicities(&merged_multimap);

    //    log::info!(
    //        "colorizing graph based on haplotype sequences of {}",
    //        params.haplotype_set
    //    );
    // intermediate data structure that maps edges of the variation graph to edges in the flow
    // graph
    log::debug!(
        "total multiplicity flow multimap: {}",
        flow_multi.values().sum::<usize>()
    );
    log::debug!(
        "total multiplicity in hap multimap: {}",
        hap_multimap.values().sum::<usize>()
    );
    log::debug!(
        "total multiplicity in merged multimap: {}",
        merged_multimap.values().sum::<usize>()
    );

    let mut out = io::BufWriter::new(std::io::stdout());
    write_lp(&g, &haplotypes, &flow_multi, &mut out)?;
    out.flush()?;

    log::info!("done");
    Ok(())
}
