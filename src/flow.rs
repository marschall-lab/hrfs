/* crate use */
use handlegraph::handle::Handle;
use quick_csv::Csv;
use regex::Regex;
use rustc_hash::{FxHashMap, FxHashSet};

/* standard use */
use std::fs;
use std::io;
use std::str::{self, FromStr};

/* private use */
use crate::{graph::*, *};

pub struct Flow {
    pub nodes: FxHashMap<(Extremity, Direction), usize>,
    // directed edges from out -> in
    pub edges: FxHashMap<Extremity, FxHashSet<(Extremity, usize)>>,
    pub sources: FxHashSet<Extremity>,
    pub sinks: FxHashSet<Extremity>,
}

impl Flow {
    fn identify_sources_sinks(
        edges: &FxHashMap<Extremity, FxHashSet<(Extremity, usize)>>,
    ) -> (FxHashSet<Extremity>, FxHashSet<Extremity>) {
        let mut sinks: FxHashSet<Extremity> = FxHashSet::default();
        let mut sources: FxHashSet<Extremity> = FxHashSet::default();
        let mut incoming: FxHashMap<Extremity, bool> = FxHashMap::default();

        for (v, adjacent_nodes) in edges.iter() {
            let mut has_nonzero_weight = false;
            for (u, weight) in adjacent_nodes.iter() {
                let w = Extremity {
                    id: u.id,
                    etype: match u.etype {
                        ExtremityType::Head => ExtremityType::Tail,
                        ExtremityType::Tail => ExtremityType::Head,
                    },
                };
                if weight > &0 {
                    has_nonzero_weight = true;
                    if match edges.get(&w) {
                        None => true,
                        Some(es) => es.len() == 0,
                    } {
                        sinks.insert(w.clone());
                    }
                    incoming.insert(u.clone(), true);
                }
            }

            let u = Extremity {
                id: v.id,
                etype: match v.etype {
                    ExtremityType::Head => ExtremityType::Tail,
                    ExtremityType::Tail => ExtremityType::Head,
                },
            };
            if has_nonzero_weight {
                incoming.entry(u).or_insert(false);
            }
        }

        sources.extend(incoming.iter().filter_map(|(v, &has_incoming)| {
            if !has_incoming {
                Some(v.clone())
            } else {
                None
            }
        }));

        (sources, sinks)
    }

    pub fn log_sources_sinks(&self) {
        log::info!(
            "identified sources: {}\nidentified sinks: {}",
            self.sources
                .iter()
                .map(|v| v.to_string())
                .collect::<Vec<String>>()
                .join(","),
            self.sinks
                .iter()
                .map(|v| v.to_string())
                .collect::<Vec<String>>()
                .join(",")
        )
    }

    pub fn new(
        nodes: FxHashMap<(Extremity, Direction), usize>,
        edges: FxHashMap<Extremity, FxHashSet<(Extremity, usize)>>,
    ) -> Flow {
        let (sources, sinks) = Flow::identify_sources_sinks(&edges);
        Flow {
            nodes,
            edges,
            sources,
            sinks,
        }
    }
}

pub fn read_flow(file: &String) -> Result<Flow, io::Error> {
    let mut nodes: FxHashMap<(Extremity, Direction), usize> = FxHashMap::default();
    let mut edges: FxHashMap<Extremity, FxHashSet<(Extremity, usize)>> = FxHashMap::default();

    let pat_node = Regex::new(r"^(i|o)(\d+)(h|t)$").unwrap();
    let pat_edge = Regex::new(r"^(i|o)(\d+)(h|t)_(i|o)(\d+)(h|t)$").unwrap();

    log::info!("loading flow solution {}", &file);
    let mut bf = io::BufReader::new(fs::File::open(&file)?);
    let reader = Csv::from_reader(&mut bf)
        .delimiter(b' ')
        .flexible(true)
        .has_header(true);

    for row in reader.into_iter() {
        let row = row.unwrap();
        let mut row_it = row.bytes_columns();
        let var = str::from_utf8(row_it.next().unwrap()).unwrap();
        let val = usize::from_str(str::from_utf8(row_it.next().unwrap()).unwrap()).unwrap();

        if let Some(m) = pat_edge.captures(var) {
            let u = Extremity {
                id: usize::from_str(&m[2]).unwrap(),
                etype: match &m[3] {
                    "h" => ExtremityType::Head,
                    _ => ExtremityType::Tail,
                },
            };
            assert!(m[1] != m[4], "Invalid flow: {}", var);
            assert!(
                &m[1] == "o" && &m[4] == "i",
                "Flow edges must be directed from o (out) to i (in), but \"{}\" isn't",
                var
            );
            let v = Extremity {
                id: usize::from_str(&m[5]).unwrap(),
                etype: match &m[6] {
                    "h" => ExtremityType::Head,
                    _ => ExtremityType::Tail,
                },
            };
            edges
                .entry(u)
                .or_insert(FxHashSet::default())
                .insert((v, val));
        } else if let Some(m) = pat_node.captures(var) {
            let v = Extremity {
                id: usize::from_str(&m[2]).unwrap(),
                etype: match &m[3] {
                    "h" => ExtremityType::Head,
                    _ => ExtremityType::Tail,
                },
            };
            let f = match &m[1] {
                "i" => Direction::In,
                _ => Direction::Out,
            };
            nodes.insert((v, f), val);
        }
    }

    Ok(Flow::new(nodes, edges))
}

pub fn extract_random_walk_from_flow(
    edges: &mut FxHashMap<Extremity, FxHashSet<(Extremity, usize)>>,
    start: &Extremity,
) -> Vec<Handle> {
    let mut res: Vec<Handle> = Vec::new();

    let mut v = start.clone();

    res.push(Handle::pack(start.id, start.etype == ExtremityType::Tail));

    while edges.contains_key(&v) {
        let mut neighbors = edges.get_mut(&v).unwrap();
        let (mut u, mut w) = pop(&mut neighbors);
        while w < 1 {
            let x = pop(&mut neighbors);
            u = x.0;
            w = x.1;
        }
        if w > 1 {
            neighbors.insert((u.clone(), w - 1));
        } else if neighbors.is_empty() {
            edges.remove(&v);
        }

        v = Extremity {
            id: u.id,
            etype: match u.etype {
                ExtremityType::Head => ExtremityType::Tail,
                ExtremityType::Tail => ExtremityType::Head,
            },
        };
        res.push(Handle::pack(v.id, v.etype == ExtremityType::Tail));
        // stop if start is revisted
        if &v == start {
            break;
        }
    }

    res
}
