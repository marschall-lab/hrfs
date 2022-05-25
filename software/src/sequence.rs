/* standard use */
use std::io::{self, Write};
use std::str;
use std::str::FromStr;

/* crate use */
use handlegraph::{
    handle::{Edge, Handle},
    handlegraph::*,
    hashgraph::HashGraph,
    pathhandlegraph::PathId,
};
use log::info;
use rustc_hash::FxHashSet;

pub fn v2extstr(v: &Handle) -> String {
    format!(
        "{}{}",
        v.unpack_number(),
        if v.is_reverse() { "t" } else { "h" }
    )
}

pub fn v2str(v: &Handle) -> String {
    format!(
        "{}{}",
        if v.is_reverse() { '<' } else { '>' },
        v.unpack_number()
    )
}

/* FIXME */
pub fn v2seq(v: &Vec<Handle>, sep: &str) -> String {
    v.iter().map(v2str).collect::<Vec<String>>().join(sep)
}
pub fn hv2seq(v: &FxHashSet<Handle>, sep: &str) -> String {
    v.iter().map(v2str).collect::<Vec<String>>().join(sep)
}

pub fn push_walk_step(walk: &mut Vec<Handle>, step: &Vec<u8>) -> Result<(), String> {
    let sid = usize::from_str(str::from_utf8(&step[1..]).unwrap()).unwrap();
    let is_rev = match step[0] {
        b'>' => Ok(false),
        b'<' => Ok(true),
        _ => Err(format!(
            "unknown orientation '{}' of segment {}",
            step[0], sid
        )),
    };
    if is_rev.is_ok() {
        walk.push(Handle::pack(sid, is_rev.unwrap()));
        Ok(())
    } else {
        Err(is_rev.err().unwrap())
    }
}

pub fn reverse_seq(seq: &Vec<Handle>) -> Vec<Handle> {
    let mut rev = seq.clone();
    rev.reverse();
    for i in 0..rev.len() {
        rev[i] = rev[i].flip();
    }
    rev
}

pub fn parse_walk(line: &String) -> Result<Vec<Handle>, String> {
    let line = line[line.find('\t').unwrap() + 1..].as_bytes().to_vec();

    let mut walk: Vec<Handle> = Vec::new();
    let mut step: Vec<u8> = Vec::new();
    for c in line {
        if (c == b'>' || c == b'<') && !step.is_empty() {
            push_walk_step(&mut walk, &step)?;
            step.clear();
        }
        step.push(c);
    }
    if !step.is_empty() {
        push_walk_step(&mut walk, &step)?;
    }
    Ok(walk)
}

pub fn write_founders<W: io::Write>(
    f: &Vec<Vec<Handle>>,
    out: &mut io::BufWriter<W>,
) -> Result<(), io::Error> {
    log::info!("writing founder sequences to file");
    f.iter()
        .enumerate()
        .try_for_each(|(i, c)| writeln!(out, "founder_seq{}\t{}", i, v2seq(c, "")))
}

pub fn write_gfa_header<W: io::Write>(out: &mut io::BufWriter<W>) -> Result<(), io::Error> {
    writeln!(out, "H\tVN:Z:1.0")
}

pub fn write_gfa_segment<W: io::Write>(
    v: &Handle,
    g: &HashGraph,
    out: &mut io::BufWriter<W>,
) -> Result<(), io::Error> {
    writeln!(
        out,
        "S\t{}\t{}",
        v.unpack_number(),
        String::from_utf8(g.sequence_vec(*v)).unwrap()
    )
}

pub fn write_gfa_link<W: io::Write>(
    u: &Handle,
    v: &Handle,
    out: &mut io::BufWriter<W>,
) -> Result<(), io::Error> {
    writeln!(
        out,
        "L\t{}\t{}\t{}\t{}\t0M",
        u.unpack_number(),
        if u.is_reverse() { '-' } else { '+' },
        v.unpack_number(),
        if v.is_reverse() { '-' } else { '+' }
    )
}

pub fn write_gfa_path<W: io::Write>(
    id: &PathId,
    g: &HashGraph,
    out: &mut io::BufWriter<W>,
) -> Result<(), io::Error> {
    let p = g.get_path(id).unwrap();
    writeln!(
        out,
        "P\t{}\t{}\t*",
        str::from_utf8(&p.name).unwrap(),
        p.nodes
            .iter()
            .map(|v| format!(
                "{}{}",
                v.unpack_number(),
                if v.is_reverse() { '-' } else { '+' }
            ))
            .collect::<Vec<String>>()
            .join(",")
    )
}

pub fn write_subset_gfa<W: io::Write>(
    graph: &HashGraph,
    subgraph_nodes: &FxHashSet<Handle>,
    subgraph_edges: &FxHashSet<Edge>,
    paths: &Vec<PathId>,
    out: &mut io::BufWriter<W>,
) -> Result<(), io::Error> {
    info!("writing subset gfa file");
    write_gfa_header(out)?;
    for v in subgraph_nodes.iter() {
        write_gfa_segment(&v, &graph, out)?;
    }
    for Edge(mut u, mut v) in graph.edges() {
        if subgraph_edges.contains(&Edge(u, v)) {
            if u.is_reverse() && v.is_reverse() {
                let w = u.flip();
                u = v.flip();
                v = w;
            }
            write_gfa_link(&u, &v, out)?;
        }
    }
    for path_id in paths.iter() {
        write_gfa_path(&path_id, &graph, out)?;
    }
    Ok(())
}

pub fn write_gfa<W: io::Write>(g: &HashGraph, out: &mut io::BufWriter<W>) -> Result<(), io::Error> {
    info!("writing gfa file");
    write_gfa_header(out)?;
    for v in g.handles() {
        write_gfa_segment(&v, &g, out)?;
    }
    for Edge(u, v) in g.edges() {
        write_gfa_link(&u, &v, out)?;
    }
    for p in g.paths.iter() {
        write_gfa_path(p.0, &g, out)?;
    }
    Ok(())
}
