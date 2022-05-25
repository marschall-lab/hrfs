/* standard use */
use std::convert::TryInto;
use std::fs;
use std::io::{self, BufRead, BufReader, Seek, Write};
use std::str::FromStr;

/* crate use */
use clap::Parser;
use itertools::Itertools;
use petgraph::graphmap::DiGraphMap;
use petgraph::Outgoing;
use regex::Regex;
use rustc_hash::{FxHashMap, FxHashSet};

/* private use */
use founderset as ff;

#[derive(clap::Parser, Debug)]
#[clap(
    version = "0.1",
    author = "Konstantinn Bonnet <bonnetk@hhu.de>, Daniel Doerr <daniel.doerr@hhu.de>",
    about = "Extract founder sequences from minimization solution"
)]
pub struct Command {
    #[clap(index = 1, help = "ilp solution", required = true)]
    pub sol: String,
}

fn parse_founder_sequences<R: BufRead>(
    f: &mut R,
) -> Result<DiGraphMap<ff::Node, ff::EdgeType>, io::Error> {
    log::info!("reading founder sequences");

    let mut g: DiGraphMap<ff::Node, ff::EdgeType> = DiGraphMap::new();

    let re = Regex::new(r"^x_(i|o)(\d+)(h|t)(\d+)_(i|o)(\d+)(h|t)(\d+)\s(.*)").expect("regex");
    f.lines().try_for_each(|i| {
        if let Err(e) = i {
            Err(e)
        } else {
            let i = &i.unwrap();
            if let Some(m) = re.captures(&i) {
                if f64::from_str(&m[9]).unwrap().round() != 0.0 {
                    let u = ff::str2node(&m[1], &m[2], &m[3], &m[4])
                        .expect(format!("expr {}: invalid node", i).as_str());
                    let v = ff::str2node(&m[5], &m[6], &m[7], &m[8])
                        .expect(format!("expr {}: invalid node", i).as_str());
                    g.add_node(u);
                    g.add_node(v);
                    g.add_edge(
                        u,
                        v,
                        if u.node == v.node {
                            ff::EdgeType::Dashed
                        } else {
                            ff::EdgeType::Solid
                        },
                    );
                }
            }
            Ok(())
        }
    })?;
    Ok(g)
}

fn parse_haplotype_assignments<R: BufRead>(
    f: &mut R,
    g: &DiGraphMap<ff::Node, ff::EdgeType>,
) -> Result<(FxHashMap<ff::Node, usize>, FxHashSet<ff::Node>), io::Error> {
    log::info!("reading colors");

    let mut cmap: FxHashMap<ff::Node, usize> = FxHashMap::default();
    let mut switch = FxHashSet::default();

    let re = Regex::new(r"^c_(i)(\d+)(h|t)(\d+)_(\d+)_(\d+)\s(.*)").expect("regex");
    let re2 = Regex::new(r"^t_(i|o)(\d+)(h|t)(\d+)_(i|o)(\d+)(h|t)(\d+)_(\d+)_(\d+)\s(.*)")
        .expect("regex");
    f.lines().try_for_each(|i| {
        if let Err(e) = i {
            Err(e)
        } else {
            let i = &i.unwrap();
            if let Some(m) = re.captures(&i) {
                let v = ff::str2node(&m[1], &m[2], &m[3], &m[4])
                    .expect(format!("expr {}: invalid node", i).as_str());
                let hap = usize::from_str(&m[5]).unwrap();
                let c = f64::from_str(&m[7]).unwrap().round() as usize;
                if c == 1 {
                    assert!(
                        !cmap.contains_key(&v),
                        "attempt to overwrite color {} of node {} with {}\n{}",
                        cmap.get(&v).unwrap(),
                        v,
                        hap,
                        i
                    );
                    cmap.insert(v, hap);
                }
                if cmap.contains_key(&v) {
                    log::debug!("{} has color {}", v, cmap.get(&v).unwrap());
                }
            } else if let Some(m) = re2.captures(&i) {
                if f64::from_str(&m[11]).unwrap().round() == 1.0 {
                    log::debug!("serpent {}", i);
                    let u = ff::str2node(&m[1], &m[2], &m[3], &m[4])
                        .expect(format!("expr {}: invalid node", i).as_str());
                    let v = ff::str2node(&m[5], &m[6], &m[7], &m[8])
                        .expect(format!("expr {}: invalid node", i).as_str());
                    // get all t_uvci = 1 to later detect t_uvci = 0 when x_uv = 1
                    if g.contains_edge(u, v) {
                        log::debug!("capturing switch at {} = 1 on matched dashed edge", u);
                        switch.insert(u);
                    }
                }
            }
            Ok(())
        }
    })?;
    Ok((cmap, switch))
}

fn walk_next(g: &DiGraphMap<ff::Node, ff::EdgeType>, u: ff::Node) -> Option<ff::Node> {
    match g.edges_directed(u, Outgoing).nth(0) {
        None => None,
        Some((_, v, _)) => Some(v),
    }
}

fn walk_sol(
    g: DiGraphMap<ff::Node, ff::EdgeType>,
    color: FxHashMap<ff::Node, usize>,
    switch: FxHashSet<ff::Node>,
) -> Vec<Vec<(u64, bool, bool, usize)>> {
    log::info!("walking founder sequences");

    ff::source_nodes(&g)
        .iter()
        .map(|s| {
            let mut f = Vec::new();
            let mut n = *s;
            let mut d = false;
            log::info!("new founder starting at source node {}", n);
            f.push((
                n.node,
                d,
                !switch.contains(&n),
                *color.get(&n).expect(&format!("node {} has no color!", n)),
            ));
            let mut seq = "walk".to_owned();
            while let Some(u) = walk_next(&g, n) {
                seq = format!("{}: {}⇒{}", seq, n, u);
                if n.direction == ff::Direction::Out && u.direction == ff::Direction::In {
                    if n.etype == u.etype {
                        d = !d;
                    }
                    f.push((
                        u.node,
                        d,
                        !switch.contains(&u),
                        *color.get(&u).expect(&format!("node {} has no color!", u)),
                    ))
                };
                n = u;
            }
            log::debug!("{}", seq);
            log::debug!(
                "founder: {}",
                f.iter()
                    .map(|(u, d, s, _)| format!(
                        "{}{}{}",
                        if *s { "|" } else { "" },
                        if *d { "<" } else { ">" },
                        u
                    ))
                    .join("")
            );
            f
        })
        .collect::<Vec<Vec<(u64, bool, bool, usize)>>>()
}

fn write_founders<W: io::Write>(
    fs: Vec<Vec<(u64, bool, bool, usize)>>,
    out: &mut io::BufWriter<W>,
) -> Result<(), io::Error> {
    log::info!("writing final haplotype-minimized founders");

    fs.iter().enumerate().try_for_each(|(i, f)| {
        let name = format!("founder_seq{}", i);
        writeln!(
            out,
            "{}\t{}",
            name,
            f.iter()
                .map(|(u, d, _, _)| format!("{}{}", if *d { "<" } else { ">" }, u))
                .join("")
        )?;
        writeln!(
            out,
            "{}\t {}",
            String::from_utf8(vec![b' '; name.len()]).unwrap(),
            f.iter()
                .chain(std::iter::once(f.iter().last().unwrap()))
                .tuple_windows()
                .map(|((u, _, _, _), (_, _, s, _))| format!(
                    "{}{}",
                    String::from_utf8(vec![b' '; u.to_string().len()]).unwrap(),
                    if *s { "|".to_string() } else { " ".to_string() }
                ))
                .join("")
        )?;
        writeln!(
            out,
            "{}\t{}",
            String::from_utf8(vec![b' '; name.len()]).unwrap(),
            f.iter()
                .chain(std::iter::once(f.iter().last().unwrap()))
                .tuple_windows()
                .enumerate()
                .map(|(i, ((_, _, s, _), (_, _, _, c)))| {
                    let s = if i == 0 || *s {
                        c.to_string()
                    } else {
                        ".".to_owned()
                    };
                    let w = std::cmp::max(2i32 - (s.len() as i32), 0);
                    format!(
                        "{}{}",
                        s,
                        String::from_utf8(vec![b'.'; w.try_into().unwrap()]).unwrap()
                    )
                })
                .join("")
        )
    })
}

fn main() -> Result<(), io::Error> {
    env_logger::init();
    let params = Command::parse();

    let f = fs::File::open(&params.sol).expect("can't open sol file");
    let mut sol = BufReader::new(f);
    let g = parse_founder_sequences(&mut sol)?;
    sol.rewind()?;

    let (cols, switch) = parse_haplotype_assignments(&mut sol, &g)?;
    let fs = walk_sol(g, cols, switch);

    let mut out = io::BufWriter::new(std::io::stdout());
    write_founders(fs, &mut out)?;
    out.flush()?;

    log::info!("done");
    Ok(())
}
