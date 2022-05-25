/* standard use */
use std::cmp::min;
use std::fs;
use std::io::{self, BufRead, Write};
use std::iter::FromIterator;

/* crate use */
use clap::Parser;
use handlegraph::handle::Handle;
use itertools::Itertools;
use rand::distributions::{Distribution, Uniform};
use rustc_hash::{FxHashMap, FxHashSet};

/* private use */
use founderset as ff;

#[derive(clap::Parser, Debug)]
#[clap(
    version = "0.1",
    author = "Konstantinn Bonnet <bonnetk@hhu.de>, Daniel Doerr <daniel.doerr@hhu.de>",
    about = "Minimize the number of recombinations on a set of fixed founder sequences by randomized coloring"
)]
pub struct Command {
    #[clap(index = 1, required = true, help = "founder sequences")]
    pub founder_set: String,

    #[clap(index = 2, required = true, help = "haplotype sequences")]
    pub haplotype_set: String,

    #[clap(
        short = 'n',
        long = "repeats",
        help = "Number of repeats for random color assignment",
        default_value = "10000"
    )]
    pub repeats: usize,
}

fn read_founderseqs<R: io::Read>(
    data: io::BufReader<R>,
) -> Result<Vec<(String, Vec<Handle>)>, String> {
    let mut seqs: Vec<(String, Vec<Handle>)> = Vec::new();

    for line_op in data.lines() {
        if let Ok(line) = line_op {
            seqs.push((
                line[..line.find('\t').unwrap()].to_string(),
                ff::parse_walk(&line)?,
            ));
        }
    }
    Ok(seqs)
}

fn haplotype_to_adj_map(
    haps: &Vec<(String, Vec<Handle>)>,
) -> FxHashMap<(Handle, Handle), FxHashSet<(usize, usize, bool)>> {
    let mut res: FxHashMap<(Handle, Handle), FxHashSet<(usize, usize, bool)>> =
        FxHashMap::default();
    for (x, (_, s)) in haps.into_iter().enumerate() {
        let n = s.len() - 1;
        s.into_iter()
            .tuple_windows()
            .enumerate()
            .for_each(|(i, (u, v))| {
                res.entry((v.flip(), u.flip()))
                    .or_insert(FxHashSet::default())
                    .insert((x, n - i, true));
                res.entry((*u, *v))
                    .or_insert(FxHashSet::default())
                    .insert((x, i, false));
            })
    }
    res
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

fn color_and_count(
    seq: &Vec<Handle>,
    adjs: &FxHashMap<(Handle, Handle), FxHashSet<(usize, usize, bool)>>,
    repeats: usize,
) -> usize {
    let mut rng = rand::thread_rng();

    // construct data structure for efficient access for random selection
    let adjs_vec: FxHashMap<(Handle, Handle), (Vec<(usize, usize, bool)>, Uniform<_>)> =
        FxHashMap::from_iter(adjs.iter().map(|((u, v), s)| {
            (
                (*u, *v),
                (s.iter().cloned().collect(), Uniform::from(0..s.len())),
            )
        }));

    let mut cur_min = usize::MAX;
    for _ in 0..repeats {
        let mut c = 0;
        let mut cur: Option<(usize, usize, bool)> = None;
        for (&u, &v) in seq.into_iter().tuple_windows() {
            let (vs, r) = adjs_vec.get(&(u, v)).expect(&format!(
                "oops, adjacency {}{} not contained in map, string is not colorable!",
                ff::v2str(&u),
                ff::v2str(&v)
            ));
            cur = Some(match cur {
                None => vs[r.sample(&mut rng)],
                Some((x, i, o)) => {
                    if vs.contains(&(x, i + 1, o)) {
                        (x, i + 1, o)
                    } else {
                        c += 1;
                        vs[r.sample(&mut rng)]
                    }
                }
            });
        }
        cur_min = min(cur_min, c);
    }
    cur_min
}

fn main() -> Result<(), std::io::Error> {
    env_logger::init();
    // initialize command line parser & parse command line arguments
    let params = Command::parse();

    log::info!(
        "performing random color assignment with {} repeats",
        params.repeats
    );

    log::info!("loading founder sequences from {}", params.founder_set);
    let founder_data = io::BufReader::new(fs::File::open(&params.founder_set)?);
    let founder_seqs = read_founderseqs(founder_data).unwrap();
    log::info!("parsed {} founders", founder_seqs.len());

    log::info!("loading haplotype sequences from {}", params.haplotype_set);
    let hap_data = io::BufReader::new(fs::File::open(&params.haplotype_set)?);
    let haplotypes = read_haplotypes(hap_data).unwrap();
    log::info!("parsed {} haplotypes", haplotypes.len());
    let hap_adjs = haplotype_to_adj_map(&haplotypes);
    log::debug!(
        "haplotype adjacency map: {}",
        hap_adjs
            .iter()
            .map(|((u, v), c)| format!(
                "{}{}: {}",
                ff::v2str(u),
                ff::v2str(v),
                c.iter()
                    .map(|(x, i, o)| format!("{}:{}{}", x, if *o { "+" } else { "-" }, i))
                    .join(", ")
            ))
            .collect::<Vec<String>>()
            .join("\n")
    );

    let mut out = io::BufWriter::new(std::io::stdout());

    let mut c = 0;
    for (name, s) in founder_seqs.iter() {
        let sc = color_and_count(s, &hap_adjs, params.repeats.clone());
        log::info!("{} has {} recombinations", name, &sc);
        c += sc;
    }
    log::info!("total #recombinations: {}", c);
    writeln!(out, "{}", c);
    out.flush()?;

    log::info!("done");
    Ok(())
}
