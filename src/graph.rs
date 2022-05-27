/* standard use */
use std::fmt;
use std::str::FromStr;
/* crate use */
use petgraph::{graphmap::DiGraphMap, Incoming, Outgoing};
use rustc_hash::FxHashSet;

#[derive(Clone, Copy, Debug, PartialOrd, Ord, PartialEq, Eq, Hash)]
pub enum ExtremityType {
    Head,
    Tail,
}

impl fmt::Display for ExtremityType {
    fn fmt(&self, fmt: &mut fmt::Formatter) -> fmt::Result {
        fmt.write_str(match &self {
            ExtremityType::Head => "h",
            ExtremityType::Tail => "t",
        })?;
        Ok(())
    }
}

pub fn str2ext(m: &str) -> Option<ExtremityType> {
    let c = m.chars().nth(0).unwrap();
    match c {
        't' => Some(ExtremityType::Tail),
        'h' => Some(ExtremityType::Head),
        _ => None,
    }
}

#[derive(Clone, Debug, PartialEq, Hash)]
pub struct Extremity {
    pub id: usize,
    pub etype: ExtremityType,
}

impl fmt::Display for Extremity {
    fn fmt(&self, fmt: &mut fmt::Formatter) -> fmt::Result {
        fmt.write_fmt(format_args!("{}{}", self.id, self.etype))?;
        Ok(())
    }
}

impl Eq for Extremity {}

#[derive(Clone, Copy, Debug, PartialOrd, Ord, PartialEq, Eq, Hash)]
pub enum Direction {
    In,
    Out,
}

impl fmt::Display for Direction {
    fn fmt(&self, fmt: &mut fmt::Formatter) -> fmt::Result {
        fmt.write_str(match &self {
            Direction::In => "i",
            Direction::Out => "o",
        })?;
        Ok(())
    }
}

pub fn str2dir(m: &str) -> Option<Direction> {
    let c = m.chars().nth(0).unwrap();
    match c {
        'i' => Some(Direction::In),
        'o' => Some(Direction::Out),
        _ => None,
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum EdgeType {
    Dashed,
    Solid,
}

impl fmt::Display for EdgeType {
    fn fmt(&self, fmt: &mut fmt::Formatter) -> fmt::Result {
        fmt.write_str(match &self {
            EdgeType::Dashed => "dashed",
            EdgeType::Solid => "solid",
        })?;
        Ok(())
    }
}
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Debug)]
pub struct Node {
    pub node: u64,
    pub etype: ExtremityType,
    pub direction: Direction,
    pub id: usize,
}

impl fmt::Display for Node {
    fn fmt(&self, fmt: &mut fmt::Formatter) -> fmt::Result {
        fmt.write_str(&format!(
            "{}{}{}{}",
            &self.direction, &self.node, &self.etype, &self.id
        ))
    }
}

pub fn str2node(d: &str, n: &str, e: &str, id: &str) -> Option<Node> {
    Some(Node {
        direction: str2dir(d)?,
        node: u64::from_str(n).unwrap(),
        etype: str2ext(e)?,
        id: usize::from_str(id).unwrap(),
    })
}

pub fn source_nodes(g: &DiGraphMap<Node, EdgeType>) -> FxHashSet<Node> {
    let mut s: FxHashSet<Node> = FxHashSet::default();
    s.extend(g.nodes().filter_map(|v| {
        if g.neighbors_directed(v, Incoming).count() == 0 {
            Some(v)
        } else {
            None
        }
    }));
    s
}

pub fn term_nodes(g: &DiGraphMap<Node, EdgeType>) -> (FxHashSet<Node>, FxHashSet<Node>) {
    let mut src = FxHashSet::default();
    let mut snk = FxHashSet::default();
    g.nodes().for_each(|v| {
        if g.neighbors_directed(v, Incoming).count() == 0 {
            src.insert(v);
        } else if g.neighbors_directed(v, Outgoing).count() == 0 {
            snk.insert(v);
        }
    });
    (src, snk)
}
