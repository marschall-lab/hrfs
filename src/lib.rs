pub mod flow;
pub mod graph;
pub mod sequence;

/* crate use */
use rustc_hash::FxHashSet;

/* private use */
pub use crate::{flow::*, graph::*, sequence::*};

// copied from da internet
// split off an arbitrary element from a (non-empty) set
pub fn pop<T>(set: &mut FxHashSet<T>) -> T
where
    T: Eq + Clone + std::hash::Hash,
{
    let elt = set.iter().next().cloned().unwrap();
    set.remove(&elt);
    elt
}
