//! Dörfler (fixed-energy) marking with DETERMINISTIC tie-breaking:
//! cells sorted by |indicator| descending, cell key ascending on ties;
//! the smallest prefix whose mass reaches θ·total is marked. Two runs
//! over the same indicators mark bitwise-identically (P2).

use std::collections::BTreeMap;

pub(crate) fn indicator_order(
    indicators: &BTreeMap<(u32, u32, u32), f64>,
) -> Vec<((u32, u32, u32), f64)> {
    let mut order: Vec<_> = indicators
        .iter()
        .map(|(&cell, &indicator)| (cell, indicator.abs()))
        .collect();
    order.sort_by(|a, b| {
        b.1.partial_cmp(&a.1)
            .expect("finite indicators")
            .then(a.0.cmp(&b.0))
    });
    order
}

/// Mark the smallest prefix carrying `theta` of the indicator mass.
#[must_use]
pub fn dorfler(indicators: &BTreeMap<(u32, u32, u32), f64>, theta: f64) -> Vec<(u32, u32, u32)> {
    let total: f64 = indicators.values().map(|v| v.abs()).sum();
    if total <= 0.0 {
        return Vec::new();
    }
    let mut marked = Vec::new();
    let mut mass = 0.0;
    for (c, v) in indicator_order(indicators) {
        if mass >= theta * total {
            break;
        }
        mass += v;
        marked.push(c);
    }
    marked
}
