//! Physical endpoint identity and closure, shared by evolution and diagnostics.

use std::collections::BTreeMap;

use super::FarField;
use crate::filament3d::{FilamentWake, Refusal};

pub(super) fn node_key(point: [f64; 3]) -> [u64; 3] {
    // Callers admit finite geometry. Signed zeros denote one physical point;
    // every other coordinate retains its exact identity.
    point.map(|v| if v == 0.0 { 0 } else { v.to_bits() })
}

/// Exact retained far sources precede the near sources, reproducing the
/// physical source order before aggregation. No images or far expansions.
pub(super) fn visit_physical_segments(
    near: &FilamentWake,
    far: &FarField,
    mut visit: impl FnMut([f64; 3], [f64; 3], f64),
) {
    for cell in &far.cells {
        for &(a, b, gamma) in &cell.segs {
            visit(a, b, gamma);
        }
    }
    near.visit_segments_with_core_growth(0.0, |s| visit(s.start, s.end, s.circulation));
}

#[derive(Default)]
struct Balance {
    sum: f64,
    correction: f64,
    magnitude: f64,
}

impl Balance {
    fn add(&mut self, value: f64) {
        let next = self.sum + value;
        self.correction += if self.sum.abs() >= value.abs() {
            (self.sum - next) + value
        } else {
            (value - next) + self.sum
        };
        self.sum = next;
        self.magnitude += value.abs();
    }

    fn closed(&self) -> bool {
        self.magnitude.is_finite()
            && (self.sum + self.correction).abs() <= 64.0 * f64::EPSILON * self.magnitude
    }
}

/// Caller admission bounds the source/cell counts and validates geometry.
/// Refusal context and cancellation ownership stay with the caller.
pub(super) fn circulation_is_closed(
    near: &FilamentWake,
    far: &FarField,
    mut checkpoint: impl FnMut() -> Result<(), Refusal>,
) -> Result<bool, Refusal> {
    let mut balances = BTreeMap::<[u64; 3], Balance>::new();
    let mut count = 0_usize;
    let mut failure = None;
    visit_physical_segments(near, far, |a, b, gamma| {
        if failure.is_some() {
            return;
        }
        if count % 256 == 0 {
            if let Err(error) = checkpoint() {
                failure = Some(error);
                return;
            }
        }
        count += 1;
        balances.entry(node_key(a)).or_default().add(-gamma);
        balances.entry(node_key(b)).or_default().add(gamma);
    });
    if let Some(error) = failure {
        return Err(error);
    }
    Ok(balances.values().all(Balance::closed))
}
