//! Invariant-preserving near→mid coarsening (bead
//! wf-root-guzez.5.18.2, E4.7-ii). Pairs of ADJACENT shed rows merge
//! into one mid-wake row under rules that preserve, per invariant
//! (each with its own executed oracle in the battery):
//!
//!   1. Kelvin closure — the merged wake is still a closed quad
//!      lattice (structural; the per-cell oracle re-runs on it),
//!   2. connectivity — n_stations + 1 nodes per row, loops closed,
//!   3. hydrodynamic impulse (streamwise-integrated spanwise
//!      circulation): γ_new = (γ_a·Δx_a + γ_b·Δx_b)/Δx_new — the
//!      per-station Σ γ·Δx is preserved EXACTLY by construction,
//!   4. first spatial moment — the merged row sits at the
//!      circulation-weighted centroid of the pair,
//!   5. core second moment — parallel-axis bookkeeping into
//!      `core2_m2` (retained, not discarded),
//!   6. symmetry — a y-symmetric wake coarsens to a y-symmetric wake.
//!
//! The mixed-norm error metric (near-probe velocity RMS delta +
//! far-probe delta) is REPORTED per coarsen — never forced to vanish.

use crate::filament3d::{FilamentWake, Refusal, ShedRow};

/// The reported coarsening error metric (mixed norm).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CoarsenMetric {
    /// RMS velocity delta over the near probe ring [m/s].
    pub near_rms_mps: f64,
    /// Velocity delta magnitude at the far probe [m/s].
    pub far_delta_mps: f64,
    /// Rows before → after.
    pub rows_before: usize,
    /// Rows after.
    pub rows_after: usize,
}

/// Circulation-weighted streamwise centroid shift and weights for one
/// row pair (span-aggregate: rows share their node line).
fn pair_weights(a: &ShedRow, b: &ShedRow) -> (f64, f64) {
    let wa: f64 = a.gamma.iter().map(|g| g.abs()).sum::<f64>().max(1e-12);
    let wb: f64 = b.gamma.iter().map(|g| g.abs()).sum::<f64>().max(1e-12);
    (wa, wb)
}

/// Coarsen the OLDEST `2*k` rows into `k` merged rows (mid-wake grows
/// from the old end; the fresh near-wake stays fully resolved).
///
/// # Errors
/// `coarsen-invalid` (k = 0, or fewer than 2k rows — cap AND cap+1 by
/// construction: exactly 2k admits, 2k−1 refuses).
pub fn coarsen_oldest(wake: &mut FilamentWake, k: usize) -> Result<CoarsenMetric, Refusal> {
    if k == 0 || wake.rows.len() < 2 * k {
        return Err(Refusal {
            code: "coarsen-invalid",
            message: format!("k {k} with {} rows", wake.rows.len()),
            ranked_repairs: vec!["coarsen at most rows/2 pairs from the old end".into()],
        });
    }
    // Probes for the REPORTED metric: a near ring around the newest
    // row's mid-span + one far point upstream.
    let n = wake.cert.n_stations;
    let newest = wake.rows.last().expect("non-empty").nodes[n / 2];
    let probes: Vec<[f64; 3]> = (0..8)
        .map(|i| {
            let ang = core::f64::consts::TAU * i as f64 / 8.0;
            [
                newest[0] + 0.5 * ang.cos(),
                newest[1],
                newest[2] + 0.5 * ang.sin(),
            ]
        })
        .collect();
    let far = [newest[0] + 20.0, newest[1], newest[2]];
    let before: Vec<[f64; 3]> = probes.iter().map(|p| wake.induced_velocity(*p)).collect();
    let far_before = wake.induced_velocity(far);
    let rows_before = wake.rows.len();

    // Pass 1: merged node lines (circulation-weighted centroids) so
    // every pair's REPRESENTED spacing is known before γ rescaling.
    let mut merged_nodes: Vec<Vec<[f64; 3]>> = Vec::with_capacity(k);
    let mut pair_w: Vec<(f64, f64)> = Vec::with_capacity(k);
    for pair in 0..k {
        let a = &wake.rows[2 * pair];
        let b = &wake.rows[2 * pair + 1];
        let (wa, wb) = pair_weights(a, b);
        let wt = wa + wb;
        pair_w.push((wa, wb));
        merged_nodes.push(
            a.nodes
                .iter()
                .zip(b.nodes.iter())
                .map(|(na, nb)| {
                    [
                        (wa * na[0] + wb * nb[0]) / wt,
                        (wa * na[1] + wb * nb[1]) / wt,
                        (wa * na[2] + wb * nb[2]) / wt,
                    ]
                })
                .collect(),
        );
    }
    // Span-aggregate x of an element (row or merged line).
    let mid = n / 2;
    let x_of = |nodes: &Vec<[f64; 3]>| nodes[mid][0];
    // Pass 2: γ rescaled so Σ γ·Δx is EXACT per station (invariant 3),
    // with Δx measured to the NEXT element in the post-merge sequence.
    let mut merged: Vec<ShedRow> = Vec::with_capacity(k);
    for pair in 0..k {
        let a = &wake.rows[2 * pair];
        let b = &wake.rows[2 * pair + 1];
        let (wa, wb) = pair_w[pair];
        let wt = wa + wb;
        // Old spacings: a→b, b→(next old element).
        let xa = x_of(&a.nodes);
        let xb = x_of(&b.nodes);
        let x_next_old = if 2 * pair + 2 < wake.rows.len() {
            x_of(&wake.rows[2 * pair + 2].nodes)
        } else {
            x_of(&wake.line_nodes)
        };
        // New spacing: merged pair → next merged line (or first tail
        // row / lifting line).
        let x_new = x_of(&merged_nodes[pair]);
        let x_new_next = if pair + 1 < k {
            x_of(&merged_nodes[pair + 1])
        } else if 2 * k < wake.rows.len() {
            x_of(&wake.rows[2 * k].nodes)
        } else {
            x_of(&wake.line_nodes)
        };
        let dx_a = xb - xa;
        let dx_b = x_next_old - xb;
        let dx_new = x_new_next - x_new;
        if dx_new.abs() < 1e-12 {
            return Err(Refusal {
                code: "coarsen-degenerate",
                message: "merged spacing collapsed to zero".into(),
                ranked_repairs: vec!["the wake rows are co-located; do not coarsen here".into()],
            });
        }
        let gamma: Vec<f64> = a
            .gamma
            .iter()
            .zip(b.gamma.iter())
            .map(|(ga, gb)| (ga * dx_a + gb * dx_b) / dx_new)
            .collect();
        // Invariant 5: parallel-axis core bookkeeping (streamwise).
        let core2 = (wa * (a.core2_m2 + (xa - x_new) * (xa - x_new))
            + wb * (b.core2_m2 + (xb - x_new) * (xb - x_new)))
            / wt;
        merged.push(ShedRow {
            nodes: merged_nodes[pair].clone(),
            gamma,
            core2_m2: core2,
        });
    }
    let tail = wake.rows.split_off(2 * k);
    wake.rows = merged;
    wake.rows.extend(tail);

    let after: Vec<[f64; 3]> = probes.iter().map(|p| wake.induced_velocity(*p)).collect();
    let far_after = wake.induced_velocity(far);
    let mut sum2 = 0.0;
    for (b, a) in before.iter().zip(after.iter()) {
        for c in 0..3 {
            sum2 += (b[c] - a[c]) * (b[c] - a[c]);
        }
    }
    let near_rms = (sum2 / probes.len() as f64).sqrt();
    let fd = ((far_before[0] - far_after[0]).powi(2)
        + (far_before[1] - far_after[1]).powi(2)
        + (far_before[2] - far_after[2]).powi(2))
    .sqrt();
    Ok(CoarsenMetric {
        near_rms_mps: near_rms,
        far_delta_mps: fd,
        rows_before,
        rows_after: wake.rows.len(),
    })
}

/// Per-station streamwise-integrated circulation (the impulse-class
/// invariant the battery holds exact): Σ over rows of γ[s]·Δx_row,
/// where Δx_row is the represented spacing (distance to the next
/// row's node line, span-aggregate).
#[must_use]
pub fn station_impulse(wake: &FilamentWake, s: usize) -> f64 {
    let n = wake.cert.n_stations;
    let mut total = 0.0;
    for (i, row) in wake.rows.iter().enumerate() {
        let x_here = row.nodes[n / 2][0];
        let x_next = if i + 1 < wake.rows.len() {
            wake.rows[i + 1].nodes[n / 2][0]
        } else {
            wake.line_nodes[n / 2][0]
        };
        total += row.gamma[s] * (x_next - x_here);
    }
    total
}

/// The FORBIDDEN naive decimation (battery falsifier ONLY — doc(hidden)
/// driver): drop every second row without re-weighting. Violates the
/// impulse invariant, which the battery proves by executing it.
#[doc(hidden)]
pub fn naive_decimate(wake: &mut FilamentWake) {
    let mut keep = Vec::new();
    for (i, row) in wake.rows.drain(..).enumerate() {
        if i % 2 == 0 {
            keep.push(row);
        }
    }
    wake.rows = keep;
}
