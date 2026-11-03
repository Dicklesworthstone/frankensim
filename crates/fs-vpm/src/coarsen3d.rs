//! Invariant-preserving near→mid coarsening (bead
//! wf-root-guzez.5.18.2, E4.7-ii). Pairs of ADJACENT shed rows merge
//! into one mid-wake row under rules that preserve, per invariant
//! (each with its own executed oracle in the battery):
//!
//!   1. Kelvin closure — the merged wake is still a closed quad
//!      lattice (structural; the per-cell oracle re-runs on it),
//!   2. connectivity — n_stations + 1 nodes per row, loops closed,
//!   3. full vector linear AND angular hydrodynamic impulse, separately
//!      retained per station to floating-point roundoff. A minimum-change
//!      circulation correction may redistribute between retained rows;
//!      unrepresentable candidate geometries refuse without changing the wake,
//!   4. oldest-node geometry when every pair exactly retains both moments
//!      (including a constant-circulation straight sheet); otherwise the
//!      pair's absolute-circulation-weighted node geometry,
//!   5. the squared core smoothing scale includes all three components
//!      of geometric separation (not a resolved viscous second moment),
//!   6. symmetry — a y-symmetric wake coarsens to a y-symmetric wake.
//!
//! The mixed-norm error metric (near-probe velocity RMS delta +
//! far-probe delta) is REPORTED per coarsen — never forced to vanish.

use crate::filament3d::{FilamentWake, Refusal, ShedRow};

mod impulse;

/// Maximum relative residual for EACH of linear and angular impulse, divided
/// by the station's original unsigned panel moments of that kind. Angular
/// moments use a common local origin. This is a floating-point allowance, not
/// a continuum error estimate; it has no dimensional floor.
pub const COARSEN_IMPULSE_RELATIVE_ROUNDOFF: f64 = 512.0 * f64::EPSILON;

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
    /// At least one station used the six-moment circulation correction,
    /// redistributing circulation among its retained rows to preserve both
    /// linear and angular impulse.
    pub used_stationwise_correction: bool,
}

/// Absolute-circulation weights for one
/// row pair (span-aggregate: rows share their node line).
fn pair_weights(a: &ShedRow, b: &ShedRow) -> (f64, f64) {
    let scale = a
        .gamma
        .iter()
        .chain(&b.gamma)
        .map(|g| g.abs())
        .fold(0.0_f64, f64::max);
    if scale == 0.0 {
        return (1.0, 1.0);
    }
    let wa = a.gamma.iter().map(|g| g.abs() / scale).sum();
    let wb = b.gamma.iter().map(|g| g.abs() / scale).sum();
    (wa, wb)
}

/// Coarsen the OLDEST `2*k` rows into `k` merged rows (mid-wake grows
/// from the old end; the fresh near-wake stays fully resolved).
///
/// # Errors
/// Invalid pair counts or source geometry/cores/ages, an attached far-field
/// interface, or the bounded pairwise/weighted-station constructions cannot
/// retain all six linear/angular impulse components. This is not a search
/// over arbitrary replacement geometries. Nonfinite candidate metrics refuse.
/// Every refusal leaves the wake unchanged.
pub fn coarsen_oldest(wake: &mut FilamentWake, k: usize) -> Result<CoarsenMetric, Refusal> {
    let count = k
        .checked_mul(2)
        .filter(|&count| count > 0 && count <= wake.rows.len());
    let Some(count) = count else {
        return Err(Refusal {
            code: "coarsen-invalid",
            message: format!("k {k} with {} rows", wake.rows.len()),
            ranked_repairs: vec!["coarsen at most rows/2 pairs from the old end".into()],
        });
    };
    if wake.downstream_gamma.is_some() {
        return Err(Refusal {
            code: "coarsen-interface-attached",
            message: "the oldest row is attached to retained far-field segments".into(),
            ranked_repairs: vec!["coarsen before transferring rows to the far field".into()],
        });
    }
    wake.admit_geometry()?;
    if wake
        .rows
        .iter()
        .any(|row| !row.core2_m2.is_finite() || row.core2_m2 < 0.0)
    {
        return Err(Refusal {
            code: "coarsen-core-invalid",
            message: "source rows need finite nonnegative squared cores".into(),
            ranked_repairs: vec!["retain finite nonnegative source cores".into()],
        });
    }
    if wake
        .rows
        .iter()
        .any(|row| !row.age_ticks.is_finite() || row.age_ticks < 0.0)
    {
        return Err(Refusal {
            code: "coarsen-age-invalid",
            message: "the source wake age is invalid".into(),
            ranked_repairs: vec!["retain finite nonnegative source ages".into()],
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

    // Try the oldest boundaries first: constant-strength straight filament
    // subdivisions telescope exactly there, including a single retained row.
    // Every pair/station must pass both moment gates before this geometry wins.
    let oldest_gamma = impulse::pairwise_oldest(wake, k)?;
    // All candidate node lines precede the station solve, so it uses the
    // actual successor geometry of every retained quad.
    let mut merged_nodes: Vec<Vec<[f64; 3]>> = Vec::with_capacity(k);
    let mut pair_w: Vec<(f64, f64)> = Vec::with_capacity(k);
    for pair in 0..k {
        let a = &wake.rows[2 * pair];
        let b = &wake.rows[2 * pair + 1];
        let (wa, wb) = pair_weights(a, b);
        let wt = wa + wb;
        pair_w.push((wa, wb));
        merged_nodes.push(if oldest_gamma.is_some() {
            a.nodes.clone()
        } else {
            a.nodes
                .iter()
                .zip(b.nodes.iter())
                .map(|(na, nb)| {
                    [
                        (wa / wt) * na[0] + (wb / wt) * nb[0],
                        (wa / wt) * na[1] + (wb / wt) * nb[1],
                        (wa / wt) * na[2] + (wb / wt) * nb[2],
                    ]
                })
                .collect()
        });
    }
    // One row-wide smoothing scale, measured at the same middle node as
    // the existing reduced core model, now using its full 3-D separation.
    let mid = n / 2;
    // Build the unpublished proposal. The weighted-geometry path solves its
    // full station moments after core/age bookkeeping, before publishing rows.
    let mut merged: Vec<ShedRow> = Vec::with_capacity(k);
    for pair in 0..k {
        let a = &wake.rows[2 * pair];
        let b = &wake.rows[2 * pair + 1];
        let (wa, wb) = pair_w[pair];
        let wt = wa + wb;
        let nodes_new = &merged_nodes[pair];
        let gamma = oldest_gamma
            .as_ref()
            .map_or_else(|| vec![0.0; n], |values| values[pair].clone());
        let separation2 = |source: [f64; 3]| -> f64 {
            (0..3)
                .map(|c| (source[c] - nodes_new[mid][c]).powi(2))
                .sum()
        };
        let core2 = (wa / wt) * (a.core2_m2 + separation2(a.nodes[mid]))
            + (wb / wt) * (b.core2_m2 + separation2(b.nodes[mid]));
        if !core2.is_finite() || nodes_new.iter().flatten().any(|v| !v.is_finite()) {
            return Err(Refusal {
                code: "coarsen-unrepresentable",
                message: "merged geometry or core exceeds the numerical range".into(),
                ranked_repairs: vec![
                    "retain these rows or rescale geometry to representable units".into(),
                ],
            });
        }
        let age_ticks = (wa / wt) * a.age_ticks + (wb / wt) * b.age_ticks;
        if !a.age_ticks.is_finite()
            || a.age_ticks < 0.0
            || !b.age_ticks.is_finite()
            || b.age_ticks < 0.0
            || !age_ticks.is_finite()
        {
            return Err(Refusal {
                code: "coarsen-age-invalid",
                message: "the source or weighted merged wake age is invalid".into(),
                ranked_repairs: vec!["retain finite nonnegative source ages".into()],
            });
        }
        merged.push(ShedRow {
            nodes: merged_nodes[pair].clone(),
            gamma,
            core2_m2: core2,
            age_ticks,
        });
    }
    let used_stationwise_correction = if oldest_gamma.is_some() {
        false
    } else {
        impulse::preserve(wake, &mut merged)?
    };
    // Move untouched rows and their node/circulation buffers into the trial.
    // Keep the original prefix for rollback. No callback or fallible operation
    // observes this temporary state before the metric check below.
    let tail = wake.rows.split_off(count);
    let original_prefix = std::mem::replace(&mut wake.rows, merged);
    wake.rows.extend(tail);

    let after: Vec<[f64; 3]> = probes.iter().map(|p| wake.induced_velocity(*p)).collect();
    let far_after = wake.induced_velocity(far);
    let mut sum2 = 0.0;
    for (b, a) in before.iter().zip(after.iter()) {
        for c in 0..3 {
            sum2 += (b[c] - a[c]) * (b[c] - a[c]);
        }
    }
    let mut near_rms = (sum2 / probes.len() as f64).sqrt();
    let mut fd = ((far_before[0] - far_after[0]).powi(2)
        + (far_before[1] - far_after[1]).powi(2)
        + (far_before[2] - far_after[2]).powi(2))
    .sqrt();
    // Preserve already-finite metric bits. Extreme but finite sampled fields
    // can overflow their squared errors while the final norm remains finite.
    if (!near_rms.is_finite() || !fd.is_finite())
        && before
            .iter()
            .chain(&after)
            .flatten()
            .chain(&far_before)
            .chain(&far_after)
            .all(|value| value.is_finite())
    {
        if !near_rms.is_finite() {
            let scale = (probes.len() as f64).sqrt();
            near_rms = 0.0;
            for (b, a) in before.iter().zip(&after) {
                for c in 0..3 {
                    let difference = b[c] - a[c];
                    let scaled = if difference.is_finite() {
                        difference / scale
                    } else {
                        b[c] / scale - a[c] / scale
                    };
                    near_rms = near_rms.hypot(scaled);
                }
            }
        }
        if !fd.is_finite() {
            fd = (far_before[0] - far_after[0])
                .hypot(far_before[1] - far_after[1])
                .hypot(far_before[2] - far_after[2]);
        }
    }
    if !near_rms.is_finite() || !fd.is_finite() {
        let tail = wake.rows.split_off(k);
        wake.rows = original_prefix;
        wake.rows.extend(tail);
        return Err(unrepresentable(
            "coarsening probe velocity or error exceeds the numerical range",
        ));
    }
    Ok(CoarsenMetric {
        near_rms_mps: near_rms,
        far_delta_mps: fd,
        rows_before,
        rows_after: wake.rows.len(),
        used_stationwise_correction,
    })
}

fn unrepresentable(message: &str) -> Refusal {
    Refusal {
        code: "coarsen-unrepresentable",
        message: message.into(),
        ranked_repairs: vec![
            "retain these rows or rescale to representable geometry and circulation".into(),
        ],
    }
}

/// Legacy streamwise-integrated circulation: Σ over rows of γ[s]·Δx_row,
/// where Δx_row is the represented spacing (distance to the next
/// row's node line, span-aggregate). This is proportional to a component
/// of physical impulse for parallel uniformly spaced span lines only.
/// Use [`FilamentWake::hydrodynamic_impulse_per_density`] for full 3-D geometry.
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
