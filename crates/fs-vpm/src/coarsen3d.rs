//! Invariant-preserving near→mid coarsening (bead
//! wf-root-guzez.5.18.2, E4.7-ii). Pairs of ADJACENT shed rows merge
//! into one mid-wake row under rules that preserve, per invariant
//! (each with its own executed oracle in the battery):
//!
//!   1. Kelvin closure — the merged wake is still a closed quad
//!      lattice (structural; the per-cell oracle re-runs on it),
//!   2. connectivity — n_stations + 1 nodes per row, loops closed,
//!   3. full vector hydrodynamic impulse — each merged quad retains
//!      γ_a A_a + γ_b A_b to floating-point roundoff, where A is its
//!      oriented vector area; unrepresentable warped pairs refuse,
//!   4. merged nodes use the pair's absolute-circulation weights,
//!   5. the squared core smoothing scale includes all three components
//!      of geometric separation (not a resolved viscous second moment),
//!   6. symmetry — a y-symmetric wake coarsens to a y-symmetric wake.
//!
//! The mixed-norm error metric (near-probe velocity RMS delta +
//! far-probe delta) is REPORTED per coarsen — never forced to vanish.

use crate::filament3d::{FilamentWake, Refusal, ShedRow, ring_vector_area};

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
/// interface, or a merged quad unable to represent all three components of
/// the original impulse. Every refusal leaves the wake unchanged.
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

    // Pass 1: all merged node lines, so each candidate quad's actual
    // successor geometry is known before circulation is chosen.
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
                        (wa / wt) * na[0] + (wb / wt) * nb[0],
                        (wa / wt) * na[1] + (wb / wt) * nb[1],
                        (wa / wt) * na[2] + (wb / wt) * nb[2],
                    ]
                })
                .collect(),
        );
    }
    // One row-wide smoothing scale, measured at the same middle node as
    // the existing reduced core model, now using its full 3-D separation.
    let mid = n / 2;
    // Pass 2: solve Γ_new A_new = Γ_a A_a + Γ_b A_b for each station.
    // A scalar strength cannot match a vector in an incompatible direction;
    // such warped pairs must remain resolved rather than lose impulse.
    let mut merged: Vec<ShedRow> = Vec::with_capacity(k);
    for pair in 0..k {
        let a = &wake.rows[2 * pair];
        let b = &wake.rows[2 * pair + 1];
        let (wa, wb) = pair_w[pair];
        let wt = wa + wb;
        let next_old = wake
            .rows
            .get(2 * pair + 2)
            .map_or(&wake.line_nodes, |row| &row.nodes);
        let nodes_new = &merged_nodes[pair];
        let next_new = if pair + 1 < k {
            &merged_nodes[pair + 1]
        } else if count < wake.rows.len() {
            &wake.rows[count].nodes
        } else {
            &wake.line_nodes
        };
        let area = |old: &[[f64; 3]], new: &[[f64; 3]], s: usize| {
            ring_vector_area(new[s], new[s + 1], old[s + 1], old[s])
        };
        let mut gamma = Vec::with_capacity(n);
        for s in 0..n {
            gamma.push(merged_circulation(
                a.gamma[s],
                area(&a.nodes, &b.nodes, s),
                b.gamma[s],
                area(&b.nodes, next_old, s),
                area(nodes_new, next_new, s),
            )?);
        }
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
    let tail = wake.rows.split_off(count);
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

fn merged_circulation(
    ga: f64,
    aa: [f64; 3],
    gb: f64,
    ab: [f64; 3],
    area: [f64; 3],
) -> Result<f64, Refusal> {
    let mut pivot = 0;
    for c in 1..3 {
        if area[c].abs() > area[pivot].abs() {
            pivot = c;
        }
    }
    if area[pivot] == 0.0 {
        return Err(Refusal {
            code: "coarsen-degenerate",
            message: "merged quad has zero representable vector area".into(),
            ranked_repairs: vec![
                "retain the resolved rows or rescale underflowing geometry".into(),
            ],
        });
    }
    let ia = aa.map(|v| ga * v);
    let ib = ab.map(|v| gb * v);
    let target: [f64; 3] = std::array::from_fn(|c| ia[c] + ib[c]);
    let gamma = target[pivot] / area[pivot];
    let represented = area.map(|v| gamma * v);
    let mut scale = 0.0_f64;
    for c in 0..3 {
        scale = scale
            .max(ia[c].abs())
            .max(ib[c].abs())
            .max(represented[c].abs());
    }
    // A vector-wide scale admits rotated cancellations in near-zero components.
    // There is no unit-dependent absolute floor. If tolerance underflows to
    // zero, only an exactly represented component can pass.
    let tolerance = 128.0 * f64::EPSILON * scale;
    if !gamma.is_finite()
        || ia
            .iter()
            .chain(&ib)
            .chain(&area)
            .chain(&target)
            .chain(&represented)
            .any(|v| !v.is_finite())
        || (0..3).any(|c| (represented[c] - target[c]).abs() > tolerance)
    {
        return Err(Refusal {
            code: "coarsen-impulse-unrepresentable",
            message: "one merged circulation cannot retain the original vector impulse".into(),
            ranked_repairs: vec![
                "retain the resolved rows for this warped or ill-scaled pair".into(),
            ],
        });
    }
    Ok(gamma)
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
