//! Linear and angular impulse on the actual oriented filament panels.

use super::{COARSEN_IMPULSE_RELATIVE_ROUNDOFF, FilamentWake, Refusal, ShedRow, unrepresentable};
use crate::filament3d::{ring_angular_impulse, ring_vector_area};
use fs_exec::reduce::Compensated;
use fs_la::factor::svd_jacobi;

fn norm(v: &[f64]) -> f64 {
    v[0].hypot(v[1]).hypot(v[2])
}

fn weighted_sum(moments: &[[f64; 6]], strengths: &[f64]) -> [f64; 6] {
    let mut sum = [Compensated::zero(); 6];
    for (moment, strength) in moments.iter().zip(strengths) {
        for d in 0..6 {
            sum[d] = sum[d].accumulate(moment[d] * strength);
        }
    }
    sum.map(Compensated::value)
}

fn allowance(moments: &[[f64; 6]], strengths: &[f64]) -> [f64; 2] {
    core::array::from_fn(|block| {
        moments
            .iter()
            .zip(strengths)
            .fold(Compensated::zero(), |sum, (moment, gamma)| {
                sum.accumulate(norm(&moment[3 * block..3 * block + 3]) * gamma.abs())
            })
            .value()
            * COARSEN_IMPULSE_RELATIVE_ROUNDOFF
    })
}

fn within(actual: [f64; 6], target: [f64; 6], tolerance: [f64; 2]) -> bool {
    let residual: [f64; 6] = core::array::from_fn(|d| actual[d] - target[d]);
    (0..2).all(|block| {
        let error = norm(&residual[3 * block..3 * block + 3]);
        error.is_finite() && error <= tolerance[block]
    })
}

struct StationMoments {
    original: Vec<[f64; 6]>,
    retained: Vec<[f64; 6]>,
}

/// A common LOCAL origin is essential: preserving J there together with I
/// preserves J at any other origin, without quadratic world-offset cancellation.
/// Normalize coordinates before area/cubic products, then normalize the linear
/// and angular blocks separately. No dimensional impulse product is needed.
fn station_moments<'a>(
    wake: &'a FilamentWake,
    k: usize,
    s: usize,
    retained_nodes: impl Fn(usize) -> &'a [[f64; 3]],
) -> Result<StationMoments, Refusal> {
    let boundary = wake
        .rows
        .get(2 * k)
        .map_or(&wake.line_nodes, |row| &row.nodes);
    let origin: [f64; 3] =
        core::array::from_fn(|d| 0.5 * boundary[s][d] + 0.5 * boundary[s + 1][d]);
    let mut length = 0.0_f64;
    for nodes in wake.rows[..2 * k]
        .iter()
        .map(|row| row.nodes.as_slice())
        .chain(core::iter::once(boundary.as_slice()))
        .chain((0..k).map(&retained_nodes))
    {
        for point in &nodes[s..s + 2] {
            let r: [f64; 3] = core::array::from_fn(|d| point[d] - origin[d]);
            let distance = norm(&r);
            if !distance.is_finite() {
                return Err(unrepresentable(
                    "local coarsening extent is unrepresentable",
                ));
            }
            length = length.max(distance);
        }
    }
    let panel = |older: &[[f64; 3]], newer: &[[f64; 3]]| {
        if length == 0.0 {
            return [0.0; 6];
        }
        let [a, b, c, d] = [newer[s], newer[s + 1], older[s + 1], older[s]]
            .map(|p| core::array::from_fn(|i| (p[i] - origin[i]) / length));
        let linear = ring_vector_area(a, b, c, d);
        let angular = ring_angular_impulse(a, b, c, d, [0.0; 3]);
        [
            linear[0], linear[1], linear[2], angular[0], angular[1], angular[2],
        ]
    };
    let mut original = Vec::with_capacity(2 * k);
    for r in 0..2 * k {
        let newer = wake
            .rows
            .get(r + 1)
            .map_or(&wake.line_nodes, |row| &row.nodes);
        original.push(panel(&wake.rows[r].nodes, newer));
    }
    let mut retained = Vec::with_capacity(k);
    for r in 0..k {
        let newer = if r + 1 < k {
            retained_nodes(r + 1)
        } else {
            boundary
        };
        retained.push(panel(retained_nodes(r), newer));
    }
    for block in 0..2 {
        let start = 3 * block;
        let scale = original
            .iter()
            .chain(&retained)
            .map(|moment| norm(&moment[start..start + 3]))
            .fold(0.0, f64::max);
        if !scale.is_finite() {
            return Err(unrepresentable(
                "normalized filament moment is unrepresentable",
            ));
        }
        if scale > 0.0 {
            for moment in original.iter_mut().chain(&mut retained) {
                for value in &mut moment[start..start + 3] {
                    *value /= scale;
                }
            }
        }
    }
    Ok(StationMoments { original, retained })
}

/// Retaining each pair's oldest boundary exactly telescopes constant-strength
/// straight filaments. Accept this geometry only when EVERY pair/station's
/// actual dimensional circulation represents both source moments. Otherwise
/// the caller uses its bounded weighted-geometry candidate and station solve.
pub(super) fn pairwise_oldest(
    wake: &FilamentWake,
    k: usize,
) -> Result<Option<Vec<Vec<f64>>>, Refusal> {
    let mut gamma = vec![vec![0.0; wake.cert.n_stations]; k];
    for s in 0..wake.cert.n_stations {
        let moments = station_moments(wake, k, s, |r| &wake.rows[2 * r].nodes)?;
        for r in 0..k {
            let scale = wake.rows[2 * r].gamma[s]
                .abs()
                .max(wake.rows[2 * r + 1].gamma[s].abs());
            if scale == 0.0 {
                continue;
            }
            let strengths = [
                wake.rows[2 * r].gamma[s] / scale,
                wake.rows[2 * r + 1].gamma[s] / scale,
            ];
            let original = &moments.original[2 * r..2 * r + 2];
            let target = weighted_sum(original, &strengths);
            let retained = moments.retained[r];
            let tolerance = allowance(original, &strengths);
            // Preserve an already-exact common circulation bit for bit. Equal
            // strengths alone do not admit geometry: both moment gates still apply.
            let common = wake.rows[2 * r].gamma[s];
            if common == wake.rows[2 * r + 1].gamma[s]
                && within(retained.map(|m| m * (common / scale)), target, tolerance)
            {
                gamma[r][s] = common;
                continue;
            }
            let pivot = (0..6)
                .max_by(|&i, &j| retained[i].abs().total_cmp(&retained[j].abs()))
                .expect("six moments");
            let value = if retained[pivot] == 0.0 {
                0.0
            } else {
                target[pivot] / retained[pivot]
            };
            let dimensional = value * scale;
            if !dimensional.is_finite()
                || (value != 0.0 && dimensional == 0.0)
                || !within(
                    retained.map(|m| m * (dimensional / scale)),
                    target,
                    tolerance,
                )
            {
                return Ok(None);
            }
            gamma[r][s] = dimensional;
        }
    }
    Ok(Some(gamma))
}

pub(super) fn preserve(wake: &FilamentWake, merged: &mut [ShedRow]) -> Result<bool, Refusal> {
    let k = merged.len();
    let mut corrected = false;
    for s in 0..wake.cert.n_stations {
        let gamma_scale = wake.rows[..2 * k]
            .iter()
            .map(|row| row.gamma[s].abs())
            .fold(0.0, f64::max);
        if gamma_scale == 0.0 {
            continue;
        }
        let moments = station_moments(wake, k, s, |r| &merged[r].nodes)?;
        let original_gamma: Vec<f64> = wake.rows[..2 * k]
            .iter()
            .map(|row| row.gamma[s] / gamma_scale)
            .collect();
        let target = weighted_sum(&moments.original, &original_gamma);
        let tolerance = allowance(&moments.original, &original_gamma);

        // Begin with the closest pairwise linear-impulse representation. The
        // correction below changes only this station's retained circulations.
        let mut gamma = Vec::with_capacity(k);
        for r in 0..k {
            let retained = moments.retained[r];
            let magnitude = norm(&retained[..3]);
            let pair = weighted_sum(
                &moments.original[2 * r..2 * r + 2],
                &original_gamma[2 * r..2 * r + 2],
            );
            let value = if magnitude == 0.0 {
                0.0
            } else {
                (0..3)
                    .map(|d| pair[d] * (retained[d] / magnitude))
                    .sum::<f64>()
                    / magnitude
            };
            if !value.is_finite() {
                return Err(unrepresentable("merged circulation is unrepresentable"));
            }
            gamma.push(value);
        }
        let represented = weighted_sum(&moments.retained, &gamma);
        if !within(represented, target, tolerance) {
            corrected = true;
            let residual: [f64; 6] = core::array::from_fn(|d| target[d] - represented[d]);
            // Factor the k-by-6 transpose directly: no Gram matrix squares its
            // condition number. Padding allows rank-deficient stations, including
            // the common rank-two straight sheet. A = U Sigma V^T gives the
            // minimum Euclidean circulation change U Sigma^+ V^T residual.
            let mut matrix = vec![0.0; 6 * k.max(6)];
            for (r, moment) in moments.retained.iter().enumerate() {
                matrix[6 * r..6 * r + 6].copy_from_slice(moment);
            }
            let svd = svd_jacobi(&matrix, k.max(6), 6);
            let rank_floor = 64.0 * f64::EPSILON * svd.sigma[0];
            for mode in 0..6 {
                if svd.sigma[mode] > rank_floor {
                    let projection = (0..6)
                        .map(|d| svd.v[6 * d + mode] * residual[d])
                        .sum::<f64>()
                        / svd.sigma[mode];
                    for (r, value) in gamma.iter_mut().enumerate() {
                        *value += svd.u[6 * r + mode] * projection;
                    }
                }
            }
        }
        // Gate the ACTUAL dimensional result separately against the original
        // unsigned linear and angular moments. Rank truncation cannot relax
        // either bound, including when a signed target happens to be zero.
        for (row, value) in merged.iter_mut().zip(&mut gamma) {
            row.gamma[s] = *value * gamma_scale;
            if !row.gamma[s].is_finite() || (*value != 0.0 && row.gamma[s] == 0.0) {
                return Err(unrepresentable("corrected circulation is unrepresentable"));
            }
            *value = row.gamma[s] / gamma_scale;
        }
        let actual = weighted_sum(&moments.retained, &gamma);
        if !within(actual, target, tolerance) {
            return Err(Refusal {
                code: "coarsen-impulse-unrepresentable",
                message: format!(
                    "station {s} cannot retain both linear and angular impulse on the candidate geometries"
                ),
                ranked_repairs: vec![
                    "retain more rows or geometry spanning both station moments".into(),
                ],
            });
        }
    }
    Ok(corrected)
}
