//! Stationwise physical impulse on the actual oriented panel geometry.

use super::{
    COARSEN_IMPULSE_RELATIVE_ROUNDOFF, FilamentWake, Refusal, ShedRow, ring_vector_area,
    unrepresentable,
};
use fs_exec::reduce::Compensated;
use fs_la::factor::svd_jacobi;

fn norm(v: [f64; 3]) -> f64 {
    v[0].hypot(v[1]).hypot(v[2])
}

/// Reuse the checked diagnostic's geometry and the actual ring orientation.
fn area(older: &[[f64; 3]], newer: &[[f64; 3]], s: usize) -> Result<[f64; 3], Refusal> {
    let value = ring_vector_area(newer[s], newer[s + 1], older[s + 1], older[s]);
    if value.iter().all(|v| v.is_finite()) {
        Ok(value)
    } else {
        Err(unrepresentable("panel area is unrepresentable"))
    }
}

fn weighted_sum(areas: &[[f64; 3]], strengths: &[f64]) -> [f64; 3] {
    let mut sum = [Compensated::zero(); 3];
    for (area, strength) in areas.iter().zip(strengths) {
        for d in 0..3 {
            sum[d] = sum[d].accumulate(area[d] * strength);
        }
    }
    sum.map(Compensated::value)
}

pub(super) fn preserve(
    wake: &FilamentWake,
    merged: &mut [ShedRow],
    stationwise: &[bool],
) -> Result<(), Refusal> {
    let k = merged.len();
    for s in 0..wake.cert.n_stations {
        if !stationwise[s] {
            continue;
        }
        let gamma_scale = wake.rows[..2 * k]
            .iter()
            .map(|row| row.gamma[s].abs())
            .fold(0.0, f64::max);
        if gamma_scale == 0.0 {
            continue;
        }
        let mut original = Vec::with_capacity(2 * k);
        for r in 0..2 * k {
            let newer = wake
                .rows
                .get(r + 1)
                .map_or(&wake.line_nodes, |row| &row.nodes);
            original.push(area(&wake.rows[r].nodes, newer, s)?);
        }
        let mut retained = Vec::with_capacity(k);
        for r in 0..k {
            let newer = merged.get(r + 1).map_or_else(
                || {
                    wake.rows
                        .get(2 * k)
                        .map_or(&wake.line_nodes, |row| &row.nodes)
                },
                |row| &row.nodes,
            );
            retained.push(area(&merged[r].nodes, newer, s)?);
        }
        // Separate geometry and circulation scales keep the solve independent
        // of length/circulation units, without forming their possibly huge product.
        let area_scale = original
            .iter()
            .chain(&retained)
            .flatten()
            .map(|v| v.abs())
            .fold(0.0, f64::max);
        if area_scale == 0.0 {
            continue;
        }
        for area in original.iter_mut().chain(&mut retained) {
            *area = area.map(|v| v / area_scale);
        }
        let original_gamma: Vec<f64> = wake.rows[..2 * k]
            .iter()
            .map(|row| row.gamma[s] / gamma_scale)
            .collect();
        let target = weighted_sum(&original, &original_gamma);
        let scale = original
            .iter()
            .zip(&original_gamma)
            .fold(Compensated::zero(), |sum, (area, gamma)| {
                sum.accumulate(norm(*area) * gamma.abs())
            })
            .value();
        let tolerance = COARSEN_IMPULSE_RELATIVE_ROUNDOFF * scale;

        // Start with each pair's closest scalar representation. On the original
        // parallel planar sheets this reduces to the old spacing-weighted law.
        let mut gamma = Vec::with_capacity(k);
        for r in 0..k {
            let magnitude = norm(retained[r]);
            let value = if magnitude == 0.0 {
                0.0
            } else {
                let pair: [f64; 3] = core::array::from_fn(|d| {
                    original[2 * r][d] * original_gamma[2 * r]
                        + original[2 * r + 1][d] * original_gamma[2 * r + 1]
                });
                (0..3)
                    .map(|d| pair[d] * (retained[r][d] / magnitude))
                    .sum::<f64>()
                    / magnitude
            };
            if !value.is_finite() {
                return Err(unrepresentable("merged circulation is unrepresentable"));
            }
            gamma.push(value);
        }
        let represented = weighted_sum(&retained, &gamma);
        let residual: [f64; 3] = core::array::from_fn(|d| target[d] - represented[d]);
        if norm(residual) > tolerance {
            // Factor the k-by-3 transpose directly: no Gram matrix squares the
            // condition number. Zero padding admits one/two-row rank deficiency.
            // A=U Sigma V^T, so delta_gamma=U Sigma^+ V^T residual is the
            // minimum Euclidean circulation change within THIS station alone.
            let mut matrix = vec![0.0; 3 * k.max(3)];
            for (r, area) in retained.iter().enumerate() {
                matrix[3 * r..3 * r + 3].copy_from_slice(area);
            }
            let svd = svd_jacobi(&matrix, k.max(3), 3);
            let rank_floor = 64.0 * f64::EPSILON * svd.sigma[0];
            for mode in 0..3 {
                if svd.sigma[mode] > rank_floor {
                    let projection = (0..3)
                        .map(|d| svd.v[3 * d + mode] * residual[d])
                        .sum::<f64>()
                        / svd.sigma[mode];
                    for (r, value) in gamma.iter_mut().enumerate() {
                        *value += svd.u[3 * r + mode] * projection;
                    }
                }
            }
        }
        // Check the ACTUAL returned strengths after dimensional conversion,
        // against the original panels. Rank truncation does not relax the bound.
        for (row, value) in merged.iter_mut().zip(&mut gamma) {
            row.gamma[s] = *value * gamma_scale;
            if !row.gamma[s].is_finite() {
                return Err(unrepresentable("corrected circulation is unrepresentable"));
            }
            *value = row.gamma[s] / gamma_scale;
        }
        let actual = weighted_sum(&retained, &gamma);
        let defect = norm(core::array::from_fn(|d| actual[d] - target[d]));
        if !defect.is_finite() || defect > tolerance {
            return Err(Refusal {
                code: "coarsen-impulse-unrepresentable",
                message: format!("station {s} vector-impulse defect {defect:e} exceeds scaled allowance {tolerance:e}"),
                ranked_repairs: vec!["retain more rows or geometry that spans this station's original vector impulse".into()],
            });
        }
    }
    Ok(())
}
