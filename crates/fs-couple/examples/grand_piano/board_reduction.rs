//! Cold bridge-driven Ritz projection of an already certified FE modal slice.
//! The low eigenpairs are copied exactly. Tail certificates belong to the
//! projected pencil; the original FE certificates remain in ReductionReport.

use super::{
    BoardGeometry, ReductionReport, cubic_triangle_shape, local_shape_displacement,
    nodal_displacement, ritz,
};
use fs_math::det;
use fs_modal::{ModePair, SliceReport};
use fs_plate::PlateModel;
use std::f64::consts::TAU;

pub(super) struct ReducedBoard {
    pub modes: Vec<ModePair>,
    pub physical_damping: Vec<f64>,
    pub report: ReductionReport,
}

fn validate_source(model: &PlateModel, source: &SliceReport) -> Result<(), String> {
    if source.below_low != 0 || source.expected != source.modes.len()
        || source.below_high.checked_sub(source.below_low) != Some(source.expected)
        || source.modes.is_empty() || source.modes.len() > ritz::MAX_SOURCE_MODES {
        return Err("board reduction requires a complete certified source slice of 1..=512 modes".into());
    }
    let mut m_phi = vec![0.0; model.free];
    for (i, pair) in source.modes.iter().enumerate() {
        if pair.phi.len() != model.free || pair.phi.iter().any(|x| !x.is_finite())
            || !pair.lambda.is_finite() || pair.lambda <= 0.0
            || !pair.residual.is_finite() || pair.residual < 0.0
            || !pair.interval.0.is_finite() || pair.interval.0 <= 0.0
            || !pair.interval.1.is_finite() || pair.interval.1 < pair.interval.0
            || pair.lambda < pair.interval.0 || pair.lambda > pair.interval.1
            || (i != 0 && pair.lambda < source.modes[i - 1].lambda) {
            return Err("invalid or nonpositive certified source soundboard eigenpair".into());
        }
        model.m.spmv(&pair.phi, &mut m_phi);
        let norm: f64 = pair.phi.iter().zip(&m_phi).map(|(a, b)| a * b).sum();
        if !norm.is_finite() || (norm - 1.0).abs() > 1e-7 {
            return Err("source soundboard eigenvector is not mass normalized".into());
        }
        for other in &source.modes[..i] {
            let product: f64 = other.phi.iter().zip(&m_phi).map(|(a, b)| a * b).sum();
            if !product.is_finite() || product.abs() > 1e-7 {
                return Err("source soundboard modes are not mass orthogonal".into());
            }
        }
    }
    Ok(())
}

pub(super) fn project(
    geometry: &BoardGeometry,
    model: &PlateModel,
    source: &SliceReport,
    keys: &[u8],
    edge_cubic_mass: bool,
    options: &ritz::RitzOptions,
) -> Result<ReducedBoard, String> {
    // No unvalidated source pair can disappear through projection.
    validate_source(model, source)?;
    let lambda: Vec<_> = source.modes.iter().map(|pair| pair.lambda).collect();
    let damping: Vec<_> = lambda.iter()
        .map(|value| 2.0 * geometry.damping_ratio * det::sqrt(*value)).collect();
    let mesh = &geometry.chart.mesh;
    let mut ports = Vec::with_capacity(keys.len());
    for &key in keys {
        let site = geometry.bridge_sites.iter().find(|site| site.midi == key)
            .ok_or_else(|| format!("geometry has no bridge station for key {key}"))?;
        let tri = mesh.tris[site.triangle];
        let cubic = edge_cubic_mass.then(|| cubic_triangle_shape(mesh, tri, site.weights));
        ports.push(source.modes.iter().map(|pair| {
            if let Some(shape) = &cubic {
                local_shape_displacement(model, &pair.phi, tri, shape)
            } else {
                (0..3).map(|i| site.weights[i]
                    * nodal_displacement(model, &pair.phi, tri[i])).sum()
            }
        }).collect::<Vec<f64>>());
    }
    let basis = ritz::bridge_basis(&lambda, &damping, &ports, options)?;
    let n = source.modes.len();
    let r = basis.columns.len();
    let low = options.keep_low_modes;
    let tail = r - low;
    let mut tail_k = vec![0.0; tail * tail];
    let mut identity = vec![0.0; tail * tail];
    for i in 0..tail {
        identity[i * tail + i] = 1.0;
        for j in 0..tail {
            tail_k[i * tail + j] = basis.stiffness[(low + i) * r + low + j];
        }
    }
    // Only the tail rotates. Even repeated protected low modes retain the
    // exact nodal eigenvectors, signs, eigenvalues and original certificates.
    let tail_pairs = fs_modal::eigh_gen_dense(&tail_k, &identity, tail)
        .map_err(|e| format!("reduced soundboard tail: {e}"))?;
    let mut modes = source.modes[..low].to_vec();
    let mut columns = basis.columns[..low].to_vec();
    for pair in tail_pairs {
        if !pair.lambda.is_finite() || pair.lambda <= 0.0
            || pair.phi.iter().any(|x| !x.is_finite())
            || !pair.interval.0.is_finite() || pair.interval.0 <= 0.0
            || !pair.interval.1.is_finite() || pair.interval.1 < pair.interval.0
            || !pair.residual.is_finite() || pair.residual < 0.0 {
            return Err("invalid or nonpositive projected soundboard eigenpair".into());
        }
        let mut column = vec![0.0; n];
        for (j, &weight) in pair.phi.iter().enumerate() {
            for (out, &entry) in column[low..].iter_mut().zip(&basis.columns[low + j][low..]) {
                *out += weight * entry;
            }
        }
        // Lower the final source-modal column to physical nodal DOFs once.
        // All bridge, volume, skin and motion projections use this same phi.
        let mut phi = vec![0.0; model.free];
        for (source_pair, &weight) in source.modes.iter().zip(&column) {
            if weight != 0.0 {
                for (out, &entry) in phi.iter_mut().zip(&source_pair.phi) {
                    *out += weight * entry;
                }
            }
        }
        modes.push(ModePair {
            lambda: pair.lambda, phi, residual: pair.residual, interval: pair.interval,
        });
        columns.push(column);
    }
    // The full congruence uses the final columns, including the tail rotation.
    // Computing from diagonal source C keeps symmetry exact and each diagonal
    // a sum of nonnegative terms. Off-diagonal physical loss is retained.
    let mut physical_damping = vec![0.0; r * r];
    for i in 0..r {
        for j in 0..=i {
            let value: f64 = (0..n)
                .map(|a| damping[a] * columns[i][a] * columns[j][a]).sum();
            if !value.is_finite() || (i == j && value < 0.0) {
                return Err("projected physical board damping overflow".into());
            }
            physical_damping[i * r + j] = value;
            physical_damping[j * r + i] = value;
        }
    }
    Ok(ReducedBoard {
        modes, physical_damping,
        report: ReductionReport {
            source_modes: n,
            protected_low_modes: low,
            source_frequency_intervals_hz: source.modes.iter().map(|pair|
                (det::sqrt(pair.interval.0) / TAU, det::sqrt(pair.interval.1) / TAU)).collect(),
            sample_hz: options.sample_hz.clone(),
            snapshot_count: basis.snapshot_count,
            max_relative_snapshot_error: basis.max_relative_snapshot_error,
        },
    })
}
