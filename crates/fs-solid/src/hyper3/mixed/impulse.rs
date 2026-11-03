//! Fixed-geometry momentum transfer on the original P2/P1 constraint.
//!
//! The current volume derivative B and retained consistent mass solve
//! [M_ff B_f^T; B_f 0] [delta-v_f, lambda] = [J_f, -B v_before].
//! The internal pressure impulse is -B^T lambda. Lambda has units Pa s;
//! it is not an endpoint pressure or a force over an invented time interval.
//! Pre-existing incompatible velocities refuse before projection. Original
//! volume rows, actual rounded endpoint velocities and pressure work are checked
//! independently of the scaled native sparse LDL solve.

use fs_exec::Cx;
use fs_solver::{ResidualClaim, SolveReport};
use fs_sparse::{Coo, Csr, DirectOrdering, Inertia, LdltOptions, SymbolicLdlt};

use super::element::PreparedMixedMesh;
use super::{MixedUpError, checkpoint, inertia_error, invalid};
use crate::inertia::{ReferenceMass, SolidImpulseReport, SolidMomentum, impulse_moments};

pub use crate::inertia::SolidImpulseControls;

/// Explicit impulse scales, physical allowances and native sparse work limits.
#[derive(Clone, Copy, Debug)]
pub struct MixedImpulseControls {
    /// Physical momentum/work gates and bounded shared accounting workspace.
    /// Its relative solve tolerance admits the recomputed scaled KKT residual.
    /// This direct lane performs no Krylov iterations or defect corrections.
    pub balance: SolidImpulseControls,
    pub velocity_scale_m_s: f64,
    pub pressure_impulse_scale_pa_s: f64,
    pub work_scale_j: f64,
    pub maximum_assembly_bytes: usize,
    pub maximum_factor_bytes: usize,
    /// Bounds every weak pressure-test volume rate / total reference volume,
    /// both before and after the event [1/s]. This is not pointwise div(v).
    pub maximum_relative_volume_rate_per_s: f64,
}

/// Pressure, applied excitation and support accounting at unchanged geometry.
#[derive(Clone, Debug)]
pub struct MixedImpulseReport {
    /// Shared original-mass accounting. `applied` and `applied_work_j` include
    /// BOTH external and internal pressure impulses, separated below.
    pub momentum: SolidImpulseReport,
    pub external: SolidMomentum,
    pub pressure: SolidMomentum,
    pub external_work_j: f64,
    pub pressure_work_j: f64,
    pub pressure_nodes: Vec<usize>,
    /// Tensile multiplier lambda, in sorted corner-node order [Pa s].
    /// A declared constant gauge is zero and never overwrites stored pressure.
    pub pressure_impulses_pa_s: Vec<f64>,
    pub volume_rate_before_m3_s: Vec<f64>,
    pub volume_rate_after_m3_s: Vec<f64>,
    pub kkt_inertia: Inertia,
}

#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
pub(super) fn constrained_impulse(
    mesh: &PreparedMixedMesh,
    mass: &ReferenceMass,
    positions_m: &[[f64; 3]],
    velocity_mps: &[[f64; 3]],
    held_velocity_dofs: &[usize],
    nodal_impulses: &[(usize, f64)],
    pressure_gauge_node: Option<usize>,
    origin_m: [f64; 3],
    controls: MixedImpulseControls,
    mut endpoint_cofactor: impl FnMut(usize, usize) -> Result<[f64; 9], MixedUpError>,
    cx: &Cx<'_>,
) -> Result<(Vec<[f64; 3]>, MixedImpulseReport), MixedUpError> {
    checkpoint(cx)?;
    let nu = mass.ndofs();
    let np = mesh.pressure_nodes.len();
    if positions_m.len() != mesh.nodes_count
        || velocity_mps.len() != mesh.nodes_count
        || nu != 3 * mesh.nodes_count
        || positions_m
            .iter()
            .chain(velocity_mps)
            .flatten()
            .any(|v| !v.is_finite())
        || origin_m.iter().any(|v| !v.is_finite())
        || [
            controls.velocity_scale_m_s,
            controls.pressure_impulse_scale_pa_s,
            controls.work_scale_j,
            controls.balance.solve.tolerance,
            controls.balance.linear_tolerance_kg_m_s,
            controls.balance.angular_tolerance_kg_m2_s,
            controls.balance.work_tolerance_j,
        ]
        .iter()
        .any(|v| !v.is_finite() || *v <= 0.0)
        || controls.balance.solve.tolerance >= 1.0
        || !controls.maximum_relative_volume_rate_per_s.is_finite()
        || controls.maximum_relative_volume_rate_per_s < 0.0
        || controls.maximum_factor_bytes == 0
    {
        return Err(invalid(
            "finite complete mixed impulse fields, positive scales and valid budgets required",
        ));
    }
    // Original gradient, its transpose, KKT COO/sort/CSR and vectors. The
    // native symbolic/numeric factor has its separate payload cap.
    let bytes = mesh
        .elements
        .len()
        .checked_mul(4 * 30 * 4)
        .and_then(|entries| entries.checked_add(mass.matrix().nnz()))
        .and_then(|entries| entries.checked_add(nu.checked_add(np)?))
        .and_then(|entries| entries.checked_mul(128))
        .ok_or_else(|| invalid("mixed impulse assembly size overflow"))?;
    if bytes > controls.maximum_assembly_bytes {
        return Err(invalid("mixed impulse exceeds assembly byte cap"));
    }
    let mut held = vec![false; nu];
    for (index, &dof) in held_velocity_dofs.iter().enumerate() {
        if index % 256 == 0 {
            checkpoint(cx)?;
        }
        if dof >= nu || held[dof] {
            return Err(invalid("distinct in-range held velocity DOFs required"));
        }
        held[dof] = true;
    }
    let free: Vec<_> = (0..nu).filter(|&d| !held[d]).collect();
    if free.is_empty() {
        return Err(invalid(
            "mixed pressure impulse requires free velocity components",
        ));
    }
    let gauge = pressure_gauge_node
        .map(|node| {
            mesh.pressure_nodes
                .binary_search(&node)
                .map_err(|_| invalid("pressure impulse gauge must name an original corner node"))
        })
        .transpose()?;
    let mut external = vec![0.0; nu];
    for (index, &(dof, value)) in nodal_impulses.iter().enumerate() {
        if index % 256 == 0 {
            checkpoint(cx)?;
        }
        if dof >= nu || !value.is_finite() {
            return Err(invalid("finite in-range applied impulses required"));
        }
        external[dof] = finite(external[dof] + value)?;
    }
    let mut gradient = Coo::new(np, nu);
    for (e, element) in mesh.elements.iter().enumerate() {
        let mut local = [[0.0; 30]; 4];
        for (q, point) in element.points.iter().enumerate() {
            checkpoint(cx)?;
            let cofactor = endpoint_cofactor(e, q)?;
            for a in 1..10 {
                for c in 0..3 {
                    let coupling = (0..3).fold(0.0, |sum, j| {
                        cofactor[3 * c + j].mul_add(point.gradients[a][j], sum)
                    });
                    for p in 0..4 {
                        let value = point.weight_m3 * point.pressure_shapes[p] * coupling;
                        local[p][3 * a + c] += value;
                        // Same anchor-difference derivative as native F(u).
                        local[p][c] -= value;
                    }
                }
            }
        }
        for p in 0..4 {
            for d in 0..30 {
                let value = finite(local[p][d])?;
                if value != 0.0 {
                    gradient.push(
                        element.pressure_indices[p],
                        3 * element.connectivity[d / 3] + d % 3,
                        value,
                    );
                }
            }
        }
    }
    let gradient = gradient.assemble();
    let transpose = fs_sparse::ops::transpose(&gradient);
    if let Some(pin) = gauge {
        let component = pressure_component(mesh, pin, cx)?;
        for (index, &dof) in free.iter().enumerate() {
            if index % 128 == 0 {
                checkpoint(cx)?;
            }
            let (columns, values) = transpose.row(dof);
            let (mut sum, mut absolute) = (0.0, 0.0);
            for (&p, &value) in columns.iter().zip(values) {
                if component[p] {
                    sum += value;
                    absolute += value.abs();
                }
            }
            if finite(sum)?.abs() > 4096.0 * f64::EPSILON * finite(absolute)?.max(f64::MIN_POSITIVE)
            {
                return Err(MixedUpError::IncompatiblePressureGauge {
                    node: mesh.pressure_nodes[pin],
                });
            }
        }
    }
    let maximum_rate = finite(controls.maximum_relative_volume_rate_per_s * mesh.total_volume_m3)?;
    let before_rate = volume_rates(&gradient, velocity_mps, maximum_rate, cx)?;
    let mut reduced = vec![None; nu];
    for (row, &dof) in free.iter().enumerate() {
        reduced[dof] = Some(row);
    }
    let mut pressure_rows = vec![None; np];
    let mut n = free.len();
    for (p, row) in pressure_rows.iter_mut().enumerate() {
        if Some(p) != gauge {
            *row = Some(n);
            n += 1;
        }
    }
    let vs = controls.velocity_scale_m_s;
    let ps = controls.pressure_impulse_scale_pa_s;
    let ws = controls.work_scale_j;
    let mut matrix = Coo::new(n, n);
    let mut rhs = vec![0.0; n];
    for (row, &dof) in free.iter().enumerate() {
        if row % 128 == 0 {
            checkpoint(cx)?;
        }
        rhs[row] = finite((vs / ws) * external[dof])?;
        let (columns, values) = mass.matrix().row(dof);
        for (&column, &value) in columns.iter().zip(values) {
            if let Some(j) = reduced[column] {
                matrix.push(row, j, finite((vs / ws) * value * vs)?);
            }
        }
    }
    for p in 0..np {
        let Some(row) = pressure_rows[p] else {
            continue;
        };
        rhs[row] = finite(-(ps / ws) * before_rate[p])?;
        let (columns, values) = gradient.row(p);
        for (&dof, &value) in columns.iter().zip(values) {
            if let Some(j) = reduced[dof] {
                let entry = finite((vs / ws) * value * ps)?;
                matrix.push(row, j, entry);
                matrix.push(j, row, entry);
            }
        }
    }
    let matrix = matrix.assemble();
    checkpoint(cx)?;
    let symbolic =
        SymbolicLdlt::analyze(&matrix, DirectOrdering::Amd).map_err(MixedUpError::Direct)?;
    checkpoint(cx)?;
    let factor = symbolic
        .factor(
            &matrix,
            &LdltOptions {
                max_factor_bytes: Some(controls.maximum_factor_bytes),
                ..LdltOptions::default()
            },
        )
        .map_err(MixedUpError::Direct)?;
    let expected = Inertia {
        positive: free.len(),
        negative: n - free.len(),
    };
    if factor.inertia() != expected {
        return Err(MixedUpError::KktInertia {
            expected,
            got: factor.inertia(),
        });
    }
    checkpoint(cx)?;
    let solution = factor.solve(&rhs);
    checkpoint(cx)?;
    let mut action = vec![0.0; n];
    matrix.spmv(&solution, &mut action);
    let relative = relative_residual(&rhs, &action)?;
    let solve = SolveReport::from_claim(
        0,
        ResidualClaim::TrueEuclidean(relative),
        controls.balance.solve.tolerance,
        vec![relative],
    );
    if !solve.converged_euclidean() {
        return Err(invalid(format!(
            "mixed impulse KKT residual {relative:e} exceeds tolerance"
        )));
    }
    drop((factor, symbolic, matrix, rhs, action));
    let mut velocity = velocity_mps.to_vec();
    for (row, &dof) in free.iter().enumerate() {
        velocity[dof / 3][dof % 3] = finite(velocity[dof / 3][dof % 3] + vs * solution[row])?;
    }
    let pressure_impulses_pa_s: Vec<_> = pressure_rows
        .iter()
        .map(|row| row.map_or(0.0, |r| ps * solution[r]))
        .collect();
    let mut pressure_impulse = vec![0.0; nu];
    transpose.spmv(&pressure_impulses_pa_s, &mut pressure_impulse);
    let mut effective = Vec::with_capacity(nu);
    let (mut external_work_j, mut pressure_work_j) = (0.0, 0.0);
    for dof in 0..nu {
        if dof % 256 == 0 {
            checkpoint(cx)?;
        }
        pressure_impulse[dof] = finite(-pressure_impulse[dof])?;
        let average = 0.5 * velocity_mps.as_flattened()[dof] + 0.5 * velocity.as_flattened()[dof];
        external_work_j = external[dof].mul_add(average, external_work_j);
        pressure_work_j = pressure_impulse[dof].mul_add(average, pressure_work_j);
        effective.push((dof, finite(external[dof] + pressure_impulse[dof])?));
    }
    let pressure =
        impulse_moments(positions_m, &pressure_impulse, origin_m, cx).map_err(inertia_error)?;
    if magnitude(pressure.linear_kg_m_s) > controls.balance.linear_tolerance_kg_m_s
        || magnitude(pressure.angular_kg_m2_s) > controls.balance.angular_tolerance_kg_m2_s
        || finite(pressure_work_j)?.abs() > controls.balance.work_tolerance_j
    {
        return Err(invalid(
            "internal pressure impulse violates momentum or zero-work allowance",
        ));
    }
    let after_rate = volume_rates(&gradient, &velocity, maximum_rate, cx)?;
    let (velocity, momentum) = mass
        .finish_impulse_candidate(
            positions_m,
            velocity_mps,
            held_velocity_dofs,
            &effective,
            origin_m,
            controls.balance,
            velocity,
            solve,
            cx,
        )
        .map_err(inertia_error)?;
    let report = MixedImpulseReport {
        momentum,
        external: impulse_moments(positions_m, &external, origin_m, cx).map_err(inertia_error)?,
        pressure,
        external_work_j: finite(external_work_j)?,
        pressure_work_j,
        pressure_nodes: mesh.pressure_nodes.clone(),
        pressure_impulses_pa_s,
        volume_rate_before_m3_s: before_rate,
        volume_rate_after_m3_s: after_rate,
        kkt_inertia: expected,
    };
    checkpoint(cx)?;
    Ok((velocity, report))
}

fn volume_rates(
    b: &Csr,
    velocity: &[[f64; 3]],
    maximum: f64,
    cx: &Cx<'_>,
) -> Result<Vec<f64>, MixedUpError> {
    checkpoint(cx)?;
    let mut rates = vec![0.0; b.nrows()];
    b.spmv(velocity.as_flattened(), &mut rates);
    if rates
        .iter()
        .any(|rate| !rate.is_finite() || rate.abs() > maximum)
    {
        return Err(invalid(
            "mixed impulse velocity violates its weak volume-rate allowance",
        ));
    }
    checkpoint(cx)?;
    Ok(rates)
}

fn pressure_component(
    mesh: &PreparedMixedMesh,
    pin: usize,
    cx: &Cx<'_>,
) -> Result<Vec<bool>, MixedUpError> {
    // Same connected constant-pressure gauge as the native mixed static owner.
    let mut roots: Vec<_> = (0..mesh.pressure_nodes.len()).collect();
    fn root(roots: &[usize], mut p: usize) -> usize {
        while roots[p] != p {
            p = roots[p];
        }
        p
    }
    for element in &mesh.elements {
        checkpoint(cx)?;
        for &p in &element.pressure_indices[1..] {
            let a = root(&roots, element.pressure_indices[0]);
            let b = root(&roots, p);
            roots[a.max(b)] = a.min(b);
        }
    }
    let component = root(&roots, pin);
    Ok((0..roots.len())
        .map(|p| root(&roots, p) == component)
        .collect())
}

fn finite(value: f64) -> Result<f64, MixedUpError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(invalid("unrepresentable mixed impulse"))
    }
}

fn relative_residual(rhs: &[f64], action: &[f64]) -> Result<f64, MixedUpError> {
    // Normalize before norms, so a finite but very large RHS cannot produce
    // an infinite denominator and a spurious zero residual claim.
    let (mut bscale, mut rscale) = (0.0_f64, 0.0_f64);
    for (&b, &a) in rhs.iter().zip(action) {
        bscale = bscale.max(b.abs());
        rscale = rscale.max(finite(b - a)?.abs());
    }
    if rscale == 0.0 {
        return Ok(0.0);
    }
    let (mut bnorm, mut rnorm) = (0.0_f64, 0.0_f64);
    for (&b, &a) in rhs.iter().zip(action) {
        bnorm = bnorm.hypot(if bscale == 0.0 { 0.0 } else { b / bscale });
        rnorm = rnorm.hypot((b - a) / rscale);
    }
    finite(if bscale == 0.0 {
        rscale * rnorm
    } else {
        (rscale / bscale) * (rnorm / bnorm)
    })
}

fn magnitude(v: [f64; 3]) -> f64 {
    v[0].hypot(v[1]).hypot(v[2])
}
