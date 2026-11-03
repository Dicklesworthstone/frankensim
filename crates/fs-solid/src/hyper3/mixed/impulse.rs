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
use fs_sparse::{Coo, Csr, DirectOrdering, Inertia, LdltFactor, LdltOptions, SymbolicLdlt};

use super::element::PreparedMixedMesh;
use super::{MixedUpError, checkpoint, inertia_error, invalid};
use crate::inertia::{
    ReferenceMass, SolidImpulseReport, SolidMomentum, admit_reference_impulses, impulse_moments,
};

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
    /// Actual native KKT factor inertia for a constrained impulse. None when
    /// the original checked mass solver supplies an unconstrained impulse.
    pub kkt_inertia: Option<Inertia>,
}

/// One fixed-geometry, supported mass/pressure operator for any number of loads.
/// The caller retains the accepted geometry unchanged while contact mobility
/// and the event share this factor. Reference mass and supports are borrowed;
/// no material history, background velocity or point trace is copied here.
#[derive(Debug)]
pub(super) struct ConstrainedImpulseResponse<'a> {
    mesh: &'a PreparedMixedMesh,
    mass: &'a ReferenceMass,
    held_velocity_dofs: &'a [usize],
    free: Vec<usize>,
    gradient: Csr,
    transpose: Csr,
    pressure_rows: Vec<Option<usize>>,
    matrix: Csr,
    factor: LdltFactor,
    controls: MixedImpulseControls,
    maximum_rate: f64,
}

struct VelocityResponse {
    increment_mps: Vec<[f64; 3]>,
    pressure_impulses_pa_s: Vec<f64>,
    solve: SolveReport,
}

impl<'a> ConstrainedImpulseResponse<'a> {
    /// Assemble the original endpoint volume derivative and factor its scaled
    /// KKT system once. This is the same owner used for accepted nodal events.
    #[allow(clippy::too_many_arguments, clippy::too_many_lines)]
    pub(super) fn prepare(
        mesh: &'a PreparedMixedMesh,
        mass: &'a ReferenceMass,
        held_velocity_dofs: &'a [usize],
        pressure_gauge_node: Option<usize>,
        controls: MixedImpulseControls,
        mut endpoint_cofactor: impl FnMut(usize, usize) -> Result<[f64; 9], MixedUpError>,
        cx: &Cx<'_>,
    ) -> Result<Self, MixedUpError> {
        checkpoint(cx)?;
        let nu = mass.ndofs();
        let np = mesh.pressure_nodes.len();
        if mesh.nodes_count.checked_mul(3) != Some(nu)
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
        // Admit supports and the shared response/accounting workspace before any
        // KKT allocation. Each load later uses this same compensated admission.
        let (held, _) = admit_reference_impulses(nu, held_velocity_dofs, &[], controls.balance, cx)
            .map_err(inertia_error)?;
        let free: Vec<_> = (0..nu).filter(|&d| !held[d]).collect();
        if free.is_empty() {
            return Err(invalid(
                "mixed pressure impulse requires free velocity components",
            ));
        }
        let gauge = pressure_gauge_node
            .map(|node| {
                mesh.pressure_nodes.binary_search(&node).map_err(|_| {
                    invalid("pressure impulse gauge must name an original corner node")
                })
            })
            .transpose()?;
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
                if finite(sum)?.abs()
                    > 4096.0 * f64::EPSILON * finite(absolute)?.max(f64::MIN_POSITIVE)
                {
                    return Err(MixedUpError::IncompatiblePressureGauge {
                        node: mesh.pressure_nodes[pin],
                    });
                }
            }
        }
        let maximum_rate =
            finite(controls.maximum_relative_volume_rate_per_s * mesh.total_volume_m3)?;
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
        for (row, &dof) in free.iter().enumerate() {
            if row % 128 == 0 {
                checkpoint(cx)?;
            }
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
        Ok(Self {
            mesh,
            mass,
            held_velocity_dofs,
            free,
            gradient,
            transpose,
            pressure_rows,
            matrix,
            factor,
            controls,
            maximum_rate,
        })
    }

    pub(super) fn kkt_inertia(&self) -> Inertia {
        self.factor.inertia()
    }

    /// Validate the original weak rate without projecting or copying a state.
    pub(super) fn validate_velocity(
        &self,
        velocity_mps: &[[f64; 3]],
        cx: &Cx<'_>,
    ) -> Result<Vec<f64>, MixedUpError> {
        checkpoint(cx)?;
        if velocity_mps.len() != self.mesh.nodes_count
            || velocity_mps.iter().flatten().any(|v| !v.is_finite())
        {
            return Err(invalid("finite complete mixed impulse velocity required"));
        }
        volume_rates(&self.gradient, velocity_mps, self.maximum_rate, cx)
    }

    /// Solve a homogeneous tangent response B delta-v = 0 directly. Mobility
    /// never subtracts an accepted, possibly much larger background velocity.
    /// Held increments are zero even when the physical supports are moving.
    pub(super) fn free_velocity_response(
        &self,
        nodal_impulses: &[(usize, f64)],
        cx: &Cx<'_>,
    ) -> Result<(Vec<[f64; 3]>, SolveReport), MixedUpError> {
        checkpoint(cx)?;
        let (_, external) = admit_reference_impulses(
            self.mass.ndofs(),
            self.held_velocity_dofs,
            nodal_impulses,
            self.controls.balance,
            cx,
        )
        .map_err(inertia_error)?;
        let response = self.solve_response(&external, None, cx)?;
        self.validate_velocity(&response.increment_mps, cx)?;
        let momentum = self
            .mass
            .action(response.increment_mps.as_flattened(), cx)
            .map_err(inertia_error)?;
        let pressure = self.pressure_impulse(&response.pressure_impulses_pa_s, cx)?;
        let (mut residual_work, mut pressure_work) = (0.0, 0.0);
        for (row, &dof) in self.free.iter().enumerate() {
            if row % 256 == 0 {
                checkpoint(cx)?;
            }
            // Match the original event owner's effective nodal RHS. Resolve
            // opposing external/pressure loads before subtracting a small
            // represented momentum increment from either large term.
            let effective = finite(external[dof] + pressure[dof])?;
            let residual = finite(momentum[dof] - effective)?;
            if residual.abs() > self.controls.balance.linear_tolerance_kg_m_s {
                return Err(invalid(
                    "mixed impulse response exceeds free momentum allowance",
                ));
            }
            let average = 0.5 * response.increment_mps.as_flattened()[dof];
            residual_work = residual.mul_add(average, residual_work);
            pressure_work = pressure[dof].mul_add(average, pressure_work);
        }
        if finite(residual_work)?.abs() > self.controls.balance.work_tolerance_j
            || finite(pressure_work)?.abs() > self.controls.balance.work_tolerance_j
        {
            return Err(invalid(
                "mixed impulse response violates zero-work allowance",
            ));
        }
        checkpoint(cx)?;
        Ok((response.increment_mps, response.solve))
    }

    fn solve_response(
        &self,
        external: &[f64],
        before_rate: Option<&[f64]>,
        cx: &Cx<'_>,
    ) -> Result<VelocityResponse, MixedUpError> {
        checkpoint(cx)?;
        let controls = self.controls;
        let vs = controls.velocity_scale_m_s;
        let ps = controls.pressure_impulse_scale_pa_s;
        let ws = controls.work_scale_j;
        let n = self.matrix.nrows();
        let mut rhs = vec![0.0; n];
        for (row, &dof) in self.free.iter().enumerate() {
            if row % 256 == 0 {
                checkpoint(cx)?;
            }
            rhs[row] = finite((vs / ws) * external[dof])?;
        }
        if let Some(before_rate) = before_rate {
            for (p, row) in self.pressure_rows.iter().enumerate() {
                if let Some(row) = row {
                    rhs[*row] = finite(-(ps / ws) * before_rate[p])?;
                }
            }
        }
        let solution = self.factor.solve(&rhs);
        checkpoint(cx)?;
        let mut action = vec![0.0; n];
        self.matrix.spmv(&solution, &mut action);
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
        let mut increment_mps = vec![[0.0; 3]; self.mesh.nodes_count];
        for (row, &dof) in self.free.iter().enumerate() {
            if row % 256 == 0 {
                checkpoint(cx)?;
            }
            increment_mps[dof / 3][dof % 3] = finite(vs * solution[row])?;
        }
        let pressure_impulses_pa_s = self
            .pressure_rows
            .iter()
            .map(|row| finite(row.map_or(0.0, |r| ps * solution[r])))
            .collect::<Result<Vec<_>, _>>()?;
        checkpoint(cx)?;
        Ok(VelocityResponse {
            increment_mps,
            pressure_impulses_pa_s,
            solve,
        })
    }

    fn pressure_impulse(
        &self,
        multipliers_pa_s: &[f64],
        cx: &Cx<'_>,
    ) -> Result<Vec<f64>, MixedUpError> {
        checkpoint(cx)?;
        let mut impulse = vec![0.0; self.mass.ndofs()];
        self.transpose.spmv(multipliers_pa_s, &mut impulse);
        for (dof, value) in impulse.iter_mut().enumerate() {
            if dof % 256 == 0 {
                checkpoint(cx)?;
            }
            *value = finite(-*value)?;
        }
        Ok(impulse)
    }

    /// Apply an event through this already factored original constraint, then
    /// audit the actual rounded velocities with the existing mass owner.
    #[allow(clippy::too_many_lines)]
    pub(super) fn impulse_candidate(
        &self,
        positions_m: &[[f64; 3]],
        velocity_mps: &[[f64; 3]],
        nodal_impulses: &[(usize, f64)],
        origin_m: [f64; 3],
        cx: &Cx<'_>,
    ) -> Result<(Vec<[f64; 3]>, MixedImpulseReport), MixedUpError> {
        checkpoint(cx)?;
        if positions_m.len() != self.mesh.nodes_count
            || positions_m.iter().flatten().any(|v| !v.is_finite())
            || origin_m.iter().any(|v| !v.is_finite())
        {
            return Err(invalid(
                "finite complete mixed impulse positions and origin required",
            ));
        }
        let before_rate = self.validate_velocity(velocity_mps, cx)?;
        let nu = self.mass.ndofs();
        let controls = self.controls;
        let (_, external) = admit_reference_impulses(
            nu,
            self.held_velocity_dofs,
            nodal_impulses,
            controls.balance,
            cx,
        )
        .map_err(inertia_error)?;
        let response = self.solve_response(&external, Some(&before_rate), cx)?;
        let mut velocity = velocity_mps.to_vec();
        for (row, &dof) in self.free.iter().enumerate() {
            if row % 256 == 0 {
                checkpoint(cx)?;
            }
            let change = response.increment_mps.as_flattened()[dof];
            if change != 0.0 {
                velocity[dof / 3][dof % 3] = finite(velocity[dof / 3][dof % 3] + change)?;
            }
        }
        let pressure_impulse = self.pressure_impulse(&response.pressure_impulses_pa_s, cx)?;
        let mut effective = Vec::with_capacity(nu);
        let (mut external_work_j, mut pressure_work_j) = (0.0, 0.0);
        for dof in 0..nu {
            if dof % 256 == 0 {
                checkpoint(cx)?;
            }
            let average =
                0.5 * velocity_mps.as_flattened()[dof] + 0.5 * velocity.as_flattened()[dof];
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
        let after_rate = self.validate_velocity(&velocity, cx)?;
        let (velocity, momentum) = self
            .mass
            .finish_impulse_candidate(
                positions_m,
                velocity_mps,
                self.held_velocity_dofs,
                &effective,
                origin_m,
                controls.balance,
                velocity,
                response.solve,
                cx,
            )
            .map_err(inertia_error)?;
        let report = MixedImpulseReport {
            momentum,
            external: impulse_moments(positions_m, &external, origin_m, cx)
                .map_err(inertia_error)?,
            pressure,
            external_work_j: finite(external_work_j)?,
            pressure_work_j,
            pressure_nodes: self.mesh.pressure_nodes.clone(),
            pressure_impulses_pa_s: response.pressure_impulses_pa_s,
            volume_rate_before_m3_s: before_rate,
            volume_rate_after_m3_s: after_rate,
            kkt_inertia: Some(self.kkt_inertia()),
        };
        checkpoint(cx)?;
        Ok((velocity, report))
    }
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

#[cfg(test)]
mod tests;
