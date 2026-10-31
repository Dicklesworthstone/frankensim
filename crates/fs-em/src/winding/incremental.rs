//! Differential winding inductance at a source-admitted nonlinear DC bias.
//!
//! The magnetic owner first solves the actual operating field. Its SAME
//! constitutive Hessian K_t then solves K_t dA/di_j = load(J_unit,j), with
//! homogeneous perturbations of the original prescribed trace. The reciprocal
//! block L_ij = (dA/di_i)^T K_t (dA/di_j) is checked against integral
//! (dA/di_j) J_unit,i dV. This is d(lambda)/di, never lambda/i or 2W/i^2.
//!
//! One sparse tangent is reused for all windings; the field is never densified.
//! Only the at-most-16 terminal block is inverted. Temperature, geometry,
//! remanence and external boundary values are held fixed. This local derivative
//! is not a finite-amplitude constitutive law or a frequency-dependent impedance.
//! The operating-point/tangent-energy construction follows the
//! [ANSYS electromagnetic theory reference, section 5.5](https://ansyshelp.ansys.com/public/Views/Secured/corp/v251/en/ans_thry/thy_emg6.html).

use super::{
    MAX_WINDINGS, WindingConfig, WindingError, admit_profiles, finite, inverse, multiply,
    pairing_on_mesh,
};
use crate::planar::{
    MagneticError, PlanarMagnetostatics, gradients,
    nonlinear::{
        IsotropicMagneticMaterial, NonlinearMagneticControls, NonlinearMagneticError,
        PlanarNonlinearField, PlanarNonlinearMagnetostatics,
    },
    poll,
};
use fs_solver::{CheckedCgError, CsrOp, SolveReport, checked_cg};
use fs_sparse::{Csr, precond::IdentityPrecond};

/// Numerical budgets plus an explicit exclusion band around source corners.
#[derive(Clone, Copy, Debug)]
pub struct IncrementalWindingConfig {
    /// Existing terminal count, closure, reciprocity and inverse controls.
    /// Here max_basis_values bounds retained profile/response/operating arrays:
    /// windings*(vertices+3*cells) + 3*vertices+14*cells. The fixed-size terminal
    /// matrices and the existing solvers' separately bounded work are additional.
    pub winding: WindingConfig,
    /// Minimum distance [T] from any nonsmooth H(B) knot and the admitted
    /// maximum |B|. Zero still refuses exact corners and range endpoints.
    /// This caller-selected margin is not a bound on operating-field error.
    pub minimum_smooth_field_margin_t: f64,
}

/// Failed physical, numerical or work admission publishes no partial basis.
#[derive(Debug)]
pub enum IncrementalWindingError {
    /// Original winding closure, allocation, reciprocity or inverse refusal.
    Winding(WindingError),
    /// Original material, geometry, nonlinear/linear solve or cancellation error.
    Magnetic(NonlinearMagneticError),
    /// A piecewise source has no unique two-sided tangent at this breakpoint.
    Nondifferentiable {
        /// Original mesh triangle.
        triangle: usize,
        /// Actual operating flux-density magnitude [T].
        observed_t: f64,
        /// Source breakpoint within the declared exclusion band [T].
        breakpoint_t: f64,
    },
    /// The operating point lacks a two-sided neighborhood in its source domain.
    SourceBoundary {
        /// Original mesh triangle whose source neighborhood is unavailable.
        triangle: usize,
    },
}
impl core::fmt::Display for IncrementalWindingError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "incremental winding inductance: {self:?}")
    }
}
impl std::error::Error for IncrementalWindingError {}
impl From<WindingError> for IncrementalWindingError {
    fn from(error: WindingError) -> Self {
        Self::Winding(error)
    }
}
impl From<NonlinearMagneticError> for IncrementalWindingError {
    fn from(error: NonlinearMagneticError) -> Self {
        Self::Magnetic(error)
    }
}
impl From<MagneticError> for IncrementalWindingError {
    fn from(error: MagneticError) -> Self {
        Self::Magnetic(error.into())
    }
}
type Result<T> = core::result::Result<T, IncrementalWindingError>;

/// An infinitesimal unit-current derivative, not a realizable 1-A field.
pub struct IncrementalWindingResponse {
    potential_per_current_h_m: Vec<f64>,
    flux_density_per_current_t_a: Vec<[f64; 2]>,
    report: SolveReport,
}
impl IncrementalWindingResponse {
    /// dA_z/di [H/m], in original mesh vertex order.
    pub fn potential_per_current_h_m(&self) -> &[f64] {
        &self.potential_per_current_h_m
    }
    /// dB/di [T/A], without averaging interfaces or extrapolating a source law.
    pub fn flux_density_per_current_t_a(&self) -> &[[f64; 2]] {
        &self.flux_density_per_current_t_a
    }
    /// Recomputed true Euclidean residual for the linearized field solve.
    pub fn report(&self) -> &SolveReport {
        &self.report
    }
}

/// A complete solved operating point and its local reciprocal terminal response.
pub struct PlanarIncrementalWindingBasis<'a> {
    operating_field: PlanarNonlinearField<'a>,
    operating_currents_a: Vec<f64>,
    operating_flux_linkages_wb: Vec<f64>,
    boundary_potential_wb_m: Vec<(usize, f64)>,
    profiles_per_m2: Vec<Vec<f64>>,
    responses: Vec<IncrementalWindingResponse>,
    active_inductance_h: Vec<f64>,
    inductance_h: Vec<f64>,
    inverse_inductance_per_h: Vec<f64>,
    end_inductance_h: Vec<f64>,
    net_turns_by_component: Vec<Vec<f64>>,
    maximum_pairing_defect_h: f64,
    inverse_residual: f64,
}

impl<'a> PlanarIncrementalWindingBasis<'a> {
    /// Solve the actual fixed-state nonlinear bias before extracting d(lambda)/di.
    ///
    /// Currents excite the original signed profiles; fixed recoil magnets and
    /// nonzero prescribed A remain physical bias sources. Unit derivatives have
    /// zero trace on EXACTLY those prescribed vertices. Both operating and
    /// derivative solves must pass their original true-residual targets.
    /// Source corners and envelope edges refuse instead of assigning a misleading
    /// unique inductance. Explicit nonnegative end self-inductances remain linear.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        problem: &'a PlanarNonlinearMagnetostatics,
        profiles_per_m2: Vec<Vec<f64>>,
        boundary_potential_wb_m: &[(usize, f64)],
        operating_currents_a: Vec<f64>,
        end_inductance_h: Vec<f64>,
        initial_potential_wb_m: &[f64],
        config: IncrementalWindingConfig,
        controls: NonlinearMagneticControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self> {
        poll(cancelled)?;
        let cfg = config.winding;
        let count = profiles_per_m2.len();
        let mesh = problem.mesh();
        let vertices = mesh.vertices().len();
        let cells = mesh.faces().len();
        let values = cells
            .checked_mul(3)
            .and_then(|n| n.checked_add(vertices))
            .and_then(|n| n.checked_mul(count))
            .and_then(|n| vertices.checked_mul(3).and_then(|v| n.checked_add(v)))
            .and_then(|n| cells.checked_mul(14).and_then(|c| n.checked_add(c)));
        if count == 0
            || count > cfg.max_windings
            || cfg.max_windings > MAX_WINDINGS
            || values.is_none_or(|n| n > cfg.max_basis_values)
            || operating_currents_a.len() != count
            || operating_currents_a.iter().any(|i| !i.is_finite())
            || end_inductance_h.len() != count
            || end_inductance_h.iter().any(|l| !l.is_finite() || *l < 0.0)
            || !cfg.pairing_tolerance_h.is_finite()
            || cfg.pairing_tolerance_h < 0.0
            || !cfg.current_balance_tolerance_turns.is_finite()
            || cfg.current_balance_tolerance_turns < 0.0
            || !cfg.inverse_residual_tolerance.is_finite()
            || cfg.inverse_residual_tolerance <= 0.0
            || cfg.inverse_residual_tolerance >= 1.0
            || !config.minimum_smooth_field_margin_t.is_finite()
            || config.minimum_smooth_field_margin_t < 0.0
        {
            return Err(
                WindingError::Invalid("invalid incremental winding data or budgets").into(),
            );
        }
        let net_turns_by_component = admit_profiles(
            mesh,
            &profiles_per_m2,
            problem.thickness_m(),
            cfg.current_balance_tolerance_turns,
            cancelled,
        )?;
        let prescribed =
            PlanarMagnetostatics::admit_mesh_boundary(mesh, boundary_potential_wb_m, cancelled)?;
        let mut current_density = vec![0.0; cells];
        for (profile, &current) in profiles_per_m2.iter().zip(&operating_currents_a) {
            for (density, &unit) in current_density.iter_mut().zip(profile) {
                poll(cancelled)?;
                *density = finite(unit.mul_add(current, *density))?;
            }
        }
        let operating_field = problem.solve(
            &current_density,
            boundary_potential_wb_m,
            initial_potential_wb_m,
            controls,
            cancelled,
        )?;
        admit_smooth_operating_point(
            &operating_field,
            config.minimum_smooth_field_margin_t,
            cancelled,
        )?;
        let (_, tangent) =
            problem.magnetic_linearization(operating_field.potential_wb_m(), cancelled)?;
        let free: Vec<_> = (0..vertices).filter(|&v| prescribed[v].is_none()).collect();
        let mut index = vec![usize::MAX; vertices];
        for (i, &v) in free.iter().enumerate() {
            index[v] = i;
        }
        let mut row_ptr = vec![0];
        let mut columns = Vec::new();
        let mut values = Vec::new();
        for &v in &free {
            poll(cancelled)?;
            let (cols, vals) = tangent.row(v);
            for (&w, &value) in cols.iter().zip(vals) {
                if index[w] != usize::MAX {
                    columns.push(index[w]);
                    values.push(finite(value)?);
                }
            }
            row_ptr.push(columns.len());
        }
        let op = CsrOp::symmetric(Csr::from_parts(
            free.len(),
            free.len(),
            row_ptr,
            columns,
            values,
        ));
        let mut responses = Vec::with_capacity(count);
        let mut operating_flux_linkages_wb = Vec::with_capacity(count);
        for (winding, profile) in profiles_per_m2.iter().enumerate() {
            poll(cancelled)?;
            let mut rhs = vec![0.0; vertices];
            for (cell, face) in mesh.faces().iter().enumerate() {
                poll(cancelled)?;
                let load = finite(profile[cell] * (mesh.face_measure(cell).unwrap() / 3.0))?;
                for &v in face {
                    rhs[v as usize] = finite(rhs[v as usize] + load)?;
                }
            }
            let reduced_rhs: Vec<_> = free.iter().map(|&v| rhs[v]).collect();
            let result = checked_cg(&op, &IdentityPrecond, &reduced_rhs, controls.linear, |_| {
                poll(cancelled)
            })
            .map_err(|error| match error {
                CheckedCgError::InvalidInput(why) => MagneticError::Invalid(why),
                CheckedCgError::Interrupted(error) => error,
            })?;
            if !result.report.converged_euclidean() {
                return Err(MagneticError::Nonconvergence(result.report).into());
            }
            let mut potential = vec![0.0; vertices];
            for (&v, &value) in free.iter().zip(&result.x) {
                potential[v] = value;
            }
            let mut flux_density = Vec::with_capacity(cells);
            for face in mesh.faces() {
                poll(cancelled)?;
                let shape = gradients(face.map(|v| mesh.vertices()[v as usize]))?;
                let a0 = potential[face[0] as usize];
                let grad: [f64; 2] = [0, 1].map(|axis| {
                    (potential[face[1] as usize] - a0).mul_add(
                        shape[1][axis],
                        (potential[face[2] as usize] - a0) * shape[2][axis],
                    )
                });
                flux_density.push([finite(grad[1])?, finite(-grad[0])?]);
            }
            responses.push(IncrementalWindingResponse {
                potential_per_current_h_m: potential,
                flux_density_per_current_t_a: flux_density,
                report: result.report,
            });
            let active =
                pairing_on_mesh(mesh, operating_field.potential_wb_m(), profile, cancelled)?;
            operating_flux_linkages_wb.push(finite(
                end_inductance_h[winding].mul_add(operating_currents_a[winding], active),
            )?);
        }
        let mut active = vec![0.0; count * count];
        let mut applied = vec![0.0; vertices];
        for j in 0..count {
            poll(cancelled)?;
            tangent.spmv(responses[j].potential_per_current_h_m(), &mut applied);
            for i in 0..=j {
                let mut gram = 0.0;
                for (&value, &effort) in responses[i]
                    .potential_per_current_h_m()
                    .iter()
                    .zip(&applied)
                {
                    poll(cancelled)?;
                    gram = finite(value.mul_add(effort, gram))?;
                }
                active[i * count + j] = gram;
                active[j * count + i] = gram;
            }
            if active[j * count + j] <= 0.0 {
                return Err(WindingError::Invalid(
                    "winding has no resolved differential magnetic energy",
                )
                .into());
            }
        }
        let mut maximum_pairing_defect_h = 0.0_f64;
        for (j, response) in responses.iter().enumerate() {
            for (i, profile) in profiles_per_m2.iter().enumerate() {
                let flux = pairing_on_mesh(
                    mesh,
                    response.potential_per_current_h_m(),
                    profile,
                    cancelled,
                )?;
                let defect = finite(flux - active[i * count + j])?.abs();
                if defect > cfg.pairing_tolerance_h {
                    return Err(WindingError::Pairing {
                        source: j,
                        receiver: i,
                        flux_h: flux,
                        energy_h: active[i * count + j],
                    }
                    .into());
                }
                maximum_pairing_defect_h = maximum_pairing_defect_h.max(defect);
            }
        }
        let mut inductance_h = active.clone();
        for i in 0..count {
            inductance_h[i * count + i] = finite(active[i * count + i] + end_inductance_h[i])?;
        }
        let (inverse_inductance_per_h, inverse_residual) = inverse(
            &inductance_h,
            count,
            cfg.inverse_residual_tolerance,
            cancelled,
        )?;
        poll(cancelled)?;
        Ok(Self {
            operating_field,
            operating_currents_a,
            operating_flux_linkages_wb,
            boundary_potential_wb_m: boundary_potential_wb_m.to_vec(),
            profiles_per_m2,
            responses,
            active_inductance_h: active,
            inductance_h,
            inverse_inductance_per_h,
            end_inductance_h,
            net_turns_by_component,
            maximum_pairing_defect_h,
            inverse_residual,
        })
    }

    /// Converged DC field with the original geometry, material laws and report.
    pub fn operating_field(&self) -> &PlanarNonlinearField<'a> {
        &self.operating_field
    }
    /// Signed DC currents [A] in the original winding order.
    pub fn operating_currents_a(&self) -> &[f64] {
        &self.operating_currents_a
    }
    /// Actual source-paired flux, including remanent/boundary bias and end flux.
    pub fn operating_flux_linkages_wb(&self) -> &[f64] {
        &self.operating_flux_linkages_wb
    }
    /// Original prescribed vertex trace, held fixed during differentiation.
    pub fn boundary_potential_wb_m(&self) -> &[(usize, f64)] {
        &self.boundary_potential_wb_m
    }
    /// Original signed turn/current densities [1/m2], in winding/cell order.
    pub fn profiles_per_m2(&self) -> &[Vec<f64>] {
        &self.profiles_per_m2
    }
    /// Derivative fields and true-residual reports, in winding order.
    pub fn unit_responses(&self) -> &[IncrementalWindingResponse] {
        &self.responses
    }
    /// Number of terminal-current coordinates.
    pub fn winding_count(&self) -> usize {
        self.responses.len()
    }
    /// Symmetric row-major active-section d(lambda_i)/d(i_j) [H].
    pub fn active_inductance_h(&self) -> &[f64] {
        &self.active_inductance_h
    }
    /// Differential inductance with the explicit linear end-storage diagonals.
    pub fn inductance_h(&self) -> &[f64] {
        &self.inductance_h
    }
    /// Admitted symmetric inverse of the full differential block [1/H].
    pub fn inverse_inductance_per_h(&self) -> &[f64] {
        &self.inverse_inductance_per_h
    }
    /// Caller-declared end/leakage self-inductances [H], not resolved fields.
    pub fn end_inductance_h(&self) -> &[f64] {
        &self.end_inductance_h
    }
    /// Signed current-closure residual [turns] per winding and mesh component.
    pub fn net_turns_by_component(&self) -> &[Vec<f64>] {
        &self.net_turns_by_component
    }
    /// Largest source-pairing versus Hessian-energy discrepancy [H].
    pub fn maximum_pairing_defect_h(&self) -> f64 {
        self.maximum_pairing_defect_h
    }
    /// Measured infinity-norm residual of L_diff times its inverse minus I.
    pub fn inverse_residual(&self) -> f64 {
        self.inverse_residual
    }

    /// Local tangent delta_lambda = L_diff delta_i; not a finite-state update.
    pub fn linearized_flux_increment_wb(&self, current_increment_a: &[f64]) -> Result<Vec<f64>> {
        Ok(multiply(
            &self.inductance_h,
            current_increment_a,
            self.winding_count(),
        )?)
    }
    /// Inverse local tangent, retaining all self and mutual winding terms.
    pub fn linearized_current_increment_a(&self, flux_increment_wb: &[f64]) -> Result<Vec<f64>> {
        Ok(multiply(
            &self.inverse_inductance_per_h,
            flux_increment_wb,
            self.winding_count(),
        )?)
    }
    /// Inductive voltage drop L_diff di/dt at this fixed operating point [V].
    /// No resistive, moving-boundary, thermal or hysteretic voltage is implied.
    pub fn voltage_from_current_rate_v(&self, current_rate_a_s: &[f64]) -> Result<Vec<f64>> {
        self.linearized_flux_increment_wb(current_rate_a_s)
    }
}

fn admit_smooth_operating_point(
    field: &PlanarNonlinearField<'_>,
    margin_t: f64,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<()> {
    for (triangle, material) in field.problem().materials().iter().enumerate() {
        poll(cancelled)?;
        match material {
            IsotropicMagneticMaterial::Reversible(law) => {
                let b = field.magnetic_flux_density_t()[triangle];
                let magnitude = b[0].hypot(b[1]);
                if law.maximum_flux_density_t() - magnitude <= margin_t {
                    return Err(IncrementalWindingError::SourceBoundary { triangle });
                }
                // Search only source breakpoints inside the exclusion band;
                // long measured curves must not cost a full scan per cell.
                let knots = &law.source_path().knots;
                let first = knots
                    .partition_point(|(b, _)| *b < magnitude - margin_t)
                    .max(1);
                let end = knots
                    .partition_point(|(b, _)| *b <= magnitude + margin_t)
                    .min(knots.len() - 1);
                for knot in first..end {
                    let left = (knots[knot].1.evidence.value.value
                        - knots[knot - 1].1.evidence.value.value)
                        / (knots[knot].0 - knots[knot - 1].0);
                    let right = (knots[knot + 1].1.evidence.value.value
                        - knots[knot].1.evidence.value.value)
                        / (knots[knot + 1].0 - knots[knot].0);
                    if left != right {
                        return Err(IncrementalWindingError::Nondifferentiable {
                            triangle,
                            observed_t: magnitude,
                            breakpoint_t: knots[knot].0,
                        });
                    }
                }
            }
            IsotropicMagneticMaterial::Recoil(law) => {
                let h = field.magnetic_field_a_m()[triangle];
                let direction = law.direction();
                let parallel = direction[0].mul_add(h[0], direction[1] * h[1]);
                let transverse = direction[0].mul_add(h[1], -direction[1] * h[0]).abs();
                let source = law.material();
                let range = source.parallel_range_a_m();
                if parallel <= range[0]
                    || parallel >= range[1]
                    || transverse >= source.maximum_transverse_a_m()
                {
                    return Err(IncrementalWindingError::SourceBoundary { triangle });
                }
            }
            IsotropicMagneticMaterial::Linear(_) => {}
        }
    }
    Ok(())
}
