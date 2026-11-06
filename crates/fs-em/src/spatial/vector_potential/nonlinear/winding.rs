//! Local differential inductance of source-resolved 3D Whitney windings.
//!
//! First solve the actual nonlinear bias, then reuse its one sparse magnetic
//! Hessian and one coexact perturbation topology for every signed unit current.
//! L_ij = d(lambda_i)/di_j = delta_A_i^T K_t delta_A_j is independently checked
//! against the original Whitney current pairing J_unit_i^T delta_A_j.
//! Neither lambda/i nor 2W/i^2 is a differential inductance at saturation.
//!
//! Perturbation traces are zero on the original prescribed edges. Temperature,
//! geometry, remanence and physical boundary traces stay fixed. Derivative
//! fields are not realizable one-ampere magnet states; only the solved bias
//! must lie in a two-sided smooth source neighborhood. No finite-amplitude
//! circuit evolution, hysteresis, eddy loss or frequency response is implied.

use super::super::{
    VectorMagneticBoundary, VectorPotentialFace, VectorTopology, apply, compatible_solve,
    gauge_error, max_difference, project, winding::SpatialWindingConfig,
};
use super::{
    NonlinearVectorControls, NonlinearVectorError, SpatialNonlinearVectorField,
    SpatialNonlinearVectorMagnetostatics,
};
use crate::spatial::{
    SpatialMagneticError, cross, dot, finite, invalid, nonlinear::SpatialNonlinearMaterial, norm,
    poll, sub,
};
use crate::winding::end::{EndWindingError, EndWindingInductance};
use fs_feec::{
    assembly::stiffness,
    weighted::{CellWeight, WeightedError, weighted_mass_matrix},
};
use fs_qty::Dims;
use fs_solver::SolveReport;
use fs_sparse::Csr;

/// Existing winding bounds plus a declared exclusion band around source corners.
#[derive(Clone, Copy, Debug)]
pub struct SpatialIncrementalWindingConfig {
    /// Here max_basis_values bounds retained profile, response and bias arrays:
    /// n*(2*edges+faces+6*cells) + 3*edges+4*faces+24*cells+2*n.
    /// Small terminal matrices and the original bounded field/gauge work are
    /// additional. The physical terminal count remains at most sixteen.
    pub winding: SpatialWindingConfig,
    /// Distance [T] from nonsmooth H(B) knots and the maximum admitted |B|.
    /// Zero still excludes exact corners/endpoints. Recoil envelopes require
    /// a strict interior operating point. This is not a field-error bound.
    pub minimum_smooth_field_margin_t: f64,
}

/// Failed admission or numerical work publishes no partial terminal basis.
#[derive(Debug)]
pub enum SpatialIncrementalWindingError {
    Magnetic(NonlinearVectorError),
    End(EndWindingError),
    Pairing {
        source: usize,
        receiver: usize,
        flux_h: f64,
        energy_h: f64,
    },
    Nondifferentiable {
        cell: usize,
        observed_t: f64,
        breakpoint_t: f64,
    },
    SourceBoundary {
        cell: usize,
    },
}
impl From<NonlinearVectorError> for SpatialIncrementalWindingError {
    fn from(error: NonlinearVectorError) -> Self {
        Self::Magnetic(error)
    }
}
impl From<SpatialMagneticError> for SpatialIncrementalWindingError {
    fn from(error: SpatialMagneticError) -> Self {
        Self::Magnetic(error.into())
    }
}
impl From<EndWindingError> for SpatialIncrementalWindingError {
    fn from(error: EndWindingError) -> Self {
        Self::End(error)
    }
}
impl core::fmt::Display for SpatialIncrementalWindingError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "spatial incremental winding: {self:?}")
    }
}
impl std::error::Error for SpatialIncrementalWindingError {}
type Result<T> = core::result::Result<T, SpatialIncrementalWindingError>;

/// An infinitesimal response; its material domain is the solved operating point.
pub struct SpatialIncrementalWindingResponse {
    potential_per_current_h: Vec<f64>,
    face_flux_per_current_h: Vec<f64>,
    flux_density_per_current_t_a: Vec<[f64; 3]>,
    report: Option<SolveReport>,
    /// Maximum original free-edge tangent-equation residual, normalized to 1 A.
    pub maximum_equation_residual: f64,
    /// Full curl-kernel compatibility defect of the original unit-current load.
    pub source_compatibility_defect: f64,
    /// Worst gauge, div-curl or recovered normal-flux defect [H].
    pub maximum_flux_defect_h: f64,
}
impl SpatialIncrementalWindingResponse {
    /// d(edge-integral A)/di [H], in original mesh edge order.
    pub fn potential_per_current_h(&self) -> &[f64] {
        &self.potential_per_current_h
    }
    /// d(face-integral B)/di [H], in original mesh face order.
    pub fn face_flux_per_current_h(&self) -> &[f64] {
        &self.face_flux_per_current_h
    }
    /// dB/di [T/A], without averaging across material interfaces.
    pub fn flux_density_per_current_t_a(&self) -> &[[f64; 3]] {
        &self.flux_density_per_current_t_a
    }
    /// Native tangent solve report; None only for an empty free-edge space.
    /// The full original-source residual after gauging is checked separately.
    pub fn report(&self) -> Option<&SolveReport> {
        self.report.as_ref()
    }
}

/// Solved nonlinear bias and its local reciprocal self/mutual response.
/// The basis does not invert L: dependent winding coordinates remain valid
/// field observations and do not acquire invented leakage or an inverse.
pub struct SpatialIncrementalWindingBasis<'m> {
    operating_field: SpatialNonlinearVectorField<'m>,
    operating_currents_a: Vec<f64>,
    operating_flux_linkages_wb: Vec<f64>,
    boundary: VectorMagneticBoundary,
    profiles_per_m2: Vec<Vec<[f64; 3]>>,
    loads_per_a: Vec<Vec<f64>>,
    responses: Vec<SpatialIncrementalWindingResponse>,
    active_inductance_h: Vec<f64>,
    inductance_h: Vec<f64>,
    end_inductance: EndWindingInductance,
    maximum_pairing_defect_h: f64,
}
impl<'m> SpatialIncrementalWindingBasis<'m> {
    /// Solve the actual current/trace/remanence bias before extracting d(lambda)/di.
    /// Profiles include turns and orientation and must close on every original
    /// face, including prescribed exterior traces. External series self-energy
    /// is explicit; zero neglects unmeshed ends. All field and gauge allowances
    /// are the native controls evaluated on a unit-current perturbation.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        problem: &'m SpatialNonlinearVectorMagnetostatics,
        profiles_per_m2: Vec<Vec<[f64; 3]>>,
        boundary: VectorMagneticBoundary,
        operating_currents_a: Vec<f64>,
        external_inductance_h: Vec<f64>,
        initial_potential_wb: &[f64],
        config: SpatialIncrementalWindingConfig,
        controls: NonlinearVectorControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self> {
        poll(cancelled)?;
        Self::new_with_end_modes(
            problem,
            profiles_per_m2,
            boundary,
            operating_currents_a,
            EndWindingInductance::new(external_inductance_h, vec![])?,
            initial_potential_wb,
            config,
            controls,
            cancelled,
        )
    }

    /// Extract the local differential field response with explicit passive
    /// end-region self and shared mutual flux. The authored end model must
    /// exclude flux already resolved by the field. Its original factors are
    /// retained for operating linkages and current-rate response, including
    /// cancelling modes that cannot be recovered from rounded matrix entries.
    #[allow(clippy::too_many_arguments)]
    pub fn new_with_end_modes(
        problem: &'m SpatialNonlinearVectorMagnetostatics,
        profiles_per_m2: Vec<Vec<[f64; 3]>>,
        boundary: VectorMagneticBoundary,
        operating_currents_a: Vec<f64>,
        end_inductance: EndWindingInductance,
        initial_potential_wb: &[f64],
        config: SpatialIncrementalWindingConfig,
        controls: NonlinearVectorControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self> {
        poll(cancelled)?;
        let n = profiles_per_m2.len();
        let mesh = problem.mesh();
        let edges = mesh.edges.len();
        let cells = mesh.tets.len();
        let faces = mesh.faces.len();
        let values = edges
            .checked_mul(2)
            .and_then(|v| v.checked_add(faces))
            .and_then(|v| cells.checked_mul(6).and_then(|c| v.checked_add(c)))
            .and_then(|v| v.checked_mul(n))
            .and_then(|v| edges.checked_mul(3).and_then(|e| v.checked_add(e)))
            .and_then(|v| faces.checked_mul(4).and_then(|f| v.checked_add(f)))
            .and_then(|v| cells.checked_mul(24).and_then(|c| v.checked_add(c)))
            .and_then(|v| n.checked_mul(2).and_then(|c| v.checked_add(c)));
        if n == 0
            || n > config.winding.max_windings
            || config.winding.max_windings > 16
            || values.is_none_or(|v| v > config.winding.max_basis_values)
            || operating_currents_a.len() != n
            || operating_currents_a.iter().any(|v| !v.is_finite())
            || end_inductance.winding_count() != n
            || boundary.potential.len() > faces
            || !config.winding.pairing_tolerance_h.is_finite()
            || config.winding.pairing_tolerance_h < 0.0
            || !config.minimum_smooth_field_margin_t.is_finite()
            || config.minimum_smooth_field_margin_t < 0.0
            || profiles_per_m2.iter().any(|p| {
                p.len() != cells
                    || p.iter().flatten().any(|v| !v.is_finite())
                    || p.iter().flatten().all(|&v| v == 0.0)
            })
        {
            return Err(invalid("bounded complete incremental winding data required").into());
        }
        let mut loads_per_a = Vec::with_capacity(n);
        let mut current = vec![[0.0; 3]; cells];
        for (profile, &amperes) in profiles_per_m2.iter().zip(&operating_currents_a) {
            poll(cancelled)?;
            // Closed unit profiles, even when their operating current is zero.
            for face in &problem.geometry.faces {
                poll(cancelled)?;
                let jump = if face.cells.len() == 2 {
                    std::array::from_fn(|k| profile[face.cells[0]][k] - profile[face.cells[1]][k])
                } else {
                    profile[face.cells[0]]
                };
                if finite(face.area_m2 * dot(jump, face.outward_normal))?.abs()
                    > controls.field.current_tolerance_a
                {
                    return Err(invalid(
                        "unit winding current has an unresolved return or discontinuity",
                    )
                    .into());
                }
            }
            let mut load = vec![0.0; edges];
            for (cell, vectors) in problem.cells.iter().enumerate() {
                poll(cancelled)?;
                for k in 0..3 {
                    current[cell][k] = finite(amperes.mul_add(profile[cell][k], current[cell][k]))?;
                }
                for local in 0..6 {
                    let edge = vectors.edges[local];
                    load[edge] = finite(load[edge] + dot(profile[cell], vectors.integrals[local]))?;
                }
            }
            loads_per_a.push(load);
        }
        let operating_field = problem.solve(
            &current,
            &boundary,
            initial_potential_wb,
            controls,
            cancelled,
        )?;
        admit_smooth_operating_point(
            &operating_field,
            config.minimum_smooth_field_margin_t,
            cancelled,
        )?;
        // The original physical trace remains in the bias. Only its derivative
        // is zero, on exactly the same edges; remanence is not a perturbation RHS.
        let homogeneous = VectorMagneticBoundary {
            potential: boundary
                .potential
                .iter()
                .map(|trace| VectorPotentialFace {
                    face: trace.face,
                    edge_integrals_wb: [0.0; 3],
                })
                .collect(),
        };
        let sources = VectorTopology {
            geometry: &problem.geometry,
            cells: &problem.cells,
        }
        .prepare_sources(
            &vec![[0.0; 3]; cells],
            &homogeneous,
            controls.field,
            cancelled,
        )?;
        let mut weights = Vec::with_capacity(cells);
        for (cell, &b) in operating_field.flux_density_t().iter().enumerate() {
            poll(cancelled)?;
            weights.push(CellWeight::Tensor(problem.response(cell, b)?.tangent));
        }
        let mass = weighted_mass_matrix(
            mesh,
            &problem.geometry.geometry,
            2,
            &weights,
            Dims([-1, -1, 2, 0, 2, 0]),
            problem.limits,
            cancelled,
        )
        .map_err(|error| {
            if error == WeightedError::Cancelled {
                SpatialMagneticError::Cancelled
            } else {
                SpatialMagneticError::Assembly(error)
            }
        })?;
        let full = stiffness(&problem.curl, mass.matrix());
        let mut index = vec![usize::MAX; edges];
        for (i, &edge) in sources.free.iter().enumerate() {
            index[edge] = i;
        }
        let mut rows = vec![0];
        let mut cols = Vec::new();
        let mut vals = Vec::new();
        for &edge in &sources.free {
            poll(cancelled)?;
            let (c, v) = full.row(edge);
            for (&column, &value) in c.iter().zip(v) {
                if index[column] != usize::MAX {
                    cols.push(index[column]);
                    vals.push(finite(value)?);
                }
            }
            rows.push(cols.len());
        }
        let reduced = Csr::from_parts(sources.free.len(), sources.free.len(), rows, cols, vals);
        let mut responses = Vec::with_capacity(n);
        let mut operating_flux_linkages_wb =
            end_inductance.flux_linkages_wb(&operating_currents_a)?;
        for (winding, load) in loads_per_a.iter().enumerate() {
            poll(cancelled)?;
            let rhs: Vec<_> = sources.free.iter().map(|&edge| load[edge]).collect();
            let riesz: Vec<_> = rhs
                .iter()
                .zip(sources.projector.metric())
                .map(|(j, metric)| finite(j / metric))
                .collect::<core::result::Result<_, _>>()?;
            let (_, projected_rhs, _) =
                project(&sources.projector, &riesz, controls.field.gauge, cancelled)?;
            let source_compatibility_defect = max_difference(&rhs, &projected_rhs)?;
            if source_compatibility_defect > controls.field.current_tolerance_a {
                return Err(invalid(
                    "unit winding current is incompatible with the complete curl kernel",
                )
                .into());
            }
            let (raw, report) =
                compatible_solve(&reduced, &projected_rhs, controls.field.field, cancelled)?;
            let (gauged, _, _) =
                project(&sources.projector, &raw, controls.field.gauge, cancelled)?;
            let before = sources
                .projector
                .rhs(&raw, cancelled)
                .map_err(gauge_error)?;
            let after = sources
                .projector
                .rhs(&gauged, cancelled)
                .map_err(gauge_error)?;
            let mut maximum_flux_defect_h = max_difference(&before, &after)?;
            let mut potential = vec![0.0; edges];
            for (&edge, &value) in sources.free.iter().zip(&gauged) {
                potential[edge] = value;
            }
            let equation = apply(&full, &potential, cancelled)?;
            let mut maximum_equation_residual = 0.0_f64;
            for &edge in &sources.free {
                poll(cancelled)?;
                maximum_equation_residual =
                    maximum_equation_residual.max(finite(equation[edge] - load[edge])?.abs());
            }
            if maximum_equation_residual > controls.field.current_tolerance_a {
                return Err(invalid(
                    "gauged tangent response exceeds its original current allowance",
                )
                .into());
            }
            let face_flux = apply(&problem.curl, &potential, cancelled)?;
            for row in mesh.d2().rows {
                poll(cancelled)?;
                let divergence = finite(
                    row.iter()
                        .map(|&(face, sign)| f64::from(sign) * face_flux[face])
                        .sum(),
                )?;
                maximum_flux_defect_h = maximum_flux_defect_h.max(divergence.abs());
            }
            let mut flux_density = Vec::with_capacity(cells);
            for vectors in &problem.cells {
                poll(cancelled)?;
                let mut b = [0.0; 3];
                for local in 0..6 {
                    for k in 0..3 {
                        b[k] = finite(
                            potential[vectors.edges[local]].mul_add(vectors.curls[local][k], b[k]),
                        )?;
                    }
                }
                flux_density.push(b);
            }
            for (i, face) in problem.geometry.faces.iter().enumerate() {
                poll(cancelled)?;
                let [a, b, c] = mesh.faces[i].map(|v| problem.positions_m()[v as usize]);
                let sign = if dot(cross(sub(b, a), sub(c, a)), face.outward_normal) > 0.0 {
                    1.0
                } else {
                    -1.0
                };
                for &cell in &face.cells {
                    let flux = finite(face.area_m2 * dot(flux_density[cell], face.outward_normal))?;
                    maximum_flux_defect_h =
                        maximum_flux_defect_h.max(finite(flux - sign * face_flux[i])?.abs());
                }
            }
            if maximum_flux_defect_h > controls.field.flux_tolerance_wb {
                return Err(invalid(
                    "unit tangent field violates the native gauge or div-curl allowance",
                )
                .into());
            }
            let active = pairing(load, operating_field.potential_wb(), cancelled)?;
            operating_flux_linkages_wb[winding] =
                finite(active + operating_flux_linkages_wb[winding])?;
            responses.push(SpatialIncrementalWindingResponse {
                potential_per_current_h: potential,
                face_flux_per_current_h: face_flux,
                flux_density_per_current_t_a: flux_density,
                report,
                maximum_equation_residual,
                source_compatibility_defect,
                maximum_flux_defect_h,
            });
        }
        let mut active_inductance_h = vec![0.0; n * n];
        for source in 0..n {
            let applied = apply(
                &full,
                responses[source].potential_per_current_h(),
                cancelled,
            )?;
            for receiver in 0..=source {
                let gram = pairing(
                    responses[receiver].potential_per_current_h(),
                    &applied,
                    cancelled,
                )?;
                active_inductance_h[receiver * n + source] = gram;
                active_inductance_h[source * n + receiver] = gram;
            }
        }
        let mut maximum_pairing_defect_h = 0.0_f64;
        for (source, response) in responses.iter().enumerate() {
            for (receiver, load) in loads_per_a.iter().enumerate() {
                let flux = pairing(load, response.potential_per_current_h(), cancelled)?;
                let energy = active_inductance_h[receiver * n + source];
                let defect = finite(flux - energy)?.abs();
                if defect > config.winding.pairing_tolerance_h {
                    return Err(SpatialIncrementalWindingError::Pairing {
                        source,
                        receiver,
                        flux_h: flux,
                        energy_h: energy,
                    });
                }
                maximum_pairing_defect_h = maximum_pairing_defect_h.max(defect);
            }
        }
        let mut inductance_h = active_inductance_h.clone();
        for (total, &end) in inductance_h.iter_mut().zip(end_inductance.matrix_h()) {
            *total = finite(*total + end)?;
        }
        for i in 0..n {
            if active_inductance_h[i * n + i] < 0.0 || inductance_h[i * n + i] <= 0.0 {
                return Err(
                    invalid("positive self inductance without invented leakage required").into(),
                );
            }
        }
        poll(cancelled)?;
        Ok(Self {
            operating_field,
            operating_currents_a,
            operating_flux_linkages_wb,
            boundary,
            profiles_per_m2,
            loads_per_a,
            responses,
            active_inductance_h,
            inductance_h,
            end_inductance,
            maximum_pairing_defect_h,
        })
    }

    pub fn operating_field(&self) -> &SpatialNonlinearVectorField<'m> {
        &self.operating_field
    }
    pub fn operating_currents_a(&self) -> &[f64] {
        &self.operating_currents_a
    }
    /// Actual total linkage, including fixed remanence, trace bias and end flux.
    pub fn operating_flux_linkages_wb(&self) -> &[f64] {
        &self.operating_flux_linkages_wb
    }
    pub fn boundary(&self) -> &VectorMagneticBoundary {
        &self.boundary
    }
    pub fn profiles_per_m2(&self) -> &[Vec<[f64; 3]>] {
        &self.profiles_per_m2
    }
    pub fn unit_responses(&self) -> &[SpatialIncrementalWindingResponse] {
        &self.responses
    }
    pub fn winding_count(&self) -> usize {
        self.responses.len()
    }
    /// Row-major differential block from the actual meshed nonlinear field [H].
    pub fn active_inductance_h(&self) -> &[f64] {
        &self.active_inductance_h
    }
    /// Full differential block including explicit end self and mutual flux [H].
    pub fn inductance_h(&self) -> &[f64] {
        &self.inductance_h
    }
    /// Independent diagonal end contribution, excluding shared-mode self terms.
    pub fn external_inductance_h(&self) -> &[f64] {
        self.end_inductance.diagonal_h()
    }
    /// Fixed passive end storage in the original terminal ordering.
    pub fn end_inductance_model(&self) -> &EndWindingInductance {
        &self.end_inductance
    }
    /// Full symmetric end self/mutual matrix [H], before adding meshed flux.
    pub fn end_inductance_matrix_h(&self) -> &[f64] {
        self.end_inductance.matrix_h()
    }
    pub fn maximum_pairing_defect_h(&self) -> f64 {
        self.maximum_pairing_defect_h
    }

    /// Original Whitney pairing for a complete A; excludes external series flux.
    /// This observes a supplied cochain and does not admit it as a physical field.
    pub fn pair_potential_wb(
        &self,
        potential_wb: &[f64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Vec<f64>> {
        poll(cancelled)?;
        if potential_wb.len() != self.operating_field.problem().mesh().edges.len()
            || potential_wb.iter().any(|v| !v.is_finite())
        {
            return Err(invalid("complete finite original edge potential required").into());
        }
        self.loads_per_a
            .iter()
            .map(|load| pairing(load, potential_wb, cancelled))
            .collect()
    }

    /// Local tangent delta_lambda = L_diff delta_i, without the fixed bias flux.
    pub fn linearized_flux_increment_wb(&self, current_increment_a: &[f64]) -> Result<Vec<f64>> {
        let n = self.winding_count();
        if current_increment_a.len() != n || current_increment_a.iter().any(|v| !v.is_finite()) {
            return Err(
                invalid("one finite current increment per original winding required").into(),
            );
        }
        let end_flux = self.end_inductance.flux_linkages_wb(current_increment_a)?;
        self.active_inductance_h
            .chunks_exact(n)
            .zip(end_flux)
            .map(|(row, end)| {
                let mut value = 0.0;
                for (&l, &i) in row.iter().zip(current_increment_a) {
                    value = finite(l.mul_add(i, value))?;
                }
                Ok(finite(value + end)?)
            })
            .collect()
    }

    /// Inductive voltage [V] for a local current rate [A/s] at this held bias.
    /// A fixed remanent or boundary flux never creates an extra voltage term.
    pub fn voltage_from_current_rate_v(&self, current_rate_a_s: &[f64]) -> Result<Vec<f64>> {
        self.linearized_flux_increment_wb(current_rate_a_s)
    }
}

fn pairing(a: &[f64], b: &[f64], cancelled: &mut impl FnMut() -> bool) -> Result<f64> {
    let mut value = 0.0;
    for (&a, &b) in a.iter().zip(b) {
        poll(cancelled)?;
        value = finite(a.mul_add(b, value))?;
    }
    Ok(value)
}

fn admit_smooth_operating_point(
    field: &SpatialNonlinearVectorField<'_>,
    margin_t: f64,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<()> {
    for (cell, material) in field.problem().materials().iter().enumerate() {
        poll(cancelled)?;
        match material {
            SpatialNonlinearMaterial::Reversible(law) => {
                let magnitude = finite(norm(field.flux_density_t()[cell]))?;
                if law.maximum_flux_density_t() - magnitude <= margin_t {
                    return Err(SpatialIncrementalWindingError::SourceBoundary { cell });
                }
                // A radial law is smooth at B=0 because its first segment is
                // linear through the origin. Only interior slope changes matter.
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
                        return Err(SpatialIncrementalWindingError::Nondifferentiable {
                            cell,
                            observed_t: magnitude,
                            breakpoint_t: knots[knot].0,
                        });
                    }
                }
            }
            SpatialNonlinearMaterial::Recoil(law) => {
                let h = field.field_strength_a_m()[cell];
                let direction = law.direction();
                let parallel = finite(dot(h, direction))?;
                let transverse = finite(norm(std::array::from_fn(|k| {
                    h[k] - parallel * direction[k]
                })))?;
                let source = law.material();
                let range = source.parallel_range_a_m();
                if parallel <= range[0]
                    || parallel >= range[1]
                    || transverse >= source.maximum_transverse_a_m()
                {
                    return Err(SpatialIncrementalWindingError::SourceBoundary { cell });
                }
            }
            SpatialNonlinearMaterial::Linear(_) => {}
        }
    }
    Ok(())
}
