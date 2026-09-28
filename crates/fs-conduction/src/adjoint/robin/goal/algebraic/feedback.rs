//! Affine Robin-reference feedback on the actual retained linear FEM operator.
//! The references vary with wall means; they are not frozen at an iterate.

mod solve;

use std::collections::BTreeMap;
use fs_sparse::Csr;
use fs_solver::goal::feedback::{
    FeedbackResidualLimits, FeedbackResidualReport, enclose_affine_feedback_error_with_inverse,
    enclose_affine_feedback_error_with_schur as enclose_affine_feedback_error,
    enclose_affine_feedback_error_with_spectral,
    enclose_affine_feedback_error_with_sparse_inverse,
    enclose_affine_feedback_error_with_spectral_inverse,
};
use fs_solver::goal::inverse::spectral::{
    SpectralError, SpectralInverseLimits, SpectralPreparation, prepare_spectral_inverse,
};
use super::{ConductionError, Cx, LinearGoalAnalyzer, bounded_solve, invalid, map_enclosure, poll};
use super::super::super::{RobinPort, bind_ports};

/// Additional explicit resources for lowering and checking a feedback law.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct RobinFeedbackAnalysisConfig {
    /// Full coupled-residual limits, including the repeated response checks.
    pub residual: FeedbackResidualLimits,
    /// Total inner CG work shared by ALL response columns. Zero selects the
    /// cheaper norm-only gain bound. Does not expand the original solid budget.
    pub max_response_iterations: usize,
    /// Conservative face/port/node work estimate checked before lowering.
    pub max_lowering_entries: usize,
}

/// A cached linear solid plus an explicit affine law `r = offset + D wall`.
/// The stored coefficient construction is ordinary binary64 assembly; the
/// subsequent residual/inverse bounds concern those exact stored coefficients.
/// No assembly-rounding, nonlinear, radiation or physical-model bound follows.
pub struct LinearRobinFeedbackAnalyzer<'m> {
    solid: LinearGoalAnalyzer<'m>,
    ports: Vec<RobinPort>,
    injection: Csr,
    feedback: Csr,
    offset: Vec<f64>,
    responses: Option<Vec<Vec<f64>>>,
    response_iterations: usize,
    config: RobinFeedbackAnalysisConfig,
    spectral_preparation: Option<SpectralPreparation>,
}

/// Regional maximum bound for the stored *coupled* solid/reference system.
#[derive(Debug, Clone, PartialEq)]
pub struct LinearRobinMaximumAnalysis {
    nominal_k: f64,
    interval_k: Option<[f64; 2]>,
    half_width_k: Option<f64>,
    free_vertices: usize,
    response_iterations: usize,
    coupled: FeedbackResidualReport,
}
impl LinearRobinMaximumAnalysis {
    /// Maximum of the supplied field in the selected region.
    #[must_use]
    pub const fn nominal_k(&self) -> f64 { self.nominal_k }
    /// Enclosure of the exact stored coupled system's regional maximum.
    #[must_use]
    pub const fn interval_k(&self) -> Option<[f64; 2]> { self.interval_k }
    /// Absolute error allowance including residual and response-solve errors.
    #[must_use]
    pub const fn algebraic_half_width_k(&self) -> Option<f64> { self.half_width_k }
    /// Number of selected free vertices; prescribed vertices stay exact.
    #[must_use]
    pub const fn free_vertices(&self) -> usize { self.free_vertices }
    /// Total response preparation work; assessments never repeat those solves.
    #[must_use]
    pub const fn response_iterations(&self) -> usize { self.response_iterations }
    /// Full coupled residual, contraction and inverse evidence.
    #[must_use]
    pub const fn coupled(&self) -> &FeedbackResidualReport { &self.coupled }
    /// A missing maximum bound never passes a goal tolerance.
    #[must_use]
    pub fn meets_absolute_tolerance(&self, tolerance: f64) -> bool {
        tolerance.is_finite() && tolerance > 0.0
            && self.half_width_k.is_some_and(|error| error <= tolerance)
    }
}

fn checked(value: f64) -> Result<f64, ConductionError> {
    if value.is_finite() { Ok(value) } else { Err(invalid("nonfinite affine Robin coefficient assembly")) }
}
fn zeros(cx: &Cx<'_>, n: usize) -> Result<Vec<f64>, ConductionError> {
    let mut v = Vec::new();
    v.try_reserve_exact(n).map_err(|_| invalid("affine Robin allocation refused"))?;
    for i in 0..n { if i % 512 == 0 { poll(cx, i)?; } v.push(0.0); }
    Ok(v)
}
fn csr(
    cx: &Cx<'_>, columns: usize, rows: &[BTreeMap<usize, f64>],
) -> Result<Csr, ConductionError> {
    let mut ptr = Vec::with_capacity(rows.len() + 1);
    let mut indices = Vec::new();
    let mut values = Vec::new();
    ptr.push(0);
    for (i, row) in rows.iter().enumerate() {
        poll(cx, i)?;
        for (&column, &value) in row {
            poll(cx, indices.len())?;
            indices.push(column); values.push(checked(value)?);
        }
        ptr.push(indices.len());
    }
    Csr::try_from_parts_with_checkpoint(rows.len(), columns, ptr, indices, values, || poll(cx, 0))?
        .ok_or_else(|| invalid("noncanonical affine Robin transfer"))
}
fn accumulate(map: &mut BTreeMap<usize, f64>, key: usize, value: f64) -> Result<(), ConductionError> {
    if value != 0.0 {
        let entry = map.entry(key).or_insert(0.0);
        *entry = checked(*entry + value)?;
    }
    Ok(())
}

impl<'m> LinearGoalAnalyzer<'m> {
    /// Bind `r_j = reference_offset_k[j] + sum_k D[j,k] wall_mean_k(T)`.
    /// D is row-major p x p; regions define the EXACT coordinate order.
    /// Other boundaries and fixed matching contacts stay in the retained solid.
    /// Prescribed-node wall contributions are included in the affine offset,
    /// never dropped by reducing the unknowns.
    ///
    /// Response columns solve the existing A x = B[:,j] using a single shared
    /// iteration cap. Finite incomplete columns remain proposals: the outward
    /// checker includes their full residual error before admitting a gain.
    /// With zero response budget, the norm-only small-gain route stays usable.
    ///
    /// # Errors
    /// Refuses malformed/duplicate/unknown/nonuniform ports, wrong matrix or
    /// vector shapes, nonfinite assembly, exceeded work/storage limits and
    /// cancellation. A failed sufficient small-gain check is not an input error.
    pub fn with_robin_feedback(
        self, cx: &Cx<'_>, regions: &[&str], reference_offset_k: &[f64],
        reference_wall_matrix: &[f64], config: RobinFeedbackAnalysisConfig,
    ) -> Result<LinearRobinFeedbackAnalyzer<'m>, ConductionError> {
        poll(cx, 0)?;
        let p = regions.len();
        let n = self.dofs().n();
        let full_n = self.problem.mesh.vertex_count();
        let square = p.checked_mul(p).ok_or_else(|| invalid("feedback port geometry overflow"))?;
        if p == 0 || p > config.residual.max_ports
            || reference_offset_k.len() != p || reference_wall_matrix.len() != square
        { return Err(invalid("affine Robin feedback requires p offsets and a bounded p by p matrix")); }
        for (i, &v) in reference_offset_k.iter().chain(reference_wall_matrix).enumerate() {
            if i % 512 == 0 { poll(cx, i)?; } checked(v)?;
        }
        // Bound port binding, lowering dense port mixing, and sparse response
        // storage before any maps, port-face copies or repeated solves.
        let lower_work = self.problem.mesh.boundary().len().checked_mul(3)
            .and_then(|v| v.checked_mul(p + 1))
            .and_then(|v| full_n.checked_mul(p).and_then(|w| v.checked_add(w)))
            .and_then(|v| v.checked_add(square))
            .ok_or_else(|| invalid("affine Robin lowering work overflow"))?;
        if lower_work > config.max_lowering_entries
            || n > config.residual.solid.max_rows
            || self.response.matrix.nnz() > config.residual.solid.max_nonzeros
        { return Err(invalid("affine Robin lowering/solid work limit exceeded")); }
        let response_entries = n.checked_mul(p).ok_or_else(|| invalid("response storage overflow"))?;
        if config.max_response_iterations > 0 && response_entries > config.residual.max_response_entries {
            return Err(invalid("affine Robin response storage limit exceeded"));
        }
        let passes = if config.max_response_iterations > 0 { p + 1 } else { 1 };
        let verify_work = n.checked_add(self.response.matrix.nnz()).and_then(|v| v.checked_mul(passes))
            .ok_or_else(|| invalid("affine Robin verification work overflow"))?;
        // A prepared maximum may already own a checked-inverse proposal
        // because consistent contact prevents comparison dominance. Retain
        // it, and admit its verification before any new response work.
        let verify_work = if self.inverse_columns.is_some() {
            let dense = n.checked_mul(n).ok_or_else(|| invalid("solid inverse storage overflow"))?;
            let inverse_work = n.checked_add(self.response.matrix.nnz())
                .and_then(|pass| n.checked_mul(pass))
                .ok_or_else(|| invalid("solid inverse verification overflow"))?;
            if n > 256 || dense > config.residual.solid.max_nonzeros
                || inverse_work > config.residual.solid.max_nonzeros
            { return Err(invalid("affine Robin solid inverse limit exceeded")); }
            verify_work.checked_add(dense).and_then(|work| work.checked_add(inverse_work))
                .ok_or_else(|| invalid("affine Robin combined verification overflow"))?
        } else { verify_work };
        if verify_work > config.residual.max_verification_entries {
            return Err(invalid("affine Robin verification work limit exceeded"));
        }
        let ports = bind_ports(cx, self.problem, regions)?;
        let mut means = Vec::with_capacity(p);
        let mut b_rows: Vec<BTreeMap<usize, f64>> = Vec::new();
        b_rows.try_reserve_exact(n).map_err(|_| invalid("Robin injection allocation refused"))?;
        for i in 0..n { if i % 512 == 0 { poll(cx, i)?; } b_rows.push(BTreeMap::new()); }
        let mut entries = 0_usize;
        for (j, port) in ports.iter().enumerate() {
            let mut mean = BTreeMap::new();
            for (vertices, area) in &port.faces {
                poll(cx, j)?;
                let weight = checked((area / 3.0) / port.area_m2)?;
                let load = checked(port.htc_w_m2_k * (area / 3.0))?;
                for &vertex in vertices {
                    accumulate(&mut mean, vertex, weight)?;
                    if let Some(i) = self.dofs().slot_of(vertex) {
                        if load != 0.0 && !b_rows[i].contains_key(&j) {
                            entries = entries.checked_add(1).ok_or_else(|| invalid("Robin injection storage overflow"))?;
                            if entries > config.residual.max_transfer_nonzeros {
                                return Err(invalid("affine Robin injection storage limit exceeded"));
                            }
                        }
                        accumulate(&mut b_rows[i], j, load)?;
                    }
                }
            }
            means.push(mean);
        }
        let mut c_rows = Vec::with_capacity(p);
        let mut offsets = zeros(cx, p)?;
        for j in 0..p {
            poll(cx, j)?;
            offsets[j] = checked(reference_offset_k[j] - ports[j].reference_k)?;
            let mut row = BTreeMap::new();
            for (k, mean) in means.iter().enumerate() {
                let coefficient = reference_wall_matrix[j * p + k];
                if coefficient == 0.0 { continue; }
                for (&vertex, &weight) in mean {
                    poll(cx, vertex)?;
                    let value = checked(coefficient * weight)?;
                    if let Some(i) = self.dofs().slot_of(vertex) {
                        if value != 0.0 && !row.contains_key(&i) {
                            entries = entries.checked_add(1).ok_or_else(|| invalid("feedback transfer storage overflow"))?;
                            if entries > config.residual.max_transfer_nonzeros {
                                return Err(invalid("affine Robin feedback storage limit exceeded"));
                            }
                        }
                        accumulate(&mut row, i, value)?;
                    } else {
                        offsets[j] = checked(value.mul_add(self.dofs().prescribed()[vertex], offsets[j]))?;
                    }
                }
            }
            c_rows.push(row);
        }
        let injection = csr(cx, p, &b_rows)?;
        let feedback = csr(cx, n, &c_rows)?;
        let mut response_iterations = 0;
        let responses = if config.max_response_iterations == 0 { None } else {
            let mut responses = Vec::with_capacity(p);
            for j in 0..p {
                poll(cx, response_iterations)?;
                let mut rhs = zeros(cx, n)?;
                for (i, slot) in rhs.iter_mut().enumerate() {
                    if i % 512 == 0 { poll(cx, i)?; }
                    *slot = injection.get(i, j);
                }
                let remaining = config.max_response_iterations - response_iterations;
                let linear = super::LinearConfig {
                    max_iterations: remaining.min(self.response.linear.max_iterations),
                    ..self.response.linear
                };
                let (column, _, iterations) = bounded_solve(cx, &self.response.matrix, &rhs, linear, false)?;
                response_iterations += iterations;
                responses.push(column);
            }
            Some(responses)
        };
        poll(cx, response_iterations)?;
        Ok(LinearRobinFeedbackAnalyzer {
            solid: self, ports, injection, feedback, offset: offsets, responses,
            response_iterations, config, spectral_preparation: None,
        })
    }
}

impl LinearRobinFeedbackAnalyzer<'_> {
    /// Prepare a sparse solid-inverse proof once for repeated field assessments.
    ///
    /// The factor is built and checked against this analyzer's immutable stored
    /// operator, not a caller-supplied matrix or stability constant. Preparation
    /// has its own explicit work/storage allowance; its system limits cannot
    /// expand the original solid envelope. Later assessments spend only their
    /// existing residual/response-verification allowance, never factor again.
    /// Inspect `spectral_preparation()` for the actual paid work and stop reason.
    /// A failed sufficient proof is retained without a certificate; ordinary
    /// assessment then keeps its prior comparison/dense-inverse behavior.
    /// Neither success nor failure changes the physical field or re-solves a
    /// response column. This opt-in builder may be called once per analyzer.
    ///
    /// # Errors
    /// Invalid proof inputs, original system limits, allocation and cancellation
    /// refuse. Sparse work/storage exhaustion returns an analyzer with an
    /// explicit unsuccessful preparation, not fabricated inverse authority.
    pub fn with_spectral_inverse(
        mut self, cx: &Cx<'_>, initial_shift: f64, mut limits: SpectralInverseLimits,
    ) -> Result<Self, ConductionError> {
        poll(cx, 0)?;
        if self.spectral_preparation.is_some() {
            return Err(invalid("spectral preparation is already retained by this analyzer"));
        }
        limits.system.max_rows = limits.system.max_rows.min(self.config.residual.solid.max_rows);
        limits.system.max_nonzeros = limits.system.max_nonzeros.min(self.config.residual.solid.max_nonzeros);
        let preparation = prepare_spectral_inverse(
            &self.solid.response.matrix, initial_shift, limits, || cx.checkpoint().is_ok(),
        ).map_err(|error| match error {
            SpectralError::Residual(error) => map_enclosure(error),
            error => invalid(&error.to_string()),
        })?;
        poll(cx, preparation.work_entries)?;
        self.spectral_preparation = Some(preparation);
        Ok(self)
    }

    /// One-time proof preparation, including failed attempts and consumed work.
    /// The certificate, when present, owns its exact matrix and is read-only.
    /// Per-assessment diagnostics separately report zero *new* preparation work.
    #[must_use]
    pub fn spectral_preparation(&self) -> Option<&SpectralPreparation> {
        self.spectral_preparation.as_ref()
    }

    /// Bound production Robin port geometry and coefficients, in declared order.
    #[must_use]
    pub fn ports(&self) -> &[RobinPort] { &self.ports }
    /// Work spent preparing all response columns together.
    #[must_use]
    pub const fn response_iterations(&self) -> usize { self.response_iterations }
    /// Original solid stability work, including any inverse proposals.
    /// Assessments never repeat the solves that produced these columns.
    #[must_use]
    pub const fn stability_iterations(&self) -> usize { self.solid.stability_iterations }
    /// Borrowed solid inverse proposals, not unchecked inverse authority.
    /// The complete coupled evaluator rechecks these against its own matrix.
    #[must_use]
    pub fn solid_inverse_columns(&self) -> Option<&[Vec<f64>]> { self.solid.inverse_columns() }
    /// Read-only stored transfer for independent numerical replay: A x = b+B(d+C x).
    #[must_use]
    pub fn stored_system(&self) -> (&Csr, &[f64], &Csr, &Csr, &[f64]) {
        (&self.solid.response.matrix, &self.solid.rhs, &self.injection, &self.feedback, &self.offset)
    }
    /// Assess the full coupled equation at another field without another solve.
    /// Region selection, material support and prescribed values use the existing
    /// maximum analyzer's admission. Active-node relocation is included.
    /// If whole-state contraction is inconclusive, checked response errors may
    /// establish port-Schur dominance within the same verification budget.
    /// Neither condition is a claim about fixed-point or dynamical stability.
    ///
    /// # Errors
    /// Existing field/selection/refusal rules, numerical range and cancellation.
    pub fn analyze_maximum(
        &self, cx: &Cx<'_>, temperature: &[f64], vertices: &[usize],
    ) -> Result<LinearRobinMaximumAnalysis, ConductionError> {
        self.analyze_maximum_inner(cx, temperature, vertices, None)
    }

    /// Assess the same complete feedback system with a bounded sparse checked
    /// Gram-inverse fallback when comparison dominance and a retained dense
    /// inverse are unavailable. No new primal/response solve is performed.
    ///
    /// The sparse limits bound fill and work, not measured peak RSS. Every
    /// shift is independently checked; the complete feedback inverse must
    /// still be established. Existing maximum assessment is unchanged unless
    /// this explicit route is requested.
    ///
    /// # Errors
    /// The usual field/selection errors plus sparse proof admission, resource,
    /// arithmetic and cancellation refusals. Missing evidence is never zero.
    pub fn analyze_maximum_with_spectral_fallback(
        &self, cx: &Cx<'_>, temperature: &[f64], vertices: &[usize],
        limits: SpectralInverseLimits,
    ) -> Result<LinearRobinMaximumAnalysis, ConductionError> {
        self.analyze_maximum_inner(cx, temperature, vertices, Some(limits))
    }

    fn analyze_maximum_inner(
        &self, cx: &Cx<'_>, temperature: &[f64], vertices: &[usize],
        spectral: Option<SpectralInverseLimits>,
    ) -> Result<LinearRobinMaximumAnalysis, ConductionError> {
        // Validate the region and field, but do not first certify a frozen
        // solid. That would repeat dense inverse verification and discard its
        // result. The complete coupled equation owns the ONE numerical pass.
        let (fixed, moving, free_vertices) = selected_maximum(&self.solid, cx, temperature, vertices)?;
        let free = self.solid.dofs().gather(temperature);
        let coupled = if let Some(certificate) = self.spectral_preparation.as_ref()
            .and_then(|preparation| preparation.certificate.as_ref())
        {
            enclose_affine_feedback_error_with_spectral_inverse(
                certificate, &self.solid.rhs, &free,
                &self.injection, &self.feedback, &self.offset,
                self.responses.as_deref(), self.solid.stability_scaling.as_deref(),
                self.config.residual, || cx.checkpoint().is_ok(),
            )
        } else { match (self.solid.inverse_columns(), self.solid.sparse_inverse()) {
            (Some(columns), _) => enclose_affine_feedback_error_with_inverse(
                &self.solid.response.matrix, &self.solid.rhs, &free,
                &self.injection, &self.feedback, &self.offset,
                self.responses.as_deref(), self.solid.stability_scaling.as_deref(),
                columns, self.config.residual, || cx.checkpoint().is_ok(),
            ),
            // Interval LDL evidence is cheaper than the shifted-Gram fallback
            // on large obtuse meshes; it runs before any spectral preparation.
            (None, Some(evidence)) => enclose_affine_feedback_error_with_sparse_inverse(
                &self.solid.response.matrix, &self.solid.rhs, &free,
                &self.injection, &self.feedback, &self.offset,
                self.responses.as_deref(), self.solid.stability_scaling.as_deref(),
                evidence, self.config.residual, || cx.checkpoint().is_ok(),
            ),
            (None, None) => match spectral {
                Some(limits) => match enclose_affine_feedback_error_with_spectral(
                    &self.solid.response.matrix, &self.solid.rhs, &free,
                    &self.injection, &self.feedback, &self.offset,
                    self.responses.as_deref(), self.solid.stability_scaling.as_deref(),
                    self.config.residual, limits, || cx.checkpoint().is_ok(),
                ) {
                    Ok(report) => Ok(report),
                    Err(SpectralError::Residual(error)) => Err(error),
                    Err(error) => return Err(invalid(&error.to_string())),
                },
                None => enclose_affine_feedback_error(
                    &self.solid.response.matrix, &self.solid.rhs, &free,
                    &self.injection, &self.feedback, &self.offset,
                    self.responses.as_deref(), self.solid.stability_scaling.as_deref(),
                    self.config.residual, || cx.checkpoint().is_ok(),
                ),
            },
        } }.map_err(map_enclosure)?;
        let error = if free_vertices == 0 { Some(0.0) } else { coupled.state_error_infinity_upper() };
        let nominal = fixed.max(moving);
        let interval = error.and_then(|error| {
            if free_vertices == 0 || error == 0.0 { return Some([nominal, nominal]); }
            let lo = fixed.max(fs_math::next_down(moving - error));
            let hi = fixed.max(fs_math::next_up(moving + error));
            (lo.is_finite() && hi.is_finite()).then_some([lo, hi])
        });
        poll(cx, 0)?;
        Ok(LinearRobinMaximumAnalysis {
            nominal_k: nominal, interval_k: interval, half_width_k: interval.and(error),
            free_vertices, response_iterations: self.response_iterations, coupled,
        })
    }
}

// Selection-only admission: no residual, inverse, primal or response solve.
// Canonical mesh order and normalized zeros preserve set-valued max semantics.
fn selected_maximum(
    solid: &LinearGoalAnalyzer<'_>, cx: &Cx<'_>, temperature: &[f64], vertices: &[usize],
) -> Result<(f64, f64, usize), ConductionError> {
    poll(cx, 0)?;
    let n = solid.problem.mesh.vertex_count();
    if vertices.is_empty() || vertices.len() > n {
        return Err(invalid("a regional maximum needs a nonempty vertex set within the mesh"));
    }
    let mut selected = Vec::new();
    selected.try_reserve_exact(n).map_err(|_| invalid("regional maximum selection allocation refused"))?;
    for i in 0..n { if i % 512 == 0 { poll(cx, i)?; } selected.push(false); }
    for (i, &vertex) in vertices.iter().enumerate() {
        if i % 512 == 0 { poll(cx, i)?; }
        if vertex >= n || selected[vertex] {
            return Err(invalid("regional maximum vertices must be distinct and in range"));
        }
        selected[vertex] = true;
    }
    super::validate_field(cx, solid.problem, solid.dofs(), temperature)?;
    let (mut fixed, mut moving) = (f64::NEG_INFINITY, f64::NEG_INFINITY);
    let mut free_vertices = 0;
    for (vertex, included) in selected.into_iter().enumerate() {
        if vertex % 512 == 0 { poll(cx, vertex)?; }
        if !included { continue; }
        let value = if temperature[vertex] == 0.0 { 0.0 } else { temperature[vertex] };
        if solid.dofs().slot_of(vertex).is_some() {
            moving = moving.max(value); free_vertices += 1;
        } else { fixed = fixed.max(value); }
    }
    poll(cx, 0)?;
    Ok((fixed, moving, free_vertices))
}
