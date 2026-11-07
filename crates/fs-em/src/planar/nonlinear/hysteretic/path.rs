//! Adaptive continuation of an explicitly authored affine magnetic load path.
//!
//! Every comparison uses the original frozen-history field solve: one full
//! increment and two sequential halves. Only the two-half history is retained.
//! Polarization, spatial pinning heat and supplied work are resolved separately;
//! the endpoint numerical remainder is never physical heat. These are empirical
//! refinement indicators, not certified global bounds or a reversal detector.
//! Callers supply separate segments at waveform turns and source events.
//!
//! This is rate-independent magnetostatics. It introduces no clock, eddy-current
//! evolution, thermal step, changing material law or extrapolated internal state.

use super::{
    MagneticError, NonlinearMagneticError, PlanarHystereticCandidate, PlanarHystereticControls,
    PlanarHystereticMagnetostatics, PlanarHystereticState, PlanarHystereticStep,
    PlanarMagnetostatics,
};
use fs_material::field_storage::hysteresis::vector::VectorHysteresisError;

/// Independent accuracy indicators and finite work across the complete segment.
#[derive(Clone, Copy, Debug)]
pub struct MagneticPathControls {
    /// Original field and return-map controls; no tolerance is relaxed.
    pub field: PlanarHystereticControls,
    pub relative_tolerance: f64,
    pub absolute_polarization_tolerance_t: f64,
    pub absolute_heat_tolerance_j: f64,
    pub absolute_work_tolerance_j: f64,
    /// Acceptance gate for cumulative physical accounting, not an error estimate.
    pub energy_tolerance_j: f64,
    /// Includes coarse, half and rejected comparisons. Each solve additionally
    /// retains its original Newton, Krylov, backtracking and constitutive caps.
    pub maximum_field_solves: usize,
    /// Total Newton updates across all attempted field solves. Zero admits only
    /// initially stationary solves, including fully prescribed fields.
    pub maximum_newton_iterations: usize,
    /// At most 52, so every attempted dyadic fraction remains distinct.
    pub maximum_refinement_depth: usize,
}

#[derive(Debug)]
pub enum MagneticPathError {
    Invalid(&'static str),
    Field(NonlinearMagneticError),
    Cancelled,
    FieldSolveBudget,
    NewtonBudget,
    /// The authored starting current does not support the accepted field on
    /// free coordinates. No initialization impulse or history update is made.
    InitialEquilibrium {
        residual_a_m: f64,
        tolerance_a_m: f64,
    },
    RefinementLimit {
        fraction: [f64; 2],
        depth: usize,
        maximum_defect_ratio: f64,
    },
    /// Native range/globalization refusal survived every permitted subdivision.
    FieldRefinementLimit {
        fraction: [f64; 2],
        depth: usize,
        error: NonlinearMagneticError,
    },
}
impl From<NonlinearMagneticError> for MagneticPathError {
    fn from(error: NonlinearMagneticError) -> Self {
        match error {
            NonlinearMagneticError::Magnetic(MagneticError::Cancelled) => Self::Cancelled,
            other => Self::Field(other),
        }
    }
}
impl From<MagneticError> for MagneticPathError {
    fn from(error: MagneticError) -> Self {
        NonlinearMagneticError::from(error).into()
    }
}
impl core::fmt::Display for MagneticPathError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "magnetic loading path: {self:?}")
    }
}
impl std::error::Error for MagneticPathError {}
type Result<T> = core::result::Result<T, MagneticPathError>;

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MagneticPathSummary {
    /// All attempted native solves, including failed and discarded trials.
    pub field_solves: usize,
    pub failed_field_solves: usize,
    pub newton_iterations: usize,
    pub linear_iterations: usize,
    pub backtracks: usize,
    pub accepted_substeps: usize,
    pub rejected_intervals: usize,
    pub maximum_refinement_depth: usize,
    /// Sum of maximum branch defects over accepted comparisons.
    pub estimated_polarization_error_t: f64,
    /// L1 per-cell physical heat defects, so spatial errors cannot cancel.
    pub estimated_pinning_heat_error_j: f64,
    /// Sum of max(work defect, retained endpoint numerical remainder).
    pub estimated_supplied_work_error_j: f64,
    pub maximum_accepted_defect_ratio: f64,
}

pub struct MagneticPathReport {
    /// Cumulative accepted physical transfers. The constitutive flux residuals
    /// and `nonlinear` report describe the final native solve; work counts for
    /// the entire comparison history live in `summary`.
    pub account: PlanarHystereticStep,
    pub summary: MagneticPathSummary,
}

pub struct MagneticPathCandidate<'a> {
    state: PlanarHystereticState<'a>,
    report: MagneticPathReport,
}
impl<'a> MagneticPathCandidate<'a> {
    pub fn state(&self) -> &PlanarHystereticState<'a> {
        &self.state
    }
    pub fn report(&self) -> &MagneticPathReport {
        &self.report
    }
    pub fn into_parts(self) -> (PlanarHystereticState<'a>, MagneticPathReport) {
        (self.state, self.report)
    }
    /// Move cumulative native heat/work into an existing consumer without a
    /// field solve, receipt clone, second deposition or history reset.
    pub fn into_magnetic_candidate(self) -> (PlanarHystereticCandidate<'a>, MagneticPathSummary) {
        (
            PlanarHystereticCandidate {
                state: self.state,
                step: self.report.account,
            },
            self.report.summary,
        )
    }
}

fn invalid(reason: &'static str) -> MagneticPathError {
    MagneticPathError::Invalid(reason)
}
fn finite(value: f64) -> Result<f64> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(invalid(
            "magnetic loading path arithmetic is unrepresentable",
        ))
    }
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<()> {
    if cancelled() {
        Err(MagneticPathError::Cancelled)
    } else {
        Ok(())
    }
}
fn add(sum: &mut f64, correction: &mut f64, value: f64) -> Result<()> {
    let next = finite(*sum + value)?;
    *correction = finite(
        *correction
            + if sum.abs() >= value.abs() {
                (*sum - next) + value
            } else {
                (value - next) + *sum
            },
    )?;
    *sum = next;
    Ok(())
}
fn count(total: &mut usize, extra: usize) -> Result<()> {
    *total = total
        .checked_add(extra)
        .ok_or_else(|| invalid("path work count overflowed"))?;
    Ok(())
}
fn distance(a: [f64; 3], b: [f64; 3]) -> Result<f64> {
    finite((a[0] - b[0]).hypot(a[1] - b[1]).hypot(a[2] - b[2]))
}
fn ratio(error: f64, tolerance: f64) -> f64 {
    if error == 0.0 { 0.0 } else { error / tolerance }
}
fn interpolate(a: f64, b: f64, fraction: f64) -> Result<f64> {
    // Convex weights avoid an overflowing b-a. Dyadic fractions are exact.
    if fraction == 1.0 {
        Ok(b)
    } else {
        finite((1.0 - fraction).mul_add(a, fraction * b))
    }
}

struct Loading<'s> {
    initial_current: &'s [f64],
    target_current: &'s [f64],
    boundary: Vec<(usize, f64, f64, f64, f64)>,
}
impl<'s> Loading<'s> {
    #[allow(clippy::too_many_arguments)]
    fn new(
        problem: &PlanarHystereticMagnetostatics,
        previous: &PlanarHystereticState<'_>,
        initial_current: &'s [f64],
        target_current: &'s [f64],
        boundary: &[(usize, f64)],
        controls: PlanarHystereticControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self> {
        if !std::ptr::eq(problem, previous.problem())
            || initial_current.len() != problem.materials().len()
            || target_current.len() != initial_current.len()
        {
            return Err(invalid(
                "original history and complete starting/target currents required",
            ));
        }
        let prescribed =
            PlanarMagnetostatics::admit_mesh_boundary(problem.mesh(), boundary, cancelled)?;
        let references =
            PlanarMagnetostatics::mesh_gauge_references(problem.mesh(), &prescribed, cancelled)?;
        let mut load = vec![0.0; problem.mesh().vertices().len()];
        for (cell, face) in problem.mesh().faces().iter().enumerate() {
            poll(cancelled)?;
            finite(target_current[cell])?;
            let value = finite(initial_current[cell] * (problem.volumes_m3[cell] / 3.0))?;
            for &vertex in face {
                load[vertex as usize] = finite(load[vertex as usize] + value)?;
            }
        }
        let effort = problem.magnetic_load(previous, cancelled)?;
        let (mut residual, mut reference) = (0.0_f64, 0.0_f64);
        let mut traces = Vec::with_capacity(boundary.len());
        for (v, value) in prescribed.iter().enumerate() {
            poll(cancelled)?;
            if let Some(target) = *value {
                let anchor = references[v];
                let gauge = prescribed[anchor].expect("admitted component gauge");
                traces.push((
                    v,
                    finite(previous.potential_wb_m()[v] - previous.potential_wb_m()[anchor])?,
                    finite(target - gauge)?,
                    gauge,
                    target,
                ));
            } else {
                residual = finite(residual.hypot(finite(effort[v] - load[v])?))?;
                reference = finite(reference.hypot(load[v]))?;
            }
        }
        // This checks the retained constitutive H directly. Calling a zero
        // increment solve here could erase an inconsistent initialization.
        let tolerance = controls
            .nonlinear
            .absolute_residual_tolerance_a_m
            .max(finite(
                controls.nonlinear.relative_residual_tolerance * reference,
            )?);
        if residual > tolerance {
            return Err(MagneticPathError::InitialEquilibrium {
                residual_a_m: residual,
                tolerance_a_m: tolerance,
            });
        }
        Ok(Self {
            initial_current,
            target_current,
            boundary: traces,
        })
    }
    fn at(
        &self,
        fraction: f64,
        current: &mut [f64],
        boundary: &mut Vec<(usize, f64)>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<()> {
        for (i, out) in current.iter_mut().enumerate() {
            poll(cancelled)?;
            *out = interpolate(self.initial_current[i], self.target_current[i], fraction)?;
        }
        boundary.clear();
        for &(v, a, b, gauge, target) in &self.boundary {
            poll(cancelled)?;
            // Interpolate physical differences per connected component. Install
            // the authored final gauge in every native solve, so no post-solve
            // shift can make retained curl(A) disagree with the actual field.
            let value = if fraction == 1.0 {
                target
            } else {
                finite(interpolate(a, b, fraction)? + gauge)?
            };
            boundary.push((v, value));
        }
        Ok(())
    }
}

struct Work {
    controls: MagneticPathControls,
    solves: usize,
    failed: usize,
    newton: usize,
    linear: usize,
    backtracks: usize,
}
impl Work {
    fn solve<'a>(
        &mut self,
        previous: &PlanarHystereticState<'a>,
        current: &[f64],
        boundary: &[(usize, f64)],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<PlanarHystereticCandidate<'a>> {
        poll(cancelled)?;
        if self.solves == self.controls.maximum_field_solves {
            return Err(MagneticPathError::FieldSolveBudget);
        }
        let remaining = self.controls.maximum_newton_iterations - self.newton;
        let mut controls = self.controls.field;
        controls.nonlinear.max_newton_iterations =
            controls.nonlinear.max_newton_iterations.min(remaining);
        self.solves += 1;
        let result = previous
            .problem()
            .candidate(previous, current, boundary, controls, cancelled);
        let candidate = match result {
            Ok(candidate) => candidate,
            Err(error) => {
                count(&mut self.failed, 1)?;
                // Range errors escaping the native solve arise at its INITIAL
                // evaluation; trial range failures become a bounded nonlinear
                // report. Count that failed work before any continuation retry.
                if let NonlinearMagneticError::Nonconvergence(ref report) = error {
                    count(&mut self.newton, report.newton_iterations)?;
                    count(&mut self.linear, report.linear_iterations)?;
                    count(&mut self.backtracks, report.backtracks)?;
                    if report.newton_iterations == remaining
                        && remaining <= self.controls.field.nonlinear.max_newton_iterations
                    {
                        return Err(MagneticPathError::NewtonBudget);
                    }
                }
                return Err(error.into());
            }
        };
        let report = &candidate.step().nonlinear;
        count(&mut self.newton, report.newton_iterations)?;
        count(&mut self.linear, report.linear_iterations)?;
        count(&mut self.backtracks, report.backtracks)?;
        Ok(candidate)
    }
}

fn subdivide_failure(
    error: MagneticPathError,
    interval: (f64, f64, usize),
    maximum_depth: usize,
    pending: &mut Vec<(f64, f64, usize)>,
    rejected: &mut usize,
) -> Result<()> {
    let MagneticPathError::Field(error) = error else {
        return Err(error);
    };
    if !matches!(
        &error,
        NonlinearMagneticError::FieldRange { .. }
            | NonlinearMagneticError::History {
                error: VectorHysteresisError::FluxRange { .. },
                ..
            }
            | NonlinearMagneticError::Nonconvergence(_)
    ) {
        // Invalid inputs, cancellation, source failures and constitutive/CG
        // work refusals must never be recast as a physical refinement problem.
        return Err(MagneticPathError::Field(error));
    }
    let (a, b, depth) = interval;
    if depth == maximum_depth {
        return Err(MagneticPathError::FieldRefinementLimit {
            fraction: [a, b],
            depth,
            error,
        });
    }
    count(rejected, 1)?;
    let middle = f64::midpoint(a, b);
    pending.push((middle, b, depth + 1));
    pending.push((a, middle, depth + 1));
    Ok(())
}

fn defects(
    start: &PlanarHystereticState<'_>,
    coarse: &PlanarHystereticCandidate<'_>,
    half: &PlanarHystereticCandidate<'_>,
    fine: &PlanarHystereticCandidate<'_>,
    fraction: f64,
    c: MagneticPathControls,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<([f64; 3], f64)> {
    let (mut polarization, mut maximum_ratio) = (0.0_f64, 0.0_f64);
    for ((old, a), b) in start
        .cell_states()
        .iter()
        .zip(coarse.state().cell_states())
        .zip(fine.state().cell_states())
    {
        for ((old, a), b) in old
            .polarizations_t()
            .iter()
            .zip(a.polarizations_t())
            .zip(b.polarizations_t())
        {
            poll(cancelled)?;
            let error = distance(*a, *b)?;
            let scale = distance(*a, *old)?.max(distance(*b, *old)?);
            let allowance = finite(
                fraction * c.absolute_polarization_tolerance_t + c.relative_tolerance * scale,
            )?;
            polarization = polarization.max(error);
            maximum_ratio = maximum_ratio.max(ratio(error, allowance));
        }
    }
    let heat = finite(half.step().irreversible_heat_j + fine.step().irreversible_heat_j)?;
    let (mut heat_error, mut correction) = (0.0, 0.0);
    for ((a, first), second) in coarse
        .step()
        .cells
        .iter()
        .zip(&half.step().cells)
        .zip(&fine.step().cells)
    {
        poll(cancelled)?;
        let error = finite(
            (a.irreversible_heat_j
                - finite(first.irreversible_heat_j + second.irreversible_heat_j)?)
            .abs(),
        )?;
        add(&mut heat_error, &mut correction, error)?;
    }
    heat_error = finite(heat_error + correction)?;
    let heat_allowance = finite(
        fraction * c.absolute_heat_tolerance_j
            + c.relative_tolerance * heat.max(coarse.step().irreversible_heat_j),
    )?;
    let work = finite(half.step().supplied_work_j + fine.step().supplied_work_j)?;
    let numerical =
        finite(half.step().numerical_dissipation_j + fine.step().numerical_dissipation_j)?;
    // Native D_num is formed stably from increments; it prevents a large bias
    // work from hiding an unresolved endpoint quadrature loss by cancellation.
    let work_error = finite((coarse.step().supplied_work_j - work).abs())?.max(numerical);
    let work_scale = finite(half.step().supplied_work_j.abs() + fine.step().supplied_work_j.abs())?
        .max(coarse.step().supplied_work_j.abs());
    let work_allowance =
        finite(fraction * c.absolute_work_tolerance_j + c.relative_tolerance * work_scale)?;
    maximum_ratio = maximum_ratio
        .max(ratio(heat_error, heat_allowance))
        .max(ratio(work_error, work_allowance));
    Ok(([polarization, heat_error, work_error], maximum_ratio))
}

struct Account {
    cells: Vec<[f64; 6]>,
    corrections: Vec<[f64; 6]>,
    total: [f64; 9],
    correction: [f64; 9],
    last: Option<PlanarHystereticStep>,
}
impl Account {
    fn new(cells: usize) -> Self {
        Self {
            cells: vec![[0.0; 6]; cells],
            corrections: vec![[0.0; 6]; cells],
            total: [0.0; 9],
            correction: [0.0; 9],
            last: None,
        }
    }
    fn push(
        &mut self,
        step: PlanarHystereticStep,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<()> {
        for (i, cell) in step.cells.iter().enumerate() {
            poll(cancelled)?;
            for (k, value) in [
                cell.stored_energy_change_j,
                cell.irreversible_heat_j,
                cell.numerical_dissipation_j,
                cell.field_work_j,
                cell.pinning_heat_j,
                cell.viscous_heat_j,
            ]
            .into_iter()
            .enumerate()
            {
                add(&mut self.cells[i][k], &mut self.corrections[i][k], value)?;
            }
        }
        for (i, value) in [
            step.stored_energy_change_j,
            step.irreversible_heat_j,
            step.numerical_dissipation_j,
            step.impressed_current_work_j,
            step.boundary_work_j,
            step.supplied_work_j,
            step.field_work_j,
            step.pinning_heat_j,
            step.viscous_heat_j,
        ]
        .into_iter()
        .enumerate()
        {
            add(&mut self.total[i], &mut self.correction[i], value)?;
        }
        self.last = Some(step);
        Ok(())
    }
    fn finish(
        mut self,
        tolerance: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<PlanarHystereticStep> {
        let mut result = self
            .last
            .ok_or_else(|| invalid("magnetic path has no accepted substeps"))?;
        let (mut local_defect, mut correction) = (0.0, 0.0);
        for (i, cell) in result.cells.iter_mut().enumerate() {
            poll(cancelled)?;
            for k in 0..6 {
                self.cells[i][k] = finite(self.cells[i][k] + self.corrections[i][k])?;
            }
            let [stored, heat, numerical, work, pinning, viscous] = self.cells[i];
            cell.stored_energy_change_j = stored;
            cell.irreversible_heat_j = heat;
            cell.pinning_heat_j = pinning;
            cell.viscous_heat_j = viscous;
            cell.numerical_dissipation_j = numerical;
            cell.field_work_j = work;
            cell.energy_residual_j = finite(work - stored - heat - numerical)?;
            add(
                &mut local_defect,
                &mut correction,
                cell.energy_residual_j.abs(),
            )?;
        }
        for i in 0..9 {
            self.total[i] = finite(self.total[i] + self.correction[i])?;
        }
        let [
            stored,
            heat,
            numerical,
            impressed,
            boundary,
            supplied,
            work,
            pinning,
            viscous,
        ] = self.total;
        result.stored_energy_change_j = stored;
        result.irreversible_heat_j = heat;
        result.pinning_heat_j = pinning;
        result.viscous_heat_j = viscous;
        result.numerical_dissipation_j = numerical;
        result.impressed_current_work_j = impressed;
        result.boundary_work_j = boundary;
        result.supplied_work_j = supplied;
        result.field_work_j = work;
        result.quadratic_work_j = 0.0;
        result.energy_residual_j = finite(supplied - stored - heat - numerical)?;
        if result.energy_residual_j.abs() > tolerance
            || finite(local_defect + correction)? > tolerance
            || finite(supplied - work)?.abs() > tolerance
            || finite(supplied - impressed - boundary)?.abs() > tolerance
        {
            return Err(invalid(
                "cumulative magnetic storage, heat, numerical loss and source work do not close",
            ));
        }
        Ok(result)
    }
}

impl PlanarHystereticMagnetostatics {
    /// Resolve one affine current-density and exterior-potential segment from
    /// the accepted magnetic history. Starting impressed current is explicit:
    /// field/memory alone does not determine its values at prescribed vertices.
    /// Its native weak equilibrium is checked before any path solve.
    ///
    /// The target boundary set stays fixed over this segment; starting trace
    /// values come from accepted A. Interpolation removes each component's
    /// independent gauge, and every native solve installs the authored target
    /// gauge. Internal states are never shifted or averaged after acceptance.
    /// Absolute allowances scale with interval fraction and relative allowances
    /// with local changes. Waveform turns require separate authored segments.
    #[allow(clippy::too_many_arguments)]
    pub fn affine_path_candidate<'a>(
        &'a self,
        previous: &PlanarHystereticState<'a>,
        initial_current_density_a_m2: &[f64],
        target_current_density_a_m2: &[f64],
        target_boundary_potential_wb_m: &[(usize, f64)],
        controls: MagneticPathControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<MagneticPathCandidate<'a>> {
        poll(cancelled)?;
        let absolute = [
            controls.absolute_polarization_tolerance_t,
            controls.absolute_heat_tolerance_j,
            controls.absolute_work_tolerance_j,
        ];
        let field = controls.field.nonlinear;
        if !controls.relative_tolerance.is_finite()
            || controls.relative_tolerance < 0.0
            || absolute.iter().any(|v| {
                !v.is_finite() || *v < 0.0 || (*v == 0.0 && controls.relative_tolerance == 0.0)
            })
            || !controls.energy_tolerance_j.is_finite()
            || controls.energy_tolerance_j < 0.0
            || controls.maximum_field_solves == 0
            || controls.maximum_refinement_depth > 52
            || !field.absolute_residual_tolerance_a_m.is_finite()
            || field.absolute_residual_tolerance_a_m <= 0.0
            || !field.relative_residual_tolerance.is_finite()
            || !(0.0..1.0).contains(&field.relative_residual_tolerance)
        {
            return Err(invalid(
                "finite path tolerances, positive solve cap, native residual controls and depth at most 52 required",
            ));
        }
        let loading = Loading::new(
            self,
            previous,
            initial_current_density_a_m2,
            target_current_density_a_m2,
            target_boundary_potential_wb_m,
            controls.field,
            cancelled,
        )?;
        let mut currents = initial_current_density_a_m2.to_vec();
        let mut boundary = Vec::with_capacity(target_boundary_potential_wb_m.len());
        let mut current = previous.clone();
        let mut work = Work {
            controls,
            solves: 0,
            failed: 0,
            newton: 0,
            linear: 0,
            backtracks: 0,
        };
        let mut account = Account::new(self.materials().len());
        // Only interval descriptors are stacked. A rejected first half is the
        // left child's coarse solve from the SAME history; right children must
        // restart from the newly refined left history.
        let mut pending = vec![(0.0, 1.0, 0usize)];
        let mut cached = None;
        let (mut accepted, mut rejected, mut deepest) = (0usize, 0usize, 0usize);
        let (mut errors, mut corrections) = ([0.0; 3], [0.0; 3]);
        let mut maximum_ratio = 0.0_f64;
        while let Some((a, b, depth)) = pending.pop() {
            poll(cancelled)?;
            deepest = deepest.max(depth);
            let midpoint = f64::midpoint(a, b);
            loading.at(b, &mut currents, &mut boundary, cancelled)?;
            let coarse = match cached.take() {
                Some(candidate) => candidate,
                None => match work.solve(&current, &currents, &boundary, cancelled) {
                    Ok(candidate) => candidate,
                    Err(error) => {
                        subdivide_failure(
                            error,
                            (a, b, depth),
                            controls.maximum_refinement_depth,
                            &mut pending,
                            &mut rejected,
                        )?;
                        continue;
                    }
                },
            };
            loading.at(midpoint, &mut currents, &mut boundary, cancelled)?;
            let half = match work.solve(&current, &currents, &boundary, cancelled) {
                Ok(candidate) => candidate,
                Err(error) => {
                    subdivide_failure(
                        error,
                        (a, b, depth),
                        controls.maximum_refinement_depth,
                        &mut pending,
                        &mut rejected,
                    )?;
                    continue;
                }
            };
            loading.at(b, &mut currents, &mut boundary, cancelled)?;
            let fine = match work.solve(half.state(), &currents, &boundary, cancelled) {
                Ok(candidate) => candidate,
                Err(error) => {
                    subdivide_failure(
                        error,
                        (a, b, depth),
                        controls.maximum_refinement_depth,
                        &mut pending,
                        &mut rejected,
                    )?;
                    cached = Some(half);
                    continue;
                }
            };
            let (defect, ratio) =
                defects(&current, &coarse, &half, &fine, b - a, controls, cancelled)?;
            if ratio > 1.0 {
                if depth == controls.maximum_refinement_depth {
                    return Err(MagneticPathError::RefinementLimit {
                        fraction: [a, b],
                        depth,
                        maximum_defect_ratio: ratio,
                    });
                }
                count(&mut rejected, 1)?;
                pending.push((midpoint, b, depth + 1));
                pending.push((a, midpoint, depth + 1));
                cached = Some(half);
            } else {
                let (_, first) = half.into_parts();
                let (next, second) = fine.into_parts();
                account.push(first, cancelled)?;
                account.push(second, cancelled)?;
                for i in 0..3 {
                    add(&mut errors[i], &mut corrections[i], defect[i])?;
                }
                maximum_ratio = maximum_ratio.max(ratio);
                count(&mut accepted, 2)?;
                current = next;
            }
        }
        for i in 0..3 {
            errors[i] = finite(errors[i] + corrections[i])?;
        }
        let report = MagneticPathReport {
            account: account.finish(controls.energy_tolerance_j, cancelled)?,
            summary: MagneticPathSummary {
                field_solves: work.solves,
                failed_field_solves: work.failed,
                newton_iterations: work.newton,
                linear_iterations: work.linear,
                backtracks: work.backtracks,
                accepted_substeps: accepted,
                rejected_intervals: rejected,
                maximum_refinement_depth: deepest,
                estimated_polarization_error_t: errors[0],
                estimated_pinning_heat_error_j: errors[1],
                estimated_supplied_work_error_j: errors[2],
                maximum_accepted_defect_ratio: maximum_ratio,
            },
        };
        poll(cancelled)?;
        Ok(MagneticPathCandidate {
            state: current,
            report,
        })
    }
}
