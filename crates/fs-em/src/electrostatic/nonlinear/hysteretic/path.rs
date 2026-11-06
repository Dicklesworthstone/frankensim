//! Adaptive continuation along one explicitly authored affine loading segment.
//!
//! A full native field increment is compared with two sequential half increments.
//! Only the two-half history can be accepted; rejected trials never become memory.
//! Branch polarization, physical switching heat and electrical work have separate
//! tolerances. Work refinement also resolves the native endpoint numerical loss,
//! which is accumulated separately and is never converted into physical heat.
//!
//! The observed full/half defects are estimates, not certified global errors or
//! an event detector. Supply separate segments at waveform turns, control-mode
//! changes and source events, and check tolerance refinement. An affine boundary
//! path can still drive nonlinear interior trajectories. This is rate-independent
//! quasistatic continuation, not a time integrator or a thermal update.

use super::super::{
    ElectrodeControl, NonlinearElectrostaticControls, NonlinearElectrostaticError, add,
};
use super::{
    HystereticElectrostaticCandidate, HystereticElectrostaticStep, PlanarHystereticElectrostatics,
};
use crate::electrostatic::ElectrostaticError;

/// Accuracy estimates and bounded work for one complete affine loading segment.
#[derive(Clone, Copy, Debug)]
pub struct HystereticPathControls {
    /// Original native nonlinear/CG controls for each field solve.
    pub field: NonlinearElectrostaticControls,
    pub relative_tolerance: f64,
    /// Absolute branch-polarization allowance across the whole segment [C/m²].
    pub absolute_polarization_tolerance_c_m2: f64,
    /// Absolute physical switching-heat allowance across the whole segment [J].
    pub absolute_heat_tolerance_j: f64,
    /// Absolute endpoint source-work allowance across the whole segment [J].
    pub absolute_work_tolerance_j: f64,
    /// Physical balance gate on the accumulated accepted path, not an error estimate.
    pub energy_tolerance_j: f64,
    /// Includes coarse, half and rejected trials; never reset on subdivision.
    pub maximum_field_solves: usize,
    /// Across every attempted native solve, in addition to each solve's own cap.
    pub maximum_operator_applications: usize,
    /// Dyadic subdivisions; at most 52 keeps all attempted fractions distinct.
    pub maximum_refinement_depth: usize,
}

#[derive(Debug)]
pub enum HystereticPathError {
    Invalid(&'static str),
    Field(NonlinearElectrostaticError),
    Cancelled,
    FieldSolveBudget,
    OperatorBudget,
    RefinementLimit {
        fraction: [f64; 2],
        depth: usize,
        maximum_defect_ratio: f64,
    },
}
impl From<NonlinearElectrostaticError> for HystereticPathError {
    fn from(error: NonlinearElectrostaticError) -> Self {
        match error {
            NonlinearElectrostaticError::Electric(ElectrostaticError::Cancelled) => Self::Cancelled,
            other => Self::Field(other),
        }
    }
}
impl core::fmt::Display for HystereticPathError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "hysteretic electric loading path: {self:?}")
    }
}
impl std::error::Error for HystereticPathError {}
type Result<T> = core::result::Result<T, HystereticPathError>;

/// Cumulative physical work/heat and observed discretization defects.
pub struct HystereticPathReport {
    /// Sums ONLY accepted substeps, including every cell's physical heat and
    /// numerical loss. Charge residuals and `nonlinear` describe the final solve.
    pub account: HystereticElectrostaticStep,
    pub summary: HystereticPathSummary,
}

/// Work counts and empirical refinement estimates, independent of cell receipts.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct HystereticPathSummary {
    pub field_solves: usize,
    pub operator_applications: usize,
    /// Two physical substeps are retained for every accepted comparison.
    pub accepted_substeps: usize,
    pub rejected_intervals: usize,
    pub maximum_refinement_depth: usize,
    /// Sum of the maximum observed branch-polarization defect per accepted interval.
    pub estimated_polarization_error_c_m2: f64,
    /// Sum of absolute per-cell heat defects; spatial heating errors cannot cancel.
    pub estimated_switching_heat_error_j: f64,
    /// Sum of max(source-work defect, retained endpoint numerical loss).
    /// The numerical-loss indicator guards precision; it is not a proved bound
    /// on the error against a continuous loading trajectory.
    pub estimated_source_work_error_j: f64,
    pub maximum_accepted_defect_ratio: f64,
}

/// Complete path result; refusing any trial leaves the original field unchanged.
pub struct HystereticPathCandidate {
    field: PlanarHystereticElectrostatics,
    report: HystereticPathReport,
}
impl HystereticPathCandidate {
    pub fn field(&self) -> &PlanarHystereticElectrostatics {
        &self.field
    }
    pub fn report(&self) -> &HystereticPathReport {
        &self.report
    }
    pub fn into_parts(self) -> (PlanarHystereticElectrostatics, HystereticPathReport) {
        (self.field, self.report)
    }
    /// Pass the complete accepted path to the original finite-enthalpy consumer.
    /// The native account carries cumulative per-cell physical heat and work;
    /// no field solve, receipt clone, extra deposition or history reset occurs.
    pub fn into_electrostatic_candidate(
        self,
    ) -> (HystereticElectrostaticCandidate, HystereticPathSummary) {
        (
            HystereticElectrostaticCandidate {
                field: self.field,
                step: self.report.account,
            },
            self.report.summary,
        )
    }
}

fn invalid(reason: &'static str) -> HystereticPathError {
    HystereticPathError::Invalid(reason)
}
fn finite(value: f64) -> Result<f64> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(invalid(
            "electric loading path arithmetic is unrepresentable",
        ))
    }
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<()> {
    if cancelled() {
        Err(HystereticPathError::Cancelled)
    } else {
        Ok(())
    }
}
fn norm(v: [f64; 3]) -> f64 {
    v[0].hypot(v[1]).hypot(v[2])
}
fn difference(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    std::array::from_fn(|i| a[i] - b[i])
}
fn ratio(error: f64, tolerance: f64) -> f64 {
    if error == 0.0 { 0.0 } else { error / tolerance }
}
fn interpolate(a: f64, b: f64, fraction: f64) -> Result<f64> {
    if fraction == 1.0 {
        return Ok(b);
    }
    if fraction == 0.0 {
        return Ok(a);
    }
    // Convex weights avoid overflow in b-a; voltages have already lost their
    // irrelevant common bias. Dyadic fractions and their complements are exact.
    finite((1.0 - fraction).mul_add(a, fraction * b))
}

struct Loading {
    start: Vec<ElectrodeControl>,
    end: Vec<ElectrodeControl>,
    rho_start: Vec<f64>,
    rho_end: Vec<f64>,
    final_reference_voltage_v: f64,
}
impl Loading {
    fn new(
        field: &PlanarHystereticElectrostatics,
        rho: &[f64],
        boundary: &[ElectrodeControl],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self> {
        if rho.len() != field.free_charge_density_c_m3.len()
            || boundary.len() != field.terminal_charge_c.len()
        {
            return Err(invalid(
                "complete volume-charge and electrode endpoints required",
            ));
        }
        let (reference, reference_v) = boundary
            .iter()
            .enumerate()
            .find_map(|(i, c)| {
                if let ElectrodeControl::Voltage(v) = *c {
                    Some((i, v))
                } else {
                    None
                }
            })
            .ok_or_else(|| invalid("an affine loading path needs a voltage reference"))?;
        finite(reference_v)?;
        let old_reference = field.terminal_voltage_v[reference];
        let mut start = Vec::with_capacity(boundary.len());
        let mut end = Vec::with_capacity(boundary.len());
        for (i, control) in boundary.iter().enumerate() {
            poll(cancelled)?;
            match *control {
                ElectrodeControl::Voltage(v) => {
                    finite(v)?;
                    start.push(ElectrodeControl::Voltage(finite(
                        field.terminal_voltage_v[i] - old_reference,
                    )?));
                    end.push(ElectrodeControl::Voltage(finite(v - reference_v)?));
                }
                ElectrodeControl::Charge(q) => {
                    start.push(ElectrodeControl::Charge(field.terminal_charge_c[i]));
                    end.push(ElectrodeControl::Charge(finite(q)?));
                }
            }
        }
        for value in rho {
            poll(cancelled)?;
            finite(*value)?;
        }
        Ok(Self {
            start,
            end,
            rho_start: field.free_charge_density_c_m3.clone(),
            rho_end: rho.to_vec(),
            final_reference_voltage_v: reference_v,
        })
    }
    fn at(
        &self,
        fraction: f64,
        boundary: &mut [ElectrodeControl],
        rho: &mut [f64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<()> {
        for (i, out) in boundary.iter_mut().enumerate() {
            poll(cancelled)?;
            *out = match (self.start[i], self.end[i]) {
                (ElectrodeControl::Voltage(a), ElectrodeControl::Voltage(b)) => {
                    ElectrodeControl::Voltage(interpolate(a, b, fraction)?)
                }
                (ElectrodeControl::Charge(a), ElectrodeControl::Charge(b)) => {
                    ElectrodeControl::Charge(interpolate(a, b, fraction)?)
                }
                _ => return Err(invalid("electrode modes changed inside an affine segment")),
            };
        }
        for (i, out) in rho.iter_mut().enumerate() {
            poll(cancelled)?;
            *out = interpolate(self.rho_start[i], self.rho_end[i], fraction)?;
        }
        Ok(())
    }
}

struct Work {
    controls: HystereticPathControls,
    field_solves: usize,
    operator_applications: usize,
}
impl Work {
    fn solve(
        &mut self,
        field: &PlanarHystereticElectrostatics,
        rho: &[f64],
        boundary: &[ElectrodeControl],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<HystereticElectrostaticCandidate> {
        poll(cancelled)?;
        if self.field_solves == self.controls.maximum_field_solves {
            return Err(HystereticPathError::FieldSolveBudget);
        }
        let remaining = self.controls.maximum_operator_applications - self.operator_applications;
        if remaining == 0 {
            return Err(HystereticPathError::OperatorBudget);
        }
        let mut controls = self.controls.field;
        controls.max_operator_applications = controls.max_operator_applications.min(remaining);
        self.field_solves += 1;
        let candidate = field
            .solve_candidate(rho, boundary, controls, cancelled)
            .map_err(|error| {
                if matches!(error, NonlinearElectrostaticError::OperatorBudget)
                    && remaining <= self.controls.field.max_operator_applications
                {
                    HystereticPathError::OperatorBudget
                } else {
                    error.into()
                }
            })?;
        let used = candidate.step().nonlinear.operator_applications;
        if used > remaining {
            return Err(HystereticPathError::OperatorBudget);
        }
        self.operator_applications += used;
        Ok(candidate)
    }
}

fn defects(
    start: &PlanarHystereticElectrostatics,
    coarse: &HystereticElectrostaticCandidate,
    half: &HystereticElectrostaticCandidate,
    fine: &HystereticElectrostaticCandidate,
    fraction: f64,
    c: HystereticPathControls,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<([f64; 3], f64)> {
    let mut polarization: f64 = 0.0;
    let mut maximum_ratio: f64 = 0.0;
    for cell in 0..start.problem.materials().len() {
        let old = start
            .cell_history(cell)
            .ok_or_else(|| invalid("missing accepted polarization"))?;
        let a = coarse
            .field()
            .cell_history(cell)
            .ok_or_else(|| invalid("missing coarse polarization"))?;
        let b = fine
            .field()
            .cell_history(cell)
            .ok_or_else(|| invalid("missing fine polarization"))?;
        for ((old, a), b) in old
            .polarizations_c_m2()
            .iter()
            .zip(a.polarizations_c_m2())
            .zip(b.polarizations_c_m2())
        {
            poll(cancelled)?;
            let error = finite(norm(difference(*a, *b)))?;
            let scale = finite(norm(difference(*a, *old)).max(norm(difference(*b, *old))))?;
            let allowance = finite(
                fraction * c.absolute_polarization_tolerance_c_m2 + c.relative_tolerance * scale,
            )?;
            polarization = polarization.max(error);
            maximum_ratio = maximum_ratio.max(ratio(error, allowance));
        }
    }
    let heat = finite(half.step().switching_heat_j + fine.step().switching_heat_j)?;
    let (mut heat_error, mut heat_correction) = (0.0, 0.0);
    for ((a, first), second) in coarse
        .step()
        .cells
        .iter()
        .zip(&half.step().cells)
        .zip(&fine.step().cells)
    {
        poll(cancelled)?;
        let local = finite(
            (a.switching_heat_j - finite(first.switching_heat_j + second.switching_heat_j)?).abs(),
        )?;
        add(&mut heat_error, &mut heat_correction, local)?;
    }
    heat_error = finite(heat_error + heat_correction)?;
    let heat_allowance = finite(
        fraction * c.absolute_heat_tolerance_j
            + c.relative_tolerance * coarse.step().switching_heat_j.max(heat),
    )?;
    let source_work =
        finite(half.step().endpoint_source_work_j + fine.step().endpoint_source_work_j)?;
    // A small full/half work difference can disappear beneath large bias work.
    // The native nonnegative remainder is formed stably from increments, so it
    // separately prevents accepting an unresolved endpoint discretization loss.
    let numerical =
        finite(half.step().numerical_dissipation_j + fine.step().numerical_dissipation_j)?;
    let work_error =
        finite((coarse.step().endpoint_source_work_j - source_work).abs())?.max(numerical);
    let work_scale = finite(
        half.step().endpoint_source_work_j.abs() + fine.step().endpoint_source_work_j.abs(),
    )?
    .max(coarse.step().endpoint_source_work_j.abs());
    let work_allowance =
        finite(fraction * c.absolute_work_tolerance_j + c.relative_tolerance * work_scale)?;
    maximum_ratio = maximum_ratio
        .max(ratio(heat_error, heat_allowance))
        .max(ratio(work_error, work_allowance));
    Ok(([polarization, heat_error, work_error], maximum_ratio))
}

struct Account {
    cells: Vec<[f64; 4]>,
    correction: Vec<[f64; 4]>,
    total: [f64; 5],
    total_correction: [f64; 5],
    last: Option<HystereticElectrostaticStep>,
}
impl Account {
    fn new(cells: usize) -> Self {
        Self {
            cells: vec![[0.0; 4]; cells],
            correction: vec![[0.0; 4]; cells],
            total: [0.0; 5],
            total_correction: [0.0; 5],
            last: None,
        }
    }
    fn push(
        &mut self,
        step: HystereticElectrostaticStep,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<()> {
        for (i, cell) in step.cells.iter().enumerate() {
            poll(cancelled)?;
            let values = [
                cell.stored_energy_change_j,
                cell.switching_heat_j,
                cell.numerical_dissipation_j,
                cell.endpoint_work_j,
            ];
            for (k, value) in values.into_iter().enumerate() {
                add(&mut self.cells[i][k], &mut self.correction[i][k], value)?;
            }
        }
        for (i, value) in [
            step.stored_change_j,
            step.switching_heat_j,
            step.numerical_dissipation_j,
            step.endpoint_field_work_j,
            step.endpoint_source_work_j,
        ]
        .into_iter()
        .enumerate()
        {
            add(&mut self.total[i], &mut self.total_correction[i], value)?;
        }
        self.last = Some(step);
        Ok(())
    }
    fn finish(
        mut self,
        stored_before_j: f64,
        tolerance: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<HystereticElectrostaticStep> {
        let mut result = self
            .last
            .ok_or_else(|| invalid("loading path has no accepted substeps"))?;
        for (i, cell) in result.cells.iter_mut().enumerate() {
            poll(cancelled)?;
            for k in 0..4 {
                self.cells[i][k] = finite(self.cells[i][k] + self.correction[i][k])?;
            }
            let [stored, heat, numerical, work] = self.cells[i];
            cell.stored_energy_change_j = stored;
            cell.switching_heat_j = heat;
            cell.numerical_dissipation_j = numerical;
            cell.endpoint_work_j = work;
            cell.energy_residual_j = finite(work - stored - heat - numerical)?;
            if cell.energy_residual_j.abs() > tolerance {
                return Err(invalid("accumulated cell switching work did not close"));
            }
        }
        for i in 0..5 {
            self.total[i] = finite(self.total[i] + self.total_correction[i])?;
        }
        let [stored, heat, numerical, field_work, source_work] = self.total;
        result.stored_before_j = stored_before_j;
        result.stored_change_j = stored;
        result.switching_heat_j = heat;
        result.numerical_dissipation_j = numerical;
        result.endpoint_field_work_j = field_work;
        result.endpoint_source_work_j = source_work;
        result.energy_residual_j = finite(source_work - stored - heat - numerical)?;
        result.field_source_work_residual_j = finite(source_work - field_work)?;
        if result.energy_residual_j.abs() > tolerance
            || result.field_source_work_residual_j.abs() > tolerance
        {
            return Err(invalid(
                "accumulated electrode work, storage, heat and numerical loss did not close",
            ));
        }
        Ok(result)
    }
}

impl PlanarHystereticElectrostatics {
    /// Continue from the accepted field along one affine electrode/volume-charge
    /// segment. Each target electrode's control mode is fixed over the segment;
    /// its initial value is the actually accepted voltage or charge. Voltages
    /// are interpolated relative to the first voltage-controlled electrode,
    /// preserving small fields beneath large common biases. The returned field
    /// restores the target's authored reference voltage.
    ///
    /// Absolute error allowances are scaled by interval fraction; relative ones
    /// scale local polarization increments and physical heat/work. The estimates
    /// are summed over accepted intervals and remain empirical discretization
    /// estimates. No extrapolated or averaged polarization is published.
    pub fn solve_affine_path_candidate(
        &self,
        target_free_charge_density_c_m3: &[f64],
        target_controls: &[ElectrodeControl],
        controls: HystereticPathControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<HystereticPathCandidate> {
        poll(cancelled)?;
        let absolute = [
            controls.absolute_polarization_tolerance_c_m2,
            controls.absolute_heat_tolerance_j,
            controls.absolute_work_tolerance_j,
        ];
        if !controls.relative_tolerance.is_finite()
            || controls.relative_tolerance < 0.0
            || absolute.iter().any(|v| {
                !v.is_finite() || *v < 0.0 || (*v == 0.0 && controls.relative_tolerance == 0.0)
            })
            || !controls.energy_tolerance_j.is_finite()
            || controls.energy_tolerance_j < 0.0
            || controls.maximum_field_solves == 0
            || controls.maximum_operator_applications == 0
            || controls.maximum_refinement_depth > 52
        {
            return Err(invalid(
                "finite path tolerances, positive work caps and depth at most 52 required",
            ));
        }
        let loading = Loading::new(
            self,
            target_free_charge_density_c_m3,
            target_controls,
            cancelled,
        )?;
        let mut boundary = loading.start.clone();
        let mut rho = loading.rho_start.clone();
        let mut current = self.clone();
        let mut work = Work {
            controls,
            field_solves: 0,
            operator_applications: 0,
        };
        let mut account = Account::new(self.problem.materials().len());
        // Only interval descriptors are stacked. Reuse a rejected interval's
        // first half as the left child's coarse solve; its starting history is
        // identical. The right child must solve again from the refined left.
        let mut pending = vec![(0.0, 1.0, 0usize)];
        let mut cached_coarse = None;
        let (mut accepted, mut rejected, mut deepest) = (0usize, 0usize, 0usize);
        let (mut estimates, mut correction) = ([0.0; 3], [0.0; 3]);
        let mut maximum_ratio: f64 = 0.0;
        while let Some((a, b, depth)) = pending.pop() {
            poll(cancelled)?;
            deepest = deepest.max(depth);
            let midpoint = f64::midpoint(a, b);
            loading.at(b, &mut boundary, &mut rho, cancelled)?;
            let coarse = match cached_coarse.take() {
                Some(value) => value,
                None => work.solve(&current, &rho, &boundary, cancelled)?,
            };
            loading.at(midpoint, &mut boundary, &mut rho, cancelled)?;
            let half = work.solve(&current, &rho, &boundary, cancelled)?;
            loading.at(b, &mut boundary, &mut rho, cancelled)?;
            let fine = work.solve(half.field(), &rho, &boundary, cancelled)?;
            let (errors, ratio) =
                defects(&current, &coarse, &half, &fine, b - a, controls, cancelled)?;
            if ratio > 1.0 {
                if depth == controls.maximum_refinement_depth {
                    return Err(HystereticPathError::RefinementLimit {
                        fraction: [a, b],
                        depth,
                        maximum_defect_ratio: ratio,
                    });
                }
                rejected += 1;
                pending.push((midpoint, b, depth + 1));
                pending.push((a, midpoint, depth + 1));
                cached_coarse = Some(half);
            } else {
                poll(cancelled)?;
                let (_, half_step) = half.into_parts();
                let (next, fine_step) = fine.into_parts();
                account.push(half_step, cancelled)?;
                account.push(fine_step, cancelled)?;
                for i in 0..3 {
                    add(&mut estimates[i], &mut correction[i], errors[i])?;
                }
                maximum_ratio = maximum_ratio.max(ratio);
                accepted = accepted
                    .checked_add(2)
                    .ok_or_else(|| invalid("path substep count overflowed"))?;
                current = next;
            }
        }
        for i in 0..3 {
            estimates[i] = finite(estimates[i] + correction[i])?;
        }
        let report = HystereticPathReport {
            account: account.finish(
                self.stored_energy_j,
                controls.energy_tolerance_j,
                cancelled,
            )?,
            summary: HystereticPathSummary {
                field_solves: work.field_solves,
                operator_applications: work.operator_applications,
                accepted_substeps: accepted,
                rejected_intervals: rejected,
                maximum_refinement_depth: deepest,
                estimated_polarization_error_c_m2: estimates[0],
                estimated_switching_heat_error_j: estimates[1],
                estimated_source_work_error_j: estimates[2],
                maximum_accepted_defect_ratio: maximum_ratio,
            },
        };
        current.reference_voltage_v = loading.final_reference_voltage_v;
        poll(cancelled)?;
        Ok(HystereticPathCandidate {
            field: current,
            report,
        })
    }
}
