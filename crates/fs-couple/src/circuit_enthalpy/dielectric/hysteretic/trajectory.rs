//! Adaptive thermal chronology along one authored affine electrical segment.
//!
//! Full and two-half trials use the original coupled endpoint solve. Only the
//! two-half histories can be accepted; neither state nor physical heat is
//! extrapolated. Electrical memory/work and per-cell thermal energy have
//! independent error indicators. This resolves switching heat during cooling,
//! including holds with no electrical change. It remains first-order backward
//! Euler, with empirical local estimates, not a certified global time bound.
//! Split at waveform turns, control changes and thermal/source events.

use super::*;

/// Independent accuracy and global work limits for one complete time segment.
#[derive(Clone, Copy, Debug)]
pub struct HystereticThermalPathControls {
    pub thermal: CircuitEnthalpySolve,
    /// Its native field policy is used for each coupled endpoint. Its P/Q/work
    /// tolerances, depth, solve and operator limits apply across the WHOLE path.
    pub path: HystereticPathControls,
    /// L1 allowance on cell-energy step-doubling defects [J], duration weighted.
    pub absolute_thermal_tolerance_j: f64,
    /// Scales local cell-energy changes, never absolute/reference enthalpy.
    pub relative_thermal_tolerance: f64,
    /// All native thermal driver evaluations, including rejected steps/stencils.
    pub maximum_thermal_evaluations: usize,
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct HystereticThermalPathSummary {
    pub attempted_steps: usize,
    pub field_operator_applications: usize,
    pub thermal_evaluations: usize,
    pub accepted_steps: usize,
    pub rejected_intervals: usize,
    pub maximum_refinement_depth: usize,
    pub estimated_polarization_error_c_m2: f64,
    pub estimated_switching_heat_error_j: f64,
    pub estimated_source_work_error_j: f64,
    pub estimated_thermal_energy_error_j: f64,
    pub maximum_accepted_defect_ratio: f64,
}

/// Both complete histories and accepted-only cumulative physical accounts.
pub struct HystereticThermalPathCandidate<'t, 'n, 'a, S = MaterialNetworkSample> {
    thermal: DepositionSimulation<'t>,
    field: ThermalHystereticElectrostatics<'n, 'a>,
    report: CircuitEnthalpyReport<S, HystereticElectrostaticStep>,
    summary: HystereticThermalPathSummary,
}
impl<'t, 'n, 'a, S> HystereticThermalPathCandidate<'t, 'n, 'a, S> {
    pub fn thermal(&self) -> &DepositionSimulation<'t> {
        &self.thermal
    }
    pub fn field(&self) -> &ThermalHystereticElectrostatics<'n, 'a> {
        &self.field
    }
    /// All energy/cell accounts span the complete accepted path. `source`,
    /// nonlinear reports and charge residuals are the FINAL substep's original
    /// sample. `maximum_cell_residual_j` is the maximum over accepted substeps;
    /// `circuit_trials` includes all attempted thermal evaluations.
    pub fn report(&self) -> &CircuitEnthalpyReport<S, HystereticElectrostaticStep> {
        &self.report
    }
    pub fn summary(&self) -> &HystereticThermalPathSummary {
        &self.summary
    }
    pub fn into_parts(
        self,
    ) -> (
        DepositionSimulation<'t>,
        ThermalHystereticElectrostatics<'n, 'a>,
        CircuitEnthalpyReport<S, HystereticElectrostaticStep>,
        HystereticThermalPathSummary,
    ) {
        (self.thermal, self.field, self.report, self.summary)
    }
}

#[derive(Clone, Copy, Default)]
struct Sum(f64, f64);
impl Sum {
    fn push<E>(&mut self, value: f64) -> Result<(), Error<E>> {
        let mut pair = (self.0, self.1);
        add(&mut pair, value)?;
        (self.0, self.1) = pair;
        Ok(())
    }
    fn value<E>(self) -> Result<f64, Error<E>> {
        finite(self.0 + self.1)
    }
}
fn sum<E>(values: impl IntoIterator<Item = f64>) -> Result<f64, Error<E>> {
    let mut result = Sum::default();
    for value in values {
        result.push(value)?;
    }
    result.value()
}
fn norm(v: [f64; 3]) -> f64 {
    v[0].hypot(v[1]).hypot(v[2])
}
fn difference(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    std::array::from_fn(|i| a[i] - b[i])
}
fn ratio(error: f64, allowance: f64) -> f64 {
    if error == 0.0 { 0.0 } else { error / allowance }
}
fn interpolate<E>(a: f64, b: f64, fraction: f64) -> Result<f64, Error<E>> {
    if fraction == 1.0 {
        Ok(b)
    } else {
        finite((1.0 - fraction).mul_add(a, fraction * b))
    }
}

struct Loading<'s> {
    field: &'s PlanarHystereticElectrostatics,
    rho: &'s [f64],
    boundary: &'s [ElectrodeControl],
    reference: usize,
    reference_v: f64,
}
impl<'s> Loading<'s> {
    fn new<E>(
        field: &'s PlanarHystereticElectrostatics,
        rho: &'s [f64],
        boundary: &'s [ElectrodeControl],
        cancel: &RefCell<impl FnMut() -> bool>,
    ) -> Result<Self, Error<E>> {
        if rho.len() != field.free_charge_density_c_m3().len()
            || boundary.len() != field.terminal_charge_c().len()
        {
            return Err(invalid(
                "complete affine charge and electrode endpoints required",
            ));
        }
        let (reference, reference_v) = boundary
            .iter()
            .enumerate()
            .find_map(|(i, b)| {
                if let ElectrodeControl::Voltage(v) = b {
                    Some((i, *v))
                } else {
                    None
                }
            })
            .ok_or_else(|| invalid("a voltage reference is required for thermal continuation"))?;
        for b in boundary {
            check(cancel)?;
            let (ElectrodeControl::Voltage(v) | ElectrodeControl::Charge(v)) = b;
            finite::<E>(*v)?;
        }
        for &v in rho {
            check(cancel)?;
            finite::<E>(v)?;
        }
        Ok(Self {
            field,
            rho,
            boundary,
            reference,
            reference_v,
        })
    }
    fn at<E>(
        &self,
        fraction: f64,
        rho: &mut [f64],
        boundary: &mut [ElectrodeControl],
        cancel: &RefCell<impl FnMut() -> bool>,
    ) -> Result<(), Error<E>> {
        for (i, b) in boundary.iter_mut().enumerate() {
            check(cancel)?;
            // Intermediate solves use a zero reference. Restore the authored
            // final gauge through the original solve, without rounding every
            // small intermediate voltage beneath a potentially large bias.
            *b = if fraction == 1.0 {
                self.boundary[i]
            } else {
                match self.boundary[i] {
                    ElectrodeControl::Voltage(v) => ElectrodeControl::Voltage(interpolate(
                        finite(
                            self.field.terminal_voltage_v()[i]
                                - self.field.terminal_voltage_v()[self.reference],
                        )?,
                        finite(v - self.reference_v)?,
                        fraction,
                    )?),
                    ElectrodeControl::Charge(q) => ElectrodeControl::Charge(interpolate(
                        self.field.terminal_charge_c()[i],
                        q,
                        fraction,
                    )?),
                }
            };
        }
        for (i, out) in rho.iter_mut().enumerate() {
            check(cancel)?;
            *out = interpolate(
                self.field.free_charge_density_c_m3()[i],
                self.rho[i],
                fraction,
            )?;
        }
        Ok(())
    }
}

fn validate<E>(controls: HystereticThermalPathControls) -> Result<(), Error<E>> {
    let p = controls.path;
    if !p.relative_tolerance.is_finite()
        || p.relative_tolerance < 0.0
        || [
            p.absolute_polarization_tolerance_c_m2,
            p.absolute_heat_tolerance_j,
            p.absolute_work_tolerance_j,
        ]
        .iter()
        .any(|v| !v.is_finite() || *v < 0.0 || (*v == 0.0 && p.relative_tolerance == 0.0))
        || !p.energy_tolerance_j.is_finite()
        || p.energy_tolerance_j < 0.0
        || p.maximum_field_solves == 0
        || p.maximum_operator_applications == 0
        || p.maximum_refinement_depth > 52
        || !controls.absolute_thermal_tolerance_j.is_finite()
        || controls.absolute_thermal_tolerance_j < 0.0
        || !controls.relative_thermal_tolerance.is_finite()
        || controls.relative_thermal_tolerance < 0.0
        || (controls.absolute_thermal_tolerance_j == 0.0
            && controls.relative_thermal_tolerance == 0.0)
        || controls.maximum_thermal_evaluations == 0
    {
        return Err(invalid(
            "finite coupled path tolerances, positive work caps and depth at most 52 required",
        ));
    }
    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn defects<S, E>(
    start: &ThermalHystereticElectrostatics<'_, '_>,
    initial: &DepositionSimulation<'_>,
    coarse: &HystereticDielectricEnthalpyCandidate<'_, '_, '_, S>,
    half: &HystereticDielectricEnthalpyCandidate<'_, '_, '_, S>,
    fine: &HystereticDielectricEnthalpyCandidate<'_, '_, '_, S>,
    fraction: f64,
    controls: HystereticThermalPathControls,
    cancel: &RefCell<impl FnMut() -> bool>,
) -> Result<([f64; 4], f64), Error<E>> {
    let p = controls.path;
    let (mut polarization, mut maximum): (f64, f64) = (0.0, 0.0);
    for cell in 0..start.field.problem().materials().len() {
        let old = start
            .field
            .cell_history(cell)
            .ok_or_else(|| invalid("missing accepted polarization"))?;
        let a = coarse
            .field
            .field
            .cell_history(cell)
            .ok_or_else(|| invalid("missing coarse polarization"))?;
        let b = fine
            .field
            .field
            .cell_history(cell)
            .ok_or_else(|| invalid("missing fine polarization"))?;
        for ((old, a), b) in old
            .polarizations_c_m2()
            .iter()
            .zip(a.polarizations_c_m2())
            .zip(b.polarizations_c_m2())
        {
            check(cancel)?;
            let error = finite(norm(difference(*a, *b)))?;
            let scale = finite(norm(difference(*a, *old)).max(norm(difference(*b, *old))))?;
            polarization = polarization.max(error);
            maximum = maximum.max(ratio(
                error,
                finite(
                    fraction * p.absolute_polarization_tolerance_c_m2
                        + p.relative_tolerance * scale,
                )?,
            ));
        }
    }
    let (a, first, b) = (
        &coarse.report.electrical,
        &half.report.electrical,
        &fine.report.electrical,
    );
    let mut heat_error = Sum::default();
    for ((a, first), b) in a.cells.iter().zip(&first.cells).zip(&b.cells) {
        check(cancel)?;
        heat_error.push(
            finite(a.switching_heat_j - sum([first.switching_heat_j, b.switching_heat_j])?)?.abs(),
        )?;
    }
    let heat_error = heat_error.value()?;
    let heat_scale = a
        .switching_heat_j
        .max(sum([first.switching_heat_j, b.switching_heat_j])?);
    maximum = maximum.max(ratio(
        heat_error,
        finite(fraction * p.absolute_heat_tolerance_j + p.relative_tolerance * heat_scale)?,
    ));
    let numerical = sum([first.numerical_dissipation_j, b.numerical_dissipation_j])?;
    let work_error = sum([
        a.endpoint_source_work_j,
        -first.endpoint_source_work_j,
        -b.endpoint_source_work_j,
    ])?
    .abs()
    .max(numerical);
    let work_scale = a.endpoint_source_work_j.abs().max(sum([
        first.endpoint_source_work_j.abs(),
        b.endpoint_source_work_j.abs(),
    ])?);
    maximum = maximum.max(ratio(
        work_error,
        finite(fraction * p.absolute_work_tolerance_j + p.relative_tolerance * work_scale)?,
    ));
    let (mut thermal_error, mut thermal_scale) = (Sum::default(), Sum::default());
    if initial.cells().len() != coarse.thermal.cells().len()
        || initial.cells().len() != fine.thermal.cells().len()
    {
        return Err(invalid("thermal topology changed inside an affine segment"));
    }
    for (&id, old) in initial.cells() {
        check(cancel)?;
        let a = coarse
            .thermal
            .cells()
            .get(&id)
            .ok_or_else(|| invalid("coarse thermal cell missing"))?;
        let b = fine
            .thermal
            .cells()
            .get(&id)
            .ok_or_else(|| invalid("fine thermal cell missing"))?;
        if a.mass_kg != old.mass_kg || b.mass_kg != old.mass_kg {
            return Err(invalid("cell inventory changed inside an affine segment"));
        }
        thermal_error.push(finite(old.mass_kg * (a.specific_h_j_kg - b.specific_h_j_kg))?.abs())?;
        thermal_scale.push(finite(
            old.mass_kg
                * (a.specific_h_j_kg - old.specific_h_j_kg)
                    .abs()
                    .max((b.specific_h_j_kg - old.specific_h_j_kg).abs()),
        )?)?;
    }
    let thermal_error = thermal_error.value()?;
    maximum = maximum.max(ratio(
        thermal_error,
        finite(
            fraction * controls.absolute_thermal_tolerance_j
                + controls.relative_thermal_tolerance * thermal_scale.value::<E>()?,
        )?,
    ));
    Ok((
        [polarization, heat_error, work_error, thermal_error],
        maximum,
    ))
}

struct Account<S> {
    cells: Vec<[Sum; 4]>,
    totals: [Sum; 11],
    maximum_cell_residual_j: f64,
    last: Option<CircuitEnthalpyReport<S, HystereticElectrostaticStep>>,
}
impl<S> Account<S> {
    fn new(cells: usize) -> Self {
        Self {
            cells: vec![[Sum::default(); 4]; cells],
            totals: [Sum::default(); 11],
            maximum_cell_residual_j: 0.0,
            last: None,
        }
    }
    fn push<E>(
        &mut self,
        report: CircuitEnthalpyReport<S, HystereticElectrostaticStep>,
        cancel: &RefCell<impl FnMut() -> bool>,
    ) -> Result<(), Error<E>> {
        for (i, cell) in report.electrical.cells.iter().enumerate() {
            check(cancel)?;
            for (sum, value) in self.cells[i].iter_mut().zip([
                cell.stored_energy_change_j,
                cell.switching_heat_j,
                cell.numerical_dissipation_j,
                cell.endpoint_work_j,
            ]) {
                sum.push(value)?;
            }
        }
        let e = &report.electrical;
        for (sum, value) in self.totals.iter_mut().zip([
            e.stored_change_j,
            e.switching_heat_j,
            e.numerical_dissipation_j,
            e.endpoint_field_work_j,
            e.endpoint_source_work_j,
            report.thermal.stored_energy_change_j,
            report.thermal.external_heat_j,
            report.deposited_joule_j,
            report.non_electrical_heat_j,
            report.thermal.residual_j,
            report.total_energy_residual_j,
        ]) {
            sum.push(value)?;
        }
        self.maximum_cell_residual_j = self
            .maximum_cell_residual_j
            .max(report.maximum_cell_residual_j);
        self.last = Some(report);
        Ok(())
    }
    fn finish<E>(
        self,
        stored_before_j: f64,
        controls: HystereticThermalPathControls,
        summary: HystereticThermalPathSummary,
        cancel: &RefCell<impl FnMut() -> bool>,
    ) -> Result<CircuitEnthalpyReport<S, HystereticElectrostaticStep>, Error<E>> {
        let mut report = self
            .last
            .ok_or_else(|| invalid("no accepted coupled substeps"))?;
        for (cell, values) in report.electrical.cells.iter_mut().zip(self.cells) {
            check(cancel)?;
            let [w, q, d, work] = values;
            cell.stored_energy_change_j = w.value()?;
            cell.switching_heat_j = q.value()?;
            cell.numerical_dissipation_j = d.value()?;
            cell.endpoint_work_j = work.value()?;
            cell.energy_residual_j = sum([
                cell.endpoint_work_j,
                -cell.stored_energy_change_j,
                -cell.switching_heat_j,
                -cell.numerical_dissipation_j,
            ])?;
            if cell.energy_residual_j.abs() > controls.path.energy_tolerance_j {
                return Err(invalid("cumulative cell switching energy did not close"));
            }
        }
        let [
            w,
            q,
            d,
            work,
            supplied,
            thermal,
            external,
            deposited,
            other,
            thermal_residual,
            total_residual,
        ] = self.totals;
        let e = &mut report.electrical;
        e.stored_before_j = stored_before_j;
        e.stored_change_j = w.value()?;
        e.switching_heat_j = q.value()?;
        e.numerical_dissipation_j = d.value()?;
        e.endpoint_field_work_j = work.value()?;
        e.endpoint_source_work_j = supplied.value()?;
        e.energy_residual_j = sum([
            e.endpoint_source_work_j,
            -e.stored_change_j,
            -e.switching_heat_j,
            -e.numerical_dissipation_j,
        ])?;
        e.field_source_work_residual_j =
            finite(e.endpoint_source_work_j - e.endpoint_field_work_j)?;
        report.thermal.stored_energy_change_j = thermal.value()?;
        report.thermal.external_heat_j = external.value()?;
        // Retain the original cell-relative audit; subtracting rounded aggregate
        // energies can erase small heat beside large reference enthalpies.
        report.thermal.residual_j = thermal_residual.value()?;
        report.thermal.substeps = summary.accepted_steps;
        report.deposited_joule_j = deposited.value()?;
        report.non_electrical_heat_j = other.value()?;
        report.total_energy_residual_j = total_residual.value()?;
        report.circuit_trials = summary.thermal_evaluations;
        report.maximum_cell_residual_j = self.maximum_cell_residual_j;
        if e.energy_residual_j.abs() > controls.path.energy_tolerance_j
            || e.field_source_work_residual_j.abs() > controls.path.energy_tolerance_j
            || (report.deposited_joule_j - e.switching_heat_j).abs()
                > controls.thermal.total_energy_tolerance_j
            || report.thermal.residual_j.abs() > controls.thermal.thermal_energy_tolerance_j
            || report.total_energy_residual_j.abs() > controls.thermal.total_energy_tolerance_j
        {
            return Err(invalid(
                "cumulative electrical and thermal energy did not close",
            ));
        }
        Ok(report)
    }
}

impl<'n, 'a> ThermalHystereticElectrostatics<'n, 'a> {
    /// Advance both histories along one affine charge/electrode segment in time.
    /// Target control modes are fixed; initial values are the accepted physical
    /// voltages/charges. The candidate owns already-deposited physical heat.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_affine_candidate<'t>(
        &self,
        thermal: &DepositionSimulation<'t>,
        target_s: f64,
        rho: &[f64],
        boundary: &[ElectrodeControl],
        controls: HystereticThermalPathControls,
        cancelled: impl FnMut() -> bool,
    ) -> Result<HystereticThermalPathCandidate<'t, 'n, 'a>, Error> {
        let cancel = RefCell::new(cancelled);
        self.advance_resolved(
            thermal,
            target_s,
            rho,
            boundary,
            controls,
            |field, state, target, rho, boundary, heat, electric, guess| {
                field.step_candidate(state, target, rho, boundary, heat, electric, guess, || {
                    (*cancel.borrow_mut())()
                })
            },
            &cancel,
        )
    }

    /// The same chronology through the original complete HeatAssembly. Every
    /// attempted endpoint retains its original source, domain and closed-energy
    /// gates, including radiation. The final report's source is the last sample.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_affine_candidate_with_heat_assembly<'t, A: HeatAssembly + ?Sized>(
        &self,
        cx: &Cx<'_>,
        assembly: &A,
        thermal: &DepositionSimulation<'t>,
        target_s: f64,
        rho: &[f64],
        boundary: &[ElectrodeControl],
        controls: HystereticThermalPathControls,
    ) -> Result<HystereticThermalPathCandidate<'t, 'n, 'a, A::Sample>, Error<A::Error>> {
        let cancel = RefCell::new(|| cx.checkpoint().is_err());
        self.advance_resolved(
            thermal,
            target_s,
            rho,
            boundary,
            controls,
            |field, state, target, rho, boundary, heat, electric, guess| {
                field.step_candidate_with_heat_assembly(
                    cx, assembly, state, target, rho, boundary, heat, electric, guess,
                )
            },
            &cancel,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn advance_resolved<'t, S, E>(
        &self,
        initial: &DepositionSimulation<'t>,
        target_s: f64,
        target_rho: &[f64],
        target_boundary: &[ElectrodeControl],
        controls: HystereticThermalPathControls,
        mut step: impl FnMut(
            &Self,
            &DepositionSimulation<'t>,
            f64,
            &[f64],
            &[ElectrodeControl],
            CircuitEnthalpySolve,
            NonlinearElectrostaticControls,
            &[f64],
        )
            -> Result<HystereticDielectricEnthalpyCandidate<'t, 'n, 'a, S>, Error<E>>,
        cancel: &RefCell<impl FnMut() -> bool>,
    ) -> Result<HystereticThermalPathCandidate<'t, 'n, 'a, S>, Error<E>> {
        check(cancel)?;
        validate(controls)?;
        if initial.time_s().to_bits() != self.time_s.to_bits()
            || !target_s.is_finite()
            || target_s <= self.time_s
        {
            return Err(invalid(
                "matching accepted clocks and a later finite target are required",
            ));
        }
        let duration = finite(target_s - self.time_s)?;
        let loading = Loading::new(&self.field, target_rho, target_boundary, cancel)?;
        let mut rho = target_rho.to_vec();
        let mut boundary = target_boundary.to_vec();
        let guess = vec![0.0; initial.cells().len()];
        admit_driver_unknowns(initial, controls.thermal, &guess)?;
        let mut summary = HystereticThermalPathSummary::default();
        let mut account = Account::new(self.field.problem().materials().len());
        let mut estimates = [Sum::default(); 4];
        let mut pending = vec![(0.0, 1.0, 0usize)];
        let mut current: Option<(DepositionSimulation<'t>, Self)> = None;
        let mut cached_coarse = None;
        while let Some((a, b, depth)) = pending.pop() {
            check(cancel)?;
            summary.maximum_refinement_depth = summary.maximum_refinement_depth.max(depth);
            let midpoint = f64::midpoint(a, b);
            let (thermal, field) = current
                .as_ref()
                .map(|(t, f)| (t, f))
                .unwrap_or((initial, self));
            let middle_s = finite(duration.mul_add(midpoint, self.time_s))?;
            let end_s = if b == 1.0 {
                target_s
            } else {
                finite(duration.mul_add(b, self.time_s))?
            };
            if middle_s <= thermal.time_s() || middle_s >= end_s {
                return Err(invalid(
                    "coupled time subdivision is no longer representable",
                ));
            }
            let mut trial = |field: &Self,
                             state: &DepositionSimulation<'t>,
                             time,
                             fraction,
                             width| {
                check(cancel)?;
                if summary.attempted_steps == controls.path.maximum_field_solves {
                    return Err(path(HystereticPathError::FieldSolveBudget));
                }
                let ops = controls.path.maximum_operator_applications
                    - summary.field_operator_applications;
                if ops == 0 {
                    return Err(path(HystereticPathError::OperatorBudget));
                }
                let evaluations =
                    controls.maximum_thermal_evaluations - summary.thermal_evaluations;
                if evaluations == 0 {
                    return Err(CircuitEnthalpyError::CircuitBudget);
                }
                let mut heat = controls.thermal;
                heat.maximum_circuit_trials = heat.maximum_circuit_trials.min(evaluations);
                heat.cell_energy_tolerance_j *= width;
                heat.thermal_energy_tolerance_j *= width;
                heat.total_energy_tolerance_j *= width;
                heat.newton.absolute_tolerance *= width;
                let mut electric = controls.path.field;
                electric.max_operator_applications = electric.max_operator_applications.min(ops);
                electric.energy_tolerance_j = electric
                    .energy_tolerance_j
                    .min(controls.path.energy_tolerance_j * width);
                loading.at(fraction, &mut rho, &mut boundary, cancel)?;
                summary.attempted_steps += 1;
                let result = step(field, state, time, &rho, &boundary, heat, electric, &guess)
                    .map_err(|error| match error {
                        CircuitEnthalpyError::Source(HystereticDielectricThermalError::Field(
                            NonlinearElectrostaticError::OperatorBudget,
                        )) if ops <= controls.path.field.max_operator_applications => {
                            path(HystereticPathError::OperatorBudget)
                        }
                        error => error,
                    })?;
                let used = result.report.electrical.nonlinear.operator_applications;
                if used > ops {
                    return Err(path(HystereticPathError::OperatorBudget));
                }
                if result.report.circuit_trials > evaluations {
                    return Err(CircuitEnthalpyError::CircuitBudget);
                }
                summary.field_operator_applications += used;
                summary.thermal_evaluations += result.report.circuit_trials;
                Ok(result)
            };
            let coarse = match cached_coarse.take() {
                Some(candidate) => candidate,
                None => trial(field, thermal, end_s, b, b - a)?,
            };
            let half = trial(field, thermal, middle_s, midpoint, midpoint - a)?;
            let fine = trial(&half.field, &half.thermal, end_s, b, b - midpoint)?;
            let (errors, ratio) = defects(
                field,
                thermal,
                &coarse,
                &half,
                &fine,
                b - a,
                controls,
                cancel,
            )?;
            if ratio > 1.0 {
                if depth == controls.path.maximum_refinement_depth {
                    return Err(path(HystereticPathError::RefinementLimit {
                        fraction: [a, b],
                        depth,
                        maximum_defect_ratio: ratio,
                    }));
                }
                summary.rejected_intervals += 1;
                pending.push((midpoint, b, depth + 1));
                pending.push((a, midpoint, depth + 1));
                cached_coarse = Some(half);
            } else {
                check(cancel)?;
                let (_, _, first, _) = half.into_parts();
                let (thermal, field, second, _) = fine.into_parts();
                account.push(first, cancel)?;
                account.push(second, cancel)?;
                for (sum, error) in estimates.iter_mut().zip(errors) {
                    sum.push(error)?;
                }
                summary.maximum_accepted_defect_ratio =
                    summary.maximum_accepted_defect_ratio.max(ratio);
                summary.accepted_steps = summary
                    .accepted_steps
                    .checked_add(2)
                    .ok_or_else(|| invalid("coupled accepted-step count overflowed"))?;
                current = Some((thermal, field));
            }
        }
        summary.estimated_polarization_error_c_m2 = estimates[0].value()?;
        summary.estimated_switching_heat_error_j = estimates[1].value()?;
        summary.estimated_source_work_error_j = estimates[2].value()?;
        summary.estimated_thermal_energy_error_j = estimates[3].value()?;
        let report = account.finish(self.field.stored_energy_j(), controls, summary, cancel)?;
        let (thermal, field) = current.ok_or_else(|| invalid("no accepted coupled endpoint"))?;
        check(cancel)?;
        Ok(HystereticThermalPathCandidate {
            thermal,
            field,
            report,
            summary,
        })
    }
}
