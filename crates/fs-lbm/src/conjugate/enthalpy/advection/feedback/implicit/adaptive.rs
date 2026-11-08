//! Local time-error control for the original thermoviscous BE endpoint solve.
//! Coarse and fine trials share numerical work, never accepted physical state.

use super::*;
use fs_conduction::enthalpy::conduction::implicit::adaptive::ConductionTimeOrder;

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AdaptiveImplicitFlowControls {
    /// Original nonlinear, finite-difference and per-interval energy controls.
    pub step: ImplicitFlowControls,
    pub initial_dt_s: f64,
    /// Floor for a complete comparison; the two fine intervals are shorter.
    pub minimum_dt_s: f64,
    pub maximum_dt_s: f64,
    pub maximum_trials: usize,
    /// Mandatory original-h accuracy, including isothermal latent storage.
    pub enthalpy_absolute_tolerance_j_kg: f64,
    pub temperature_absolute_tolerance_k: Option<f64>,
    /// Scales changes from each comparison's start, never absolute h datums.
    pub relative_tolerance: f64,
    /// Whole-advance energy allowance, including the sum of absolute accepted
    /// algebraic energy residuals. It does not multiply by the substep count.
    pub energy_tolerance_j: f64,
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct ImplicitFlowTimeError {
    pub maximum_enthalpy_difference_j_kg: f64,
    pub maximum_temperature_difference_k: Option<f64>,
    /// Corrected fine local error / (absolute + relative * physical change).
    pub normalized: f64,
}

#[derive(Debug)]
pub struct AdaptiveImplicitFlowReport<F, T = FluidTransportState> {
    pub thermal: EnergyAudit,
    /// The original work/physical-loss owner, integrated only on fine steps.
    pub mechanical: ViscousHeatingAudit,
    pub deposited_heat_j: f64,
    pub net_stream_enthalpy_j: f64,
    pub other_thermal_input_j: f64,
    pub combined_residual_j: f64,
    pub absolute_accepted_energy_residual_j: f64,
    /// Includes every coarse, fine, rejected and finite-difference flow solve.
    pub momentum_solves: usize,
    pub trials: usize,
    pub accepted_intervals: usize,
    pub accuracy_rejections: usize,
    pub maximum_accepted_normalized_error: f64,
    pub last_error: ImplicitFlowTimeError,
    pub suggested_dt_s: f64,
    /// Final original coupled endpoint. Its one-step accounts are already
    /// included above; they are not additional physical work or heat.
    pub last_step: ImplicitFlowReport<F, T>,
}

#[derive(Debug)]
pub enum AdaptiveImplicitFlowError {
    /// Original source, momentum, nonlinear, budget and cancellation failures.
    Solve(ImplicitFlowError),
    Invalid(&'static str),
    MinimumStep {
        time_s: f64,
        dt_s: f64,
    },
    TimeAccuracy {
        time_s: f64,
        dt_s: f64,
        error: ImplicitFlowTimeError,
    },
    EnergyBalance {
        residual_j: f64,
        absolute_accepted_residual_j: f64,
    },
}
impl From<ImplicitFlowError> for AdaptiveImplicitFlowError {
    fn from(error: ImplicitFlowError) -> Self {
        Self::Solve(error)
    }
}
impl From<DepositionError> for AdaptiveImplicitFlowError {
    fn from(error: DepositionError) -> Self {
        Self::Solve(error.into())
    }
}
impl core::fmt::Display for AdaptiveImplicitFlowError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "adaptive implicit material flow: {self:?}")
    }
}
impl std::error::Error for AdaptiveImplicitFlowError {}

// Compensate physical changes and signed heat channels, not absolute formation
// energies. This small account has no source, stepping or acceptance authority.
#[derive(Clone, Copy, Default)]
struct Sum {
    sum: f64,
    correction: f64,
}
impl Sum {
    fn add(&mut self, value: f64) -> Result<(), DepositionError> {
        let value = finite(value)?;
        let next = finite(self.sum + value)?;
        let lost = if self.sum.abs() >= value.abs() {
            (self.sum - next) + value
        } else {
            (value - next) + self.sum
        };
        self.correction = finite(self.correction + lost)?;
        self.sum = next;
        Ok(())
    }
    fn value(self) -> Result<f64, DepositionError> {
        finite(self.sum + self.correction)
    }
}

fn time_error(
    initial: &DepositionSimulation<'_>,
    coarse: &DepositionSimulation<'_>,
    fine: &DepositionSimulation<'_>,
    correction: f64,
    controls: AdaptiveImplicitFlowControls,
    checkpoint: &mut impl FnMut() -> Result<(), DepositionError>,
) -> Result<ImplicitFlowTimeError, DepositionError> {
    let mut error = ImplicitFlowTimeError {
        maximum_temperature_difference_k: controls.temperature_absolute_tolerance_k.map(|_| 0.0),
        ..ImplicitFlowTimeError::default()
    };
    let component = |old: f64, a: f64, b: f64, absolute: f64| {
        let difference = finite(b - a)?.abs();
        let change = finite(a - old)?.abs().max(finite(b - old)?.abs());
        let allowance = finite(absolute + finite(controls.relative_tolerance * change)?)?;
        let normalized = (difference / allowance) * correction;
        if normalized.is_nan() || normalized < 0.0 || allowance <= 0.0 {
            return Err(DepositionError::Arithmetic);
        }
        Ok((difference, normalized))
    };
    for (index, (&id, cell)) in initial.cells().iter().enumerate() {
        if index % 256 == 0 {
            checkpoint()?;
        }
        let (difference, normalized) = component(
            cell.specific_h_j_kg,
            coarse.cells()[&id].specific_h_j_kg,
            fine.cells()[&id].specific_h_j_kg,
            controls.enthalpy_absolute_tolerance_j_kg,
        )?;
        error.maximum_enthalpy_difference_j_kg =
            error.maximum_enthalpy_difference_j_kg.max(difference);
        error.normalized = error.normalized.max(normalized);
        if let Some(tolerance) = controls.temperature_absolute_tolerance_k {
            let (difference, normalized) = component(
                initial.temperature_k(id)?,
                coarse.temperature_k(id)?,
                fine.temperature_k(id)?,
                tolerance,
            )?;
            error.maximum_temperature_difference_k = Some(
                error
                    .maximum_temperature_difference_k
                    .unwrap_or(0.0)
                    .max(difference),
            );
            error.normalized = error.normalized.max(normalized);
        }
    }
    checkpoint()?;
    Ok(error)
}

impl<F: Clone, R, T: FlowTransportSnapshot> MaterialFlowFeedback<'_, F, R, (), T>
where
    R: Clone + FnMut(usize, &DepositedCell, f64) -> Result<F, MaterialError>,
{
    /// Adapt the existing quasisteady thermoviscous BE solve in physical time.
    ///
    /// Each comparison is one full native step versus two native fine steps.
    /// Original material/rheology, momentum, dissipation, fluid conduction and
    /// stream enthalpy are resolved anew at every trial endpoint. Only the two
    /// fine states and their physical accounts are accepted. The same immutable
    /// blueprint is used throughout; no trial becomes a new material history.
    ///
    /// Mandatory h error is datum-invariant and detects latent-energy changes
    /// at fixed temperature. Optional T error queries the original calorics.
    /// The first-order fine correction is one for equal halves, or
    /// (a²+b²)/(2ab) for actual representable spans. This is a local estimate,
    /// not a certified global/spatial error bound. Nonlinear/source failures
    /// remain typed refusals rather than being relabeled as accuracy retries.
    ///
    /// Retained momentum, buoyancy and wrapped histories retain the one-step
    /// owner's exclusions. The whole positive interval must contain no native
    /// events; target-time events remain pending. The native accepted-substeps
    /// budget counts both fine steps. All speculative work shares ONE original
    /// max_solves budget. Failure publishes no thermal/flow prefix. Resume by
    /// passing the accepted thermal state to this same blueprint.
    pub fn advance_implicit_adaptive_to<'t>(
        &self,
        thermal: &DepositionSimulation<'t>,
        target_s: f64,
        controls: AdaptiveImplicitFlowControls,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<
        (DepositionSimulation<'t>, AdaptiveImplicitFlowReport<F, T>),
        AdaptiveImplicitFlowError,
    > {
        check_gate(self.gate, &mut checkpoint)?;
        let positive = |x: f64| x.is_finite() && x > 0.0;
        if !target_s.is_finite()
            || target_s <= thermal.time_s()
            || !positive(controls.minimum_dt_s)
            || !positive(controls.initial_dt_s)
            || !positive(controls.maximum_dt_s)
            || controls.minimum_dt_s > controls.initial_dt_s
            || controls.initial_dt_s > controls.maximum_dt_s
            || !positive(controls.enthalpy_absolute_tolerance_j_kg)
            || controls
                .temperature_absolute_tolerance_k
                .is_some_and(|x| !positive(x))
            || !controls.relative_tolerance.is_finite()
            || !(0.0..1.0).contains(&controls.relative_tolerance)
            || !controls.energy_tolerance_j.is_finite()
            || controls.energy_tolerance_j < 0.0
            || thermal.cells().len() > controls.step.max_unknowns
        {
            return Err(AdaptiveImplicitFlowError::Invalid(
                "invalid adaptive implicit-flow target, controls or cell budget",
            ));
        }
        // Refuse an interior event before any momentum work, even when a short
        // first proposal would not yet encounter it. No hidden event handling.
        BackwardEulerStep::new(thermal, target_s, || check_gate(self.gate, &mut checkpoint))?;
        let initial = vec![0.0; thermal.cells().len()];
        controls.step.admit(initial.len(), &initial)?;
        let mut mechanical = self
            .viscous_heating
            .clone()
            .ok_or(ImplicitFlowError::Invalid(
                "implicit material flow needs the existing stationary Stokes heating owner",
            ))?;
        mechanical.begin();
        let work = Cell::new(0);
        let mut next = thermal.clone();
        let mut suggested = controls.initial_dt_s;
        let (mut trials, mut accepted, mut rejected) = (0usize, 0usize, 0usize);
        let mut maximum_error: f64 = 0.0;
        // External heat, physical loss, stream enthalpy, other supply, |residual|.
        let mut heat = [Sum::default(); 5];
        let (last_step, last_error) = loop {
            check_gate(self.gate, &mut checkpoint)?;
            if trials >= controls.maximum_trials {
                return Err(DepositionError::Budget("adaptive implicit-flow trials").into());
            }
            if thermal
                .budget()
                .max_substeps
                .saturating_sub(mechanical.audit.substeps)
                < 2
            {
                return Err(
                    DepositionError::Budget("adaptive implicit-flow accepted intervals").into(),
                );
            }
            let start = next.time_s();
            let remaining = target_s - start;
            let proposed = suggested.min(remaining);
            let mut end = if proposed == remaining {
                target_s
            } else {
                finite(start + proposed)?.min(target_s)
            };
            if end - start > controls.maximum_dt_s {
                end = end.next_down();
            }
            let middle = start.midpoint(end);
            let (dt, a, b) = (end - start, middle - start, end - middle);
            if middle <= start || middle >= end || (dt < controls.minimum_dt_s && end != target_s) {
                return Err(AdaptiveImplicitFlowError::MinimumStep {
                    time_s: start,
                    dt_s: dt,
                });
            }
            let correction = ConductionTimeOrder::First
                .fine_error_correction(a, b)
                .ok_or(DepositionError::Arithmetic)?;
            trials += 1;
            let (coarse, _) = self.step_implicit_with_work(
                &next,
                end,
                controls.step,
                &initial,
                &work,
                &mut checkpoint,
            )?;
            let (half, first) = self.step_implicit_with_work(
                &next,
                middle,
                controls.step,
                &initial,
                &work,
                &mut checkpoint,
            )?;
            let (fine, last) = self.step_implicit_with_work(
                &half,
                end,
                controls.step,
                &initial,
                &work,
                &mut checkpoint,
            )?;
            let error = time_error(&next, &coarse, &fine, correction, controls, &mut || {
                check_gate(self.gate, &mut checkpoint)
            })?;
            if error.normalized > 1.0 {
                rejected += 1;
                if dt <= controls.minimum_dt_s {
                    return Err(AdaptiveImplicitFlowError::TimeAccuracy {
                        time_s: start,
                        dt_s: dt,
                        error,
                    });
                }
                let reduction = ConductionTimeOrder::First
                    .rejected_reduction_factor(error.normalized)
                    .ok_or(DepositionError::Arithmetic)?;
                suggested = (dt * reduction).max(controls.minimum_dt_s);
                if suggested >= dt {
                    return Err(AdaptiveImplicitFlowError::MinimumStep {
                        time_s: start,
                        dt_s: dt,
                    });
                }
                continue;
            }
            for (report, start, span, end) in [(&first, start, a, middle), (&last, middle, b, end)]
            {
                mechanical.accept(
                    report
                        .endpoint
                        .viscous_power
                        .as_ref()
                        .ok_or(ImplicitFlowError::Invalid(
                            "accepted implicit endpoint lacks its original physical power",
                        ))?,
                    start,
                    span,
                    end,
                )?;
                for (sum, value) in heat.iter_mut().zip([
                    report.thermal.external_heat_j,
                    report.deposited_heat_j,
                    report.net_stream_enthalpy_j,
                    report.other_thermal_input_j,
                    report.thermal.residual_j.abs(),
                ]) {
                    sum.add(value)?;
                }
            }
            let absolute = heat[4].value()?;
            if absolute > controls.energy_tolerance_j {
                return Err(AdaptiveImplicitFlowError::EnergyBalance {
                    residual_j: 0.0,
                    absolute_accepted_residual_j: absolute,
                });
            }
            accepted += 1;
            maximum_error = maximum_error.max(error.normalized);
            let growth = ConductionTimeOrder::First
                .accepted_growth_factor(error.normalized)
                .ok_or(DepositionError::Arithmetic)?;
            suggested = (dt * growth).clamp(controls.minimum_dt_s, controls.maximum_dt_s);
            next = fine;
            if next.time_s() == target_s {
                break (last, error);
            }
        };
        let mut storage = Sum::default();
        for (index, (&id, cell)) in thermal.cells().iter().enumerate() {
            if index % 256 == 0 {
                check_gate(self.gate, &mut checkpoint)?;
            }
            storage.add(finite(
                cell.mass_kg * finite(next.cells()[&id].specific_h_j_kg - cell.specific_h_j_kg)?,
            )?)?;
        }
        let stored = storage.value()?;
        let external = heat[0].value()?;
        let mut balance = storage;
        balance.add(-heat[0].sum)?;
        balance.add(-heat[0].correction)?;
        let residual = balance.value()?;
        let deposited = heat[1].value()?;
        let stream = heat[2].value()?;
        let other = heat[3].value()?;
        let absolute = heat[4].value()?;
        let mut combined = storage;
        for value in [
            heat[2].sum,
            heat[2].correction,
            heat[3].sum,
            heat[3].correction,
            mechanical.audit.supplied_wall_work_j,
            mechanical.audit.supplied_pressure_work_j,
            mechanical.audit.supplied_body_work_j,
        ] {
            combined.add(-value)?;
        }
        let combined = combined.value()?;
        let power_allowance =
            finite(mechanical.controls.power_tolerance_w * mechanical.audit.duration_s)?;
        if residual.abs() > controls.energy_tolerance_j
            || combined.abs() > finite(controls.energy_tolerance_j + power_allowance)?
            || (deposited - mechanical.audit.generated_heat_j).abs() > power_allowance
        {
            return Err(AdaptiveImplicitFlowError::EnergyBalance {
                residual_j: if residual.abs() > controls.energy_tolerance_j {
                    residual
                } else {
                    combined
                },
                absolute_accepted_residual_j: absolute,
            });
        }
        check_gate(self.gate, &mut checkpoint)?;
        Ok((
            next,
            AdaptiveImplicitFlowReport {
                thermal: EnergyAudit {
                    incoming_mass_kg: 0.0,
                    incoming_energy_j: 0.0,
                    external_heat_j: external,
                    stored_energy_change_j: stored,
                    residual_j: residual,
                    substeps: mechanical.audit.substeps,
                    events: 0,
                },
                mechanical: mechanical.audit,
                deposited_heat_j: deposited,
                net_stream_enthalpy_j: stream,
                other_thermal_input_j: other,
                combined_residual_j: combined,
                absolute_accepted_energy_residual_j: absolute,
                momentum_solves: work.get(),
                trials,
                accepted_intervals: accepted,
                accuracy_rejections: rejected,
                maximum_accepted_normalized_error: maximum_error,
                last_error,
                suggested_dt_s: suggested,
                last_step,
            },
        ))
    }
}
