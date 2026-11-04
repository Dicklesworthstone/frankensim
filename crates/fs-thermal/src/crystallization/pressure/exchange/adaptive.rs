//! Local time-error control around the original conservative contact solve.
//!
//! The two half steps are the physical trajectory. The full comparison step
//! spends work but supplies no heat, phase history or pressure work. This is
//! separate from each parcel's adaptive kinetics on a prescribed heat path.

use fs_math::dd::Dd;

use super::super::{
    CrystallizationError, CrystallizationState, check_cancelled, finite, invalid, nonnegative,
    parcel::IsobaricCrystallizationParcel, positive, product, sum,
};
use super::{
    IsobaricExchangeControls, IsobaricExchangeSideControls, IsobaricExchangeStep,
    exchange_heat_impl,
};

/// Explicit local accuracy and whole-call work/energy limits.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct IsobaricExchangeAdaptivity {
    pub initial_dt_s: f64,
    pub maximum_dt_s: f64,
    /// Bounds full comparison spans. Accepted half steps can be shorter; an
    /// exact final remainder below this floor can pass its first comparison.
    pub minimum_dt_s: f64,
    /// Includes rejected comparisons; at most three original contact solves
    /// per comparison, each with the caller's separate evaluation ceiling.
    pub maximum_trials: usize,
    pub temperature_tolerance_k: f64,
    pub specific_enthalpy_tolerance_j_kg: f64,
    pub crystalline_mass_fraction_tolerance: f64,
    pub exposure_tolerance: f64,
    /// Scales CHANGES from the accepted state, never an arbitrary caloric datum
    /// or absolute Kelvin. Zero requests purely absolute local tolerances.
    pub relative_tolerance: f64,
    /// Sum of noncancelling accepted-step defects and final account rounding.
    pub energy_tolerance_j: f64,
}

impl Default for IsobaricExchangeAdaptivity {
    fn default() -> Self {
        Self {
            initial_dt_s: 1.0,
            maximum_dt_s: 1.0,
            minimum_dt_s: 1e-9,
            maximum_trials: 10_000,
            temperature_tolerance_k: 1e-3,
            specific_enthalpy_tolerance_j_kg: 1.0,
            crystalline_mass_fraction_tolerance: 1e-6,
            exposure_tolerance: 1e-6,
            relative_tolerance: 1e-3,
            energy_tolerance_j: 1e-4,
        }
    }
}

impl IsobaricExchangeAdaptivity {
    pub fn validate(self) -> Result<(), CrystallizationError> {
        if !positive(self.minimum_dt_s)
            || !positive(self.maximum_dt_s)
            || !positive(self.initial_dt_s)
            || self.minimum_dt_s > self.initial_dt_s
            || self.initial_dt_s > self.maximum_dt_s
            || self.maximum_trials == 0
            || self.maximum_trials > usize::MAX / 3
            || !positive(self.temperature_tolerance_k)
            || !positive(self.specific_enthalpy_tolerance_j_kg)
            || !positive(self.crystalline_mass_fraction_tolerance)
            || self.crystalline_mass_fraction_tolerance > 1.0
            || !positive(self.exposure_tolerance)
            || !nonnegative(self.relative_tolerance)
            || self.relative_tolerance >= 1.0
            || !nonnegative(self.energy_tolerance_j)
        {
            return Err(invalid("invalid adaptive isobaric contact controls"));
        }
        Ok(())
    }
}

/// Original physical increments over the accepted trajectory, in joules.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct IsobaricExchangeBodyAccount {
    pub supplied_heat_j: f64,
    pub enthalpy_change_j: f64,
    pub internal_energy_change_j: f64,
    pub boundary_work_j: f64,
    pub energy_residual_j: f64,
}

/// Local coarse/fine discrepancies; these are not certified global bounds.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct IsobaricExchangeTimeError {
    pub temperature_k: f64,
    pub specific_enthalpy_j_kg: f64,
    pub crystalline_mass_fraction: f64,
    pub exposure: f64,
    pub normalized: f64,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AdaptiveIsobaricExchangeStep {
    pub initial_time_s: f64,
    pub final_time_s: f64,
    pub heat_left_to_right_j: f64,
    pub left: IsobaricExchangeBodyAccount,
    pub right: IsobaricExchangeBodyAccount,
    pub energy_residual_j: f64,
    pub absolute_energy_residual_j: f64,
    /// Sum of absolute contact-equation residuals on retained half steps.
    /// This is a constitutive solve error, separate from energy conservation.
    pub absolute_contact_residual_j: f64,
    pub accepted_substeps: usize,
    pub attempted_trials: usize,
    pub rejected_trials: usize,
    /// Includes successful coarse and accuracy-rejected comparison work.
    pub evaluated_candidates: usize,
    pub maximum_accepted_error: IsobaricExchangeTimeError,
    pub suggested_dt_s: f64,
}

/// Evolve a finite contact to an exact endpoint, retaining the original parcel
/// sources, remelting selections, masses, reference volumes and pressure loads.
///
/// The contact's supplied constant total conductance is held over this call.
/// The explicit heat bracket must remain valid for every candidate; it is not
/// extrapolated or silently clipped to material domains. Original source,
/// bracket, nonlinear and cancellation errors propagate without replacement.
/// Only excessive local time discrepancies request smaller comparison spans.
///
/// Temperature, enthalpy, crystal MASS fraction and unsaturated exposure all
/// participate in acceptance. Equal plateau temperatures alone cannot hide an
/// unresolved latent inventory. No extrapolated state is published. The caller
/// must still resolve source features through maximum_dt_s and study global
/// refinement; these local estimates certify neither phase-event timing nor
/// an arbitrary kinetic source's unseen narrow temperature features.
#[allow(clippy::too_many_arguments, clippy::type_complexity)]
pub fn advance_exchange_to<'ll, 'ls, 'rl, 'rs>(
    left: &IsobaricCrystallizationParcel<'ll, 'ls>,
    right: &IsobaricCrystallizationParcel<'rl, 'rs>,
    conductance_w_k: f64,
    target_s: f64,
    left_controls: IsobaricExchangeSideControls,
    right_controls: IsobaricExchangeSideControls,
    exchange: IsobaricExchangeControls,
    adaptivity: IsobaricExchangeAdaptivity,
    cancelled: impl FnMut() -> bool,
) -> Result<
    (
        IsobaricCrystallizationParcel<'ll, 'ls>,
        IsobaricCrystallizationParcel<'rl, 'rs>,
        AdaptiveIsobaricExchangeStep,
    ),
    CrystallizationError,
> {
    advance_exchange_to_impl(
        left,
        right,
        |_, _| Ok(conductance_w_k),
        target_s,
        left_controls,
        right_controls,
        exchange,
        adaptivity,
        Some(conductance_w_k),
        cancelled,
    )
}

/// Adaptive contact with an original conductance query at both proposed
/// endpoint states. The same source is used by full and retained half steps.
///
/// The query may read temperature, kinetic mass fraction, volume and original
/// source identities from each side. It must be deterministic and immutable in
/// its physical meaning: comparisons and rejected candidates must not advance
/// source history. No conductivity mixture rule or contact geometry is inferred.
/// Each evaluated value must be finite and nonnegative. Source errors propagate
/// unchanged, and the same time-error/work/energy semantics as
/// [`advance_exchange_to`] apply to the complete coupled nonlinear evolution.
#[allow(clippy::too_many_arguments, clippy::type_complexity)]
pub fn advance_exchange_to_with_conductance<'ll, 'ls, 'rl, 'rs>(
    left: &IsobaricCrystallizationParcel<'ll, 'ls>,
    right: &IsobaricCrystallizationParcel<'rl, 'rs>,
    conductance: impl FnMut(
        &IsobaricCrystallizationParcel<'ll, 'ls>,
        &IsobaricCrystallizationParcel<'rl, 'rs>,
    ) -> Result<f64, CrystallizationError>,
    target_s: f64,
    left_controls: IsobaricExchangeSideControls,
    right_controls: IsobaricExchangeSideControls,
    exchange: IsobaricExchangeControls,
    adaptivity: IsobaricExchangeAdaptivity,
    cancelled: impl FnMut() -> bool,
) -> Result<
    (
        IsobaricCrystallizationParcel<'ll, 'ls>,
        IsobaricCrystallizationParcel<'rl, 'rs>,
        AdaptiveIsobaricExchangeStep,
    ),
    CrystallizationError,
> {
    advance_exchange_to_impl(
        left,
        right,
        conductance,
        target_s,
        left_controls,
        right_controls,
        exchange,
        adaptivity,
        None,
        cancelled,
    )
}

#[allow(clippy::too_many_arguments, clippy::type_complexity)]
fn advance_exchange_to_impl<'ll, 'ls, 'rl, 'rs>(
    left: &IsobaricCrystallizationParcel<'ll, 'ls>,
    right: &IsobaricCrystallizationParcel<'rl, 'rs>,
    mut conductance: impl FnMut(
        &IsobaricCrystallizationParcel<'ll, 'ls>,
        &IsobaricCrystallizationParcel<'rl, 'rs>,
    ) -> Result<f64, CrystallizationError>,
    target_s: f64,
    left_controls: IsobaricExchangeSideControls,
    right_controls: IsobaricExchangeSideControls,
    exchange: IsobaricExchangeControls,
    adaptivity: IsobaricExchangeAdaptivity,
    fixed_conductance_w_k: Option<f64>,
    mut cancelled: impl FnMut() -> bool,
) -> Result<
    (
        IsobaricCrystallizationParcel<'ll, 'ls>,
        IsobaricCrystallizationParcel<'rl, 'rs>,
        AdaptiveIsobaricExchangeStep,
    ),
    CrystallizationError,
> {
    check_cancelled(&mut cancelled)?;
    adaptivity.validate()?;
    exchange.validate()?;
    left_controls.validate()?;
    right_controls.validate()?;
    if fixed_conductance_w_k.is_some_and(|value| !nonnegative(value))
        || !target_s.is_finite()
        || target_s <= left.time_s()
        || left.time_s() != right.time_s()
        || left_controls.caloric.caloric.temperature_tolerance_k
            > adaptivity.temperature_tolerance_k
        || right_controls.caloric.caloric.temperature_tolerance_k
            > adaptivity.temperature_tolerance_k
        || exchange.maximum_evaluations > usize::MAX / (3 * adaptivity.maximum_trials)
    {
        return Err(invalid(
            "invalid adaptive contact time, tolerance or work bound",
        ));
    }
    let mut report = AdaptiveIsobaricExchangeStep {
        initial_time_s: left.time_s(),
        final_time_s: target_s,
        heat_left_to_right_j: 0.0,
        left: IsobaricExchangeBodyAccount::default(),
        right: IsobaricExchangeBodyAccount::default(),
        energy_residual_j: 0.0,
        absolute_energy_residual_j: 0.0,
        absolute_contact_residual_j: 0.0,
        accepted_substeps: 0,
        attempted_trials: 0,
        rejected_trials: 0,
        evaluated_candidates: 0,
        maximum_accepted_error: IsobaricExchangeTimeError::default(),
        suggested_dt_s: adaptivity.initial_dt_s,
    };
    let mut current_left = left.clone();
    let mut current_right = right.clone();
    let mut totals = Totals::default();
    let mut requested_dt = adaptivity.initial_dt_s;
    while current_left.time_s() < target_s {
        check_cancelled(&mut cancelled)?;
        if report.attempted_trials == adaptivity.maximum_trials {
            return Err(CrystallizationError::NonConverged);
        }
        let time = current_left.time_s();
        let remaining = finite(target_s - time)?;
        let endpoint = if requested_dt >= remaining {
            target_s
        } else {
            finite(time + requested_dt)?
        };
        let dt = finite(endpoint - time)?;
        let midpoint = finite(time + 0.5 * dt)?;
        if midpoint <= time || midpoint >= endpoint {
            return Err(CrystallizationError::NonConverged);
        }
        report.attempted_trials += 1;
        let (coarse_left, coarse_right, coarse) = exchange_heat_impl(
            &current_left,
            &current_right,
            &mut conductance,
            dt,
            left_controls,
            right_controls,
            exchange,
            fixed_conductance_w_k,
            &mut cancelled,
        )?;
        let (middle_left, middle_right, first) = exchange_heat_impl(
            &current_left,
            &current_right,
            &mut conductance,
            midpoint - time,
            left_controls,
            right_controls,
            exchange,
            fixed_conductance_w_k,
            &mut cancelled,
        )?;
        let (fine_left, fine_right, second) = exchange_heat_impl(
            &middle_left,
            &middle_right,
            &mut conductance,
            endpoint - midpoint,
            left_controls,
            right_controls,
            exchange,
            fixed_conductance_w_k,
            &mut cancelled,
        )?;
        if coarse_left.time_s() != endpoint || fine_left.time_s() != endpoint {
            return Err(invalid(
                "adaptive contact could not represent the requested endpoint",
            ));
        }
        report.evaluated_candidates +=
            coarse.evaluated_candidates + first.evaluated_candidates + second.evaluated_candidates;
        // Unequal representable halves need the first-order correction for
        // (a²+b²)/(2ab). Equal halves give exactly one, without extrapolation.
        let a = first.duration_s / dt;
        let b = second.duration_s / dt;
        let correction = finite((a * a + b * b) / (2.0 * a * b))?;
        let mut error = IsobaricExchangeTimeError::default();
        compare(
            &mut error,
            &current_left,
            coarse_left.state(),
            fine_left.state(),
            adaptivity,
            correction,
        )?;
        compare(
            &mut error,
            &current_right,
            coarse_right.state(),
            fine_right.state(),
            adaptivity,
            correction,
        )?;
        if error.normalized > 1.0 {
            report.rejected_trials += 1;
            let half = 0.5 * dt;
            if half < adaptivity.minimum_dt_s || !positive(half) {
                return Err(CrystallizationError::NonConverged);
            }
            requested_dt = half;
            continue;
        }
        totals.include(first)?;
        totals.include(second)?;
        if totals.defects.to_f64() > adaptivity.energy_tolerance_j {
            return Err(invalid(
                "adaptive contact cumulative energy tolerance exceeded",
            ));
        }
        report.accepted_substeps += 2;
        report.maximum_accepted_error.include(error);
        current_left = fine_left;
        current_right = fine_right;
        requested_dt = if error.normalized <= 0.125 {
            if dt >= 0.5 * adaptivity.maximum_dt_s {
                adaptivity.maximum_dt_s
            } else {
                2.0 * dt
            }
        } else {
            dt.min(adaptivity.maximum_dt_s)
        };
        requested_dt = requested_dt.max(adaptivity.minimum_dt_s);
        report.suggested_dt_s = requested_dt;
    }
    totals.finish(left, right, &current_left, &current_right, &mut report)?;
    if report.absolute_energy_residual_j > adaptivity.energy_tolerance_j {
        return Err(invalid("adaptive contact final energy tolerance exceeded"));
    }
    check_cancelled(&mut cancelled)?;
    Ok((current_left, current_right, report))
}

fn compare(
    error: &mut IsobaricExchangeTimeError,
    initial: &IsobaricCrystallizationParcel<'_, '_>,
    coarse: CrystallizationState,
    fine: CrystallizationState,
    controls: IsobaricExchangeAdaptivity,
    correction: f64,
) -> Result<(), CrystallizationError> {
    let old = initial.state();
    let quantities = [
        (
            old.temperature_k(),
            coarse.temperature_k(),
            fine.temperature_k(),
            controls.temperature_tolerance_k,
        ),
        (
            old.specific_enthalpy_j_kg(),
            coarse.specific_enthalpy_j_kg(),
            fine.specific_enthalpy_j_kg(),
            controls.specific_enthalpy_tolerance_j_kg,
        ),
        (
            old.crystalline_mass_fraction(),
            coarse.crystalline_mass_fraction(),
            fine.crystalline_mass_fraction(),
            controls.crystalline_mass_fraction_tolerance,
        ),
        (
            old.exposure(),
            coarse.exposure(),
            fine.exposure(),
            controls.exposure_tolerance,
        ),
    ];
    let mut differences = [0.0; 4];
    for (i, (old, coarse_value, fine_value, absolute)) in quantities.into_iter().enumerate() {
        if i == 3
            && coarse.crystalline_mass_fraction()
                == initial.law().config().maximum_crystalline_mass_fraction
            && fine.crystalline_mass_fraction()
                == initial.law().config().maximum_crystalline_mass_fraction
        {
            continue;
        }
        let difference = finite(fine_value - coarse_value)?.abs();
        let change = finite(coarse_value - old)?
            .abs()
            .max(finite(fine_value - old)?.abs());
        let allowance = finite(absolute + controls.relative_tolerance * change)?;
        let ratio = (difference / allowance) * correction;
        if ratio.is_nan() {
            return Err(invalid("invalid adaptive contact discrepancy"));
        }
        differences[i] = difference;
        error.normalized = error.normalized.max(ratio);
    }
    error.temperature_k = error.temperature_k.max(differences[0]);
    error.specific_enthalpy_j_kg = error.specific_enthalpy_j_kg.max(differences[1]);
    error.crystalline_mass_fraction = error.crystalline_mass_fraction.max(differences[2]);
    error.exposure = error.exposure.max(differences[3]);
    Ok(())
}

impl IsobaricExchangeTimeError {
    fn include(&mut self, other: Self) {
        self.temperature_k = self.temperature_k.max(other.temperature_k);
        self.specific_enthalpy_j_kg = self
            .specific_enthalpy_j_kg
            .max(other.specific_enthalpy_j_kg);
        self.crystalline_mass_fraction = self
            .crystalline_mass_fraction
            .max(other.crystalline_mass_fraction);
        self.exposure = self.exposure.max(other.exposure);
        self.normalized = self.normalized.max(other.normalized);
    }
}

#[derive(Default)]
struct Totals {
    heat: Dd,
    enthalpy: [Dd; 2],
    internal: [Dd; 2],
    work: [Dd; 2],
    defects: Dd,
    contact: Dd,
}
impl Totals {
    fn include(&mut self, step: IsobaricExchangeStep) -> Result<(), CrystallizationError> {
        self.heat = self.heat + Dd::from_f64(step.heat_left_to_right_j);
        self.defects = self.defects + Dd::from_f64(step.absolute_energy_residual_j);
        self.contact = self.contact + Dd::from_f64(step.residual_j.abs());
        for (i, body) in [step.left, step.right].into_iter().enumerate() {
            self.enthalpy[i] = self.enthalpy[i] + Dd::from_f64(body.enthalpy_change_j);
            self.internal[i] = self.internal[i] + Dd::from_f64(body.internal_energy_change_j);
            self.work[i] = self.work[i] + Dd::from_f64(body.boundary_work_j);
        }
        for value in [self.heat, self.defects, self.contact]
            .into_iter()
            .chain(self.enthalpy)
            .chain(self.internal)
            .chain(self.work)
        {
            finite(value.hi)?;
            finite(value.lo)?;
        }
        Ok(())
    }

    fn finish(
        &mut self,
        left: &IsobaricCrystallizationParcel<'_, '_>,
        right: &IsobaricCrystallizationParcel<'_, '_>,
        final_left: &IsobaricCrystallizationParcel<'_, '_>,
        final_right: &IsobaricCrystallizationParcel<'_, '_>,
        report: &mut AdaptiveIsobaricExchangeStep,
    ) -> Result<(), CrystallizationError> {
        report.heat_left_to_right_j = self.heat.to_f64();
        report.absolute_contact_residual_j = self.contact.to_f64();
        let mut bodies = [IsobaricExchangeBodyAccount::default(); 2];
        for (i, (initial, final_body, sign)) in
            [(left, final_left, -1.0), (right, final_right, 1.0)]
                .into_iter()
                .enumerate()
        {
            let body = &mut bodies[i];
            body.supplied_heat_j = sign * report.heat_left_to_right_j;
            body.enthalpy_change_j = product(
                initial.mass_kg(),
                finite(
                    final_body.state().specific_enthalpy_j_kg()
                        - initial.state().specific_enthalpy_j_kg(),
                )?,
            )?;
            body.internal_energy_change_j = self.internal[i].to_f64();
            body.boundary_work_j = product(
                initial.pressure_pa(),
                finite(final_body.volume_m3() - initial.volume_m3())?,
            )?;
            body.energy_residual_j = sum(&[
                body.internal_energy_change_j,
                body.boundary_work_j,
                -body.supplied_heat_j,
            ])?;
            // Charge final endpoint/telescope and public-total rounding instead
            // of hiding it behind cancellation of two large formation datums.
            self.defects = self.defects
                + (Dd::from_f64(body.enthalpy_change_j) - self.enthalpy[i]).abs()
                + (Dd::from_f64(body.boundary_work_j) - self.work[i]).abs()
                + Dd::from_f64(self.internal[i].lo.abs())
                + Dd::from_f64(self.heat.lo.abs());
        }
        report.left = bodies[0];
        report.right = bodies[1];
        report.energy_residual_j = sum(&[
            bodies[0].internal_energy_change_j,
            bodies[1].internal_energy_change_j,
            bodies[0].boundary_work_j,
            bodies[1].boundary_work_j,
        ])?;
        report.absolute_energy_residual_j = finite(self.defects.to_f64())?;
        // The published whole-system residual must itself fit the report even
        // if final scalar summation incurred another rounding error.
        report.absolute_energy_residual_j = report
            .absolute_energy_residual_j
            .max(report.energy_residual_j.abs())
            .max(sum(&[
                bodies[0].energy_residual_j.abs(),
                bodies[1].energy_residual_j.abs(),
            ])?);
        Ok(())
    }
}

#[cfg(test)]
#[path = "adaptive/tests.rs"]
mod tests;
