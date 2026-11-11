//! Local time refinement for simultaneous reaction, heat exchange and volume work.
//!
//! One full interval and two half intervals use the same original thermal
//! interval kernel and the same prescribed logarithmic volume path. Only the
//! fine gas AND finite-body endpoints may enter the accepted trajectory. Raw
//! discrepancies control gas temperature, every original species, gas work,
//! exchanged heat and actual finite-body enthalpy; body temperature alone would
//! miss latent storage. These comparisons do not certify a global time error.

use super::super::adaptive::{
    AdaptiveMaterialCompressionError, AdaptiveMaterialCompressionSolve, compare, error_root,
    scaled_difference,
};
use super::super::{
    CompressionInterval, CompressionTotals, MAX_COMPRESSION_SUBSTEPS, MaterialCompressionError,
    accumulate, advance_interval_with, positive,
};
use super::{
    ClosedGasParcel, MaterialThermalCompressionBodyStep, MaterialThermalCompressionReactor,
    MaterialThermalCompressionStep, ThermalInterval,
};
use crate::chemistry::heat::{HeatBath, material::FiniteHeatBody};
use crate::material::isochoric::MaterialIsochoricReactor;
use fs_material::thermochemical_gas::energy::entropy::GasEntropySolve;
use fs_math::dd::Dd;
use fs_qty::{Time, Volume};

pub use super::super::adaptive::MAX_COMPRESSION_EVALUATIONS;

/// Native calorics and explicit local full/fine comparison allowances.
///
/// The prepared reactor retains its original kinetics controls. Positive
/// absolute allowances are required even when a particular channel is inactive.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AdaptiveMaterialThermalCompressionSolve {
    pub volume: GasEntropySolve,
    pub absolute_temperature_tolerance_k: f64,
    /// Allowance for EACH original species, including initially absent products.
    pub absolute_amount_tolerance_mol: f64,
    pub absolute_work_tolerance_j: f64,
    pub absolute_heat_tolerance_j: f64,
    /// Actual mass times source-specific enthalpy discrepancy [J]. This controls
    /// latent phase inventory even when both body temperatures are identical.
    /// The corresponding comparison is zero for a held-temperature bath.
    pub absolute_body_enthalpy_tolerance_j: f64,
    /// In [0,1). Scales interval CHANGES, never absolute formation/enthalpy
    /// datums, carrier inventories or work/heat accumulated before the trial.
    pub relative_tolerance: f64,
    pub initial_step_s: f64,
    /// An exact final or prescribed-forcing-knot remainder may be smaller; a
    /// failed comparison there cannot reduce again below this minimum.
    pub minimum_step_s: f64,
    pub maximum_step_s: f64,
    /// Accepted full/fine comparisons, 1..=4096; each retains two intervals.
    /// A positive-duration fixed-volume call consumes one accepted interval.
    pub maximum_intervals: usize,
    /// ALL original interval probes, including rejected/coarse work, 3..=12288.
    /// Fixed volume needs one probe after admission. Native chemistry retains
    /// its own per-call work bound; total chemical work is separately reported.
    pub maximum_evaluations: usize,
}

/// Accepted original physical accounts and actual adaptive work.
///
/// The default physical payload describes a held bath. The finite-body alias
/// retains the original body endpoint and its distinct energy/work accounts.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AdaptiveMaterialThermalCompressionStep<P = MaterialThermalCompressionStep> {
    pub physical: P,
    pub accepted_intervals: usize,
    pub rejected_intervals: usize,
    pub evaluations: usize,
    /// Gas T, maximum original-species amount, gas work, exchanged heat, body H.
    /// Raw full/fine differences divided by their declared allowances, with no
    /// Richardson divisor or extrapolation. Accepted ratios are at most one.
    pub maximum_accepted_error_ratios: [f64; 5],
    pub minimum_accepted_step_s: f64,
    pub maximum_accepted_step_s: f64,
    pub suggested_step_s: f64,
    /// Chemical work across ALL successful probes, including coarse/rejected
    /// trajectories; physical counts retain accepted fine trajectories only.
    pub total_chemical_accepted_steps: usize,
    pub total_chemical_rejected_steps: usize,
}

pub type AdaptiveMaterialThermalCompressionBodyStep<'m> =
    AdaptiveMaterialThermalCompressionStep<MaterialThermalCompressionBodyStep<'m>>;

/// Any refusal preserves the complete caller parcel and supplied original body.
#[derive(Clone, Debug, PartialEq)]
pub enum AdaptiveMaterialThermalCompressionError {
    Invalid(&'static str),
    /// Original source, kinetics, heat/body or volume refusal. Physical failures
    /// are terminal; only an actual comparison failure requests smaller time.
    Physical(MaterialCompressionError),
    Accuracy {
        elapsed_time_s: f64,
        step_s: f64,
        error_ratios: [f64; 5],
    },
    WorkLimit {
        elapsed_time_s: f64,
        accepted_intervals: usize,
        rejected_intervals: usize,
        evaluations: usize,
    },
    TimeResolution {
        elapsed_time_s: f64,
        target_time_s: f64,
    },
    Cancelled,
}

type Error = AdaptiveMaterialThermalCompressionError;

impl core::fmt::Display for AdaptiveMaterialThermalCompressionError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "adaptive reacting thermal compression: {self:?}")
    }
}
impl std::error::Error for AdaptiveMaterialThermalCompressionError {}
impl From<MaterialCompressionError> for AdaptiveMaterialThermalCompressionError {
    fn from(error: MaterialCompressionError) -> Self {
        match error {
            MaterialCompressionError::Cancelled => Self::Cancelled,
            other => Self::Physical(other),
        }
    }
}
// A kinetics refusal reaches this type as the compression owner would report
// it: cancellation stays cancellation, anything else is a physical refusal.
impl From<crate::material::MaterialKineticsError> for AdaptiveMaterialThermalCompressionError {
    fn from(error: crate::material::MaterialKineticsError) -> Self {
        MaterialCompressionError::from(error).into()
    }
}
// The original comparison helpers currently return only Invalid or Cancelled.
// Preserve every native error shape if their internal scope later expands.
impl From<AdaptiveMaterialCompressionError> for AdaptiveMaterialThermalCompressionError {
    fn from(error: AdaptiveMaterialCompressionError) -> Self {
        match error {
            AdaptiveMaterialCompressionError::Invalid(message) => Self::Invalid(message),
            AdaptiveMaterialCompressionError::Physical(error) => Self::from(error),
            AdaptiveMaterialCompressionError::Cancelled => Self::Cancelled,
            AdaptiveMaterialCompressionError::Accuracy {
                elapsed_time_s,
                step_s,
                error_ratios,
            } => Self::Accuracy {
                elapsed_time_s,
                step_s,
                error_ratios: [error_ratios[0], error_ratios[1], error_ratios[2], 0.0, 0.0],
            },
            AdaptiveMaterialCompressionError::WorkLimit {
                elapsed_time_s,
                accepted_intervals,
                rejected_intervals,
                evaluations,
            } => Self::WorkLimit {
                elapsed_time_s,
                accepted_intervals,
                rejected_intervals,
                evaluations,
            },
            AdaptiveMaterialCompressionError::TimeResolution {
                elapsed_time_s,
                target_time_s,
            } => Self::TimeResolution {
                elapsed_time_s,
                target_time_s,
            },
        }
    }
}

impl AdaptiveMaterialThermalCompressionSolve {
    fn validate(self) -> Result<(), Error> {
        if ![
            self.absolute_temperature_tolerance_k,
            self.absolute_amount_tolerance_mol,
            self.absolute_work_tolerance_j,
            self.absolute_heat_tolerance_j,
            self.absolute_body_enthalpy_tolerance_j,
            self.initial_step_s,
            self.minimum_step_s,
            self.maximum_step_s,
        ]
        .into_iter()
        .all(positive)
            || !self.relative_tolerance.is_finite()
            || !(0.0..1.0).contains(&self.relative_tolerance)
            || self.minimum_step_s > self.initial_step_s
            || self.initial_step_s > self.maximum_step_s
            || !(1..=MAX_COMPRESSION_SUBSTEPS).contains(&self.maximum_intervals)
            || !(3..=MAX_COMPRESSION_EVALUATIONS).contains(&self.maximum_evaluations)
        {
            return Err(Error::Invalid(
                "positive channel tolerances, 0<=relative<1, ordered positive time bounds and bounded work required",
            ));
        }
        Ok(())
    }
}

impl MaterialThermalCompressionReactor {
    /// Adapt the compression/reaction/heat composition against a held bath.
    ///
    /// Native chemistry still controls its own local error and source/entropy
    /// gates. Outer tolerances below that accuracy may refuse at the minimum
    /// step. Fixed volume and zero duration retain one original reactor call.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_with_heat_bath_adaptive(
        &mut self,
        parcel: &mut ClosedGasParcel,
        final_volume: Volume,
        duration: Time,
        bath: HeatBath,
        controls: AdaptiveMaterialThermalCompressionSolve,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<AdaptiveMaterialThermalCompressionStep, Error> {
        let (candidate, step) = advance(
            self,
            parcel,
            final_volume,
            duration,
            bath,
            controls,
            &mut cancelled,
        )?;
        poll(&mut cancelled)?;
        *parcel = candidate;
        Ok(step)
    }

    /// Adapt both gas and finite-body endpoints, including original latent H.
    ///
    /// Gas heat cancels the body's original enthalpy transfer once. Body
    /// pressure-volume work, gas work and native residuals retain their separate
    /// accounts. A rejected probe or refused stroke publishes neither endpoint.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_with_finite_body_adaptive<'m>(
        &mut self,
        parcel: &mut ClosedGasParcel,
        final_volume: Volume,
        duration: Time,
        body: FiniteHeatBody<'m>,
        controls: AdaptiveMaterialThermalCompressionSolve,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<AdaptiveMaterialThermalCompressionBodyStep<'m>, Error> {
        let (candidate, step) = advance(
            self,
            parcel,
            final_volume,
            duration,
            body,
            controls,
            &mut cancelled,
        )?;
        poll(&mut cancelled)?;
        *parcel = candidate;
        Ok(step)
    }
}

// The two thermal partners differ only in their existing native reaction call
// and physical report. The controller and volume/chemistry kernel are shared.
pub(super) trait Partner: Copy {
    type Physical: Copy;

    fn react<F: FnMut() -> bool>(
        self,
        reactor: &mut MaterialIsochoricReactor,
        parcel: &mut ClosedGasParcel,
        duration: Time,
        clock: IntervalTime,
        cancelled: &mut F,
    ) -> Result<(Self, ThermalInterval, [f64; 5]), MaterialCompressionError>;

    /// Constant partners keep the original one-call fixed-volume route.
    fn is_time_constant(self) -> bool {
        true
    }

    /// The first prescribed forcing knot strictly after the accepted time.
    /// Native constant partners have no interior time discontinuities.
    fn next_boundary(self, _after_s: f64, duration_s: f64) -> f64 {
        duration_s
    }

    /// Original specific h and invariant mass, without an absolute H product.
    fn body_inventory(self) -> Option<[f64; 2]> {
        None
    }

    fn finish(self, gas: MaterialThermalCompressionStep, body: [f64; 5]) -> Self::Physical;
}

/// Actual physical endpoints of one original interval, also used for forcing.
#[derive(Clone, Copy)]
pub(super) struct IntervalTime {
    pub(super) start_s: f64,
    pub(super) end_s: f64,
}

impl IntervalTime {
    pub(super) fn midpoint_s(self) -> f64 {
        self.start_s + 0.5 * (self.end_s - self.start_s)
    }
}

impl Partner for HeatBath {
    type Physical = MaterialThermalCompressionStep;

    fn react<F: FnMut() -> bool>(
        self,
        reactor: &mut MaterialIsochoricReactor,
        parcel: &mut ClosedGasParcel,
        duration: Time,
        _clock: IntervalTime,
        cancelled: &mut F,
    ) -> Result<(Self, ThermalInterval, [f64; 5]), MaterialCompressionError> {
        let step = reactor.advance_with_heat_bath(parcel, duration, self, cancelled)?;
        let entropy = -step.chemistry.balance().entropy_transfer_j_k;
        Ok((self, ThermalInterval::from_step(&step, entropy), [0.0; 5]))
    }

    fn finish(self, gas: MaterialThermalCompressionStep, _: [f64; 5]) -> Self::Physical {
        gas
    }
}

impl<'m> Partner for FiniteHeatBody<'m> {
    type Physical = MaterialThermalCompressionBodyStep<'m>;

    fn react<F: FnMut() -> bool>(
        self,
        reactor: &mut MaterialIsochoricReactor,
        parcel: &mut ClosedGasParcel,
        duration: Time,
        _clock: IntervalTime,
        cancelled: &mut F,
    ) -> Result<(Self, ThermalInterval, [f64; 5]), MaterialCompressionError> {
        let step = reactor.advance_with_finite_body(parcel, duration, self, cancelled)?;
        Ok((
            step.body,
            ThermalInterval::from_step(&step.gas, step.body_entropy_change_j_k),
            [
                step.body_report.enthalpy_change_j,
                step.body_report.internal_energy_change_j,
                step.body_report.boundary_work_j,
                step.combined_heat_residual_j,
                step.combined_energy_residual_j,
            ],
        ))
    }

    fn body_inventory(self) -> Option<[f64; 2]> {
        Some([self.state().specific_enthalpy_j_kg(), self.mass().value()])
    }

    fn finish(self, gas: MaterialThermalCompressionStep, body: [f64; 5]) -> Self::Physical {
        MaterialThermalCompressionBodyStep {
            gas,
            body: self,
            body_enthalpy_change_j: body[0],
            body_internal_energy_change_j: body[1],
            body_boundary_work_j: body[2],
            body_entropy_change_j_k: gas.thermal_partner_entropy_change_j_k,
            combined_heat_residual_j: body[3],
            combined_energy_residual_j: body[4],
        }
    }
}

/// Accepted physical accounts. Rejected/coarse reductions never enter this sum.
#[derive(Clone, Copy, Default)]
pub(super) struct Totals {
    gas: CompressionTotals,
    // Heat, counterpart entropy, native gas-plus-counterpart production.
    thermal: [(f64, f64); 3],
    // Original body H, U, outward work and both coupled first-law residuals.
    body: [(f64, f64); 5],
}

impl Totals {
    pub(super) fn interval<B: Partner, F: FnMut() -> bool>(
        &mut self,
        reactor: &mut MaterialIsochoricReactor,
        candidate: &mut ClosedGasParcel,
        partner: &mut B,
        interval: CompressionInterval,
        clock: IntervalTime,
        volume: GasEntropySolve,
        cancelled: &mut F,
    ) -> Result<(), MaterialCompressionError> {
        let Self { gas, thermal, body } = self;
        let entropy = advance_interval_with(
            candidate,
            interval,
            volume,
            gas,
            cancelled,
            |candidate, duration, cancelled| {
                let (next, step, body_change) =
                    partner.react(reactor, candidate, duration, clock, cancelled)?;
                for (total, value) in thermal.iter_mut().zip([
                    step.heat_input_j,
                    step.partner_entropy_change_j_k,
                    step.entropy_production_j_k,
                ]) {
                    accumulate(total, value);
                }
                for (total, value) in body.iter_mut().zip(body_change) {
                    accumulate(total, value);
                }
                *partner = next;
                Ok(step.chemical)
            },
        )?;
        for change in entropy {
            accumulate(&mut thermal[2], change);
        }
        Ok(())
    }

    fn merge(&mut self, other: Self) -> Result<(), Error> {
        self.gas.merge(other.gas)?;
        for (into, from) in self
            .thermal
            .iter_mut()
            .zip(other.thermal)
            .chain(self.body.iter_mut().zip(other.body))
        {
            accumulate(into, from.0);
            accumulate(into, -from.1);
        }
        Ok(())
    }

    pub(super) fn finish<B: Partner>(
        self,
        partner: B,
        initial_volume_m3: f64,
        final_volume_m3: f64,
        duration: Time,
    ) -> Result<B::Physical, MaterialCompressionError> {
        let gas = self
            .gas
            .finish(initial_volume_m3, final_volume_m3, duration)?;
        if self
            .thermal
            .iter()
            .chain(&self.body)
            .any(|value| !value.0.is_finite())
        {
            return Err(MaterialCompressionError::Invalid(
                "accumulated thermal/body work, heat or entropy is unrepresentable",
            ));
        }
        Ok(partner.finish(
            MaterialThermalCompressionStep {
                elapsed_time: duration,
                initial_volume_m3,
                final_volume_m3,
                boundary_work_on_gas_j: gas.boundary_work_on_gas_j,
                heat_input_j: self.thermal[0].0,
                internal_energy_change_j: gas.internal_energy_change_j,
                energy_residual_j: gas.energy_residual_j,
                gas_entropy_change_j_k: gas.entropy_change_j_k,
                thermal_partner_entropy_change_j_k: self.thermal[1].0,
                entropy_production_j_k: self.thermal[2].0,
                substeps: gas.substeps,
                chemical_accepted_steps: gas.chemical_accepted_steps,
                chemical_rejected_steps: gas.chemical_rejected_steps,
            },
            self.body.map(|value| value.0),
        ))
    }
}

#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
pub(super) fn advance<B: Partner, F: FnMut() -> bool>(
    reactor: &mut MaterialThermalCompressionReactor,
    parcel: &ClosedGasParcel,
    final_volume: Volume,
    duration: Time,
    mut partner: B,
    controls: AdaptiveMaterialThermalCompressionSolve,
    cancelled: &mut F,
) -> Result<
    (
        ClosedGasParcel,
        AdaptiveMaterialThermalCompressionStep<B::Physical>,
    ),
    Error,
> {
    poll(cancelled)?;
    controls.validate()?;
    let initial_volume = parcel.volume_m3();
    let final_volume_m3 = final_volume.value();
    let duration_s = duration.value();
    if !positive(final_volume_m3) || !duration_s.is_finite() || duration_s < 0.0 {
        return Err(Error::Invalid(
            "positive finite volume and nonnegative finite duration required",
        ));
    }
    let moving = final_volume_m3 != initial_volume;
    let ratio = final_volume_m3 / initial_volume;
    if !positive(ratio) || moving && !positive(duration_s) {
        return Err(Error::Invalid(
            "moving strokes require positive time and a representable volume ratio",
        ));
    }
    let mut candidate = parcel.clone();
    candidate
        .set_volume_isentropic(initial_volume, controls.volume, &mut *cancelled)
        .map_err(MaterialCompressionError::from)?;
    let mut totals = Totals::default();
    if !moving && partner.is_time_constant() {
        totals.interval(
            &mut reactor.reactor,
            &mut candidate,
            &mut partner,
            CompressionInterval {
                duration_s,
                midpoint_volume_m3: initial_volume,
                final_volume_m3,
            },
            IntervalTime {
                start_s: 0.0,
                end_s: duration_s,
            },
            controls.volume,
            cancelled,
        )?;
        let physical = totals.finish(partner, initial_volume, final_volume_m3, duration)?;
        return Ok((
            candidate,
            AdaptiveMaterialThermalCompressionStep {
                physical,
                accepted_intervals: usize::from(duration_s > 0.0),
                rejected_intervals: 0,
                evaluations: 1,
                maximum_accepted_error_ratios: [0.0; 5],
                minimum_accepted_step_s: duration_s,
                maximum_accepted_step_s: duration_s,
                suggested_step_s: controls.initial_step_s,
                total_chemical_accepted_steps: totals.gas.chemical_accepted,
                total_chemical_rejected_steps: totals.gas.chemical_rejected,
            },
        ));
    }

    // Bind canonical keys once, without cloning source cards or rebuilding
    // keyed maps during retries. The prepared reactor remains mutably reusable.
    let mut species = Vec::with_capacity(reactor.bank().thermodynamics().components().len());
    for source in reactor.bank().thermodynamics().components() {
        poll(cancelled)?;
        species.push(source.species().clone());
    }
    let kinetics = reactor.kinetics_controls();
    let comparison = AdaptiveMaterialCompressionSolve {
        kinetics,
        volume: controls.volume,
        absolute_temperature_tolerance_k: controls.absolute_temperature_tolerance_k,
        absolute_amount_tolerance_mol: controls.absolute_amount_tolerance_mol,
        absolute_work_tolerance_j: controls.absolute_work_tolerance_j,
        relative_tolerance: controls.relative_tolerance,
        initial_step_s: controls.initial_step_s,
        minimum_step_s: controls.minimum_step_s,
        maximum_step_s: controls.maximum_step_s,
        maximum_intervals: controls.maximum_intervals,
        maximum_evaluations: controls.maximum_evaluations,
    };
    let log_ratio = fs_math::det::ln(ratio);
    let volume_at = |time: f64| {
        if time == duration_s {
            final_volume_m3
        } else {
            initial_volume * fs_math::det::exp(log_ratio * (time / duration_s))
        }
    };
    let (mut elapsed, mut step) = (0.0, controls.initial_step_s);
    let (mut accepted, mut rejected, mut evaluations) = (0, 0, 0);
    let mut chemical = [0_usize; 2];
    let mut maxima = [0.0_f64; 5];
    let (mut minimum, mut maximum) = (f64::INFINITY, 0.0_f64);
    while elapsed < duration_s {
        poll(cancelled)?;
        if accepted == controls.maximum_intervals || evaluations + 3 > controls.maximum_evaluations
        {
            return Err(Error::WorkLimit {
                elapsed_time_s: elapsed,
                accepted_intervals: accepted,
                rejected_intervals: rejected,
                evaluations,
            });
        }
        let boundary = partner.next_boundary(elapsed, duration_s);
        let end = if step >= boundary - elapsed {
            boundary
        } else {
            elapsed + step
        };
        let span = end - elapsed;
        let middle = elapsed + 0.5 * span;
        let first_midpoint = elapsed + 0.5 * (middle - elapsed);
        let second_midpoint = middle + 0.5 * (end - middle);
        if !(elapsed < first_midpoint
            && first_midpoint < middle
            && middle < second_midpoint
            && second_midpoint < end)
        {
            return Err(Error::TimeResolution {
                elapsed_time_s: elapsed,
                target_time_s: end,
            });
        }
        let endpoint = volume_at(end);
        let midpoint = volume_at(middle);
        let mut coarse = candidate.clone();
        let mut coarse_partner = partner;
        let mut coarse_totals = Totals::default();
        evaluations += 1;
        coarse_totals.interval(
            &mut reactor.reactor,
            &mut coarse,
            &mut coarse_partner,
            CompressionInterval {
                duration_s: span,
                midpoint_volume_m3: midpoint,
                final_volume_m3: endpoint,
            },
            IntervalTime {
                start_s: elapsed,
                end_s: end,
            },
            controls.volume,
            cancelled,
        )?;
        let mut fine = candidate.clone();
        let mut fine_partner = partner;
        let mut fine_totals = Totals::default();
        evaluations += 1;
        fine_totals.interval(
            &mut reactor.reactor,
            &mut fine,
            &mut fine_partner,
            CompressionInterval {
                duration_s: middle - elapsed,
                midpoint_volume_m3: volume_at(first_midpoint),
                final_volume_m3: midpoint,
            },
            IntervalTime {
                start_s: elapsed,
                end_s: middle,
            },
            controls.volume,
            cancelled,
        )?;
        evaluations += 1;
        fine_totals.interval(
            &mut reactor.reactor,
            &mut fine,
            &mut fine_partner,
            CompressionInterval {
                duration_s: end - middle,
                midpoint_volume_m3: volume_at(second_midpoint),
                final_volume_m3: endpoint,
            },
            IntervalTime {
                start_s: middle,
                end_s: end,
            },
            controls.volume,
            cancelled,
        )?;
        // Validate every physical account before comparing or accepting a state.
        coarse_totals.finish(
            coarse_partner,
            candidate.volume_m3(),
            endpoint,
            Time::new(span),
        )?;
        fine_totals.finish(
            fine_partner,
            candidate.volume_m3(),
            endpoint,
            Time::new(span),
        )?;
        for (total, counts) in chemical.iter_mut().zip([
            [
                coarse_totals.gas.chemical_accepted,
                fine_totals.gas.chemical_accepted,
            ],
            [
                coarse_totals.gas.chemical_rejected,
                fine_totals.gas.chemical_rejected,
            ],
        ]) {
            *total = counts[0]
                .checked_add(counts[1])
                .and_then(|n| total.checked_add(n))
                .ok_or(Error::Invalid("total chemical work count overflow"))?;
        }
        let original = compare(
            &species,
            [&candidate, &coarse, &fine],
            [coarse_totals.gas.work.0, fine_totals.gas.work.0],
            comparison,
            cancelled,
        )?;
        let errors = [
            original[0],
            original[1],
            original[2],
            scaled_difference(
                0.0,
                coarse_totals.thermal[0].0,
                fine_totals.thermal[0].0,
                controls.absolute_heat_tolerance_j,
                controls.relative_tolerance,
            )?,
            body_error(partner, coarse_partner, fine_partner, controls)?,
        ];
        let error = errors.into_iter().fold(0.0_f64, f64::max);
        if error <= 1.0 {
            totals.merge(fine_totals)?;
            candidate = fine;
            partner = fine_partner;
            elapsed = end;
            accepted += 1;
            minimum = minimum.min(span);
            maximum = maximum.max(span);
            for i in 0..5 {
                maxima[i] = maxima[i].max(errors[i]);
            }
            let factor = if error == 0.0 {
                2.0
            } else {
                (0.9 / error_root(error, kinetics.chemistry.method)).clamp(0.5, 2.0)
            };
            step = if span >= controls.maximum_step_s / factor {
                controls.maximum_step_s
            } else {
                (span * factor).max(controls.minimum_step_s)
            };
        } else {
            rejected += 1;
            let factor = (0.8 / error_root(error, kinetics.chemistry.method)).clamp(0.1, 0.5);
            let next = (span * factor).max(controls.minimum_step_s);
            if span <= controls.minimum_step_s || !positive(next) || next >= span {
                return Err(Error::Accuracy {
                    elapsed_time_s: elapsed,
                    step_s: span,
                    error_ratios: errors,
                });
            }
            step = next;
        }
    }
    let physical = totals.finish(partner, initial_volume, final_volume_m3, duration)?;
    Ok((
        candidate,
        AdaptiveMaterialThermalCompressionStep {
            physical,
            accepted_intervals: accepted,
            rejected_intervals: rejected,
            evaluations,
            maximum_accepted_error_ratios: maxima,
            minimum_accepted_step_s: minimum,
            maximum_accepted_step_s: maximum,
            suggested_step_s: step,
            total_chemical_accepted_steps: chemical[0],
            total_chemical_rejected_steps: chemical[1],
        },
    ))
}

// Subtract represented specific enthalpies before multiplying by original mass.
// Neither a large enthalpy datum nor an absolute H product sets the error scale.
fn body_error<B: Partner>(
    initial: B,
    coarse: B,
    fine: B,
    controls: AdaptiveMaterialThermalCompressionSolve,
) -> Result<f64, Error> {
    let Some([initial_h, mass]) = initial.body_inventory() else {
        return Ok(0.0);
    };
    let [coarse_h, _] = coarse
        .body_inventory()
        .ok_or(Error::Invalid("missing coarse body state"))?;
    let [fine_h, _] = fine
        .body_inventory()
        .ok_or(Error::Invalid("missing fine body state"))?;
    let change = ((Dd::from_f64(coarse_h) - Dd::from_f64(initial_h)) * Dd::from_f64(mass))
        .to_f64()
        .abs()
        .max(
            ((Dd::from_f64(fine_h) - Dd::from_f64(initial_h)) * Dd::from_f64(mass))
                .to_f64()
                .abs(),
        );
    let difference = ((Dd::from_f64(fine_h) - Dd::from_f64(coarse_h)) * Dd::from_f64(mass))
        .to_f64()
        .abs();
    let scale = controls.absolute_body_enthalpy_tolerance_j + controls.relative_tolerance * change;
    if !positive(scale) || !difference.is_finite() {
        return Err(Error::Invalid(
            "body enthalpy comparison is unrepresentable",
        ));
    }
    Ok(difference / scale)
}

fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), Error> {
    if cancelled() {
        Err(Error::Cancelled)
    } else {
        Ok(())
    }
}
