//! Conserved internal-energy/volume inversion using the original PT/PH owners.
//!
//! For each pressure, the existing PH solve imposes `h = u_target + p*v_target`.
//! The outer volume residual therefore closes both conserved coordinates at its
//! root. Away from that root this curve need not be monotone; only a declared
//! sign bracket is used. Finite phase-screen evidence is not upgraded to a
//! global entropy maximum or a unique unrestricted U,V equilibrium.

use super::{PhaseEnthalpyControls, ResolvedThermochemicalPhaseState};
use crate::thermochemical_gas::nonideal::ThermochemicalPengRobinsonModel;
use crate::thermochemical_gas::nonideal::flash::{
    CaloricSample, E,
    enthalpy::{EnthalpyWork, Target, spacing},
    failure, poll, positive,
};
use fs_qty::{Pressure, Temperature};
use fs_thermochem::eos::peng_robinson::phase_screen::{
    PengRobinsonPhaseScreenOptionsV1, PengRobinsonPhaseScreenV1,
};

/// Maximum pressure inversions in one conserved-energy/volume solve.
pub const MAX_ENERGY_VOLUME_PRESSURE_EVALUATIONS: usize = 128;
/// Maximum original PT selections across all nested pressure inversions.
pub const MAX_ENERGY_VOLUME_PT_EVALUATIONS: usize = 16_384;

/// Explicit physical closure tolerances, source interval and nested work limits.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PhaseEnergyVolumeControls {
    /// Increasing positive pressure endpoints [Pa]. Their solved volume
    /// residuals must enclose zero from positive to negative.
    pub pressure_bracket: [Pressure; 2],
    /// Absolute mass-specific volume residual tolerance [m³/kg].
    pub specific_volume_tolerance_m3_kg: f64,
    /// Absolute mass-specific internal-energy residual tolerance [J/kg].
    pub internal_energy_tolerance_j_kg: f64,
    /// Accepted pressure correction or remaining bracket tolerance [Pa].
    pub pressure_tolerance_pa: f64,
    /// Whole-call number of pressure/enthalpy inversions, including endpoints.
    pub maximum_pressure_evaluations: usize,
    /// Whole-call original PT work, including the reference and all inner solves.
    pub maximum_pt_evaluations: usize,
    /// Original phase policy, source temperatures and per-pressure PH work cap.
    /// Its enthalpy tolerance is tightened to reserve the energy closure budget.
    pub phase: PhaseEnthalpyControls,
}
impl PhaseEnergyVolumeControls {
    /// A caller-declared pressure bracket and the original phase/temperature policy.
    #[must_use]
    pub fn new(pressure_bracket: [Pressure; 2], phase: PhaseEnthalpyControls) -> Self {
        Self {
            pressure_bracket,
            specific_volume_tolerance_m3_kg: 1e-12,
            internal_energy_tolerance_j_kg: 1e-5,
            pressure_tolerance_pa: 1e-3,
            maximum_pressure_evaluations: 64,
            maximum_pt_evaluations: 4096,
            phase,
        }
    }

    /// Validate all nested limits before source evaluations or allocations.
    pub fn validate(self) -> Result<(), E> {
        self.phase.validate()?;
        let [lo, hi] = self.pressure_bracket.map(|p| p.value());
        if !positive(lo)
            || !positive(hi)
            || lo >= hi
            || !positive(self.specific_volume_tolerance_m3_kg)
            || !positive(self.internal_energy_tolerance_j_kg)
            || !positive(self.pressure_tolerance_pa)
            || self.maximum_pressure_evaluations == 0
            || self.maximum_pressure_evaluations > MAX_ENERGY_VOLUME_PRESSURE_EVALUATIONS
            || self.maximum_pt_evaluations == 0
            || self.maximum_pt_evaluations > MAX_ENERGY_VOLUME_PT_EVALUATIONS
        {
            return Err(E::Invalid("invalid bounded energy/volume flash controls"));
        }
        Ok(())
    }

    /// Conservative total cubic EOS ceiling, distinct from actual PT calls.
    pub fn maximum_eos_evaluation_bound(self) -> Result<usize, E> {
        self.validate()?;
        self.phase
            .phase_options()
            .maximum_eos_evaluation_bound()
            .map_err(E::Phase)?
            .checked_mul(self.maximum_pt_evaluations)
            .ok_or(E::Invalid("energy/volume EOS work bound overflow"))
    }
}

/// An admitted U,V endpoint with independently reconstructed conserved residuals.
#[derive(Clone, Debug, PartialEq)]
pub struct PhaseEnergyVolumeSolution {
    state: ResolvedThermochemicalPhaseState,
    energy_residual_j_kg: f64,
    volume_residual_m3_kg: f64,
    accepted_increment_j_kg: Option<f64>,
    entropy_change_j_kg_k: Option<f64>,
    pressure_evaluations: usize,
    pt_evaluations: usize,
    eos_evaluation_upper_bound: usize,
}
impl PhaseEnergyVolumeSolution {
    /// Complete original-source state and its finite phase-selection evidence.
    #[must_use]
    pub const fn state(&self) -> &ResolvedThermochemicalPhaseState {
        &self.state
    }
    /// Move the accepted state without another phase/source evaluation.
    #[must_use]
    pub fn into_state(self) -> ResolvedThermochemicalPhaseState {
        self.state
    }
    /// Direct original-law internal-energy coordinate minus target [J/kg].
    #[must_use]
    pub const fn energy_residual_j_kg(&self) -> f64 {
        self.energy_residual_j_kg
    }
    /// Original phase volume per feed mass minus target [m³/kg].
    #[must_use]
    pub const fn volume_residual_m3_kg(&self) -> f64 {
        self.volume_residual_m3_kg
    }
    /// Source-integrated actual Delta-u, present only for the increment API.
    #[must_use]
    pub const fn accepted_increment_j_kg(&self) -> Option<f64> {
        self.accepted_increment_j_kg
    }
    /// Original reference-invariant Delta-s, present only for the increment API.
    /// A process consumer must impose its own heat/work/entropy convention.
    #[must_use]
    pub const fn entropy_change_j_kg_k(&self) -> Option<f64> {
        self.entropy_change_j_kg_k
    }
    /// Actual attempted pressure/enthalpy inversions, including endpoints.
    #[must_use]
    pub const fn pressure_evaluations(&self) -> usize {
        self.pressure_evaluations
    }
    /// Actual original PT calls across all inversions, including the reference.
    #[must_use]
    pub const fn pt_evaluations(&self) -> usize {
        self.pt_evaluations
    }
    /// Observed screen work plus conservative admitted split work ceilings.
    #[must_use]
    pub const fn eos_evaluation_upper_bound(&self) -> usize {
        self.eos_evaluation_upper_bound
    }
}

pub(super) type Sample = CaloricSample<PengRobinsonPhaseScreenV1>;

enum EnergyTarget {
    Absolute(f64),
    Increment { initial: Sample, requested: f64 },
}
impl EnergyTarget {
    fn requested(&self) -> f64 {
        match self {
            Self::Absolute(value) => *value,
            Self::Increment { requested, .. } => *requested,
        }
    }
    fn tolerance(&self, configured: f64) -> Result<f64, E> {
        let value = match self {
            Self::Increment { requested, .. } if *requested != 0.0 => {
                configured.min(0.25 * requested.abs())
            }
            _ => configured,
        };
        if positive(value * 0.125) {
            Ok(value)
        } else {
            Err(E::Invalid(
                "energy increment is below represented solve resolution",
            ))
        }
    }
}

#[derive(Default)]
pub(super) struct Work {
    pub pressure: usize,
    pub pt: usize,
    pub eos_bound: usize,
}
pub(super) struct PressureSample {
    pub sample: Sample,
    // The public owner retains either internal-energy or entropy units.
    pub coordinate: f64,
    pub volume_residual: f64,
    pub entropy_change: Option<f64>,
}

// One pressure bracket engine serves both original U,V and S,V constraints.
// Each physical owner supplies and independently checks its caloric coordinate.
pub(super) struct VolumeRootControls {
    pub pressure_bracket: [Pressure; 2],
    pub volume_tolerance: f64,
    pub pressure_tolerance: f64,
    pub coordinate_target: f64,
    pub coordinate_tolerance: f64,
}

impl ThermochemicalPengRobinsonModel {
    /// Recover original-source pressure, temperature and phases at specified
    /// mass-specific internal energy and volume. The whole declared rectangle
    /// must be supported by the source and admit every sampled PT state; failed
    /// stability searches and phase-amount ambiguity refuse without fallback.
    /// Absolute reference spacing must resolve the requested energy tolerance.
    pub fn phase_equilibrium_internal_energy_volume(
        &self,
        specific_volume_m3_kg: f64,
        specific_internal_energy_j_kg: f64,
        controls: PhaseEnergyVolumeControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PhaseEnergyVolumeSolution, E> {
        self.admit_energy_volume(specific_volume_m3_kg, controls, None, &mut cancelled)?;
        if !specific_internal_energy_j_kg.is_finite()
            || spacing(specific_internal_energy_j_kg) > controls.internal_energy_tolerance_j_kg
        {
            return Err(E::Invalid(
                "absolute internal energy cannot resolve tolerance; use source increments",
            ));
        }
        self.solve_energy_volume(
            specific_volume_m3_kg,
            EnergyTarget::Absolute(specific_internal_energy_j_kg),
            controls,
            Work::default(),
            &mut cancelled,
        )
    }

    /// Change original-source internal energy and optionally volume at fixed
    /// species inventory. A fixed-volume vessel supplies its unchanged volume
    /// per mass. No boundary-work or process-entropy law is inferred if volume
    /// changes. Ideal species increments are integrated before departure
    /// differences, so formation and entropy datums cancel before rounding.
    /// Zero increment at identical volume returns an in-bracket reference
    /// exactly after one PT call.
    #[allow(clippy::too_many_arguments)]
    pub fn phase_equilibrium_internal_energy_volume_increment(
        &self,
        initial_temperature: Temperature,
        initial_pressure: Pressure,
        final_specific_volume_m3_kg: f64,
        change_j_kg: f64,
        controls: PhaseEnergyVolumeControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PhaseEnergyVolumeSolution, E> {
        self.admit_energy_volume(
            final_specific_volume_m3_kg,
            controls,
            Some(initial_temperature),
            &mut cancelled,
        )?;
        if !positive(initial_pressure.value()) || !change_j_kg.is_finite() {
            return Err(E::Invalid(
                "finite energy increment and positive reference pressure required",
            ));
        }
        let mut work = Work::default();
        let initial = self.energy_volume_sample(
            initial_temperature,
            initial_pressure,
            controls,
            &mut work,
            &mut cancelled,
        )?;
        self.solve_energy_volume(
            final_specific_volume_m3_kg,
            EnergyTarget::Increment {
                initial,
                requested: change_j_kg,
            },
            controls,
            work,
            &mut cancelled,
        )
    }

    fn admit_energy_volume(
        &self,
        volume: f64,
        controls: PhaseEnergyVolumeControls,
        initial: Option<Temperature>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), E> {
        controls.validate()?;
        if !positive(volume) || spacing(volume) > controls.specific_volume_tolerance_m3_kg {
            return Err(E::Invalid(
                "positive representable specific volume required",
            ));
        }
        self.admit_enthalpy_flash(
            controls.pressure_bracket[0],
            controls.phase.enthalpy,
            initial,
            cancelled,
        )
    }

    fn energy_volume_sample(
        &self,
        temperature: Temperature,
        pressure: Pressure,
        controls: PhaseEnergyVolumeControls,
        work: &mut Work,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Sample, E> {
        self.volume_caloric_sample(
            temperature,
            pressure,
            controls.phase.phase_options(),
            controls.maximum_pt_evaluations,
            work,
            cancelled,
        )
    }

    pub(super) fn volume_caloric_sample(
        &self,
        temperature: Temperature,
        pressure: Pressure,
        phase_options: PengRobinsonPhaseScreenOptionsV1,
        maximum_pt_evaluations: usize,
        work: &mut Work,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Sample, E> {
        poll(cancelled)?;
        if work.pt >= maximum_pt_evaluations {
            return Err(E::WorkLimit {
                flash_evaluations: work.pt,
            });
        }
        work.pt += 1;
        let sample = self.caloric_phase_sample(temperature, pressure, phase_options, cancelled)?;
        work.eos_bound = work
            .eos_bound
            .checked_add(sample.eos_evaluation_upper_bound)
            .ok_or(E::Invalid("caloric/volume EOS work accounting overflow"))?;
        Ok(sample)
    }

    fn internal_energy_coordinate(&self, sample: &Sample, target: &EnergyTarget) -> Result<f64, E> {
        let value = match target {
            EnergyTarget::Absolute(_) => sample.specific_u,
            EnergyTarget::Increment { initial, .. } => {
                let ideal = self
                    .source()
                    .thermodynamics()
                    .specific_internal_energy_change(initial.temperature, sample.temperature)
                    .map_err(|error| failure("original flash internal-energy increment", error))?
                    .value();
                ideal
                    + (sample.residual_u - initial.residual_u)
                        / self.mechanics().molar_mass().value()
            }
        };
        if value.is_finite() {
            Ok(value)
        } else {
            Err(E::Invalid(
                "energy/volume internal-energy coordinate overflow",
            ))
        }
    }

    #[allow(clippy::too_many_arguments)]
    fn energy_volume_pressure(
        &self,
        pressure: f64,
        volume: f64,
        target: &EnergyTarget,
        controls: PhaseEnergyVolumeControls,
        energy_tolerance: f64,
        work: &mut Work,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<PressureSample, E> {
        poll(cancelled)?;
        if work.pressure >= controls.maximum_pressure_evaluations {
            return Err(E::PressureWorkLimit {
                pressure_evaluations: work.pressure,
            });
        }
        work.pressure += 1;
        let (enthalpy_target, requested) = match target {
            EnergyTarget::Absolute(u) => {
                let h = u + pressure * volume;
                (Target::Absolute(h), h)
            }
            EnergyTarget::Increment { initial, requested } => {
                let initial_volume = initial.molar_volume / self.mechanics().molar_mass().value();
                // Difference of p*v without subtracting two large full products.
                let pv_change = pressure * (volume - initial_volume)
                    + (pressure - initial.pressure.value()) * initial_volume;
                let h = requested + pv_change;
                (
                    Target::Increment {
                        initial: initial.clone(),
                        requested: h,
                    },
                    h,
                )
            }
        };
        let mut inner = controls.phase.enthalpy;
        inner.enthalpy_tolerance_j_kg = inner.enthalpy_tolerance_j_kg.min(0.125 * energy_tolerance);
        if !requested.is_finite() || spacing(requested) > inner.enthalpy_tolerance_j_kg {
            return Err(E::Invalid(
                "energy/volume enthalpy coordinate cannot resolve tolerance",
            ));
        }
        let outcome = self.solve_flash_enthalpy(
            Pressure::new(pressure),
            enthalpy_target,
            inner,
            EnthalpyWork::default(),
            &mut |t, p, c| self.energy_volume_sample(t, p, controls, work, c),
            cancelled,
        );
        // Even a per-pressure PH refusal reports the complete UV call's PT
        // count, including its original reference and earlier pressure solves.
        let outcome = outcome.map_err(|error| match error {
            E::WorkLimit { .. } => E::WorkLimit {
                flash_evaluations: work.pt,
            },
            other => other,
        })?;
        let energy_coordinate = self.internal_energy_coordinate(&outcome.sample, target)?;
        let actual_volume = outcome.sample.molar_volume / self.mechanics().molar_mass().value();
        if !positive(actual_volume)
            || matches!(target, EnergyTarget::Absolute(_))
                && spacing(energy_coordinate) > energy_tolerance
        {
            return Err(E::Invalid(
                "energy/volume original-state resolution is insufficient",
            ));
        }
        Ok(PressureSample {
            sample: outcome.sample,
            coordinate: energy_coordinate,
            volume_residual: actual_volume - volume,
            entropy_change: outcome.entropy_change_j_kg_k,
        })
    }

    fn finish_energy_volume(
        &self,
        accepted: PressureSample,
        target: &EnergyTarget,
        work: Work,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<PhaseEnergyVolumeSolution, E> {
        let state = self.publish_caloric_phase(accepted.sample, cancelled)?;
        poll(cancelled)?;
        Ok(PhaseEnergyVolumeSolution {
            state,
            energy_residual_j_kg: accepted.coordinate - target.requested(),
            volume_residual_m3_kg: accepted.volume_residual,
            accepted_increment_j_kg: matches!(target, EnergyTarget::Increment { .. })
                .then_some(accepted.coordinate),
            entropy_change_j_kg_k: accepted.entropy_change,
            pressure_evaluations: work.pressure,
            pt_evaluations: work.pt,
            eos_evaluation_upper_bound: work.eos_bound,
        })
    }

    fn solve_energy_volume(
        &self,
        volume: f64,
        target: EnergyTarget,
        controls: PhaseEnergyVolumeControls,
        mut work: Work,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<PhaseEnergyVolumeSolution, E> {
        let energy_tolerance = target.tolerance(controls.internal_energy_tolerance_j_kg)?;
        if let EnergyTarget::Increment { initial, requested } = &target {
            if *requested == 0.0
                && volume == initial.molar_volume / self.mechanics().molar_mass().value()
                && initial.temperature.value()
                    >= controls.phase.enthalpy.temperature_bracket[0].value()
                && initial.temperature.value()
                    <= controls.phase.enthalpy.temperature_bracket[1].value()
                && initial.pressure.value() >= controls.pressure_bracket[0].value()
                && initial.pressure.value() <= controls.pressure_bracket[1].value()
            {
                return self.finish_energy_volume(
                    PressureSample {
                        sample: initial.clone(),
                        coordinate: 0.0,
                        volume_residual: 0.0,
                        entropy_change: Some(0.0),
                    },
                    &target,
                    work,
                    cancelled,
                );
            }
        }
        // Reserve energy room for the outer volume residual, independently of
        // the final direct U check. This is a numerical safeguard, not a Cp.
        let volume_tolerance = controls
            .specific_volume_tolerance_m3_kg
            .min(0.25 * energy_tolerance / controls.pressure_bracket[1].value());
        if !positive(volume_tolerance) {
            return Err(E::Invalid("energy/volume closure tolerance underflow"));
        }
        let accepted = solve_volume_bracket(
            VolumeRootControls {
                pressure_bracket: controls.pressure_bracket,
                volume_tolerance,
                pressure_tolerance: controls.pressure_tolerance_pa,
                coordinate_target: target.requested(),
                coordinate_tolerance: energy_tolerance,
            },
            &mut |pressure, cancelled| {
                self.energy_volume_pressure(
                    pressure,
                    volume,
                    &target,
                    controls,
                    energy_tolerance,
                    &mut work,
                    cancelled,
                )
            },
            cancelled,
        )?;
        self.finish_energy_volume(accepted, &target, work, cancelled)
    }
}

pub(super) fn solve_volume_bracket<C: FnMut() -> bool>(
    controls: VolumeRootControls,
    evaluate: &mut impl FnMut(f64, &mut C) -> Result<PressureSample, E>,
    cancelled: &mut C,
) -> Result<PressureSample, E> {
    let mut lower = evaluate(controls.pressure_bracket[0].value(), cancelled)?;
    let mut upper = evaluate(controls.pressure_bracket[1].value(), cancelled)?;
    if lower.volume_residual < 0.0 || upper.volume_residual > 0.0 {
        return Err(E::OutsideVolumeBracket);
    }
    let mut last_pressure = f64::NAN;
    let mut last_step = f64::INFINITY;
    let mut iteration = 0;
    loop {
        poll(cancelled)?;
        let lo = lower.sample.pressure.value();
        let hi = upper.sample.pressure.value();
        let lower_best = lower.volume_residual.abs() <= upper.volume_residual.abs();
        let best = if lower_best { &lower } else { &upper };
        let best_pressure = best.sample.pressure.value();
        let coordinate_residual = best.coordinate - controls.coordinate_target;
        if best.volume_residual.abs() <= controls.volume_tolerance
            && coordinate_residual.abs() <= controls.coordinate_tolerance
            && (best.volume_residual == 0.0
                || hi - lo <= controls.pressure_tolerance
                || best_pressure == last_pressure && last_step <= controls.pressure_tolerance)
        {
            return Ok(if lower_best { lower } else { upper });
        }
        let fraction = lower.volume_residual / (lower.volume_residual - upper.volume_residual);
        let secant = lo + fraction * (hi - lo);
        let pressure = if iteration % 4 != 3 && secant.is_finite() && secant > lo && secant < hi {
            secant
        } else {
            lo + 0.5 * (hi - lo)
        };
        if pressure <= lo || pressure >= hi {
            return Err(E::UnresolvedPressure);
        }
        last_pressure = pressure;
        last_step = (pressure - best_pressure).abs();
        let candidate = evaluate(pressure, cancelled)?;
        if candidate.volume_residual >= 0.0 {
            lower = candidate;
        } else {
            upper = candidate;
        }
        iteration += 1;
    }
}
