//! Caloric flash with one shared bounded temperature root engine.
//!
//! A safeguarded secant solve uses the actual original PT equilibrium calorics;
//! its secant is numerical and is never exposed as a mixed heat capacity.

use super::{CaloricSample, E, ResolvedThermochemicalFlashState, poll, positive};
use crate::thermochemical_gas::{failure, nonideal::ThermochemicalPengRobinsonModel};
use fs_qty::{Pressure, Temperature};
use fs_thermochem::eos::peng_robinson::flash::{
    MAX_PENG_ROBINSON_FLASH_ITERATIONS_V1, MAX_RACHFORD_RICE_ITERATIONS_V1,
    PengRobinsonFlashOptionsV1,
};
use fs_thermochem::kinetics::{KineticsErrorV1, admit_caloric_interval};

/// Maximum number of original PT calls in one enthalpy inversion.
pub const MAX_ENTHALPY_FLASH_EVALUATIONS: usize = 128;

/// Bounded physical tolerances and complete nested phase-solve work policy.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FlashEnthalpyControls {
    /// Increasing original-source temperatures. The strict two-phase API requires
    /// separated endpoints when searching; the screened API admits either phase.
    /// An exact fixed point needs its reference inside the declared bracket.
    pub temperature_bracket: [Temperature; 2],
    /// Absolute specific-enthalpy residual tolerance [J/kg].
    pub enthalpy_tolerance_j_kg: f64,
    /// Temperature correction or bracket tolerance [K].
    pub temperature_tolerance_k: f64,
    /// Whole-call PT work, including endpoints and an increment reference state.
    pub maximum_flash_evaluations: usize,
    /// Original PT phase-equilibrium controls, unchanged by this energy owner.
    pub flash_options: PengRobinsonFlashOptionsV1,
}
impl FlashEnthalpyControls {
    /// Default physical tolerances with an explicit caller-supplied phase bracket.
    #[must_use]
    pub fn new(temperature_bracket: [Temperature; 2]) -> Self {
        Self {
            temperature_bracket,
            enthalpy_tolerance_j_kg: 1e-5,
            temperature_tolerance_k: 1e-8,
            maximum_flash_evaluations: 64,
            flash_options: PengRobinsonFlashOptionsV1 {
                log_k_tolerance: 1e-12,
                ..Default::default()
            },
        }
    }
    /// Admit the numerical limits; species/phase domains are checked by their owners.
    pub fn validate(self) -> Result<(), E> {
        let phase = self.flash_options;
        if !CaloricRootControls::from(self).is_valid()
            || phase.maximum_iterations == 0
            || phase.maximum_iterations > MAX_PENG_ROBINSON_FLASH_ITERATIONS_V1
            || phase.rachford_rice_iterations == 0
            || phase.rachford_rice_iterations > MAX_RACHFORD_RICE_ITERATIONS_V1
            || !positive(phase.log_k_tolerance)
        {
            return Err(E::Invalid("invalid bounded enthalpy flash controls"));
        }
        Ok(())
    }
    /// Worst-case cubic EOS work of the entire call, not an observed iteration count.
    /// The original PT owner permits at most `6*maximum_iterations+2` EOS calls.
    pub fn maximum_eos_evaluation_bound(self) -> Result<usize, E> {
        self.validate()?;
        Ok(self.maximum_flash_evaluations * (6 * self.flash_options.maximum_iterations + 2))
    }
}

/// Accepted energy-constrained phase state and independently evaluated closure.
#[derive(Clone, Debug, PartialEq)]
pub struct FlashEnthalpySolution {
    state: ResolvedThermochemicalFlashState,
    residual_j_kg: f64,
    accepted_increment_j_kg: Option<f64>,
    entropy_change_j_kg_k: Option<f64>,
    flash_evaluations: usize,
    eos_evaluation_upper_bound: usize,
}
impl FlashEnthalpySolution {
    /// Complete original-source two-phase endpoint.
    #[must_use]
    pub const fn state(&self) -> &ResolvedThermochemicalFlashState {
        &self.state
    }
    /// Move the accepted state without rebuilding source receipts.
    #[must_use]
    pub fn into_state(self) -> ResolvedThermochemicalFlashState {
        self.state
    }
    /// Evaluated enthalpy coordinate minus the requested target [J/kg].
    #[must_use]
    pub const fn residual_j_kg(&self) -> f64 {
        self.residual_j_kg
    }
    /// Source-integrated actual Delta h, present only for the increment API.
    #[must_use]
    pub const fn accepted_increment_j_kg(&self) -> Option<f64> {
        self.accepted_increment_j_kg
    }
    /// Original-reference-invariant Delta s, present only for the increment API.
    /// An application imposing a process entropy inequality must check this value.
    #[must_use]
    pub const fn entropy_change_j_kg_k(&self) -> Option<f64> {
        self.entropy_change_j_kg_k
    }
    /// Actual original PT flash calls, including reference and bracket endpoints.
    #[must_use]
    pub const fn flash_evaluations(&self) -> usize {
        self.flash_evaluations
    }
    /// Conservative EOS-work bound for those calls, not actual EOS evaluations.
    #[must_use]
    pub const fn eos_evaluation_upper_bound(&self) -> usize {
        self.eos_evaluation_upper_bound
    }
}

pub(super) enum Target<P> {
    Absolute(f64),
    Increment {
        initial: CaloricSample<P>,
        requested: f64,
    },
}
impl<P> Target<P> {
    pub fn requested(&self) -> f64 {
        match self {
            Self::Absolute(value) => *value,
            Self::Increment { requested, .. } => *requested,
        }
    }
}

pub(super) fn spacing(value: f64) -> f64 {
    (value.next_up() - value)
        .abs()
        .max((value - value.next_down()).abs())
}

#[derive(Clone, Copy, Default)]
pub(super) struct EnthalpyWork {
    pub calls: usize,
    pub eos_upper_bound: usize,
}

// Both caloric constraints use the same safeguarded temperature iteration and
// work ledger. The coordinate retains its physical units at the public boundary.
#[derive(Clone, Copy)]
pub(super) enum CaloricCoordinate {
    Enthalpy,
    Entropy,
}

#[derive(Clone, Copy)]
pub(super) struct CaloricRootControls {
    pub temperature_bracket: [Temperature; 2],
    pub residual_tolerance: f64,
    pub temperature_tolerance_k: f64,
    pub maximum_pt_evaluations: usize,
}
impl CaloricRootControls {
    pub fn is_valid(self) -> bool {
        let [lower, upper] = self.temperature_bracket.map(|t| t.value());
        positive(lower)
            && positive(upper)
            && lower < upper
            && positive(self.residual_tolerance)
            && positive(self.temperature_tolerance_k)
            && self.maximum_pt_evaluations > 0
            && self.maximum_pt_evaluations <= MAX_ENTHALPY_FLASH_EVALUATIONS
    }
}
impl From<FlashEnthalpyControls> for CaloricRootControls {
    fn from(value: FlashEnthalpyControls) -> Self {
        Self {
            temperature_bracket: value.temperature_bracket,
            residual_tolerance: value.enthalpy_tolerance_j_kg,
            temperature_tolerance_k: value.temperature_tolerance_k,
            maximum_pt_evaluations: value.maximum_flash_evaluations,
        }
    }
}

pub(super) struct CaloricRoot<P> {
    pub sample: CaloricSample<P>,
    pub coordinate: f64,
    pub work: EnthalpyWork,
}

pub(super) struct EnthalpyOutcome<P> {
    pub sample: CaloricSample<P>,
    pub residual_j_kg: f64,
    pub accepted_increment_j_kg: Option<f64>,
    pub entropy_change_j_kg_k: Option<f64>,
    pub flash_evaluations: usize,
    pub eos_evaluation_upper_bound: usize,
}

impl ThermochemicalPengRobinsonModel {
    /// Recover a two-phase temperature/composition at fixed pressure and full
    /// specific enthalpy. Single-phase or unresolved bracket endpoints refuse;
    /// absence of a split is not treated as a global stability certificate.
    /// Use the increment API when absolute reference spacing exceeds tolerance.
    pub fn two_phase_flash_pressure_enthalpy(
        &self,
        pressure: Pressure,
        specific_enthalpy_j_kg: f64,
        controls: FlashEnthalpyControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<FlashEnthalpySolution, E> {
        self.admit_enthalpy_flash(pressure, controls, None, &mut cancelled)?;
        if !specific_enthalpy_j_kg.is_finite()
            || spacing(specific_enthalpy_j_kg) > controls.enthalpy_tolerance_j_kg
        {
            return Err(E::Invalid(
                "absolute enthalpy cannot resolve tolerance; use source increments",
            ));
        }
        let outcome = self.solve_flash_enthalpy(
            pressure,
            Target::Absolute(specific_enthalpy_j_kg),
            controls,
            EnthalpyWork::default(),
            &mut |t, p, c| self.caloric_flash_sample(t, p, controls.flash_options, c),
            &mut cancelled,
        )?;
        self.publish_enthalpy_flash(outcome, &mut cancelled)
    }

    /// Change pressure and/or enthalpy while retaining every feed species amount.
    /// Both reference and final states are resolved by the original two-phase PT
    /// owner. Delta h integrates the original NASA laws before adding departure
    /// differences, so species formation datums cancel before rounding. Delta h=0
    /// describes an adiabatic throttle only when kinetic/potential changes are
    /// negligible; this API imposes no process entropy or flow-rate model.
    /// At unchanged pressure and zero increment, a reference inside the declared
    /// bracket is returned exactly after one PT call, without endpoint searches.
    #[allow(clippy::too_many_arguments)]
    pub fn two_phase_flash_pressure_enthalpy_increment(
        &self,
        initial_temperature: Temperature,
        initial_pressure: Pressure,
        final_pressure: Pressure,
        change_j_kg: f64,
        controls: FlashEnthalpyControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<FlashEnthalpySolution, E> {
        self.admit_enthalpy_flash(
            final_pressure,
            controls,
            Some(initial_temperature),
            &mut cancelled,
        )?;
        if !positive(initial_pressure.value()) || !change_j_kg.is_finite() {
            return Err(E::Invalid(
                "finite enthalpy increment and positive reference pressure required",
            ));
        }
        let mut work = EnthalpyWork::default();
        let initial = self.enthalpy_flash_sample(
            initial_temperature,
            initial_pressure,
            controls,
            &mut work,
            &mut |t, p, c| self.caloric_flash_sample(t, p, controls.flash_options, c),
            &mut cancelled,
        )?;
        let outcome = self.solve_flash_enthalpy(
            final_pressure,
            Target::Increment {
                initial,
                requested: change_j_kg,
            },
            controls,
            work,
            &mut |t, p, c| self.caloric_flash_sample(t, p, controls.flash_options, c),
            &mut cancelled,
        )?;
        self.publish_enthalpy_flash(outcome, &mut cancelled)
    }

    pub(super) fn admit_enthalpy_flash(
        &self,
        pressure: Pressure,
        controls: FlashEnthalpyControls,
        initial: Option<Temperature>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), E> {
        controls.validate()?;
        self.admit_caloric_flash(pressure, controls.temperature_bracket, initial, cancelled)
    }

    pub(super) fn admit_caloric_flash(
        &self,
        pressure: Pressure,
        temperature_bracket: [Temperature; 2],
        initial: Option<Temperature>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), E> {
        poll(cancelled)?;
        if !positive(pressure.value()) || initial.is_some_and(|t| !positive(t.value())) {
            return Err(E::Invalid(
                "positive flash pressure and reference temperature required",
            ));
        }
        let [mut lower, mut upper] = temperature_bracket;
        if let Some(t) = initial {
            lower = Temperature::new(lower.value().min(t.value()));
            upper = Temperature::new(upper.value().max(t.value()));
        }
        for species in self.source().thermodynamics().components() {
            admit_caloric_interval(species, lower, upper, &mut *cancelled).map_err(|error| {
                match error {
                    KineticsErrorV1::Cancelled => E::Cancelled,
                    other => {
                        E::Material(failure("flash continuous original caloric interval", other))
                    }
                }
            })?;
        }
        poll(cancelled)
    }

    pub(super) fn enthalpy_flash_sample<P, C: FnMut() -> bool>(
        &self,
        t: Temperature,
        p: Pressure,
        controls: FlashEnthalpyControls,
        work: &mut EnthalpyWork,
        evaluate: &mut impl FnMut(Temperature, Pressure, &mut C) -> Result<CaloricSample<P>, E>,
        cancelled: &mut C,
    ) -> Result<CaloricSample<P>, E> {
        self.bounded_caloric_sample(
            t,
            p,
            controls.maximum_flash_evaluations,
            work,
            evaluate,
            cancelled,
        )
    }

    pub(super) fn bounded_caloric_sample<P, C: FnMut() -> bool>(
        &self,
        t: Temperature,
        p: Pressure,
        maximum_pt_evaluations: usize,
        work: &mut EnthalpyWork,
        evaluate: &mut impl FnMut(Temperature, Pressure, &mut C) -> Result<CaloricSample<P>, E>,
        cancelled: &mut C,
    ) -> Result<CaloricSample<P>, E> {
        poll(cancelled)?;
        if work.calls >= maximum_pt_evaluations {
            return Err(E::WorkLimit {
                flash_evaluations: work.calls,
            });
        }
        work.calls += 1;
        let sample = evaluate(t, p, cancelled)?;
        work.eos_upper_bound += sample.eos_evaluation_upper_bound;
        Ok(sample)
    }

    pub(super) fn caloric_enthalpy_change<P>(
        &self,
        initial: &CaloricSample<P>,
        sample: &CaloricSample<P>,
    ) -> Result<f64, E> {
        let ideal = self
            .source()
            .thermodynamics()
            .specific_enthalpy_change(initial.temperature, sample.temperature)
            .map_err(|error| failure("original flash enthalpy increment", error))?
            .value();
        let change = ideal
            + (sample.residual_h - initial.residual_h) / self.mechanics().molar_mass().value();
        if change.is_finite() {
            Ok(change)
        } else {
            Err(E::Invalid("flash enthalpy increment overflow"))
        }
    }

    pub(super) fn caloric_internal_energy_change<P>(
        &self,
        initial: &CaloricSample<P>,
        sample: &CaloricSample<P>,
    ) -> Result<f64, E> {
        let ideal = self
            .source()
            .thermodynamics()
            .specific_internal_energy_change(initial.temperature, sample.temperature)
            .map_err(|error| failure("original phase internal-energy increment", error))?
            .value();
        let change = ideal
            + (sample.residual_u - initial.residual_u) / self.mechanics().molar_mass().value();
        if change.is_finite() {
            Ok(change)
        } else {
            Err(E::Invalid(
                "original phase internal-energy increment overflow",
            ))
        }
    }

    pub(super) fn caloric_entropy_change<P>(
        &self,
        initial: &CaloricSample<P>,
        sample: &CaloricSample<P>,
    ) -> Result<f64, E> {
        let ideal = self
            .source()
            .thermodynamics()
            .specific_entropy_change(
                initial.temperature,
                initial.pressure,
                sample.temperature,
                sample.pressure,
            )
            .map_err(|error| failure("original flash entropy increment", error))?
            .value();
        let change = ideal
            + (sample.entropy_correction - initial.entropy_correction)
                / self.mechanics().molar_mass().value();
        if change.is_finite() {
            Ok(change)
        } else {
            Err(E::Invalid("flash entropy increment overflow"))
        }
    }

    fn caloric_coordinate<P>(
        &self,
        sample: &CaloricSample<P>,
        target: &Target<P>,
        coordinate: CaloricCoordinate,
    ) -> Result<f64, E> {
        let value = match (coordinate, target) {
            (CaloricCoordinate::Enthalpy, Target::Absolute(_)) => sample.specific_h,
            (CaloricCoordinate::Entropy, Target::Absolute(_)) => sample.specific_s,
            (CaloricCoordinate::Enthalpy, Target::Increment { initial, .. }) => {
                self.caloric_enthalpy_change(initial, sample)?
            }
            (CaloricCoordinate::Entropy, Target::Increment { initial, .. }) => {
                self.caloric_entropy_change(initial, sample)?
            }
        };
        if value.is_finite() {
            Ok(value)
        } else {
            Err(E::Invalid("flash caloric coordinate overflow"))
        }
    }

    fn finish_enthalpy_flash<P>(
        &self,
        sample: CaloricSample<P>,
        target: &Target<P>,
        coordinate: f64,
        work: EnthalpyWork,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<EnthalpyOutcome<P>, E> {
        let entropy_change = if let Target::Increment { initial, .. } = target {
            Some(self.caloric_entropy_change(initial, &sample)?)
        } else {
            None
        };
        let increment = matches!(target, Target::Increment { .. }).then_some(coordinate);
        poll(cancelled)?;
        Ok(EnthalpyOutcome {
            sample,
            residual_j_kg: coordinate - target.requested(),
            accepted_increment_j_kg: increment,
            entropy_change_j_kg_k: entropy_change,
            flash_evaluations: work.calls,
            eos_evaluation_upper_bound: work.eos_upper_bound,
        })
    }

    fn publish_enthalpy_flash(
        &self,
        outcome: EnthalpyOutcome<fs_thermochem::eos::peng_robinson::flash::PengRobinsonFlashV1>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<FlashEnthalpySolution, E> {
        let state = self.publish_caloric_flash(outcome.sample, cancelled)?;
        poll(cancelled)?;
        Ok(FlashEnthalpySolution {
            state,
            residual_j_kg: outcome.residual_j_kg,
            accepted_increment_j_kg: outcome.accepted_increment_j_kg,
            entropy_change_j_kg_k: outcome.entropy_change_j_kg_k,
            flash_evaluations: outcome.flash_evaluations,
            eos_evaluation_upper_bound: outcome.eos_evaluation_upper_bound,
        })
    }

    pub(super) fn solve_flash_enthalpy<P: Clone, C: FnMut() -> bool>(
        &self,
        pressure: Pressure,
        target: Target<P>,
        controls: FlashEnthalpyControls,
        work: EnthalpyWork,
        evaluate: &mut impl FnMut(Temperature, Pressure, &mut C) -> Result<CaloricSample<P>, E>,
        cancelled: &mut C,
    ) -> Result<EnthalpyOutcome<P>, E> {
        let root = self.solve_flash_coordinate(
            pressure,
            &target,
            CaloricCoordinate::Enthalpy,
            controls.into(),
            work,
            evaluate,
            cancelled,
        )?;
        self.finish_enthalpy_flash(root.sample, &target, root.coordinate, root.work, cancelled)
    }

    #[allow(clippy::too_many_arguments)]
    pub(super) fn solve_flash_coordinate<P: Clone, C: FnMut() -> bool>(
        &self,
        pressure: Pressure,
        target: &Target<P>,
        coordinate: CaloricCoordinate,
        controls: CaloricRootControls,
        mut work: EnthalpyWork,
        evaluate: &mut impl FnMut(Temperature, Pressure, &mut C) -> Result<CaloricSample<P>, E>,
        cancelled: &mut C,
    ) -> Result<CaloricRoot<P>, E> {
        if let Target::Increment { initial, requested } = target {
            let temperature = initial.temperature.value();
            if *requested == 0.0
                && initial.pressure == pressure
                && temperature >= controls.temperature_bracket[0].value()
                && temperature <= controls.temperature_bracket[1].value()
            {
                return Ok(CaloricRoot {
                    sample: initial.clone(),
                    coordinate: 0.0,
                    work,
                });
            }
        }
        let mut lower = self.bounded_caloric_sample(
            controls.temperature_bracket[0],
            pressure,
            controls.maximum_pt_evaluations,
            &mut work,
            evaluate,
            cancelled,
        )?;
        let mut upper = self.bounded_caloric_sample(
            controls.temperature_bracket[1],
            pressure,
            controls.maximum_pt_evaluations,
            &mut work,
            evaluate,
            cancelled,
        )?;
        let mut lower_value = self.caloric_coordinate(&lower, target, coordinate)?;
        let mut upper_value = self.caloric_coordinate(&upper, target, coordinate)?;
        if matches!(target, Target::Absolute(_))
            && [lower_value, upper_value]
                .into_iter()
                .any(|value| spacing(value) > controls.residual_tolerance)
        {
            return Err(E::Invalid(
                "absolute phase calorics cannot resolve tolerance; use source increments",
            ));
        }
        if !positive(upper_value - lower_value) {
            return Err(E::Invalid(
                "resolved phase caloric coordinate must increase across the source interval",
            ));
        }
        let requested = target.requested();
        if requested < lower_value || requested > upper_value {
            return Err(match coordinate {
                CaloricCoordinate::Enthalpy => E::OutsideEnthalpyBracket,
                CaloricCoordinate::Entropy => E::OutsideEntropyBracket,
            });
        }
        let mut last_candidate_temperature = f64::NAN;
        let mut last_temperature_step = f64::INFINITY;
        let mut iteration = 0;
        loop {
            poll(cancelled)?;
            let lo = lower.temperature.value();
            let hi = upper.temperature.value();
            let low_residual = lower_value - requested;
            let high_residual = upper_value - requested;
            let best_is_lower = low_residual.abs() <= high_residual.abs();
            let (best_t, best_value, best_residual) = if best_is_lower {
                (lo, lower_value, low_residual)
            } else {
                (hi, upper_value, high_residual)
            };
            if best_residual.abs() <= controls.residual_tolerance
                && (best_residual == 0.0
                    || hi - lo <= controls.temperature_tolerance_k
                    || (best_t == last_candidate_temperature
                        && last_temperature_step <= controls.temperature_tolerance_k))
            {
                let sample = if best_is_lower { lower } else { upper };
                return Ok(CaloricRoot {
                    sample,
                    coordinate: best_value,
                    work,
                });
            }
            // Strict-interior secants retain fast smooth-branch convergence;
            // every fourth step bisects to prevent a stale endpoint stalling.
            let fraction = (requested - lower_value) / (upper_value - lower_value);
            let secant = lo + fraction * (hi - lo);
            let t = if iteration % 4 != 3 && secant.is_finite() && secant > lo && secant < hi {
                secant
            } else {
                lo + 0.5 * (hi - lo)
            };
            if t <= lo || t >= hi {
                return Err(E::UnresolvedTemperature);
            }
            last_candidate_temperature = t;
            last_temperature_step = (t - best_t).abs();
            let candidate = self.bounded_caloric_sample(
                Temperature::new(t),
                pressure,
                controls.maximum_pt_evaluations,
                &mut work,
                evaluate,
                cancelled,
            )?;
            let value = self.caloric_coordinate(&candidate, target, coordinate)?;
            if value < lower_value - controls.residual_tolerance
                || value > upper_value + controls.residual_tolerance
            {
                return Err(E::Invalid(
                    "original admitted phase caloric branch is not monotone",
                ));
            }
            if value <= requested {
                lower = candidate;
                lower_value = value;
            } else {
                upper = candidate;
                upper_value = value;
            }
            iteration += 1;
        }
    }
}
