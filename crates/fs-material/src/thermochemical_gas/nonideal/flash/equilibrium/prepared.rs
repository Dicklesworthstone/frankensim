//! Repeated original-source PH increments with one immutable inlet selection.
//!
//! A prepared path borrows its original mixture and fixes the inlet, source
//! bracket and phase policy. It retains a default outlet pressure; individual
//! queries may supply another positive pressure without preparing the inlet
//! again. Every endpoint uses the native caloric root and source-integrated
//! enthalpy/entropy increments.

use super::{
    CaloricSample, E, EnthalpyWork, PengRobinsonPhaseScreenV1, PhaseEnthalpyControls,
    PhaseEnthalpySolution, Target, ThermochemicalPengRobinsonModel, poll, positive,
};
use fs_qty::{Pressure, Temperature};

/// Immutable original inlet and policy for repeated PH queries.
///
/// Construction performs one PT selection. Charge its preparation counters once
/// to a surrounding operation; each returned solution counts only new work.
/// No state can be reused with a different mixture, phase policy or source
/// bracket. The original finite stability evidence is retained, not strengthened.
#[derive(Clone, Debug)]
pub struct PreparedPhaseEnthalpyPath<'m> {
    model: &'m ThermochemicalPengRobinsonModel,
    initial: CaloricSample<PengRobinsonPhaseScreenV1>,
    final_pressure: Pressure,
    controls: PhaseEnthalpyControls,
    preparation_pt_evaluations: usize,
    preparation_eos_evaluation_upper_bound: usize,
}

impl ThermochemicalPengRobinsonModel {
    /// Select one original inlet for repeated source-integrated PH increments.
    ///
    /// The supplied controls fix the caloric bracket, tolerances and phase
    /// policy. Their PT limit applies separately to each subsequent query;
    /// preparation work is reported separately so callers can budget the whole
    /// operation without charging a cached inlet repeatedly.
    pub fn prepare_phase_enthalpy_path(
        &self,
        initial_temperature: Temperature,
        initial_pressure: Pressure,
        final_pressure: Pressure,
        controls: PhaseEnthalpyControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PreparedPhaseEnthalpyPath<'_>, E> {
        controls.validate()?;
        self.admit_enthalpy_flash(
            final_pressure,
            controls.enthalpy,
            Some(initial_temperature),
            &mut cancelled,
        )?;
        if !positive(initial_pressure.value()) {
            return Err(E::Invalid(
                "positive prepared PH reference pressure required",
            ));
        }
        let mut work = EnthalpyWork::default();
        let initial = self.enthalpy_flash_sample(
            initial_temperature,
            initial_pressure,
            controls.enthalpy,
            &mut work,
            &mut |t, p, c| self.caloric_phase_sample(t, p, controls.phase_options(), c),
            &mut cancelled,
        )?;
        poll(&mut cancelled)?;
        Ok(PreparedPhaseEnthalpyPath {
            model: self,
            initial,
            final_pressure,
            controls,
            preparation_pt_evaluations: work.calls,
            preparation_eos_evaluation_upper_bound: work.eos_upper_bound,
        })
    }
}

impl PreparedPhaseEnthalpyPath<'_> {
    /// PT selections performed when preparing this immutable inlet.
    #[must_use]
    pub const fn preparation_pt_evaluations(&self) -> usize {
        self.preparation_pt_evaluations
    }

    /// Preparation's original EOS work ceiling, not a measured EOS count.
    #[must_use]
    pub const fn preparation_eos_evaluation_upper_bound(&self) -> usize {
        self.preparation_eos_evaluation_upper_bound
    }

    /// Resolve a source-integrated increment from the same original inlet.
    ///
    /// The allowance bounds new PT selections and is capped by the prepared
    /// controls. Returned counters exclude preparation. Zero remaining allowance
    /// still permits exact zero heat at unchanged pressure when the inlet lies
    /// inside the bracket: that path performs no new PT selection.
    ///
    /// Failure or cancellation leaves this immutable path reusable. No caloric
    /// chart, heat-capacity approximation or absolute-energy subtraction is used.
    pub fn solve_increment(
        &self,
        change_j_kg: f64,
        remaining_pt_evaluations: usize,
        cancelled: impl FnMut() -> bool,
    ) -> Result<PhaseEnthalpySolution, E> {
        self.solve_increment_at_pressure(
            self.final_pressure, change_j_kg, remaining_pt_evaluations, cancelled,
        )
    }

    /// Resolve an original-inlet enthalpy increment at a new positive pressure.
    ///
    /// Pressure changes do not replace or chain the reference state. This lets
    /// Rayleigh/Hugoniot and other pressure paths reuse one original phase
    /// selection. The immutable caloric interval was admitted at preparation;
    /// the native PT owner admits every new pressure/temperature endpoint.
    /// Counters and the remaining-work allowance cover only new PT work.
    ///
    /// Zero heat still requires a solve if pressure differs from the inlet.
    /// Neither a process entropy law nor pressure-flow work is added here.
    pub fn solve_increment_at_pressure(
        &self,
        final_pressure: Pressure,
        change_j_kg: f64,
        remaining_pt_evaluations: usize,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PhaseEnthalpySolution, E> {
        poll(&mut cancelled)?;
        if !positive(final_pressure.value()) || !change_j_kg.is_finite() {
            return Err(E::Invalid(
                "positive prepared PH pressure and finite enthalpy increment required",
            ));
        }
        let temperature = self.initial.temperature.value();
        let bracket = self.controls.enthalpy.temperature_bracket;
        let exact_rest = change_j_kg == 0.0
            && self.initial.pressure == final_pressure
            && temperature >= bracket[0].value()
            && temperature <= bracket[1].value();
        if remaining_pt_evaluations == 0 && !exact_rest {
            return Err(E::WorkLimit {
                flash_evaluations: 0,
            });
        }
        let mut controls = self.controls.enthalpy;
        // Keep valid scalar controls for the zero-query shortcut. Every other
        // call receives at most the actual remaining external work allowance.
        controls.maximum_flash_evaluations = controls
            .maximum_flash_evaluations
            .min(remaining_pt_evaluations.max(1));
        let outcome = self.model.solve_flash_enthalpy(
            final_pressure,
            Target::Increment {
                initial: self.initial.clone(),
                requested: change_j_kg,
            },
            controls,
            EnthalpyWork::default(),
            &mut |t, p, c| {
                self.model
                    .caloric_phase_sample(t, p, self.controls.phase_options(), c)
            },
            &mut cancelled,
        )?;
        self.model.publish_phase_enthalpy(outcome, &mut cancelled)
    }
}
