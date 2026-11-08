//! Repeated original-source PS increments anchored to one immutable inlet.
//!
//! Pressure may vary between queries. Every increment retains the same original
//! inventory, source datums, phase policy and temperature bracket.

use super::super::CaloricSample;
use super::{
    CaloricCoordinate, E, EnthalpyWork, PengRobinsonPhaseScreenV1, PhaseEntropyControls,
    PhaseEntropySolution, Target, ThermochemicalPengRobinsonModel, poll,
};
use fs_qty::{Pressure, Temperature};

/// One prepared original inlet for pressure/entropy paths and isentropic devices.
///
/// Preparation performs one original PT selection. Each returned solution counts
/// only subsequent work; surrounding devices charge preparation once. This is
/// an immutable source reference, not an interpolated caloric surrogate.
#[derive(Clone, Debug)]
pub struct PreparedPhaseEntropyPath<'m> {
    model: &'m ThermochemicalPengRobinsonModel,
    initial: CaloricSample<PengRobinsonPhaseScreenV1>,
    controls: PhaseEntropyControls,
    preparation_pt_evaluations: usize,
    preparation_eos_evaluation_upper_bound: usize,
}

impl ThermochemicalPengRobinsonModel {
    /// Prepare one inlet for repeated original entropy increments at varying p.
    pub fn prepare_phase_entropy_path(
        &self,
        initial_temperature: Temperature,
        initial_pressure: Pressure,
        controls: PhaseEntropyControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PreparedPhaseEntropyPath<'_>, E> {
        controls.validate()?;
        self.admit_caloric_flash(
            initial_pressure,
            controls.temperature_bracket,
            Some(initial_temperature),
            &mut cancelled,
        )?;
        let mut work = EnthalpyWork::default();
        let initial = self.bounded_caloric_sample(
            initial_temperature,
            initial_pressure,
            controls.maximum_pt_evaluations,
            &mut work,
            &mut |t, p, c| self.caloric_phase_sample(t, p, controls.phase_options(), c),
            &mut cancelled,
        )?;
        poll(&mut cancelled)?;
        Ok(PreparedPhaseEntropyPath {
            model: self,
            initial,
            controls,
            preparation_pt_evaluations: work.calls,
            preparation_eos_evaluation_upper_bound: work.eos_upper_bound,
        })
    }
}

impl PreparedPhaseEntropyPath<'_> {
    /// Original PT selections used to prepare the immutable inlet.
    #[must_use]
    pub const fn preparation_pt_evaluations(&self) -> usize {
        self.preparation_pt_evaluations
    }

    /// Original preparation EOS ceiling; not an observed EOS-call count.
    #[must_use]
    pub const fn preparation_eos_evaluation_upper_bound(&self) -> usize {
        self.preparation_eos_evaluation_upper_bound
    }

    /// Recover an endpoint from the original inlet, never from the previous query.
    ///
    /// Zero entropy increment traces the original equilibrium isentrope. The
    /// pressure must be positive; the fixed source bracket and phase policy are
    /// admitted by the native owner. All energy differences cancel species
    /// reference datums before rounding. Returned PT counts exclude preparation.
    ///
    /// An exact zero increment at the inlet pressure needs no additional PT
    /// selection if the inlet lies in the bracket, even with zero remaining work.
    /// Refusal and cancellation leave this path reusable.
    pub fn solve_increment(
        &self,
        pressure: Pressure,
        change_j_kg_k: f64,
        remaining_pt_evaluations: usize,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PhaseEntropySolution, E> {
        poll(&mut cancelled)?;
        if !change_j_kg_k.is_finite() {
            return Err(E::Invalid("finite prepared PS entropy increment required"));
        }
        self.model.admit_caloric_flash(
            pressure,
            self.controls.temperature_bracket,
            Some(self.initial.temperature),
            &mut cancelled,
        )?;
        let t = self.initial.temperature.value();
        let [lo, hi] = self.controls.temperature_bracket;
        let exact_rest = change_j_kg_k == 0.0
            && pressure == self.initial.pressure
            && t >= lo.value()
            && t <= hi.value();
        if remaining_pt_evaluations == 0 && !exact_rest {
            return Err(E::WorkLimit {
                flash_evaluations: 0,
            });
        }
        let mut scalar = self.controls.increment_root_controls(change_j_kg_k)?;
        scalar.maximum_pt_evaluations = scalar
            .maximum_pt_evaluations
            .min(remaining_pt_evaluations.max(1));
        let target = Target::Increment {
            initial: self.initial.clone(),
            requested: change_j_kg_k,
        };
        let root = self.model.solve_flash_coordinate(
            pressure,
            &target,
            CaloricCoordinate::Entropy,
            scalar,
            EnthalpyWork::default(),
            &mut |t, p, c| {
                self.model
                    .caloric_phase_sample(t, p, self.controls.phase_options(), c)
            },
            &mut cancelled,
        )?;
        self.model
            .publish_phase_entropy(root, &target, &mut cancelled)
    }
}
