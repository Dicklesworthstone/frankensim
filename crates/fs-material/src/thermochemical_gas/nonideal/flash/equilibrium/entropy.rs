//! Entropy-constrained equilibrium at a specified pressure.
//!
//! The original screened PT owner resolves every trial, including phase amounts
//! and phase mixing entropy. The shared caloric temperature engine in `enthalpy`
//! closes the constraint without a frozen or invented mixture heat capacity.
//! An increment of zero is an isentropic endpoint constraint; a consumer supplies
//! the process assumptions needed to interpret its energy changes as work.

use super::super::enthalpy::{
    CaloricCoordinate, CaloricRoot, CaloricRootControls, EnthalpyWork,
    MAX_ENTHALPY_FLASH_EVALUATIONS, Target, spacing,
};
use super::super::{E, ThermochemicalPengRobinsonModel, poll, positive};
use super::{PhaseEnthalpyControls, ResolvedThermochemicalPhaseState};
use fs_qty::{Pressure, Temperature};
use fs_thermochem::eos::peng_robinson::{
    flash::PengRobinsonFlashOptionsV1,
    phase_screen::{PengRobinsonPhaseScreenOptionsV1, PengRobinsonPhaseScreenV1},
    stability::PengRobinsonStabilityOptionsV1,
};

/// Entropy and volume closure for reversible compression and expansion.
pub mod volume;

/// Shared hard ceiling on original PT selections in a single caloric inversion.
pub const MAX_ENTROPY_FLASH_EVALUATIONS: usize = MAX_ENTHALPY_FLASH_EVALUATIONS;

/// Physical entropy/temperature tolerances and complete nested phase work policy.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PhaseEntropyControls {
    /// Increasing source temperatures; both endpoints must be admitted.
    pub temperature_bracket: [Temperature; 2],
    /// Absolute specific-entropy residual tolerance [J/(kg K)].
    pub entropy_tolerance_j_kg_k: f64,
    /// Temperature correction or enclosing-bracket tolerance [K].
    pub temperature_tolerance_k: f64,
    /// Whole-call PT work, including the increment reference and both endpoints.
    pub maximum_pt_evaluations: usize,
    /// Original composition/fugacity closure, unchanged by the caloric root.
    pub flash_options: PengRobinsonFlashOptionsV1,
    /// Both original stationary searches must converge for every sampled state.
    pub stability: PengRobinsonStabilityOptionsV1,
}
impl PhaseEntropyControls {
    /// Default physical tolerances with an explicit original-source bracket.
    #[must_use]
    pub fn new(temperature_bracket: [Temperature; 2]) -> Self {
        let phase = PhaseEnthalpyControls::new(temperature_bracket);
        Self {
            temperature_bracket,
            entropy_tolerance_j_kg_k: 1e-8,
            temperature_tolerance_k: 1e-8,
            maximum_pt_evaluations: 64,
            flash_options: phase.enthalpy.flash_options,
            stability: phase.stability,
        }
    }

    fn root_controls(self) -> CaloricRootControls {
        CaloricRootControls {
            temperature_bracket: self.temperature_bracket,
            residual_tolerance: self.entropy_tolerance_j_kg_k,
            temperature_tolerance_k: self.temperature_tolerance_k,
            maximum_pt_evaluations: self.maximum_pt_evaluations,
        }
    }

    fn increment_root_controls(self, change: f64) -> Result<CaloricRootControls, E> {
        let mut scalar = self.root_controls();
        if change != 0.0 {
            scalar.residual_tolerance = scalar.residual_tolerance.min(0.25 * change.abs());
            if !positive(scalar.residual_tolerance) {
                return Err(E::Invalid(
                    "entropy increment cannot resolve a positive tolerance",
                ));
            }
        }
        Ok(scalar)
    }

    /// The exact bounded phase policy used for the reference and every trial.
    #[must_use]
    pub const fn phase_options(self) -> PengRobinsonPhaseScreenOptionsV1 {
        PengRobinsonPhaseScreenOptionsV1 {
            stability: self.stability,
            flash: self.flash_options,
        }
    }

    /// Admit numerical limits before any phase or source work.
    pub fn validate(self) -> Result<(), E> {
        if !self.root_controls().is_valid() {
            return Err(E::Invalid("invalid bounded entropy flash controls"));
        }
        self.phase_options().validate().map_err(E::Phase)
    }

    /// Conservative cubic-EOS ceiling, distinct from the observed PT-call count.
    pub fn maximum_eos_evaluation_bound(self) -> Result<usize, E> {
        self.validate()?;
        Ok(self.maximum_pt_evaluations
            * self
                .phase_options()
                .maximum_eos_evaluation_bound()
                .map_err(E::Phase)?)
    }
}

/// Accepted phase state with original-source entropy closure and energy changes.
#[derive(Clone, Debug, PartialEq)]
pub struct PhaseEntropySolution {
    state: ResolvedThermochemicalPhaseState,
    residual_j_kg_k: f64,
    accepted_increment_j_kg_k: Option<f64>,
    enthalpy_change_j_kg: Option<f64>,
    internal_energy_change_j_kg: Option<f64>,
    pt_evaluations: usize,
    eos_evaluation_upper_bound: usize,
}
impl PhaseEntropySolution {
    /// Complete accepted phase state and retained original source receipts.
    #[must_use]
    pub const fn state(&self) -> &ResolvedThermochemicalPhaseState {
        &self.state
    }
    /// Move the accepted state without another phase or source evaluation.
    #[must_use]
    pub fn into_state(self) -> ResolvedThermochemicalPhaseState {
        self.state
    }
    /// Evaluated entropy coordinate minus its requested target [J/(kg K)].
    #[must_use]
    pub const fn residual_j_kg_k(&self) -> f64 {
        self.residual_j_kg_k
    }
    /// Original-source Delta s, present only for the increment API.
    #[must_use]
    pub const fn accepted_increment_j_kg_k(&self) -> Option<f64> {
        self.accepted_increment_j_kg_k
    }
    /// Original-source Delta h, present only for the increment API.
    /// Formation datums cancel before rounding; no process-work model is implied.
    #[must_use]
    pub const fn enthalpy_change_j_kg(&self) -> Option<f64> {
        self.enthalpy_change_j_kg
    }
    /// Independently integrated original-source Delta u, only for the increment API.
    #[must_use]
    pub const fn internal_energy_change_j_kg(&self) -> Option<f64> {
        self.internal_energy_change_j_kg
    }
    /// Actual PT selections, including reference and bracket endpoint work.
    #[must_use]
    pub const fn pt_evaluations(&self) -> usize {
        self.pt_evaluations
    }
    /// Sum of original screen EOS work and any split-work ceilings.
    /// This is an upper bound rather than an observed EOS evaluation count.
    #[must_use]
    pub const fn eos_evaluation_upper_bound(&self) -> usize {
        self.eos_evaluation_upper_bound
    }
}

impl ThermochemicalPengRobinsonModel {
    /// Recover temperature and phase amounts from fixed pressure and full entropy.
    /// Every trial retains the original finite stability evidence; this API does
    /// not claim a global Gibbs minimum or resolve indeterminate phase amounts.
    /// Use source increments if entropy reference spacing exceeds the tolerance.
    pub fn phase_equilibrium_pressure_entropy(
        &self,
        pressure: Pressure,
        specific_entropy_j_kg_k: f64,
        controls: PhaseEntropyControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PhaseEntropySolution, E> {
        controls.validate()?;
        self.admit_caloric_flash(pressure, controls.temperature_bracket, None, &mut cancelled)?;
        if !specific_entropy_j_kg_k.is_finite()
            || spacing(specific_entropy_j_kg_k) > controls.entropy_tolerance_j_kg_k
        {
            return Err(E::Invalid(
                "absolute entropy cannot resolve tolerance; use source increments",
            ));
        }
        let target = Target::Absolute(specific_entropy_j_kg_k);
        let root = self.solve_flash_coordinate(
            pressure,
            &target,
            CaloricCoordinate::Entropy,
            controls.root_controls(),
            EnthalpyWork::default(),
            &mut |t, p, c| self.caloric_phase_sample(t, p, controls.phase_options(), c),
            &mut cancelled,
        )?;
        self.publish_phase_entropy(root, &target, &mut cancelled)
    }

    /// Change pressure and entropy while retaining every original feed amount.
    /// Delta s integrates the original NASA entropy and pressure terms, then
    /// adds the original phase-redistribution/departure difference. Independent
    /// species entropy and formation datums cancel before floating-point rounding.
    ///
    /// At zero Delta s and changed pressure this solves the reversible adiabatic
    /// endpoint, including any admitted phase transition. At unchanged pressure,
    /// zero increment, and a reference inside the bracket, one original PT call
    /// returns that exact state. Nonzero increments reserve at least three quarters
    /// of their magnitude against entropy residual error, even with a loose caller
    /// tolerance. No kinetics, heat transfer, flow rate, or mixed Cp is inferred.
    #[allow(clippy::too_many_arguments)]
    pub fn phase_equilibrium_pressure_entropy_increment(
        &self,
        initial_temperature: Temperature,
        initial_pressure: Pressure,
        final_pressure: Pressure,
        change_j_kg_k: f64,
        controls: PhaseEntropyControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PhaseEntropySolution, E> {
        controls.validate()?;
        self.admit_caloric_flash(
            final_pressure,
            controls.temperature_bracket,
            Some(initial_temperature),
            &mut cancelled,
        )?;
        if !positive(initial_pressure.value()) || !change_j_kg_k.is_finite() {
            return Err(E::Invalid(
                "finite entropy increment and positive reference pressure required",
            ));
        }
        let scalar = controls.increment_root_controls(change_j_kg_k)?;
        let mut work = EnthalpyWork::default();
        let initial = self.bounded_caloric_sample(
            initial_temperature,
            initial_pressure,
            controls.maximum_pt_evaluations,
            &mut work,
            &mut |t, p, c| self.caloric_phase_sample(t, p, controls.phase_options(), c),
            &mut cancelled,
        )?;
        let target = Target::Increment {
            initial,
            requested: change_j_kg_k,
        };
        let root = self.solve_flash_coordinate(
            final_pressure,
            &target,
            CaloricCoordinate::Entropy,
            scalar,
            work,
            &mut |t, p, c| self.caloric_phase_sample(t, p, controls.phase_options(), c),
            &mut cancelled,
        )?;
        self.publish_phase_entropy(root, &target, &mut cancelled)
    }

    fn publish_phase_entropy(
        &self,
        root: CaloricRoot<PengRobinsonPhaseScreenV1>,
        target: &Target<PengRobinsonPhaseScreenV1>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<PhaseEntropySolution, E> {
        poll(cancelled)?;
        let (increment, enthalpy_change, internal_energy_change) = match target {
            Target::Absolute(_) => (None, None, None),
            Target::Increment { initial, .. } => (
                Some(root.coordinate),
                Some(self.caloric_enthalpy_change(initial, &root.sample)?),
                Some(self.caloric_internal_energy_change(initial, &root.sample)?),
            ),
        };
        let residual = root.coordinate - target.requested();
        let state = self.publish_caloric_phase(root.sample, cancelled)?;
        poll(cancelled)?;
        Ok(PhaseEntropySolution {
            state,
            residual_j_kg_k: residual,
            accepted_increment_j_kg_k: increment,
            enthalpy_change_j_kg: enthalpy_change,
            internal_energy_change_j_kg: internal_energy_change,
            pt_evaluations: root.work.calls,
            eos_evaluation_upper_bound: root.work.eos_upper_bound,
        })
    }
}
