//! Source-bound enthalpy paths across screened homogeneous and two-phase states.
//!
//! Both original tangent-plane searches must converge. Their finite evidence is
//! retained without claiming global stability; a witnessed instability must
//! resolve the original split. The same caloric and scalar enthalpy owners serve
//! this API and the stricter two-phase API. No mixed Cp or transport is invented.

use super::enthalpy::{EnthalpyOutcome, EnthalpyWork, FlashEnthalpyControls, Target, spacing};
use super::{
    CaloricSample, E, FrozenIdealGasMixtureEvaluationV1, PengRobinsonMixtureV1,
    ThermochemicalPengRobinsonModel, failure, poll, positive, residual_entropy,
};
use fs_qty::{Pressure, Temperature};
use fs_thermochem::eos::peng_robinson::{
    flash::PengRobinsonFlashErrorV1,
    phase_screen::{PengRobinsonPhaseScreenOptionsV1, PengRobinsonPhaseScreenV1},
    stability::PengRobinsonStabilityOptionsV1,
};

/// Conserved internal-energy and volume closure over the same phase owner.
pub mod volume;

/// Energy-root controls plus the explicit, bounded phase stability screen.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PhaseEnthalpyControls {
    /// Shared original-source interval, residual tolerances and PT work cap.
    pub enthalpy: FlashEnthalpyControls,
    /// Both stationary searches must converge for each PT state.
    pub stability: PengRobinsonStabilityOptionsV1,
}
impl PhaseEnthalpyControls {
    /// A source temperature bracket that may span liquid, coexistence and vapor.
    #[must_use]
    pub fn new(temperature_bracket: [Temperature; 2]) -> Self {
        Self {
            enthalpy: FlashEnthalpyControls::new(temperature_bracket),
            stability: PengRobinsonStabilityOptionsV1 {
                stationarity_tolerance: 1e-12,
                tangent_plane_tolerance: 1e-12,
                ..Default::default()
            },
        }
    }
    /// The exact phase policy applied at every sampled temperature.
    #[must_use]
    pub const fn phase_options(self) -> PengRobinsonPhaseScreenOptionsV1 {
        PengRobinsonPhaseScreenOptionsV1 {
            stability: self.stability,
            flash: self.enthalpy.flash_options,
        }
    }
    /// Validate the scalar and both nested phase policies before doing work.
    pub fn validate(self) -> Result<(), E> {
        self.enthalpy.validate()?;
        self.phase_options().validate().map_err(E::Phase)
    }
    /// Conservative total cubic EOS ceiling, distinct from observed PT calls.
    pub fn maximum_eos_evaluation_bound(self) -> Result<usize, E> {
        self.validate()?;
        Ok(self.enthalpy.maximum_flash_evaluations
            * self
                .phase_options()
                .maximum_eos_evaluation_bound()
                .map_err(E::Phase)?)
    }
}

/// An admitted caloric state with the exact source law and finite phase evidence.
/// Homogeneous states have no inferred vapor amount or liquid/vapor phase label.
#[derive(Clone, Debug, PartialEq)]
pub struct ResolvedThermochemicalPhaseState {
    law: PengRobinsonMixtureV1,
    reference: FrozenIdealGasMixtureEvaluationV1,
    sample: CaloricSample<PengRobinsonPhaseScreenV1>,
}
impl ResolvedThermochemicalPhaseState {
    /// Exact retained source critical parameters, interactions and feed inventory.
    #[must_use]
    pub const fn constitutive_law(&self) -> &PengRobinsonMixtureV1 {
        &self.law
    }
    /// Complete original NASA feed evaluation and reference convention receipts.
    #[must_use]
    pub const fn reference(&self) -> &FrozenIdealGasMixtureEvaluationV1 {
        &self.reference
    }
    /// Screened feed root or original separated phases, with both search receipts.
    #[must_use]
    pub const fn phase(&self) -> &PengRobinsonPhaseScreenV1 {
        &self.sample.phase
    }
    /// Absolute temperature [K].
    #[must_use]
    pub fn temperature_k(&self) -> f64 {
        self.sample.temperature.value()
    }
    /// Common pressure [Pa].
    #[must_use]
    pub fn pressure_pa(&self) -> f64 {
        self.sample.pressure.value()
    }
    /// Original or additive phase molar volume per exact feed mole [m³/mol].
    #[must_use]
    pub const fn molar_volume_m3_mol(&self) -> f64 {
        self.sample.molar_volume
    }
    /// Exact retained feed mass divided by its volume [kg/m³].
    #[must_use]
    pub fn density_kg_m3(&self) -> f64 {
        self.law.molar_mass().value() / self.sample.molar_volume
    }
    /// Liquid mass amount, present only for an admitted split.
    #[must_use]
    pub fn liquid_mass_fraction(&self) -> Option<f64> {
        self.phase()
            .two_phase()
            .map(|_| self.sample.liquid_mass_fraction)
    }
    /// Vapor mass amount, present only for an admitted split.
    #[must_use]
    pub fn vapor_mass_fraction(&self) -> Option<f64> {
        self.phase()
            .two_phase()
            .map(|_| self.sample.vapor_mass_fraction)
    }
    /// Vapor volume amount, present only for an admitted split.
    #[must_use]
    pub fn vapor_volume_fraction(&self) -> Option<f64> {
        self.phase()
            .two_phase()
            .map(|_| self.sample.vapor_volume_fraction)
    }
    /// Full source-reference enthalpy [J/kg].
    #[must_use]
    pub const fn specific_enthalpy_j_kg(&self) -> f64 {
        self.sample.specific_h
    }
    /// Full source-reference internal energy [J/kg].
    #[must_use]
    pub const fn specific_internal_energy_j_kg(&self) -> f64 {
        self.sample.specific_u
    }
    /// Full entropy with the actual phase mixing correction [J/(kg K)].
    #[must_use]
    pub const fn specific_entropy_j_kg_k(&self) -> f64 {
        self.sample.specific_s
    }
}

/// Accepted homogeneous or separated endpoint of one conserved-enthalpy solve.
#[derive(Clone, Debug, PartialEq)]
pub struct PhaseEnthalpySolution {
    state: ResolvedThermochemicalPhaseState,
    residual_j_kg: f64,
    accepted_increment_j_kg: Option<f64>,
    entropy_change_j_kg_k: Option<f64>,
    pt_evaluations: usize,
    eos_evaluation_upper_bound: usize,
}
impl PhaseEnthalpySolution {
    /// Complete accepted phase state and all original source receipts.
    #[must_use]
    pub const fn state(&self) -> &ResolvedThermochemicalPhaseState {
        &self.state
    }
    /// Move the state without another source or phase evaluation.
    #[must_use]
    pub fn into_state(self) -> ResolvedThermochemicalPhaseState {
        self.state
    }
    /// Independently evaluated enthalpy coordinate minus requested target [J/kg].
    #[must_use]
    pub const fn residual_j_kg(&self) -> f64 {
        self.residual_j_kg
    }
    /// Source-integrated actual enthalpy increment, only for the increment API.
    #[must_use]
    pub const fn accepted_increment_j_kg(&self) -> Option<f64> {
        self.accepted_increment_j_kg
    }
    /// Datum-invariant actual entropy change, only for the increment API.
    #[must_use]
    pub const fn entropy_change_j_kg_k(&self) -> Option<f64> {
        self.entropy_change_j_kg_k
    }
    /// Actual attempted PT selections, including the reference and endpoints.
    #[must_use]
    pub const fn pt_evaluations(&self) -> usize {
        self.pt_evaluations
    }
    /// Sum of observed screen EOS work plus optional split work ceilings.
    /// This is an upper bound, not an observed total EOS count.
    #[must_use]
    pub const fn eos_evaluation_upper_bound(&self) -> usize {
        self.eos_evaluation_upper_bound
    }
}

impl ThermochemicalPengRobinsonModel {
    /// Resolve homogeneous or separated PT calorics using the original two-search
    /// screen. No absence-of-split fallback or global-minimum claim is made.
    pub fn evaluate_phase_equilibrium(
        &self,
        temperature: Temperature,
        pressure: Pressure,
        options: PengRobinsonPhaseScreenOptionsV1,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<ResolvedThermochemicalPhaseState, E> {
        let sample = self.caloric_phase_sample(temperature, pressure, options, &mut cancelled)?;
        self.publish_caloric_phase(sample, &mut cancelled)
    }

    /// Recover temperature and phase amount from fixed pressure and full enthalpy.
    /// The bracket may cross the mixture bubble/dew boundaries. Every sampled
    /// state must pass the original phase and caloric admissions; unresolved or
    /// amount-indeterminate states refuse. This does not find a global minimum.
    pub fn phase_equilibrium_pressure_enthalpy(
        &self,
        pressure: Pressure,
        specific_enthalpy_j_kg: f64,
        controls: PhaseEnthalpyControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PhaseEnthalpySolution, E> {
        controls.validate()?;
        self.admit_enthalpy_flash(pressure, controls.enthalpy, None, &mut cancelled)?;
        if !specific_enthalpy_j_kg.is_finite()
            || spacing(specific_enthalpy_j_kg) > controls.enthalpy.enthalpy_tolerance_j_kg
        {
            return Err(E::Invalid(
                "absolute enthalpy cannot resolve tolerance; use source increments",
            ));
        }
        let outcome = self.solve_flash_enthalpy(
            pressure,
            Target::Absolute(specific_enthalpy_j_kg),
            controls.enthalpy,
            EnthalpyWork::default(),
            &mut |t, p, c| self.caloric_phase_sample(t, p, controls.phase_options(), c),
            &mut cancelled,
        )?;
        self.publish_phase_enthalpy(outcome, &mut cancelled)
    }

    /// Change enthalpy and/or pressure while retaining the exact feed inventory.
    /// Ideal source increments are integrated before EOS departure differences;
    /// even large independent species formation/entropy datums cancel before
    /// rounding. Heating or cooling may cross both mixture phase boundaries.
    /// This is a caloric constraint, not a flow or process-entropy model.
    /// At unchanged pressure and zero increment, a reference inside the bracket
    /// returns exactly after one PT selection.
    #[allow(clippy::too_many_arguments)]
    pub fn phase_equilibrium_pressure_enthalpy_increment(
        &self,
        initial_temperature: Temperature,
        initial_pressure: Pressure,
        final_pressure: Pressure,
        change_j_kg: f64,
        controls: PhaseEnthalpyControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PhaseEnthalpySolution, E> {
        controls.validate()?;
        self.admit_enthalpy_flash(
            final_pressure,
            controls.enthalpy,
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
            controls.enthalpy,
            &mut work,
            &mut |t, p, c| self.caloric_phase_sample(t, p, controls.phase_options(), c),
            &mut cancelled,
        )?;
        let outcome = self.solve_flash_enthalpy(
            final_pressure,
            Target::Increment {
                initial,
                requested: change_j_kg,
            },
            controls.enthalpy,
            work,
            &mut |t, p, c| self.caloric_phase_sample(t, p, controls.phase_options(), c),
            &mut cancelled,
        )?;
        self.publish_phase_enthalpy(outcome, &mut cancelled)
    }

    fn caloric_phase_sample(
        &self,
        temperature: Temperature,
        pressure: Pressure,
        options: PengRobinsonPhaseScreenOptionsV1,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<CaloricSample<PengRobinsonPhaseScreenV1>, E> {
        poll(cancelled)?;
        let reference = self
            .source()
            .thermodynamics()
            .evaluate_properties(temperature, pressure)
            .map_err(|error| failure("screened phase feed thermodynamics", error))?;
        let phase = self
            .mechanics()
            .phase_screen_pt_cancellable(temperature, pressure, options, &mut *cancelled)
            .map_err(|error| match error {
                PengRobinsonFlashErrorV1::Cancelled => E::Cancelled,
                other => E::Phase(other),
            })?;
        let work = phase.eos_evaluation_upper_bound();
        if let Some(split) = phase.two_phase() {
            let values =
                self.two_phase_caloric_sample(temperature, pressure, split, reference, cancelled)?;
            return Ok(values.with_phase(phase, work));
        }
        let root = phase.stability().feed();
        let mass = self.mechanics().molar_mass().value();
        if root
            .caloric_mechanics(reference.cv().value(), mass)
            .is_none()
        {
            return Err(E::Invalid(
                "screened homogeneous phase requires stable original caloric mechanics",
            ));
        }
        let molar_volume = root.molar_volume_m3_per_mol();
        let residual_h = root.residual_enthalpy_j_mol();
        let residual_u = root.residual_internal_energy_j_mol();
        let entropy_correction =
            residual_entropy(root, self.mechanics().mole_fractions(), temperature.value());
        let specific_h = (reference.h().value() + residual_h) / mass;
        let specific_u = (reference.u().value() + residual_u) / mass;
        let specific_s = (reference.s().value() + entropy_correction) / mass;
        if [
            residual_h,
            residual_u,
            entropy_correction,
            specific_h,
            specific_u,
            specific_s,
        ]
        .into_iter()
        .any(|value| !value.is_finite())
            || !positive(molar_volume)
            || !positive(mass / molar_volume)
        {
            return Err(E::Invalid(
                "screened homogeneous calorics are unrepresentable",
            ));
        }
        poll(cancelled)?;
        Ok(CaloricSample {
            phase,
            temperature,
            pressure,
            eos_evaluation_upper_bound: work,
            molar_volume,
            // Private placeholders are never exposed as homogeneous phase amounts.
            liquid_mass_fraction: 0.0,
            vapor_mass_fraction: 0.0,
            vapor_volume_fraction: 0.0,
            residual_h,
            residual_u,
            entropy_correction,
            specific_h,
            specific_u,
            specific_s,
        })
    }

    fn publish_caloric_phase(
        &self,
        sample: CaloricSample<PengRobinsonPhaseScreenV1>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<ResolvedThermochemicalPhaseState, E> {
        poll(cancelled)?;
        let reference = self
            .source()
            .thermodynamics()
            .evaluate(sample.temperature, sample.pressure)
            .map_err(|error| failure("accepted screened phase source receipt", error))?;
        let state = ResolvedThermochemicalPhaseState {
            law: self.mechanics().clone(),
            reference,
            sample,
        };
        poll(cancelled)?;
        Ok(state)
    }

    fn publish_phase_enthalpy(
        &self,
        outcome: EnthalpyOutcome<PengRobinsonPhaseScreenV1>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<PhaseEnthalpySolution, E> {
        let state = self.publish_caloric_phase(outcome.sample, cancelled)?;
        poll(cancelled)?;
        Ok(PhaseEnthalpySolution {
            state,
            residual_j_kg: outcome.residual_j_kg,
            accepted_increment_j_kg: outcome.accepted_increment_j_kg,
            entropy_change_j_kg_k: outcome.entropy_change_j_kg_k,
            pt_evaluations: outcome.flash_evaluations,
            eos_evaluation_upper_bound: outcome.eos_evaluation_upper_bound,
        })
    }
}
