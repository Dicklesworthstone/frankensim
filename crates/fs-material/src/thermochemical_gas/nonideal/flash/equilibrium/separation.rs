//! Reusable original liquid/vapor sources from an already admitted phase split.
//!
//! Routing equilibrated phases does not add latent heat or remove ideal mixing
//! entropy a second time. Species references cancel before the original phase
//! departures are compared; actual species amounts remain a separate audit.

use super::super::{
    ThermochemicalFlashError, ThermochemicalPengRobinsonModel, failure, positive, residual_entropy,
    x_log_x,
};
use super::ResolvedThermochemicalPhaseState;
use crate::thermochemical_gas::energy::GasEnergyError;
use crate::thermochemical_gas::nonideal::energy::parcel::transfer::mixing::model_for_source;
use fs_math::eft::two_sum;
use fs_qty::{Amount, Pressure, Temperature};
use fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R;
use fs_thermochem::eos::peng_robinson::{PengRobinsonPhaseRootV1, PengRobinsonPhaseV1};
use fs_thermochem::mixture::FrozenIdealGasMixtureEvaluationV1;

/// The original source or phase could not supply both physical outlets.
#[derive(Clone, Debug, PartialEq)]
pub enum PhaseSeparationError {
    /// Invalid input or unrepresentable phase calorics.
    Invalid(&'static str),
    /// The supplied state has no admitted nontrivial liquid/vapor split.
    NoSeparatedPhases,
    /// Original species, caloric, EOS or source-identity refusal.
    Source(ThermochemicalFlashError),
    /// The original source-preserving composition/catalog join refused.
    Catalog(GasEnergyError),
    /// Cancellation before publishing either outlet.
    Cancelled,
}
impl core::fmt::Display for PhaseSeparationError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "original equilibrium phase separation: {self:?}")
    }
}
impl std::error::Error for PhaseSeparationError {}
type Error = PhaseSeparationError;

/// One original homogeneous phase with its own reusable source composition.
///
/// This is a caloric state, not a dilute-gas transport evaluation. Its branch
/// comes from the admitted parent split; no new global-stability claim is made.
#[derive(Clone, Debug, PartialEq)]
pub struct SeparatedMixturePhase {
    model: ThermochemicalPengRobinsonModel,
    reference: FrozenIdealGasMixtureEvaluationV1,
    root: PengRobinsonPhaseRootV1,
    phase: PengRobinsonPhaseV1,
    specific_enthalpy_j_kg: f64,
    specific_internal_energy_j_kg: f64,
    specific_entropy_j_kg_k: f64,
}
impl SeparatedMixturePhase {
    /// Original cards, species and PR pair data bound to this outlet composition.
    pub const fn model(&self) -> &ThermochemicalPengRobinsonModel {
        &self.model
    }
    /// Actual original NASA source receipts for this phase composition.
    pub const fn reference(&self) -> &FrozenIdealGasMixtureEvaluationV1 {
        &self.reference
    }
    /// Original EOS root, recomputed after canonical composition closure.
    pub const fn root(&self) -> &PengRobinsonPhaseRootV1 {
        &self.root
    }
    /// Explicit liquid or vapor branch inherited from the original split.
    pub const fn phase(&self) -> PengRobinsonPhaseV1 {
        self.phase
    }
    /// Unchanged separator temperature [K].
    pub fn temperature_k(&self) -> f64 {
        self.reference.temperature().value()
    }
    /// Unchanged separator pressure [Pa].
    pub fn pressure_pa(&self) -> f64 {
        self.reference.pressure().value()
    }
    /// Actual outlet density, using its own phase molar mass [kg/m³].
    pub fn density_kg_m3(&self) -> f64 {
        self.root.mass_density_kg_per_m3()
    }
    /// Original phase enthalpy, including pressure-flow work once [J/kg].
    pub const fn specific_enthalpy_j_kg(&self) -> f64 {
        self.specific_enthalpy_j_kg
    }
    /// Original phase internal energy [J/kg].
    pub const fn specific_internal_energy_j_kg(&self) -> f64 {
        self.specific_internal_energy_j_kg
    }
    /// Original phase entropy, including its actual mixing entropy [J/(kg K)].
    pub const fn specific_entropy_j_kg_k(&self) -> f64 {
        self.specific_entropy_j_kg_k
    }
}

/// Both reusable phase sources and original, reference-independent accounts.
#[derive(Clone, Debug, PartialEq)]
pub struct PhaseSeparationSolution {
    /// Unchanged original inlet, including its finite phase-screen evidence.
    pub inlet: ResolvedThermochemicalPhaseState,
    /// Liquid and vapor outlets, in that order throughout this report.
    pub outlets: [SeparatedMixturePhase; 2],
    /// Original molar amounts per feed mole: [1-beta, beta].
    pub mole_fractions: [f64; 2],
    /// Each phase mole amount times its own molar mass divided by feed molar mass.
    /// Numerical total/species mass defects are audited by the flow consumer.
    pub mass_fractions: [f64; 2],
    /// Outgoing-minus-incoming enthalpy per feed mass, after species references
    /// cancel on the conserved ledger [J/kg]. No correction heat is applied.
    pub enthalpy_residual_j_kg: f64,
    /// Outgoing-minus-incoming entropy per feed mass [J/(kg K)], unclipped.
    /// The inlet already includes the original phase mixing contribution.
    pub entropy_residual_j_kg_k: f64,
    /// Two bounded cubic EOS evaluations after the retained inlet selection.
    pub additional_eos_evaluations: usize,
}

impl ThermochemicalPengRobinsonModel {
    /// Separate an admitted inlet without repeating its stability/flash solve.
    ///
    /// The original mechanical law AND complete caloric receipts must match.
    /// Source composition closure retains every positive species and its cards;
    /// each outlet then uses its own original PR root and NASA calorics. No gas
    /// transport data are required or extrapolated, including for a liquid.
    ///
    /// This is routing at unchanged T,p, without pressure loss, heat, shaft work,
    /// reaction or phase kinetics. A homogeneous or amount-indeterminate inlet
    /// cannot be assigned a guessed vapor quality and refuses. Physical flow
    /// tolerances and actual total/species mass checks belong to the consumer.
    /// Absolute h/s retain source datums; reference-safe residuals assume the
    /// intended species ledger, with its floating-point defects audited separately.
    pub fn separate_phase_equilibrium(
        &self,
        inlet: ResolvedThermochemicalPhaseState,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PhaseSeparationSolution, Error> {
        poll(&mut cancelled)?;
        let split = inlet.phase().two_phase().ok_or(Error::NoSeparatedPhases)?;
        if inlet.constitutive_law() != self.mechanics() {
            return Err(Error::Invalid(
                "separator source differs from the admitted inlet",
            ));
        }
        let temperature = Temperature::new(inlet.temperature_k());
        let pressure = Pressure::new(inlet.pressure_pa());
        let original = self
            .source()
            .thermodynamics()
            .evaluate(temperature, pressure)
            .map_err(|e| source_error("separator inlet source identity", e))?;
        if inlet.reference() != &original {
            return Err(Error::Invalid(
                "separator caloric source differs from the admitted inlet",
            ));
        }
        let mole_fractions = [1.0 - split.vapor_fraction_mol(), split.vapor_fraction_mol()];
        let outlets = [
            self.separated_outlet(
                split.liquid_mole_fractions(),
                temperature,
                pressure,
                PengRobinsonPhaseV1::Liquid,
                &mut cancelled,
            )?,
            self.separated_outlet(
                split.vapor_mole_fractions(),
                temperature,
                pressure,
                PengRobinsonPhaseV1::Vapor,
                &mut cancelled,
            )?,
        ];
        let mass = self.mechanics().molar_mass().value();
        let mass_fractions = std::array::from_fn(|i| {
            mole_fractions[i] * (outlets[i].model.mechanics().molar_mass().value() / mass)
        });
        // Every species is at the SAME temperature and pressure before/after
        // routing. Its ideal h/s reference cancels before phase sums are formed.
        // Recomputed PR departures and phase mixing retain closure roundoff.
        let residual_h = sum([
            mole_fractions[0] * outlets[0].root.residual_enthalpy_j_mol(),
            mole_fractions[1] * outlets[1].root.residual_enthalpy_j_mol(),
            -inlet.sample.residual_h,
        ]) / mass;
        let outlet_entropy = std::array::from_fn::<_, 2, _>(|i| {
            let x = outlets[i].model.mechanics().mole_fractions();
            residual_entropy(&outlets[i].root, x, temperature.value()) - R * x_log_x(x)
        });
        let residual_s = sum([
            mole_fractions[0] * outlet_entropy[0],
            mole_fractions[1] * outlet_entropy[1],
            R * x_log_x(self.mechanics().mole_fractions()),
            -inlet.sample.entropy_correction,
        ]) / mass;
        if !residual_h.is_finite()
            || !residual_s.is_finite()
            || mass_fractions.iter().any(|x| !positive(*x) || *x >= 1.0)
        {
            return Err(Error::Invalid(
                "separated phase amounts or calorics are unrepresentable",
            ));
        }
        poll(&mut cancelled)?;
        Ok(PhaseSeparationSolution {
            inlet,
            outlets,
            mole_fractions,
            mass_fractions,
            enthalpy_residual_j_kg: residual_h,
            entropy_residual_j_kg_k: residual_s,
            additional_eos_evaluations: 2,
        })
    }

    fn separated_outlet(
        &self,
        fractions: &[f64],
        temperature: Temperature,
        pressure: Pressure,
        phase: PengRobinsonPhaseV1,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<SeparatedMixturePhase, Error> {
        poll(cancelled)?;
        let amounts: Vec<_> = fractions.iter().map(|&x| Amount::new(x)).collect();
        let source = self
            .source()
            .with_mole_amounts(&amounts)
            .map_err(|e| Error::Source(e.into()))?;
        let model = model_for_source(&[self], source).map_err(Error::Catalog)?;
        let reference = model
            .source()
            .thermodynamics()
            .evaluate(temperature, pressure)
            .map_err(|e| source_error("separated phase calorics", e))?;
        poll(cancelled)?;
        let root = model
            .mechanics()
            .evaluate_pt(temperature, pressure)
            .map_err(|e| source_error("separated phase EOS", e))?
            .select(phase)
            .clone();
        let properties = reference.properties();
        let mass = model.mechanics().molar_mass().value();
        let h = (properties.h().value() + root.residual_enthalpy_j_mol()) / mass;
        let u = (properties.u().value() + root.residual_internal_energy_j_mol()) / mass;
        let s = (properties.s().value()
            + residual_entropy(
                &root,
                model.mechanics().mole_fractions(),
                temperature.value(),
            ))
            / mass;
        if root
            .caloric_mechanics(properties.cv().value(), mass)
            .is_none()
            || !positive(root.mass_density_kg_per_m3())
            || [h, u, s].iter().any(|x| !x.is_finite())
        {
            return Err(Error::Invalid(
                "separated phase needs stable finite original calorics",
            ));
        }
        poll(cancelled)?;
        Ok(SeparatedMixturePhase {
            model,
            reference,
            root,
            phase,
            specific_enthalpy_j_kg: h,
            specific_internal_energy_j_kg: u,
            specific_entropy_j_kg_k: s,
        })
    }
}

fn source_error(stage: &'static str, error: impl core::fmt::Display) -> Error {
    Error::Source(failure(stage, error).into())
}
fn sum(values: impl IntoIterator<Item = f64>) -> f64 {
    let (mut total, mut correction) = (0.0, 0.0);
    for value in values {
        let (next, error) = two_sum(total, value);
        total = next;
        correction += error;
    }
    total + correction
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), Error> {
    if cancelled() {
        Err(Error::Cancelled)
    } else {
        Ok(())
    }
}
