//! A closed reacting parcel driven by accounted heat, work and volume.
//!
//! The caller supplies Q into the gas, W by the gas and the final volume.
//! The prepared native closure solves Delta U = Q-W at fixed original elements;
//! chemistry, temperature and pressure all respond. It does not invent a piston
//! pressure path, reaction rate or heat-transfer coefficient. An optional single
//! reservoir audit checks the actual system-plus-reservoir entropy increment.

use super::{
    FreezeVolumeControls, GasIsochoricEquilibriumChamber, GasIsochoricEquilibriumError,
    ThermochemicalGasModel,
};
use fs_matdb::QueryPoint;
use fs_math::eft::two_sum;
use fs_qty::{Amount, Energy, Temperature, Volume};
use fs_thermochem::SpeciesId;
use fs_thermochem::equilibrium::composition::{
    CompositionEquilibriumErrorV1, EquilibriumExpansionErrorV1, IdealGasEquilibriumV1,
    volume::{
        EquilibriumInternalEnergySolveV1, EquilibriumVolumeErrorV1, IdealGasIsochoricEquilibriumV1,
        energy_path::{EquilibriumEnergyVolumePathV1, EquilibriumEnergyVolumeSampleV1},
    },
};

/// Fixed original source path, work bounds and extensive conservation allowances.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DrivenEquilibriumGasSolve {
    /// Native TV/composition/temperature controls and source-relative U tolerance.
    pub equilibrium: EquilibriumInternalEnergySolveV1,
    /// Original feed/product/material mass discrepancy limit [kg].
    pub mass_tolerance_kg: f64,
    /// Material-freezing, per-step and complete-reference energy allowance [J].
    /// This allowance is not multiplied by the number of updates.
    pub energy_tolerance_j: f64,
    /// Absolute material-minus-native entropy freezing allowance [J/K].
    /// Also bounds a numerical entropy decrease for an adiabatic (Q=0) update.
    pub entropy_tolerance_j_k: f64,
}

/// Original source balances for one accepted heat/work/volume change.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DrivenEquilibriumGasStep {
    /// Exact accepted volume minus previous volume [m3].
    pub volume_change_m3: f64,
    /// Heat supplied INTO the gas [J].
    pub supplied_heat_j: f64,
    /// Work done BY the gas, removing energy [J].
    pub boundary_work_by_gas_j: f64,
    /// Nearest representable Q-W for this step [J].
    pub requested_internal_energy_change_j: f64,
    /// Source-integrated material energy change, including freeze defects [J].
    pub stored_energy_change_j: f64,
    /// Compensated Delta U + W - Q for the actual supplied accounts [J].
    pub first_law_residual_j: f64,
    /// Source-integrated material U-U_birth [J], never a running energy sum.
    pub internal_energy_change_since_reference_j: f64,
    /// Original source U-U_birth minus cumulative supplied Q-W [J].
    pub reference_energy_residual_j: f64,
    /// Difference of the two rounded absolute material storage totals [J].
    /// Informational: at a large formation datum this can lose small heat that
    /// the direct source-integrated energy balance above retains.
    pub rounded_absolute_storage_change_j: f64,
    /// Actual material entropy change [J/K], including mixing and chemistry.
    pub entropy_change_j_k: f64,
    /// Material S-S_birth on the original source reference [J/K].
    pub entropy_change_since_reference_j_k: f64,
    /// All native pressure trials, including bracket and predictor trials.
    pub pressure_iterations: usize,
    /// All native composition iterations, including discarded trials.
    pub composition_iterations: usize,
    /// Native interior temperature iterations.
    pub temperature_iterations: usize,
}

/// An accepted update coupled to one explicit constant-temperature reservoir.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DrivenEquilibriumGasReservoirStep {
    /// Original energy/volume/chemistry balances.
    pub equilibrium: DrivenEquilibriumGasStep,
    /// Entropy lost by the reservoir, -Q/T_reservoir [J/K].
    pub reservoir_entropy_change_j_k: f64,
    /// Gas plus reservoir entropy change [J/K], never clipped to zero.
    pub entropy_production_j_k: f64,
}

/// Failure leaves every accepted physical state and cumulative account unchanged.
#[derive(Clone, Debug, PartialEq)]
pub enum DrivenEquilibriumGasError {
    /// Original native source, chemistry, volume or energy refusal.
    Chemistry(EquilibriumVolumeErrorV1),
    /// Material/source/transport/freezing refusal.
    Material(GasIsochoricEquilibriumError),
    /// Invalid or unrepresentable input or extensive state.
    Invalid(&'static str),
    /// Original step or reference first-law residual exceeded its allowance.
    Balance(DrivenEquilibriumGasStep),
    /// A zero-heat closed update would decrease actual gas entropy beyond its allowance.
    AdiabaticEntropy(DrivenEquilibriumGasStep),
    /// Actual gas-plus-reservoir entropy fell below the declared tolerance.
    Entropy(DrivenEquilibriumGasReservoirStep),
    /// Cancellation before publication.
    Cancelled,
}
impl core::fmt::Display for DrivenEquilibriumGasError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "driven reacting equilibrium parcel: {self:?}")
    }
}
impl std::error::Error for DrivenEquilibriumGasError {}
impl From<EquilibriumVolumeErrorV1> for DrivenEquilibriumGasError {
    fn from(error: EquilibriumVolumeErrorV1) -> Self {
        match error {
            EquilibriumVolumeErrorV1::Equilibrium(CompositionEquilibriumErrorV1::Cancelled)
            | EquilibriumVolumeErrorV1::Response(EquilibriumExpansionErrorV1::Equilibrium(
                CompositionEquilibriumErrorV1::Cancelled,
            )) => Self::Cancelled,
            other => Self::Chemistry(other),
        }
    }
}
impl From<GasIsochoricEquilibriumError> for DrivenEquilibriumGasError {
    fn from(error: GasIsochoricEquilibriumError) -> Self {
        match error {
            GasIsochoricEquilibriumError::Cancelled => Self::Cancelled,
            other => Self::Material(other),
        }
    }
}
type Error = DrivenEquilibriumGasError;
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), Error> {
    if cancelled() {
        Err(Error::Cancelled)
    } else {
        Ok(())
    }
}

// Two-word original input account. This retains small Q/W terms before forming
// the representable native target; physical residuals use both words.
#[derive(Clone, Copy, Default)]
struct Account {
    sum: f64,
    tail: f64,
}
impl Account {
    fn add(&mut self, value: f64) {
        let (s, e) = two_sum(self.sum, value);
        let (sum, tail) = two_sum(s, self.tail + e);
        self.sum = sum;
        self.tail = tail;
    }
    fn value(self) -> f64 {
        self.sum + self.tail
    }
    fn residual(self, stored: f64) -> f64 {
        let (difference, tail) = two_sum(stored, -self.sum);
        difference + (tail - self.tail)
    }
}

/// One accepted material state, immutable original inventory and prepared UV work.
///
/// The numerical workspace may change on a refused query, but the accepted
/// chamber, entropy/energy reference and supplied-input account do not. This
/// represents the infinitely fast ideal-gas chemical-equilibrium limit. Frozen
/// transport products and the separate reacting response are both explicit.
pub struct DrivenEquilibriumGasParcel<'a> {
    bank: &'a ThermochemicalGasModel,
    chemistry: &'a IdealGasEquilibriumV1,
    feed: Vec<(SpeciesId, Amount)>,
    context: QueryPoint,
    path: EquilibriumEnergyVolumePathV1<'a>,
    controls: DrivenEquilibriumGasSolve,
    chamber: GasIsochoricEquilibriumChamber,
    sample: EquilibriumEnergyVolumeSampleV1,
    inputs: Account,
    reference_freeze_energy_j: f64,
    reference_freeze_entropy_j_k: f64,
    stored_since_reference_j: f64,
    entropy_since_reference_j_k: f64,
}
impl core::fmt::Debug for DrivenEquilibriumGasParcel<'_> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        f.debug_struct("DrivenEquilibriumGasParcel")
            .field("chamber", &self.chamber)
            .field("controls", &self.controls)
            .field("energy_input_since_reference_j", &self.inputs.value())
            .finish_non_exhaustive()
    }
}

impl<'a> DrivenEquilibriumGasParcel<'a> {
    /// Establish a reacting equilibrium at original T,V and feed elements.
    /// The complete immutable bank must include all admissible products.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        bank: &'a ThermochemicalGasModel,
        chemistry: &'a IdealGasEquilibriumV1,
        context: &QueryPoint,
        feed: &[(SpeciesId, Amount)],
        temperature: Temperature,
        volume: Volume,
        controls: DrivenEquilibriumGasSolve,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, Error> {
        poll(&mut cancelled)?;
        if !controls.energy_tolerance_j.is_finite()
            || controls.energy_tolerance_j <= 0.0
            || !controls.entropy_tolerance_j_k.is_finite()
            || controls.entropy_tolerance_j_k <= 0.0
        {
            return Err(Error::Invalid(
                "positive finite material energy/entropy allowances required",
            ));
        }
        bank.admit_equilibrium_bank(chemistry, controls.mass_tolerance_kg, &mut cancelled)
            .map_err(GasIsochoricEquilibriumError::from)?;
        let mut path = EquilibriumEnergyVolumePathV1::new(
            chemistry,
            temperature,
            volume,
            feed,
            controls.equilibrium,
            &mut cancelled,
        )?;
        let sample = path.sample(volume, Energy::new(0.0), &mut cancelled)?;
        let native = path.publish_current(&mut cancelled)?;
        let (chamber, reference_freeze_energy_j, freeze_entropy) = freeze(
            bank,
            chemistry,
            context,
            feed,
            native,
            path.feed_mass_kg(),
            controls,
            &mut cancelled,
        )?;
        poll(&mut cancelled)?;
        Ok(Self {
            bank,
            chemistry,
            feed: feed.to_vec(),
            context: context.clone(),
            path,
            controls,
            chamber,
            sample,
            inputs: Account::default(),
            reference_freeze_energy_j,
            reference_freeze_entropy_j_k: freeze_entropy,
            stored_since_reference_j: 0.0,
            entropy_since_reference_j_k: 0.0,
        })
    }
    /// Source-bearing accepted chemical endpoint and common frozen transport view.
    pub const fn chamber(&self) -> &GasIsochoricEquilibriumChamber {
        &self.chamber
    }
    /// Accepted reacting thermodynamics and native relative-energy residual.
    pub const fn sample(&self) -> EquilibriumEnergyVolumeSampleV1 {
        self.sample
    }
    /// Immutable original species/element inventory.
    pub fn feed(&self) -> &[(SpeciesId, Amount)] {
        &self.feed
    }
    /// Exact original feed mass [kg].
    pub fn mass_kg(&self) -> f64 {
        self.chamber.parcel().mass_kg()
    }
    /// Exact accepted physical volume [m3].
    pub fn volume_m3(&self) -> f64 {
        self.chamber.parcel().volume_m3()
    }
    /// Fixed source/conservation/work controls.
    pub const fn controls(&self) -> DrivenEquilibriumGasSolve {
        self.controls
    }
    /// Nearest representable cumulative supplied Q-W since the reference [J].
    pub fn energy_input_since_reference_j(&self) -> f64 {
        self.inputs.value()
    }

    /// Apply supplied heat/work at the prescribed final volume atomically.
    /// No pressure trajectory or work approximation is inferred from the endpoints.
    /// Zero heat requires nondecreasing gas entropy within the fixed allowance.
    /// For nonzero heat, use `advance_with_reservoir` when the reservoir is known;
    /// an unspecified thermal source does not imply a second-law audit.
    pub fn advance(
        &mut self,
        volume: Volume,
        heat: Energy,
        work_by_gas: Energy,
        cancelled: impl FnMut() -> bool,
    ) -> Result<DrivenEquilibriumGasStep, Error> {
        self.advance_inner(volume, heat, work_by_gas, None, cancelled)
    }

    /// Re-equilibrate a rigid chamber after supplied heat.
    pub fn add_heat(
        &mut self,
        heat: Energy,
        cancelled: impl FnMut() -> bool,
    ) -> Result<DrivenEquilibriumGasStep, Error> {
        self.advance(
            Volume::new(self.volume_m3()),
            heat,
            Energy::new(0.0),
            cancelled,
        )
    }

    /// Include one constant-temperature reservoir in the entropy balance.
    ///
    /// Reversible mechanical work carries no entropy. Actual Delta S_gas-Q/T
    /// must be at least minus the explicit numerical allowance, or nothing is
    /// published. This audit does not supply a thermal conductance or time law.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_with_reservoir(
        &mut self,
        volume: Volume,
        heat: Energy,
        work_by_gas: Energy,
        reservoir_temperature: Temperature,
        entropy_tolerance_j_k: f64,
        cancelled: impl FnMut() -> bool,
    ) -> Result<DrivenEquilibriumGasReservoirStep, Error> {
        if !reservoir_temperature.value().is_finite()
            || reservoir_temperature.value() <= 0.0
            || !entropy_tolerance_j_k.is_finite()
            || entropy_tolerance_j_k < 0.0
        {
            return Err(Error::Invalid(
                "positive finite reservoir temperature and nonnegative entropy allowance required",
            ));
        }
        let step = self.advance_inner(
            volume,
            heat,
            work_by_gas,
            Some((reservoir_temperature.value(), entropy_tolerance_j_k)),
            cancelled,
        )?;
        Ok(reservoir_report(step, reservoir_temperature.value()))
    }

    fn advance_inner(
        &mut self,
        volume: Volume,
        heat: Energy,
        work: Energy,
        reservoir: Option<(f64, f64)>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<DrivenEquilibriumGasStep, Error> {
        poll(&mut cancelled)?;
        if !volume.value().is_finite()
            || volume.value() <= 0.0
            || !heat.value().is_finite()
            || !work.value().is_finite()
        {
            return Err(Error::Invalid(
                "positive finite volume and finite heat/work required",
            ));
        }
        let mut inputs = self.inputs;
        inputs.add(heat.value());
        inputs.add(-work.value());
        let target = inputs.value();
        let requested = heat.value() - work.value();
        if !target.is_finite()
            || !requested.is_finite()
            || !inputs.residual(target).is_finite()
            || inputs.residual(target).abs() > self.controls.energy_tolerance_j
        {
            return Err(Error::Invalid(
                "unrepresentable original supplied energy balance",
            ));
        }
        if volume.value() == self.volume_m3() && heat.value() == 0.0 && work.value() == 0.0 {
            return Ok(DrivenEquilibriumGasStep {
                volume_change_m3: 0.0,
                supplied_heat_j: 0.0,
                boundary_work_by_gas_j: 0.0,
                requested_internal_energy_change_j: 0.0,
                stored_energy_change_j: 0.0,
                first_law_residual_j: 0.0,
                internal_energy_change_since_reference_j: self.stored_since_reference_j,
                reference_energy_residual_j: inputs.residual(self.stored_since_reference_j),
                rounded_absolute_storage_change_j: 0.0,
                entropy_change_j_k: 0.0,
                entropy_change_since_reference_j_k: self.entropy_since_reference_j_k,
                pressure_iterations: 0,
                composition_iterations: 0,
                temperature_iterations: 0,
            });
        }
        let sample = self
            .path
            .sample(volume, Energy::new(target), &mut cancelled)?;
        let native = self.path.publish_current(&mut cancelled)?;
        let (candidate, freeze_energy, freeze_entropy) = freeze(
            self.bank,
            self.chemistry,
            &self.context,
            &self.feed,
            native,
            self.mass_kg(),
            self.controls,
            &mut cancelled,
        )?;
        let stored = sample.internal_energy_change_j
            + (freeze_energy - self.reference_freeze_energy_j);
        let entropy =
            sample.entropy_change_j_k + (freeze_entropy - self.reference_freeze_entropy_j_k);
        let change = stored - self.stored_since_reference_j;
        let mut step_inputs = Account::default();
        step_inputs.add(heat.value());
        step_inputs.add(-work.value());
        let report = DrivenEquilibriumGasStep {
            volume_change_m3: volume.value() - self.volume_m3(),
            supplied_heat_j: heat.value(),
            boundary_work_by_gas_j: work.value(),
            requested_internal_energy_change_j: requested,
            stored_energy_change_j: change,
            first_law_residual_j: step_inputs.residual(change),
            internal_energy_change_since_reference_j: stored,
            reference_energy_residual_j: inputs.residual(stored),
            rounded_absolute_storage_change_j: candidate.balance().parcel_internal_energy_j
                - self.chamber.balance().parcel_internal_energy_j,
            entropy_change_j_k: entropy - self.entropy_since_reference_j_k,
            entropy_change_since_reference_j_k: entropy,
            pressure_iterations: sample.pressure_iterations,
            composition_iterations: sample.composition_iterations,
            temperature_iterations: sample.temperature_iterations,
        };
        if [
            change,
            stored,
            entropy,
            report.first_law_residual_j,
            report.reference_energy_residual_j,
            report.entropy_change_j_k,
            report.rounded_absolute_storage_change_j,
        ]
        .iter()
        .any(|x| !x.is_finite())
            || report.first_law_residual_j.abs() > self.controls.energy_tolerance_j
            || report.reference_energy_residual_j.abs() > self.controls.energy_tolerance_j
            || candidate.parcel().mass_kg() != self.mass_kg()
        {
            return Err(Error::Balance(report));
        }
        if heat.value() == 0.0
            && report.entropy_change_j_k < -self.controls.entropy_tolerance_j_k
        {
            return Err(Error::AdiabaticEntropy(report));
        }
        if let Some((temperature, allowance)) = reservoir {
            let audited = reservoir_report(report, temperature);
            if !audited.entropy_production_j_k.is_finite()
                || audited.entropy_production_j_k < -allowance
            {
                return Err(Error::Entropy(audited));
            }
        }
        poll(&mut cancelled)?;
        self.chamber = candidate;
        self.sample = sample;
        self.inputs = inputs;
        self.stored_since_reference_j = stored;
        self.entropy_since_reference_j_k = entropy;
        Ok(report)
    }
}

fn reservoir_report(
    step: DrivenEquilibriumGasStep,
    temperature: f64,
) -> DrivenEquilibriumGasReservoirStep {
    let reservoir = -step.supplied_heat_j / temperature;
    DrivenEquilibriumGasReservoirStep {
        equilibrium: step,
        reservoir_entropy_change_j_k: reservoir,
        entropy_production_j_k: step.entropy_change_j_k + reservoir,
    }
}

#[allow(clippy::too_many_arguments)]
fn freeze(
    bank: &ThermochemicalGasModel,
    chemistry: &IdealGasEquilibriumV1,
    context: &QueryPoint,
    feed: &[(SpeciesId, Amount)],
    state: IdealGasIsochoricEquilibriumV1,
    mass: f64,
    controls: DrivenEquilibriumGasSolve,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(GasIsochoricEquilibriumChamber, f64, f64), Error> {
    let chamber = bank.freeze_volume_equilibrium(
        chemistry,
        context,
        feed,
        state,
        Some(mass),
        FreezeVolumeControls {
            mass_tolerance_kg: controls.mass_tolerance_kg,
            energy_tolerance_j: controls.energy_tolerance_j,
            relative_volume_tolerance: controls.equilibrium.volume.relative_volume_tolerance,
        },
        &mut *cancelled,
    )?;
    let properties = chamber.parcel().state().thermodynamics().properties();
    let entropy = (mass / properties.molar_mass().value()) * properties.s().value();
    let defect = entropy - chamber.equilibrium().thermodynamic().entropy_j_k();
    if !entropy.is_finite() || !defect.is_finite() || defect.abs() > controls.entropy_tolerance_j_k
    {
        return Err(Error::Invalid(
            "material entropy freezing exceeded the declared allowance",
        ));
    }
    // Reconstruct the amount-only freezing energy difference directly, before
    // absolute storage totals round away small heat. Both states use the same
    // species sources and temperature, so only represented amounts can differ.
    let chemical = chamber.equilibrium().equilibrium();
    let mixture = chamber.parcel().model().thermodynamics();
    let total = mass / properties.molar_mass().value();
    let mut energy = Account::default();
    for (model, &fraction) in mixture.components().iter().zip(mixture.mole_composition().fractions()) {
        poll(cancelled)?;
        let product = chemical.species().iter().find(|s| &s.species == model.species())
            .ok_or(Error::Invalid("frozen product lacks its original equilibrium species"))?;
        let source = model.evaluate_properties(chemical.temperature()).map_err(|source| {
            EquilibriumVolumeErrorV1::from(CompositionEquilibriumErrorV1::SpeciesEvaluation {
                species: model.species().clone(), source,
            })
        })?;
        energy.add((total * fraction - product.amount.value()) * source.u().value());
    }
    let energy_defect = energy.value();
    if !energy_defect.is_finite() || energy_defect.abs() > controls.energy_tolerance_j {
        return Err(Error::Invalid("source-relative material freezing exceeded the energy allowance"));
    }
    poll(cancelled)?;
    Ok((chamber, energy_defect, defect))
}

#[cfg(test)]
mod tests;
