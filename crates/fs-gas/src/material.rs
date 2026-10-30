//! Finite-rate sealed reactors using the common source-bearing gas parcel.
//!
//! The L3 kinetics owner chooses the composition after an actual elapsed time.
//! The existing material owner then reconstructs one state at the parcel's
//! original mass, volume and total internal energy. Transport and acoustics use
//! those same products and sources. No additional reaction heat is applied.

/// Shared material-card viscosity/conductivity in conservative gas transport.
#[cfg(feature = "experimental-euler")]
pub mod transport;

use crate::chemistry::{
    IsochoricKineticsError, IsochoricKineticsResult, IsochoricKineticsSolve,
    advance_adiabatic_energy,
};
use fs_material::MaterialError;
use fs_material::thermochemical_gas::ThermochemicalGasModel;
use fs_material::thermochemical_gas::energy::parcel::ClosedGasParcel;
use fs_material::thermochemical_gas::energy::{GasEnergyError, GasEnergySolve};
use fs_qty::{Amount, Energy, Temperature, Time, Volume};
use fs_thermochem::SpeciesId;
use fs_thermochem::kinetics::{IdealGasKineticsV1, KineticsErrorV1};

/// Explicit time-integration, material-inversion and cross-owner tolerances.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MaterialKineticsSolve {
    /// Adaptive finite-rate integration with one conserved total internal energy.
    pub chemistry: IsochoricKineticsSolve,
    /// Original material energy inversion; its tolerance is in J/kg.
    pub material_energy: GasEnergySolve,
    /// Maximum feed/product/density mass discrepancy [kg].
    pub mass_tolerance_kg: f64,
    /// Maximum relative temperature or pressure change while resolving the
    /// material endpoint. This bounds source/roundoff reconciliation, not a
    /// permitted physical heating or expansion process.
    pub relative_state_tolerance: f64,
}

/// Actual balances after material resolution, in addition to reactor balances.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MaterialKineticsBalance {
    /// Original feed reconstructed from parcel mass and molar composition [kg].
    pub feed_mass_kg: f64,
    /// Feed mass minus exactly preserved parcel mass [kg].
    pub feed_mass_residual_kg: f64,
    /// Integrated species mass minus exactly preserved parcel mass [kg].
    pub chemical_mass_residual_kg: f64,
    /// Final resolved density times volume minus preserved parcel mass [kg].
    pub density_mass_residual_kg: f64,
    /// Actual final material internal energy [J], including formation energy.
    pub stored_internal_energy_j: f64,
    /// Actual stored internal energy minus the original parcel energy [J].
    pub energy_residual_j: f64,
    /// Material energy minus the finite-rate chemistry endpoint energy [J].
    pub freeze_energy_residual_j: f64,
    /// Relative temperature change at material resolution.
    pub relative_temperature_correction: f64,
    /// Relative pressure change at material resolution.
    pub relative_pressure_correction: f64,
    /// Largest absolute mole-fraction correction at the exact-simplex boundary.
    pub maximum_mole_fraction_correction: f64,
}

/// An accepted finite-time reaction update and its actual material balances.
#[derive(Clone, Debug, PartialEq)]
pub struct MaterialKineticsStep {
    /// Original mechanism, source receipts, integrated species and work counts.
    pub chemistry: IsochoricKineticsResult,
    /// Cross-owner checks on the state actually published to the parcel.
    pub balance: MaterialKineticsBalance,
}

/// Atomic refusal: the original parcel remains unchanged.
#[derive(Clone, Debug, PartialEq)]
pub enum MaterialKineticsError {
    /// Invalid controls or unrepresentable extensive state.
    Invalid(&'static str),
    /// The complete bank or active feed does not retain the original species law.
    SourceMismatch(SpeciesId),
    /// Finite-rate chemistry refused the state or exhausted its work budget.
    Chemistry(IsochoricKineticsError),
    /// Material source selection or composition admission failed.
    Material(MaterialError),
    /// The original material energy inversion refused the endpoint.
    Energy(GasEnergyError),
    /// Actual material balances exceeded the requested tolerances.
    Balance(MaterialKineticsBalance),
    /// Cancellation before the complete update was published.
    Cancelled,
}

impl core::fmt::Display for MaterialKineticsError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "source-bearing gas reactor: {self:?}")
    }
}
impl std::error::Error for MaterialKineticsError {}
impl From<MaterialError> for MaterialKineticsError {
    fn from(error: MaterialError) -> Self {
        Self::Material(error)
    }
}
impl From<GasEnergyError> for MaterialKineticsError {
    fn from(error: GasEnergyError) -> Self {
        match error {
            GasEnergyError::Cancelled => Self::Cancelled,
            other => Self::Energy(other),
        }
    }
}
impl From<IsochoricKineticsError> for MaterialKineticsError {
    fn from(error: IsochoricKineticsError) -> Self {
        match error {
            IsochoricKineticsError::Cancelled
            | IsochoricKineticsError::Chemistry(KineticsErrorV1::Cancelled) => Self::Cancelled,
            other => Self::Chemistry(other),
        }
    }
}

/// Advance an existing rigid, adiabatic gas parcel through finite-rate chemistry.
///
/// The mechanism and complete product bank must share the exact original NASA-9
/// models. Active feed transport cards, water identity and source packs must
/// also match. Initially absent products may appear through the admitted
/// reactions. Successful publication preserves mass and volume bit-for-bit;
/// failed kinetics, source resolution, balances or cancellation preserve every
/// field of the original parcel. Zero duration is an exact parcel no-op after
/// admission. Material acoustics are frozen-composition at the new state.
#[allow(clippy::too_many_arguments)]
pub fn advance_closed_parcel_adiabatic(
    parcel: &mut ClosedGasParcel,
    bank: &ThermochemicalGasModel,
    mechanism: &IdealGasKineticsV1,
    duration: Time,
    controls: MaterialKineticsSolve,
    mut cancelled: impl FnMut() -> bool,
) -> Result<MaterialKineticsStep, MaterialKineticsError> {
    poll(&mut cancelled)?;
    controls.material_energy.validate()?;
    if !positive(controls.mass_tolerance_kg)
        || !positive(controls.relative_state_tolerance)
        || controls.relative_state_tolerance >= 1.0
    {
        return Err(MaterialKineticsError::Invalid(
            "positive finite mass and relative state tolerances below one required",
        ));
    }
    admit_sources(parcel, bank, mechanism)?;
    let mass = parcel.mass_kg();
    let volume = parcel.volume_m3();
    let amount = mass
        / parcel
            .state()
            .thermodynamics()
            .properties()
            .molar_mass()
            .value();
    let target = mass * parcel.state().specific_internal_energy_j_kg();
    let material_tolerance = mass * controls.material_energy.energy_tolerance_j_kg;
    if !positive(amount) || !target.is_finite() || !positive(material_tolerance) {
        return Err(MaterialKineticsError::Invalid(
            "unrepresentable parcel amount, energy or tolerance",
        ));
    }
    let feed: Vec<_> = parcel
        .model()
        .thermodynamics()
        .components()
        .iter()
        .zip(
            parcel
                .model()
                .thermodynamics()
                .mole_composition()
                .fractions(),
        )
        .map(|(species, &fraction)| (species.species().clone(), Amount::new(amount * fraction)))
        .collect();
    if feed.iter().any(|(_, amount)| !positive(amount.value())) {
        return Err(MaterialKineticsError::Invalid(
            "positive feed species amount is unrepresentable",
        ));
    }
    let feed_mass = compensated_sum(feed.iter().map(|(id, amount)| {
        let model = mechanism
            .thermodynamics()
            .models()
            .iter()
            .find(|model| model.species() == id)
            .expect("source admission checked active species");
        amount.value() * model.molar_mass().value()
    }));
    if !feed_mass.is_finite() || (feed_mass - mass).abs() > controls.mass_tolerance_kg {
        return Err(MaterialKineticsError::Invalid(
            "reconstructed feed does not preserve parcel mass",
        ));
    }
    let chemistry = advance_adiabatic_energy(
        mechanism,
        &feed,
        Energy::new(target),
        Volume::new(volume),
        Temperature::new(parcel.state().state().temperature),
        duration,
        controls.chemistry,
        &mut cancelled,
    )?;
    poll(&mut cancelled)?;
    let mut candidate = parcel.clone();
    if duration.value() != 0.0 {
        let products = bank.with_mole_amounts(chemistry.amounts())?;
        candidate.replace_model_conserving_mass_energy(
            products,
            controls.material_energy,
            &mut cancelled,
        )?;
    }
    let total = compensated_sum(chemistry.amounts().iter().map(|n| n.value()));
    if !positive(total) {
        return Err(MaterialKineticsError::Invalid(
            "unrepresentable final amount",
        ));
    }
    let product_models = candidate.model().thermodynamics().components();
    let product_fractions = candidate
        .model()
        .thermodynamics()
        .mole_composition()
        .fractions();
    let correction = mechanism
        .thermodynamics()
        .models()
        .iter()
        .zip(chemistry.amounts())
        .map(|(model, amount)| {
            let fraction = product_models
                .iter()
                .position(|p| p.species() == model.species())
                .map_or(0.0, |i| product_fractions[i]);
            (fraction - amount.value() / total).abs()
        })
        .fold(0.0_f64, f64::max);
    let state = candidate.state();
    let energy = mass * state.specific_internal_energy_j_kg();
    let balance = MaterialKineticsBalance {
        feed_mass_kg: feed_mass,
        feed_mass_residual_kg: feed_mass - mass,
        chemical_mass_residual_kg: chemistry.mass().value() - mass,
        density_mass_residual_kg: state.state().density * volume - mass,
        stored_internal_energy_j: energy,
        energy_residual_j: energy - target,
        freeze_energy_residual_j: energy - chemistry.internal_energy().value(),
        relative_temperature_correction: state.state().temperature
            / chemistry.temperature().value()
            - 1.0,
        relative_pressure_correction: state.state().pressure / chemistry.pressure().value() - 1.0,
        maximum_mole_fraction_correction: correction,
    };
    let freezing_tolerance = material_tolerance + controls.chemistry.energy_tolerance_j;
    if !positive(freezing_tolerance)
        || ![
            balance.chemical_mass_residual_kg,
            balance.density_mass_residual_kg,
            balance.stored_internal_energy_j,
            balance.energy_residual_j,
            balance.freeze_energy_residual_j,
            balance.relative_temperature_correction,
            balance.relative_pressure_correction,
            balance.maximum_mole_fraction_correction,
        ]
        .iter()
        .all(|value| value.is_finite())
        || balance.chemical_mass_residual_kg.abs() > controls.mass_tolerance_kg
        || balance.density_mass_residual_kg.abs() > controls.mass_tolerance_kg
        || balance.energy_residual_j.abs() > material_tolerance
        || balance.freeze_energy_residual_j.abs() > freezing_tolerance
        || balance.relative_temperature_correction.abs() > controls.relative_state_tolerance
        || balance.relative_pressure_correction.abs() > controls.relative_state_tolerance
        || balance.maximum_mole_fraction_correction
            > 64.0 * f64::EPSILON * chemistry.amounts().len() as f64
    {
        return Err(MaterialKineticsError::Balance(balance));
    }
    poll(&mut cancelled)?;
    *parcel = candidate;
    Ok(MaterialKineticsStep { chemistry, balance })
}

fn admit_sources(
    parcel: &ClosedGasParcel,
    bank: &ThermochemicalGasModel,
    mechanism: &IdealGasKineticsV1,
) -> Result<(), MaterialKineticsError> {
    let models = bank.thermodynamics().components();
    let chemical = mechanism.thermodynamics().models();
    if models.len() != chemical.len() {
        return Err(MaterialKineticsError::Invalid(
            "material bank must cover the complete mechanism",
        ));
    }
    for (material, source) in models.iter().zip(chemical) {
        if material != source {
            return Err(MaterialKineticsError::SourceMismatch(
                source.species().clone(),
            ));
        }
    }
    let active = parcel.model().thermodynamics().components();
    for (i, source) in active.iter().enumerate() {
        let Some(j) = models.iter().position(|m| m.species() == source.species()) else {
            return Err(MaterialKineticsError::SourceMismatch(
                source.species().clone(),
            ));
        };
        if source != &models[j]
            || parcel.model().transport_materials()[i] != bank.transport_materials()[j]
        {
            return Err(MaterialKineticsError::SourceMismatch(
                source.species().clone(),
            ));
        }
    }
    let active_water = bank
        .water_species()
        .filter(|water| active.iter().any(|source| source.species() == *water));
    if active_water != parcel.model().water_species() {
        return Err(MaterialKineticsError::Invalid(
            "active feed water identities differ",
        ));
    }
    let packs: Vec<_> = bank
        .source_packs()
        .iter()
        .filter(|pack| {
            active
                .iter()
                .any(|source| source.species() == &pack.species)
        })
        .collect();
    if packs.len() != parcel.model().source_packs().len()
        || packs
            .iter()
            .zip(parcel.model().source_packs())
            .any(|(a, b)| *a != b)
    {
        return Err(MaterialKineticsError::Invalid(
            "active feed source-pack identities differ",
        ));
    }
    Ok(())
}

fn positive(value: f64) -> bool {
    value.is_finite() && value > 0.0
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), MaterialKineticsError> {
    if cancelled() {
        Err(MaterialKineticsError::Cancelled)
    } else {
        Ok(())
    }
}
fn compensated_sum(values: impl Iterator<Item = f64>) -> f64 {
    let (mut sum, mut error) = (0.0, 0.0);
    for value in values {
        let adjusted = value - error;
        let next = sum + adjusted;
        error = (next - sum) - adjusted;
        sum = next;
    }
    sum
}

#[cfg(test)]
mod tests;
