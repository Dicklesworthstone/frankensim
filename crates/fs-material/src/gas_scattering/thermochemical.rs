//! Molecular optics attached to an already accepted NASA-9 gas evaluation.
//!
//! Composition, molar mass, density, heat capacities and transport remain the
//! original thermochemical owner's outputs. Only molecular cross sections and
//! Lorentz–Lorenz polarizabilities are combined here. A gas-field sample never
//! passes through a calorically-perfect `GasSpec` or another gas law.

use super::{GasScatteringError, MolecularOpticsSum, MolecularRayleigh, RayleighGasState};
use crate::fluid_transport::FluidTransport;
use crate::gas_mixture::MAX_GAS_MIXTURE_COMPONENTS;
use crate::thermochemical_gas::ResolvedThermochemicalGasState;
use fs_blake3::DomainHasher;
use fs_thermochem::SpeciesId;

/// Explicit optical reference assigned to an existing canonical species.
/// The assignment is caller authority; no species name supplies optical data.
#[derive(Clone, Debug, PartialEq)]
pub struct SpeciesMolecularRayleigh {
    /// Identity from the original thermochemistry species owner.
    pub species: SpeciesId,
    /// Original reference refractivity, depolarization and source identity.
    pub optical: MolecularRayleigh,
}

/// Original resolved thermochemical state and its scattering/refractivity view.
/// Retains every source/phase admission result without reconstructing the gas.
#[derive(Clone, Debug, PartialEq)]
pub struct ThermochemicalRayleighState {
    resolved: ResolvedThermochemicalGasState,
    optical_sources: Vec<SpeciesMolecularRayleigh>,
    rayleigh: RayleighGasState,
}

impl ThermochemicalRayleighState {
    /// Attach one explicit optical row per active accepted species. Input rows
    /// may be unordered; missing, duplicate or foreign species refuse atomically.
    /// Work/storage are bounded by the shared 128-species material limit.
    ///
    /// Uses the accepted canonical **mole** fractions and evaluated molar mass,
    /// not a new basis conversion or a molar mass recovered from another gas
    /// constant. `AxialGasField::state_at` results can be passed directly here.
    /// The original thermochemical source domain and phase admission survive.
    ///
    /// The model identity binds species, component masses, composition, optical
    /// sources, original transport cards and source packs independently of T/p.
    /// Selected NASA-9 regions and evaluated transport receipts belong to the
    /// full state identity: crossing a source temperature-region boundary does
    /// not change the declared molecular optics. A fixed-model render grid
    /// still refuses differing composition or optical references between cells.
    ///
    /// Molecular polarizability remains the explicitly constant visible-band
    /// model; temperature-dependent heat capacities do not establish optical
    /// dispersion, Raman scattering, absorption or thermal emission.
    pub fn try_new(
        resolved: ResolvedThermochemicalGasState,
        optical_sources: &[SpeciesMolecularRayleigh],
    ) -> Result<Self, GasScatteringError> {
        resolved
            .state()
            .validate_coherence()
            .map_err(GasScatteringError::Gas)?;
        let receipt = resolved.thermodynamics().receipt();
        let count = receipt.species().len();
        if count == 0
            || count > MAX_GAS_MIXTURE_COMPONENTS
            || optical_sources.len() != count
            || receipt.mole_fraction_bits().len() != count
            || receipt.component_receipts().len() != count
            || resolved.components().len() != count
        {
            return Err(GasScatteringError::Components);
        }
        let mut optical_sources = optical_sources.to_vec();
        optical_sources.sort_by(|a, b| a.species.cmp(&b.species));
        if optical_sources
            .iter()
            .zip(receipt.species())
            .any(|(row, species)| &row.species != species)
        {
            return Err(GasScatteringError::Components);
        }
        let mut optics = MolecularOpticsSum::default();
        let mut hash = DomainHasher::new("org.frankensim.material.thermochemical-rayleigh.v1");
        words(
            &mut hash,
            &[
                count as u64,
                receipt.molar_mass_bits(),
                resolved.state().water_mole_fraction.to_bits(),
            ],
        );
        for (index, row) in optical_sources.iter().enumerate() {
            let fraction_bits = receipt.mole_fraction_bits()[index];
            optics.add(f64::from_bits(fraction_bits), row.optical)?;
            text(&mut hash, row.species.as_str());
            words(
                &mut hash,
                &[
                    fraction_bits,
                    receipt.component_receipts()[index].molar_mass_bits(),
                ],
            );
            row.optical.hash_into(&mut hash);
            hash.update(resolved.components()[index].card_identity().as_bytes());
        }
        words(&mut hash, &[resolved.source_packs().len() as u64]);
        for pack in resolved.source_packs() {
            text(&mut hash, pack.species.as_str());
            hash.update(pack.species_pack.as_bytes());
            hash.update(pack.nasa9_pack.as_bytes());
            hash.update(pack.transport_pack.as_bytes());
        }
        let model_identity = hash.finalize();
        hash_evaluation(&resolved, &mut hash)?;
        let rayleigh = optics.finish(
            *resolved.state(),
            resolved.thermodynamics().properties().molar_mass().value(),
            model_identity,
            hash,
        )?;
        Ok(Self {
            resolved,
            optical_sources,
            rayleigh,
        })
    }

    /// Original resolved gas, including all thermochemical and material receipts.
    #[must_use]
    pub const fn resolved(&self) -> &ResolvedThermochemicalGasState {
        &self.resolved
    }

    /// Original optical assignments, paired with canonical active species.
    #[must_use]
    pub fn optical_sources(&self) -> &[SpeciesMolecularRayleigh] {
        &self.optical_sources
    }

    /// Same accepted numerical gas with exact molecular spectrum/phase. This
    /// view is compatible with the native borrowed gas-cell rendering adapter.
    #[must_use]
    pub const fn rayleigh(&self) -> &RayleighGasState {
        &self.rayleigh
    }
}

// Bind the original accepted evaluation, including caloric-reference changes
// that need not alter rho/cp or any optical coefficient. Existing material
// bundle identities already cover their complete queries and source receipts.
fn hash_evaluation(
    state: &ResolvedThermochemicalGasState,
    hash: &mut DomainHasher,
) -> Result<(), GasScatteringError> {
    let mixture = state.thermodynamics().receipt();
    words(
        hash,
        &[
            mixture.evaluator_version() as u64,
            mixture.gas_constant_bits(),
            mixture.temperature_bits(),
            mixture.pressure_bits(),
            mixture.reference_pressure_bits(),
            mixture.phase() as u64,
            mixture.eos() as u64,
            mixture.declared_basis() as u64,
            mixture.declared_fraction_sum_bits(),
            mixture.mole_fraction_sum_bits(),
            mixture.x_log_x_bits(),
        ],
    );
    text(hash, mixture.fs_math_version());
    text(hash, mixture.fs_qty_version());
    text(hash, mixture.elemental_reference().as_str());
    words(hash, mixture.declared_fraction_bits());
    for component in mixture.component_receipts() {
        words(
            hash,
            &[
                component.evaluator_version() as u64,
                component.gas_constant_bits(),
                component.temperature_bits(),
                component.reference_pressure_bits(),
                component.phase() as u64,
                component.eos() as u64,
                component.region_index() as u64,
                component.region_min_bits(),
                component.region_max_bits(),
            ],
        );
        text(hash, component.fs_math_version());
        text(hash, component.elemental_reference().as_str());
        words(hash, component.coefficient_bits());
        hash.update(component.card_identity());
        hash.update(&[u8::from(component.caloric_reference().is_some())]);
        if let Some(reference) = component.caloric_reference() {
            words(
                hash,
                &[
                    reference.evaluator_version() as u64,
                    reference.anchor_temperature().value().to_bits(),
                    reference.anchor_region() as u64,
                    reference.maximum_enthalpy_correction().value().to_bits(),
                    reference.maximum_entropy_correction().value().to_bits(),
                    reference.regions().len() as u64,
                ],
            );
            for region in reference.regions() {
                hash.update(region.card_identity());
                words(
                    hash,
                    &[
                        region.enthalpy_offset().value().to_bits(),
                        region.entropy_offset().value().to_bits(),
                    ],
                );
            }
        }
    }
    for material in state.components() {
        hash.update(material.identity().as_bytes());
    }
    let eos = state.equation_of_state().receipt();
    text(hash, eos.fs_qty_version());
    words(
        hash,
        &[
            eos.evaluator_version() as u64,
            eos.eos() as u64,
            eos.gas_constant_bits(),
            eos.molar_mass_bits(),
            eos.temperature_bits(),
            eos.pressure_bits(),
            eos.mass_density_bits(),
            eos.molar_volume_bits(),
            eos.specific_gas_constant_bits(),
            eos.compressibility_factor_bits(),
        ],
    );
    let viscosity = state.viscosity().receipt();
    text(hash, viscosity.fs_math_version());
    text(hash, viscosity.fs_qty_version());
    words(
        hash,
        &[
            viscosity.evaluator_version() as u64,
            viscosity.rule() as u64,
            viscosity.temperature_bits(),
            viscosity.declared_basis() as u64,
            viscosity.declared_fraction_sum_bits(),
            viscosity.mole_fraction_sum_bits(),
            viscosity.mixture_viscosity_bits(),
        ],
    );
    for species in viscosity.species() {
        text(hash, species.as_str());
    }
    for values in [
        viscosity.molar_mass_bits(),
        viscosity.pure_viscosity_bits(),
        viscosity.declared_fraction_bits(),
        viscosity.mole_fraction_bits(),
        viscosity.denominator_bits(),
    ] {
        words(hash, values);
    }
    let properties = state.thermodynamics().properties();
    let transport = state.transport_state().map_err(GasScatteringError::Gas)?;
    for value in [
        properties.cp().value(),
        properties.cv().value(),
        properties.h().value(),
        properties.u().value(),
        properties.s().value(),
        properties.g().value(),
        properties.mass_specific_cp().value(),
        properties.mass_specific_cv().value(),
        properties.mass_specific_gas_constant().value(),
        state.specific_internal_energy_j_kg(),
        state.specific_enthalpy_j_kg(),
        transport.density_kg_m3(),
        transport.specific_heat_j_kg_k(),
        transport.thermal_conductivity_w_m_k(),
        transport.dynamic_viscosity_pa_s(),
        transport.kinematic_viscosity_m2_s(),
        transport.thermal_diffusivity_m2_s(),
        transport.prandtl(),
    ] {
        hash.update(&value.to_bits().to_le_bytes());
    }
    Ok(())
}

fn words(hash: &mut DomainHasher, values: &[u64]) {
    hash.update(&(values.len() as u64).to_le_bytes());
    for value in values {
        hash.update(&value.to_le_bytes());
    }
}
fn text(hash: &mut DomainHasher, value: &str) {
    hash.update(&(value.len() as u64).to_le_bytes());
    hash.update(value.as_bytes());
}

#[cfg(test)]
mod tests;
