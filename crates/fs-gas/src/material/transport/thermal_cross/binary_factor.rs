//! Binary thermal-diffusion factors feeding the original reciprocal face law.
//!
//! The source supplies the local dimensionless factor alpha_12 defined by
//! grad ln(x_2/x_1) = alpha_12 grad ln(T) at zero species flux and constant p.
//! Species order is explicit; reversing it reverses alpha's sign. An end-to-end
//! experimental ratio is usable only if its source supports this local law.
//!
//! In the existing mass frame this gives
//! H_1-H_2 = alpha_12 R T (x_1/M_2 + x_2/M_1).
//! Set H_2=0 as a deterministic kinetic-heat gauge. The existing face owner
//! supplies both Soret force and reciprocal Dufour heat; no new flux kernel,
//! caloric approximation or inferred transport coefficient is introduced.
//! Van-Brunt, Farrell and Monroe, arXiv:2109.05082, section 6.3.

use super::{ThermalConductivityConvention, face_query, reference_state};
use super::super::{material_input, poll};
use crate::multicomponent::dissipation::{
    DissipationError, TransportCoefficients, TransportLaw, validate_binary,
};
use fs_matdb::{MaterialCard, PropertyKey, QueryPoint};
use fs_material::state_point::{
    MaterialPropertySelection, ResolvedMaterialStatePoint, ScalarAdmissibility,
    ScalarPropertyRequirement, resolve_material_state_point,
};
use fs_qty::{Dims, Pressure, QuantitySpec, Temperature};
use fs_thermochem::equilibrium::composition::IdealGasEquilibriumV1;
use fs_thermochem::{
    MassSpecificEnergyQuantityV1, SpeciesId, UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R,
};

/// One oriented pair's local, dimensionless thermal-diffusion factor.
#[derive(Clone, Debug, PartialEq)]
pub struct BinaryThermalDiffusionFactorMaterial {
    /// [species 1, species 2] in grad ln(x_2/x_1) = alpha_12 grad ln(T).
    /// This order is physical and is never inferred from names or molar masses.
    pub species: [SpeciesId; 2],
    pub material: MaterialCard,
    /// Explicit finite signed dimensionless factor, not a Soret coefficient [1/K].
    pub property: PropertyKey,
    /// Select once at the reference query, then retain that exact claim.
    pub selection: MaterialPropertySelection,
}

/// One complete source query and its derived canonical kinetic heats.
#[derive(Clone, Debug, PartialEq)]
pub struct ResolvedBinaryThermalDiffusionFactor {
    temperature: Temperature,
    pressure: Pressure,
    requested_fractions: [f64; 2],
    mole_fractions: [f64; 2],
    factor: f64,
    heats: [MassSpecificEnergyQuantityV1; 2],
    material: ResolvedMaterialStatePoint,
}
impl ResolvedBinaryThermalDiffusionFactor {
    pub const fn temperature(&self) -> Temperature {
        self.temperature
    }
    pub const fn pressure(&self) -> Pressure {
        self.pressure
    }
    /// Mole fractions in the underlying canonical caloric-species order.
    pub const fn mole_fractions(&self) -> &[f64; 2] {
        &self.mole_fractions
    }
    /// Original signed alpha_12, using the source assignment's pair orientation.
    pub const fn factor(&self) -> f64 {
        self.factor
    }
    /// Derived H_i [J/kg] in canonical order; the source's second species has H=0.
    /// These kinetic coefficients are separate from partial thermodynamic enthalpy.
    pub const fn heats(&self) -> &[MassSpecificEnergyQuantityV1; 2] {
        &self.heats
    }
    /// Original selected factor claim and actual query coordinates.
    pub const fn material(&self) -> &ResolvedMaterialStatePoint {
        &self.material
    }
    fn matches(&self, t: Temperature, p: Pressure, x: &[f64]) -> bool {
        self.temperature == t && self.pressure == p && self.requested_fractions.as_slice() == x
    }
}

/// Source-qualified binary Soret/Dufour transport for all original gas-face owners.
///
/// The underlying law owns diffusivity, viscosity and conductivity. Its
/// conductivity must use the explicitly declared zero-diffusive-species-flux
/// convention. The factor alone does not supply diffusivity or conductivity.
///
/// Exactly two canonical caloric species are required, including either species
/// when its current concentration is zero. No trace is added and no binary data
/// are extrapolated to a third species. All source T,p,composition restrictions
/// are enforced at each query. The existing numerical owner handles reference
/// frames, reciprocal heat, conservation and finite-step entropy admission.
///
/// Successful repeated evaluate/heats queries at one face reuse the same factor
/// resolution. Failure preserves output buffers and the last complete factor
/// receipt. Source resolution allocates; the conversion itself is constant work.
pub struct MaterialBinarySoretDufourTransport<L: TransportLaw> {
    law: L,
    convention: ThermalConductivityConvention,
    source: BinaryThermalDiffusionFactorMaterial,
    composition_axes: Vec<(SpeciesId, String)>,
    reference_query: QueryPoint,
    molar_masses_kg_mol: [f64; 2],
    first_species: usize,
    resolved: ResolvedBinaryThermalDiffusionFactor,
    binary_workspace: [f64; 4],
}
impl<L: TransportLaw> MaterialBinarySoretDufourTransport<L> {
    /// Bind an oriented local factor, exact binary bank and typed mole-fraction axes.
    pub fn new(
        law: L,
        convention: ThermalConductivityConvention,
        source: BinaryThermalDiffusionFactorMaterial,
        mut composition_axes: Vec<(SpeciesId, String)>,
        reference_query: QueryPoint,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, DissipationError> {
        poll(&mut cancelled)?;
        let models = law.source().models();
        if models.len() != 2
            || source.species[0] == source.species[1]
            || source.species.iter().any(|id| !models.iter().any(|m| m.species() == id))
            || source.property.quantity() != QuantitySpec::dimensional(Dims::NONE)
        {
            return Err(DissipationError::Invalid(
                "binary thermal factor requires two original species and a dimensionless property",
            ));
        }
        let first_species = if &source.species[0] == models[0].species() { 0 } else { 1 };
        let molar_masses_kg_mol = [
            models[0].molar_mass().value(), models[1].molar_mass().value(),
        ];
        let (temperature, pressure, fractions) =
            reference_state(law.source(), &mut composition_axes, &reference_query)?;
        let (query, canonical) = face_query(
            &reference_query, &composition_axes, temperature, pressure, &fractions,
        )?;
        let requirement = ScalarPropertyRequirement::try_with_key(
            &source.property, ScalarAdmissibility::Finite,
        ).map_err(|e| material_input("thermal diffusion factor requirement", e))?;
        let material = resolve_material_state_point(
            &source.material, &query, &[requirement], source.selection.clone(),
        ).map_err(|e| material_input("thermal diffusion factor source", e))?;
        let factor = factor(&source, &material)?;
        let heats = kinetic_heats(factor, temperature, &canonical,
            molar_masses_kg_mol, first_species)?;
        poll(&mut cancelled)?;
        Ok(Self {
            law, convention, source, composition_axes, reference_query,
            molar_masses_kg_mol, first_species,
            resolved: ResolvedBinaryThermalDiffusionFactor {
                temperature, pressure,
                requested_fractions: [fractions[0], fractions[1]],
                mole_fractions: [canonical[0], canonical[1]],
                factor, heats, material,
            },
            binary_workspace: [0.0; 4],
        })
    }

    pub const fn law(&self) -> &L {
        &self.law
    }
    pub const fn conductivity_convention(&self) -> ThermalConductivityConvention {
        self.convention
    }
    /// Original assignment, including the physical sign/orientation convention.
    pub const fn factor_source(&self) -> &BinaryThermalDiffusionFactorMaterial {
        &self.source
    }
    pub fn composition_axes(&self) -> &[(SpeciesId, String)] {
        &self.composition_axes
    }
    pub const fn reference_query(&self) -> &QueryPoint {
        &self.reference_query
    }
    pub const fn last_resolved(&self) -> &ResolvedBinaryThermalDiffusionFactor {
        &self.resolved
    }

    fn resolve(
        &self,
        temperature: Temperature,
        pressure: Pressure,
        fractions: &[f64],
        cancelled: &mut dyn FnMut() -> bool,
    ) -> Result<ResolvedBinaryThermalDiffusionFactor, DissipationError> {
        poll(cancelled)?;
        let (query, canonical) = face_query(
            &self.reference_query, &self.composition_axes, temperature, pressure, fractions,
        )?;
        let material = self.resolved.material.at_query_point(&self.source.material, &query)
            .map_err(|e| material_input("thermal diffusion factor refresh", e))?;
        let factor = factor(&self.source, &material)?;
        let heats = kinetic_heats(factor, temperature, &canonical,
            self.molar_masses_kg_mol, self.first_species)?;
        poll(cancelled)?;
        Ok(ResolvedBinaryThermalDiffusionFactor {
            temperature, pressure,
            requested_fractions: [fractions[0], fractions[1]],
            mole_fractions: [canonical[0], canonical[1]],
            factor, heats, material,
        })
    }
}
impl<L: TransportLaw> TransportLaw for MaterialBinarySoretDufourTransport<L> {
    fn source(&self) -> &IdealGasEquilibriumV1 {
        self.law.source()
    }
    fn evaluate(
        &mut self,
        temperature: Temperature,
        pressure: Pressure,
        fractions: &[f64],
        binary_diffusivities: &mut [f64],
        cancelled: &mut dyn FnMut() -> bool,
    ) -> Result<TransportCoefficients, DissipationError> {
        poll(cancelled)?;
        if binary_diffusivities.len() != 4 {
            return Err(DissipationError::Invalid("binary thermal factor output axis"));
        }
        let candidate = if self.resolved.matches(temperature, pressure, fractions) {
            None
        } else {
            Some(self.resolve(temperature, pressure, fractions, cancelled)?)
        };
        let resolved = candidate.as_ref().unwrap_or(&self.resolved);
        let coefficients = self.law.evaluate(
            temperature, pressure, &resolved.mole_fractions,
            &mut self.binary_workspace, cancelled,
        )?;
        validate_binary(&self.binary_workspace, 2)?;
        poll(cancelled)?;
        binary_diffusivities.copy_from_slice(&self.binary_workspace);
        if let Some(candidate) = candidate {
            self.resolved = candidate;
        }
        Ok(coefficients)
    }
    fn heats_of_transport(
        &mut self,
        temperature: Temperature,
        pressure: Pressure,
        fractions: &[f64],
        output_j_kg: &mut [f64],
        cancelled: &mut dyn FnMut() -> bool,
    ) -> Result<bool, DissipationError> {
        poll(cancelled)?;
        if output_j_kg.len() != 2 {
            return Err(DissipationError::Invalid("binary thermal factor heat output axis"));
        }
        let candidate = if self.resolved.matches(temperature, pressure, fractions) {
            None
        } else {
            Some(self.resolve(temperature, pressure, fractions, cancelled)?)
        };
        let resolved = candidate.as_ref().unwrap_or(&self.resolved);
        poll(cancelled)?;
        for (output, heat) in output_j_kg.iter_mut().zip(&resolved.heats) {
            *output = heat.value();
        }
        if let Some(candidate) = candidate {
            self.resolved = candidate;
        }
        Ok(true)
    }
}

fn factor(
    source: &BinaryThermalDiffusionFactorMaterial,
    resolved: &ResolvedMaterialStatePoint,
) -> Result<f64, DissipationError> {
    let property = resolved.property_by_key(&source.property)
        .ok_or(DissipationError::Invalid("missing binary thermal diffusion factor"))?;
    if !property.has_only_fixed_external_conditions(Some("temperature"), false) {
        return Err(DissipationError::Invalid(
            "thermal diffusion factor requires static local T,p,composition and external-condition support",
        ));
    }
    Ok(property.value_si())
}

fn kinetic_heats(
    factor: f64,
    temperature: Temperature,
    x: &[f64],
    masses: [f64; 2],
    first: usize,
) -> Result<[MassSpecificEnergyQuantityV1; 2], DissipationError> {
    let second = 1 - first;
    // mean(M)/(M_1 M_2) = x_1/M_2 + x_2/M_1 avoids a molar-mass product.
    let difference = (factor * R) * temperature.value()
        * (x[first] / masses[second] + x[second] / masses[first]);
    if !difference.is_finite() || (factor != 0.0 && difference == 0.0) {
        return Err(DissipationError::Invalid(
            "thermal diffusion factor heat conversion is unrepresentable",
        ));
    }
    let mut heats = [
        MassSpecificEnergyQuantityV1::new(0.0), MassSpecificEnergyQuantityV1::new(0.0),
    ];
    heats[first] = MassSpecificEnergyQuantityV1::new(difference);
    Ok(heats)
}
