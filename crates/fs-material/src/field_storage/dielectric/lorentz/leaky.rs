//! Original Lorentz storage with independently sourced DC Ohmic conduction.
//!
//! The unchanged inertial polarization owner supplies D, P, polarization
//! current, restoring/kinetic storage and pole friction. The independent
//! conductor supplies J = sigma E and sigma |E|^2 heat, with no added storage.
//! DC conductivity must exclude loss already represented by Lorentz damping.
//! This binding retains one complete source bundle and makes no second source
//! query or coefficient fit when an existing dielectric is adopted.
//!
//! The field owner enforces total-current continuity, including both Delta D
//! and integrated conduction. It must not conserve dielectric free charge in
//! a leaking cell. Coefficients and geometry remain fixed; this binding adds
//! no carrier inertia, transport equation, temperature refresh or induction.

use crate::conductor::ELECTRICAL_CONDUCTIVITY_DIMS;

use super::super::leaky::dc_conductivity;
use super::{
    FieldStorageError, LorentzPoleProperties, MaterialCard, MaterialPropertySelection, QueryPoint,
    ResolvedLorentzMaterial, ResolvedMaterialStatePoint, ScalarAdmissibility,
    ScalarPropertyRequirement, requirements, resolve_material_state_point,
};

/// One original resonant dielectric and its independent static DC conductor.
/// Zero conductivity is an explicitly sourced perfect-insulator limit.
#[derive(Clone, Debug, PartialEq)]
pub struct ResolvedLeakyLorentzMaterial {
    dielectric: ResolvedLorentzMaterial,
    conductivity_property: String,
    conductivity_s_m: f64,
}

impl ResolvedLeakyLorentzMaterial {
    /// Adopt a distinct DC conductivity already retained in the dielectric's
    /// complete material snapshot. All pole states and source identities stay
    /// with the original Lorentz owner; no bundle is cloned or re-resolved.
    ///
    /// Missing/aliased roles, incorrect SI or semantic quantities, negative
    /// conductivity, and frequency/field/strain/stress dependence refuse.
    pub fn try_from_dielectric(
        dielectric: ResolvedLorentzMaterial,
        conductivity_property: &str,
    ) -> Result<Self, FieldStorageError> {
        let conductivity_s_m = dc_conductivity(
            dielectric.material(),
            conductivity_property,
            core::iter::once(dielectric.high_frequency_property()).chain(
                dielectric.pole_properties().iter().flat_map(|pole| {
                    [
                        pole.permittivity_increment.as_str(),
                        pole.angular_frequency.as_str(),
                        pole.damping_rate.as_str(),
                    ]
                }),
            ),
        )?;
        Ok(Self {
            dielectric,
            conductivity_property: conductivity_property.to_owned(),
            conductivity_s_m,
        })
    }

    /// Bind inertial polarization and conduction in one complete original
    /// multiphysics snapshot, preserving unrelated properties and receipts.
    pub fn try_from_material(
        material: ResolvedMaterialStatePoint,
        high_frequency_property: &str,
        pole_properties: &[LorentzPoleProperties],
        conductivity_property: &str,
    ) -> Result<Self, FieldStorageError> {
        Self::try_from_dielectric(
            ResolvedLorentzMaterial::try_from_material(
                material,
                high_frequency_property,
                pole_properties,
            )?,
            conductivity_property,
        )
    }

    /// Resolve the original storage, inertia, damping and DC conductor roles
    /// together. Angular natural frequencies retain their semantic quantity;
    /// damping and conductivity independently permit exact zero.
    pub fn resolve(
        card: &MaterialCard,
        point: &QueryPoint,
        selection: MaterialPropertySelection,
        high_frequency_property: &str,
        pole_properties: &[LorentzPoleProperties],
        conductivity_property: &str,
    ) -> Result<Self, FieldStorageError> {
        let mut requirements = requirements(high_frequency_property, pole_properties)?;
        requirements.push(ScalarPropertyRequirement::try_new(
            conductivity_property,
            ELECTRICAL_CONDUCTIVITY_DIMS,
            ScalarAdmissibility::NonNegative,
        )?);
        Self::try_from_material(
            resolve_material_state_point(card, point, &requirements, selection)?,
            high_frequency_property,
            pole_properties,
            conductivity_property,
        )
    }

    /// Complete original snapshot, including all retained multiphysics roles.
    pub const fn material(&self) -> &ResolvedMaterialStatePoint {
        self.dielectric.material()
    }

    /// Unchanged resonant storage and physical pole-friction law.
    pub const fn dielectric(&self) -> &ResolvedLorentzMaterial {
        &self.dielectric
    }

    /// Exact original source name for DC Ohmic leakage.
    pub fn conductivity_property(&self) -> &str {
        &self.conductivity_property
    }

    /// Nonnegative absolute DC conductivity [S/m], excluding pole friction.
    pub const fn conductivity_s_m(&self) -> f64 {
        self.conductivity_s_m
    }

    /// Move the original source owner and admitted conductor into the field
    /// solver. The field supplies conduction current and heat independently
    /// of the Lorentz midpoint recovery; conductivity is never permittivity.
    pub fn into_parts(self) -> (ResolvedLorentzMaterial, String, f64) {
        (
            self.dielectric,
            self.conductivity_property,
            self.conductivity_s_m,
        )
    }
}
