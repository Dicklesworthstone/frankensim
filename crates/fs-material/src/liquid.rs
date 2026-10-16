//! Source-resolved single-phase liquid thermal and transport tuple.
//!
//! The caller's material card and complete query point define the liquid-phase
//! conditions. This adapter resolves exactly the four required bulk claims and
//! derives kinematic viscosity, thermal diffusivity, and Prandtl number. It
//! supplies no equation of state, phase boundary, concentration conversion,
//! chemistry-name default, corrosion model, or experimental-accuracy claim.

use fs_matdb::{MaterialCard, QueryPoint};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Density, Dims, QuantitySpec};

use crate::state_point::{
    DENSITY_PROPERTY, MaterialPropertySelection, MaterialStatePointError,
    ResolvedMaterialStatePoint, SPECIFIC_HEAT_CAPACITY_DIMS, SPECIFIC_HEAT_CAPACITY_PROPERTY,
    ScalarAdmissibility, ScalarPropertyRequirement, THERMAL_CONDUCTIVITY_PROPERTY,
    resolve_material_state_point,
};

/// Canonical source property for dynamic viscosity [Pa s].
pub const DYNAMIC_VISCOSITY_PROPERTY: &str = "dynamic_viscosity";
/// Dynamic-viscosity dimensions `[L, M, t, T, I, N]` for Pa s.
pub const DYNAMIC_VISCOSITY_DIMS: Dims = Dims([-1, 1, -1, 0, 0, 0]);
/// Canonical source property for kinematic viscosity [m²/s].
pub const KINEMATIC_VISCOSITY_PROPERTY: &str = "kinematic_viscosity";
/// Kinematic-viscosity dimensions `[L, M, t, T, I, N]` for m²/s.
pub const KINEMATIC_VISCOSITY_DIMS: Dims = Dims([2, 0, -1, 0, 0, 0]);

/// One source-resolved liquid transport state at its complete query point.
///
/// The retained material bundle contains all four selected property-use
/// receipts and the exact query identity. It does not certify that the source
/// conditions remain single phase beyond the claims' own validity envelopes.
#[derive(Clone, Debug, PartialEq)]
pub struct ResolvedLiquidState {
    material: ResolvedMaterialStatePoint,
    density_kg_m3: f64,
    specific_heat_j_kg_k: f64,
    thermal_conductivity_w_m_k: f64,
    dynamic_viscosity_pa_s: f64,
    kinematic_viscosity_m2_s: f64,
    thermal_diffusivity_m2_s: f64,
    prandtl: f64,
}

impl ResolvedLiquidState {
    /// The complete evidence-bearing source bundle and its four receipts.
    #[must_use]
    pub const fn material(&self) -> &ResolvedMaterialStatePoint {
        &self.material
    }

    /// Density [kg/m³].
    #[must_use]
    pub const fn density_kg_m3(&self) -> f64 {
        self.density_kg_m3
    }

    /// Specific heat capacity [J/(kg K)].
    #[must_use]
    pub const fn specific_heat_j_kg_k(&self) -> f64 {
        self.specific_heat_j_kg_k
    }

    /// Thermal conductivity [W/(m K)].
    #[must_use]
    pub const fn thermal_conductivity_w_m_k(&self) -> f64 {
        self.thermal_conductivity_w_m_k
    }

    /// Dynamic viscosity [Pa s].
    #[must_use]
    pub const fn dynamic_viscosity_pa_s(&self) -> f64 {
        self.dynamic_viscosity_pa_s
    }

    /// Kinematic viscosity `ν = μ / ρ` [m²/s].
    #[must_use]
    pub const fn kinematic_viscosity_m2_s(&self) -> f64 {
        self.kinematic_viscosity_m2_s
    }

    /// Thermal diffusivity `α = k / (ρ cp)` [m²/s].
    #[must_use]
    pub const fn thermal_diffusivity_m2_s(&self) -> f64 {
        self.thermal_diffusivity_m2_s
    }

    /// Prandtl number `Pr = μ cp / k`.
    #[must_use]
    pub const fn prandtl(&self) -> f64 {
        self.prandtl
    }
}

/// Resolve one single-phase liquid transport tuple from four sourced claims.
///
/// The card owns phase, composition, pressure, temperature, and any other
/// validity conditions. Every required claim must be positive at the same
/// query point; missing, ambiguous, wrong-unit, pinned, or out-of-domain
/// claims refuse through [`resolve_material_state_point`].
///
/// # Errors
/// Refuses incomplete or inadmissible source bundles and non-finite or
/// non-positive derived transport quantities.
pub fn resolve_liquid_state(
    card: &MaterialCard,
    point: &QueryPoint,
    selection: MaterialPropertySelection,
) -> Result<ResolvedLiquidState, MaterialStatePointError> {
    resolve_liquid_state_with_viscosity(
        card,
        point,
        selection,
        DYNAMIC_VISCOSITY_PROPERTY,
        DYNAMIC_VISCOSITY_DIMS,
        false,
    )
}

/// Resolve one liquid transport tuple when the source reports kinematic
/// viscosity rather than dynamic viscosity.
///
/// This is a separate admission path: it requires the canonical
/// `kinematic_viscosity` claim at the same complete point as density, heat
/// capacity, and thermal conductivity, then derives `mu = rho nu`. It never
/// falls back to, or infers, a dynamic-viscosity claim.
///
/// # Errors
/// Refuses the same source failures as [`resolve_liquid_state`] and an
/// unrepresentable dynamic-viscosity conversion.
pub fn resolve_liquid_state_from_kinematic(
    card: &MaterialCard,
    point: &QueryPoint,
    selection: MaterialPropertySelection,
) -> Result<ResolvedLiquidState, MaterialStatePointError> {
    resolve_liquid_state_with_viscosity(
        card,
        point,
        selection,
        KINEMATIC_VISCOSITY_PROPERTY,
        KINEMATIC_VISCOSITY_DIMS,
        true,
    )
}

fn resolve_liquid_state_with_viscosity(
    card: &MaterialCard,
    point: &QueryPoint,
    selection: MaterialPropertySelection,
    viscosity_property: &str,
    viscosity_dims: Dims,
    viscosity_is_kinematic: bool,
) -> Result<ResolvedLiquidState, MaterialStatePointError> {
    let requirements = [
        ScalarPropertyRequirement::try_new(
            DENSITY_PROPERTY,
            Density::DIMS,
            ScalarAdmissibility::StrictlyPositive,
        )?,
        ScalarPropertyRequirement::try_new(
            SPECIFIC_HEAT_CAPACITY_PROPERTY,
            SPECIFIC_HEAT_CAPACITY_DIMS,
            ScalarAdmissibility::StrictlyPositive,
        )?,
        ScalarPropertyRequirement::try_with_quantity(
            THERMAL_CONDUCTIVITY_PROPERTY,
            QuantitySpec::semantic(SemanticType::new(
                QuantityKind::ThermalConductivity,
                ValueForm::Static,
            )),
            ScalarAdmissibility::StrictlyPositive,
        )?,
        ScalarPropertyRequirement::try_new(
            viscosity_property,
            viscosity_dims,
            ScalarAdmissibility::StrictlyPositive,
        )?,
    ];
    let material = resolve_material_state_point(card, point, &requirements, selection)?;
    let value = |name| {
        material
            .property(name)
            .expect("required liquid property was resolved")
            .value_si()
    };
    let density_kg_m3 = value(DENSITY_PROPERTY);
    let specific_heat_j_kg_k = value(SPECIFIC_HEAT_CAPACITY_PROPERTY);
    let thermal_conductivity_w_m_k = value(THERMAL_CONDUCTIVITY_PROPERTY);
    let source_viscosity = value(viscosity_property);
    let dynamic_viscosity_pa_s = if viscosity_is_kinematic {
        dynamic_viscosity_from_kinematic(density_kg_m3, source_viscosity)?
    } else {
        source_viscosity
    };
    let (kinematic_viscosity_m2_s, thermal_diffusivity_m2_s, prandtl) = derived_transport(
        density_kg_m3,
        specific_heat_j_kg_k,
        thermal_conductivity_w_m_k,
        dynamic_viscosity_pa_s,
    )?;
    Ok(ResolvedLiquidState {
        material,
        density_kg_m3,
        specific_heat_j_kg_k,
        thermal_conductivity_w_m_k,
        dynamic_viscosity_pa_s,
        kinematic_viscosity_m2_s,
        thermal_diffusivity_m2_s,
        prandtl,
    })
}

fn dynamic_viscosity_from_kinematic(
    density_kg_m3: f64,
    kinematic_viscosity_m2_s: f64,
) -> Result<f64, MaterialStatePointError> {
    let dynamic_viscosity_pa_s = density_kg_m3 * kinematic_viscosity_m2_s;
    if !(dynamic_viscosity_pa_s.is_finite() && dynamic_viscosity_pa_s > 0.0) {
        return Err(MaterialStatePointError::InvalidDerived {
            quantity: "liquid dynamic viscosity from kinematic viscosity",
        });
    }
    Ok(dynamic_viscosity_pa_s)
}

fn derived_transport(
    density_kg_m3: f64,
    specific_heat_j_kg_k: f64,
    thermal_conductivity_w_m_k: f64,
    dynamic_viscosity_pa_s: f64,
) -> Result<(f64, f64, f64), MaterialStatePointError> {
    if [
        density_kg_m3,
        specific_heat_j_kg_k,
        thermal_conductivity_w_m_k,
        dynamic_viscosity_pa_s,
    ]
    .into_iter()
    .any(|value| !value.is_finite() || value <= 0.0)
    {
        return Err(MaterialStatePointError::InvalidDerived {
            quantity: "finite positive liquid source properties",
        });
    }
    let kinematic_viscosity_m2_s = dynamic_viscosity_pa_s / density_kg_m3;
    let thermal_diffusivity_m2_s =
        (thermal_conductivity_w_m_k / density_kg_m3) / specific_heat_j_kg_k;
    let prandtl = (dynamic_viscosity_pa_s / thermal_conductivity_w_m_k) * specific_heat_j_kg_k;
    if [kinematic_viscosity_m2_s, thermal_diffusivity_m2_s, prandtl]
        .into_iter()
        .any(|value| !value.is_finite() || value <= 0.0)
    {
        return Err(MaterialStatePointError::InvalidDerived {
            quantity: "finite positive liquid transport quantity",
        });
    }
    Ok((kinematic_viscosity_m2_s, thermal_diffusivity_m2_s, prandtl))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn g0_derived_transport_matches_liquid_relations() {
        let (nu, alpha, prandtl) = derived_transport(998.2, 4_182.0, 0.598, 1.002e-3).unwrap();
        assert!((nu - 1.003_806_852_334_201_6e-6).abs() < 1e-20);
        assert!((alpha - 1.432_516_358_234_875_4e-7).abs() < 1e-21);
        assert!((prandtl - 7.007_297_658_862_877_5).abs() < 1e-13);
    }

    #[test]
    fn g0_derived_transport_refuses_nonphysical_or_unrepresentable_output() {
        assert!(derived_transport(0.0, 1.0, 1.0, 1.0).is_err());
        assert!(derived_transport(f64::MAX, f64::MAX, f64::MIN_POSITIVE, 1.0).is_err());
    }

    #[test]
    fn g0_kinematic_conversion_requires_representable_dynamic_viscosity() {
        assert_eq!(
            dynamic_viscosity_from_kinematic(1_000.0, 1.0e-6).unwrap(),
            1.0e-3
        );
        assert!(dynamic_viscosity_from_kinematic(f64::MAX, f64::MAX).is_err());
    }
}
