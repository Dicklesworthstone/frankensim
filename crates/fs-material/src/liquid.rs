//! Source-resolved single-phase liquid thermal and transport tuple.
//!
//! The caller's material card and complete query point define the liquid-phase
//! conditions. This adapter resolves four required bulk claims and, when
//! requested, an independent isobaric volumetric-expansion claim. It derives
//! kinematic viscosity, thermal diffusivity, and Prandtl number, but supplies
//! no equation of state, phase boundary, concentration conversion,
//! chemistry-name default, corrosion model, or experimental-accuracy claim.

pub mod surface;
mod thermal;
pub use thermal::ResolvedLiquidExpansion;

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
/// Canonical source property for instantaneous isobaric volumetric expansion.
/// This is `beta = -(1/rho) (d rho/dT)_p` at fixed composition, not a linear
/// expansion coefficient or a mean coefficient over a temperature interval.
pub const VOLUMETRIC_EXPANSION_PROPERTY: &str = "volumetric_thermal_expansion_coefficient";
/// Isobaric volumetric-expansion dimensions [1/K].
pub const VOLUMETRIC_EXPANSION_DIMS: Dims = Dims([0, 0, 0, -1, 0, 0]);

/// Exact source keys assigned to density, heat capacity, conductivity and viscosity.
/// Source vocabulary is preserved; selecting roles never rewrites material claims.
#[derive(Clone, Debug, PartialEq)]
pub struct LiquidPropertyRequirements {
    properties: [ScalarPropertyRequirement; 4],
    viscosity_basis: crate::liquid_temperature::LiquidViscosityBasis,
}

impl LiquidPropertyRequirements {
    /// Assign four distinct source keys to their physical liquid transport roles.
    /// Conductivity requires its static semantic type; the other roles require
    /// dimension-only scalars. Numerical admissibility remains source-explicit.
    pub fn try_new(
        density: ScalarPropertyRequirement,
        specific_heat: ScalarPropertyRequirement,
        thermal_conductivity: ScalarPropertyRequirement,
        viscosity: ScalarPropertyRequirement,
        viscosity_basis: crate::liquid_temperature::LiquidViscosityBasis,
    ) -> Result<Self, MaterialStatePointError> {
        let viscosity_dims = match viscosity_basis {
            crate::liquid_temperature::LiquidViscosityBasis::Dynamic => DYNAMIC_VISCOSITY_DIMS,
            crate::liquid_temperature::LiquidViscosityBasis::Kinematic => KINEMATIC_VISCOSITY_DIMS,
        };
        let properties = [density, specific_heat, thermal_conductivity, viscosity];
        let quantities = [
            QuantitySpec::dimensional(Density::DIMS),
            QuantitySpec::dimensional(SPECIFIC_HEAT_CAPACITY_DIMS),
            QuantitySpec::semantic(SemanticType::new(
                QuantityKind::ThermalConductivity,
                ValueForm::Static,
            )),
            QuantitySpec::dimensional(viscosity_dims),
        ];
        for (index, (property, quantity)) in properties.iter().zip(quantities).enumerate() {
            if property.quantity() != quantity
                || property.hardness_test().is_some()
                || property.elastic_component().is_some()
                || property.strain_component().is_some()
                || property.stress_component().is_some()
                || properties[..index]
                    .iter()
                    .any(|prior| prior.name() == property.name())
            {
                return Err(MaterialStatePointError::InvalidRequirement {
                    property: property.name().to_owned(),
                    reason: "liquid transport requires distinct exact scalar keys of the declared physical quantities",
                });
            }
        }
        Ok(Self {
            properties,
            viscosity_basis,
        })
    }

    /// Canonical source names with explicitly selected dynamic or kinematic viscosity.
    pub fn canonical(
        viscosity_basis: crate::liquid_temperature::LiquidViscosityBasis,
    ) -> Result<Self, MaterialStatePointError> {
        let (name, dims) = match viscosity_basis {
            crate::liquid_temperature::LiquidViscosityBasis::Dynamic => {
                (DYNAMIC_VISCOSITY_PROPERTY, DYNAMIC_VISCOSITY_DIMS)
            }
            crate::liquid_temperature::LiquidViscosityBasis::Kinematic => {
                (KINEMATIC_VISCOSITY_PROPERTY, KINEMATIC_VISCOSITY_DIMS)
            }
        };
        let [density, specific_heat, conductivity, viscosity] = liquid_requirements(name, dims)?;
        Self::try_new(
            density,
            specific_heat,
            conductivity,
            viscosity,
            viscosity_basis,
        )
    }

    /// Original requirements in density, heat capacity, conductivity, viscosity order.
    pub const fn properties(&self) -> &[ScalarPropertyRequirement; 4] {
        &self.properties
    }

    /// Original source viscosity basis; there is no fallback to another claim.
    pub const fn viscosity_basis(&self) -> crate::liquid_temperature::LiquidViscosityBasis {
        self.viscosity_basis
    }
}

/// One source-resolved liquid transport state at its complete query point.
///
/// The retained material bundle contains every selected property-use
/// receipt and the exact query identity. It does not certify that the source
/// conditions remain single phase beyond the claims' own validity envelopes.
#[derive(Clone, Debug, PartialEq)]
pub struct ResolvedLiquidState {
    material: ResolvedMaterialStatePoint,
    density_expansion: Option<ResolvedLiquidExpansion>,
    requirements: LiquidPropertyRequirements,
    density_kg_m3: f64,
    specific_heat_j_kg_k: f64,
    thermal_conductivity_w_m_k: f64,
    dynamic_viscosity_pa_s: f64,
    kinematic_viscosity_m2_s: f64,
    thermal_diffusivity_m2_s: f64,
    prandtl: f64,
    volumetric_expansion_per_k: Option<f64>,
}

impl ResolvedLiquidState {
    /// Adopt an already resolved complete bundle without selecting new claims.
    ///
    /// All original properties, query meanings and receipts are retained. The
    /// four liquid requirements must match this consumer exactly. An optional
    /// volumetric-expansion claim must match its finite signed requirement too;
    /// additional coupled properties refresh with them at every later query.
    ///
    /// # Errors
    /// Refuses missing or incompatible requirements and invalid derived transport.
    pub fn try_from_material(
        material: ResolvedMaterialStatePoint,
        viscosity_basis: crate::liquid_temperature::LiquidViscosityBasis,
    ) -> Result<Self, MaterialStatePointError> {
        Self::try_from_material_with_requirements(
            material,
            &LiquidPropertyRequirements::canonical(viscosity_basis)?,
        )
    }

    /// Adopt the retained bundle with explicit original source keys and roles.
    pub fn try_from_material_with_requirements(
        material: ResolvedMaterialStatePoint,
        requirements: &LiquidPropertyRequirements,
    ) -> Result<Self, MaterialStatePointError> {
        for required in requirements.properties() {
            if material
                .property(required.name())
                .is_none_or(|property| property.requirement() != required)
            {
                return Err(MaterialStatePointError::InvalidRequirement {
                    property: required.name().to_owned(),
                    reason: "retained liquid bundle must contain each exact transport requirement",
                });
            }
        }
        liquid_from_resolved_material(material, requirements)
    }

    /// Refresh every transport property at a new complete physical query while
    /// retaining the original material card, source claims and viscosity basis.
    /// Temperature, pressure, shear rate and composition coordinates may evolve
    /// together, provided their names/conventions and source support are retained.
    /// In particular, a kinematic-viscosity source is converted with the NEW
    /// density. No derived transport coefficient remains at the reference state.
    ///
    /// # Errors
    /// Refuses any source-bundle refresh or derived-transport failure atomically.
    /// A retained density-expansion path also requires its original nonthermal
    /// context and a temperature with a unique tangent inside that path.
    /// This evaluates a constitutive state; it does not supply a fluid solver,
    /// phase stability, or a path-dependent rheological history.
    pub fn at_query_point(
        &self,
        card: &MaterialCard,
        point: &QueryPoint,
    ) -> Result<Self, MaterialStatePointError> {
        let mut updated = liquid_from_resolved_material(
            self.material.at_query_point(card, point)?,
            &self.requirements,
        )?;
        if let Some(expansion) = &self.density_expansion {
            updated.density_expansion =
                Some(expansion.at_material(&updated.material, updated.density_kg_m3)?);
        }
        Ok(updated)
    }

    /// The complete evidence-bearing source bundle and all retained receipts.
    #[must_use]
    pub const fn material(&self) -> &ResolvedMaterialStatePoint {
        &self.material
    }

    /// Exact original source keys, physical roles and viscosity basis.
    pub const fn requirements(&self) -> &LiquidPropertyRequirements {
        &self.requirements
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

    /// Source-resolved isobaric volumetric expansion [1/K] at this same point.
    /// An explicitly admitted density path supplies its analytic local tangent;
    /// otherwise the independently selected expansion claim supplies the value.
    /// Missing source information is `None`, distinct from known zero or
    /// negative expansion. An isolated density value supplies no derivative.
    #[must_use]
    pub const fn volumetric_expansion_per_k(&self) -> Option<f64> {
        match &self.density_expansion {
            Some(expansion) => Some(expansion.volumetric_expansion_per_k()),
            None => self.volumetric_expansion_per_k,
        }
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
    resolve_liquid_state_with_requirements(
        card,
        point,
        &LiquidPropertyRequirements::canonical(
            crate::liquid_temperature::LiquidViscosityBasis::Dynamic,
        )?,
        selection,
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
    resolve_liquid_state_with_requirements(
        card,
        point,
        &LiquidPropertyRequirements::canonical(
            crate::liquid_temperature::LiquidViscosityBasis::Kinematic,
        )?,
        selection,
    )
}

/// Resolve transport and isobaric volumetric expansion from five sourced claims.
///
/// Every property is selected at the SAME complete query and retained in the
/// original liquid bundle. The explicitly selected viscosity basis never falls
/// back to the other basis. Expansion may be positive, zero or negative; it is
/// neither inferred from density samples nor replaced by ideal-gas `1/T`.
/// This admits the declared source laws, not an independent thermodynamic
/// consistency or liquid-phase-stability certification of those laws.
///
/// # Errors
/// Refuses a missing, ambiguous, wrong-quantity or unsupported expansion claim
/// as well as every failure of the four-property transport resolution.
pub fn resolve_liquid_state_with_expansion(
    card: &MaterialCard,
    point: &QueryPoint,
    viscosity_basis: crate::liquid_temperature::LiquidViscosityBasis,
    selection: MaterialPropertySelection,
) -> Result<ResolvedLiquidState, MaterialStatePointError> {
    let liquid = LiquidPropertyRequirements::canonical(viscosity_basis)?;
    let mut requirements = liquid.properties().to_vec();
    requirements.push(expansion_requirement()?);
    let material = resolve_material_state_point(card, point, &requirements, selection)?;
    liquid_from_resolved_material(material, &liquid)
}

/// Resolve transport using explicitly selected original source keys.
/// All claims are queried together and retained without aliases or source mutation.
pub fn resolve_liquid_state_with_requirements(
    card: &MaterialCard,
    point: &QueryPoint,
    requirements: &LiquidPropertyRequirements,
    selection: MaterialPropertySelection,
) -> Result<ResolvedLiquidState, MaterialStatePointError> {
    let material = resolve_material_state_point(card, point, requirements.properties(), selection)?;
    liquid_from_resolved_material(material, requirements)
}

fn liquid_requirements(
    viscosity_property: &str,
    viscosity_dims: Dims,
) -> Result<[ScalarPropertyRequirement; 4], MaterialStatePointError> {
    Ok([
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
    ])
}

// Only checked liquid admission and exact-source refresh enter here, so all
// four required complete keys have already passed the same consumer checks.
fn liquid_from_resolved_material(
    material: ResolvedMaterialStatePoint,
    requirements: &LiquidPropertyRequirements,
) -> Result<ResolvedLiquidState, MaterialStatePointError> {
    let volumetric_expansion_per_k = sourced_expansion(&material)?;
    let viscosity_is_kinematic = requirements.viscosity_basis()
        == crate::liquid_temperature::LiquidViscosityBasis::Kinematic;
    let value = |name| {
        material
            .property(name)
            .expect("required liquid property was resolved")
            .value_si()
    };
    let density_kg_m3 = value(requirements.properties()[0].name());
    let specific_heat_j_kg_k = value(requirements.properties()[1].name());
    let thermal_conductivity_w_m_k = value(requirements.properties()[2].name());
    let source_viscosity = value(requirements.properties()[3].name());
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
        density_expansion: None,
        requirements: requirements.clone(),
        density_kg_m3,
        specific_heat_j_kg_k,
        thermal_conductivity_w_m_k,
        dynamic_viscosity_pa_s,
        kinematic_viscosity_m2_s,
        thermal_diffusivity_m2_s,
        prandtl,
        volumetric_expansion_per_k,
    })
}

fn expansion_requirement() -> Result<ScalarPropertyRequirement, MaterialStatePointError> {
    ScalarPropertyRequirement::try_new(
        VOLUMETRIC_EXPANSION_PROPERTY,
        VOLUMETRIC_EXPANSION_DIMS,
        ScalarAdmissibility::Finite,
    )
}

fn sourced_expansion(
    material: &ResolvedMaterialStatePoint,
) -> Result<Option<f64>, MaterialStatePointError> {
    let mut candidates = material
        .properties()
        .iter()
        .filter(|property| property.requirement().name() == VOLUMETRIC_EXPANSION_PROPERTY);
    let Some(property) = candidates.next() else {
        return Ok(None);
    };
    if candidates.next().is_some() || property.requirement() != &expansion_requirement()? {
        return Err(MaterialStatePointError::InvalidRequirement {
            property: VOLUMETRIC_EXPANSION_PROPERTY.to_owned(),
            reason: "retained liquid expansion must be one exact finite isobaric volumetric coefficient",
        });
    }
    Ok(Some(property.value_si()))
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
