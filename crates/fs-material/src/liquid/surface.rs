//! Source-resolved isothermal liquid momentum and free-surface coefficients.
//!
//! A liquid flow need not invent heat capacity or conductivity when it does
//! not evolve heat. Density, the explicitly chosen viscosity basis, and only
//! an active surface-tension law are resolved at one complete source point.
//! Surface tension remains conditional on that source's interface environment;
//! this adapter supplies no gas, purity, surfactant or contact-angle defaults.

use fs_matdb::{MaterialCard, QueryPoint};
use fs_qty::{Density, Dims, Frequency, QuantitySpec};
use fs_qty::semantic::QuantityKind;

use super::{
    DYNAMIC_VISCOSITY_DIMS, DYNAMIC_VISCOSITY_PROPERTY, KINEMATIC_VISCOSITY_DIMS,
    KINEMATIC_VISCOSITY_PROPERTY, dynamic_viscosity_from_kinematic,
};
use crate::liquid_temperature::LiquidViscosityBasis;
use crate::state_point::{
    DENSITY_PROPERTY, MaterialPropertySelection, MaterialStatePointError,
    ResolvedMaterialStatePoint, ScalarAdmissibility, ScalarPropertyRequirement,
    resolve_material_state_point,
};

/// Canonical source property for liquid surface tension [N/m].
pub const SURFACE_TENSION_PROPERTY: &str = "surface_tension";
/// Surface-tension dimensions [L, M, t, T, I, N] for N/m.
pub const SURFACE_TENSION_DIMS: Dims = Dims([0, 1, -2, 0, 0, 0]);

/// The scenario's explicit free-surface approximation.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum LiquidCapillarity {
    /// Omit capillary stress. Missing surface data is not interpreted as zero.
    Neglected,
    /// Use a positive sourced coefficient, frozen at the complete query point.
    /// Spatially varying temperature or composition needs a different law.
    Sourced,
}

/// Frozen source state for an isothermal, single-phase liquid free surface.
///
/// Every used claim and source condition remains in `material()`. This bundle
/// does not infer a phase diagram, an ambient gas law, wetting, or surface-energy
/// evolution. A coefficient reported for an unknown surface environment stays
/// unknown; it is not promoted to an arbitrary liquid/gas pairing.
#[derive(Clone, Debug, PartialEq)]
pub struct ResolvedLiquidSurfaceState {
    material: ResolvedMaterialStatePoint,
    viscosity_basis: LiquidViscosityBasis,
    density_kg_m3: f64,
    dynamic_viscosity_pa_s: f64,
    kinematic_viscosity_m2_s: f64,
    surface_tension_n_m: Option<f64>,
}

impl ResolvedLiquidSurfaceState {
    /// Adopt the SAME immutable bundle already resolved for other physics.
    ///
    /// Extra properties, original selections and all receipts are retained.
    /// The requested viscosity basis is the independent source; another viscosity
    /// observation is not reconciled with it or turned into a derived claim.
    /// Only validity-wide constants or positive absolute-temperature curves are
    /// admitted. Other response curves and frequency/rate/strain-restricted
    /// coefficients require a constitutive owner that evolves those coordinates.
    /// Fixed pressure, composition and interface envelopes remain source context.
    ///
    /// # Errors
    /// Refuses missing/ambiguous required properties, incorrect quantities,
    /// nonpositive values, unsupported source dependence and viscosity overflow.
    pub fn try_from_material(
        material: ResolvedMaterialStatePoint,
        viscosity_basis: LiquidViscosityBasis,
        capillarity: LiquidCapillarity,
    ) -> Result<Self, MaterialStatePointError> {
        let (viscosity_property, viscosity_dims) = viscosity_requirement(viscosity_basis);
        let density_kg_m3 = coefficient(&material, DENSITY_PROPERTY, Density::DIMS)?;
        let source_viscosity = coefficient(&material, viscosity_property, viscosity_dims)?;
        let (dynamic_viscosity_pa_s, kinematic_viscosity_m2_s) = match viscosity_basis {
            LiquidViscosityBasis::Dynamic => (source_viscosity, source_viscosity / density_kg_m3),
            LiquidViscosityBasis::Kinematic => (
                dynamic_viscosity_from_kinematic(density_kg_m3, source_viscosity)?, source_viscosity,
            ),
        };
        if !kinematic_viscosity_m2_s.is_finite() || kinematic_viscosity_m2_s <= 0.0 {
            return Err(MaterialStatePointError::InvalidDerived {
                quantity: "liquid kinematic viscosity from dynamic viscosity",
            });
        }
        let surface_tension_n_m = match capillarity {
            LiquidCapillarity::Neglected => None,
            LiquidCapillarity::Sourced => Some(coefficient(
                &material, SURFACE_TENSION_PROPERTY, SURFACE_TENSION_DIMS,
            )?),
        };
        Ok(Self {
            material, viscosity_basis, density_kg_m3, dynamic_viscosity_pa_s,
            kinematic_viscosity_m2_s, surface_tension_n_m,
        })
    }

    /// Complete source bundle, including original viscosity and interface context.
    #[must_use]
    pub const fn material(&self) -> &ResolvedMaterialStatePoint { &self.material }

    /// Source-reported viscosity convention; no fallback or synthetic claim.
    #[must_use]
    pub const fn viscosity_basis(&self) -> LiquidViscosityBasis { self.viscosity_basis }

    /// Frozen liquid density [kg/m³].
    #[must_use]
    pub const fn density_kg_m3(&self) -> f64 { self.density_kg_m3 }

    /// Dynamic viscosity [Pa s], derived with this same density for a ν source.
    #[must_use]
    pub const fn dynamic_viscosity_pa_s(&self) -> f64 { self.dynamic_viscosity_pa_s }

    /// Kinematic viscosity [m²/s] used by the momentum operator.
    #[must_use]
    pub const fn kinematic_viscosity_m2_s(&self) -> f64 { self.kinematic_viscosity_m2_s }

    /// Sourced constant surface tension [N/m], or an explicitly neglected law.
    #[must_use]
    pub const fn surface_tension_n_m(&self) -> Option<f64> { self.surface_tension_n_m }

    /// The chosen capillary approximation, retained independently of material name.
    #[must_use]
    pub const fn capillarity(&self) -> LiquidCapillarity {
        if self.surface_tension_n_m.is_some() {
            LiquidCapillarity::Sourced
        } else {
            LiquidCapillarity::Neglected
        }
    }
}

/// Resolve only the coefficients the isothermal momentum/surface operator uses.
///
/// All required claims share the same full query and explicit selection policy.
/// Source-supported pressure, temperature, composition and interface conditions
/// are retained without defaults. Choosing `Neglected` requires no σ claim;
/// choosing `Sourced` refuses missing, ambiguous or nonpositive σ.
///
/// # Errors
/// Refuses the complete bundle on source/domain/unit/selection errors and on
/// unrepresentable viscosity conversion. Thermal properties are never requested.
pub fn resolve_liquid_surface_state(
    card: &MaterialCard,
    point: &QueryPoint,
    viscosity_basis: LiquidViscosityBasis,
    capillarity: LiquidCapillarity,
    selection: MaterialPropertySelection,
) -> Result<ResolvedLiquidSurfaceState, MaterialStatePointError> {
    let (viscosity_property, viscosity_dims) = viscosity_requirement(viscosity_basis);
    let mut requirements = vec![
        ScalarPropertyRequirement::try_new(
            DENSITY_PROPERTY, Density::DIMS, ScalarAdmissibility::StrictlyPositive,
        )?,
        ScalarPropertyRequirement::try_new(
            viscosity_property, viscosity_dims, ScalarAdmissibility::StrictlyPositive,
        )?,
    ];
    if capillarity == LiquidCapillarity::Sourced {
        requirements.push(ScalarPropertyRequirement::try_new(
            SURFACE_TENSION_PROPERTY, SURFACE_TENSION_DIMS, ScalarAdmissibility::StrictlyPositive,
        )?);
    }
    let material = resolve_material_state_point(card, point, &requirements, selection)?;
    ResolvedLiquidSurfaceState::try_from_material(material, viscosity_basis, capillarity)
}

fn viscosity_requirement(basis: LiquidViscosityBasis) -> (&'static str, Dims) {
    match basis {
        LiquidViscosityBasis::Dynamic => (DYNAMIC_VISCOSITY_PROPERTY, DYNAMIC_VISCOSITY_DIMS),
        LiquidViscosityBasis::Kinematic => (KINEMATIC_VISCOSITY_PROPERTY, KINEMATIC_VISCOSITY_DIMS),
    }
}

fn coefficient(
    material: &ResolvedMaterialStatePoint,
    name: &'static str,
    dims: Dims,
) -> Result<f64, MaterialStatePointError> {
    let property = material.property(name)
        .ok_or(MaterialStatePointError::InvalidDerived { quantity: name })?;
    if property.requirement().quantity() != QuantitySpec::dimensional(dims)
        || !property.is_constant_at_fixed_temperature()
    {
        return Err(MaterialStatePointError::InvalidRequirement {
            property: name.to_owned(),
            reason: "isothermal liquid needs its exact scalar quantity and a constant or absolute-temperature-only source",
        });
    }
    // A positive sample within a shear/frequency/strain band does not certify
    // an instantaneous Newtonian coefficient outside that band during flow.
    for (axis, _) in &property.answer().receipt.query_point {
        if property.answer().evidence.model.validity.bound(axis).is_none() { continue; }
        let quantity = property.answer().receipt.axis_quantities.get(axis).copied();
        let kind = quantity.and_then(QuantitySpec::semantic_type).map(|ty| ty.kind());
        let dynamic_kind = matches!(kind, Some(
            QuantityKind::Frequency(_) | QuantityKind::AngularVelocity(_)
                | QuantityKind::Strain { .. } | QuantityKind::Stress
        ));
        let rate_dimension = quantity.is_some_and(|q| q.dims() == Frequency::DIMS);
        let legacy_dynamic_axis = kind.is_none() && [
            "omega", "frequency", "strain", "epsilon", "eps", "shear_strain",
            "strain_rate", "shear_rate", "plastic_strain", "equivalent_plastic_strain", "stress",
        ].contains(&axis.as_str());
        if dynamic_kind || rate_dimension || legacy_dynamic_axis {
            return Err(MaterialStatePointError::InvalidRequirement {
                property: name.to_owned(),
                reason: "frequency, rate, strain or stress bounded liquid data needs an evolving constitutive owner",
            });
        }
    }
    let value = property.value_si();
    if !value.is_finite() || value <= 0.0 {
        return Err(MaterialStatePointError::OutsideConsumerDomain {
            property: name.to_owned(), value, admissibility: ScalarAdmissibility::StrictlyPositive,
        });
    }
    Ok(value)
}
