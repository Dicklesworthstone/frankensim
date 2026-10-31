//! Source-resolved isothermal liquid momentum and free-surface coefficients.
//!
//! A liquid flow need not invent heat capacity or conductivity when it does
//! not evolve heat. Density, the explicitly chosen viscosity basis, and only
//! an active surface-tension law are resolved at one complete source point.
//! Surface tension remains conditional on that source's interface environment;
//! this adapter supplies no gas, purity, surfactant or contact-angle defaults.

/// Reversible adhesion work driving the same source-bound equilibrium wall.
pub mod adhesion;
/// Shared sourced contact-line friction and advancing/receding hysteresis.
pub mod contact_line;

use fs_matdb::{InterfaceSystemCard, MaterialCard, QueryPoint, SurfaceSpec, SystemContext};
use fs_qty::semantic::QuantityKind;
use fs_qty::{Density, Dims, Frequency, QuantitySpec};

use super::{
    DYNAMIC_VISCOSITY_DIMS, DYNAMIC_VISCOSITY_PROPERTY, KINEMATIC_VISCOSITY_DIMS,
    KINEMATIC_VISCOSITY_PROPERTY, dynamic_viscosity_from_kinematic,
};
use crate::liquid_temperature::LiquidViscosityBasis;
use crate::state_point::{
    DENSITY_PROPERTY, MaterialPropertySelection, MaterialStatePointError,
    ResolvedInterfaceStatePoint, ResolvedMaterialStatePoint, ResolvedScalarProperty,
    ScalarAdmissibility, ScalarPropertyRequirement, resolve_interface_state_point,
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

/// Source state for a single-phase liquid free surface.
///
/// Every used claim and source condition remains in `material()`. This bundle
/// does not infer a phase diagram, an ambient gas law, wetting, or surface-energy
/// evolution. A coefficient reported for an unknown surface environment stays
/// unknown; it is not promoted to an arbitrary liquid/gas pairing.
/// The ordinary constructor admits frozen coefficients. The explicit rheology
/// constructor instead requires an evolving temperature/shear owner; consumers
/// must refuse that bundle when they cannot supply that evolution.
#[derive(Clone, Debug, PartialEq)]
pub struct ResolvedLiquidSurfaceState {
    material: ResolvedMaterialStatePoint,
    viscosity_basis: LiquidViscosityBasis,
    density_kg_m3: f64,
    dynamic_viscosity_pa_s: f64,
    kinematic_viscosity_m2_s: f64,
    surface_tension_n_m: Option<f64>,
    rheology_axes: Option<(String, String)>,
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
                dynamic_viscosity_from_kinematic(density_kg_m3, source_viscosity)?,
                source_viscosity,
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
                &material,
                SURFACE_TENSION_PROPERTY,
                SURFACE_TENSION_DIMS,
            )?),
        };
        Ok(Self {
            material,
            viscosity_basis,
            density_kg_m3,
            dynamic_viscosity_pa_s,
            kinematic_viscosity_m2_s,
            surface_tension_n_m,
            rheology_axes: None,
        })
    }

    /// Admit a source for an owner that evolves temperature and physical shear
    /// rate `sqrt(2 D:D)` together. This does NOT grant frozen isothermal use:
    /// [`Self::rheology_axes`] remains set and must be consumed by a runtime
    /// constitutive solve. Capillarity is explicitly neglected on this path.
    ///
    /// Only these two response axes may vary. The rate must have the explicit
    /// dimensional quantity 1/s, not cyclic/angular frequency; the temperature
    /// must be positive absolute kelvin. Other dynamic applicability refuses.
    /// Every original property and selected claim remains in the bundle.
    pub fn try_for_shear_rate_evolution(
        material: ResolvedMaterialStatePoint,
        viscosity_basis: LiquidViscosityBasis,
        temperature_axis: &str,
        shear_rate_axis: &str,
    ) -> Result<Self, MaterialStatePointError> {
        Self::try_for_shear_rate_evolution_with_capillarity(
            material,
            viscosity_basis,
            LiquidCapillarity::Neglected,
            temperature_axis,
            shear_rate_axis,
        )
    }

    /// Admit T/shear evolution while retaining the scenario's capillary choice.
    ///
    /// All active coefficients come from the same complete source bundle.
    /// The receiving runtime must own every admitted response and explicitly
    /// refuse any capillary dependence it cannot evolve. This does not grant
    /// frozen use or supply surface-energy, Marangoni or wetting dynamics.
    ///
    /// # Errors
    /// Refuses the same unsupported axes and invalid momentum coefficients as
    /// the noncapillary entry point, and missing or invalid active surface data.
    pub fn try_for_shear_rate_evolution_with_capillarity(
        material: ResolvedMaterialStatePoint,
        viscosity_basis: LiquidViscosityBasis,
        capillarity: LiquidCapillarity,
        temperature_axis: &str,
        shear_rate_axis: &str,
    ) -> Result<Self, MaterialStatePointError> {
        let quantities = material.axis_quantities();
        let temperature_quantity = quantities.get(temperature_axis).copied();
        let temperature_ok = temperature_quantity
            == Some(QuantitySpec::semantic(fs_qty::semantic::SemanticType::new(
                QuantityKind::AbsoluteTemperature,
                fs_qty::semantic::ValueForm::Static,
            )))
            || temperature_axis == "T"
                && (temperature_quantity.is_none()
                    || temperature_quantity
                        == Some(QuantitySpec::dimensional(fs_qty::Temperature::DIMS)));
        if temperature_axis == shear_rate_axis
            || !temperature_ok
            || !material
                .query_point()
                .iter()
                .find(|(axis, _)| axis == temperature_axis)
                .is_some_and(|(_, t)| t.is_finite() && *t > 0.0)
            || quantities.get(shear_rate_axis) != Some(&QuantitySpec::dimensional(Frequency::DIMS))
            || !material
                .query_point()
                .iter()
                .find(|(axis, _)| axis == shear_rate_axis)
                .is_some_and(|(_, rate)| rate.is_finite() && *rate >= 0.0)
        {
            return Err(MaterialStatePointError::InvalidDerived {
                quantity: "liquid rheology requires named absolute temperature and physical shear-rate axes",
            });
        }
        for property in material.properties() {
            if property.source_axes().iter().any(|(axis, dims)| {
                !(axis == temperature_axis && *dims == fs_qty::Temperature::DIMS
                    || axis == shear_rate_axis && *dims == Frequency::DIMS)
            }) {
                return Err(MaterialStatePointError::InvalidRequirement {
                    property: property.requirement().name().to_owned(),
                    reason: "liquid rheology evolves only its declared temperature and shear rate",
                });
            }
            require_no_unowned_dynamic_axes(property, Some(shear_rate_axis))?;
        }
        let (name, dims) = viscosity_requirement(viscosity_basis);
        let value = |name: &'static str, dims| -> Result<f64, MaterialStatePointError> {
            let property = material
                .property(name)
                .ok_or(MaterialStatePointError::InvalidDerived { quantity: name })?;
            if property.requirement().quantity() != QuantitySpec::dimensional(dims)
                || !property.value_si().is_finite()
                || property.value_si() <= 0.0
            {
                return Err(MaterialStatePointError::InvalidDerived { quantity: name });
            }
            Ok(property.value_si())
        };
        let density_kg_m3 = value(DENSITY_PROPERTY, Density::DIMS)?;
        let source_viscosity = value(name, dims)?;
        let (dynamic_viscosity_pa_s, kinematic_viscosity_m2_s) = match viscosity_basis {
            LiquidViscosityBasis::Dynamic => (source_viscosity, source_viscosity / density_kg_m3),
            LiquidViscosityBasis::Kinematic => (
                dynamic_viscosity_from_kinematic(density_kg_m3, source_viscosity)?,
                source_viscosity,
            ),
        };
        if !kinematic_viscosity_m2_s.is_finite() || kinematic_viscosity_m2_s <= 0.0 {
            return Err(MaterialStatePointError::InvalidDerived {
                quantity: "liquid rheological kinematic viscosity",
            });
        }
        let surface_tension_n_m = match capillarity {
            LiquidCapillarity::Neglected => None,
            LiquidCapillarity::Sourced => {
                Some(value(SURFACE_TENSION_PROPERTY, SURFACE_TENSION_DIMS)?)
            }
        };
        Ok(Self {
            material,
            viscosity_basis,
            density_kg_m3,
            dynamic_viscosity_pa_s,
            kinematic_viscosity_m2_s,
            surface_tension_n_m,
            rheology_axes: Some((temperature_axis.to_owned(), shear_rate_axis.to_owned())),
        })
    }

    /// Required evolving coordinates `(temperature, physical shear rate)`.
    /// `Some` is a prohibition on frozen use, not a grant of static authority.
    #[must_use]
    pub fn rheology_axes(&self) -> Option<(&str, &str)> {
        self.rheology_axes
            .as_ref()
            .map(|(temperature, rate)| (temperature.as_str(), rate.as_str()))
    }

    /// Complete source bundle, including original viscosity and interface context.
    #[must_use]
    pub const fn material(&self) -> &ResolvedMaterialStatePoint {
        &self.material
    }

    /// Source-reported viscosity convention; no fallback or synthetic claim.
    #[must_use]
    pub const fn viscosity_basis(&self) -> LiquidViscosityBasis {
        self.viscosity_basis
    }

    /// Liquid density at this reference source point [kg/m³].
    #[must_use]
    pub const fn density_kg_m3(&self) -> f64 {
        self.density_kg_m3
    }

    /// Dynamic viscosity [Pa s], derived with this same density for a ν source.
    #[must_use]
    pub const fn dynamic_viscosity_pa_s(&self) -> f64 {
        self.dynamic_viscosity_pa_s
    }

    /// Kinematic viscosity [m²/s] used by the momentum operator.
    #[must_use]
    pub const fn kinematic_viscosity_m2_s(&self) -> f64 {
        self.kinematic_viscosity_m2_s
    }

    /// Sourced constant surface tension [N/m], or an explicitly neglected law.
    #[must_use]
    pub const fn surface_tension_n_m(&self) -> Option<f64> {
        self.surface_tension_n_m
    }

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
            DENSITY_PROPERTY,
            Density::DIMS,
            ScalarAdmissibility::StrictlyPositive,
        )?,
        ScalarPropertyRequirement::try_new(
            viscosity_property,
            viscosity_dims,
            ScalarAdmissibility::StrictlyPositive,
        )?,
    ];
    if capillarity == LiquidCapillarity::Sourced {
        requirements.push(ScalarPropertyRequirement::try_new(
            SURFACE_TENSION_PROPERTY,
            SURFACE_TENSION_DIMS,
            ScalarAdmissibility::StrictlyPositive,
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
    let property = material
        .property(name)
        .ok_or(MaterialStatePointError::InvalidDerived { quantity: name })?;
    require_frozen_coefficient(property, name, dims)?;
    let value = property.value_si();
    if !value.is_finite() || value <= 0.0 {
        return Err(MaterialStatePointError::OutsideConsumerDomain {
            property: name.to_owned(),
            value,
            admissibility: ScalarAdmissibility::StrictlyPositive,
        });
    }
    Ok(value)
}

// Momentum, capillarity and equilibrium wetting share the same rule for a
// frozen instantaneous coefficient; a valid dynamic sample is insufficient.
fn require_frozen_coefficient(
    property: &ResolvedScalarProperty,
    name: &'static str,
    dims: Dims,
) -> Result<(), MaterialStatePointError> {
    if property.requirement().quantity() != QuantitySpec::dimensional(dims) {
        return Err(MaterialStatePointError::InvalidRequirement {
            property: name.to_owned(),
            reason: "frozen liquid coefficient needs its exact scalar quantity and a constant or absolute-temperature-only source",
        });
    }
    require_frozen_coefficient_support(property)
}

/// Admit support for an instantaneous coefficient evaluated at fixed temperature.
///
/// Constants and absolute-temperature-only curves are allowed. Source bounds in
/// frequency, strain, stress, rate or flow coordinates need a constitutive owner
/// that actually evolves those coordinates; they cannot be silently frozen.
/// Pressure, composition and other fixed external conditions remain context.
///
/// This checks SUPPORT ONLY. The caller must separately require the exact
/// quantity, value range and material/state assignment. In particular, semantic
/// thermal-conductivity quantities are retained without conversion to a merely
/// dimensional key. This does not establish source accuracy or full-span cover.
///
/// # Errors
/// Refuses a non-temperature response curve or an unevolved dynamic source bound.
pub fn require_frozen_coefficient_support(
    property: &ResolvedScalarProperty,
) -> Result<(), MaterialStatePointError> {
    let name = property.requirement().name();
    if !property.is_constant_at_fixed_temperature() {
        return Err(MaterialStatePointError::InvalidRequirement {
            property: name.to_owned(),
            reason: "instantaneous coefficient needs a scalar constant or absolute-temperature-only source",
        });
    }
    require_no_unowned_dynamic_axes(property, None)
}

fn require_no_unowned_dynamic_axes(
    property: &ResolvedScalarProperty,
    shear_rate_axis: Option<&str>,
) -> Result<(), MaterialStatePointError> {
    let name = property.requirement().name();
    for (axis, _) in &property.answer().receipt.query_point {
        if shear_rate_axis == Some(axis.as_str()) {
            continue;
        }
        if property
            .answer()
            .evidence
            .model
            .validity
            .bound(axis)
            .is_none()
        {
            continue;
        }
        let quantity = property.answer().receipt.axis_quantities.get(axis).copied();
        let kind = quantity
            .and_then(QuantitySpec::semantic_type)
            .map(|ty| ty.kind());
        let dynamic_kind = matches!(
            kind,
            Some(
                QuantityKind::Frequency(_)
                    | QuantityKind::AngularVelocity(_)
                    | QuantityKind::Strain { .. }
                    | QuantityKind::Stress
            )
        );
        let dynamic_dimension = quantity
            .is_some_and(|q| q.dims() == Frequency::DIMS || q.dims() == Dims([1, 0, -1, 0, 0, 0]));
        let legacy_dynamic_axis = kind.is_none()
            && [
                "omega",
                "frequency",
                "strain",
                "epsilon",
                "eps",
                "shear_strain",
                "strain_rate",
                "shear_rate",
                "plastic_strain",
                "equivalent_plastic_strain",
                "stress",
                "contact_line_speed",
                "velocity",
                "capillary_number",
            ]
            .contains(&axis.as_str());
        if dynamic_kind || dynamic_dimension || legacy_dynamic_axis {
            return Err(MaterialStatePointError::InvalidRequirement {
                property: name.to_owned(),
                reason: "frequency, rate, strain, stress or flow bounded data needs an evolving constitutive owner",
            });
        }
    }
    Ok(())
}

/// Canonical equilibrium contact angle through the liquid, in radians.
pub const EQUILIBRIUM_CONTACT_ANGLE_PROPERTY: &str = "equilibrium-contact-angle";

/// One sourced equilibrium substrate/liquid response in its declared context.
///
/// Surface A is the substrate; surface B is the liquid. The caller assigns the
/// named environment to the physical scene: retaining that context does not
/// resolve gas composition or establish equivalence with a separate gas EOS.
/// This response supplies ONLY the angle. Liquid/gas tension remains in the
/// receiving `ResolvedLiquidSurfaceState`, including its source conditions.
#[derive(Clone, Debug, PartialEq)]
pub struct ResolvedLiquidWettingState {
    resolved: ResolvedInterfaceStatePoint,
    surfaces: [SurfaceSpec; 2],
    context: SystemContext,
    contact_angle_radians: f64,
    young_dupre: Option<adhesion::YoungDupreSource>,
}

impl ResolvedLiquidWettingState {
    /// Adopt the same complete source bundle used by other interface laws.
    /// The original card, ordered surfaces and equilibrium-angle restrictions
    /// are checked without selecting another claim or discarding other data.
    pub fn try_from_resolved(
        card: &InterfaceSystemCard,
        resolved: ResolvedInterfaceStatePoint,
    ) -> Result<Self, MaterialStatePointError> {
        if resolved.card_identity() != card.content_hash() {
            return Err(MaterialStatePointError::InvalidDerived {
                quantity: "wetting interface requires its original ordered surface card",
            });
        }
        let angle = resolved
            .property(EQUILIBRIUM_CONTACT_ANGLE_PROPERTY)
            .ok_or_else(|| MaterialStatePointError::InvalidRequirement {
                property: EQUILIBRIUM_CONTACT_ANGLE_PROPERTY.into(),
                reason: "missing or ambiguous equilibrium contact angle",
            })?;
        require_frozen_coefficient(angle, EQUILIBRIUM_CONTACT_ANGLE_PROPERTY, Dims([0; 6]))?;
        let contact_angle_radians = angle.value_si();
        if !(contact_angle_radians > 0.0 && contact_angle_radians < core::f64::consts::PI) {
            return Err(MaterialStatePointError::InvalidDerived {
                quantity: "equilibrium contact angle must lie strictly between zero and pi",
            });
        }
        Ok(Self {
            resolved,
            surfaces: [card.surface_a().clone(), card.surface_b().clone()],
            context: card.context().clone(),
            contact_angle_radians,
            young_dupre: None,
        })
    }

    /// Refresh the entire source bundle and rebuild the equilibrium wall law.
    /// The original phase/process assignments remain fixed; this does not
    /// evolve wetting hysteresis or a moving contact line.
    pub fn at_query_point(
        &self,
        card: &InterfaceSystemCard,
        point: &QueryPoint,
    ) -> Result<Self, MaterialStatePointError> {
        if self.young_dupre.is_some() {
            return Err(MaterialStatePointError::InvalidDerived {
                quantity: "adhesion-derived angle refresh also requires the current sourced liquid tension",
            });
        }
        Self::try_from_resolved(card, self.resolved.at_query_point(card, point)?)
    }

    /// The original selected interface claim and its complete usage receipt.
    #[must_use]
    pub const fn resolved(&self) -> &ResolvedInterfaceStatePoint {
        &self.resolved
    }

    /// Ordered substrate/liquid identities, including their preparation frames.
    #[must_use]
    pub const fn surfaces(&self) -> &[SurfaceSpec; 2] {
        &self.surfaces
    }

    /// Source-reported environment, medium, third body and history.
    #[must_use]
    pub const fn context(&self) -> &SystemContext {
        &self.context
    }

    /// Equilibrium angle through liquid, strictly between zero and pi.
    /// `angle_origin` distinguishes a source claim from a Young--Dupre prediction.
    #[must_use]
    pub const fn contact_angle_radians(&self) -> f64 {
        self.contact_angle_radians
    }

    /// Require the declared liquid identity and SAME complete source point.
    ///
    /// This checks the surface's phase/process/revision identity and all query
    /// values/descriptors. SurfaceSpec contains no bulk-card content hash; this
    /// check does not claim stronger binding than that shared identity contract.
    pub fn require_liquid_state(
        &self,
        liquid: &ResolvedLiquidSurfaceState,
    ) -> Result<(), MaterialStatePointError> {
        let bulk = liquid.material();
        if &self.surfaces[1].material != bulk.material()
            || self.resolved.query_point() != bulk.query_point()
            || self.resolved.axis_quantities() != bulk.axis_quantities()
        {
            return Err(MaterialStatePointError::InvalidDerived {
                quantity: "wetting interface and liquid material/state assignment mismatch",
            });
        }
        if let Some(source) = &self.young_dupre {
            source.require_liquid(liquid)?;
        }
        Ok(())
    }
}

/// Resolve only the equilibrium angle needed by a geometric wall closure.
///
/// Missing equilibrium data never falls back to advancing/receding endpoints
/// or a neutral wall. Constants or absolute-temperature-only sources may be
/// frozen at this point; dynamic rate/strain/frequency support is not silently
/// broadened into a static response. Complete interface identity and context
/// remain accessible. This supplies no hysteresis, moving-contact-line law,
/// thermocapillary traction or experimental accuracy claim.
///
/// # Errors
/// Refuses source/selection/unit/domain failures, angles outside (0, pi), and
/// source dependence requiring an evolving constitutive owner.
pub fn resolve_liquid_wetting_state(
    card: &InterfaceSystemCard,
    point: &QueryPoint,
    selection: MaterialPropertySelection,
) -> Result<ResolvedLiquidWettingState, MaterialStatePointError> {
    let requirement = ScalarPropertyRequirement::try_new(
        EQUILIBRIUM_CONTACT_ANGLE_PROPERTY,
        Dims([0; 6]),
        ScalarAdmissibility::OpenInterval {
            lower: 0.0,
            upper: core::f64::consts::PI,
        },
    )?;
    let resolved = resolve_interface_state_point(card, point, &[requirement], selection)?;
    ResolvedLiquidWettingState::try_from_resolved(card, resolved)
}
