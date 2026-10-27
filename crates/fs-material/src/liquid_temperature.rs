//! Temperature updates of one source-resolved liquid, shared by physics owners.
//!
//! The initial resolution selects the material claims once. Subsequent states
//! change only the declared absolute-temperature coordinate and pin those same
//! claims. Crossing their support refuses; a new temperature cannot silently
//! choose another formulation, phase, pressure, composition or source table.

use fs_matdb::{MaterialCard, QueryPoint};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{QuantitySpec, Temperature};

use crate::liquid::{
    ResolvedLiquidState, resolve_liquid_state, resolve_liquid_state_from_kinematic,
};
use crate::state_point::{
    MaterialPropertySelection, MaterialStatePointError, ResolvedMaterialStatePoint,
};

/// Which viscosity quantity the source actually reports. No implicit fallback.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum LiquidViscosityBasis {
    /// Dynamic viscosity [Pa s].
    Dynamic,
    /// Kinematic viscosity [m2/s], converted using density at the SAME state.
    Kinematic,
}

impl LiquidViscosityBasis {
    fn resolve(
        self, card: &MaterialCard, point: &QueryPoint, selection: MaterialPropertySelection,
    ) -> Result<ResolvedLiquidState, MaterialStatePointError> {
        match self {
            Self::Dynamic => resolve_liquid_state(card, point, selection),
            Self::Kinematic => resolve_liquid_state_from_kinematic(card, point, selection),
        }
    }

    #[cfg(test)]
    fn property(self) -> &'static str {
        match self {
            Self::Dynamic => crate::liquid::DYNAMIC_VISCOSITY_PROPERTY,
            Self::Kinematic => crate::liquid::KINEMATIC_VISCOSITY_PROPERTY,
        }
    }
}

/// A fixed-context liquid law, not a frozen set of transport coefficients.
///
/// Owns its initial query and resolution and borrows the immutable source card.
/// The original source authority survives each evaluation. This does not assert
/// liquid phase stability, extrapolate, integrate energy, or evolve pressure.
#[derive(Clone, Debug)]
pub struct LiquidTemperatureLaw<'a> {
    card: &'a MaterialCard,
    reference_point: QueryPoint,
    reference_state: ResolvedLiquidState,
    temperature_axis: String,
    viscosity_basis: LiquidViscosityBasis,
}

impl<'a> LiquidTemperatureLaw<'a> {
    /// Resolve once and bind the selected four claims for all later temperatures.
    /// Typed temperature axes must be absolute kelvin, not a difference. Legacy
    /// untyped `T` is accepted; any other axis requires its explicit descriptor.
    ///
    /// # Errors
    /// Refuses a missing/mistyped temperature and all source-resolution errors.
    pub fn try_new(
        card: &'a MaterialCard, point: &QueryPoint, temperature_axis: &str,
        viscosity_basis: LiquidViscosityBasis, selection: MaterialPropertySelection,
    ) -> Result<Self, MaterialStatePointError> {
        validate_temperature(point, temperature_axis)?;
        let reference_state = viscosity_basis.resolve(card, point, selection)?;
        Ok(Self {
            card, reference_point: point.clone(), reference_state,
            temperature_axis: temperature_axis.to_owned(), viscosity_basis,
        })
    }

    /// Bind an existing complete source bundle to its original immutable card.
    ///
    /// This preserves optional coupled claims along with rho, Cp, conductivity
    /// and the explicitly selected viscosity basis. Every later temperature
    /// refresh pins the complete original bundle; no second material epoch is
    /// introduced for additional coefficients such as surface tension.
    ///
    /// # Errors
    /// Refuses a different card, invalid temperature coordinate or incompatible
    /// liquid requirements in the supplied bundle.
    pub fn try_from_material(
        card: &'a MaterialCard,
        material: ResolvedMaterialStatePoint,
        temperature_axis: &str,
        viscosity_basis: LiquidViscosityBasis,
    ) -> Result<Self, MaterialStatePointError> {
        if material.card_identity() != card.content_hash() {
            return Err(invalid("liquid temperature law requires the original material card"));
        }
        let reference_point = material.query();
        validate_temperature(&reference_point, temperature_axis)?;
        let reference_state = ResolvedLiquidState::try_from_material(material, viscosity_basis)?;
        Ok(Self {
            card, reference_point, reference_state,
            temperature_axis: temperature_axis.to_owned(), viscosity_basis,
        })
    }

    /// Evaluate all transport properties at a new temperature, atomically.
    /// Kinematic-source conversion uses the NEW density, not the reference one.
    /// Every nonthermal coordinate and its descriptor stays exactly unchanged.
    ///
    /// # Errors
    /// Refuses invalid temperature, out-of-support pins and inadmissible results.
    pub fn at_temperature(&self, temperature: Temperature)
        -> Result<ResolvedLiquidState, MaterialStatePointError>
    {
        if !temperature.value().is_finite() || temperature.value() <= 0.0 {
            return Err(invalid("liquid evaluation temperature must be positive and finite"));
        }
        let point = self.reference_point.clone().with(&self.temperature_axis, temperature.value())
            .map_err(|_| invalid("liquid temperature violates its retained quantity convention"))?;
        self.reference_state.at_query_point(self.card, &point)
    }

    /// The immutable card owning every evaluated state.
    #[must_use]
    pub const fn card(&self) -> &'a MaterialCard { self.card }
    /// Initial selected source state, before temperature updates.
    #[must_use]
    pub const fn reference_state(&self) -> &ResolvedLiquidState { &self.reference_state }
    /// Complete initial query, including all nonthermal conditions.
    #[must_use]
    pub const fn reference_point(&self) -> &QueryPoint { &self.reference_point }
    /// Explicit temperature coordinate used by both material and caloric owners.
    #[must_use]
    pub fn temperature_axis(&self) -> &str { &self.temperature_axis }
    /// Original source viscosity convention.
    #[must_use]
    pub const fn viscosity_basis(&self) -> LiquidViscosityBasis { self.viscosity_basis }
}

fn validate_temperature(
    point: &QueryPoint,
    temperature_axis: &str,
) -> Result<(), MaterialStatePointError> {
    let absolute = QuantitySpec::semantic(SemanticType::new(
        QuantityKind::AbsoluteTemperature, ValueForm::Static,
    ));
    let descriptor_ok = match point.axis_quantities().get(temperature_axis) {
        Some(q) => *q == absolute || *q == QuantitySpec::dimensional(Temperature::DIMS),
        None => temperature_axis == "T",
    };
    if !descriptor_ok || !point.axes().get(temperature_axis)
        .is_some_and(|t| t.is_finite() && *t > 0.0)
    {
        return Err(invalid("liquid temperature law needs an explicit positive absolute-temperature coordinate"));
    }
    Ok(())
}

fn invalid(quantity: &'static str) -> MaterialStatePointError {
    MaterialStatePointError::InvalidDerived { quantity }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fluid_transport::FluidTransport;
    use crate::liquid::{DYNAMIC_VISCOSITY_PROPERTY, KINEMATIC_VISCOSITY_PROPERTY};
    use crate::state_point::{
        DENSITY_PROPERTY, SPECIFIC_HEAT_CAPACITY_PROPERTY, THERMAL_CONDUCTIVITY_PROPERTY,
    };
    use fs_evidence::ValidityDomain;
    use fs_matdb::{
        ClaimSet, InterpolationPolicy, MaterialStateId, PropertyClaim, PropertyKey, PropertyValue,
        Provenance, UncertaintyModel,
    };
    use fs_qty::Dims;

    fn card(basis: LiquidViscosityBasis, split_viscosity: bool) -> MaterialCard {
        let mut claims = ClaimSet::new();
        for (name, dims, values) in [
            (DENSITY_PROPERTY, Dims([-3,1,0,0,0,0]), [1000.0, 800.0]),
            (SPECIFIC_HEAT_CAPACITY_PROPERTY, Dims([2,0,-2,-1,0,0]), [2000.0, 4000.0]),
            (THERMAL_CONDUCTIVITY_PROPERTY, Dims([1,1,-3,-1,0,0]), [0.5, 1.0]),
            (basis.property(), if basis == LiquidViscosityBasis::Dynamic {
                Dims([-1,1,-1,0,0,0]) } else { Dims([2,0,-1,0,0,0]) },
                if basis == LiquidViscosityBasis::Dynamic { [0.002, 0.001] } else { [2e-6, 1e-6] }),
        ] {
            let quantity = if name == THERMAL_CONDUCTIVITY_PROPERTY {
                QuantitySpec::semantic(SemanticType::new(QuantityKind::ThermalConductivity, ValueForm::Static))
            } else { QuantitySpec::dimensional(dims) };
            let mut claim = PropertyClaim {
                key: PropertyKey::with_quantity(name, quantity),
                value: PropertyValue::Curve { abscissa: "T".into(), abscissa_dims: Temperature::DIMS,
                    knots: vec![(300.0, values[0]), (400.0, values[1])], dims },
                validity: ValidityDomain::unconstrained()
                    .with_quantity(
                        "T",
                        QuantitySpec::semantic(SemanticType::new(
                            QuantityKind::AbsoluteTemperature,
                            ValueForm::Static,
                        )),
                        300.0,
                        400.0,
                    )
                    .with("pressure", 100000.0, 100000.0),
                uncertainty: UncertaintyModel::Unstated, interpolation: InterpolationPolicy::LinearInside,
                observations: vec![], provenance: Provenance {
                    source: "synthetic liquid law, not measurements".into(), license: "test-only".into(), artifact: None },
            };
            if split_viscosity && name == basis.property() {
                claim.validity = claim.validity.with("T", 300.0, 350.0);
                claims.insert_claim(claim.clone()).unwrap();
                claim.validity = claim.validity.with("T", 351.0, 400.0);
            }
            claims.insert_claim(claim).unwrap();
        }
        MaterialCard::assemble(MaterialStateId { chemistry: "synthetic liquid".into(), phase: "liquid".into(),
            process: "numerical fixture".into(), revision: 0 }, claims, vec![]).unwrap()
    }
    fn point() -> QueryPoint {
        QueryPoint::new().with_quantity("T", QuantitySpec::semantic(SemanticType::new(
            QuantityKind::AbsoluteTemperature, ValueForm::Static)), 300.0).unwrap()
            .with("pressure", 100000.0).unwrap()
    }

    #[test]
    fn temperature_updates_transport_and_preserves_claims_and_context() {
        let card = card(LiquidViscosityBasis::Dynamic, false);
        let law = LiquidTemperatureLaw::try_new(&card, &point(), "T", LiquidViscosityBasis::Dynamic,
            MaterialPropertySelection::SingleClaimOnly).unwrap();
        let hot = law.at_temperature(Temperature::new(400.0)).unwrap();
        assert_eq!(
            hot,
            law.reference_state()
                .at_query_point(
                    &card,
                    &law.reference_state()
                        .material()
                        .query()
                        .with("T", 400.0)
                        .unwrap()
                )
                .unwrap()
        );
        assert_eq!(hot.density_kg_m3(), 800.0);
        assert_eq!(hot.specific_heat_j_kg_k(), 4000.0);
        assert_eq!(hot.dynamic_viscosity_pa_s(), 0.001);
        assert_eq!(hot.thermal_conductivity_w_m_k(), 1.0);
        assert_eq!(hot.transport_state().unwrap().prandtl(), 4.0);
        for name in [DENSITY_PROPERTY, SPECIFIC_HEAT_CAPACITY_PROPERTY, THERMAL_CONDUCTIVITY_PROPERTY, DYNAMIC_VISCOSITY_PROPERTY] {
            let receipt = &hot.material().property(name).unwrap().answer().receipt;
            assert_eq!(receipt.selected, law.reference_state().material().property(name).unwrap().answer().receipt.selected);
            assert_eq!(receipt.axis_quantities, *point().axis_quantities());
            assert!(receipt.query_point.contains(&("pressure".into(), 100000.0)));
        }
        assert_eq!(hot, law.at_temperature(Temperature::new(400.0)).unwrap());
        assert_eq!(law.reference_point(), &point());
    }

    #[test]
    fn kinematic_basis_uses_the_updated_density() {
        let card = card(LiquidViscosityBasis::Kinematic, false);
        let law = LiquidTemperatureLaw::try_new(&card, &point(), "T", LiquidViscosityBasis::Kinematic,
            MaterialPropertySelection::SingleClaimOnly).unwrap();
        let hot = law.at_temperature(Temperature::new(400.0)).unwrap();
        assert!((hot.dynamic_viscosity_pa_s() - 0.0008).abs() < 1e-18);
        assert!(hot.material().property(KINEMATIC_VISCOSITY_PROPERTY).is_some());
        assert!(hot.material().property(DYNAMIC_VISCOSITY_PROPERTY).is_none());
    }

    #[test]
    fn updates_do_not_silently_switch_claims_at_a_source_boundary() {
        let card = card(LiquidViscosityBasis::Dynamic, true);
        let law = LiquidTemperatureLaw::try_new(&card, &point(), "T", LiquidViscosityBasis::Dynamic,
            MaterialPropertySelection::SingleClaimOnly).unwrap();
        assert!(resolve_liquid_state(&card, &point().with("T", 375.0).unwrap(),
            MaterialPropertySelection::SingleClaimOnly).is_ok());
        assert!(law.at_temperature(Temperature::new(375.0)).is_err());
        assert!(law.at_temperature(Temperature::new(325.0)).is_ok());
    }

    #[test]
    fn invalid_temperature_semantics_and_source_domain_refuse() {
        let card = card(LiquidViscosityBasis::Dynamic, false);
        let bad = point().with_quantity("T", QuantitySpec::semantic(SemanticType::new(
            QuantityKind::TemperatureDifference, ValueForm::Static)), 300.0).unwrap();
        assert!(LiquidTemperatureLaw::try_new(&card, &bad, "T", LiquidViscosityBasis::Dynamic,
            MaterialPropertySelection::SingleClaimOnly).is_err());
        let law = LiquidTemperatureLaw::try_new(&card, &point(), "T", LiquidViscosityBasis::Dynamic,
            MaterialPropertySelection::SingleClaimOnly).unwrap();
        for t in [0.0, -1.0, f64::NAN, f64::INFINITY, 299.0, 401.0] {
            assert!(law.at_temperature(Temperature::new(t)).is_err());
        }
    }
}
