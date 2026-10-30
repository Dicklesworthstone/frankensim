//! Source-preserving thermal expansion on the original liquid transport owner.

use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::fluid_transport::{FluidTransport, ThermalFluidTransport};
use fs_material::liquid::{
    DYNAMIC_VISCOSITY_DIMS, DYNAMIC_VISCOSITY_PROPERTY, KINEMATIC_VISCOSITY_DIMS,
    KINEMATIC_VISCOSITY_PROPERTY, ResolvedLiquidState, VOLUMETRIC_EXPANSION_DIMS,
    VOLUMETRIC_EXPANSION_PROPERTY, resolve_liquid_state, resolve_liquid_state_with_expansion,
};
use fs_material::liquid_temperature::{LiquidTemperatureLaw, LiquidViscosityBasis};
use fs_material::state_point::{
    DENSITY_PROPERTY, MaterialPropertySelection, SPECIFIC_HEAT_CAPACITY_DIMS,
    SPECIFIC_HEAT_CAPACITY_PROPERTY, ScalarAdmissibility, ScalarPropertyRequirement,
    THERMAL_CONDUCTIVITY_DIMS, THERMAL_CONDUCTIVITY_PROPERTY, resolve_material_state_point,
};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Density, Dims, QuantitySpec, Temperature};

fn absolute_temperature() -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(
        QuantityKind::AbsoluteTemperature,
        ValueForm::Static,
    ))
}

fn point(temperature: f64) -> QueryPoint {
    QueryPoint::new()
        .with_quantity("T", absolute_temperature(), temperature)
        .unwrap()
        .with("pressure", 100_000.0)
        .unwrap()
        .with("solute_mass_fraction", 0.2)
        .unwrap()
}

fn card(
    basis: LiquidViscosityBasis,
    density: [f64; 2],
    expansion: Option<[f64; 2]>,
    expansion_dims: Dims,
    expansion_maximum_temperature: f64,
) -> MaterialCard {
    let (viscosity_name, viscosity_dims, viscosity) = match basis {
        LiquidViscosityBasis::Dynamic => (
            DYNAMIC_VISCOSITY_PROPERTY,
            DYNAMIC_VISCOSITY_DIMS,
            [0.001, 0.0006],
        ),
        LiquidViscosityBasis::Kinematic => (
            KINEMATIC_VISCOSITY_PROPERTY,
            KINEMATIC_VISCOSITY_DIMS,
            [1.0e-6, 0.6e-6],
        ),
    };
    let mut data = vec![
        (DENSITY_PROPERTY, Density::DIMS, density),
        (
            SPECIFIC_HEAT_CAPACITY_PROPERTY,
            SPECIFIC_HEAT_CAPACITY_DIMS,
            [4000.0, 4200.0],
        ),
        (
            THERMAL_CONDUCTIVITY_PROPERTY,
            THERMAL_CONDUCTIVITY_DIMS,
            [0.6, 0.65],
        ),
        (viscosity_name, viscosity_dims, viscosity),
    ];
    if let Some(values) = expansion {
        data.push((VOLUMETRIC_EXPANSION_PROPERTY, expansion_dims, values));
    }
    let mut claims = ClaimSet::new();
    for (name, dims, values) in data {
        let quantity = if name == THERMAL_CONDUCTIVITY_PROPERTY {
            QuantitySpec::semantic(SemanticType::new(
                QuantityKind::ThermalConductivity,
                ValueForm::Static,
            ))
        } else {
            QuantitySpec::dimensional(dims)
        };
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::with_quantity(name, quantity),
                value: PropertyValue::Curve {
                    abscissa: "T".into(),
                    abscissa_dims: Temperature::DIMS,
                    knots: vec![(300.0, values[0]), (400.0, values[1])],
                    dims,
                },
                validity: ValidityDomain::unconstrained()
                    .with_quantity(
                        "T",
                        absolute_temperature(),
                        300.0,
                        if name == VOLUMETRIC_EXPANSION_PROPERTY {
                            expansion_maximum_temperature
                        } else {
                            400.0
                        },
                    )
                    .with("pressure", 100_000.0, 100_000.0)
                    .with("solute_mass_fraction", 0.2, 0.2),
                interpolation: InterpolationPolicy::LinearInside,
                uncertainty: UncertaintyModel::Unstated,
                observations: vec![],
                provenance: Provenance {
                    source: "synthetic liquid expansion fixture, not measurements".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic liquid".into(),
            phase: "liquid".into(),
            process: "numerical fixture".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}

#[test]
fn expansion_temperature_updates_share_the_original_transport_and_sources() {
    for basis in [
        LiquidViscosityBasis::Dynamic,
        LiquidViscosityBasis::Kinematic,
    ] {
        // Endpoint beta values agree with -(d rho/dT)/rho of this density
        // curve. Independent interpolation between endpoints remains the
        // declared source law, not a thermodynamic-consistency certification.
        let card = card(
            basis,
            [1000.0, 980.0],
            Some([0.2 / 1000.0, 0.2 / 980.0]),
            VOLUMETRIC_EXPANSION_DIMS,
            400.0,
        );
        let law = LiquidTemperatureLaw::try_new_with_expansion(
            &card,
            &point(300.0),
            "T",
            basis,
            MaterialPropertySelection::SingleClaimOnly,
        )
        .unwrap();
        let hot = law.at_temperature(Temperature::new(400.0)).unwrap();
        let thermal = hot.thermal_transport_state().unwrap();
        assert_eq!(hot.density_kg_m3(), 980.0);
        assert_eq!(hot.specific_heat_j_kg_k(), 4200.0);
        assert_eq!(hot.thermal_conductivity_w_m_k(), 0.65);
        assert_eq!(
            hot.dynamic_viscosity_pa_s(),
            match basis {
                LiquidViscosityBasis::Dynamic => 0.0006,
                LiquidViscosityBasis::Kinematic => 980.0 * 0.6e-6,
            }
        );
        assert_eq!(thermal.transport(), hot.transport_state().unwrap());
        assert_eq!(thermal.temperature_k(), 400.0);
        assert_eq!(thermal.volumetric_expansion_per_k(), Some(0.2 / 980.0));
        assert_ne!(
            hot.volumetric_expansion_per_k(),
            law.reference_state().volumetric_expansion_per_k()
        );
        assert_eq!(hot.material().card_identity(), card.content_hash());
        assert_eq!(hot.material().properties().len(), 5);
        for property in hot.material().properties() {
            let receipt = &property.answer().receipt;
            assert_eq!(
                receipt.selected,
                law.reference_state()
                    .material()
                    .property(property.requirement().name())
                    .unwrap()
                    .answer()
                    .receipt
                    .selected
            );
            assert_eq!(receipt.query_point, hot.material().query_point());
            assert_eq!(receipt.axis_quantities, *point(400.0).axis_quantities());
        }
        assert_eq!(hot.material().query(), point(400.0));
        assert_eq!(
            ResolvedLiquidState::try_from_material(hot.material().clone(), basis).unwrap(),
            hot
        );
        let rebound =
            LiquidTemperatureLaw::try_from_material(&card, hot.material().clone(), "T", basis)
                .unwrap();
        assert_eq!(
            &rebound.at_temperature(Temperature::new(300.0)).unwrap(),
            law.reference_state()
        );
    }
}

#[test]
fn negative_zero_and_unknown_expansion_remain_distinct() {
    for (density, beta) in [([1000.0, 1020.0], -0.2 / 1000.0), ([1000.0; 2], 0.0)] {
        let card = card(
            LiquidViscosityBasis::Dynamic,
            density,
            Some([beta, beta * density[0] / density[1]]),
            VOLUMETRIC_EXPANSION_DIMS,
            400.0,
        );
        let state = resolve_liquid_state_with_expansion(
            &card,
            &point(300.0),
            LiquidViscosityBasis::Dynamic,
            MaterialPropertySelection::SingleClaimOnly,
        )
        .unwrap();
        assert_eq!(
            state
                .thermal_transport_state()
                .unwrap()
                .volumetric_expansion_per_k(),
            Some(beta)
        );
        // Selecting only four properties is still an explicit missing-beta
        // path, even when a separate source claim exists on the card.
        let four = resolve_liquid_state(
            &card,
            &point(300.0),
            MaterialPropertySelection::SingleClaimOnly,
        )
        .unwrap();
        assert_eq!(
            four.thermal_transport_state()
                .unwrap()
                .volumetric_expansion_per_k(),
            None
        );
    }
    let missing = card(
        LiquidViscosityBasis::Dynamic,
        [1000.0, 980.0],
        None,
        VOLUMETRIC_EXPANSION_DIMS,
        400.0,
    );
    assert!(
        resolve_liquid_state_with_expansion(
            &missing,
            &point(300.0),
            LiquidViscosityBasis::Dynamic,
            MaterialPropertySelection::SingleClaimOnly,
        )
        .is_err()
    );
}

#[test]
fn expansion_support_limits_the_whole_refresh_atomically() {
    let card = card(
        LiquidViscosityBasis::Dynamic,
        [1000.0, 980.0],
        Some([0.2 / 1000.0, 0.2 / 980.0]),
        VOLUMETRIC_EXPANSION_DIMS,
        350.0,
    );
    let state = resolve_liquid_state_with_expansion(
        &card,
        &point(300.0),
        LiquidViscosityBasis::Dynamic,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let unchanged = state.clone();
    assert!(state.at_query_point(&card, &point(400.0)).is_err());
    assert_eq!(state, unchanged);
    assert!(
        resolve_liquid_state(
            &card,
            &point(400.0),
            MaterialPropertySelection::SingleClaimOnly,
        )
        .is_ok()
    );
}

#[test]
fn adopting_a_complete_bundle_does_not_erase_mistyped_expansion() {
    let wrong_dims = Temperature::DIMS;
    let card = card(
        LiquidViscosityBasis::Dynamic,
        [1000.0, 980.0],
        Some([1.0, 2.0]),
        wrong_dims,
        400.0,
    );
    assert!(
        resolve_liquid_state_with_expansion(
            &card,
            &point(300.0),
            LiquidViscosityBasis::Dynamic,
            MaterialPropertySelection::SingleClaimOnly,
        )
        .is_err()
    );
    let base = resolve_liquid_state(
        &card,
        &point(300.0),
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let mut requirements: Vec<_> = base
        .material()
        .properties()
        .iter()
        .map(|property| property.requirement().clone())
        .collect();
    requirements.push(
        ScalarPropertyRequirement::try_new(
            VOLUMETRIC_EXPANSION_PROPERTY,
            wrong_dims,
            ScalarAdmissibility::Finite,
        )
        .unwrap(),
    );
    let malformed = resolve_material_state_point(
        &card,
        &point(300.0),
        &requirements,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    assert!(
        ResolvedLiquidState::try_from_material(malformed, LiquidViscosityBasis::Dynamic).is_err()
    );
}
