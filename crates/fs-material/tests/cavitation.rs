//! G0/G3 shared cavitation coefficients. Manufactured cards below are numerical
//! fixtures; the original IF97 comparison checks source preservation.

use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::cavitation::{
    CavitationProperties, CavitationPropertyRequirements, CavitationState, SaturatedVaporState,
    SaturatedWaterCavitation, resolve_cavitation_state,
};
use fs_material::liquid::surface::SURFACE_TENSION_DIMS;
use fs_material::liquid::{DYNAMIC_VISCOSITY_DIMS, KINEMATIC_VISCOSITY_DIMS};
use fs_material::liquid_temperature::LiquidViscosityBasis;
use fs_material::state_point::{
    MaterialPropertySelection, ScalarAdmissibility, ScalarPropertyRequirement,
};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Density, Dims, Pressure, QuantitySpec, Temperature, Velocity};
use fs_thermochem::water::if97::If97SaturationPair;

fn quantity(kind: QuantityKind) -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(kind, ValueForm::Static))
}
fn point(t: f64) -> QueryPoint {
    QueryPoint::new()
        .with_quantity(
            "temperature",
            quantity(QuantityKind::AbsoluteTemperature),
            t,
        )
        .unwrap()
        .with_quantity("pressure", quantity(QuantityKind::Pressure), 100_000.0)
        .unwrap()
}
fn requirement(name: &str, dims: Dims) -> ScalarPropertyRequirement {
    ScalarPropertyRequirement::try_new(name, dims, ScalarAdmissibility::StrictlyPositive).unwrap()
}
fn fixture(
    basis: LiquidViscosityBasis,
    with_vapor: bool,
) -> (MaterialCard, CavitationPropertyRequirements) {
    let viscosity_dims = match basis {
        LiquidViscosityBasis::Dynamic => DYNAMIC_VISCOSITY_DIMS,
        LiquidViscosityBasis::Kinematic => KINEMATIC_VISCOSITY_DIMS,
    };
    let viscosity = match basis {
        LiquidViscosityBasis::Dynamic => [0.001, 0.0005],
        LiquidViscosityBasis::Kinematic => [1e-6, 0.5e-6],
    };
    let mut rows = vec![
        ("source-rho", Density::DIMS, [1000.0, 900.0]),
        ("source-viscosity", viscosity_dims, viscosity),
        ("source-sigma", SURFACE_TENSION_DIMS, [0.07, 0.06]),
        (
            "source-equilibrium-pressure",
            Pressure::DIMS,
            [2000.0, 20_000.0],
        ),
        ("source-c", Velocity::DIMS, [1500.0, 1600.0]),
    ];
    if with_vapor {
        rows.push(("source-vapor-density", Density::DIMS, [0.02, 0.2]));
        rows.push(("source-latent-h", Dims([2, 0, -2, 0, 0, 0]), [2.5e6, 2.2e6]));
    }
    let requirements = CavitationPropertyRequirements::try_new(
        std::array::from_fn(|i| requirement(rows[i].0, rows[i].1)),
        with_vapor.then(|| std::array::from_fn(|i| requirement(rows[i + 5].0, rows[i + 5].1))),
        basis,
    )
    .unwrap();
    let mut claims = ClaimSet::new();
    for (name, dims, values) in rows {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Curve {
                    abscissa: "temperature".into(),
                    abscissa_dims: Temperature::DIMS,
                    knots: vec![(300.0, values[0]), (400.0, values[1])],
                    dims,
                },
                validity: ValidityDomain::unconstrained()
                    .with_quantity(
                        "temperature",
                        quantity(QuantityKind::AbsoluteTemperature),
                        300.0,
                        400.0,
                    )
                    .with_quantity(
                        "pressure",
                        quantity(QuantityKind::Pressure),
                        100_000.0,
                        100_000.0,
                    ),
                interpolation: InterpolationPolicy::LinearInside,
                uncertainty: UncertaintyModel::Unstated,
                observations: vec![],
                provenance: Provenance {
                    source: "manufactured cavitation coefficients".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "manufactured fluid".into(),
            phase: "liquid with equilibrium vapor data".into(),
            process: "numerical fixture".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    (card, requirements)
}

#[test]
fn numerical_limits_are_explicit_and_missing_vapor_is_not_invented() {
    let state = CavitationState::try_new(300.0, 1000.0, 0.0, 0.0, 0.0, 1500.0).unwrap();
    assert_eq!(state.saturated_vapor(), None);
    assert_eq!(
        (&state as &dyn CavitationProperties)
            .cavitation_state()
            .unwrap(),
        state
    );
    for bad in [f64::NAN, f64::INFINITY, -1.0] {
        assert!(CavitationState::try_new(300.0, 1000.0, bad, 0.0, 0.0, 1500.0).is_err());
        assert!(CavitationState::try_new(300.0, 1000.0, 0.0, bad, 0.0, 1500.0).is_err());
        assert!(CavitationState::try_new(300.0, 1000.0, 0.0, 0.0, bad, 1500.0).is_err());
    }
    assert!(CavitationState::try_new(0.0, 1000.0, 0.0, 0.0, 0.0, 1500.0).is_err());
    assert!(CavitationState::try_new(300.0, 0.0, 0.0, 0.0, 0.0, 1500.0).is_err());
    assert!(CavitationState::try_new(300.0, 1000.0, 0.0, 0.0, 0.0, 0.0).is_err());
    assert!(
        state
            .with_saturated_vapor(SaturatedVaporState::try_new(1000.0, 2.4e6).unwrap())
            .is_err()
    );
    assert!(SaturatedVaporState::try_new(0.02, 0.0).is_err());
}

#[test]
fn original_water_saturation_and_transport_are_retained_without_another_fit() {
    let mut previous_pressure = 0.0;
    for t in [293.15, 333.15, 373.15] {
        let source = If97SaturationPair::at_temperature(Temperature::new(t)).unwrap();
        let water = SaturatedWaterCavitation::try_new(source, 0.07).unwrap();
        let state = water.cavitation_state().unwrap();
        assert_eq!(*water.saturation(), source);
        assert_eq!(state.temperature_k(), t);
        assert_eq!(state.saturation_pressure_pa(), source.pressure().value());
        assert_eq!(
            state.liquid_density_kg_m3(),
            source.liquid().density().value()
        );
        assert_eq!(
            state.liquid_sound_speed_m_s(),
            source.liquid().sound_speed().value()
        );
        assert_eq!(
            state.dynamic_viscosity_pa_s(),
            water.liquid_transport().dynamic_viscosity().value()
        );
        assert_eq!(
            state.saturated_vapor().unwrap().density_kg_m3(),
            source.vapor().density().value()
        );
        assert_eq!(
            state.saturated_vapor().unwrap().latent_enthalpy_j_kg(),
            source.latent_enthalpy().value()
        );
        assert!(state.saturation_pressure_pa() > previous_pressure);
        previous_pressure = state.saturation_pressure_pa();
    }
}

#[test]
fn all_source_roles_refresh_together_with_original_keys_and_viscosity_basis() {
    for basis in [
        LiquidViscosityBasis::Dynamic,
        LiquidViscosityBasis::Kinematic,
    ] {
        let (card, requirements) = fixture(basis, true);
        let initial = resolve_cavitation_state(
            &card,
            &point(300.0),
            &requirements,
            "temperature",
            MaterialPropertySelection::SingleClaimOnly,
        )
        .unwrap();
        let updated = initial.at_query_point(&card, &point(400.0)).unwrap();
        assert_eq!(
            updated.material().card_identity(),
            initial.material().card_identity()
        );
        assert_eq!(updated.material().properties().len(), 7);
        assert_eq!(updated.requirements(), &requirements);
        assert_eq!(updated.state().liquid_density_kg_m3(), 900.0);
        let expected_mu = if basis == LiquidViscosityBasis::Kinematic {
            900.0 * 0.5e-6
        } else {
            0.0005
        };
        assert_eq!(updated.state().dynamic_viscosity_pa_s(), expected_mu);
        assert_eq!(updated.state().saturation_pressure_pa(), 20_000.0);
        assert_eq!(updated.state().surface_tension_n_m(), 0.06);
        assert_eq!(updated.state().liquid_sound_speed_m_s(), 1600.0);
        assert_eq!(
            updated
                .state()
                .saturated_vapor()
                .unwrap()
                .latent_enthalpy_j_kg(),
            2.2e6
        );
        assert_eq!(initial.state().liquid_density_kg_m3(), 1000.0);
        assert!(initial.at_query_point(&card, &point(401.0)).is_err());
        assert!(
            updated
                .material()
                .property("specific_heat_capacity")
                .is_none()
        );
        assert!(
            updated
                .material()
                .property("thermal_conductivity")
                .is_none()
        );
    }
}

#[test]
fn optional_calorics_and_temperature_roles_are_not_guessed() {
    let (card, requirements) = fixture(LiquidViscosityBasis::Dynamic, false);
    let state = resolve_cavitation_state(
        &card,
        &point(350.0),
        &requirements,
        "temperature",
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    assert_eq!(state.state().saturated_vapor(), None);
    assert_eq!(state.material().properties().len(), 5);
    assert!(
        resolve_cavitation_state(
            &card,
            &point(350.0),
            &requirements,
            "pressure",
            MaterialPropertySelection::SingleClaimOnly,
        )
        .is_err()
    );
    let wrong = requirement("pressure-as-sound", Pressure::DIMS);
    let mut roles = requirements.liquid().clone();
    roles[4] = wrong;
    assert!(
        CavitationPropertyRequirements::try_new(roles, None, LiquidViscosityBasis::Dynamic)
            .is_err()
    );
    let mut roles = requirements.liquid().clone();
    roles[4] = requirement(roles[0].name(), Velocity::DIMS);
    assert!(
        CavitationPropertyRequirements::try_new(roles, None, LiquidViscosityBasis::Dynamic)
            .is_err()
    );
}
