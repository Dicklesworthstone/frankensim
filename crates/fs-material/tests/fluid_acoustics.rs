//! G0/G3: a common acoustic state preserves each original source and query.
//! Manufactured source coefficients below are numerical fixtures, not measurements.

use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::fluid_acoustics::{
    FluidAcousticParameter, FluidAcousticState, FluidAcoustics, resolve_acoustic_liquid_state,
    resolve_acoustic_liquid_state_with_requirements, resolve_fluid_acoustic_state,
};
use fs_material::fluid_transport::FluidTransport;
use fs_material::gas::{GasSpec, GasState};
use fs_material::liquid::{
    DYNAMIC_VISCOSITY_DIMS, KINEMATIC_VISCOSITY_DIMS, LiquidPropertyRequirements,
    resolve_liquid_state_with_requirements,
};
use fs_material::liquid_temperature::LiquidViscosityBasis;
use fs_material::state_point::{
    MaterialPropertySelection, SPECIFIC_HEAT_CAPACITY_DIMS, ScalarAdmissibility,
    ScalarPropertyRequirement,
};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Density, Dims, Pressure, QuantitySpec, Temperature, Velocity};
use fs_thermochem::water::if97::region3::If97DenseState;
use fs_thermochem::water::if97::regions::If97FluidState;
use fs_thermochem::water::if97::{If97Region, If97State};

fn near(actual: f64, expected: f64) {
    assert!(
        (actual / expected - 1.0).abs() < 3e-13,
        "{actual} != {expected}"
    );
}

fn temperature_quantity() -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(
        QuantityKind::AbsoluteTemperature,
        ValueForm::Static,
    ))
}

fn conductivity_quantity() -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(
        QuantityKind::ThermalConductivity,
        ValueForm::Static,
    ))
}

fn point(t: f64) -> QueryPoint {
    QueryPoint::new()
        .with_quantity("T", temperature_quantity(), t)
        .unwrap()
        .with("pressure", 100_000.0)
        .unwrap()
        .with("phase-liquid", 1.0)
        .unwrap()
}

fn requirement(name: &str, dims: Dims) -> ScalarPropertyRequirement {
    ScalarPropertyRequirement::try_new(name, dims, ScalarAdmissibility::StrictlyPositive).unwrap()
}

fn source_card(data: Vec<(&str, QuantitySpec, [f64; 2])>) -> MaterialCard {
    let mut claims = ClaimSet::new();
    for (name, quantity, values) in data {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::with_quantity(name, quantity),
                value: PropertyValue::Curve {
                    abscissa: "T".into(),
                    abscissa_dims: Temperature::DIMS,
                    knots: vec![(300.0, values[0]), (400.0, values[1])],
                    dims: quantity.dims(),
                },
                validity: ValidityDomain::unconstrained()
                    .with_quantity("T", temperature_quantity(), 300.0, 400.0)
                    .with("pressure", 100_000.0, 100_000.0)
                    .with("phase-liquid", 1.0, 1.0),
                interpolation: InterpolationPolicy::LinearInside,
                uncertainty: UncertaintyModel::Unstated,
                observations: vec![],
                provenance: Provenance {
                    source: "manufactured fluid acoustic fixture".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "manufactured fluid".into(),
            phase: "liquid".into(),
            process: "numerical fixture".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}

fn liquid_fixture(
    basis: LiquidViscosityBasis,
    custom_names: bool,
) -> (MaterialCard, LiquidPropertyRequirements) {
    let names = if custom_names {
        [
            "mass-density",
            "specific-heat-capacity",
            "thermal-conductivity",
            "source-viscosity",
        ]
    } else {
        [
            "density",
            "specific_heat_capacity",
            "thermal_conductivity",
            if basis == LiquidViscosityBasis::Dynamic {
                "dynamic_viscosity"
            } else {
                "kinematic_viscosity"
            },
        ]
    };
    let (viscosity_dims, viscosity) = match basis {
        LiquidViscosityBasis::Dynamic => (DYNAMIC_VISCOSITY_DIMS, [0.001, 0.0005]),
        LiquidViscosityBasis::Kinematic => (KINEMATIC_VISCOSITY_DIMS, [1e-6, 0.5e-6]),
    };
    let requirements = LiquidPropertyRequirements::try_new(
        requirement(names[0], Density::DIMS),
        requirement(names[1], SPECIFIC_HEAT_CAPACITY_DIMS),
        ScalarPropertyRequirement::try_with_quantity(
            names[2],
            conductivity_quantity(),
            ScalarAdmissibility::StrictlyPositive,
        )
        .unwrap(),
        requirement(names[3], viscosity_dims),
        basis,
    )
    .unwrap();
    let card = source_card(vec![
        (
            names[0],
            QuantitySpec::dimensional(Density::DIMS),
            [1000.0, 960.0],
        ),
        (
            names[1],
            QuantitySpec::dimensional(SPECIFIC_HEAT_CAPACITY_DIMS),
            [4000.0, 4200.0],
        ),
        (names[2], conductivity_quantity(), [0.6, 0.65]),
        (
            names[3],
            QuantitySpec::dimensional(viscosity_dims),
            viscosity,
        ),
        (
            "sound_speed",
            QuantitySpec::dimensional(Velocity::DIMS),
            [1500.0, 1600.0],
        ),
    ]);
    (card, requirements)
}

#[test]
fn gas_pressure_scaling_and_coherence_are_shared_with_original_owner() {
    let a = GasState::try_new(&GasSpec::dry_air_ussa1976(), 300.0, 100_000.0).unwrap();
    let b = GasState::try_new(&GasSpec::dry_air_ussa1976(), 300.0, 200_000.0).unwrap();
    let av = a.acoustic_state().unwrap();
    let bv = b.acoustic_state().unwrap();
    assert_eq!(av.sound_speed_m_s(), bv.sound_speed_m_s());
    near(
        bv.adiabatic_bulk_modulus_pa(),
        2.0 * av.adiabatic_bulk_modulus_pa(),
    );
    near(
        bv.characteristic_impedance_pa_s_m(),
        2.0 * av.characteristic_impedance_pa_s_m(),
    );
    near(
        bv.adiabatic_compressibility_per_pa(),
        av.adiabatic_compressibility_per_pa() / 2.0,
    );
    let borrowed: &dyn FluidAcoustics = &a;
    assert_eq!((&borrowed).acoustic_state().unwrap(), av);
    let mut bad = a.clone();
    bad.sound_speed *= 1.01;
    assert!(bad.acoustic_state().is_err());
    let mut bad = a;
    bad.prandtl *= 1.01;
    assert!(bad.acoustic_state().is_err());
}

#[test]
fn if97_liquid_steam_and_dense_states_use_original_thermodynamic_derivatives() {
    for (p, t, region) in [
        (3e6, 300.0, If97Region::Region1),
        (3500.0, 700.0, If97Region::Region2),
    ] {
        let water =
            If97State::at_pressure_temperature(Pressure::new(p), Temperature::new(t), region)
                .unwrap();
        let acoustic = water.acoustic_state().unwrap();
        assert_eq!(
            acoustic.density_kg_m3().to_bits(),
            water.density().value().to_bits()
        );
        assert_eq!(
            acoustic.sound_speed_m_s().to_bits(),
            water.sound_speed().value().to_bits()
        );
        near(
            acoustic.adiabatic_compressibility_per_pa(),
            water.isothermal_compressibility().value() * water.cv().value() / water.cp().value(),
        );
        assert_eq!(
            If97FluidState::Gibbs(water).acoustic_state().unwrap(),
            acoustic
        );
    }
    let dense =
        If97DenseState::at_density_temperature(Density::new(500.0), Temperature::new(650.0))
            .unwrap();
    let acoustic = dense.acoustic_state().unwrap();
    assert_eq!(
        acoustic.sound_speed_m_s().to_bits(),
        dense.sound_speed().value().to_bits()
    );
    near(
        acoustic.adiabatic_compressibility_per_pa(),
        dense.isothermal_compressibility().value() * dense.cv().value() / dense.cp().value(),
    );
}

#[test]
fn source_acoustics_keep_complete_query_and_original_claims_without_transport_defaults() {
    let card = source_card(vec![
        (
            "rho-original",
            QuantitySpec::dimensional(Density::DIMS),
            [1000.0, 960.0],
        ),
        (
            "measured-wave-speed",
            QuantitySpec::dimensional(Velocity::DIMS),
            [1500.0, 1600.0],
        ),
    ]);
    let density = requirement("rho-original", Density::DIMS);
    let parameter =
        FluidAcousticParameter::SoundSpeed(requirement("measured-wave-speed", Velocity::DIMS));
    let cold = resolve_fluid_acoustic_state(
        &card,
        &point(300.0),
        &density,
        &parameter,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let hot = cold.at_query_point(&card, &point(400.0)).unwrap();
    assert_eq!(hot.material().card_identity(), card.content_hash());
    assert_eq!(hot.material().query(), point(400.0));
    assert_eq!(hot.material().properties().len(), 2);
    assert_ne!(cold.identity(), hot.identity());
    for property in hot.material().properties() {
        let previous = cold
            .material()
            .property(property.requirement().name())
            .unwrap();
        assert_eq!(
            property.answer().receipt.selected,
            previous.answer().receipt.selected
        );
        assert_eq!(
            property.answer().receipt.query_point,
            hot.material().query_point()
        );
    }
    assert_eq!(hot.state().sound_speed_m_s(), 1600.0);
    assert!(cold.at_query_point(&card, &point(401.0)).is_err());
    assert!(
        cold.at_query_point(&card, &point(300.0).with("pressure", 200_000.0).unwrap())
            .is_err()
    );
}

#[test]
fn explicit_adiabatic_modulus_role_preserves_source_and_rejects_wrong_quantity() {
    let card = source_card(vec![
        (
            "density",
            QuantitySpec::dimensional(Density::DIMS),
            [1000.0; 2],
        ),
        (
            "bulk-law",
            QuantitySpec::dimensional(Pressure::DIMS),
            [2.25e9; 2],
        ),
    ]);
    let density = requirement("density", Density::DIMS);
    let modulus = requirement("bulk-law", Pressure::DIMS);
    let acoustic = resolve_fluid_acoustic_state(
        &card,
        &point(300.0),
        &density,
        &FluidAcousticParameter::AdiabaticBulkModulus(modulus.clone()),
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    assert_eq!(acoustic.state().adiabatic_bulk_modulus_pa(), 2.25e9);
    assert_eq!(acoustic.state().sound_speed_m_s(), 1500.0);
    assert!(
        resolve_fluid_acoustic_state(
            &card,
            &point(300.0),
            &density,
            &FluidAcousticParameter::SoundSpeed(modulus),
            MaterialPropertySelection::SingleClaimOnly
        )
        .is_err()
    );
    let pressure = ScalarPropertyRequirement::try_with_quantity(
        "bulk-law",
        QuantitySpec::semantic(SemanticType::new(QuantityKind::Pressure, ValueForm::Static)),
        ScalarAdmissibility::StrictlyPositive,
    )
    .unwrap();
    assert!(
        resolve_fluid_acoustic_state(
            &card,
            &point(300.0),
            &density,
            &FluidAcousticParameter::AdiabaticBulkModulus(pressure),
            MaterialPropertySelection::SingleClaimOnly
        )
        .is_err()
    );
}

#[test]
fn liquid_transport_and_acoustics_refresh_atomically_for_both_viscosity_bases() {
    for basis in [
        LiquidViscosityBasis::Dynamic,
        LiquidViscosityBasis::Kinematic,
    ] {
        let (card, _) = liquid_fixture(basis, false);
        let parameter = FluidAcousticParameter::sound_speed().unwrap();
        let cold = resolve_acoustic_liquid_state(
            &card,
            &point(300.0),
            &parameter,
            basis,
            MaterialPropertySelection::SingleClaimOnly,
        )
        .unwrap();
        let hot = cold.at_query_point(&card, &point(400.0)).unwrap();
        let acoustic = hot.acoustic_state().unwrap();
        let transport = hot.transport_state().unwrap();
        assert_eq!(
            acoustic.density_kg_m3().to_bits(),
            transport.density_kg_m3().to_bits()
        );
        assert_eq!(acoustic.density_kg_m3(), 960.0);
        assert_eq!(acoustic.sound_speed_m_s(), 1600.0);
        assert_eq!(hot.material().properties().len(), 5);
        assert_eq!(hot.material().query(), point(400.0));
        near(
            transport.dynamic_viscosity_pa_s(),
            if basis == LiquidViscosityBasis::Dynamic {
                0.0005
            } else {
                960.0 * 0.5e-6
            },
        );
        assert!(cold.at_query_point(&card, &point(401.0)).is_err());
        assert!(
            resolve_acoustic_liquid_state(
                &card,
                &point(300.0),
                &FluidAcousticParameter::adiabatic_bulk_modulus().unwrap(),
                basis,
                MaterialPropertySelection::SingleClaimOnly
            )
            .is_err()
        );
    }
}

#[test]
fn exact_source_vocabulary_requires_no_aliases_or_claim_rewriting() {
    let (card, requirements) = liquid_fixture(LiquidViscosityBasis::Dynamic, true);
    let original_hash = card.content_hash();
    let liquid = resolve_liquid_state_with_requirements(
        &card,
        &point(300.0),
        &requirements,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let combined = resolve_acoustic_liquid_state_with_requirements(
        &card,
        &point(300.0),
        &FluidAcousticParameter::sound_speed().unwrap(),
        &requirements,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    assert_eq!(
        liquid.transport_state().unwrap(),
        combined.transport_state().unwrap()
    );
    assert_eq!(combined.material().card_identity(), original_hash);
    assert_eq!(card.content_hash(), original_hash);
    assert_eq!(combined.liquid_requirements(), &requirements);
    assert!(
        combined
            .material()
            .property("specific-heat-capacity")
            .is_some()
    );
    assert!(
        combined
            .material()
            .property("specific_heat_capacity")
            .is_none()
    );
    let expanding = liquid
        .with_density_expansion(
            &card,
            "T",
            Temperature::new(300.0),
            Temperature::new(400.0),
        )
        .unwrap();
    let hot = expanding.at_query_point(&card, &point(400.0)).unwrap();
    near(hot.density_kg_m3(), 960.0);
    near(hot.density_expansion().unwrap().density_slope_kg_m3_k(), -0.4);
    near(hot.volumetric_expansion_per_k().unwrap(), 0.4 / 960.0);
    assert_eq!(hot.material().card_identity(), original_hash);
    assert!(hot.material().property("mass-density").is_some());
    assert!(hot.material().property("density").is_none());
    let p = requirements.properties();
    assert!(
        LiquidPropertyRequirements::try_new(
            p[0].clone(),
            requirement("bad-cp", Pressure::DIMS),
            p[2].clone(),
            p[3].clone(),
            LiquidViscosityBasis::Dynamic
        )
        .is_err()
    );
}

#[test]
fn invalid_or_unrepresentable_acoustic_coefficients_refuse() {
    for v in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(FluidAcousticState::try_new(v, 1500.0).is_err());
        assert!(FluidAcousticState::try_new(1000.0, v).is_err());
        assert!(FluidAcousticState::from_adiabatic_bulk_modulus(v, 2e9).is_err());
        assert!(FluidAcousticState::from_adiabatic_bulk_modulus(1000.0, v).is_err());
    }
    assert!(FluidAcousticState::try_new(f64::MAX, 2.0).is_err());
    assert!(FluidAcousticState::try_new(f64::MIN_POSITIVE, f64::MIN_POSITIVE).is_err());
    assert!(FluidAcousticState::from_adiabatic_bulk_modulus(f64::MIN_POSITIVE, f64::MAX).is_err());
}
