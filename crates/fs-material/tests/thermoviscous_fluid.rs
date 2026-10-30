//! Coherent caloric/transport binding; numerical fixtures are not measurements.

use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::fluid_acoustics::thermoviscous::{
    ResolvedThermoviscousLiquidState, ThermoviscousFluid, ThermoviscousFluidState,
    resolve_thermoviscous_liquid_state_with_requirements,
};
use fs_material::fluid_acoustics::{FluidAcousticParameter, FluidAcousticState, FluidAcoustics};
use fs_material::fluid_transport::{FluidTransport, FluidTransportState, ThermalFluidTransport};
use fs_material::gas::{GasSpec, GasState};
use fs_material::liquid::{DYNAMIC_VISCOSITY_DIMS, LiquidPropertyRequirements};
use fs_material::liquid_temperature::LiquidViscosityBasis;
use fs_material::state_point::{
    MaterialPropertySelection, SPECIFIC_HEAT_CAPACITY_DIMS, ScalarAdmissibility,
    ScalarPropertyRequirement,
};
use fs_qty::semantic::{Composition, CompositionBasis, QuantityKind, SemanticType, ValueForm};
use fs_qty::{Density, Dims, QuantitySpec, Temperature, Velocity};

fn semantic(kind: QuantityKind) -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(kind, ValueForm::Static))
}
fn point(temperature: f64) -> QueryPoint {
    QueryPoint::new()
        .with_quantity(
            "T",
            semantic(QuantityKind::AbsoluteTemperature),
            temperature,
        )
        .unwrap()
        .with("pressure", 100_000.0)
        .unwrap()
        .with("phase-liquid", 1.0)
        .unwrap()
}
fn requirement(name: &str, dims: Dims) -> ScalarPropertyRequirement {
    ScalarPropertyRequirement::try_new(name, dims, ScalarAdmissibility::StrictlyPositive).unwrap()
}

#[test]
fn gas_composition_and_pressure_share_original_caloric_and_transport_state() {
    let composition = Composition::new(CompositionBasis::MoleFraction, vec![0.98, 0.02]).unwrap();
    let specs = [GasSpec::dry_air_ussa1976(), GasSpec::water_vapor_nist()];
    let a = GasState::try_new_mixture(&specs, &composition, 300.0, 100_000.0, Some(1)).unwrap();
    let b = GasState::try_new_mixture(&specs, &composition, 300.0, 200_000.0, Some(1)).unwrap();
    let av = a.thermoviscous_state().unwrap();
    let bv = b.thermoviscous_state().unwrap();
    assert_eq!(av.heat_capacity_ratio().to_bits(), a.gamma.to_bits());
    assert_eq!(av.transport(), a.transport_state().unwrap());
    assert_eq!(av.acoustic(), a.acoustic_state().unwrap());
    assert_eq!(
        av.acoustic().sound_speed_m_s(),
        bv.acoustic().sound_speed_m_s()
    );
    assert_eq!(
        av.transport().kinematic_viscosity_m2_s() / 2.0,
        bv.transport().kinematic_viscosity_m2_s()
    );
    assert_eq!(
        av.transport().thermal_diffusivity_m2_s() / 2.0,
        bv.transport().thermal_diffusivity_m2_s()
    );
    let borrowed: &dyn ThermoviscousFluid = &a;
    assert_eq!((&borrowed).thermoviscous_state().unwrap(), av);
    let mut stale = a;
    stale.gamma *= 1.01;
    assert!(stale.thermoviscous_state().is_err());
}

#[test]
fn no_thermal_coupling_is_explicit_and_incompatible_states_refuse() {
    let acoustic = FluidAcousticState::try_new(1000.0, 1500.0).unwrap();
    let transport = FluidTransportState::try_new(1000.0, 4000.0, 0.6, 0.001).unwrap();
    let zero = ThermoviscousFluidState::try_new(acoustic, transport, 4000.0).unwrap();
    assert_eq!(zero.heat_capacity_ratio() - 1.0, 0.0);
    assert_eq!(zero.thermoviscous_state().unwrap(), zero);
    for cv in [
        4001.0,
        0.0,
        -1.0,
        f64::NAN,
        f64::INFINITY,
        f64::MIN_POSITIVE,
    ] {
        assert!(ThermoviscousFluidState::try_new(acoustic, transport, cv).is_err());
    }
    for ratio in [0.99, 0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(
            ThermoviscousFluidState::from_heat_capacity_ratio(acoustic, transport, ratio).is_err()
        );
    }
    let other_density = FluidTransportState::try_new(999.0, 4000.0, 0.6, 0.001).unwrap();
    assert!(ThermoviscousFluidState::try_new(acoustic, other_density, 3900.0).is_err());
}

fn liquid_fixture() -> (
    MaterialCard,
    LiquidPropertyRequirements,
    ScalarPropertyRequirement,
) {
    let cp_quantity = QuantitySpec::dimensional(SPECIFIC_HEAT_CAPACITY_DIMS);
    let conductivity = semantic(QuantityKind::ThermalConductivity);
    let mut claims = ClaimSet::new();
    for (name, quantity, values) in [
        (
            "rho-source",
            QuantitySpec::dimensional(Density::DIMS),
            [1000.0, 960.0],
        ),
        ("cp-source", cp_quantity, [4000.0, 4200.0]),
        ("cv-source", cp_quantity, [3900.0, 3800.0]),
        ("conductivity-source", conductivity, [0.6, 0.65]),
        (
            "viscosity-source",
            QuantitySpec::dimensional(DYNAMIC_VISCOSITY_DIMS),
            [0.001, 0.0005],
        ),
        (
            "wave-source",
            QuantitySpec::dimensional(Velocity::DIMS),
            [1500.0, 1600.0],
        ),
    ] {
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
                    .with_quantity(
                        "T",
                        semantic(QuantityKind::AbsoluteTemperature),
                        300.0,
                        400.0,
                    )
                    .with("pressure", 100_000.0, 100_000.0)
                    .with("phase-liquid", 1.0, 1.0),
                interpolation: InterpolationPolicy::LinearInside,
                uncertainty: UncertaintyModel::Unstated,
                observations: vec![],
                provenance: Provenance {
                    source: "manufactured thermoviscous liquid".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "manufactured fluid".into(),
            phase: "liquid".into(),
            process: "numerical fixture".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    let liquid = LiquidPropertyRequirements::try_new(
        requirement("rho-source", Density::DIMS),
        requirement("cp-source", SPECIFIC_HEAT_CAPACITY_DIMS),
        ScalarPropertyRequirement::try_with_quantity(
            "conductivity-source",
            conductivity,
            ScalarAdmissibility::StrictlyPositive,
        )
        .unwrap(),
        requirement("viscosity-source", DYNAMIC_VISCOSITY_DIMS),
        LiquidViscosityBasis::Dynamic,
    )
    .unwrap();
    (
        card,
        liquid,
        requirement("cv-source", SPECIFIC_HEAT_CAPACITY_DIMS),
    )
}

#[test]
fn liquid_refresh_preserves_original_six_claims_and_independent_caloric_roles() {
    let (card, liquid, cv) = liquid_fixture();
    let parameter = FluidAcousticParameter::SoundSpeed(requirement("wave-source", Velocity::DIMS));
    let cold = resolve_thermoviscous_liquid_state_with_requirements(
        &card,
        &point(300.0),
        &parameter,
        &liquid,
        &cv,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let original = cold.clone();
    let hot = cold.at_query_point(&card, &point(400.0)).unwrap();
    assert_eq!(cold.state().heat_capacity_ratio(), 4000.0 / 3900.0);
    assert_eq!(hot.state().heat_capacity_ratio(), 4200.0 / 3800.0);
    assert_eq!(hot.state().acoustic().sound_speed_m_s(), 1600.0);
    assert_eq!(hot.state().transport().dynamic_viscosity_pa_s(), 0.0005);
    assert_eq!(hot.material().properties().len(), 6);
    assert_eq!(hot.material().query(), point(400.0));
    assert_eq!(hot.material().card_identity(), card.content_hash());
    assert_eq!(hot.acoustic_liquid().liquid_requirements(), &liquid);
    assert_eq!(hot.cv_requirement(), &cv);
    assert_ne!(cold.identity(), hot.identity());
    for property in hot.material().properties() {
        let previous = cold
            .material()
            .property(property.requirement().name())
            .unwrap();
        assert_eq!(property.requirement(), previous.requirement());
        assert_eq!(
            property.answer().receipt.selected,
            previous.answer().receipt.selected
        );
        assert_eq!(
            property.answer().receipt.query_point,
            hot.material().query_point()
        );
        card.claims()
            .verify_receipt(&property.answer().receipt)
            .unwrap();
    }
    assert_eq!(
        hot.thermal_transport_state().unwrap(),
        hot.acoustic_liquid().thermal_transport_state().unwrap()
    );
    assert!(cold.at_query_point(&card, &point(401.0)).is_err());
    assert!(
        cold.at_query_point(&card, &point(300.0).with("pressure", 200_000.0).unwrap())
            .is_err()
    );
    assert_eq!(cold, original);
    // Assigning the isobaric source to the isochoric role would silently turn
    // off thermal attenuation; even equal-valued roles need distinct claims.
    assert!(
        ResolvedThermoviscousLiquidState::try_from_acoustic(
            cold.acoustic_liquid().clone(),
            &liquid.properties()[1]
        )
        .is_err()
    );
    assert!(
        ResolvedThermoviscousLiquidState::try_from_acoustic(
            cold.acoustic_liquid().clone(),
            &requirement("absent-cv", SPECIFIC_HEAT_CAPACITY_DIMS)
        )
        .is_err()
    );
}
