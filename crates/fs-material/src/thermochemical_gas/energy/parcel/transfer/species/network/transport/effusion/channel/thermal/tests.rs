//! Original material capacity, finite storage and shared channel transactions.
use super::super::super::tests::{controls, energy, gas, interface, near};
use super::*;
use crate::thermal::transport::volumetric::{
    PreparedVolumetricThermalTransport, VOLUMETRIC_HEAT_CAPACITY_DIMS, VolumetricThermalPathQuery,
};
use crate::thermochemical_gas::energy::parcel::transfer::species::effusion::channel::{
    ClausingControls, CylindricalMolecularGeometry,
};
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSelection, ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim,
    PropertyKey, PropertyValue, Provenance, QueryPoint, SelectionPolicy, UncertaintyModel,
};
use fs_qty::{
    Dims, QuantitySpec, Temperature,
    semantic::{QuantityKind, SemanticType, ValueForm},
};

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
}
/// Explicitly authored isochoric C_vol=A+B*T; no measurement claim.
pub(in super::super) fn material(
    a: f64,
    b: f64,
    span: [f64; 2],
) -> PreparedVolumetricThermalTransport {
    material_with_conductivity(a, b, span, [1.0; 2])
}

/// Independently authored affine conductivity on the same supported interval.
pub(super) fn material_with_conductivity(
    a: f64,
    b: f64,
    span: [f64; 2],
    k: [f64; 2],
) -> PreparedVolumetricThermalTransport {
    let capacity = PropertyKey::new(
        "isochoric-volumetric-capacity",
        VOLUMETRIC_HEAT_CAPACITY_DIMS,
    );
    let conductivity = PropertyKey::with_quantity(
        "thermal-conductivity",
        quantity(QuantityKind::ThermalConductivity),
    );
    let mut claims = ClaimSet::new();
    for (key, values, dims) in [
        (
            capacity.clone(),
            span.map(|t| (t, a + b * t)),
            VOLUMETRIC_HEAT_CAPACITY_DIMS,
        ),
        (
            conductivity.clone(),
            [(span[0], k[0]), (span[1], k[1])],
            Dims([1, 1, -3, -1, 0, 0]),
        ),
    ] {
        claims
            .insert_claim(PropertyClaim {
                key,
                value: PropertyValue::Curve {
                    abscissa: "temperature".into(),
                    abscissa_dims: Temperature::DIMS,
                    knots: values.to_vec(),
                    dims,
                },
                validity: ValidityDomain::unconstrained().with_quantity(
                    "temperature",
                    quantity(QuantityKind::AbsoluteTemperature),
                    span[0],
                    span[1],
                ),
                interpolation: InterpolationPolicy::LinearInside,
                uncertainty: UncertaintyModel::Unstated,
                observations: vec![],
                provenance: Provenance {
                    source: "synthetic reversible isochoric wall fixture".into(),
                    license: "CC0-1.0".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "authored finite wall".into(),
            phase: "solid".into(),
            process: "fixed volume and composition".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    PreparedVolumetricThermalTransport::try_new(
        &card,
        VolumetricThermalPathQuery {
            volumetric_heat_capacity: &capacity,
            conductivity: &conductivity,
            temperature_axis: "temperature",
            lower: &point(span[0]),
            upper: &point(span[1]),
            selections: [ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly); 2],
        },
    )
    .unwrap()
}

pub(super) fn channel(length: f64) -> DiffuseMolecularChannel {
    DiffuseMolecularChannel::try_new(
        interface(0.0),
        CylindricalMolecularGeometry::try_new(
            1e-8,
            length,
            1e-4,
            "synthetic diffuse finite-wall channel",
        )
        .unwrap(),
    )
    .unwrap()
}

pub(super) fn link(channel: &DiffuseMolecularChannel) -> FiniteChannelNetworkLink<'_> {
    FiniteChannelNetworkLink {
        id: 7,
        parcels: [0, 1],
        channel,
        area_m2: 1e-6,
        wall: 9,
    }
}

pub(super) fn endpoints(gases: &BTreeMap<u64, ClosedGasParcel>) -> BTreeMap<u64, GasEnergySolve> {
    gases.keys().map(|&id| (id, energy())).collect()
}

#[test]
fn repeated_collision_accommodation_drives_both_gases_and_one_finite_wall() {
    use fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R;
    let source = material(1e6, 0.0, [200.0, 1800.0]);
    let law = FixedVolumeThermalLaw::try_new(&source, 1000.0, 200.0)
        .unwrap()
        .with_sensible_entropy();
    let walls = BTreeMap::from([(
        9,
        FixedVolumeChannelWall::try_new_isochoric(&law, 1e-6, 900.0).unwrap(),
    )]);
    let path = channel(2e-8)
        .with_thermal_accommodation(
            0.25,
            ClausingControls {
                initial_rings: 8,
                maximum_rings: 512,
                maximum_iterations: 256,
                absolute_tolerance: 1e-6,
                relative_tolerance: 1e-4,
            },
            "synthetic common per-collision energy accommodation",
            || false,
        )
        .unwrap();
    let gases = BTreeMap::from([
        (0, gas(400.0, 100_000.0, 0.0)),
        (1, gas(400.0, 100_000.0, 0.0)),
    ]);
    let accepted = step_candidate(
        &gases,
        &walls,
        GasTransportPaths::default(),
        &[link(&path)],
        0.001,
        &endpoints(&gases),
        controls(),
        || false,
    )
    .unwrap();
    let model = path.thermal_accommodation().unwrap();
    let probability = model.reflected_heat_probability() + model.transmitted_heat_probability();
    let incident =
        0.001 * 1e-6 * 100_000.0 / (2.0 * core::f64::consts::PI * 0.028 * R * 400.0).sqrt();
    let heat = incident * 2.0 * R * 500.0 * probability;
    for id in [0, 1] {
        let n = gases[&id].mass_kg() / 0.028;
        near(
            accepted.parcels()[&id].mass_kg(),
            gases[&id].mass_kg(),
            1e-18,
        );
        near(
            accepted.parcels()[&id].state().state().temperature,
            400.0 + heat / (1.5 * R * n),
            3e-8,
        );
    }
    near(
        accepted.walls()[&9].temperature_k(),
        900.0 - 2.0 * heat,
        1e-10,
    );
    near(accepted.step().wall_energy_change_j, -2.0 * heat, 1e-9);
    assert!(accepted.step().energy_residual_j.abs() <= controls().exchange.energy_tolerance_j);
    assert!(accepted.step().total_entropy_change_j_k > 0.0);
}

#[test]
fn nonlinear_finite_hot_and_cold_walls_close_actual_energy_and_entropy() {
    let source = material(1e6, 1000.0, [200.0, 1800.0]);
    let path = channel(1e-8);
    let gases = BTreeMap::from([
        (0, gas(400.0, 100_000.0, 0.0)),
        (1, gas(400.0, 100_000.0, 0.0)),
    ]);
    let original = gases.clone();
    for t0 in [250.0, 900.0] {
        let mut previous = None;
        for (density, datum) in [(1000.0, 200.0), (8000.0, 1200.0)] {
            let law = FixedVolumeThermalLaw::try_new(&source, density, datum)
                .unwrap()
                .with_sensible_entropy();
            let walls = BTreeMap::from([(
                9,
                FixedVolumeChannelWall::try_new_isochoric(&law, 1e-6, t0).unwrap(),
            )]);
            let result = step_candidate(
                &gases,
                &walls,
                GasTransportPaths::default(),
                &[link(&path)],
                0.001,
                &endpoints(&gases),
                controls(),
                || false,
            )
            .unwrap();
            let t1 = result.walls()[&9].temperature_k();
            let row = &result.step().walls[0];
            let expected_energy = 1e-6 * (t1 - t0) * (1e6 + 500.0 * (t0 + t1));
            let expected_entropy = 1e-6 * (1e6 * ((t1 - t0) / t0).ln_1p() + 1000.0 * (t1 - t0));
            near(row.energy_change_j, expected_energy, 2e-10);
            near(row.entropy_change_j_k, expected_entropy, 2e-12);
            near(row.energy_change_j, row.heat_into_wall_j, 1e-9);
            assert!(row.entropy_change_j_k < row.heat_into_wall_j / t0);
            assert!(result.step().total_entropy_change_j_k > 0.0);
            assert!(
                result.step().energy_residual_j.abs() <= controls().exchange.energy_tolerance_j
            );
            assert_eq!(result.walls()[&9].volume_m3(), walls[&9].volume_m3());
            assert_eq!(result.walls()[&9].mass_kg(), walls[&9].mass_kg());
            assert_eq!(result.walls()[&9].law().identity(), law.identity());
            near(walls[&9].temperature_k(), t0, 1e-10);
            if t0 > 400.0 {
                assert!(t1 < t0);
            } else {
                assert!(t1 > t0);
            }
            if let Some(old) = previous {
                near(t1, old, 1e-10);
            }
            previous = Some(t1);
        }
    }
    assert_eq!(gases, original);
}

#[test]
fn opposing_channels_pool_one_shared_wall_before_the_only_reconstruction() {
    let source = material(1000.0, 0.0, [499.0, 501.0]);
    let law = FixedVolumeThermalLaw::try_new(&source, 1000.0, 500.0)
        .unwrap()
        .with_sensible_entropy();
    let walls = BTreeMap::from([(
        9,
        FixedVolumeChannelWall::try_new_isochoric(&law, 1e-6, 500.0).unwrap(),
    )]);
    let path = channel(1e-8);
    // Equal original molecular incident rates at unequal gas temperatures.
    let p_hot = 200_000.0;
    let gases = BTreeMap::from([
        (0, gas(800.0, p_hot, 0.0)),
        (1, gas(800.0, p_hot, 0.0)),
        (2, gas(200.0, 100_000.0, 0.0)),
        (3, gas(200.0, 100_000.0, 0.0)),
    ]);
    let a = link(&path);
    let b = FiniteChannelNetworkLink {
        id: 8,
        parcels: [2, 3],
        ..a
    };
    let run = |links: &[FiniteChannelNetworkLink<'_>]| {
        step_candidate(
            &gases,
            &walls,
            GasTransportPaths::default(),
            links,
            0.001,
            &endpoints(&gases),
            controls(),
            || false,
        )
        .unwrap()
    };
    let first = run(&[a, b]);
    let reversed = run(&[b, a]);
    assert_eq!(first.parcels(), reversed.parcels());
    assert_eq!(first.step(), reversed.step());
    assert_eq!(first.step().walls.len(), 1);
    near(first.walls()[&9].temperature_k(), 500.0, 1e-8);
    for stream in &first.step().channels {
        let heat = stream.wall.heat_into_wall_j;
        assert!(
            law.temperature_k(heat / walls[&9].mass_kg()).is_err(),
            "each unpooled heat would leave the narrow source interval"
        );
    }
    assert!(first.step().total_entropy_change_j_k > 0.0);
}

#[test]
fn stiff_wall_refusal_preserves_both_finite_inventories() {
    let source = material(1e6, 1000.0, [200.0, 1800.0]);
    let law = FixedVolumeThermalLaw::try_new(&source, 1000.0, 200.0)
        .unwrap()
        .with_sensible_entropy();
    let walls = BTreeMap::from([(
        9,
        FixedVolumeChannelWall::try_new_isochoric(&law, 1e-11, 900.0).unwrap(),
    )]);
    let path = channel(1e-8);
    let gases = BTreeMap::from([
        (0, gas(400.0, 100_000.0, 0.0)),
        (1, gas(400.0, 100_000.0, 0.0)),
    ]);
    let original = gases.clone();
    let old_energy = walls[&9].specific_energy_j_kg();
    assert!(matches!(
        step_candidate(
            &gases,
            &walls,
            GasTransportPaths::default(),
            &[link(&path)],
            0.001,
            &endpoints(&gases),
            controls(),
            || false
        ),
        Err(FiniteChannelError::WallDomain { .. } | FiniteChannelError::WallOvershoot { .. })
    ));
    assert_eq!(gases, original);
    assert_eq!(walls[&9].specific_energy_j_kg(), old_energy);
    assert_eq!(walls[&9].temperature_k(), 900.0);
}

#[test]
fn cancelled_late_candidate_cannot_publish_only_the_gas_endpoint() {
    let source = material(1e6, 0.0, [200.0, 1800.0]);
    let law = FixedVolumeThermalLaw::try_new(&source, 1000.0, 200.0)
        .unwrap()
        .with_sensible_entropy();
    let walls = BTreeMap::from([(
        9,
        FixedVolumeChannelWall::try_new_isochoric(&law, 1e-6, 600.0).unwrap(),
    )]);
    let path = channel(1e-8);
    let gases = BTreeMap::from([
        (0, gas(400.0, 100_000.0, 0.0)),
        (1, gas(400.0, 100_000.0, 0.0)),
    ]);
    let original = gases.clone();
    let mut calls = 0;
    step_candidate(
        &gases,
        &walls,
        GasTransportPaths::default(),
        &[link(&path)],
        0.001,
        &endpoints(&gases),
        controls(),
        || {
            calls += 1;
            false
        },
    )
    .unwrap();
    let mut second_calls = 0;
    assert!(matches!(
        step_candidate(
            &gases,
            &walls,
            GasTransportPaths::default(),
            &[link(&path)],
            0.001,
            &endpoints(&gases),
            controls(),
            || {
                second_calls += 1;
                second_calls == calls
            }
        ),
        Err(FiniteChannelError::Cancelled)
    ));
    assert_eq!(gases, original);
    assert_eq!(walls[&9].temperature_k(), 600.0);
}

#[test]
fn ballistic_channel_keeps_the_finite_wall_out_of_its_heat_and_domain_queries() {
    let source = material(1e6, 0.0, [1800.0, 2000.0]);
    let law = FixedVolumeThermalLaw::try_new(&source, 1000.0, 1800.0)
        .unwrap()
        .with_sensible_entropy();
    let walls = BTreeMap::from([(
        9,
        FixedVolumeChannelWall::try_new_isochoric(&law, 1e-6, 1900.0).unwrap(),
    )]);
    let path = channel(0.0);
    let gases = BTreeMap::from([
        (0, gas(600.0, 150_000.0, 0.0)),
        (1, gas(300.0, 100_000.0, 0.0)),
    ]);
    let result = step_candidate(
        &gases,
        &walls,
        GasTransportPaths::default(),
        &[link(&path)],
        0.001,
        &endpoints(&gases),
        controls(),
        || false,
    )
    .unwrap();
    assert_eq!(
        result.walls()[&9].specific_energy_j_kg(),
        walls[&9].specific_energy_j_kg()
    );
    assert_eq!(result.step().heat_into_walls_j, 0.0);
    assert_eq!(result.step().wall_entropy_change_j_k, 0.0);
}

#[test]
fn unavailable_material_entropy_and_missing_wall_are_explicit_refusals() {
    let source = material(1e6, 0.0, [200.0, 1800.0]);
    let law = FixedVolumeThermalLaw::try_new(&source, 1000.0, 200.0).unwrap();
    assert!(matches!(
        FixedVolumeChannelWall::try_new_isochoric(&law, 1e-6, 400.0),
        Err(FiniteChannelError::Invalid(_))
    ));
    let path = channel(1e-8);
    let gases = BTreeMap::from([
        (0, gas(400.0, 100_000.0, 0.0)),
        (1, gas(400.0, 100_000.0, 0.0)),
    ]);
    assert!(matches!(
        step_candidate(
            &gases,
            &BTreeMap::new(),
            GasTransportPaths::default(),
            &[link(&path)],
            0.001,
            &endpoints(&gases),
            controls(),
            || false
        ),
        Err(FiniteChannelError::Invalid(_))
    ));
}

fn gas_with_caloric_span(name: &str, temperature_k: f64, span: [f64; 2]) -> ClosedGasParcel {
    use crate::thermochemical_gas::{
        ThermochemicalGasModel,
        tests::{point as gas_point, species_in_domain, transport_in_domain},
    };
    use fs_thermochem::{Composition, CompositionBasis, mixture::FrozenIdealGasMixtureModelV1};
    let model = ThermochemicalGasModel::new(
        FrozenIdealGasMixtureModelV1::new(
            vec![species_in_domain(name, 0.028, 2.5, 0.0, 0.0, span)],
            Composition::new(CompositionBasis::MoleFraction, vec![1.0]).unwrap(),
        )
        .unwrap(),
        vec![transport_in_domain(
            name,
            1.7e-5,
            0.025,
            [200.0, 2000.0],
            [1.0, 1e7],
        )],
        None,
    )
    .unwrap();
    ClosedGasParcel::new(model, gas_point(temperature_k, 100_000.0), 0.001).unwrap()
}

fn channel_for_catalog(
    catalog: &crate::thermochemical_gas::ThermochemicalGasModel,
    length_m: f64,
) -> DiffuseMolecularChannel {
    use crate::thermochemical_gas::energy::parcel::transfer::species::effusion::{
        EffusiveInterface, EffusiveSpeciesTransmission,
    };
    let mouth = interface(0.0);
    let transmissions = catalog
        .thermodynamics()
        .components()
        .iter()
        .map(|species| {
            EffusiveSpeciesTransmission::try_new(
                species.species().clone(),
                1.0,
                "synthetic source-specific channel",
            )
            .unwrap()
        })
        .collect();
    DiffuseMolecularChannel::try_new(
        EffusiveInterface::try_new(
            catalog.clone(),
            transmissions,
            mouth.support(),
            mouth.regime().clone(),
        )
        .unwrap(),
        CylindricalMolecularGeometry::try_new(
            1e-8,
            length_m,
            1e-4,
            "synthetic source-specific channel geometry",
        )
        .unwrap(),
    )
    .unwrap()
}

#[test]
fn species_specific_accommodation_drives_closed_and_vacuum_finite_walls() {
    use crate::thermochemical_gas::energy::parcel::transfer::species::effusion::channel::SpeciesThermalAccommodation;
    use crate::thermochemical_gas::{ThermochemicalGasModel, tests::point as gas_point};
    use fs_qty::chemistry::SpeciesId;
    use fs_thermochem::{
        Composition, CompositionBasis, UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R,
        mixture::FrozenIdealGasMixtureModelV1,
    };
    let narrow = gas_with_caloric_span("A", 400.0, [200.0, 600.0]);
    let wide = gas_with_caloric_span("B", 400.0, [200.0, 1800.0]);
    let catalog = ThermochemicalGasModel::new(
        FrozenIdealGasMixtureModelV1::new(
            vec![
                narrow.model().thermodynamics().components()[0].clone(),
                wide.model().thermodynamics().components()[0].clone(),
            ],
            Composition::new(CompositionBasis::MoleFraction, vec![0.5, 0.5]).unwrap(),
        )
        .unwrap(),
        vec![
            narrow.model().transport_materials()[0].clone(),
            wide.model().transport_materials()[0].clone(),
        ],
        None,
    )
    .unwrap();
    let mixture =
        ClosedGasParcel::new(catalog.clone(), gas_point(400.0, 100_000.0), 0.001).unwrap();
    let path = channel_for_catalog(&catalog, 1e-8)
        .with_species_thermal_accommodation(
            [("A", 0.0), ("B", 1.0)]
                .into_iter()
                .map(|(name, alpha)| {
                    SpeciesThermalAccommodation::try_new(
                        SpeciesId::new(name).unwrap(),
                        alpha,
                        "synthetic species-specific energy accommodation",
                    )
                    .unwrap()
                })
                .collect(),
            ClausingControls::default(),
            || false,
        )
        .unwrap();
    let source = material(1e6, 0.0, [200.0, 1800.0]);
    let law = FixedVolumeThermalLaw::try_new(&source, 1000.0, 200.0)
        .unwrap()
        .with_sensible_entropy();
    let walls = BTreeMap::from([(
        9,
        FixedVolumeChannelWall::try_new_isochoric(&law, 1e-6, 900.0).unwrap(),
    )]);
    let gases = BTreeMap::from([(0, mixture.clone()), (1, mixture.clone())]);
    let mut solve = endpoints(&gases);
    for control in solve.values_mut() {
        control.maximum_temperature_k = 600.0;
    }
    let closed = step_candidate(
        &gases,
        &walls,
        GasTransportPaths::default(),
        &[link(&path)],
        0.001,
        &solve,
        controls(),
        || false,
    )
    .unwrap();
    // A and B have equal molar masses and cp=5R/2. Half of the incident
    // molecules are elastic A; only B transfers heat to or from this wall.
    let incident_b =
        0.001 * 1e-6 * 50_000.0 / (2.0 * core::f64::consts::PI * 0.028 * R * 400.0).sqrt();
    let tau = path.geometry().transmission_probability();
    let returned = incident_b * path.geometry().reflected_probability() * 2.0 * R * 500.0;
    let transmitted = incident_b * path.geometry().wall_transmitted_probability() * 2.0 * R * 500.0;
    let n0 = mixture.mass_kg() / 0.028;
    for id in [0, 1] {
        near(
            closed.parcels()[&id].state().state().temperature,
            400.0 + (returned + transmitted) / (1.5 * R * n0),
            3e-8,
        );
    }
    near(
        closed.step().wall_energy_change_j,
        -2.0 * (returned + transmitted),
        1e-9,
    );
    near(
        closed.walls()[&9].temperature_k(),
        900.0 - 2.0 * (returned + transmitted),
        1e-9,
    );
    assert!(closed.walls()[&9].temperature_k() > 600.0);
    assert!(closed.step().energy_residual_j.abs() <= controls().exchange.energy_tolerance_j);
    let reflected = &closed.step().channels[0].wall.reflected_species;
    assert_eq!(reflected.len(), 2);
    assert_eq!(reflected[0].mass_kg, reflected[1].mass_kg);
    assert!(reflected[0].mass_kg[0] > 0.0);

    let vacuum_gases = BTreeMap::from([(0, mixture)]);
    let vacuum_solve = BTreeMap::from([(0, solve[&0])]);
    let escaped = vacuum::step_candidate(
        &vacuum_gases,
        &walls,
        GasTransportPaths::default(),
        &[vacuum::FiniteVacuumChannelLink {
            id: 7,
            parcel: 0,
            channel: &path,
            area_m2: 1e-6,
            wall: 9,
        }],
        0.001,
        &vacuum_solve,
        controls(),
        || false,
    )
    .unwrap();
    let n_out = 2.0 * incident_b * tau;
    near(
        escaped.parcels()[&0].state().state().temperature,
        (1.5 * R * n0 * 400.0 - n_out * 2.0 * R * 400.0 + returned) / (1.5 * R * (n0 - n_out)),
        3e-8,
    );
    near(
        escaped.step().gas.audit.escaped_energy_j,
        n_out * 2.0 * R * 400.0 + transmitted,
        1e-11,
    );
    near(
        escaped.step().wall_energy_change_j,
        -returned - transmitted,
        1e-9,
    );
    assert!(escaped.walls()[&9].temperature_k() > 600.0);
    assert!(escaped.step().energy_residual_j.abs() <= controls().exchange.energy_tolerance_j);
    let escaped_species = &escaped.step().gas.outlets[0].escaped_species;
    assert_eq!(escaped_species.len(), 2);
    near(
        escaped_species[0].mass_kg,
        escaped_species[1].mass_kg,
        1e-18,
    );
    assert_eq!(walls[&9].temperature_k(), 900.0);
}

#[test]
fn final_wall_temperature_must_fit_every_actually_exposed_original_caloric_source() {
    let source = material(1000.0, 0.0, [200.0, 1800.0]);
    let law = FixedVolumeThermalLaw::try_new(&source, 1000.0, 200.0)
        .unwrap()
        .with_sensible_entropy();
    let walls = BTreeMap::from([(
        9,
        FixedVolumeChannelWall::try_new_isochoric(&law, 3e-7, 500.0).unwrap(),
    )]);
    for cold_upper in [1800.0, 600.0] {
        let cold = gas_with_caloric_span("A", 400.0, [200.0, cold_upper]);
        let hot = gas_with_caloric_span("B", 1500.0, [200.0, 1800.0]);
        let cold_path = channel_for_catalog(cold.model(), 1e-8);
        let hot_path = channel_for_catalog(hot.model(), 1e-8);
        cold_path
            .admit_wall_for_parcel(&cold, 500.0, || false)
            .unwrap();
        hot_path
            .admit_wall_for_parcel(&hot, 500.0, || false)
            .unwrap();
        let gases = BTreeMap::from([(0, cold.clone()), (1, cold), (2, hot.clone()), (3, hot)]);
        let links = [
            link(&cold_path),
            FiniteChannelNetworkLink {
                id: 8,
                parcels: [2, 3],
                channel: &hot_path,
                area_m2: 1e-6,
                wall: 9,
            },
        ];
        // Cold A must have a valid full caloric bracket even when its source
        // stops at 600 K. Hot B still needs its original 1800 K upper bound.
        let mut solve = endpoints(&gases);
        for id in [0, 1] {
            solve.get_mut(&id).unwrap().maximum_temperature_k = 600.0;
        }
        let result = step_candidate(
            &gases,
            &walls,
            GasTransportPaths::default(),
            &links,
            0.001,
            &solve,
            controls(),
            || false,
        );
        if cold_upper == 1800.0 {
            let accepted = result.unwrap();
            assert!(accepted.walls()[&9].temperature_k() > 600.0);
            assert!(accepted.walls()[&9].temperature_k() < 800.0);
            assert!(accepted.parcels()[&0].state().state().temperature < 600.0);
        } else {
            assert!(matches!(
                result,
                Err(FiniteChannelError::Transport(GasTransportError::Effusion(
                    _
                )))
            ));
        }
        assert_eq!(walls[&9].temperature_k(), 500.0);
    }
}

#[test]
fn endpoint_wall_admission_skips_absent_species_and_inactive_wall_paths() {
    use crate::thermochemical_gas::ThermochemicalGasModel;
    use fs_thermochem::{Composition, CompositionBasis, mixture::FrozenIdealGasMixtureModelV1};
    let narrow = gas_with_caloric_span("A", 400.0, [200.0, 600.0]);
    let wide = gas_with_caloric_span("B", 1500.0, [200.0, 1800.0]);
    let catalog = ThermochemicalGasModel::new(
        FrozenIdealGasMixtureModelV1::new(
            vec![
                narrow.model().thermodynamics().components()[0].clone(),
                wide.model().thermodynamics().components()[0].clone(),
            ],
            Composition::new(CompositionBasis::MoleFraction, vec![0.5, 0.5]).unwrap(),
        )
        .unwrap(),
        vec![
            narrow.model().transport_materials()[0].clone(),
            wide.model().transport_materials()[0].clone(),
        ],
        None,
    )
    .unwrap();
    let active = channel_for_catalog(&catalog, 1e-8);
    active
        .admit_wall_for_parcel(&wide, 650.0, || false)
        .unwrap();
    assert!(
        active
            .admit_wall_for_parcel(&narrow, 650.0, || false)
            .is_err()
    );
    assert!(matches!(
        active.admit_wall_for_parcel(&wide, 650.0, || true),
        Err(crate::thermochemical_gas::energy::GasEnergyError::Cancelled)
    ));
    let ballistic = channel_for_catalog(&catalog, 0.0);
    ballistic
        .admit_wall_for_parcel(&narrow, f64::NAN, || false)
        .unwrap();

    let source = material(1e6, 0.0, [200.0, 1800.0]);
    let law = FixedVolumeThermalLaw::try_new(&source, 1000.0, 200.0)
        .unwrap()
        .with_sensible_entropy();
    let walls = BTreeMap::from([(
        9,
        FixedVolumeChannelWall::try_new_isochoric(&law, 1e-6, 900.0).unwrap(),
    )]);
    let gases = BTreeMap::from([(0, narrow.clone()), (1, narrow)]);
    let closed = FiniteChannelNetworkLink {
        area_m2: 0.0,
        ..link(&active)
    };
    let mut solve = endpoints(&gases);
    for control in solve.values_mut() {
        control.maximum_temperature_k = 600.0;
    }
    let result = step_candidate(
        &gases,
        &walls,
        GasTransportPaths::default(),
        &[closed],
        0.001,
        &solve,
        controls(),
        || false,
    )
    .unwrap();
    assert_eq!(result.parcels(), &gases);
    assert_eq!(
        result.walls()[&9].specific_energy_j_kg(),
        walls[&9].specific_energy_j_kg()
    );
    assert_eq!(result.step().heat_into_walls_j, 0.0);
}
