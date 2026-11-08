//! G1/G3/G4: original material resistance, nonlinear/latent storage and rollback.

use super::*;
use fs_evidence::ValidityDomain;
use fs_gas::dynamics::nonideal::SharedWallLink;
use fs_gas::multicomponent::network::wall_exchange::{
    MaterialWallConduction, PhaseWallConductionSource, WallConductionEvaluation,
};
use fs_gas::multicomponent::wall_exchange::{WallCaloricConstraint, WallExchangeError};
use fs_matdb::{
    ClaimSelection, ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim,
    PropertyKey, PropertyValue, Provenance, QueryPoint, SelectionPolicy, UncertaintyModel,
};
use fs_material::phase::{EnthalpyPhaseKnot, EquilibriumEnthalpyPhaseCurve};
use fs_material::state_point::{MaterialPropertySelection, THERMAL_CONDUCTIVITY_DIMS};
use fs_material::thermal::conduction::ReservoirLayerGeometry;
use fs_material::thermal::source::ThermalMaterialPoint;
use fs_material::thermal::transport::storage::{
    HeatStorageAmount, HeatStorageComponentQuery, PreparedHeatStorage,
};
use fs_material::thermal::transport::{PreparedThermalTransport, ThermalTransportSourceQuery};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Area, Dims, QuantitySpec};
use std::sync::Arc;

fn heat_controls() -> WallControl {
    WallControl {
        relative_energy_tolerance: 1e-12,
        ..wall_controls()
    }
}

fn fixed_interval() -> StepControl {
    StepControl {
        maximum_halvings: 0,
        ..controls()
    }
}

fn temperature_spec() -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(
        QuantityKind::AbsoluteTemperature,
        ValueForm::Static,
    ))
}

// Manufactured original card: k(T)=a+bT, Cv(T)=c+dT, one kilogram.
// These are analytic test laws, not measured material property claims.
fn temperature_source(
    conductivity: [f64; 2],
    capacity: [f64; 2],
    span: [f64; 2],
) -> (Arc<PreparedThermalTransport>, MaterialWallCaloric<'static>) {
    let keys = [
        PropertyKey::new("density", Dims([-3, 1, 0, 0, 0, 0])),
        PropertyKey::new("constant_volume_specific_heat", Dims([2, 0, -2, -1, 0, 0])),
        PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS),
    ];
    let [lo, hi] = span;
    let mut claims = ClaimSet::new();
    for (key, [a, b]) in keys.iter().zip([[1000.0, 0.0], capacity, conductivity]) {
        claims
            .insert_claim(PropertyClaim {
                key: key.clone(),
                value: PropertyValue::Curve {
                    abscissa: "T".into(),
                    abscissa_dims: Temperature::DIMS,
                    knots: vec![(lo, a + b * lo), (hi, a + b * hi)],
                    dims: key.dims(),
                },
                validity: ValidityDomain::unconstrained().with_quantity(
                    "T",
                    temperature_spec(),
                    lo,
                    hi,
                ),
                interpolation: InterpolationPolicy::LinearInside,
                uncertainty: UncertaintyModel::Unstated,
                observations: vec![],
                provenance: Provenance {
                    source: "manufactured real-gas wall conduction".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic wall".into(),
            phase: "fixed solid".into(),
            process: "manufactured".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    let point = |t| {
        QueryPoint::new()
            .with_quantity("T", temperature_spec(), t)
            .unwrap()
    };
    let (lower, upper) = (point(lo), point(hi));
    let selection = ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly);
    let transport = PreparedThermalTransport::try_new(
        &card,
        ThermalTransportSourceQuery {
            density: &keys[0],
            specific_heat: &keys[1],
            conductivity: &keys[2],
            temperature_axis: "T",
            lower: &lower,
            upper: &upper,
            selections: [selection; 3],
        },
    )
    .unwrap();
    let storage = PreparedHeatStorage::try_new(&[HeatStorageComponentQuery {
        id: 0,
        card: &card,
        capacity: &keys[1],
        amount: HeatStorageAmount::MassKg(1.0),
        temperature_axis: "T",
        lower: &lower,
        upper: &upper,
        selection,
    }])
    .unwrap();
    (
        Arc::new(transport),
        MaterialWallCaloric::new(
            Arc::new(storage),
            WallCaloricConstraint::IsochoricInternalEnergy,
        ),
    )
}

// Cp=500 J/(kg K), Tm=550 K, latent heat=100 kJ/kg, rho=750 kg/m³.
// The separately authored conductivity is k(phi)=10+190 phi W/(m K).
fn phase_source(maximum_fraction: f64, anchor: f64) -> Arc<PhaseWallConductionSource> {
    let key = PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS);
    let mut claims = ClaimSet::new();
    let claim = |end: f64, scale: f64| PropertyClaim {
        key: key.clone(),
        value: PropertyValue::Curve {
            abscissa: "phi".into(),
            abscissa_dims: Dims::NONE,
            knots: vec![(0.0, scale * 10.0), (end, scale * (10.0 + 190.0 * end))],
            dims: key.dims(),
        },
        validity: ValidityDomain::unconstrained().with_quantity(
            "phi",
            QuantitySpec::dimensional(Dims::NONE),
            0.0,
            end,
        ),
        interpolation: InterpolationPolicy::LinearInside,
        uncertainty: UncertaintyModel::Unstated,
        observations: vec![],
        provenance: Provenance {
            source: "manufactured phase conductivity".into(),
            license: "test-only".into(),
            artifact: None,
        },
    };
    let selected = claims.insert_claim(claim(maximum_fraction, 1.0)).unwrap();
    if maximum_fraction < 1.0 {
        // A broad alternative must not replace the original pinned claim when
        // a later heat-transfer stage leaves its admitted phase interval.
        claims.insert_claim(claim(1.0, 10.0)).unwrap();
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic phase wall".into(),
            phase: "solid-liquid".into(),
            process: "manufactured".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    let curve = Arc::new(
        EquilibriumEnthalpyPhaseCurve::try_new(
            card.content_hash(),
            [
                (0.0, 300.0, 0.0),
                (125_000.0, 550.0, 0.0),
                (225_000.0, 550.0, 1.0),
                (450_000.0, 1000.0, 1.0),
            ]
            .map(|(h, t, phi)| EnthalpyPhaseKnot {
                specific_enthalpy_j_kg: h,
                temperature_k: t,
                liquid_mass_fraction: phi,
                bulk_density_kg_m3: 750.0,
            })
            .to_vec(),
        )
        .unwrap(),
    );
    let context = QueryPoint::new()
        .with_quantity("T", temperature_spec(), 550.0)
        .unwrap()
        .with_quantity("phi", QuantitySpec::dimensional(Dims::NONE), 0.0)
        .unwrap();
    Arc::new(
        PhaseWallConductionSource::new(
            curve,
            ThermalMaterialPoint {
                material: &card,
                context: &context,
                temperature_axis: "T",
                liquid_fraction_axis: Some("phi"),
            },
            &key,
            MaterialPropertySelection::PinnedByProperty(vec![(key.name().to_owned(), selected)]),
            anchor,
        )
        .unwrap(),
    )
}

fn phase_input(source: &Arc<PhaseWallConductionSource>, h: f64) -> RealGasWallContact<'static> {
    let caloric = MaterialWallCaloric::from_phase_curve(
        source.curve().clone(),
        wall_cell(h),
        [0.0, 450_000.0],
    )
    .unwrap();
    RealGasWallContact {
        contact: caloric.initial_contact(0.0, WallDrag::None).unwrap(),
        heat_storage: Some(caloric),
    }
}

fn region(input: RealGasWallContact<'static>, cell: usize, film: f64) -> SharedWall<'static> {
    let mut wall = shared_wall(input, &[film]);
    wall.ports[0].gas = SharedWallGas::Branch { branch: 0, cell };
    wall
}

fn temperature(walls: &RealGasWallOperator<'_>, index: usize) -> f64 {
    let storage = &walls.shared_wall(index).unwrap().storage;
    match &storage.heat_storage {
        Some(source) => source.temperature_k(storage.thermal).unwrap(),
        None => match storage.thermal {
            WallThermal::Lumped { temperature_k, .. } => temperature_k,
            _ => panic!("finite inventory required"),
        },
    }
}

fn phase_enthalpy(walls: &RealGasWallOperator<'_>, index: usize) -> f64 {
    let WallThermal::Caloric {
        mass_kg,
        specific_energy_j_kg,
        ..
    } = walls.shared_wall(index).unwrap().storage.thermal
    else {
        panic!("original phase inventory required");
    };
    assert_eq!(mass_kg, 0.75);
    specific_energy_j_kg
}

#[test]
fn integrated_profile_matches_pr_heat_and_one_nonlinear_inventory_with_solid_only_support() {
    let dt = 1e-6;
    let (initial_gas, final_gas, final_wall, face) = (400.0, 420.0, 600.0, 500.0);
    for minimum in [450.0, 510.0] {
        let mut flow = uniform(gas(), initial_gas, 0.0);
        let mass = 20.0 * flow.mesh().volumes()[0];
        let q = mass * (internal_energy(20.0, final_gas) - internal_energy(20.0, initial_gas));
        // Manufacture the original endpoint Fourier integral and film balance.
        // Cwall(T)=T for one kg gives exactly three q in its initial storage.
        let initial_wall = (final_wall * final_wall + 6.0 * q).sqrt();
        let film = q / (dt * (face - final_gas));
        let length = dt * 0.5 * (final_wall * final_wall - face * face) / q;
        let (transport, caloric) = temperature_source([0.0, 1.0], [0.0, 1.0], [minimum, 1000.0]);
        let layer = MaterialWallConduction::new(transport, Length::new(length), Area::new(1.0))
            .unwrap()
            .with_integrated_profile()
            .unwrap();
        let input = RealGasWallContact {
            contact: caloric.contact(initial_wall, 0.0, WallDrag::None).unwrap(),
            heat_storage: Some(caloric),
        };
        let mut walls =
            RealGasWallOperator::new(&flow, vec![contact(WallThermal::Adiabatic); 3], || false)
                .unwrap()
                .with_shared_walls(
                    vec![shared_wall(input, &[film; 3])],
                    SharedWallControl::default(),
                    || false,
                )
                .unwrap()
                .with_shared_wall_conduction(vec![vec![Some(layer.clone()); 3]], || false)
                .unwrap();
        let old = (
            flow.states().to_vec(),
            flow.primitives().to_vec(),
            flow.time(),
            walls.shared_wall(0).unwrap().clone(),
        );
        let outcome = flow.step_with_walls(
            Time::new(dt),
            fixed_interval(),
            &mut walls,
            heat_controls(),
            || false,
        );
        if minimum > face {
            // Both retained wall states are supported, but the actual solid
            // face is not. An endpoint-only k(Twall) shortcut would miss this.
            assert!(initial_wall > minimum && final_wall > minimum);
            assert!(matches!(outcome, Err(RealGasFlowError::StepRetries { .. })));
            assert_eq!(
                (
                    flow.states().to_vec(),
                    flow.primitives().to_vec(),
                    flow.time(),
                    walls.shared_wall(0).unwrap().clone()
                ),
                old
            );
            assert_eq!(walls.last_report(), None);
            assert_eq!(walls.last_shared_wall_report(), None);
            continue;
        }
        let report = outcome.unwrap();
        assert_eq!(report.elapsed, Time::new(dt));
        assert!(layer.source().unwrap().at_temperature(final_gas).is_err());
        assert!(layer.source().unwrap().at_temperature(face).is_ok());
        near(temperature(&walls, 0), final_wall, 8e-8);
        for (i, primitive) in flow.primitives().iter().enumerate() {
            near(primitive.temperature, final_gas, 8e-8);
            assert_eq!(flow.states()[i].density, old.0[i].density);
            assert_eq!(flow.states()[i].momentum, old.0[i].momentum);
            let actual_heat =
                flow.mesh().volumes()[i] * (flow.states()[i].energy - old.0[i].energy);
            let actual_face = primitive.temperature + actual_heat / (dt * film);
            near(actual_heat, q, 2e-5);
            near(actual_face, face, 8e-8);
            near(
                actual_heat,
                dt * (final_wall * final_wall - actual_face * actual_face) / (2.0 * length),
                3e-5,
            );
        }
        let endpoint_sampled = 1.0 / (1.0 / film + length / final_wall);
        assert!((dt * endpoint_sampled * (final_wall - final_gas) - q).abs() > 100.0);
        near(
            walls.last_report().unwrap().wall_storage_change_j,
            0.5 * (final_wall * final_wall - initial_wall * initial_wall),
            3e-5,
        );
        near(report.wall_heat, 3.0 * q, 3e-5);
        assert_eq!(walls.shared_wall(0).unwrap().ports, old.3.ports);
        assert_eq!(
            walls.shared_wall_conduction(0).unwrap(),
            vec![Some(layer.clone()); 3]
        );
        check_wall_balance(&flow, &walls, report);
    }
}

#[test]
fn implicit_phase_conductivity_changes_heat_on_one_latent_inventory_in_both_directions() {
    let (dt, film, length, area) = (1e-5, 1e8, 1e-5, 0.25);
    for (h0, gas0) in [(125_000.0, 800.0), (180_000.0, 450.0)] {
        let source = phase_source(1.0, h0);
        let layer = MaterialWallConduction::from_phase_source(
            source.clone(),
            Length::new(length),
            Area::new(area),
        )
        .unwrap()
        .with_evaluation(WallConductionEvaluation::ImplicitEndpoint);
        let mut flow = uniform(gas(), gas0, 0.0);
        let mass = 20.0 * flow.mesh().volumes()[0];
        let mut walls =
            RealGasWallOperator::new(&flow, vec![contact(WallThermal::Adiabatic); 3], || false)
                .unwrap()
                .with_shared_walls(
                    vec![shared_wall(phase_input(&source, h0), &[film; 3])],
                    SharedWallControl::default(),
                    || false,
                )
                .unwrap()
                .with_shared_wall_conduction(vec![vec![Some(layer); 3]], || false)
                .unwrap();
        // Independent PR inversion and phase inventory in the physical heat
        // coordinate. The bracket remains within the original latent plateau.
        let gas_at = |q: f64| {
            root(gas0.min(550.0), gas0.max(550.0), |t| {
                mass * (internal_energy(20.0, t) - internal_energy(20.0, gas0)) - q
            })
        };
        let bound = if gas0 > 550.0 {
            -0.75 * (225_000.0 - h0) / 3.0
        } else {
            0.75 * (h0 - 125_000.0) / 3.0
        };
        let q = root(bound.min(0.0), bound.max(0.0), |q| {
            let h = h0 - 3.0 * q / 0.75;
            let k = 10.0 + 190.0 * (h - 125_000.0) / 100_000.0;
            let g = 1.0 / (1.0 / film + length / (area * k));
            q - dt * g * (550.0 - gas_at(q))
        });
        let expected_h = h0 - 3.0 * q / 0.75;
        let expected_gas = gas_at(q);
        assert!((125_000.0..225_000.0).contains(&expected_h));
        let initial_k = 10.0 + 190.0 * (h0 - 125_000.0) / 100_000.0;
        let initial_g = 1.0 / (1.0 / film + length / (area * initial_k));
        let lagged = root(gas0.min(550.0), gas0.max(550.0), |t| {
            mass * (internal_energy(20.0, t) - internal_energy(20.0, gas0))
                - dt * initial_g * (550.0 - t)
        });
        assert!((lagged - expected_gas).abs() > 0.1);
        if gas0 > 550.0 {
            assert!(
                q.abs() > (dt * initial_g * (550.0 - gas0)).abs(),
                "a frozen initial-flux bracket excludes this endpoint"
            );
        }
        let report = flow
            .step_with_walls(
                Time::new(dt),
                fixed_interval(),
                &mut walls,
                heat_controls(),
                || false,
            )
            .unwrap();
        assert_eq!(report.elapsed, Time::new(dt));
        near(phase_enthalpy(&walls, 0), expected_h, 2e-5);
        assert_eq!(temperature(&walls, 0), 550.0);
        for primitive in flow.primitives() {
            near(primitive.temperature, expected_gas, 5e-8);
        }
        near(report.wall_heat, 3.0 * q, 2e-5);
        near(
            walls.last_report().unwrap().wall_storage_change_j,
            -3.0 * q,
            2e-5,
        );
        let accepted = walls.shared_wall(0).unwrap();
        assert!(Arc::ptr_eq(
            accepted
                .storage
                .heat_storage
                .as_ref()
                .unwrap()
                .phase_curve()
                .unwrap(),
            source.curve()
        ));
        assert!(
            accepted
                .ports
                .iter()
                .all(|port| port.conductance_w_per_k == film)
        );
        check_wall_balance(&flow, &walls, report);
    }
}

#[test]
fn two_sided_radial_material_links_refresh_original_capacity_and_phase_conductivity() {
    let (dt, interface) = (1e-5, 1e8);
    let source = phase_source(1.0, 155_000.0);
    let (transport, caloric) = temperature_source([100.0, 0.1], [0.0, 1.0], [300.0, 1000.0]);
    let hot = RealGasWallContact {
        contact: caloric.contact(750.0, 0.0, WallDrag::None).unwrap(),
        heat_storage: Some(caloric),
    };
    let layers = [
        Some(
            MaterialWallConduction::from_geometry(
                transport,
                ReservoirLayerGeometry::cylindrical_face(0.1, 0.100001, 0.5).unwrap(),
            )
            .unwrap(),
        ),
        Some(
            MaterialWallConduction::from_phase_source_with_geometry(
                source.clone(),
                ReservoirLayerGeometry::spherical_face(0.2, 0.200001)
                    .unwrap()
                    .with_area_fraction(0.25)
                    .unwrap(),
            )
            .unwrap()
            .with_evaluation(WallConductionEvaluation::ImplicitEndpoint),
        ),
    ];
    let mut flow = uniform(gas(), 550.0, 0.0);
    let mut walls =
        RealGasWallOperator::new(&flow, vec![contact(WallThermal::Adiabatic); 3], || false)
            .unwrap()
            .with_shared_walls(
                vec![
                    region(hot, 0, 0.0),
                    region(phase_input(&source, 155_000.0), 1, 0.0),
                ],
                SharedWallControl::default(),
                || false,
            )
            .unwrap()
            .with_shared_wall_links(
                vec![SharedWallLink {
                    walls: [0, 1],
                    conductance_w_per_k: interface,
                }],
                || false,
            )
            .unwrap()
            .with_shared_wall_link_conduction(vec![layers.clone()], || false)
            .unwrap();
    let initial_gas = (flow.states().to_vec(), flow.primitives().to_vec());
    let (mut ta, mut h) = (750.0_f64, 155_000.0);
    for step in 1..=2 {
        // Independent cylindrical and spherical resistance integrals. The
        // first k is lagged at this stage, the second is implicit in latent h.
        let ra = ((0.100001_f64 - 0.1) / 0.1).ln_1p()
            / (core::f64::consts::TAU * 0.5 * (100.0 + 0.1 * ta));
        let q = root(0.0, 0.75 * (225_000.0 - h), |q| {
            let hot = (ta * ta - 2.0 * q).sqrt();
            let cold_h = h + q / 0.75;
            let kb = 10.0 + 190.0 * (cold_h - 125_000.0) / 100_000.0;
            let rb = (1.0 / 0.2 - 1.0 / 0.200001) / (4.0 * core::f64::consts::PI * 0.25 * kb);
            let g = 1.0 / (1.0 / interface + ra + rb);
            q - dt * g * (hot - 550.0)
        });
        ta = (ta * ta - 2.0 * q).sqrt();
        h += q / 0.75;
        let report = flow
            .step_with_walls(
                Time::new(dt),
                fixed_interval(),
                &mut walls,
                heat_controls(),
                || false,
            )
            .unwrap();
        near(temperature(&walls, 0), ta, 8e-8);
        near(phase_enthalpy(&walls, 1), h, 3e-5);
        assert_eq!(temperature(&walls, 1), 550.0);
        assert_eq!(
            (flow.states().to_vec(), flow.primitives().to_vec()),
            initial_gas
        );
        assert_eq!(report.wall_heat, 0.0);
        near(
            walls.last_report().unwrap().wall_storage_change_j,
            0.0,
            2e-7,
        );
        near(
            walls.last_report().unwrap().gas_wall_storage_residual_j,
            0.0,
            3e-7,
        );
        let receipt = walls.shared_wall_link_reports().unwrap()[0];
        near(receipt.heat_a_to_b_j, q, 2e-5);
        near(receipt.wall_energy_change_j[0], -q, 2e-5);
        near(receipt.wall_energy_change_j[1], q, 2e-5);
        near(receipt.energy_residual_j, 0.0, 2e-7);
        assert_eq!(walls.shared_wall_link_conduction(0).unwrap(), &layers);
        assert_eq!(
            walls.shared_wall_links().unwrap()[0].conductance_w_per_k,
            interface
        );
        near(
            0.5 * (ta * ta - 750.0 * 750.0) + 0.75 * (h - 155_000.0),
            0.0,
            3e-5,
        );
        assert_eq!(flow.time(), Time::new(f64::from(step) * dt));
    }
}

#[test]
fn later_links_cannot_publish_a_wall_outside_an_earlier_pinned_conductivity_domain() {
    for on_link in [false, true] {
        for evaluation in [
            WallConductionEvaluation::InitialState,
            WallConductionEvaluation::ImplicitEndpoint,
        ] {
            let source = phase_source(0.001, 125_000.0);
            let layer = MaterialWallConduction::from_phase_source(
                source.clone(),
                Length::new(1e-5),
                Area::new(0.25),
            )
            .unwrap()
            .with_evaluation(evaluation);
            let mut flow = uniform(gas(), 550.0, 0.0);
            let room = |t| {
                contact(WallThermal::Lumped {
                    temperature_k: t,
                    heat_capacity_j_per_k: 400.0,
                    conductance_w_per_k: 0.0,
                })
            };
            let mut walls =
                RealGasWallOperator::new(&flow, vec![contact(WallThermal::Adiabatic); 3], || false)
                    .unwrap()
                    .with_shared_walls(
                        vec![
                            region(
                                phase_input(&source, 125_000.0),
                                0,
                                if on_link { 0.0 } else { 1e8 },
                            ),
                            region(room(550.0), 1, 0.0),
                            region(room(800.0), 2, 0.0),
                        ],
                        SharedWallControl::default(),
                        || false,
                    )
                    .unwrap()
                    .with_shared_wall_links(
                        vec![
                            SharedWallLink {
                                walls: [0, 1],
                                conductance_w_per_k: 1e8,
                            },
                            SharedWallLink {
                                walls: [2, 0],
                                conductance_w_per_k: 1e8,
                            },
                        ],
                        || false,
                    )
                    .unwrap();
            walls = if on_link {
                walls
                    .with_shared_wall_link_conduction(
                        vec![[Some(layer), None], [None, None]],
                        || false,
                    )
                    .unwrap()
            } else {
                walls
                    .with_shared_wall_conduction(
                        vec![vec![Some(layer)], vec![None], vec![None]],
                        || false,
                    )
                    .unwrap()
            };
            // A short first interval is supported. In the next interval the
            // gas/first-link stages stay at 550 K, then the last link advances
            // phi beyond the retained claim, although a broad alternative exists.
            flow.step_with_walls(
                Time::new(1e-9),
                fixed_interval(),
                &mut walls,
                heat_controls(),
                || false,
            )
            .unwrap();
            assert!(phase_enthalpy(&walls, 0) > 125_000.0);
            assert!(phase_enthalpy(&walls, 0) < 125_100.0);
            assert!(source.at_specific_enthalpy(125_200.0).is_err());
            let saved = (
                flow.states().to_vec(),
                flow.primitives().to_vec(),
                flow.time(),
                (0..3)
                    .map(|i| walls.shared_wall(i).unwrap().clone())
                    .collect::<Vec<_>>(),
                walls.last_report(),
                walls.last_shared_wall_report(),
                walls.shared_wall_link_reports().unwrap().to_vec(),
            );
            let error = flow
                .step_with_walls(
                    Time::new(1e-6),
                    fixed_interval(),
                    &mut walls,
                    heat_controls(),
                    || false,
                )
                .unwrap_err();
            let RealGasFlowError::StepRetries { last, .. } = error else {
                panic!("{error:?}");
            };
            assert!(matches!(*last, RealGasFlowError::Thermal(ref e)
                if matches!(**e, WallExchangeError::PhaseConduction(_))));
            assert_eq!(
                (
                    flow.states().to_vec(),
                    flow.primitives().to_vec(),
                    flow.time(),
                    (0..3)
                        .map(|i| walls.shared_wall(i).unwrap().clone())
                        .collect::<Vec<_>>(),
                    walls.last_report(),
                    walls.last_shared_wall_report(),
                    walls.shared_wall_link_reports().unwrap().to_vec()
                ),
                saved
            );
        }
    }
}
