//! Original variable-density liquid sources meet through native finite mixing.
use super::*;
use fs_couple::liquid_transfer::transfer_interface_mass_streams;
use fs_lbm::material::free_surface::enthalpy::FreeSurfaceEnthalpyError;
use fs_material::liquid::surface::ResolvedLiquidSurfaceState;
use fs_material::phase::{HeatCapacitySourceQuery, ResolvedHeatCapacityCurve};
use fs_material::state_point::{
    ScalarAdmissibility, ScalarPropertyRequirement, resolve_material_state_point,
};

const GRID: f64 = 1.0 / 1024.0;
const PACKET: f64 = 3.0 / 33_554_432.0;

#[derive(Clone, Copy, PartialEq)]
enum Context {
    Same,
    Pressure,
    OptionalSet,
    OptionalClaim,
}

fn source_point(temperature: f64, pressure: f64) -> QueryPoint {
    QueryPoint::new()
        .with("T", temperature)
        .unwrap()
        .with("pressure", pressure)
        .unwrap()
}

fn with_thermal_group(
    plateau: bool,
    datum: f64,
    velocities: [[f64; 3]; 3],
    context: Context,
    run: impl for<'s, 'm> FnOnce(
        &mut MaterialFreeSurfaceEnthalpy3<'s, 'm>,
        &mut MaterialFreeSurfaceEnthalpy3<'s, 'm>,
        &mut MaterialFreeSurfaceEnthalpy3<'s, 'm>,
        &CancelGate,
    ),
) {
    let mut claims = ClaimSet::new();
    for (name, dims, values) in [
        (
            "density",
            Dims([-3, 1, 0, 0, 0, 0]),
            [1024.0, if plateau { 1024.0 } else { 896.0 }, 768.0],
        ),
        (
            "specific_heat_capacity",
            Dims([2, 0, -2, -1, 0, 0]),
            [2048.0; 3],
        ),
        ("thermal_conductivity", Dims([1, 1, -3, -1, 0, 0]), [0.5; 3]),
        ("dynamic_viscosity", Dims([-1, 1, -1, 0, 0, 0]), [0.1; 3]),
        (
            "electrical_conductivity",
            Dims([-3, -1, 3, 0, 0, 2]),
            [2.0; 3],
        ),
    ] {
        let quantity = if name == "thermal_conductivity" {
            QuantitySpec::semantic(SemanticType::new(
                QuantityKind::ThermalConductivity,
                ValueForm::Static,
            ))
        } else {
            QuantitySpec::dimensional(dims)
        };
        let ranges = if name == "electrical_conductivity" && context == Context::OptionalClaim {
            vec![(300.0, 400.0), (400.0, 500.0)]
        } else {
            vec![(300.0, 500.0)]
        };
        for (lower, upper) in ranges {
            claims
                .insert_claim(PropertyClaim {
                    key: PropertyKey::with_quantity(name, quantity),
                    value: PropertyValue::Curve {
                        abscissa: "T".into(),
                        abscissa_dims: Temperature::DIMS,
                        knots: [300.0, 400.0, 500.0].into_iter().zip(values).collect(),
                        dims,
                    },
                    validity: ValidityDomain::unconstrained()
                        .with("T", lower, upper)
                        .with("pressure", 100000.0, 101000.0),
                    uncertainty: UncertaintyModel::Unstated,
                    interpolation: InterpolationPolicy::LinearInside,
                    observations: vec![],
                    provenance: Provenance {
                        source: "synthetic thermal-density mixing".into(),
                        license: "test-only".into(),
                        artifact: None,
                    },
                })
                .unwrap();
        }
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "one synthetic liquid".into(),
            phase: "liquid".into(),
            process: "density mixing fixture".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    let chart = ResolvedHeatCapacityCurve::try_new_exact(
        &card,
        HeatCapacitySourceQuery {
            heat_capacity: &PropertyKey::with_quantity(
                "specific_heat_capacity",
                QuantitySpec::dimensional(Dims([2, 0, -2, -1, 0, 0])),
            ),
            density: &PropertyKey::with_quantity(
                "density",
                QuantitySpec::dimensional(Dims([-3, 1, 0, 0, 0, 0])),
            ),
            temperature_axis: "T",
            lower: &source_point(300.0, 100000.0),
            upper: &source_point(500.0, 100000.0),
            selection: fs_matdb::SelectionPolicy::SingleClaimOnly,
        },
        SolidLiquidPhase::Liquid,
        datum,
    )
    .unwrap();
    let temperatures = [300.0, 500.0, if plateau { 350.0 } else { 400.0 }];
    let laws: [LiquidTemperatureLaw<'_>; 3] = core::array::from_fn(|index| {
        let pressure = if context == Context::Pressure && index == 1 {
            101000.0
        } else {
            100000.0
        };
        let point = source_point(temperatures[index], pressure);
        let original = LiquidTemperatureLaw::try_new(
            &card,
            &point,
            "T",
            LiquidViscosityBasis::Dynamic,
            MaterialPropertySelection::SingleClaimOnly,
        )
        .unwrap();
        let mut requirements = original
            .reference_state()
            .material()
            .properties()
            .iter()
            .map(|p| p.requirement().clone())
            .collect::<Vec<_>>();
        if context == Context::OptionalClaim || (context == Context::OptionalSet && index == 0) {
            requirements.push(
                ScalarPropertyRequirement::try_new(
                    "electrical_conductivity",
                    Dims([-3, -1, 3, 0, 0, 2]),
                    ScalarAdmissibility::StrictlyPositive,
                )
                .unwrap(),
            );
        }
        let material = resolve_material_state_point(
            &card,
            &point,
            &requirements,
            MaterialPropertySelection::SingleClaimOnly,
        )
        .unwrap();
        LiquidTemperatureLaw::try_from_material(&card, material, "T", LiquidViscosityBasis::Dynamic)
            .unwrap()
    });
    let sources: [ResolvedLiquidSurfaceState; 3] = core::array::from_fn(|index| {
        ResolvedLiquidSurfaceState::try_from_material(
            laws[index].reference_state().material().clone(),
            LiquidViscosityBasis::Dynamic,
            LiquidCapillarity::Neglected,
        )
        .unwrap()
    });
    let gate = CancelGate::new_clock_free();
    fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default()).scope(|arena| {
        let cx = fs_exec::Cx::new(
            &gate,
            arena,
            fs_exec::StreamKey {
                seed: 91,
                kernel_id: 29,
                tile: 0,
                iteration: 0,
            },
            fs_exec::Budget::INFINITE,
            fs_exec::ExecMode::Deterministic,
        );
        let make = |index: usize| {
            let dx = if index == 1 { 2.0 * GRID } else { GRID };
            let h = datum + (temperatures[index] - 300.0) * 2048.0;
            plan_free_surface_from_material(&sources[index], dx, GRID, 0.1)
                .unwrap()
                .free_surface3_from_cells(
                    MaterialFreeSurfaceDomain3 {
                        cells: [8; 3],
                        gravity_m_s2: [0.0; 3],
                        contact: ContactModel3::Neutral,
                        solid_tiles: &[],
                        max_cells: 512,
                    },
                    &[MaterialFreeSurfaceInitialCell3 {
                        coordinates: CELL,
                        fill_fraction: if index == 2 { 0.25 } else { 0.75 },
                        velocity_m_s: velocities[index],
                    }],
                    &cx,
                )
                .unwrap()
                .with_enthalpy(
                    &laws[index],
                    chart.curve(),
                    FreeSurfaceEnthalpyControls::default(),
                    |_| h,
                    &gate,
                )
                .unwrap()
        };
        let mut cold = make(0);
        let mut hot = make(1);
        let mut receiver = make(2);
        run(&mut cold, &mut hot, &mut receiver, &gate);
    });
}

fn packet(
    donor: &MaterialFreeSurfaceEnthalpy3<'_, '_>,
    receiver: &MaterialFreeSurfaceEnthalpy3<'_, '_>,
    mass_kg: f64,
) -> LiquidMassTransfer {
    let origin = |liquid: &MaterialFreeSurfaceEnthalpy3<'_, '_>| {
        let dx = liquid.flow().scaling().cell_size_m();
        CELL.map(|coordinate| 0.5 - (coordinate as f64 + 0.5) * dx)
    };
    LiquidMassTransfer {
        donor_origin_m: origin(donor),
        receiver_origin_m: origin(receiver),
        ..transfer(mass_kg)
    }
}

fn inventory(liquid: &MaterialFreeSurfaceEnthalpy3<'_, '_>) -> [f64; 6] {
    let sample = liquid
        .flow()
        .sample_at_si(CELL[0], CELL[1], CELL[2])
        .unwrap()
        .unwrap();
    let mass = liquid.flow().ledger_mass_kg().unwrap();
    let kinetic = 0.5 * mass * sample.velocity_m_s.iter().map(|v| v * v).sum::<f64>();
    [
        mass,
        liquid.total_enthalpy_j().unwrap() + kinetic,
        mass * sample.velocity_m_s[0],
        mass * sample.velocity_m_s[1],
        mass * sample.velocity_m_s[2],
        sample.pressure_perturbation_pa,
    ]
}

#[test]
fn g1_balanced_hot_cold_densities_and_native_mass_units_share_one_original_receiver() {
    with_thermal_group(
        false,
        0.0,
        [[0.0; 3]; 3],
        Context::Same,
        |cold, hot, receiver, gate| {
            assert!(!cold.shares_mass_transfer_source(receiver));
            assert!(!hot.shares_mass_transfer_source(receiver));
            assert!(cold.shares_mass_transfer_material(receiver));
            assert!(hot.shares_mass_transfer_material(receiver));
            let before = [inventory(cold), inventory(hot), inventory(receiver)];
            let transfers = [
                packet(cold, receiver, PACKET),
                packet(hot, receiver, PACKET),
            ];
            let old_h = receiver.specific_enthalpy_at(CELL).unwrap();
            let report = transfer_interface_mass_streams(
                &mut [&mut *cold, &mut *hot],
                receiver,
                &transfers,
                2,
                controls(),
                gate,
            )
            .unwrap();
            let after = [inventory(cold), inventory(hot), inventory(receiver)];
            for component in 0..5 {
                near(
                    after.iter().map(|x| x[component]).sum(),
                    before.iter().map(|x| x[component]).sum(),
                    if component == 1 { 1e-12 } else { 1e-18 },
                );
            }
            assert_eq!(receiver.specific_enthalpy_at(CELL), Some(old_h));
            assert_eq!(receiver.material_at(CELL).unwrap().density_kg_m3(), 896.0);
            assert_eq!(cold.flow().scaling().material().density_kg_m3(), 1024.0);
            assert_eq!(hot.flow().scaling().material().density_kg_m3(), 768.0);
            assert_eq!(report.receiver.mixing_heat_j, 0.0);
            assert_eq!(report.point_m, [0.5; 3]);
            assert_eq!(report.donors.len(), 2);
            for donor in &report.donors {
                assert_eq!(donor.exchange.mass_change_kg, -PACKET);
            }
            near(
                report.reference_volumes.withdrawn_m3,
                PACKET / 1024.0 + PACKET / 768.0,
                1e-24,
            );
            near(
                report.reference_volumes.received_m3,
                2.0 * PACKET / 896.0,
                1e-24,
            );
            assert!(report.reference_volumes.net_change_m3 < 0.0);
            for sample in after {
                assert_eq!(sample[5], 0.0);
            }
            assert_eq!(cold.accepted_steps(), 0);
            assert_eq!(hot.accepted_steps(), 0);
            assert_eq!(receiver.accepted_steps(), 0);
        },
    );
}

#[test]
fn g1_different_density_streams_thermalize_only_their_actual_relative_motion() {
    let velocities = [
        [0.035, -0.012, 0.003],
        [-0.025, 0.018, -0.004],
        [0.005, -0.007, 0.009],
    ];
    let mut reference_heat = None;
    for (datum, boost) in [(0.0, [0.0; 3]), (1e8, [0.01, 0.005, -0.002])] {
        let velocities = velocities.map(|u| core::array::from_fn(|a| u[a] + boost[a]));
        with_thermal_group(
            true,
            datum,
            velocities,
            Context::Same,
            |cold, hot, receiver, gate| {
                let before = [inventory(cold), inventory(hot), inventory(receiver)];
                let old_h = receiver.specific_enthalpy_at(CELL).unwrap();
                let sample = |owner: &MaterialFreeSurfaceEnthalpy3<'_, '_>| {
                    owner
                        .flow()
                        .sample_at_si(CELL[0], CELL[1], CELL[2])
                        .unwrap()
                        .unwrap()
                        .velocity_m_s
                };
                let u = [sample(cold), sample(hot), sample(receiver)];
                let masses = [2.0 * PACKET, PACKET, before[2][0]];
                let total_mass = masses.iter().sum::<f64>();
                let mut expected_heat = 0.0;
                for i in 0..3 {
                    for j in i + 1..3 {
                        expected_heat += 0.5 * masses[i] * masses[j] / total_mass
                            * (0..3).map(|a| (u[i][a] - u[j][a]).powi(2)).sum::<f64>();
                    }
                }
                let transfers = [
                    packet(cold, receiver, masses[0]),
                    packet(hot, receiver, masses[1]),
                ];
                let report = transfer_interface_mass_streams(
                    &mut [&mut *cold, &mut *hot],
                    receiver,
                    &transfers,
                    2,
                    controls(),
                    gate,
                )
                .unwrap();
                near(report.receiver.mixing_heat_j, expected_heat, 2e-23);
                assert!(report.receiver.mixing_heat_j > 0.0);
                if let Some(heat) = reference_heat {
                    near(report.receiver.mixing_heat_j, heat, 2e-23);
                } else {
                    reference_heat = Some(report.receiver.mixing_heat_j);
                }
                let expected_h = old_h
                    + (masses[0] * (datum - old_h)
                        + masses[1] * (datum + 409600.0 - old_h)
                        + expected_heat)
                        / total_mass;
                near(
                    receiver.specific_enthalpy_at(CELL).unwrap(),
                    expected_h,
                    1e-7,
                );
                assert_eq!(receiver.material_at(CELL).unwrap().density_kg_m3(), 1024.0);
                let after = [inventory(cold), inventory(hot), inventory(receiver)];
                for component in 0..5 {
                    near(
                        after.iter().map(|x| x[component]).sum(),
                        before.iter().map(|x| x[component]).sum(),
                        if component == 1 { 1e-10 } else { 1e-18 },
                    );
                }
                assert_eq!(after[2][5], before[2][5]);
                assert!(report.enthalpy_kinetic_defect_j.abs() <= controls().energy_tolerance_j);
            },
        );
    }
}

#[test]
fn g3_single_different_density_packet_matches_one_native_stream_batch() {
    with_thermal_group(
        true,
        0.0,
        [[0.0; 3], [0.02, 0.0, 0.0], [0.0; 3]],
        Context::Same,
        |_, hot, receiver, gate| {
            let mut batch_donor = hot.clone();
            let mut batch_receiver = receiver.clone();
            let transfer = packet(hot, receiver, PACKET);
            let single =
                transfer_interface_mass(hot, receiver, transfer, controls(), gate).unwrap();
            let batch = transfer_interface_mass_streams(
                &mut [&mut batch_donor],
                &mut batch_receiver,
                &[transfer],
                1,
                controls(),
                gate,
            )
            .unwrap();
            assert_eq!(single.donor, batch.donors[0]);
            assert_eq!(single.receiver, batch.receiver);
            assert_eq!(single.reference_volumes, batch.reference_volumes);
            unchanged(hot, &batch_donor);
            unchanged(receiver, &batch_receiver);
            assert!(single.receiver.mixing_heat_j > 0.0);
        },
    );
}

#[test]
fn g4_unmodeled_final_density_capacity_and_cancellation_preserve_every_original() {
    with_thermal_group(
        false,
        0.0,
        [[0.0; 3]; 3],
        Context::Same,
        |cold, hot, receiver, gate| {
            let before = [cold.clone(), hot.clone(), receiver.clone()];
            let unbalanced = [
                packet(cold, receiver, PACKET),
                packet(hot, receiver, 2.0 * PACKET),
            ];
            assert!(matches!(
                transfer_interface_mass_streams(
                    &mut [&mut *cold, &mut *hot],
                    receiver,
                    &unbalanced,
                    2,
                    controls(),
                    gate,
                ),
                Err(FreeSurfaceEnthalpyError::Invalid(
                    "free-surface heat needs pure liquid at the original fixed caloric density"
                ))
            ));
            for (owner, original) in [
                (&*cold, &before[0]),
                (&*hot, &before[1]),
                (&*receiver, &before[2]),
            ] {
                unchanged(owner, original);
            }
            let valid = [
                packet(cold, receiver, PACKET),
                packet(hot, receiver, PACKET),
            ];
            assert!(
                transfer_interface_mass_streams(
                    &mut [&mut *cold, &mut *hot],
                    receiver,
                    &valid,
                    1,
                    controls(),
                    gate,
                )
                .is_err()
            );
            let overflowing = [
                packet(cold, receiver, 6.0 * PACKET),
                packet(hot, receiver, 6.0 * PACKET),
            ];
            assert!(
                transfer_interface_mass_streams(
                    &mut [&mut *cold, &mut *hot],
                    receiver,
                    &overflowing,
                    2,
                    controls(),
                    gate,
                )
                .is_err()
            );
            let cancelled = CancelGate::new_clock_free();
            cancelled.request();
            assert!(
                transfer_interface_mass_streams(
                    &mut [&mut *cold, &mut *hot],
                    receiver,
                    &valid,
                    2,
                    controls(),
                    &cancelled,
                )
                .is_err()
            );
            for (owner, original) in [
                (&*cold, &before[0]),
                (&*hot, &before[1]),
                (&*receiver, &before[2]),
            ] {
                unchanged(owner, original);
            }
        },
    );
}

#[test]
fn g4_complete_optional_claims_and_nonthermal_context_cannot_be_discarded() {
    for context in [
        Context::Pressure,
        Context::OptionalSet,
        Context::OptionalClaim,
    ] {
        with_thermal_group(
            true,
            0.0,
            [[0.0; 3]; 3],
            context,
            |cold, hot, receiver, gate| {
                let before = [cold.clone(), hot.clone(), receiver.clone()];
                if context == Context::OptionalSet {
                    assert!(!cold.shares_mass_transfer_material(receiver));
                } else {
                    assert!(!hot.shares_mass_transfer_material(receiver));
                }
                let transfers = [
                    packet(cold, receiver, 2.0 * PACKET),
                    packet(hot, receiver, PACKET),
                ];
                assert!(
                    transfer_interface_mass_streams(
                        &mut [&mut *cold, &mut *hot],
                        receiver,
                        &transfers,
                        2,
                        controls(),
                        gate,
                    )
                    .is_err()
                );
                for (owner, original) in [
                    (&*cold, &before[0]),
                    (&*hot, &before[1]),
                    (&*receiver, &before[2]),
                ] {
                    unchanged(owner, original);
                }
            },
        );
    }
}
