//! Manufactured conservation/contact-transport checks, not an instrument fit.
#[path = "dry/temperature.rs"]
mod temperature;

use super::*;
use fs_couple::bowed_string::runtime::reciprocal::bank::multiple::{
    BowDryContactHeat, DryReceiverHeatPosition,
    performance::bow_stations::{BowStationKnot, BowStationProgram},
    stations::BowStationMotionLimits,
};
use fs_couple::render::schedule::force::coupled::contact::friction::{
    interface::DryFrictionBinding,
    multiple::layer::thermal::dry::{DryContactFootprint, DryContactHeat, DryThermalCell},
};
use fs_material::thermal::source::ThermalMaterialPoint;
use fs_tribo::{
    ApplicabilityRange as Range, DryFrictionApplicability, DryInterfaceSystemCard, InputAuthority,
    InterfaceMedium, InterfaceSystemRef,
};

const LENGTH: f64 = 0.002;
const WIDTH: f64 = 0.001;
const DEPTH: f64 = 1e-5;

fn dry_source(maximum_temperature: f64) -> DryFrictionBinding {
    DryFrictionBinding::new(
        DryInterfaceSystemCard::new(
            InterfaceSystemRef::new(
                "synthetic/string->finite-receiver",
                "synthetic/dry-history",
                "manufactured heat balance",
                InputAuthority::SyntheticFixture,
                InterfaceMedium::Dry,
            )
            .unwrap(),
            law(),
            DryFrictionApplicability::new(
                Range::new(280.0, maximum_temperature).unwrap(),
                Range::new(0.001, 1e6).unwrap(),
                Range::new(0.0, 10.0).unwrap(),
            )
            .unwrap(),
        )
        .unwrap(),
        300.0,
        AREA,
    )
    .unwrap()
}
fn dry_bank(cx: &Cx<'_>, gate: &CancelGate, maximum_temperature: f64) -> MultiBowStringBank {
    reacting_surface(joint(cx, gate, 7800.0), gate)
        .with_bow_interface(0, dry_source(maximum_temperature))
        .unwrap()
}
fn transport(card: &MaterialCard) -> ResolvedMaterialStatePoint {
    let requirements = [
        ("density", Density::DIMS),
        ("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS),
    ]
    .into_iter()
    .map(|(name, dims)| {
        ScalarPropertyRequirement::try_new(name, dims, ScalarAdmissibility::StrictlyPositive)
            .unwrap()
    })
    .collect::<Vec<_>>();
    resolve_material_state_point(
        card,
        &point(300.0),
        &requirements,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap()
}
fn id(side: usize, cell: usize) -> CellId {
    CellId((100 + 100 * side + cell) as u64)
}
fn dry_network(
    curve: &EquilibriumEnthalpyPhaseCurve,
    temperatures: [f64; 2],
    depths: [f64; 2],
) -> DepositionSimulation<'_> {
    let mut events = Vec::new();
    for side in 0..2 {
        for cell in 0..3 {
            let volume = LENGTH * WIDTH * depths[side];
            events.push(TimedEvent {
                time_s: 0.0,
                event: Event::Birth(DepositedCell {
                    id: id(side, cell),
                    segment_id: cell as u64,
                    body_id: side as u64,
                    material_id: MaterialId(0),
                    mass_kg: DENSITY * volume,
                    specific_h_j_kg: CP * (temperatures[side] - 300.0),
                    geometry: CellGeometry {
                        start_m: [cell as f64 * LENGTH, 0.0, side as f64 * DEPTH],
                        end_m: [(cell + 1) as f64 * LENGTH, 0.0, side as f64 * DEPTH],
                        volume_m3: volume,
                        surface_area_m2: 2.0 * (LENGTH * WIDTH + depths[side] * (LENGTH + WIDTH)),
                    },
                    boundary: Boundary {
                        ambient_k: 300.0,
                        htc_w_m2_k: 0.0,
                        reservoir: None,
                        source_w: 0.0,
                    },
                }),
            });
        }
    }
    DepositionSimulation::new(
        events,
        &[(MaterialId(0), curve as &dyn ThermalLaw)],
        DepositionBudget {
            max_cells: 6,
            max_links: 1,
            max_events: 8,
            max_substeps: 8,
        },
    )
    .unwrap()
    .advance_to(
        0.0,
        StepConfig {
            max_dt_s: 1e-5,
            energy_tolerance_j: 1e-13,
        },
        || Ok(()),
    )
    .unwrap()
    .0
}
fn route<'a>(
    card: &'a MaterialCard,
    curve: &'a EquilibriumEnthalpyPhaseCurve,
    bundle: &'a ResolvedMaterialStatePoint,
    context: &'a QueryPoint,
    depths: [f64; 2],
) -> BowDryContactHeat<'a> {
    BowDryContactHeat {
        heat: DryContactHeat {
            port: 0,
            surfaces: std::array::from_fn(|side| {
                (0..3)
                    .map(|cell| DryThermalCell {
                        cell: id(side, cell),
                        material: ThermalMaterialPoint {
                            material: card,
                            context,
                            temperature_axis: "T",
                            liquid_fraction_axis: None,
                        },
                        caloric: curve,
                        transport: bundle,
                        conductivity_property: "thermal_conductivity",
                        density_property: "density",
                        interval_m: [cell as f64 * LENGTH, (cell + 1) as f64 * LENGTH],
                        width_m: WIDTH,
                        half_depth_m: 0.5 * depths[side],
                    })
                    .collect()
            }),
            footprint: DryContactFootprint {
                intervals_m: [[0.002, 0.004]; 2],
                widths_m: [WIDTH; 2],
            },
        },
        string_coordinate_origin_m: -0.097,
        receiver: DryReceiverHeatPosition::Fixed,
    }
}
fn solve_controls() -> LayerThermalSolve {
    LayerThermalSolve {
        maximum_cells: 6,
        maximum_mechanical_trials: 512,
        difference_step_j: 1e-9,
        ..controls(6)
    }
}
fn guesses(heat: &DepositionSimulation<'_>) -> Vec<f64> {
    heat.cells()
        .values()
        .map(|c| c.mass_kg * c.specific_h_j_kg)
        .collect()
}
fn motion_limits() -> BowStationMotionLimits {
    BowStationMotionLimits {
        maximum_speed_m_s: 200.0,
        maximum_step_m: 0.005,
        max_update_terms: 1_048_576,
    }
}

#[test]
fn stationary_contact_conducts_without_slip_and_release_switches_off_both_fluxes() {
    with_cx(|cx, gate| {
        let (card, caloric) = coating(false, false);
        let bundle = transport(&card);
        let context = point(300.0);
        let depths = [DEPTH, 2.0 * DEPTH];
        let routes = [route(&card, caloric.curve(), &bundle, &context, depths)];
        for load in [0.0, 0.1] {
            let mut bank = dry_bank(cx, gate, 340.0);
            let forces = vec![0.0; bank.system().mode_count()];
            let mut heat = dry_network(caloric.curve(), [310.0, 290.0], depths);
            let saved = heat.cells().clone();
            let initial = total_heat(&heat);
            let frame = bank
                .step_with_dry_contact_thermal_motion(
                    &mut heat,
                    &routes,
                    &[],
                    solve_controls(),
                    &guesses_from_cells(&saved),
                    &forces,
                    &[drive(0.0, load), drive(0.0, 0.0)],
                    None,
                    None,
                    None,
                    gate,
                    other,
                )
                .unwrap();
            assert_eq!(frame.routed_dry_heat_j, 0.0);
            assert!(frame.total_energy_residual_j.abs() < 2e-12);
            assert!((total_heat(&heat) - initial).abs() < 2e-12);
            let contact = &frame.dry_contacts[0];
            if load == 0.0 {
                assert_eq!(contact.interface_temperature_k, None);
                assert!(contact.cell_heat_j.is_empty());
                assert_eq!(heat.cells(), &saved);
            } else {
                let dt = bank.system().sample_period_s();
                let conductance = 0.2 * AREA / (0.5 * (depths[0] + depths[1]));
                let capacities = depths.map(|d| DENSITY * AREA * d * CP);
                let difference =
                    20.0 / (1.0 + dt * conductance * (1.0 / capacities[0] + 1.0 / capacities[1]));
                let exchanged = dt * conductance * difference;
                assert!(
                    (heat.temperature_k(id(0, 1)).unwrap() - (310.0 - exchanged / capacities[0]))
                        .abs()
                        < 1e-7
                );
                assert!(
                    (heat.temperature_k(id(1, 1)).unwrap() - (290.0 + exchanged / capacities[1]))
                        .abs()
                        < 1e-7
                );
                assert!(contact.cell_heat_j.iter().any(|(_, q)| *q < 0.0));
                assert!(contact.cell_heat_j.iter().any(|(_, q)| *q > 0.0));
                assert_eq!(
                    bank.bow_interface(0).unwrap().temperature_kelvin(),
                    contact.interface_temperature_k.unwrap()
                );
            }
            for side in 0..2 {
                for cell in [0, 2] {
                    assert_eq!(heat.cells()[&id(side, cell)], saved[&id(side, cell)]);
                }
            }
        }
    });
}
fn guesses_from_cells(cells: &BTreeMap<CellId, DepositedCell>) -> Vec<f64> {
    cells
        .values()
        .map(|c| c.mass_kg * c.specific_h_j_kg)
        .collect()
}

#[test]
fn finite_capacity_and_half_cell_geometry_determine_heat_split_without_changing_dry_mechanics() {
    with_cx(|cx, gate| {
        let (card, caloric) = coating(false, false);
        let bundle = transport(&card);
        let context = point(300.0);
        let mut reference = dry_bank(cx, gate, 340.0);
        let forces = vec![0.0; reference.system().mode_count()];
        let drives = [drive(0.05, 0.1), drive(0.0, 0.0)];
        let expected = reference.step(&forces, &drives).unwrap().clone();
        for depths in [[DEPTH, DEPTH], [DEPTH, 2.0 * DEPTH]] {
            let mut bank = dry_bank(cx, gate, 340.0);
            let mut heat = dry_network(caloric.curve(), [300.0; 2], depths);
            let routes = [route(&card, caloric.curve(), &bundle, &context, depths)];
            let frame = bank
                .step_with_dry_contact_thermal_motion(
                    &mut heat,
                    &routes,
                    &[],
                    solve_controls(),
                    &[0.0; 6],
                    &forces,
                    &drives,
                    None,
                    None,
                    None,
                    gate,
                    other,
                )
                .unwrap();
            assert_eq!(frame.mechanics, expected);
            assert!(frame.routed_dry_heat_j > 0.0);
            assert_eq!(frame.routed_dry_heat_j, expected.ports[0].dissipation_j);
            let dt = bank.system().sample_period_s();
            let effective = depths.map(|depth| {
                let g = 0.2 * AREA / (0.5 * depth);
                let capacity = DENSITY * AREA * depth * CP;
                g * capacity / (capacity + dt * g)
            });
            let expected_left =
                frame.routed_dry_heat_j * effective[0] / (effective[0] + effective[1]);
            let left = heat.cells()[&id(0, 1)].mass_kg * heat.cells()[&id(0, 1)].specific_h_j_kg;
            assert!((left - expected_left).abs() < 2e-12);
            assert!((total_heat(&heat) - frame.routed_dry_heat_j).abs() < 2e-12);
            assert!(
                (frame.dry_contacts[0]
                    .cell_heat_j
                    .iter()
                    .map(|(_, q)| q)
                    .sum::<f64>()
                    - frame.routed_dry_heat_j)
                    .abs()
                    < 2e-12
            );
        }
    });
}

#[test]
fn moving_and_reversing_bow_visits_original_warm_cells_and_failed_source_rolls_back_geometry() {
    with_cx(|cx, gate| {
        let (card, caloric) = coating(false, false);
        let bundle = transport(&card);
        let context = point(300.0);
        let depths = [DEPTH; 2];
        let routes = [route(&card, caloric.curve(), &bundle, &context, depths)];
        let make = |temperature| {
            dry_bank(cx, gate, temperature)
                .with_bow_station_motion(motion_limits(), gate)
                .unwrap()
        };
        let mut bank = make(340.0);
        let mut heat = dry_network(caloric.curve(), [300.0; 2], depths);
        let forces = vec![0.0; bank.system().mode_count()];
        let drives = [drive(0.05, 0.1), drive(0.0, 0.0)];
        let mut previous = 0.2;
        let mut departed_heat = 0.0;
        for (step, station) in [0.202, 0.206, 0.202, 0.198].into_iter().enumerate() {
            let guess = guesses(&heat);
            let frame = bank
                .step_with_dry_contact_thermal_motion(
                    &mut heat,
                    &routes,
                    &[],
                    solve_controls(),
                    &guess,
                    &forces,
                    &drives,
                    Some(&[station, 0.3]),
                    None,
                    None,
                    gate,
                    other,
                )
                .unwrap();
            let interval = frame.dry_contacts[0].footprint.intervals_m[0];
            let centre = -0.097 + 0.25 * (station + previous);
            assert!((0.5 * (interval[0] + interval[1]) - centre).abs() < 1e-16);
            previous = station;
            if step == 0 {
                departed_heat = heat.cells()[&id(0, 1)].specific_h_j_kg;
                assert!(departed_heat > 0.0);
            } else if step <= 2 {
                assert!((heat.cells()[&id(0, 1)].specific_h_j_kg - departed_heat).abs() < 1e-7);
            } else {
                assert!(heat.cells()[&id(0, 1)].specific_h_j_kg > departed_heat);
            }
            assert!(frame.total_energy_residual_j.abs() < 2e-12);
        }
        let mut refused = make(300.0001);
        let mut expected = make(300.0001);
        let mut cold = dry_network(caloric.curve(), [300.0; 2], depths);
        let saved = (
            states(&refused),
            cold.cells().clone(),
            refused.bow_station_fractions().unwrap().to_vec(),
            refused.system().port(0).unwrap().left.clone(),
            refused.bow_interface(0).unwrap().clone(),
        );
        assert!(
            refused
                .step_with_dry_contact_thermal_motion(
                    &mut cold,
                    &routes,
                    &[],
                    solve_controls(),
                    &[0.0; 6],
                    &forces,
                    &drives,
                    Some(&[0.202, 0.3]),
                    None,
                    None,
                    gate,
                    other
                )
                .is_err()
        );
        assert_eq!(states(&refused), saved.0);
        assert_eq!(cold.cells(), &saved.1);
        assert_eq!(refused.bow_station_fractions().unwrap(), saved.2);
        assert_eq!(refused.system().port(0).unwrap().left, saved.3);
        assert_eq!(refused.bow_interface(0).unwrap(), &saved.4);
        assert_eq!(cold.time_s(), 0.0);
        assert_eq!(refused.system().samples_rendered(), 0);
        assert_eq!(
            refused
                .step_with_bow_stations_under_gate(&forces, &drives, &[0.202, 0.3], None, gate)
                .unwrap(),
            expected
                .step_with_bow_stations_under_gate(&forces, &drives, &[0.202, 0.3], None, gate)
                .unwrap()
        );

        // A coincident first strip edge can round BELOW local zero after a
        // much larger world-to-material translation. Preserve the footprint's
        // area when snapping that edge, even though its other edge is interior.
        let half = 0.0017;
        let mut zero_route = routes[0].clone();
        zero_route.string_coordinate_origin_m = half - 0.1;
        zero_route.heat.footprint.intervals_m[0] = [0.0, 2.0 * half];
        zero_route.heat.footprint.widths_m[0] = AREA / (2.0 * half);
        assert!(0.1 + zero_route.string_coordinate_origin_m - half < 0.0);
        let mut at_zero = make(340.0);
        let mut zero_heat = dry_network(caloric.curve(), [300.0; 2], depths);
        let frame = at_zero
            .step_with_dry_contact_thermal_motion(
                &mut zero_heat,
                &[zero_route.clone()],
                &[],
                solve_controls(),
                &[0.0; 6],
                &forces,
                &drives,
                None,
                None,
                None,
                gate,
                other,
            )
            .unwrap();
        assert_eq!(
            frame.dry_contacts[0].footprint.intervals_m[0],
            [0.0, 2.0 * half]
        );
        zero_route.string_coordinate_origin_m -= 1e-5;
        let guess = guesses(&zero_heat);
        assert!(
            at_zero
                .step_with_dry_contact_thermal_motion(
                    &mut zero_heat,
                    &[zero_route],
                    &[],
                    solve_controls(),
                    &guess,
                    &forces,
                    &drives,
                    None,
                    None,
                    None,
                    gate,
                    other
                )
                .is_err()
        );
        assert_eq!(at_zero.system().samples_rendered(), 1);
    });
}

#[test]
fn dry_thermal_performance_preserves_canonical_moving_pressure_and_callback_boundaries() {
    with_cx(|cx, gate| {
        let (card, caloric) = coating(false, false);
        let bundle = transport(&card);
        let context = point(300.0);
        let depths = [DEPTH; 2];
        let routes = [route(&card, caloric.curve(), &bundle, &context, depths)];
        let make = || {
            let schedule = GestureSchedule::try_new(
                701,
                vec![
                    GestureTrack {
                        id: "bow-0".into(),
                        target: GestureTarget::BowStroke { string: 0 },
                        initial: value(0.05, 0.1, 0.2),
                        events: vec![],
                    },
                    GestureTrack {
                        id: "bow-1".into(),
                        target: GestureTarget::BowStroke { string: 1 },
                        initial: value(0.0, 0.0, 0.3),
                        events: vec![],
                    },
                ],
            )
            .unwrap();
            dry_bank(cx, gate, 340.0)
                .into_scheduled(
                    schedule,
                    tracks(),
                    MultiBowPerformanceConfig {
                        samples: 48,
                        max_block: 16,
                        ..performance_config()
                    },
                )
                .unwrap()
                .with_bow_station_motion(
                    BowStationProgram {
                        knots: vec![
                            BowStationKnot {
                                sample: 24,
                                station_fractions: vec![0.204, 0.3],
                            },
                            BowStationKnot {
                                sample: 48,
                                station_fractions: vec![0.2, 0.3],
                            },
                        ],
                        limits: motion_limits(),
                        max_knots: 2,
                        max_coordinate_values: 4,
                    },
                    gate,
                )
                .unwrap()
        };
        let mut single = make();
        let mut chunks = make();
        let mut reference = make();
        let mut a = single
            .with_dry_contact_thermal_network(
                dry_network(caloric.curve(), [300.0; 2], depths),
                routes.to_vec(),
                vec![],
                solve_controls(),
                gate,
                other,
            )
            .unwrap();
        let mut b = chunks
            .with_dry_contact_thermal_network(
                dry_network(caloric.curve(), [300.0; 2], depths),
                routes.to_vec(),
                vec![],
                solve_controls(),
                gate,
                other,
            )
            .unwrap();
        let mut x = [0.0; 48];
        let mut y = x;
        let mut z = x;
        for block in x.chunks_mut(1) {
            a.block(block).unwrap();
        }
        for block in y.chunks_mut(7) {
            b.block(block).unwrap();
        }
        for block in z.chunks_mut(16) {
            reference.block(block).unwrap();
        }
        assert_eq!(x, y);
        assert_eq!(x, z);
        assert_eq!(a.thermal().cells(), b.thermal().cells());
        assert_eq!(states(a.performance().bank()), states(reference.bank()));
        assert_eq!(
            a.performance().bank().bow_station_fractions(),
            reference.bank().bow_station_fractions()
        );
        assert!(total_heat(a.thermal()) > 0.0);
    });
}

#[test]
fn accelerating_fiber_heat_uses_actual_material_midpoint_and_refusal_keeps_support_history() {
    with_cx(|cx, gate| {
        let (card, caloric) = coating(false, false);
        let bundle = transport(&card);
        let context = point(300.0);
        let depths = [DEPTH; 2];
        let make = || {
            base(cx, gate, 7800.0)
                .with_tensioned_fiber_ribbon(
                    crate::fiber_ribbon::fibers(1, 1, 1300.0, 3e9, 50.0, 40.0, 48000, gate),
                    0.4,
                    ribbon_limits(),
                    gate,
                )
                .unwrap()
                .with_fiber_support_motion(0.0, crate::fiber_ribbon::support_limits(), gate)
                .unwrap()
                .with_bow_interface(0, dry_source(340.0))
                .unwrap()
        };
        let mut bank = make();
        let mut reference = make();
        let mut heat = dry_network(caloric.curve(), [300.0; 2], depths);
        let mut binding = route(&card, caloric.curve(), &bundle, &context, depths);
        binding.receiver = DryReceiverHeatPosition::Fiber {
            coordinate_origin_m: -0.237,
        };
        let routes = [binding];
        let forces = vec![0.0; bank.system().mode_count()];
        let dt = bank.system().sample_period_s();
        for target in [0.02, -0.02, 0.01] {
            let initial = bank.fiber_support_kinematics().unwrap();
            // The support follows constant acceleration within the sample.
            // t=dt/2 is NOT the average of the two endpoint positions.
            let expected_centre = -0.237 + 0.24
                - initial.displacement_m
                - dt * (0.375 * initial.velocity_m_s + 0.125 * target);
            let drive = [drive(target, 0.001)];
            let guess = guesses(&heat);
            let frame = bank
                .step_with_dry_contact_thermal_motion(
                    &mut heat,
                    &routes,
                    &[],
                    solve_controls(),
                    &guess,
                    &forces,
                    &drive,
                    None,
                    None,
                    None,
                    gate,
                    other,
                )
                .unwrap();
            let interval = frame.dry_contacts[0].footprint.intervals_m[1];
            assert!((0.5 * (interval[0] + interval[1]) - expected_centre).abs() < 5e-17);
            let ordinary = reference.step(&forces, &drive).unwrap();
            assert_eq!(&frame.mechanics, ordinary);
            assert_eq!(
                bank.fiber_support_kinematics(),
                reference.fiber_support_kinematics()
            );
            assert!(frame.total_energy_residual_j.abs() < 2e-12);
        }
        assert!(total_heat(&heat) > 0.0);
        let saved = (
            states(&bank),
            heat.cells().clone(),
            bank.fiber_support_kinematics(),
            bank.system().last_frame().cloned(),
            bank.bow_interface(0).cloned(),
        );
        let guess = guesses(&heat);
        assert!(
            bank.step_with_dry_contact_thermal_motion(
                &mut heat,
                &routes,
                &[],
                solve_controls(),
                &guess,
                &forces,
                &[drive(0.02, 0.001)],
                None,
                None,
                None,
                gate,
                |_| -> Result<StepCoefficients<()>, &'static str> {
                    Err("late dry heat source refusal")
                }
            )
            .is_err()
        );
        assert_eq!(states(&bank), saved.0);
        assert_eq!(heat.cells(), &saved.1);
        assert_eq!(bank.fiber_support_kinematics(), saved.2);
        assert_eq!(bank.system().last_frame(), saved.3.as_ref());
        assert_eq!(bank.bow_interface(0), saved.4.as_ref());
        assert_eq!(bank.system().samples_rendered(), 3);
        assert_eq!(
            bank.step(&forces, &[drive(0.02, 0.001)]).unwrap(),
            reference.step(&forces, &[drive(0.02, 0.001)]).unwrap()
        );
    });
}
