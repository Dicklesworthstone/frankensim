//! Gas loading, molecular heat and moving walls on one original thermal tick.

use super::*;
use fs_lbm::material::free_surface::enthalpy::{
    FreeSurfaceWallGasPressure, FreeSurfaceWallSources,
};
use fs_lbm::material::free_surface::{MaterialGasPressureLink3, MaterialPlanarRigidWallMotion3};
use std::cell::Cell;

fn with_pocket<R>(
    reference: bool,
    exercise: impl FnOnce(&mut MaterialFreeSurfaceEnthalpy3<'_, '_>, &CancelGate) -> R,
) -> R {
    let card = card(0.02, 1000.0);
    let law = LiquidTemperatureLaw::try_new(
        &card,
        &point(400.0),
        "T",
        LiquidViscosityBasis::Dynamic,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let source = resolve_liquid_surface_state(
        &card,
        &point(400.0),
        LiquidViscosityBasis::Dynamic,
        LiquidCapillarity::Neglected,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let chart = caloric(&card);
    let gate = CancelGate::new_clock_free();
    let mut cells = initial(true);
    // An actual interior gas pocket leaves both complete z wall planes wet.
    cells.retain(|cell| !cell.coordinates.iter().all(|c| (3..5).contains(c)));
    for cell in &mut cells {
        cell.velocity_m_s = [0.0; 3];
    }
    with_cx(&gate, |cx| {
        let mut scaling = plan_free_surface_from_material(&source, 0.001, 0.001, 0.1).unwrap();
        if reference {
            scaling = scaling.with_pressure_reference(100000.0).unwrap();
        }
        let flow = scaling
            .free_surface3_from_cells(domain(), &cells, cx)
            .unwrap();
        let mut liquid = flow
            .with_enthalpy(
                &law,
                &chart,
                FreeSurfaceEnthalpyControls::default(),
                |_| 200000.0,
                &gate,
            )
            .unwrap();
        exercise(&mut liquid, &gate)
    })
}

fn motions() -> [MaterialPlanarRigidWallMotion3; 2] {
    [PlanarWall3::ZMin, PlanarWall3::ZMax].map(|wall| {
        let upper = wall == PlanarWall3::ZMax;
        MaterialPlanarRigidWallMotion3 {
            wall,
            origin_m: [0.004, 0.004, if upper { 0.008 } else { 0.0 }],
            velocity_at_origin_m_s: [if upper { -0.015 } else { 0.025 }, 0.004, 0.0],
            angular_speed_rad_s: if upper { -0.8 } else { 1.2 },
        }
    })
}

fn heat() -> [FreeSurfaceCellHeat; 1] {
    [FreeSurfaceCellHeat {
        coordinates: [3, 3, 0],
        energy_j: 0.001,
    }]
}

fn cell_sources<'a>(
    heat: &'a [FreeSurfaceCellHeat],
    pressure_pa: &'a dyn Fn([f64; 3]) -> f64,
) -> FreeSurfaceWallSources<'a> {
    FreeSurfaceWallSources {
        heat,
        molecular_heating: false,
        gas_pressure: Some(FreeSurfaceWallGasPressure::Cell {
            maximum_density_excursion: 0.05,
            pressure_pa,
        }),
    }
}

fn same_liquid(
    actual: &MaterialFreeSurfaceEnthalpy3<'_, '_>,
    expected: &MaterialFreeSurfaceEnthalpy3<'_, '_>,
) {
    assert_eq!(
        actual.flow().lattice().grid().state_bits(),
        expected.flow().lattice().grid().state_bits()
    );
    assert_eq!(
        actual.total_enthalpy_j().unwrap(),
        expected.total_enthalpy_j().unwrap()
    );
    assert_eq!(actual.accepted_steps(), expected.accepted_steps());
    assert_eq!(actual.physical_time_s(), expected.physical_time_s());
    for p in [[3, 3, 0], [3, 3, 2], [4, 3, 5], [3, 3, 7]] {
        assert_eq!(
            actual.specific_enthalpy_at(p),
            expected.specific_enthalpy_at(p)
        );
        assert_eq!(actual.material_at(p), expected.material_at(p));
    }
}

#[test]
fn g1_reference_pressure_retains_all_three_original_thermal_wall_paths() {
    with_pocket(true, |liquid, gate| {
        let calls = Cell::new(0);
        let reference = |_| {
            calls.set(calls.get() + 1);
            100000.0
        };
        let heat = heat();
        let sources = cell_sources(&heat, &reference);
        let velocities = field();
        let mut original_field = liquid.clone();
        let mut loaded_field = liquid.clone();
        let original = original_field
            .step_with_planar_wall_velocity_field_si(motion(&velocities), controls(), &heat, gate)
            .unwrap();
        let loaded = loaded_field
            .step_with_planar_wall_velocity_field_and_sources_si(
                motion(&velocities),
                controls(),
                sources,
                gate,
            )
            .unwrap();
        assert_eq!(loaded.flow, original);
        assert!(loaded.applied.is_empty());
        same_liquid(&loaded_field, &original_field);

        let motion = motions();
        let (original_single, original) = liquid
            .planar_rigid_wall_motion_candidate_si(motion[0], 1e-11, &heat, gate)
            .unwrap();
        let (loaded_single, loaded) = liquid
            .planar_rigid_wall_motion_candidate_with_sources_si(motion[0], 1e-11, sources, gate)
            .unwrap();
        assert_eq!(loaded.flow, original);
        assert!(loaded.applied.is_empty());
        same_liquid(&loaded_single, &original_single);

        let (original_pair, original) = liquid
            .planar_rigid_wall_pair_candidate_si(motion, 1e-11, &heat, gate)
            .unwrap();
        let (loaded_pair, loaded) = liquid
            .planar_rigid_wall_pair_candidate_with_sources_si(motion, 1e-11, sources, gate)
            .unwrap();
        assert_eq!(loaded.flow, original);
        assert!(loaded.applied.is_empty());
        same_liquid(&loaded_pair, &original_pair);
        assert!(
            calls.get() > 0,
            "real gas links must sample the pressure field"
        );
        assert_eq!(liquid.accepted_steps(), 0);
    });
}

#[test]
fn g1_cell_and_link_loads_share_wall_heat_and_one_native_mass_trace() {
    with_pocket(true, |liquid, gate| {
        let before_h = liquid.total_enthalpy_j().unwrap();
        let before_mass = liquid.flow().ledger_mass_kg().unwrap();
        let heat = heat();
        let cell_calls = Cell::new(0);
        let link_calls = Cell::new(0);
        let load = |p: [f64; 3]| 100000.0 + if p[0] < 0.004 { 0.5 } else { -0.5 };
        let cell_pressure = |p| {
            cell_calls.set(cell_calls.get() + 1);
            load(p)
        };
        let link_pressure = |link: MaterialGasPressureLink3| {
            link_calls.set(link_calls.get() + 1);
            assert!(link.link.incoming_direction > 0);
            load(link.liquid_center_m)
        };
        let mut cell = liquid.clone();
        let mut link = liquid.clone();
        let mut unloaded = liquid.clone();
        for step in 0..3 {
            let heat = if step == 0 { heat.as_slice() } else { &[] };
            let (next_cell, actual) = cell
                .planar_rigid_wall_pair_candidate_with_sources_si(
                    motions(),
                    1e-11,
                    cell_sources(heat, &cell_pressure),
                    gate,
                )
                .unwrap();
            let (next_link, expected) = link
                .planar_rigid_wall_pair_candidate_with_sources_si(
                    motions(),
                    1e-11,
                    FreeSurfaceWallSources {
                        heat,
                        molecular_heating: false,
                        gas_pressure: Some(FreeSurfaceWallGasPressure::Link {
                            maximum_density_excursion: 0.05,
                            pressure_pa: &link_pressure,
                        }),
                    },
                    gate,
                )
                .unwrap();
            unloaded = unloaded
                .planar_rigid_wall_pair_candidate_si(motions(), 1e-11, heat, gate)
                .unwrap()
                .0;
            assert_eq!(actual, expected);
            same_liquid(&next_cell, &next_link);
            assert!(actual.flow.flow.mass_routes > 0);
            assert_eq!(actual.flow.flow.viscous_heat_j, 0.0);
            near(actual.flow.flow.mass_after_kg, before_mass, 2e-12);
            near(actual.flow.flow.enthalpy_after_j, before_h + 0.001, 2e-11);
            assert!(actual.flow.walls.fluid_work_j != 0.0);
            assert_eq!(next_cell.accepted_steps(), step + 1);
            cell = next_cell;
            link = next_link;
        }
        assert!(cell_calls.get() > 0);
        assert!(link_calls.get() > cell_calls.get());
        assert_ne!(
            cell.flow().lattice().grid().state_bits(),
            unloaded.flow().lattice().grid().state_bits()
        );
        near(
            cell.total_enthalpy_j().unwrap(),
            unloaded.total_enthalpy_j().unwrap(),
            2e-11,
        );
        assert_eq!(liquid.accepted_steps(), 0);
        assert_eq!(liquid.total_enthalpy_j().unwrap(), before_h);
    });
}

#[test]
fn g3_pressure_and_moving_pair_apply_original_molecular_heat_once() {
    with_pocket(true, |liquid, gate| {
        let load = |p: [f64; 3]| 100000.0 + if p[0] < 0.004 { 0.5 } else { -0.5 };
        *liquid = liquid
            .planar_rigid_wall_pair_candidate_with_sources_si(
                motions(),
                1e-11,
                cell_sources(&[], &load),
                gate,
            )
            .unwrap()
            .0;
        let mut source_probe = liquid.clone();
        let source = source_probe.step_with_viscous_heating(&[], gate).unwrap();
        let applied = source.applied;
        let molecular: f64 = applied.iter().map(|sample| sample.heat_j).sum();
        assert!(molecular > 0.0);
        let equivalent: Vec<_> = applied
            .iter()
            .map(|sample| FreeSurfaceCellHeat {
                coordinates: sample.coordinates,
                energy_j: sample.heat_j,
            })
            .collect();
        let before_h = liquid.total_enthalpy_j().unwrap();
        let before_steps = liquid.accepted_steps();
        let (expected, reference) = liquid
            .planar_rigid_wall_pair_candidate_with_sources_si(
                motions(),
                1e-11,
                cell_sources(&equivalent, &load),
                gate,
            )
            .unwrap();
        let sources = FreeSurfaceWallSources {
            molecular_heating: true,
            ..cell_sources(&[], &load)
        };
        let (actual, report) = liquid
            .planar_rigid_wall_pair_candidate_with_sources_si(motions(), 1e-11, sources, gate)
            .unwrap();
        let (replay, replay_report) = liquid
            .planar_rigid_wall_pair_candidate_with_sources_si(motions(), 1e-11, sources, gate)
            .unwrap();
        assert_eq!(report, replay_report);
        assert_eq!(report.applied, applied);
        assert_eq!(report.flow.walls, reference.flow.walls);
        same_liquid(&actual, &expected);
        same_liquid(&actual, &replay);
        assert_eq!(report.flow.flow.external_heat_j, 0.0);
        assert_eq!(report.flow.flow.viscous_heat_j, molecular);
        near(
            actual.total_enthalpy_j().unwrap(),
            before_h + molecular,
            2e-11,
        );
        assert_eq!(actual.accepted_steps(), before_steps + 1);
        assert_eq!(liquid.accepted_steps(), before_steps);
        assert_eq!(liquid.total_enthalpy_j().unwrap(), before_h);
    });
}

#[test]
fn g4_late_pressure_source_and_cancel_refusals_leave_wall_heat_and_clock_untouched() {
    with_pocket(true, |liquid, gate| {
        let original = liquid.clone();
        let velocities = field();
        let heat = heat();
        for cancel in [false, true] {
            let attempt_gate = CancelGate::new_clock_free();
            let calls = Cell::new(0);
            let pressure = |_| {
                calls.set(calls.get() + 1);
                if calls.get() == 5 {
                    if cancel {
                        attempt_gate.request();
                    } else {
                        return 0.0;
                    }
                }
                100000.5
            };
            let result = liquid.step_with_planar_wall_velocity_field_and_sources_si(
                motion(&velocities),
                controls(),
                FreeSurfaceWallSources {
                    heat: &heat,
                    molecular_heating: true,
                    gas_pressure: Some(FreeSurfaceWallGasPressure::Link {
                        maximum_density_excursion: 0.05,
                        pressure_pa: &pressure,
                    }),
                },
                &attempt_gate,
            );
            assert!(result.is_err());
            assert!(calls.get() >= 5);
            same_liquid(liquid, &original);
        }
        let invalid_heat = [FreeSurfaceCellHeat {
            coordinates: [3, 3, 0],
            energy_j: 6.0,
        }];
        assert!(
            liquid
                .planar_rigid_wall_pair_candidate_with_sources_si(
                    motions(),
                    1e-11,
                    cell_sources(&invalid_heat, &|_| 100000.5),
                    gate,
                )
                .is_err()
        );
        same_liquid(liquid, &original);
        let mut expected = original.clone();
        let reference = expected
            .step_with_planar_wall_velocity_field_and_sources_si(
                motion(&velocities),
                controls(),
                cell_sources(&heat, &|_| 100000.5),
                gate,
            )
            .unwrap();
        let retried = liquid
            .step_with_planar_wall_velocity_field_and_sources_si(
                motion(&velocities),
                controls(),
                cell_sources(&heat, &|_| 100000.5),
                gate,
            )
            .unwrap();
        assert_eq!(retried, reference);
        same_liquid(liquid, &expected);
    });
    with_pocket(false, |liquid, gate| {
        let original = liquid.clone();
        assert!(
            liquid
                .planar_rigid_wall_motion_candidate_with_sources_si(
                    motions()[0],
                    1e-11,
                    cell_sources(&[], &|_| 100000.0),
                    gate,
                )
                .is_err()
        );
        same_liquid(liquid, &original);
    });
}
