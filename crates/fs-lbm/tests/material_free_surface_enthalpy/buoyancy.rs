//! Source density reaches the retained native force and moving-liquid state.

use super::*;
use fs_lbm::material::free_surface::enthalpy::MaterialFreeSurfaceEnthalpy3;
use fs_lbm::material::transport::boussinesq::BoussinesqConfig;

#[derive(Clone, Copy)]
struct Case {
    full: bool,
    buoyant: bool,
    gravity: [f64; 3],
    shear: f64,
    presteps: usize,
    density_limit: f64,
    controls: FreeSurfaceEnthalpyControls,
}

impl Default for Case {
    fn default() -> Self {
        Self {
            full: true,
            buoyant: true,
            gravity: [0.0, 0.0, -1.0],
            shear: 0.0,
            presteps: 0,
            density_limit: 0.03,
            controls: FreeSurfaceEnthalpyControls::default(),
        }
    }
}

fn config(case: Case) -> BoussinesqConfig<3> {
    BoussinesqConfig {
        gravity_m_s2: case.gravity,
        max_relative_density_departure: case.density_limit,
        max_relative_specific_heat_departure: 0.0,
    }
}

fn with_liquid<R>(
    card: &MaterialCard,
    case: Case,
    h: impl Fn([usize; 3]) -> f64,
    exercise: impl FnOnce(&mut MaterialFreeSurfaceEnthalpy3<'_, '_>, &CancelGate) -> R,
) -> R {
    let law = LiquidTemperatureLaw::try_new(
        card,
        &point(300.0),
        "T",
        LiquidViscosityBasis::Dynamic,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let source = resolve_liquid_surface_state(
        card,
        &point(300.0),
        LiquidViscosityBasis::Dynamic,
        LiquidCapillarity::Neglected,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let chart = caloric(card);
    let gate = CancelGate::new_clock_free();
    let mut cells = initial(case.full);
    for cell in &mut cells {
        cell.velocity_m_s = if case.full {
            [case.shear * (cell.coordinates[2] as f64 - 3.5), 0.0, 0.0]
        } else {
            [0.02, 0.0, 0.0]
        };
    }
    with_cx(&gate, |cx| {
        let mut flow = plan_free_surface_from_material(&source, 0.001, 0.001, 0.1)
            .unwrap()
            .free_surface3_from_cells(
                MaterialFreeSurfaceDomain3 {
                    gravity_m_s2: case.gravity,
                    ..domain()
                },
                &cells,
                cx,
            )
            .unwrap();
        for _ in 0..case.presteps {
            flow.step_with_gate(&gate).unwrap();
        }
        let before = flow.lattice().grid().state_bits();
        let velocity = flow.sample_at_si(3, 3, 3).unwrap().unwrap().velocity_m_s;
        let mut thermal = if case.buoyant {
            flow.with_boussinesq_enthalpy(&law, &chart, case.controls, config(case), h, &gate)
        } else {
            flow.with_enthalpy(&law, &chart, case.controls, h, &gate)
        }
        .unwrap();
        if !case.buoyant {
            assert_eq!(thermal.flow().lattice().grid().state_bits(), before);
        }
        let after = thermal
            .flow()
            .sample_at_si(3, 3, 3)
            .unwrap()
            .unwrap()
            .velocity_m_s;
        for axis in 0..3 {
            assert!((after[axis] - velocity[axis]).abs() < 2e-15);
        }
        exercise(&mut thermal, &gate)
    })
}

#[test]
fn g3_zero_density_departure_preserves_the_original_gravity_path_exactly() {
    let card = card(0.1, 1000.0);
    for full in [true, false] {
        let run = |buoyant| {
            with_liquid(
                &card,
                Case {
                    full,
                    buoyant,
                    gravity: [0.0, 0.0, -0.25],
                    shear: 0.003,
                    ..Case::default()
                },
                |p| if p[0] < 4 { 100000.0 } else { 300000.0 },
                |thermal, gate| {
                    let report = thermal.step(gate).unwrap();
                    (
                        thermal.flow().lattice().grid().state_bits(),
                        thermal.flow().sample_at_si(3, 3, 3).unwrap().unwrap(),
                        thermal.total_enthalpy_j().unwrap(),
                        report.mass_after_kg,
                    )
                },
            )
        };
        assert_eq!(run(false), run(true));
    }
}

#[test]
fn g1_density_source_and_changed_heat_drive_the_native_held_force_impulse() {
    let run = |hot_density| {
        let card = card(0.1, hot_density);
        with_liquid(
            &card,
            Case {
                presteps: 1,
                ..Case::default()
            },
            |_| 200000.0,
            |thermal, gate| {
                let mass = thermal.flow().ledger_mass_kg().unwrap();
                let energy = thermal.total_enthalpy_j().unwrap();
                thermal.step(gate).unwrap();
                let first = thermal
                    .flow()
                    .sample_at_si(3, 3, 3)
                    .unwrap()
                    .unwrap()
                    .velocity_m_s[2];
                // Attach to an already-evolved flow. The original reference
                // step contributes -0.001 m/s; the new held interval contributes
                // g*dt*rho(400 K)/rho0 before a wall signal reaches this sample.
                let rho_400 = (1000.0 + hot_density) / 2.0;
                near(first, -0.001 - 0.001 * rho_400 / 1000.0, 2e-11);
                near(
                    thermal.flow().lattice().force_density_at(3, 3, 3).unwrap()[2],
                    -0.001 * rho_400 / 1000.0,
                    2e-11,
                );
                near(thermal.total_enthalpy_j().unwrap(), energy, 2e-11);

                // Add the same specific heat to every actual native mass. Force
                // changes at the thermal split but does not add a second half kick.
                let mut heat = Vec::new();
                for z in 0..8 {
                    for y in 0..8 {
                        for x in 0..8 {
                            let cell = thermal.flow().sample_at_si(x, y, z).unwrap().unwrap();
                            heat.push(FreeSurfaceCellHeat {
                                coordinates: [x, y, z],
                                energy_j: cell.density_kg_m3 * 1e-9 * 100000.0,
                            });
                        }
                    }
                }
                thermal.step_with_heat(&heat, gate).unwrap();
                let second = thermal
                    .flow()
                    .sample_at_si(3, 3, 3)
                    .unwrap()
                    .unwrap()
                    .velocity_m_s[2];
                let rho_450 = 1000.0 + 0.75 * (hot_density - 1000.0);
                near(second - first, -0.001 * rho_450 / 1000.0, 2e-11);
                near(thermal.flow().ledger_mass_kg().unwrap(), mass, 2e-11);
                near(
                    thermal.total_enthalpy_j().unwrap() - energy,
                    heat.iter().map(|q| q.energy_j).sum(),
                    2e-11,
                );
                (first, second)
            },
        )
    };
    let fixed = run(1000.0);
    let expanding = run(980.0);
    assert!(expanding.0 > fixed.0);
    assert!(expanding.1 > fixed.1);
}

#[test]
fn g3_density_only_changes_do_not_replace_reference_inertia_or_viscosity() {
    let run = |hot_density| {
        let card = card(0.1, hot_density);
        with_liquid(
            &card,
            Case {
                gravity: [0.0; 3],
                shear: 0.003,
                ..Case::default()
            },
            |_| 200000.0,
            |thermal, gate| {
                for _ in 0..2 {
                    thermal.step(gate).unwrap();
                }
                (
                    thermal.flow().lattice().grid().state_bits(),
                    thermal.total_enthalpy_j().unwrap(),
                )
            },
        )
    };
    // No gravity means no density-force difference. Dynamic mu/rho0 must stay
    // identical even though the retained material's own mu/rho(T) is different.
    assert_eq!(run(1000.0), run(980.0));
}

#[test]
fn g1_partial_fill_nonzero_gravity_has_a_continuous_density_reference_limit() {
    let run = |hot_density| {
        let card = card(0.1, hot_density);
        with_liquid(
            &card,
            Case {
                full: false,
                gravity: [0.0, 0.0, -0.25],
                ..Case::default()
            },
            |_| 200000.0,
            |thermal, gate| {
                thermal.step(gate).unwrap();
                let mut velocities = Vec::new();
                for z in 1..=6 {
                    for y in 1..=6 {
                        for x in 1..=3 {
                            velocities.push(
                                thermal
                                    .flow()
                                    .sample_at_si(x, y, z)
                                    .unwrap()
                                    .unwrap()
                                    .velocity_m_s,
                            );
                        }
                    }
                }
                // Any newly created interface also needs a retained native
                // source frame after h extension and endpoint admission.
                for z in 0..8 {
                    for y in 0..8 {
                        for x in 0..8 {
                            if let Some(sample) = thermal.flow().sample_at_si(x, y, z).unwrap() {
                                assert!(sample.velocity_m_s.iter().all(|v| v.is_finite()));
                                assert!(
                                    thermal.flow().lattice().force_density_at(x, y, z).is_some()
                                );
                            }
                        }
                    }
                }
                velocities
            },
        )
    };
    let reference = run(1000.0);
    let error = |actual: Vec<[f64; 3]>| {
        actual
            .iter()
            .flatten()
            .zip(reference.iter().flatten())
            .map(|(a, b)| (a - b).abs())
            .fold(0.0_f64, f64::max)
    };
    let coarse = error(run(999.999));
    let fine = error(run(999.9995));
    assert!(
        coarse > 1e-13 && coarse < 1e-8,
        "no O(reference gravity) boundary jump: {coarse}"
    );
    assert!(
        fine > 0.35 * coarse && fine < 0.75 * coarse,
        "density refinement: {coarse} -> {fine}"
    );
}

#[test]
fn g4_source_bounds_cancellation_and_late_transport_refusal_preserve_the_force_epoch() {
    let card = card(0.1, 980.0);
    let case = Case {
        density_limit: 0.015,
        ..Case::default()
    };
    let clean = with_liquid(
        &card,
        case,
        |_| 200000.0,
        |thermal, gate| {
            thermal.step(gate).unwrap();
            (
                thermal.flow().lattice().grid().state_bits(),
                thermal.total_enthalpy_j().unwrap(),
                thermal.flow().lattice().force_density_at(3, 3, 3).unwrap(),
            )
        },
    );
    with_liquid(
        &card,
        case,
        |_| 200000.0,
        |thermal, gate| {
            let bits = thermal.flow().lattice().grid().state_bits();
            let energy = thermal.total_enthalpy_j().unwrap();
            let force = thermal.flow().lattice().force_density_at(3, 3, 3).unwrap();
            let state = thermal.material_at([3, 3, 3]).unwrap().clone();
            for energy_j in [0.15, 0.3] {
                // First reaches supported T=475 but exceeds density projection;
                // second leaves the original source's temperature domain entirely.
                assert!(
                    thermal
                        .step_with_heat(
                            &[FreeSurfaceCellHeat {
                                coordinates: [3, 3, 3],
                                energy_j,
                            }],
                            gate
                        )
                        .is_err()
                );
                assert_eq!(thermal.flow().lattice().grid().state_bits(), bits);
                assert_eq!(thermal.total_enthalpy_j().unwrap(), energy);
                assert_eq!(
                    thermal.flow().lattice().force_density_at(3, 3, 3).unwrap(),
                    force
                );
                assert_eq!(thermal.material_at([3, 3, 3]).unwrap(), &state);
            }
            let cancelled = CancelGate::new_clock_free();
            cancelled.request();
            assert!(thermal.step(&cancelled).is_err());
            assert_eq!(thermal.flow().lattice().grid().state_bits(), bits);
            thermal.step(gate).unwrap();
            assert_eq!(thermal.flow().lattice().grid().state_bits(), clean.0);
            assert_eq!(thermal.total_enthalpy_j().unwrap(), clean.1);
            assert_eq!(
                thermal.flow().lattice().force_density_at(3, 3, 3).unwrap(),
                clean.2
            );
        },
    );
    with_liquid(
        &card,
        Case {
            full: false,
            controls: FreeSurfaceEnthalpyControls {
                max_mass_routes: 1,
                ..case.controls
            },
            ..case
        },
        |_| 200000.0,
        |thermal, gate| {
            let bits = thermal.flow().lattice().grid().state_bits();
            let force = thermal.flow().lattice().force_density_at(3, 3, 3).unwrap();
            let energy = thermal.total_enthalpy_j().unwrap();
            assert!(matches!(
                thermal.step(gate),
                Err(FreeSurfaceEnthalpyError::Flow(_))
            ));
            assert_eq!(thermal.flow().lattice().grid().state_bits(), bits);
            assert_eq!(
                thermal.flow().lattice().force_density_at(3, 3, 3).unwrap(),
                force
            );
            assert_eq!(thermal.total_enthalpy_j().unwrap(), energy);
            assert_eq!(thermal.flow().lattice().conversions().to_fluid, 0);
        },
    );
}

#[cfg(feature = "enthalpy-coupling")]
#[test]
fn g3_conduction_and_buoyancy_use_the_same_source_temperature() {
    use fs_lbm::material::free_surface::enthalpy::conduction::{
        FreeSurfaceConductionControls, FreeSurfaceConductionFace,
    };
    let card = card(0.1, 980.0);
    with_liquid(
        &card,
        Case::default(),
        |p| if p == [3, 3, 3] { 100000.0 } else { 300000.0 },
        |thermal, gate| {
            let energy = thermal.total_enthalpy_j().unwrap();
            let report = thermal
                .step_with_conduction(
                    &[FreeSurfaceConductionFace {
                        cells: [[3, 3, 3], [4, 3, 3]],
                        area_m2: 1e-6,
                        distances_m: [0.0005; 2],
                    }],
                    &[],
                    FreeSurfaceConductionControls {
                        maximum_faces: 1,
                        maximum_substeps: 10,
                        maximum_dt_s: 0.00025,
                        energy_tolerance_j: 1e-9,
                    },
                    gate,
                )
                .unwrap();
            assert_eq!(report.thermal.substeps, 4);
            assert!(report.thermal_endpoint[0].materials[0].density_kg_m3() < 995.0);
            let state = thermal.material_at([3, 3, 3]).unwrap();
            near(
                thermal.flow().lattice().force_density_at(3, 3, 3).unwrap()[2],
                -0.001 * state.density_kg_m3() / 1000.0,
                2e-12,
            );
            near(thermal.total_enthalpy_j().unwrap(), energy, 2e-11);
        },
    );
}
