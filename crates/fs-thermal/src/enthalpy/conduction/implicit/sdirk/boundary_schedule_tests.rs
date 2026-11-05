use super::source_schedule::adaptive_controls;
use super::*;
use crate::enthalpy::conduction::implicit::adaptive::AdaptiveImplicitConductionError;
use crate::enthalpy::conduction::implicit::boundary::{GrayAmbientBoundary, GrayAmbientControls};

fn baths(ambient_k: f64, reservoir_k: Option<f64>, htc_w_m2_k: f64) -> Boundary {
    Boundary {
        ambient_k,
        htc_w_m2_k,
        ..boundary(reservoir_k)
    }
}

fn change(time_s: f64, boundary: Boundary) -> TimedEvent {
    TimedEvent {
        time_s,
        event: Event::Boundary(id(0), boundary),
    }
}

#[test]
fn scheduled_bath_temperatures_refine_on_fixed_and_adaptive_chains_and_graphs() {
    let law = Law::constant(2.0, 1.0);
    let (mut thermal, observer) = network(
        &law,
        &[400.0],
        &[1.0],
        &[],
        Some(300.0),
        vec![
            change(0.35, baths(240.0, Some(250.0), 1.0 / 6.0)),
            change(0.8, baths(320.0, Some(350.0), 1.0 / 6.0)),
        ],
    );
    thermal.cells.get_mut(&id(0)).unwrap().boundary = baths(260.0, Some(300.0), 1.0 / 6.0);
    let before = thermal.cells().clone();
    let direct = chain(&thermal, &observer, &law);
    let graph = ImplicitConductionGraph::new(
        &thermal,
        &observer,
        &[(MaterialId(0), &law)],
        ImplicitConductionLinearControls::default(),
        || Ok(()),
    )
    .unwrap();
    let mut workspace = ImplicitConductionWorkspace::new(1).unwrap();
    // C=2 J/K; reservoir G=1 W/K and the remaining 3 m2 give ambient G=0.5.
    // Compose exact exponential responses to each physical bath history.
    let mut exact = 400.0;
    for (span, ambient, reservoir) in [
        (0.35_f64, 260.0, 300.0),
        (0.45, 240.0, 250.0),
        (0.4, 320.0, 350.0),
    ] {
        let equilibrium = (0.5 * ambient + reservoir) / 1.5;
        exact = equilibrium + (exact - equilibrium) * (-1.5 * span / 2.0).exp();
    }
    for sparse in [false, true] {
        let mut errors = Vec::new();
        for dt in [0.2, 0.1, 0.05] {
            let (next, report) = if sparse {
                graph.advance_sdirk2_to(&thermal, 1.2, controls(dt), &mut workspace, || Ok(()))
            } else {
                direct.advance_sdirk2_to(&thermal, 1.2, controls(dt), &mut workspace, || Ok(()))
            }
            .unwrap();
            let t = next.temperature_k(id(0)).unwrap();
            errors.push((t - exact).abs());
            near(report.energy.external_heat_j, 2.0 * (t - 400.0), 1e-8);
            assert_eq!(report.energy.events, 2);
        }
        for pair in errors.windows(2) {
            assert!(
                (3.5..4.5).contains(&(pair[0] / pair[1])),
                "bath-event errors: {errors:?}"
            );
        }
        let (next, report) = if sparse {
            graph.advance_sdirk2_adaptive_to(
                &thermal,
                1.2,
                adaptive_controls(0.2),
                &mut workspace,
                || Ok(()),
            )
        } else {
            direct.advance_sdirk2_adaptive_to(
                &thermal,
                1.2,
                adaptive_controls(0.2),
                &mut workspace,
                || Ok(()),
            )
        }
        .unwrap();
        near(next.temperature_k(id(0)).unwrap(), exact, 4e-4);
        assert_eq!(report.solve.energy.events, 2);
        assert!(report.maximum_accepted_normalized_error <= 1.0);
    }
    assert_eq!(thermal.cells(), &before);
    assert_eq!(thermal.time_s(), 0.0);
}

#[test]
fn native_robin_switches_include_insulation_and_target_time_continuation() {
    let law = Law::constant(3.0, 1.0);
    let (mut thermal, observer) = network(
        &law,
        &[400.0],
        &[2.0],
        &[],
        None,
        vec![
            change(0.35, baths(290.0, None, 0.0)),
            change(0.8, baths(500.0, None, 1.0)),
            change(1.2, baths(200.0, None, 0.25)),
        ],
    );
    thermal.cells.get_mut(&id(0)).unwrap().boundary = baths(350.0, None, 0.5);
    let solver = chain(&thermal, &observer, &law);
    let mut workspace = ImplicitConductionWorkspace::new(1).unwrap();
    // C=6 J/K, exposed area=4 m2. Each nonzero Robin interval uses its own G.
    let cooled = 350.0 + 50.0 * amplification(2.0 * 0.35 / 6.0);
    let heated = 500.0 + (cooled - 500.0) * amplification(4.0 * 0.4 / 6.0);
    let (next, report) = solver
        .advance_sdirk2_to(&thermal, 1.2, controls(1.0), &mut workspace, || Ok(()))
        .unwrap();
    near(next.temperature_k(id(0)).unwrap(), heated, 1e-11);
    near(report.energy.external_heat_j, 6.0 * (heated - 400.0), 1e-9);
    assert_eq!(report.energy.events, 3);
    assert_eq!(next.cells()[&id(0)].boundary, baths(200.0, None, 0.25));

    let (at_switch, _) = solver
        .advance_sdirk2_to(&thermal, 0.35, controls(1.0), &mut workspace, || Ok(()))
        .unwrap();
    let (insulated, insulation) = solver
        .advance_sdirk2_to(&at_switch, 0.8, controls(1.0), &mut workspace, || Ok(()))
        .unwrap();
    near(insulated.temperature_k(id(0)).unwrap(), cooled, 1e-11);
    assert_eq!(insulation.energy.external_heat_j, 0.0);
    assert_eq!(insulation.entropy.reservoir_entropy_input_j_k, 0.0);
    assert_eq!(insulation.entropy.conduction_generation_j_k, 0.0);
    // The event at 1.2 s changes only the returned boundary. Cooling resumes
    // under its new G=1 W/K and 200 K bath in the following physical interval.
    let (continued, continuation) = solver
        .advance_sdirk2_to(&next, 1.25, controls(1.0), &mut workspace, || Ok(()))
        .unwrap();
    let expected = 200.0 + (heated - 200.0) * amplification(0.05 / 6.0);
    near(continued.temperature_k(id(0)).unwrap(), expected, 1e-11);
    near(
        continuation.energy.external_heat_j,
        6.0 * (expected - heated),
        1e-9,
    );
    assert_eq!(continuation.energy.events, 0);
    assert_eq!(thermal.time_s(), 0.0);
}

#[test]
fn bath_jumps_preserve_native_stage_heat_and_reservoir_entropy() {
    let law = Law::constant(3.0, 1.0);
    let (mut thermal, observer) = network(
        &law,
        &[600.0],
        &[2.0],
        &[],
        Some(350.0),
        vec![change(0.7, baths(500.0, Some(650.0), 1.0 / 6.0))],
    );
    thermal.cells.get_mut(&id(0)).unwrap().boundary = baths(250.0, Some(350.0), 1.0 / 6.0);
    let solver = chain(&thermal, &observer, &law);
    let mut workspace = ImplicitConductionWorkspace::new(1).unwrap();
    let (next, report) = solver
        .advance_sdirk2_to(&thermal, 1.6, controls(2.0), &mut workspace, || Ok(()))
        .unwrap();
    let mut temperature = 600.0;
    let mut heat = 0.0;
    let mut entropy_input = 0.0;
    let mut generation = 0.0;
    for (dt, ambient, reservoir) in [(0.7, 250.0, 350.0), (0.9, 500.0, 650.0)] {
        let equilibrium = (0.5 * ambient + reservoir) / 1.5;
        let first = equilibrium + (temperature - equilibrium) / (1.0 + GAMMA * dt * 1.5 / 6.0);
        let second = equilibrium + (temperature - equilibrium) * amplification(dt * 1.5 / 6.0);
        // Independently integrate each physical bath's q and q/T_bath. The jump
        // changes neither cell storage nor reservoir entropy at a single instant.
        for (t, weight) in [(first, 1.0 - GAMMA), (second, GAMMA)] {
            for (g, bath) in [(0.5, ambient), (1.0, reservoir)] {
                let q = g * (bath - t);
                heat += dt * weight * q;
                entropy_input += dt * weight * q / bath;
                generation += dt * weight * g * (bath - t).powi(2) / (bath * t);
            }
        }
        temperature = second;
    }
    let storage_entropy = 6.0 * (temperature / 600.0).ln();
    near(next.temperature_k(id(0)).unwrap(), temperature, 1e-11);
    near(report.energy.external_heat_j, heat, 1e-9);
    near(
        report.energy.stored_energy_change_j,
        6.0 * (temperature - 600.0),
        1e-9,
    );
    near(
        report.entropy.reservoir_entropy_input_j_k,
        entropy_input,
        1e-11,
    );
    near(report.entropy.conduction_generation_j_k, generation, 1e-11);
    near(
        report.entropy.thermal_storage_change_j_k,
        storage_entropy,
        1e-11,
    );
    near(
        report.entropy.temporal_defect_j_k,
        storage_entropy - entropy_input - generation,
        1e-11,
    );
    assert!(generation > 0.0);
    assert!(!report.entropy.has_unclassified_sources);
    assert_eq!(report.energy.events, 1);
    assert_eq!(report.energy.substeps, 2);
}

#[test]
fn gray_exterior_receives_scheduled_baths_and_preserves_interval_accounts() {
    let law = Law::constant(50.0, 2.0);
    let (thermal, observer) = network(
        &law,
        &[600.0],
        &[1.0],
        &[],
        None,
        vec![change(0.3, baths(450.0, None, 0.0))],
    );
    let observer = observer
        .with_ambient_faces(
            &thermal,
            &[ConductiveAmbient {
                cell_id: id(0),
                half_distance_m: 0.5,
                surface_resistance_m2_k_w: 0.02,
            }],
            || Ok(()),
        )
        .unwrap();
    let gray = GrayAmbientBoundary::new(&law, 0.8, 5.0, GrayAmbientControls::default()).unwrap();
    let geometry = ImplicitBoundaryGeometry {
        area_m2: 4.0,
        half_distance_m: 0.5,
        surface_resistance_m2_k_w: 0.02,
    };
    let cold = gray
        .evaluate(50.0 * 600.0, 300.0, geometry, &mut || Ok(()))
        .unwrap();
    let warm = gray
        .evaluate(50.0 * 600.0, 450.0, geometry, &mut || Ok(()))
        .unwrap();
    assert!(warm.outward_power_w > 0.0 && warm.outward_power_w < cold.outward_power_w);
    let surface = warm.surface_temperature_k;
    near(
        warm.outward_power_w,
        4.0 * (5.0 * (surface - 450.0)
            + 0.8 * 5.670_374_419e-8 * (surface.powi(4) - 450.0_f64.powi(4))),
        1e-7,
    );
    let solver = chain(&thermal, &observer, &law)
        .with_boundary_laws(&thermal, &[], &[(id(0), &gray)], || Ok(()))
        .unwrap();
    let mut workspace = ImplicitConductionWorkspace::new(1).unwrap();
    let (scheduled, report) = solver
        .advance_sdirk2_to(&thermal, 0.7, controls(0.1), &mut workspace, || Ok(()))
        .unwrap();

    // An independent trajectory has no scheduled event. Change its bath between
    // calls and reprepare the boundary owner at the actual new physical state.
    let mut manual = thermal.clone();
    manual.schedule.truncate(manual.cursor);
    let (mut first, first_report) = solver
        .advance_sdirk2_to(&manual, 0.3, controls(0.1), &mut workspace, || Ok(()))
        .unwrap();
    first.cells.get_mut(&id(0)).unwrap().boundary.ambient_k = 450.0;
    let prepared = chain(&first, &observer, &law)
        .with_boundary_laws(&first, &[], &[(id(0), &gray)], || Ok(()))
        .unwrap();
    let (manual, second_report) = prepared
        .advance_sdirk2_to(&first, 0.7, controls(0.1), &mut workspace, || Ok(()))
        .unwrap();
    near(
        scheduled.temperature_k(id(0)).unwrap(),
        manual.temperature_k(id(0)).unwrap(),
        1e-10,
    );
    near(
        report.energy.external_heat_j,
        first_report.energy.external_heat_j + second_report.energy.external_heat_j,
        1e-8,
    );
    near(
        report.entropy.reservoir_entropy_input_j_k,
        first_report.entropy.reservoir_entropy_input_j_k
            + second_report.entropy.reservoir_entropy_input_j_k,
        1e-10,
    );
    near(
        report.entropy.conduction_generation_j_k,
        first_report.entropy.conduction_generation_j_k
            + second_report.entropy.conduction_generation_j_k,
        1e-10,
    );
    near(
        report.energy.stored_energy_change_j,
        50.0 * (scheduled.temperature_k(id(0)).unwrap() - 600.0),
        1e-8,
    );
    assert_eq!(report.energy.events, 1);
    assert_eq!(thermal.cells()[&id(0)].boundary.ambient_k, 300.0);
}

#[test]
fn overridden_exchange_coefficients_refuse_before_any_boundary_interval() {
    let law = Law::constant(1.0, 1.0);
    for change_ambient_htc in [true, false] {
        let mut unsupported = boundary(Some(350.0));
        if change_ambient_htc {
            unsupported.htc_w_m2_k = 2.0;
        } else {
            unsupported.reservoir.as_mut().unwrap().conductance_w_k = 2.0;
        }
        let (thermal, observer) = network(
            &law,
            &[400.0],
            &[1.0],
            &[],
            Some(300.0),
            vec![
                change(0.25, baths(320.0, Some(350.0), 0.0)),
                change(1.0, unsupported),
            ],
        );
        let observer = observer
            .with_ambient_faces(
                &thermal,
                &[ConductiveAmbient {
                    cell_id: id(0),
                    half_distance_m: 1.0,
                    surface_resistance_m2_k_w: 0.1,
                }],
                || Ok(()),
            )
            .unwrap();
        let before = thermal.cells().clone();
        let solver = chain(&thermal, &observer, &law);
        let mut workspace = ImplicitConductionWorkspace::new(1).unwrap();
        let step = ImplicitConductionControls {
            maximum_residual_evaluations: 1,
            ..controls(0.1)
        };
        assert!(matches!(
            solver.advance_sdirk2_to(&thermal, 1.0, step, &mut workspace, || Ok(())),
            Err(ImplicitConductionError::Invalid(_))
        ));
        let mut adaptive = adaptive_controls(0.1);
        adaptive.step = step;
        assert!(matches!(
            solver.advance_sdirk2_adaptive_to(&thermal, 1.0, adaptive, &mut workspace, || Ok(())),
            Err(AdaptiveImplicitConductionError::Implicit(
                ImplicitConductionError::Invalid(_)
            ))
        ));
        assert_eq!(thermal.cells(), &before);
        assert_eq!(thermal.time_s(), 0.0);
    }
}
