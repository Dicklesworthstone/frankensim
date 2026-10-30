//! G1/G3/G4: true curved path and finite arrest in the sourced heat transaction.
use super::*;
use fs_couple::friction_heat::directional::trajectory::*;
use fs_math::det;
use fs_tribo::tangential::elliptic::relaxation::EllipticRelaxationControls;

#[path = "trajectory/march.rs"]
mod march;

fn trajectory_controls() -> EllipticRelaxationControls {
    EllipticRelaxationControls {
        time_tolerance_s: 1e-10,
        travel_tolerance_m: 1e-10,
        maximum_evaluations: 40_000,
        maximum_root_iterations: 96,
    }
}
fn heat_controls() -> DirectionalHeatSolve {
    let mut solve = controls();
    solve.motion.distance_bracket_m = [0.0, 2.0];
    solve
}

/// Independently eliminate two linear cells and their common massless face.
/// This is backward Euler at the given step, not a continuous-time heat oracle.
fn finite_cell_oracle(old: [f64; 2], work: f64, dt: f64) -> (f64, [f64; 2], [f64; 2]) {
    let g = [2.0, 6.0];
    let effective = [0, 1].map(|i| dt * g[i] / (1.0 + dt * g[i] / CAPACITY[i]));
    let ti = (work + effective[0] * old[0] + effective[1] * old[1]) / (effective[0] + effective[1]);
    let heat = [0, 1].map(|i| effective[i] * (ti - old[i]));
    let next = [0, 1].map(|i| old[i] + heat[i] / CAPACITY[i]);
    (ti, next, heat)
}

#[test]
fn finite_arrest_heats_the_original_cells_and_conduction_continues_at_rest() {
    let (a, ca) = material(0);
    let (b, cb) = material(1);
    let cards = [a, b];
    let curves = [ca, cb];
    let q = query();
    let card = interface(&cards, [scalar(0.2), scalar(0.4)], 400.0);
    let law = law(&card, &q, [1.0, 0.0]);
    let mut state = law.initial_state();
    let initial = state.clone();
    let temperatures = [290.0, 310.0];
    let mut heat = thermal(&curves, temperatures);
    let route = FrictionHeatRoute::from_link(&heat, 7, CONDUCTIVITY).unwrap();
    let (mut bodies, props, contact) = mechanics(false);
    let original = bodies;
    let dt = 0.8;
    let w = 101.0 / 120.0;
    let force = 2.0;
    let stop = 1.0 / (w * force);
    let distance = 0.5 * stop;
    let work = 0.5 / w;
    let (expected_ti, expected_t, expected_heat) = finite_cell_oracle(temperatures, work, dt);
    let report = apply_source_directional_trajectory_heat(
        &mut bodies,
        &props,
        &mut state,
        &law,
        contact,
        &mut heat,
        &route,
        binding(&cards, &curves, &card, &q),
        dt,
        heat_controls(),
        trajectory_controls(),
        &[0.0; 2],
        305.0,
        |_, p| transport(&cards, p),
        &CancelGate::new(),
    )
    .unwrap();
    close(
        report.mechanical.relaxation.kinetic_phase_end_s.unwrap(),
        stop,
    );
    close(report.mechanical.slip_velocity_after_m_s[0], 0.0);
    close(report.actual_distance_m, distance);
    close(report.mechanical.relaxation.slip_increment_m[0], distance);
    close(report.heat.dissipated_j, work);
    close(report.mechanical.kinetic_energy_change_j, -work);
    close(report.heat.interface_temperature_k, expected_ti);
    close(report.combined_energy_residual_j, 0.0);
    for i in 0..2 {
        close(heat.temperature_k(CELLS[i]).unwrap(), expected_t[i]);
        close(report.heat.body_heat_j[i], expected_heat[i]);
        close(report.source[i].query().axes()["T"], expected_t[i]);
        let source = &state.friction_sources()[i];
        close(source.state_coordinate("speed").unwrap(), distance / dt);
        close(source.state_coordinate("path").unwrap(), 0.25 + distance);
        for property in source.properties() {
            assert_eq!(
                property.answer().receipt.selected,
                initial.friction_sources()[i]
                    .property(property.requirement().name())
                    .unwrap()
                    .answer()
                    .receipt
                    .selected
            );
        }
    }
    assert_eq!(bodies[2], original[2]);
    assert_eq!(heat.links().len(), 1);
    let next_dt = 0.2;
    let (next_ti, next_t, next_heat) = finite_cell_oracle(expected_t, 0.0, next_dt);
    let stationary = apply_source_directional_trajectory_heat(
        &mut bodies,
        &props,
        &mut state,
        &law,
        contact,
        &mut heat,
        &route,
        binding(&cards, &curves, &card, &q),
        next_dt,
        heat_controls(),
        trajectory_controls(),
        &[0.0; 2],
        next_ti,
        |_, p| transport(&cards, p),
        &CancelGate::new(),
    )
    .unwrap();
    close(stationary.heat.dissipated_j, 0.0);
    close(stationary.actual_distance_m, 0.0);
    close(stationary.heat.interface_temperature_k, next_ti);
    close(heat.time_s(), 1.0);
    assert!(stationary.heat.body_heat_j[0] > 0.0);
    assert!(stationary.heat.body_heat_j[1] < 0.0);
    for i in 0..2 {
        close(heat.temperature_k(CELLS[i]).unwrap(), next_t[i]);
        close(stationary.heat.body_heat_j[i], next_heat[i]);
        close(
            state.friction_sources()[i]
                .state_coordinate("path")
                .unwrap(),
            0.25 + distance,
        );
    }
    close(stationary.combined_energy_residual_j, 0.0);
}

#[test]
fn original_temperature_weakening_changes_stop_time_and_path_in_the_same_phase() {
    let (a, ca) = material(0);
    let (b, cb) = material(1);
    let cards = [a, b];
    let curves = [ca, cb];
    let q = query();
    let card = interface(
        &cards,
        [
            curve("Ti", Temperature::DIMS, [(280.0, 0.22), (340.0, 0.16)]),
            scalar(0.4),
        ],
        400.0,
    );
    let law = law(&card, &q, [1.0, 0.0]);
    let mut state = law.initial_state();
    let mut heat = thermal(&curves, [300.0; 2]);
    let route = FrictionHeatRoute::from_link(&heat, 7, CONDUCTIVITY).unwrap();
    let (mut bodies, props, contact) = mechanics(false);
    let dt = 0.8;
    let w = 101.0 / 120.0;
    let work = 0.5 / w;
    let (ti, _, _) = finite_cell_oracle([300.0; 2], work, dt);
    let force = 2.0 - 0.01 * (ti - 300.0);
    let stop = 1.0 / (w * force);
    let distance = 0.5 * stop;
    let report = apply_source_directional_trajectory_heat(
        &mut bodies,
        &props,
        &mut state,
        &law,
        contact,
        &mut heat,
        &route,
        binding(&cards, &curves, &card, &q),
        dt,
        heat_controls(),
        trajectory_controls(),
        &[0.0; 2],
        ti,
        |_, p| transport(&cards, p),
        &CancelGate::new(),
    )
    .unwrap();
    assert!(stop > 0.61 && stop < dt);
    close(
        report.mechanical.relaxation.kinetic_phase_end_s.unwrap(),
        stop,
    );
    close(report.actual_distance_m, distance);
    close(report.heat.interface_temperature_k, ti);
    close(report.heat.dissipated_j, work);
    close(report.combined_energy_residual_j, 0.0);
    close(
        state.friction_sources()[0]
            .property("kinetic-x")
            .unwrap()
            .value_si(),
        force / 10.0,
    );
}

#[test]
fn curved_unsigned_source_path_and_full_tensor_impulse_are_frame_covariant() {
    let (a, ca) = material(0);
    let (b, cb) = material(1);
    let cards = [a, b];
    let curves = [ca, cb];
    let q = query();
    let card = interface(
        &cards,
        [
            curve("speed", SPEED, [(0.0, 0.1), (20.0, 0.7)]),
            curve("path", LENGTH, [(0.0, 0.2), (10.0, 0.6)]),
        ],
        400.0,
    );
    let mut reference = None;
    for rotated in [false, true] {
        let (mut bodies, props, mut contact) = mechanics(true);
        let axis = if rotated {
            contact.frame = TangentFrame::new([0.0, 0.0, 1.0], [0.6, 0.8, 0.0]).unwrap();
            [1.0, 0.0]
        } else {
            [0.6, 0.8]
        };
        let law = law(&card, &q, axis);
        let mut state = law.initial_state();
        let mut heat = thermal(&curves, [300.0; 2]);
        let route = FrictionHeatRoute::from_link(&heat, 7, CONDUCTIVITY).unwrap();
        let dt = 2.0;
        let report = apply_source_directional_trajectory_heat(
            &mut bodies,
            &props,
            &mut state,
            &law,
            contact,
            &mut heat,
            &route,
            binding(&cards, &curves, &card, &q),
            dt,
            heat_controls(),
            trajectory_controls(),
            &[0.0; 2],
            305.0,
            |_, p| transport(&cards, p),
            &CancelGate::new(),
        )
        .unwrap();
        let native = &report.mechanical.relaxation;
        assert!(native.kinetic_phase_end_s.is_some());
        assert!(report.mechanical.mobility_per_kg[0][1].abs() > 1e-4);
        assert!(
            native.distance_m - det::hypot(native.slip_increment_m[0], native.slip_increment_m[1])
                > 1e-4
        );
        assert!(report.source_motion_residual_m <= 0.0);
        assert!(report.source_motion_residual_m.abs() < 3e-12);
        for source in state.friction_sources() {
            close(
                source.state_coordinate("speed").unwrap(),
                native.distance_m / dt,
            );
            close(
                source.state_coordinate("path").unwrap(),
                0.25 + native.distance_m,
            );
            close(
                source.state_coordinate("Ta").unwrap(),
                report.phases[0].temperature_k(),
            );
            close(
                source.state_coordinate("Tb").unwrap(),
                report.phases[1].temperature_k(),
            );
        }
        let values = [
            report.mechanical.force_impulse_on_a_world_ns.x,
            report.mechanical.force_impulse_on_a_world_ns.y,
            report.actual_distance_m,
            report.heat.dissipated_j,
            report.heat.interface_temperature_k,
            native.kinetic_phase_end_s.unwrap(),
        ];
        if let Some(old) = reference {
            for (a, b) in values.into_iter().zip(old) {
                close(a, b);
            }
        } else {
            reference = Some(values);
        }
        close(report.combined_energy_residual_j, 0.0);
    }
}

#[test]
fn failed_sourced_arrest_publishes_neither_body_history_nor_heat_and_can_retry() {
    let (a, ca) = material(0);
    let (b, cb) = material(1);
    let cards = [a, b];
    let curves = [ca, cb];
    let q = query();
    let card = interface(&cards, [scalar(0.2), scalar(0.4)], 400.0);
    let law = law(&card, &q, [1.0, 0.0]);
    let initial = law.initial_state();
    let base = thermal(&curves, [300.0; 2]);
    let route = FrictionHeatRoute::from_link(&base, 7, CONDUCTIVITY).unwrap();
    let (original, props, contact) = mechanics(false);
    for mode in 0..5 {
        let mut bodies = original;
        let mut state = initial.clone();
        let mut heat = base.clone();
        let gate = CancelGate::new();
        let mut solve = heat_controls();
        if mode == 2 {
            solve.motion.distance_bracket_m = [0.0, 1e-4];
        }
        if mode == 3 {
            solve.maximum_trials = 1;
        }
        let mut material = binding(&cards, &curves, &card, &q);
        if mode == 4 {
            material.caloric.swap(0, 1);
        }
        let result = apply_source_directional_trajectory_heat(
            &mut bodies,
            &props,
            &mut state,
            &law,
            contact,
            &mut heat,
            &route,
            material,
            0.8,
            solve,
            trajectory_controls(),
            &[0.0; 2],
            305.0,
            |_, p| {
                if mode == 0 {
                    return Err(ThermalMaterialError::Invalid("transport refused"));
                }
                if mode == 1 {
                    gate.request();
                }
                transport(&cards, p)
            },
            &gate,
        );
        assert!(result.is_err(), "failure mode {mode} unexpectedly accepted");
        assert_eq!(bodies, original);
        assert_eq!(state, initial);
        assert_eq!(heat.cells(), base.cells());
        assert_eq!(heat.links(), base.links());
        assert_eq!(heat.time_s(), base.time_s());
    }
    let mut bodies = original;
    let mut state = initial.clone();
    let mut heat = base.clone();
    let report = apply_source_directional_trajectory_heat(
        &mut bodies,
        &props,
        &mut state,
        &law,
        contact,
        &mut heat,
        &route,
        binding(&cards, &curves, &card, &q),
        0.8,
        heat_controls(),
        trajectory_controls(),
        &[0.0; 2],
        305.0,
        |_, p| transport(&cards, p),
        &CancelGate::new(),
    )
    .unwrap();
    assert!(report.mechanical.relaxation.kinetic_phase_end_s.is_some());
    assert_ne!(state, initial);
    assert_eq!(bodies[2], original[2]);
    close(report.combined_energy_residual_j, 0.0);
}
