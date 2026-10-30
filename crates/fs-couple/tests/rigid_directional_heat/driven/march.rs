//! G1/G4: driven continuation, continuous source convergence and span rollback.
use super::*;
use fs_couple::friction_heat::directional::trajectory::march::{
    DirectionalMarchControls, DirectionalMarchTolerance,
};

fn march_controls(step: f64) -> DirectionalMarchControls {
    DirectionalMarchControls {
        initial_step_s: step,
        minimum_step_s: 1e-7,
        maximum_step_s: step,
        mean_speed_bracket_m_s: [0.0, 10.0],
        maximum_attempts: 500,
        maximum_accepted_phases: 800,
        tolerance: DirectionalMarchTolerance {
            linear_momentum_ns: 1e-5,
            angular_momentum_nms: 1e-6,
            cell_energy_j: 1e-5,
            distance_m: 1e-6,
            dissipated_j: 1e-5,
            relative: 1e-4,
        },
        energy_tolerance_j: 1e-7,
        source_path_tolerance_m: 1e-8,
    }
}

fn phase_controls() -> DirectionalHeatSolve {
    let mut solve = controls();
    solve.motion.distance_bracket_m = [0.0, 2.0];
    solve
}

fn relative_x(
    bodies: &[RigidBodyState],
    props: &[MassProperties],
    contact: RigidDirectionalContact,
) -> f64 {
    let a = bodies[0]
        .point_kinematics(props[0], contact.binding.arm_a_body)
        .unwrap();
    let b = bodies[1]
        .point_kinematics(props[1], contact.binding.arm_b_body)
        .unwrap();
    a.point_velocity_world.x - b.point_velocity_world.x
}

#[test]
fn driven_span_breaks_away_arrests_and_holds_between_kinetic_and_static_thresholds() {
    let (a, ca) = material(0);
    let (b, cb) = material(1);
    let cards = [a, b];
    let curves = [ca, cb];
    let q = query();
    let card = interface(&cards, [scalar(0.2), scalar(0.4)], 400.0);
    let law = law(&card, &q, [1.0, 0.0]);
    let mut history = law.initial_state();
    let mut mode = EllipticDrivenMode::Resting;
    let mut heat = thermal(&curves, [300.0; 2]);
    let route = FrictionHeatRoute::from_link(&heat, 7, CONDUCTIVITY).unwrap();
    let (mut bodies, props, contact) = resting_mechanics(0.0);
    let original = bodies;
    let mobility = 101.0 / 120.0;
    let driven_time = 0.02;
    let velocity = (10.0 - 2.0 * mobility) * driven_time;
    let driven_path = 0.5 * velocity * driven_time;
    let report = advance_source_directional_driven_heat(
        &mut bodies,
        &props,
        &mut history,
        &mut mode,
        &law,
        contact,
        drive(20.0),
        &mut heat,
        &route,
        binding(&cards, &curves, &card, &q),
        driven_time,
        phase_controls(),
        driven_controls(),
        march_controls(0.01),
        300.0,
        |_, p| transport(&cards, p),
        &CancelGate::new(),
    )
    .unwrap();
    assert_eq!(mode, EllipticDrivenMode::Sliding);
    assert!(report.phases.len() >= 4);
    assert!(
        report
            .phases
            .iter()
            .all(|phase| { phase.report.mechanical.driven.mode == EllipticDrivenMode::Sliding })
    );
    close(relative_x(&bodies, &props, contact), velocity);
    close(report.distance_m, driven_path);
    close(report.net_slip_m[0], driven_path);
    close(report.dissipated_j, 2.0 * driven_path);
    close(
        report.external_mechanical_work_j,
        90.0 * driven_time * driven_time,
    );
    close(bodies[0].linear_momentum_world().x, 18.0 * driven_time);
    close(bodies[1].linear_momentum_world().x, 2.0 * driven_time);
    for body in &bodies[..2] {
        close(body.angular_momentum_body().y, -0.2 * driven_time);
    }
    close(report.mechanical_energy_residual_j, 0.0);
    close(report.thermal_energy_residual_j, 0.0);
    close(report.combined_energy_residual_j, 0.0);

    let before_brake = bodies;
    let brake_time = 0.4;
    let stop = velocity / (2.0 * mobility - 1.0);
    let brake_path = 0.5 * velocity * stop;
    let hold_traction = 1.0 / mobility;
    let impulse = 2.0 * stop + hold_traction * (brake_time - stop);
    // During kinetic braking A's 2N drive balances its 2N friction. After
    // arrest it accelerates at (2-q)/(2kg) while contact slip remains zero.
    let initial_a_speed = before_brake[0].linear_momentum_world().x / 2.0;
    let hold_acceleration = (2.0 - hold_traction) / 2.0;
    let brake_work = 2.0 * initial_a_speed * brake_time
        + hold_acceleration * (brake_time - stop) * (brake_time - stop);
    let report = advance_source_directional_driven_heat(
        &mut bodies,
        &props,
        &mut history,
        &mut mode,
        &law,
        contact,
        drive(2.0),
        &mut heat,
        &route,
        binding(&cards, &curves, &card, &q),
        brake_time,
        phase_controls(),
        driven_controls(),
        march_controls(0.04),
        300.0,
        |_, p| transport(&cards, p),
        &CancelGate::new(),
    )
    .unwrap();
    assert_eq!(mode, EllipticDrivenMode::Resting);
    close(relative_x(&bodies, &props, contact), 0.0);
    close(report.distance_m, brake_path);
    close(report.dissipated_j, 2.0 * brake_path);
    close(report.external_mechanical_work_j, brake_work);
    close(
        bodies[0].linear_momentum_world().x,
        before_brake[0].linear_momentum_world().x + 2.0 * brake_time - impulse,
    );
    close(
        bodies[1].linear_momentum_world().x,
        before_brake[1].linear_momentum_world().x + impulse,
    );
    for i in 0..2 {
        close(
            bodies[i].angular_momentum_body().y,
            before_brake[i].angular_momentum_body().y - 0.1 * impulse,
        );
    }
    let arrest = report
        .phases
        .iter()
        .position(|phase| phase.report.mechanical.driven.first_arrest_time_s.is_some())
        .unwrap();
    let event = &report.phases[arrest];
    assert!(event.report.mechanical.driven.arrest_time_is_exact);
    close(
        event.start_s + event.report.mechanical.driven.first_arrest_time_s.unwrap(),
        driven_time + stop,
    );
    assert!(arrest + 1 < report.phases.len());
    for phase in &report.phases[arrest + 1..] {
        assert_eq!(
            phase.report.mechanical.driven.mode,
            EllipticDrivenMode::Resting
        );
        assert_eq!(phase.report.heat.dissipated_j, 0.0);
        assert_eq!(phase.report.actual_distance_m, 0.0);
    }
    close(report.combined_energy_residual_j, 0.0);
    let cumulative_path = driven_path + brake_path;
    for id in CELLS {
        close(
            heat.temperature_k(id).unwrap(),
            300.0 + 2.0 * cumulative_path / 0.08,
        );
    }

    // This drive exceeds kinetic resistance but is below static capacity.
    // Retained Resting, including any actual momentum roundoff, must hold.
    let before_hold = bodies;
    let old_temperatures = CELLS.map(|id| heat.temperature_k(id).unwrap());
    let hold_time = 0.03;
    let traction = 6.0 / mobility;
    assert!(traction > 2.0 && traction < 9.0);
    let a_speed = before_hold[0].linear_momentum_world().x / 2.0;
    let hold_work = 12.0 * (a_speed * hold_time + 0.25 * (12.0 - traction) * hold_time * hold_time);
    let temperature = heat.temperature_k(CELLS[0]).unwrap();
    let report = advance_source_directional_driven_heat(
        &mut bodies,
        &props,
        &mut history,
        &mut mode,
        &law,
        contact,
        drive(12.0),
        &mut heat,
        &route,
        binding(&cards, &curves, &card, &q),
        hold_time,
        phase_controls(),
        driven_controls(),
        march_controls(0.01),
        temperature,
        |_, p| transport(&cards, p),
        &CancelGate::new(),
    )
    .unwrap();
    assert_eq!(mode, EllipticDrivenMode::Resting);
    close(report.distance_m, 0.0);
    close(report.dissipated_j, 0.0);
    close(report.external_mechanical_work_j, hold_work);
    close(report.kinetic_energy_change_j, hold_work);
    close(report.combined_energy_residual_j, 0.0);
    for i in 0..2 {
        close(heat.temperature_k(CELLS[i]).unwrap(), old_temperatures[i]);
    }
    for phase in &report.phases {
        assert_eq!(
            phase.report.mechanical.driven.mode,
            EllipticDrivenMode::Resting
        );
        assert_eq!(phase.report.actual_distance_m, 0.0);
    }
    for source in history.friction_sources() {
        close(
            source.state_coordinate("path").unwrap(),
            0.25 + cumulative_path,
        );
        close(source.state_coordinate("speed").unwrap(), 0.0);
    }
    assert_eq!(bodies[2], original[2]);
    close(heat.time_s(), driven_time + brake_time + hold_time);
}

#[test]
fn three_driven_time_rungs_converge_to_exact_speed_dependent_friction_ode() {
    let (a, ca) = material(0);
    let (b, cb) = material(1);
    let cards = [a, b];
    let curves = [ca, cb];
    let q = query();
    let coefficient = curve("speed", SPEED, [(0.0, 0.2), (10.0, 0.7)]);
    let card = interface(&cards, [coefficient, scalar(0.4)], 400.0);
    let law = law(&card, &q, [1.0, 0.0]);
    let initial = law.initial_state();
    let base = thermal(&curves, [300.0; 2]);
    let route = FrictionHeatRoute::from_link(&base, 7, CONDUCTIVITY).unwrap();
    let (original, props, contact) = resting_mechanics(0.0);
    let mobility: f64 = 101.0 / 120.0;
    let duration: f64 = 0.2;
    // F=2+.5v and v'=10-WF: an independently integrated linear-drag ODE.
    // Its exact path also fixes friction heat and full external body work.
    let rate = 0.5 * mobility;
    let terminal = (10.0 - 2.0 * mobility) / rate;
    let expected_velocity = terminal * (1.0 - (-rate * duration).exp());
    let expected_path = terminal * duration - expected_velocity / rate;
    let expected_heat = (10.0 * expected_path - 0.5 * expected_velocity.powi(2)) / mobility;
    let expected_work =
        100.0 * duration * duration - 10.0 * (5.0 * duration * duration - expected_path) / mobility;
    let mut errors = Vec::new();
    for steps in [4, 8, 16] {
        let mut bodies = original;
        let mut history = initial.clone();
        let mut mode = EllipticDrivenMode::Resting;
        let mut heat = base.clone();
        let mut control = march_controls(duration / f64::from(steps));
        // Preserve these three grids while still evaluating each independent
        // full/two-half comparison. Mean-speed sampling is midpoint for this
        // linear source law, so this particular ODE converges at second order.
        control.tolerance = DirectionalMarchTolerance {
            linear_momentum_ns: 1.0,
            angular_momentum_nms: 1.0,
            cell_energy_j: 1.0,
            distance_m: 1.0,
            dissipated_j: 1.0,
            relative: 0.0,
        };
        let report = advance_source_directional_driven_heat(
            &mut bodies,
            &props,
            &mut history,
            &mut mode,
            &law,
            contact,
            drive(20.0),
            &mut heat,
            &route,
            binding(&cards, &curves, &card, &q),
            duration,
            phase_controls(),
            driven_controls(),
            control,
            300.0,
            |_, p| transport(&cards, p),
            &CancelGate::new(),
        )
        .unwrap();
        let velocity = relative_x(&bodies, &props, contact);
        let impulse = (10.0 * duration - velocity) / mobility;
        errors.push([
            (velocity - expected_velocity).abs(),
            (report.distance_m - expected_path).abs(),
            (report.dissipated_j - expected_heat).abs(),
            (report.external_mechanical_work_j - expected_work).abs(),
        ]);
        assert_eq!(mode, EllipticDrivenMode::Sliding);
        assert_eq!(report.rejected_attempts, 0);
        assert!(report.maximum_accepted_normalized_defect <= 1.0);
        close(
            bodies[0].linear_momentum_world().x,
            20.0 * duration - impulse,
        );
        close(bodies[1].linear_momentum_world().x, impulse);
        for body in &bodies[..2] {
            close(body.angular_momentum_body().y, -0.1 * impulse);
        }
        close(
            report.dissipated_j,
            (10.0 * report.distance_m - 0.5 * velocity * velocity) / mobility,
        );
        close(
            report.external_mechanical_work_j,
            100.0 * duration * duration
                - 10.0 * (5.0 * duration * duration - report.distance_m) / mobility,
        );
        close(report.mechanical_energy_residual_j, 0.0);
        close(report.thermal_energy_residual_j, 0.0);
        close(report.combined_energy_residual_j, 0.0);
        for id in CELLS {
            close(
                heat.temperature_k(id).unwrap(),
                300.0 + report.dissipated_j / 0.08,
            );
        }
        let last = report.phases.last().unwrap();
        let mean_speed = last.report.actual_distance_m / (last.end_s - last.start_s);
        close(
            history.friction_sources()[0]
                .property("kinetic-x")
                .unwrap()
                .value_si(),
            0.2 + 0.05 * mean_speed,
        );
        for (i, source) in history.friction_sources().iter().enumerate() {
            close(source.state_coordinate("speed").unwrap(), mean_speed);
            close(
                source.state_coordinate("path").unwrap(),
                0.25 + report.distance_m,
            );
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
                card.claims()
                    .verify_receipt(&property.answer().receipt)
                    .unwrap();
            }
        }
        assert!(report.source_path_error_m < control.source_path_tolerance_m);
        assert_eq!(heat.time_s(), duration);
        assert_eq!(bodies[2], original[2]);
    }
    for pair in errors.windows(2) {
        for metric in 0..4 {
            let ratio = pair[0][metric] / pair[1][metric];
            assert!(ratio > 3.4 && ratio < 4.6, "ODE convergence: {errors:?}");
        }
    }
    for (error, tolerance) in errors[2].into_iter().zip([1.1e-6, 2.5e-6, 4e-5, 3e-5]) {
        assert!(error < tolerance, "finest driven ODE errors: {errors:?}");
    }
}

#[test]
fn later_driven_source_cancel_and_phase_budget_failures_rollback_mode_and_all_owners() {
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
    let (original, props, contact) = resting_mechanics(0.0);
    for failure in 0..3 {
        let mut bodies = original;
        let mut history = initial.clone();
        let mut mode = EllipticDrivenMode::Resting;
        let mut heat = base.clone();
        let mut control = march_controls(0.01);
        if failure == 2 {
            control.maximum_accepted_phases = 4;
        }
        let latest_query = std::cell::Cell::new(0.0_f64);
        let gate = CancelGate::new();
        let result = advance_source_directional_driven_heat(
            &mut bodies,
            &props,
            &mut history,
            &mut mode,
            &law,
            contact,
            drive(20.0),
            &mut heat,
            &route,
            binding(&cards, &curves, &card, &q),
            0.06,
            phase_controls(),
            driven_controls(),
            control,
            300.0,
            |endpoint, phases| {
                latest_query.set(latest_query.get().max(endpoint.time_s()));
                if endpoint.time_s() > 0.025 {
                    if failure == 0 {
                        return Err(ThermalMaterialError::Invalid("later driven source refusal"));
                    }
                    if failure == 1 {
                        gate.request();
                    }
                }
                transport(&cards, phases)
            },
            &gate,
        );
        assert!(
            result.is_err(),
            "late failure {failure} unexpectedly published"
        );
        assert!(
            latest_query.get() >= 0.02,
            "must advance private driven phases"
        );
        if failure == 2 {
            assert!(matches!(
                result,
                Err(DrivenDirectionalHeatMarchError::AcceptedPhaseBudget)
            ));
        }
        assert_eq!(mode, EllipticDrivenMode::Resting);
        assert_eq!(bodies, original);
        assert_eq!(history, initial);
        assert_eq!(heat.cells(), base.cells());
        assert_eq!(heat.links(), base.links());
        assert_eq!(heat.time_s(), base.time_s());
    }
    let mut bodies = original;
    let mut history = initial;
    let mut mode = EllipticDrivenMode::Resting;
    let mut heat = base.clone();
    let duration = 0.06;
    let report = advance_source_directional_driven_heat(
        &mut bodies,
        &props,
        &mut history,
        &mut mode,
        &law,
        contact,
        drive(20.0),
        &mut heat,
        &route,
        binding(&cards, &curves, &card, &q),
        duration,
        phase_controls(),
        driven_controls(),
        march_controls(0.01),
        300.0,
        |_, p| transport(&cards, p),
        &CancelGate::new(),
    )
    .unwrap();
    assert_eq!(mode, EllipticDrivenMode::Sliding);
    assert!(report.phases.len() > 4);
    close(
        report.distance_m,
        0.5 * (10.0 - 2.0 * 101.0 / 120.0) * duration * duration,
    );
    close(report.dissipated_j, 2.0 * report.distance_m);
    close(
        report.external_mechanical_work_j,
        90.0 * duration * duration,
    );
    close(report.combined_energy_residual_j, 0.0);
    assert_eq!(bodies[2], original[2]);
    assert_eq!(heat.time_s(), duration);
}
