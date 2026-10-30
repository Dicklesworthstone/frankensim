//! G1/G3/G4: source/thermal time refinement, whole-span energy and rollback.
use super::*;
use fs_couple::friction_heat::directional::trajectory::march::*;

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
fn three_time_rungs_converge_to_independent_thermal_friction_ode() {
    let (a, ca) = material(0);
    let (b, cb) = material(1);
    let cards = [a, b];
    let curves = [ca, cb];
    let q = query();
    // The finite cell capacities are proportional to their face conductances.
    // Equal initial bulk T therefore stays equal, exactly: T=300+D/C_total.
    // mu(Ta)=.2+.01*(Ta-300), D=(1-v²)/(2W), so dv/dt=-A+Bv².
    // This independent continuous-time solution tests the time-discretization
    // axis; matching a single backward-Euler step would not establish that.
    let coefficient = curve("Ta", Temperature::DIMS, [(290.0, 0.1), (330.0, 0.5)]);
    let card = interface(&cards, [coefficient.clone(), coefficient], 400.0);
    let law = law(&card, &q, [1.0, 0.0]);
    let initial = law.initial_state();
    let base = thermal(&curves, [300.0; 2]);
    let route = FrictionHeatRoute::from_link(&base, 7, CONDUCTIVITY).unwrap();
    let (original, props, contact) = mechanics(false);
    let w: f64 = 101.0 / 120.0;
    let beta: f64 = 10.0 * 0.01 / (2.0 * 0.08);
    let alpha = 2.0 * w + beta;
    let duration = 0.2;
    let expected_v = (alpha / beta).sqrt()
        * ((beta / alpha).sqrt().atanh() - (alpha * beta).sqrt() * duration).tanh();
    let mut errors = Vec::new();
    for steps in [4, 8, 16] {
        let mut bodies = original;
        let mut state = initial.clone();
        let mut heat = base.clone();
        let mut control = march_controls(duration / f64::from(steps));
        // Keep the declared three time grids; the controller still performs
        // its independent full/two-half comparison and publishes finer states.
        control.tolerance = DirectionalMarchTolerance {
            linear_momentum_ns: 1.0,
            angular_momentum_nms: 1.0,
            cell_energy_j: 1.0,
            distance_m: 1.0,
            dissipated_j: 1.0,
            relative: 0.0,
        };
        let report = advance_source_directional_heat(
            &mut bodies,
            &props,
            &mut state,
            &law,
            contact,
            &mut heat,
            &route,
            binding(&cards, &curves, &card, &q),
            duration,
            heat_controls(),
            trajectory_controls(),
            control,
            300.0,
            |_, p| transport(&cards, p),
            &CancelGate::new(),
        )
        .unwrap();
        let v = relative_x(&bodies, &props, contact);
        errors.push((v - expected_v).abs());
        let dissipated = (1.0 - v * v) / (2.0 * w);
        close(report.dissipated_j, dissipated);
        for id in CELLS {
            close(heat.temperature_k(id).unwrap(), 300.0 + dissipated / 0.08);
        }
        assert_eq!(report.rejected_attempts, 0);
        assert!(report.maximum_accepted_normalized_defect <= 1.0);
        close(report.combined_energy_residual_j, 0.0);
        close(report.mechanical_energy_residual_j, 0.0);
        close(report.thermal_energy_residual_j, 0.0);
        assert_eq!(heat.time_s(), duration);
        assert_eq!(bodies[2], original[2]);
        for pair in report.phases.windows(2) {
            assert_eq!(pair[0].end_s, pair[1].start_s);
        }
        assert_eq!(report.phases[0].start_s, 0.0);
        assert_eq!(report.phases.last().unwrap().end_s, duration);
        assert!(report.source_path_error_m < control.source_path_tolerance_m);
        close(
            state.friction_sources()[0]
                .state_coordinate("path")
                .unwrap(),
            0.25 + report.distance_m,
        );
    }
    assert!(errors[2] < 1.5e-3, "thermal ODE errors: {errors:?}");
    for pair in errors.windows(2) {
        let ratio = pair[0] / pair[1];
        assert!(
            ratio > 1.6 && ratio < 2.4,
            "first-order thermal feedback ladder: {errors:?}"
        );
    }
}

#[test]
fn adaptive_refinement_reaches_arrest_and_keeps_conducting_without_extra_heat() {
    let (a, ca) = material(0);
    let (b, cb) = material(1);
    let cards = [a, b];
    let curves = [ca, cb];
    let q = query();
    let card = interface(&cards, [scalar(0.2), scalar(0.4)], 400.0);
    let law = law(&card, &q, [1.0, 0.0]);
    let mut state = law.initial_state();
    let mut heat = thermal(&curves, [290.0, 310.0]);
    let route = FrictionHeatRoute::from_link(&heat, 7, CONDUCTIVITY).unwrap();
    let (mut bodies, props, contact) = mechanics(false);
    let original = bodies;
    let mut control = march_controls(0.2);
    control.tolerance.cell_energy_j = 1e-3;
    let report = advance_source_directional_heat(
        &mut bodies,
        &props,
        &mut state,
        &law,
        contact,
        &mut heat,
        &route,
        binding(&cards, &curves, &card, &q),
        0.9,
        heat_controls(),
        trajectory_controls(),
        control,
        305.0,
        |_, p| transport(&cards, p),
        &CancelGate::new(),
    )
    .unwrap();
    let w = 101.0 / 120.0;
    close(relative_x(&bodies, &props, contact), 0.0);
    close(report.distance_m, 1.0 / (4.0 * w));
    close(report.dissipated_j, 0.5 / w);
    close(report.net_slip_m[0], report.distance_m);
    close(report.net_slip_m[1], 0.0);
    close(report.external_heat_j, 0.0);
    close(report.combined_energy_residual_j, 0.0);
    assert!(report.rejected_attempts > 0);
    assert!(
        report.thermal_trials
            > report
                .phases
                .iter()
                .map(|p| p.report.thermal_trials)
                .sum::<usize>()
    );
    assert!(report.maximum_accepted_normalized_defect <= 1.0);
    let arrested = report
        .phases
        .iter()
        .position(|p| {
            p.report.mechanical.relaxation.kinetic_phase_end_s.is_some()
                && p.report.mechanical.relaxation.velocity_after_m_s == [0.0; 2]
        })
        .unwrap();
    let phase = &report.phases[arrested];
    close(
        phase.start_s
            + phase
                .report
                .mechanical
                .relaxation
                .kinetic_phase_end_s
                .unwrap(),
        1.0 / (2.0 * w),
    );
    assert!(arrested + 1 < report.phases.len());
    for p in &report.phases[arrested + 1..] {
        assert!(p.report.heat.dissipated_j < 1e-20);
        assert!(p.report.actual_distance_m < 1e-20);
    }
    // Continuous two-cell exchange decays as exp(-100t); at .9 s both are
    // essentially the conserved weighted initial mean plus mechanical heat.
    let equilibrium = 305.0 + report.dissipated_j / 0.08;
    for id in CELLS {
        assert!((heat.temperature_k(id).unwrap() - equilibrium).abs() < 2e-4);
    }
    assert_eq!(bodies[2], original[2]);
    assert_eq!(heat.time_s(), 0.9);
    assert_eq!(heat.links().len(), 1);
    assert!(report.source_path_error_m <= control.source_path_tolerance_m);
}

#[test]
fn later_source_cancel_and_budget_failures_rollback_the_whole_span_and_retry() {
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
    for mode in 0..4 {
        let mut bodies = original;
        let mut state = initial.clone();
        let mut heat = base.clone();
        let mut control = march_controls(0.01);
        if mode == 2 {
            control.maximum_attempts = 2;
        }
        if mode == 3 {
            control.maximum_accepted_phases = 4;
        }
        let latest_query = std::cell::Cell::new(0.0_f64);
        let gate = CancelGate::new();
        let result = advance_source_directional_heat(
            &mut bodies,
            &props,
            &mut state,
            &law,
            contact,
            &mut heat,
            &route,
            binding(&cards, &curves, &card, &q),
            0.06,
            heat_controls(),
            trajectory_controls(),
            control,
            300.0,
            |endpoint, p| {
                latest_query.set(latest_query.get().max(endpoint.time_s()));
                if endpoint.time_s() > 0.025 {
                    if mode == 0 {
                        return Err(ThermalMaterialError::Invalid("later source refusal"));
                    }
                    if mode == 1 {
                        gate.request();
                    }
                }
                transport(&cards, p)
            },
            &gate,
        );
        assert!(result.is_err(), "mode {mode} unexpectedly published");
        assert!(
            latest_query.get() >= 0.02,
            "must exercise private accepted substeps"
        );
        if mode == 2 {
            assert!(matches!(
                result,
                Err(DirectionalHeatMarchError::AttemptBudget)
            ));
        }
        if mode == 3 {
            assert!(matches!(
                result,
                Err(DirectionalHeatMarchError::AcceptedPhaseBudget)
            ));
        }
        assert_eq!(bodies, original);
        assert_eq!(state, initial);
        assert_eq!(heat.cells(), base.cells());
        assert_eq!(heat.links(), base.links());
        assert_eq!(heat.time_s(), base.time_s());
    }
    let mut bodies = original;
    let mut state = initial.clone();
    let mut heat = base.clone();
    let report = advance_source_directional_heat(
        &mut bodies,
        &props,
        &mut state,
        &law,
        contact,
        &mut heat,
        &route,
        binding(&cards, &curves, &card, &q),
        0.06,
        heat_controls(),
        trajectory_controls(),
        march_controls(0.01),
        300.0,
        |_, p| transport(&cards, p),
        &CancelGate::new(),
    )
    .unwrap();
    close(report.combined_energy_residual_j, 0.0);
    assert_eq!(heat.time_s(), 0.06);
    assert_ne!(state, initial);
    assert_eq!(bodies[2], original[2]);
}

#[test]
fn unresolved_time_accuracy_and_incompatible_speed_domains_refuse_atomically() {
    let (a, ca) = material(0);
    let (b, cb) = material(1);
    let cards = [a, b];
    let curves = [ca, cb];
    let q = query();
    let card = interface(&cards, [scalar(0.2), scalar(0.4)], 400.0);
    let law = law(&card, &q, [1.0, 0.0]);
    let initial = law.initial_state();
    let base = thermal(&curves, [290.0, 310.0]);
    let route = FrictionHeatRoute::from_link(&base, 7, CONDUCTIVITY).unwrap();
    let (original, props, contact) = mechanics(false);
    for mode in 0..3 {
        let mut bodies = original;
        let mut state = initial.clone();
        let mut heat = base.clone();
        let mut control = march_controls(0.01);
        let mut phase = heat_controls();
        control.minimum_step_s = 0.01;
        if mode == 0 {
            control.tolerance.cell_energy_j = 1e-12;
            control.tolerance.relative = 0.0;
        } else if mode == 1 {
            control.mean_speed_bracket_m_s = [2.0, 3.0];
        } else {
            phase.motion.distance_bracket_m = [f64::NAN, f64::NAN];
        }
        let result = advance_source_directional_heat(
            &mut bodies,
            &props,
            &mut state,
            &law,
            contact,
            &mut heat,
            &route,
            binding(&cards, &curves, &card, &q),
            0.02,
            phase,
            trajectory_controls(),
            control,
            305.0,
            |_, p| transport(&cards, p),
            &CancelGate::new(),
        );
        if mode == 0 {
            assert!(matches!(
                result,
                Err(DirectionalHeatMarchError::Unresolved { .. })
            ));
        } else {
            assert!(result.is_err());
        }
        assert_eq!(bodies, original);
        assert_eq!(state, initial);
        assert_eq!(heat.cells(), base.cells());
        assert_eq!(heat.time_s(), 0.0);
    }
}
