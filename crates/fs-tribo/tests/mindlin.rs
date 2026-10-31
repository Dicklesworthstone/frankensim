//! G0/G1/G3/G5: analytic Hertz--Cattaneo--Mindlin force, memory and work tests.
//! All elastic and interface parameters below are synthetic fixtures.
use fs_math::det;
use fs_tribo::mindlin::{MAX_HISTORY_INTERVALS, MindlinDeresiewiczPoint, MindlinPointParameters};
use fs_tribo::{TangentialLaw, TangentialSolveConfig, TangentialStepInput, TriboError};

fn parameters() -> MindlinPointParameters {
    MindlinPointParameters {
        normal_force_n: 10.0,
        contact_radius_m: 1e-3,
        shear_moduli_pa: [30e9, 20e9],
        poisson_ratios: [0.25, 0.3],
        friction_coefficient: 0.6,
        maximum_history_intervals: 32,
    }
}

fn model() -> MindlinDeresiewiczPoint {
    MindlinDeresiewiczPoint::new(parameters()).unwrap()
}

fn input(law: &MindlinDeresiewiczPoint, increment: f64, compliance: f64) -> TangentialStepInput {
    TangentialStepInput {
        normal_n: law.parameters().normal_force_n,
        slip_increment_m: [increment, 0.0],
        dt_s: 1e-3,
        contact_temperature_k: Some(300.0),
        nominal_pressure_pa: None,
        compliance_m_per_n: compliance,
        solve: TangentialSolveConfig {
            max_iterations: 64,
            maximum_force_n: 10.0,
            force_absolute_tolerance_n: 1e-13,
            force_relative_tolerance: 1e-13,
        },
    }
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual={actual:e}, expected={expected:e}, tolerance={tolerance:e}"
    );
}

#[test]
fn returned_traction_satisfies_its_actual_compliant_load_line_after_reversals() {
    let law = model();
    let virgin = law.initial_state();
    let loaded = law.advance(&virgin,
        &input(&law, 0.7 * law.gross_slip_displacement_m(), 0.0)).unwrap().state;
    let reversed = law.advance(&loaded,
        &input(&law, -0.2 * law.gross_slip_displacement_m(), 0.0)).unwrap().state;
    // The returned effort drives the body, so the original law must agree at
    // free-C*RETURNED_TRACTION, not merely at the last internal root iterate.
    for (history, state) in [&virgin, &loaded, &reversed].into_iter().enumerate() {
        for ratio in [0.1, 1.0, 10.0, 100.0, 10000.0] {
            let compliance = ratio / law.initial_stiffness_n_per_m();
            for fraction in [-2.0, -0.7, -0.01, 0.01, 0.31, 1.4] {
                let free = fraction * law.gross_slip_displacement_m();
                let mut request = input(&law, free, compliance);
                request.solve.force_absolute_tolerance_n = 1e-10;
                request.solve.force_relative_tolerance = 1e-9;
                let step = law.advance(state, &request).unwrap();
                let actual = free - compliance * step.traction_n[0];
                let mut prescribed = request;
                prescribed.slip_increment_m[0] = actual;
                prescribed.compliance_m_per_n = 0.0;
                let checked = law.advance(state, &prescribed).unwrap();
                let residual = step.traction_n[0] - checked.traction_n[0];
                let tolerance = request.solve.force_absolute_tolerance_n
                    + request.solve.force_relative_tolerance
                        * step.traction_n[0].abs().max(checked.traction_n[0].abs());
                assert!(residual.abs() <= tolerance,
                    "history={history}, C*k0={ratio}, free/d={fraction}, residual={residual:e}, tolerance={tolerance:e}");
                assert!(step.dissipated_j >= 0.0);
                close(step.closure_residual_j, 0.0,
                    2e-14 * law.gross_slip_capacity_n() * law.gross_slip_displacement_m());
            }
        }
    }
}

fn backbone(law: &MindlinDeresiewiczPoint, displacement: f64) -> f64 {
    let fraction = displacement.abs() / law.gross_slip_displacement_m();
    let magnitude = if fraction >= 1.0 {
        law.gross_slip_capacity_n()
    } else {
        law.gross_slip_capacity_n() * (1.0 - det::pow(1.0 - fraction, 1.5))
    };
    displacement.signum() * magnitude
}

#[test]
fn g1_monotonic_curve_matches_johnson_7_42() {
    let law = model();
    let p = law.parameters();
    let compliance = (2.0 - p.poisson_ratios[0]) / p.shear_moduli_pa[0]
        + (2.0 - p.poisson_ratios[1]) / p.shear_moduli_pa[1];
    for fraction in [0.0, 1e-5, 0.01, 0.1, 0.5, 0.9, 0.999] {
        // Force-controlled Johnson 7.42 independently gives the imposed slip.
        let displacement = 3.0 * p.friction_coefficient * p.normal_force_n
            / (16.0 * p.contact_radius_m)
            * compliance
            * (1.0 - det::pow(1.0 - fraction, 2.0 / 3.0));
        for sign in [-1.0, 1.0] {
            let step = law
                .advance(&law.initial_state(), &input(&law, sign * displacement, 0.0))
                .unwrap();
            let force = sign * fraction * p.friction_coefficient * p.normal_force_n;
            close(
                step.state.endpoint_force_n(),
                force,
                1e-10 * force.abs().max(1e-14),
            );
            assert!(step.dissipated_j >= 0.0);
            close(
                step.closure_residual_j,
                0.0,
                2e-14 * law.gross_slip_capacity_n() * law.gross_slip_displacement_m(),
            );
        }
    }
}

#[test]
fn g1_initial_stiffness_is_8_a_gstar() {
    let law = model();
    let p = parameters();
    let exact = 8.0 * p.contact_radius_m
        / ((2.0 - p.poisson_ratios[0]) / p.shear_moduli_pa[0]
            + (2.0 - p.poisson_ratios[1]) / p.shear_moduli_pa[1]);
    close(law.initial_stiffness_n_per_m(), exact, 1e-12 * exact);
    let displacement = 1e-10 * law.gross_slip_displacement_m();
    let step = law
        .advance(&law.initial_state(), &input(&law, displacement, 0.0))
        .unwrap();
    // The nonlinear secant tends to k0; the held effort tends to k0*u/2.
    close(
        step.state.endpoint_force_n() / displacement,
        exact,
        3e-11 * exact,
    );
    close(
        step.traction_n[0] / displacement,
        0.5 * exact,
        2e-11 * exact,
    );
}

#[test]
fn g1_masing_reversal_and_nested_return_point_memory() {
    let law = model();
    let amplitude = 0.75 * law.gross_slip_displacement_m();
    let peak = law
        .advance(&law.initial_state(), &input(&law, amplitude, 0.0))
        .unwrap()
        .state;
    for fraction in [0.0, 0.1, 0.5, 1.0, 1.5, 2.0] {
        let retreat = fraction * amplitude;
        let step = law.advance(&peak, &input(&law, -retreat, 0.0)).unwrap();
        let expected = peak.endpoint_force_n() - 2.0 * backbone(&law, retreat / 2.0);
        close(
            step.state.endpoint_force_n(),
            expected,
            2e-13 * law.gross_slip_capacity_n(),
        );
    }
    let below = law
        .advance(&peak, &input(&law, -0.5 * amplitude, 0.0))
        .unwrap()
        .state;
    let nested = law
        .advance(&below, &input(&law, 0.2 * amplitude, 0.0))
        .unwrap()
        .state;
    let returned = law
        .advance(&nested, &input(&law, -0.2 * amplitude, 0.0))
        .unwrap()
        .state;
    close(returned.endpoint_force_n(), below.endpoint_force_n(), 1e-13);
    close(
        law.stored_energy_j(&returned),
        law.stored_energy_j(&below),
        1e-13 * law.gross_slip_capacity_n() * law.gross_slip_displacement_m(),
    );
    // The closed nested loop has erased its memory: future response agrees.
    for increment in [-0.4 * amplitude, 0.8 * amplitude] {
        let a = law
            .advance(&returned, &input(&law, increment, 0.0))
            .unwrap();
        let b = law.advance(&below, &input(&law, increment, 0.0)).unwrap();
        close(
            a.state.endpoint_force_n(),
            b.state.endpoint_force_n(),
            1e-13,
        );
        close(
            a.dissipated_j,
            b.dissipated_j,
            1e-13 * law.gross_slip_capacity_n() * law.gross_slip_displacement_m(),
        );
    }
}

fn cycle(law: &MindlinDeresiewiczPoint, force_fraction: f64, half_steps: usize) -> (f64, f64) {
    let amplitude =
        law.gross_slip_displacement_m() * (1.0 - det::pow(1.0 - force_fraction, 2.0 / 3.0));
    let mut state = law
        .advance(&law.initial_state(), &input(law, amplitude, 0.0))
        .unwrap()
        .state;
    let start_energy = law.stored_energy_j(&state);
    let mut heat = 0.0;
    let mut endpoint_work = 0.0;
    for sign in [-1.0, 1.0] {
        for _ in 0..half_steps {
            let increment = sign * 2.0 * amplitude / half_steps as f64;
            let old_force = state.endpoint_force_n();
            let step = law.advance(&state, &input(law, increment, 0.0)).unwrap();
            heat += step.dissipated_j;
            endpoint_work += 0.5 * (old_force + step.state.endpoint_force_n()) * increment;
            close(
                step.stored_change_j,
                law.stored_energy_j(&step.state) - law.stored_energy_j(&state),
                2e-14 * law.gross_slip_capacity_n() * law.gross_slip_displacement_m(),
            );
            close(
                step.closure_residual_j,
                0.0,
                1e-14 * law.gross_slip_capacity_n() * law.gross_slip_displacement_m(),
            );
            state = step.state;
        }
    }
    close(
        law.stored_energy_j(&state),
        start_energy,
        1e-12 * law.gross_slip_capacity_n() * law.gross_slip_displacement_m(),
    );
    (heat, endpoint_work)
}

#[test]
fn g1_masing_loop_energy_matches_closed_form_and_cubic_limit() {
    let law = model();
    let capacity = law.gross_slip_capacity_n();
    let critical = law.gross_slip_displacement_m();
    for fraction in [0.1, 0.5, 0.9] {
        // Integral of the independent Masing backbone loop, factored to
        // preserve its Q^3 small-amplitude asymptote. h is the stick-radius ratio.
        let h = det::pow(1.0 - fraction, 1.0 / 3.0);
        let one_minus_h = fraction / (1.0 + h + h * h);
        let exact = 0.8
            * capacity
            * critical
            * one_minus_h
            * one_minus_h
            * one_minus_h
            * (1.0 + 3.0 * h + h * h);
        let mut errors = Vec::new();
        for steps in [100, 200, 400] {
            let (heat, endpoint_work) = cycle(&law, fraction, steps);
            // The constitutive interval integrals are exact; no step ladder
            // is needed for them. A separate endpoint-work observer is O(dt^2).
            close(heat, exact, 2e-11 * exact);
            errors.push((endpoint_work - exact).abs());
        }
        for pair in errors.windows(2) {
            let ratio = pair[0] / pair[1];
            assert!(
                (3.9..4.1).contains(&ratio),
                "fraction={fraction}, ratio={ratio}"
            );
        }
    }
    let a = cycle(&law, 0.001, 100).0;
    let b = cycle(&law, 0.002, 100).0;
    let slope = det::ln(b / a) / det::ln(2.0);
    close(slope, 3.0, 0.002);
}

#[test]
fn g0_gross_slip_rest_preload_and_frictionless_limit() {
    let law = model();
    let d = law.gross_slip_displacement_m();
    let loaded = law
        .advance(&law.initial_state(), &input(&law, 2.0 * d, 0.0))
        .unwrap();
    close(
        loaded.state.endpoint_force_n(),
        law.gross_slip_capacity_n(),
        1e-12,
    );
    let sliding = law
        .advance(&loaded.state, &input(&law, 0.5 * d, 0.0))
        .unwrap();
    close(sliding.traction_n[0], law.gross_slip_capacity_n(), 1e-12);
    close(
        sliding.dissipated_j,
        law.gross_slip_capacity_n() * 0.5 * d,
        1e-14 * d,
    );
    close(sliding.stored_change_j, 0.0, 1e-14 * d);
    let resting = law.advance(&sliding.state, &input(&law, 0.0, 0.0)).unwrap();
    assert_eq!(resting.state, sliding.state);
    assert_eq!(resting.dissipated_j, 0.0);
    assert_eq!(resting.traction_n[0], resting.state.endpoint_force_n());

    let frictionless = MindlinDeresiewiczPoint::new(MindlinPointParameters {
        friction_coefficient: 0.0,
        ..parameters()
    })
    .unwrap();
    let step = frictionless
        .advance(
            &frictionless.initial_state(),
            &input(&frictionless, 0.1, 1.0),
        )
        .unwrap();
    assert_eq!(step.traction_n, [0.0; 2]);
    assert_eq!(step.regime, fs_tribo::FrictionRegime::Sliding);
    assert_eq!(step.slip_increment_m, [0.1, 0.0]);
    assert_eq!(step.dissipated_j, 0.0);
}

#[test]
fn g1_compliant_load_line_uses_path_work_and_releases_stored_energy() {
    let law = model();
    let d = law.gross_slip_displacement_m();
    let compliance = 0.75 / law.initial_stiffness_n_per_m();
    let initial = law.initial_state();
    let loaded = law.advance(&initial, &input(&law, d, compliance)).unwrap();
    let reference = law
        .advance(&initial, &input(&law, loaded.slip_increment_m[0], 0.0))
        .unwrap();
    close(loaded.traction_n[0], reference.traction_n[0], 1e-12);
    close(
        loaded.slip_increment_m[0],
        d - compliance * loaded.traction_n[0],
        1e-13 * d,
    );
    assert_ne!(loaded.traction_n[0], loaded.state.endpoint_force_n());
    assert!(
        loaded.iterations < 12,
        "analytic tangent should accelerate the scalar solve"
    );
    // No free driver increment: retained shear still unloads into the body.
    let relaxing = law
        .advance(&loaded.state, &input(&law, 0.0, compliance))
        .unwrap();
    assert!(relaxing.slip_increment_m[0] < 0.0);
    assert!(relaxing.stored_change_j < 0.0);
    assert!(relaxing.traction_n[0] * relaxing.slip_increment_m[0] < 0.0);
    assert!(relaxing.dissipated_j >= 0.0);
    close(
        relaxing.closure_residual_j,
        0.0,
        1e-13 * law.gross_slip_capacity_n() * d,
    );
}

#[test]
fn g0_bounded_reversal_memory_refuses_atomically_and_wipes_on_gross_slip() {
    // These synthetic inputs make d=1, so the nested reversal locations are exact.
    let law = MindlinDeresiewiczPoint::new(MindlinPointParameters {
        normal_force_n: 2.0,
        contact_radius_m: 1.0,
        shear_moduli_pa: [0.65625; 2],
        poisson_ratios: [0.25; 2],
        friction_coefficient: 0.5,
        maximum_history_intervals: 3,
    })
    .unwrap();
    assert_eq!(law.gross_slip_displacement_m(), 1.0);
    let a = law
        .advance(&law.initial_state(), &input(&law, 0.5, 0.0))
        .unwrap()
        .state;
    assert_eq!(a.history_intervals(), 2);
    let b = law.advance(&a, &input(&law, -0.25, 0.0)).unwrap().state;
    assert_eq!(b.history_intervals(), 3);
    let checkpoint = b.clone();
    assert!(matches!(
        law.advance(&b, &input(&law, 0.125, 0.0)),
        Err(TriboError::OutsideApplicability {
            field: "mindlin_history_intervals",
            ..
        })
    ));
    assert_eq!(b, checkpoint);
    let erased = law.advance(&b, &input(&law, 3.0, 0.0)).unwrap().state;
    assert_eq!(erased.history_intervals(), 1);
    close(erased.endpoint_force_n(), 1.0, 1e-14);
}

#[test]
fn g0_domain_force_iteration_and_checkpoint_refusals() {
    let law = model();
    let state = law.initial_state();
    let mut x = input(&law, law.gross_slip_displacement_m(), 0.0);
    x.normal_n = x.normal_n.next_up();
    assert!(matches!(
        law.advance(&state, &x),
        Err(TriboError::OutsideApplicability {
            field: "mindlin_fixed_normal_force_n",
            ..
        })
    ));
    x = input(&law, 0.0, 0.0);
    x.slip_increment_m[1] = 1e-20;
    assert!(law.advance(&state, &x).is_err());
    let other = MindlinDeresiewiczPoint::new(MindlinPointParameters {
        contact_radius_m: 2e-3,
        ..parameters()
    })
    .unwrap();
    assert!(
        law.advance(&other.initial_state(), &input(&law, 0.0, 0.0))
            .is_err()
    );
    assert!(law.stored_energy_j(&other.initial_state()).is_nan());
    x = input(&law, law.gross_slip_displacement_m(), 0.0);
    x.solve.maximum_force_n = 0.01;
    assert!(law.advance(&state, &x).is_err());
    x = input(
        &law,
        law.gross_slip_displacement_m(),
        1.0 / law.initial_stiffness_n_per_m(),
    );
    x.solve.max_iterations = 1;
    assert!(matches!(
        law.advance(&state, &x),
        Err(TriboError::TangentialSolve { .. })
    ));
    for value in [0, MAX_HISTORY_INTERVALS + 1] {
        assert!(
            MindlinDeresiewiczPoint::new(MindlinPointParameters {
                maximum_history_intervals: value,
                ..parameters()
            })
            .is_err()
        );
    }
    assert!(
        MindlinDeresiewiczPoint::new(MindlinPointParameters {
            maximum_history_intervals: MAX_HISTORY_INTERVALS,
            ..parameters()
        })
        .is_ok()
    );
}

#[test]
fn g3_step_partition_invariance_and_g5_checkpoint_replay() {
    let law = model();
    let d = law.gross_slip_displacement_m();
    let path = [0.7, -0.6, 0.2, -0.1, 0.3, -1.5, 2.5];
    let run = |subdivisions: usize| {
        let mut state = law.initial_state();
        let mut heat = 0.0;
        for normalized_increment in path {
            for _ in 0..subdivisions {
                let step = law
                    .advance(
                        &state,
                        &input(&law, normalized_increment * d / subdivisions as f64, 0.0),
                    )
                    .unwrap();
                heat += step.dissipated_j;
                state = step.state;
            }
        }
        (state, heat)
    };
    let whole = run(1);
    for count in [4, 16] {
        let divided = run(count);
        close(
            whole.0.endpoint_force_n(),
            divided.0.endpoint_force_n(),
            1e-13,
        );
        close(whole.1, divided.1, 1e-13 * law.gross_slip_capacity_n() * d);
        close(
            law.stored_energy_j(&whole.0),
            law.stored_energy_j(&divided.0),
            1e-13 * law.gross_slip_capacity_n() * d,
        );
    }
    assert_eq!(run(4), run(4));
    let saved = whole.0.clone();
    let a = law.advance(&whole.0, &input(&law, -0.3 * d, 0.0)).unwrap();
    let b = law.advance(&saved, &input(&law, -0.3 * d, 0.0)).unwrap();
    assert_eq!(a, b);
}
