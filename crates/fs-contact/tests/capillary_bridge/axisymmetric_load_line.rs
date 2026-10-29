//! G1/G3/G4: actual 3-D capillary work and finite-body mechanical compliance.
use super::*;
use fs_contact::capillary_bridge::axisymmetric::load_line::*;
use fs_contact::capillary_bridge::axisymmetric::*;
use fs_contact::capillary_bridge::load_line::BridgeLoadLineConfig;

fn model3d(sigma: f64) -> SourceAxisymmetricCapillaryBridge {
    let (liquid, lower) = sources(1.1, sigma, 300.0, true);
    let (_, upper) = sources(1.3, sigma, 300.0, true);
    SourceAxisymmetricCapillaryBridge::new(
        liquid,
        [lower[0].clone(), upper[1].clone()],
        VOLUME,
        AxisymmetricBranch {
            minimum_suction_parameter: 0.0,
            maximum_suction_parameter: 32.0,
            minimum_stability_margin: 0.01,
        },
        AxisymmetricLimits {
            minimum_gap_m: 1e-8,
            maximum_gap_m: 0.1,
            maximum_wetted_radius_m: 0.1,
            maximum_suction_pa: 1e8,
            maximum_abs_force_n: 10.0,
            maximum_abs_surface_energy_j: 1.0,
            maximum_abs_step_work_j: 1.0,
            gravity_m_s2: 0.0,
            maximum_bond_number: 0.01,
        },
    )
    .unwrap()
}
fn numerical3d() -> AxisymmetricNumerics {
    AxisymmetricNumerics {
        quadrature_relative_tolerance: 8e-14,
        gap_relative_tolerance: 3e-14,
        maximum_quadrature_panels: 4096,
        maximum_root_iterations: 80,
    }
}
fn controls() -> AxisymmetricLoadLineConfig {
    AxisymmetricLoadLineConfig {
        mechanical: BridgeLoadLineConfig {
            max_iterations: 32,
            max_evaluations: 160,
            max_backtracks: 16,
            gap_absolute_tolerance_m: 3e-15,
            gap_relative_tolerance: 1e-10,
            maximum_gap_increment_m: 1e-4,
            minimum_compatibility_slope: 0.05,
        },
        work: AxisymmetricWorkBudget {
            equilibrium: numerical3d(),
            absolute_tolerance_j: 1e-18,
            relative_tolerance: 1e-5,
        },
        force_absolute_tolerance_n: 1e-10,
        force_relative_tolerance: 1e-6,
    }
}

#[test]
fn g1_compliant_axisymmetric_equilibrium_has_actual_work_and_mechanical_tangent() {
    let source = model3d(0.05);
    let gate = CancelGate::new_clock_free();
    let previous = source.solve_at_suction(2.0, numerical3d(), &gate).unwrap();
    assert_eq!(source.evaluate(previous).unwrap(), *previous.profile());
    let compliance = 0.005;
    let solver = AxisymmetricBridgeLoadLine::new(&source, compliance, controls()).unwrap();
    assert!(std::ptr::eq(solver.source(), &source));
    assert_eq!(solver.config(), controls());
    assert_eq!(solver.compliance_m_per_n(), compliance);
    let free = 1.04 * previous.profile().gap_m;
    let result = solver.trial(previous, free, &gate).unwrap();
    let next = result.bridge.state();
    let gap = next.profile().gap_m;
    assert_eq!(
        result.compatibility_residual_m,
        (gap - free) - compliance * result.bridge.mean_normal_force_n
    );
    assert!(
        result.compatibility_residual_m.abs() + result.estimated_compatibility_error_m
            <= result.compatibility_tolerance_m
    );
    assert!(result.estimated_force_error_n <= result.force_error_tolerance_n);
    assert_eq!(
        result.bridge.work_on_liquid_j,
        next.profile().surface_free_energy_j - previous.profile().surface_free_energy_j
    );
    near(
        -result.bridge.mean_normal_force_n * (gap - previous.profile().gap_m),
        result.bridge.work_on_liquid_j,
        3e-15,
        1e-20,
    );
    let delta = 1e-4 * previous.profile().gap_m;
    let lo = solver.trial(previous, free - delta, &gate).unwrap();
    let hi = solver.trial(previous, free + delta, &gate).unwrap();
    near(
        (hi.bridge.state().profile().gap_m - lo.bridge.state().profile().gap_m) / (2.0 * delta),
        1.0 / result.compatibility_slope,
        2e-6,
        1.0,
    );
}

#[test]
fn g1_two_finite_bodies_exchange_axisymmetric_surface_work_and_momentum() {
    let source = model3d(0.05);
    let gate = CancelGate::new_clock_free();
    let mut state = source.solve_at_suction(2.0, numerical3d(), &gate).unwrap();
    let (ma, mb, dt) = (2e-4, 3e-4, 2e-5);
    let compliance = 0.5 * dt * dt * (1.0 / ma + 1.0 / mb);
    let solver = AxisymmetricBridgeLoadLine::new(&source, compliance, controls()).unwrap();
    let mut xa = 0.0;
    let mut xb = state.profile().gap_m;
    let mut va = 0.0;
    let mut vb = 0.0;
    let kinetic = |va: f64, vb: f64| 0.5 * ma * va * va + 0.5 * mb * vb * vb;
    let initial_gap = xb;
    let initial = state.profile().surface_free_energy_j;
    for _ in 0..200 {
        let old_kinetic = kinetic(va, vb);
        let result = solver
            .trial(state, xb - xa + dt * (vb - va), &gate)
            .unwrap();
        let force = result.bridge.mean_normal_force_n;
        xa += dt * va - 0.5 * dt * dt * force / ma;
        xb += dt * vb + 0.5 * dt * dt * force / mb;
        va -= dt * force / ma;
        vb += dt * force / mb;
        let actual = solver.advance_actual(state, xb - xa, &gate).unwrap();
        near(actual.mean_normal_force_n, force, 1e-6, 1e-7);
        near(actual.state().profile().gap_m, xb - xa, 2e-13, 0.0);
        near(ma * va + mb * vb, 0.0, 2e-13, 1e-5);
        near(
            kinetic(va, vb) - old_kinetic + actual.work_on_liquid_j,
            0.0,
            1e-8,
            1e-8,
        );
        near(
            kinetic(va, vb) + actual.state().profile().surface_free_energy_j,
            initial,
            1e-8,
            1e-8,
        );
        assert!(result.evaluations <= controls().mechanical.max_evaluations);
        assert_eq!(actual.state().profile().liquid_volume_m3, VOLUME);
        state = actual.state();
    }
    assert!(state.profile().gap_m < 0.95 * initial_gap);
    assert!(kinetic(va, vb) > 0.0);
}

#[test]
fn g3_original_force_dynamics_converge_under_time_refinement() {
    let source = model3d(0.05);
    let gate = CancelGate::new_clock_free();
    let start = source.solve_at_suction(2.0, numerical3d(), &gate).unwrap();
    let inverse_mass = 1.0 / 2e-4 + 1.0 / 3e-4;
    let total_time = 0.002;
    let force = |h| {
        source
            .solve_at_gap(h, numerical3d(), &gate)
            .unwrap()
            .profile()
            .normal_force_n
    };
    // Independent classical RK4 of the original continuous relative equation.
    let mut reference = [start.profile().gap_m, 0.0];
    let dt = total_time / 256.0;
    let rhs = |y: [f64; 2]| [y[1], inverse_mass * force(y[0])];
    let plus = |y: [f64; 2], k: [f64; 2], s: f64| [y[0] + s * k[0], y[1] + s * k[1]];
    for _ in 0..256 {
        let a = rhs(reference);
        let b = rhs(plus(reference, a, 0.5 * dt));
        let c = rhs(plus(reference, b, 0.5 * dt));
        let d = rhs(plus(reference, c, dt));
        reference = std::array::from_fn(|i| {
            reference[i] + dt * (a[i] + 2.0 * b[i] + 2.0 * c[i] + d[i]) / 6.0
        });
    }
    let mut errors = [0.0; 3];
    for (i, steps) in [20, 40, 80].into_iter().enumerate() {
        let dt = total_time / f64::from(steps);
        let solver =
            AxisymmetricBridgeLoadLine::new(&source, 0.5 * dt * dt * inverse_mass, controls())
                .unwrap();
        let mut state = start;
        let mut gap = start.profile().gap_m;
        let mut velocity = 0.0;
        for _ in 0..steps {
            let step = solver.trial(state, gap + dt * velocity, &gate).unwrap();
            let force = step.bridge.mean_normal_force_n;
            gap += dt * velocity + 0.5 * dt * dt * inverse_mass * force;
            velocity += dt * inverse_mass * force;
            state = solver.advance_actual(state, gap, &gate).unwrap().state();
        }
        errors[i] = (gap - reference[0]).abs() + total_time * (velocity - reference[1]).abs();
    }
    assert!(
        errors[0] / errors[1] > 3.5 && errors[1] / errors[2] > 3.5,
        "{errors:?}"
    );
    assert!(
        errors[0] / errors[1] < 4.5 && errors[1] / errors[2] < 4.5,
        "{errors:?}"
    );
}

#[test]
fn g1_prescribed_motion_exact_rest_and_unresolved_finite_work_stay_distinct() {
    let source = model3d(0.05);
    let gate = CancelGate::new_clock_free();
    let previous = source.solve_at_suction(2.0, numerical3d(), &gate).unwrap();
    let solver = AxisymmetricBridgeLoadLine::new(&source, 0.0, controls()).unwrap();
    let rest = solver
        .trial(previous, previous.profile().gap_m, &gate)
        .unwrap();
    assert_eq!(rest.bridge.state(), previous);
    assert_eq!(rest.bridge.work_on_liquid_j, 0.0);
    assert_eq!(
        rest.bridge.mean_normal_force_n,
        previous.profile().normal_force_n
    );
    assert_eq!(
        rest.bridge.mean_force_gap_derivative_n_per_m,
        0.5 * previous.profile().force_gap_derivative_n_per_m
    );
    assert_eq!(rest.iterations, 0);
    let target = 0.98 * previous.profile().gap_m;
    let prescribed = solver.trial(previous, target, &gate).unwrap();
    near(
        prescribed.bridge.state().profile().gap_m,
        target,
        3e-14,
        0.0,
    );
    assert!(prescribed.bridge.work_on_liquid_j < 0.0);
    assert!(matches!(
        solver.trial(previous, previous.profile().gap_m.next_up(), &gate),
        Err(AxisymmetricLoadLineError::Bridge(
            AxisymmetricError::WorkUnresolved { .. }
        ))
    ));
    assert_eq!(rest.bridge.state(), previous);
}

#[test]
fn g0_g4_source_force_accuracy_attractive_slope_and_budget_refusals_are_atomic() {
    let source = model3d(0.05);
    let other = model3d(0.06);
    let gate = CancelGate::new_clock_free();
    let previous = source.solve_at_suction(2.0, numerical3d(), &gate).unwrap();
    let h = previous.profile().gap_m;
    let solver = AxisymmetricBridgeLoadLine::new(&source, 0.001, controls()).unwrap();
    let expected = solver.trial(previous, 1.02 * h, &gate).unwrap();
    let other_state = other.solve_at_suction(2.0, numerical3d(), &gate).unwrap();
    assert_eq!(
        source.evaluate(other_state),
        Err(AxisymmetricError::StateMismatch)
    );
    assert!(matches!(
        solver.trial(other_state, h, &gate),
        Err(AxisymmetricLoadLineError::Bridge(
            AxisymmetricError::StateMismatch
        ))
    ));
    let cancelled = CancelGate::new_clock_free();
    cancelled.request();
    assert!(matches!(
        solver.trial(previous, 1.02 * h, &cancelled),
        Err(AxisymmetricLoadLineError::Cancelled)
    ));

    let mut limited = controls();
    limited.mechanical.max_evaluations = 1;
    assert!(matches!(
        AxisymmetricBridgeLoadLine::new(&source, 0.001, limited)
            .unwrap()
            .trial(previous, 1.02 * h, &gate),
        Err(AxisymmetricLoadLineError::Budget)
    ));

    let mut force_accuracy = controls();
    force_accuracy.force_absolute_tolerance_n = 1e-30;
    force_accuracy.force_relative_tolerance = 1e-30;
    assert!(matches!(
        AxisymmetricBridgeLoadLine::new(&source, 0.0, force_accuracy)
            .unwrap()
            .trial(previous, h, &gate),
        Err(AxisymmetricLoadLineError::ForceAccuracy { .. })
    ));

    let mut gap_accuracy = controls();
    gap_accuracy.mechanical.gap_absolute_tolerance_m = 1e-30;
    gap_accuracy.mechanical.gap_relative_tolerance = 1e-30;
    assert!(matches!(
        AxisymmetricBridgeLoadLine::new(&source, 0.0, gap_accuracy)
            .unwrap()
            .trial(previous, h, &gate),
        Err(AxisymmetricLoadLineError::CompatibilityAccuracy { .. })
    ));

    let derivative = previous.profile().force_gap_derivative_n_per_m;
    assert!(derivative > 0.0);
    let unstable = AxisymmetricBridgeLoadLine::new(&source, 4.0 / derivative, controls()).unwrap();
    assert!(matches!(
        unstable.trial(previous, h, &gate),
        Err(AxisymmetricLoadLineError::CompatibilitySlope { .. })
    ));

    let mut short = controls();
    short.mechanical.maximum_gap_increment_m = 1e-10;
    assert!(matches!(
        AxisymmetricBridgeLoadLine::new(&source, 0.0, short)
            .unwrap()
            .trial(previous, 1.02 * h, &gate),
        Err(AxisymmetricLoadLineError::GapIncrement)
    ));
    assert_eq!(solver.trial(previous, 1.02 * h, &gate).unwrap(), expected);
    assert_eq!(source.evaluate(previous).unwrap(), *previous.profile());
}
