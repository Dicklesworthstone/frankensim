//! G1/G3: full body compliance against independent matrix/KKT references.
use fs_tribo::layer_contact::ThirdBodyLayerContact;
use fs_tribo::tangential::matrix::{
    MatrixTangentialConfig, MatrixTangentialError, advance_with_matrix_compliance,
};
use fs_tribo::{
    FrictionLaw, FrictionRegime, TangentialLaw, TangentialSolveConfig, TangentialStepInput,
};

fn input(free: [f64; 2]) -> TangentialStepInput {
    TangentialStepInput {
        normal_n: 1.0,
        slip_increment_m: free,
        dt_s: 0.02,
        contact_temperature_k: Some(300.0),
        nominal_pressure_pa: None,
        compliance_m_per_n: 0.0,
        solve: TangentialSolveConfig {
            max_iterations: 128,
            maximum_force_n: 100.0,
            force_absolute_tolerance_n: 1e-13,
            force_relative_tolerance: 1e-12,
        },
    }
}

fn config() -> MatrixTangentialConfig {
    MatrixTangentialConfig {
        maximum_evaluations: 512,
        displacement_tolerance_m: 2e-13,
        energy_tolerance_j: 2e-12,
        matrix_relative_tolerance: 128.0 * f64::EPSILON,
    }
}

fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 3e-10 * actual.abs().max(expected.abs()).max(1.0),
        "actual {actual:.17e}, expected {expected:.17e}"
    );
}
fn close2(actual: [f64; 2], expected: [f64; 2]) {
    for i in 0..2 {
        close(actual[i], expected[i]);
    }
}
fn norm(x: [f64; 2]) -> f64 {
    x[0].hypot(x[1])
}
fn dot(a: [f64; 2], b: [f64; 2]) -> f64 {
    a[0] * b[0] + a[1] * b[1]
}
fn mul(a: [[f64; 2]; 2], x: [f64; 2]) -> [f64; 2] {
    [dot(a[0], x), dot(a[1], x)]
}
fn solve(a: [[f64; 2]; 2], b: [f64; 2]) -> [f64; 2] {
    let determinant = a[0][0] * a[1][1] - a[0][1] * a[1][0];
    [
        (a[1][1] * b[0] - a[0][1] * b[1]) / determinant,
        (a[0][0] * b[1] - a[1][0] * b[0]) / determinant,
    ]
}
fn shifted(mut a: [[f64; 2]; 2], shift: f64) -> [[f64; 2]; 2] {
    a[0][0] += shift;
    a[1][1] += shift;
    a
}

// Independent Coulomb KKT: (C + alpha I) F = u_free and |F| = mu_k N.
// A scalar multiplier root is unrelated to the production resolvent iteration.
fn coulomb_reference(
    c: [[f64; 2]; 2],
    free: [f64; 2],
    static_limit: f64,
    kinetic_limit: f64,
) -> ([f64; 2], FrictionRegime) {
    let static_force = solve(c, free);
    if static_force.iter().all(|v| v.is_finite()) && norm(static_force) <= static_limit {
        return (static_force, FrictionRegime::Sticking);
    }
    let mut lo = 0.0;
    let mut hi = norm(free) / kinetic_limit;
    for _ in 0..200 {
        let mid = 0.5 * lo + 0.5 * hi;
        if norm(solve(shifted(c, mid), free)) > kinetic_limit {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    (
        solve(shifted(c, 0.5 * lo + 0.5 * hi), free),
        FrictionRegime::Sliding,
    )
}

#[test]
fn isotropic_and_zero_matrices_preserve_the_original_native_candidate() {
    let law = FrictionLaw::Coulomb {
        static_mu: 0.6,
        kinetic_mu: 0.4,
    };
    let layer = ThirdBodyLayerContact::new(400.0, 0.8, 0.7, 300.0, 1.0).unwrap();
    for c in [0.0, 0.003] {
        for free in [[0.0, 0.0], [0.002, -0.003], [-0.04, 0.02]] {
            let original = input(free);
            let scalar = TangentialStepInput {
                compliance_m_per_n: c,
                ..original
            };
            let expected = law.advance(&law.initial_state(), &scalar).unwrap();
            let actual = advance_with_matrix_compliance(
                &law,
                &law.initial_state(),
                &original,
                [[c, 0.0], [0.0, c]],
                config(),
                || false,
            )
            .unwrap();
            assert_eq!(actual.contact, expected);
            assert_eq!(actual.law_evaluations, 1);
            assert_eq!(actual.iterations, 0);
            let expected = layer.advance(&layer.initial_state(), &scalar).unwrap();
            let actual = advance_with_matrix_compliance(
                &layer,
                &layer.initial_state(),
                &original,
                [[c, 0.0], [0.0, c]],
                config(),
                || false,
            )
            .unwrap();
            assert_eq!(actual.contact, expected);
            assert_eq!(actual.law_evaluations, 1);
        }
    }
}

#[test]
fn full_matrix_static_reaction_matches_the_unconstrained_matrix_solve() {
    let law = FrictionLaw::Coulomb {
        static_mu: 0.8,
        kinetic_mu: 0.5,
    };
    for c in [
        [[0.002, 0.0006], [0.0006, 0.001]],
        [[10.0, 0.0], [0.0, 0.001]],
    ] {
        let free = [0.0003, -0.0002];
        let expected = solve(c, free);
        assert!(norm(expected) < 0.8);
        let actual = advance_with_matrix_compliance(
            &law,
            &law.initial_state(),
            &input(free),
            c,
            config(),
            || false,
        )
        .unwrap();
        close2(actual.contact.traction_n, expected);
        assert_eq!(actual.contact.regime, FrictionRegime::Sticking);
        assert_eq!(actual.contact.slip_increment_m, [0.0; 2]);
        assert_eq!(actual.contact.dissipated_j, 0.0);
        assert!(norm(actual.contact.compatibility_residual_m) <= config().displacement_tolerance_m);
    }
    // A valid kinetic solution also exists here, but the original native
    // convention gives a feasible static reaction priority. Starting from a
    // smallest-eigenvalue scalar probe alone would wrongly publish sliding.
    let law = FrictionLaw::Coulomb {
        static_mu: 1.0,
        kinetic_mu: 0.5,
    };
    let c = [[0.1, 0.0], [0.0, 0.5]];
    let free = input([0.0, 0.4]);
    let actual =
        advance_with_matrix_compliance(&law, &law.initial_state(), &free, c, config(), || false)
            .unwrap();
    assert_eq!(actual.contact.regime, FrictionRegime::Sticking);
    close2(actual.contact.traction_n, [0.0, 0.8]);
    assert_eq!(actual.contact.slip_increment_m, [0.0; 2]);
    assert_eq!(actual.contact.dissipated_j, 0.0);
    // A force budget cannot turn a refused feasible static reaction into slip.
    let mut limited = free;
    limited.solve.maximum_force_n = 0.7;
    assert!(matches!(
        advance_with_matrix_compliance(&law, &law.initial_state(), &limited, c, config(), || false),
        Err(MatrixTangentialError::Law(_))
    ));
}

#[test]
fn unequal_static_kinetic_coefficients_slide_on_the_full_anisotropic_load_line() {
    let law = FrictionLaw::Coulomb {
        static_mu: 0.35,
        kinetic_mu: 0.2,
    };
    let c = [[0.1, 0.08], [0.08, 0.1]];
    for free in [[0.016, 0.004], [-0.016, -0.004], [0.004, 0.016]] {
        let (expected, regime) = coulomb_reference(c, free, 0.35, 0.2);
        assert_eq!(regime, FrictionRegime::Sliding);
        let actual = advance_with_matrix_compliance(
            &law,
            &law.initial_state(),
            &input(free),
            c,
            config(),
            || false,
        )
        .unwrap();
        close2(actual.contact.traction_n, expected);
        assert_eq!(actual.contact.regime, regime);
        let response = mul(c, expected);
        let slip = [free[0] - response[0], free[1] - response[1]];
        close2(actual.contact.slip_increment_m, slip);
        close(actual.contact.dissipated_j, 0.2 * norm(slip));
        assert!(actual.body_work_residual_j.abs() <= config().energy_tolerance_j);
        // The anisotropic body response rotates traction relative to free motion.
        assert!((expected[0] * free[1] - expected[1] * free[0]).abs() > 1e-4);
    }
}

#[test]
fn retained_elastic_preload_uses_both_diagonal_and_cross_axis_body_response() {
    let law = ThirdBodyLayerContact::new(400.0, 100.0, 0.6, 300.0, 1.0).unwrap();
    let loaded = law
        .advance(&law.initial_state(), &input([0.08, -0.04]))
        .unwrap()
        .state;
    let original = loaded.clone();
    let c = [[0.003, 0.001], [0.001, 0.007]];
    let actual =
        advance_with_matrix_compliance(&law, &loaded, &input([0.0; 2]), c, config(), || false)
            .unwrap();
    let a = [
        [1.0 + 200.0 * c[0][0], 200.0 * c[0][1]],
        [200.0 * c[1][0], 1.0 + 200.0 * c[1][1]],
    ];
    let expected = solve(a, [32.0, -16.0]);
    close2(actual.contact.traction_n, expected);
    assert_eq!(actual.contact.dissipated_j, 0.0);
    assert!(actual.contact.stored_change_j < 0.0);
    close(
        actual.contact.stored_change_j,
        -dot(expected, mul(c, expected)),
    );
    assert_eq!(loaded, original);
}

// Independent layer return: solve one plastic-flow multiplier after eliminating
// the full matrix elastic equation; never calls the production adapter or law.
fn layer_reference(
    c: [[f64; 2]; 2],
    free: [f64; 2],
    e: [f64; 2],
    k: f64,
    yield_force: f64,
    resistance: f64,
    dt: f64,
) -> ([f64; 2], [f64; 2], [f64; 2], f64) {
    let half = 0.5 * k;
    let a = [
        [1.0 + half * c[0][0], half * c[0][1]],
        [half * c[1][0], 1.0 + half * c[1][1]],
    ];
    let rhs = [k * e[0] + half * free[0], k * e[1] + half * free[1]];
    let elastic = solve(a, rhs);
    let flow = dt / resistance;
    let b = half * flow;
    let force = if norm(elastic) <= yield_force {
        elastic
    } else {
        let mut lo = 0.0;
        let mut hi = b;
        for _ in 0..200 {
            let mid = 0.5 * lo + 0.5 * hi;
            let magnitude = norm(solve(shifted(a, mid), rhs));
            let residual = mid - b * (1.0 - yield_force / magnitude).max(0.0);
            if residual < 0.0 {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        solve(shifted(a, 0.5 * lo + 0.5 * hi), rhs)
    };
    let magnitude = norm(force);
    let plastic_factor = if magnitude > yield_force {
        flow * (1.0 - yield_force / magnitude)
    } else {
        0.0
    };
    let plastic = [plastic_factor * force[0], plastic_factor * force[1]];
    let response = mul(c, force);
    let elastic = [
        e[0] + free[0] - response[0] - plastic[0],
        e[1] + free[1] - response[1] - plastic[1],
    ];
    let distance = norm(plastic);
    (
        force,
        elastic,
        plastic,
        yield_force * distance + resistance * distance * distance / dt,
    )
}

#[test]
fn nonproportional_layer_history_and_heat_match_an_independent_multiplier_solve() {
    let law = ThirdBodyLayerContact::new(200.0, 1.0, 0.8, 300.0, 1.0).unwrap();
    let c = [[0.01, 0.004], [0.004, 0.003]];
    let mut state = law.initial_state();
    let mut elastic = [0.0; 2];
    let mut plastic = [0.0; 2];
    let mut heat = 0.0;
    for i in 0..48 {
        let angle = 0.37 * f64::from(i);
        let free = [0.025 * angle.cos(), 0.018 * angle.sin()];
        let expected = layer_reference(c, free, elastic, 200.0, 1.0, 0.8, 0.02);
        let actual =
            advance_with_matrix_compliance(&law, &state, &input(free), c, config(), || false)
                .unwrap();
        close2(actual.contact.traction_n, expected.0);
        close2(actual.contact.state.elastic_slip_m(), expected.1);
        for axis in 0..2 {
            plastic[axis] += expected.2[axis];
        }
        close2(actual.contact.state.plastic_slip_m(), plastic);
        close(actual.contact.dissipated_j, expected.3);
        assert!(actual.body_work_residual_j.abs() <= config().energy_tolerance_j);
        heat += actual.contact.dissipated_j;
        elastic = expected.1;
        state = actual.contact.state;
    }
    assert!(heat > 0.1);
    assert!(state.plastic_slip_m().iter().all(|p| p.abs() > 1e-3));
}

fn rotate(v: [f64; 2], angle: f64) -> [f64; 2] {
    let (s, c) = angle.sin_cos();
    [c * v[0] - s * v[1], s * v[0] + c * v[1]]
}
fn rotate_matrix(matrix: [[f64; 2]; 2], angle: f64) -> [[f64; 2]; 2] {
    let first = rotate(mul(matrix, rotate([1.0, 0.0], -angle)), angle);
    let second = rotate(mul(matrix, rotate([0.0, 1.0], -angle)), angle);
    [[first[0], second[0]], [first[1], second[1]]]
}

#[test]
fn rotating_tangent_axes_transports_force_retained_strain_and_heat_together() {
    let law = ThirdBodyLayerContact::new(200.0, 1.0, 0.8, 300.0, 1.0).unwrap();
    let c = [[0.01, 0.004], [0.004, 0.003]];
    let angle = 0.713;
    let rotated = rotate_matrix(c, angle);
    let mut a = law.initial_state();
    let mut b = law.initial_state();
    for free in [[0.02, 0.03], [-0.01, 0.005], [0.0, -0.04], [0.03, 0.0]] {
        let left =
            advance_with_matrix_compliance(&law, &a, &input(free), c, config(), || false).unwrap();
        let right = advance_with_matrix_compliance(
            &law,
            &b,
            &input(rotate(free, angle)),
            rotated,
            config(),
            || false,
        )
        .unwrap();
        close2(
            right.contact.traction_n,
            rotate(left.contact.traction_n, angle),
        );
        close2(
            right.contact.state.elastic_slip_m(),
            rotate(left.contact.state.elastic_slip_m(), angle),
        );
        close2(
            right.contact.state.plastic_slip_m(),
            rotate(left.contact.state.plastic_slip_m(), angle),
        );
        close(right.contact.dissipated_j, left.contact.dissipated_j);
        close(right.contact.stored_change_j, left.contact.stored_change_j);
        a = left.contact.state;
        b = right.contact.state;
    }
}

#[test]
fn rank_one_and_frictionless_limits_preserve_unconstrained_motion() {
    let law = FrictionLaw::Coulomb {
        static_mu: 0.6,
        kinetic_mu: 0.4,
    };
    let c = [[0.01, 0.0], [0.0, 0.0]];
    let still = advance_with_matrix_compliance(
        &law,
        &law.initial_state(),
        &input([0.001, 0.0]),
        c,
        config(),
        || false,
    )
    .unwrap();
    close2(still.contact.traction_n, [0.1, 0.0]);
    assert_eq!(still.contact.regime, FrictionRegime::Sticking);
    let free = [0.001, 0.0001];
    let mut bounded = input(free);
    bounded.solve.maximum_force_n = 3.0;
    let actual =
        advance_with_matrix_compliance(&law, &law.initial_state(), &bounded, c, config(), || false)
            .unwrap();
    let (expected, _) = coulomb_reference(c, free, 0.6, 0.4);
    close2(actual.contact.traction_n, expected);
    close(actual.contact.slip_increment_m[1], free[1]);
    let zero = FrictionLaw::Coulomb {
        static_mu: 0.0,
        kinetic_mu: 0.0,
    };
    let actual = advance_with_matrix_compliance(
        &zero,
        &zero.initial_state(),
        &input(free),
        c,
        config(),
        || false,
    )
    .unwrap();
    assert_eq!(actual.contact.traction_n, [0.0; 2]);
    assert_eq!(actual.contact.slip_increment_m, free);
    assert_eq!(actual.contact.dissipated_j, 0.0);
    let layer = ThirdBodyLayerContact::new(100.0, 100.0, 1.0, 300.0, 1.0).unwrap();
    let loaded = layer
        .advance(&layer.initial_state(), &input([0.005, -0.004]))
        .unwrap()
        .state;
    let held =
        advance_with_matrix_compliance(&layer, &loaded, &input([0.0; 2]), c, config(), || false)
            .unwrap();
    close(held.contact.traction_n[1], -0.4);
    close(held.contact.state.elastic_slip_m()[1], -0.004);
}

#[test]
fn invalid_matrix_budget_source_or_cancellation_never_changes_accepted_history() {
    let law = ThirdBodyLayerContact::new(200.0, 1.0, 0.8, 300.0, 1.0).unwrap();
    let state = law
        .advance(&law.initial_state(), &input([0.02, -0.01]))
        .unwrap()
        .state;
    let before = state.clone();
    let free = input([0.01, 0.02]);
    for matrix in [
        [[-0.1, 0.0], [0.0, 1.0]],
        [[0.1, 0.2], [0.2, 0.1]],
        [[1.0, 0.1], [0.2, 1.0]],
        [[f64::NAN, 0.0], [0.0, 1.0]],
    ] {
        assert!(matches!(
            advance_with_matrix_compliance(&law, &state, &free, matrix, config(), || false),
            Err(MatrixTangentialError::Invalid(_))
        ));
    }
    let c = [[0.01, 0.004], [0.004, 0.003]];
    let budget = MatrixTangentialConfig {
        maximum_evaluations: 1,
        ..config()
    };
    assert!(matches!(
        advance_with_matrix_compliance(&law, &state, &free, c, budget, || false),
        Err(MatrixTangentialError::Unresolved { evaluations: 1, .. })
    ));
    let mut polls = 0;
    assert!(matches!(
        advance_with_matrix_compliance(&law, &state, &free, c, config(), || {
            polls += 1;
            polls == 2
        }),
        Err(MatrixTangentialError::Cancelled)
    ));
    assert_eq!(polls, 2);
    let wrong_temperature = TangentialStepInput {
        contact_temperature_k: Some(301.0),
        ..free
    };
    assert!(matches!(
        advance_with_matrix_compliance(&law, &state, &wrong_temperature, c, config(), || false),
        Err(MatrixTangentialError::Law(_))
    ));
    let double_compliance = TangentialStepInput {
        compliance_m_per_n: 0.01,
        ..free
    };
    assert!(matches!(
        advance_with_matrix_compliance(&law, &state, &double_compliance, c, config(), || false),
        Err(MatrixTangentialError::Invalid(_))
    ));
    assert_eq!(state, before);
    assert_eq!(
        advance_with_matrix_compliance(&law, &state, &free, c, config(), || false).unwrap(),
        advance_with_matrix_compliance(&law, &before, &free, c, config(), || false).unwrap()
    );
}
