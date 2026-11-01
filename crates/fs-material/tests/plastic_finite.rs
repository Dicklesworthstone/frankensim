//! Finite multiplicative kinematics, independent return-map references,
//! work-conjugate potential and all 81 consistent tangent entries.
use fs_material::plastic::finite::{
    FiniteJ2, FiniteJ2Error, FiniteJ2Limits, FiniteJ2Parameters, FiniteJ2State,
};

const IDENTITY: [f64; 9] = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0];

fn parameters() -> FiniteJ2Parameters {
    FiniteJ2Parameters {
        shear_modulus_pa: 80e9,
        bulk_modulus_pa: 160e9,
        yield_stress_pa: 250e6,
        hardening_modulus_pa: 2e9,
    }
}

fn limits() -> FiniteJ2Limits {
    FiniteJ2Limits {
        minimum_jacobian: 0.05,
        maximum_elastic_log_strain: 0.3,
        maximum_plastic_log_strain: 3.0,
        maximum_accumulated_plastic_strain: 3.0,
    }
}

fn law() -> FiniteJ2 {
    FiniteJ2::new(parameters(), limits()).unwrap()
}

fn near(actual: f64, expected: f64, absolute: f64, relative: f64) {
    assert!(
        (actual - expected).abs() <= absolute + relative * actual.abs().max(expected.abs()),
        "{actual:.16e} != {expected:.16e}; abs difference {:.8e}",
        (actual - expected).abs()
    );
}

fn coaxial(shear: f64, dilation: f64) -> [f64; 9] {
    [
        (dilation + shear).exp(),
        0.0,
        0.0,
        0.0,
        (dilation - 0.5 * shear).exp(),
        0.0,
        0.0,
        0.0,
        (dilation - 0.5 * shear).exp(),
    ]
}

fn multiply(a: &[f64; 9], b: &[f64; 9]) -> [f64; 9] {
    core::array::from_fn(|ij| {
        (0..3)
            .map(|k| a[3 * (ij / 3) + k] * b[3 * k + ij % 3])
            .sum()
    })
}

fn transpose(a: &[f64; 9]) -> [f64; 9] {
    core::array::from_fn(|i| a[3 * (i % 3) + i / 3])
}

fn determinant(a: &[f64; 9]) -> f64 {
    a[0] * (a[4] * a[8] - a[5] * a[7]) - a[1] * (a[3] * a[8] - a[5] * a[6])
        + a[2] * (a[3] * a[7] - a[4] * a[6])
}

fn metric(s: &FiniteJ2State) -> [f64; 9] {
    let c = s.inverse_plastic_cauchy_green;
    [c[0], c[3], c[5], c[3], c[1], c[4], c[5], c[4], c[2]]
}

fn rotation(angle: f64) -> [f64; 9] {
    let (s, c) = angle.sin_cos();
    [c, -s, 0.0, s, c, 0.0, 0.0, 0.0, 1.0]
}

#[test]
fn independent_coaxial_reference_fixes_flow_stress_energy_and_kirchhoff_convention() {
    for hardening in [0.0, 2e9] {
        let p = FiniteJ2Parameters {
            hardening_modulus_pa: hardening,
            ..parameters()
        };
        let material = FiniteJ2::new(p, limits()).unwrap();
        let (shear, dilation) = (0.25, 0.04);
        let f = coaxial(shear, dilation);
        let old = material.initial_state();
        let r = material.response(&f, &old).unwrap();
        let dg = (3.0 * p.shear_modulus_pa * shear - p.yield_stress_pa)
            / (3.0 * p.shear_modulus_pa + hardening);
        let ee = shear - dg;
        let j = (3.0 * dilation).exp();
        let q = p.yield_stress_pa + hardening * dg;
        near(r.plastic_increment, dg, 3e-16, 0.0);
        near(r.state.alpha, dg, 3e-16, 0.0);
        near(r.jacobian, j, 3e-16, 2e-16);
        near(r.equivalent_kirchhoff_pa, q, 1e-5, 2e-14);
        let cp = [(-2.0 * dg).exp(), dg.exp(), dg.exp(), 0.0, 0.0, 0.0];
        for (actual, expected) in r.state.inverse_plastic_cauchy_green.iter().zip(cp) {
            near(*actual, expected, 3e-15, 0.0);
        }
        let expected_tau = [
            p.bulk_modulus_pa * 3.0 * dilation + 2.0 * p.shear_modulus_pa * ee,
            p.bulk_modulus_pa * 3.0 * dilation - p.shear_modulus_pa * ee,
            p.bulk_modulus_pa * 3.0 * dilation - p.shear_modulus_pa * ee,
        ];
        for (i, expected) in [0, 4, 8].into_iter().zip(expected_tau) {
            near(r.kirchhoff_pa[i], expected, 2e-4, 2e-14);
            near(r.first_piola_pa[i], expected / f[i], 2e-4, 2e-14);
        }
        let stored = 0.5 * p.bulk_modulus_pa * (3.0 * dilation).powi(2)
            + 1.5 * p.shear_modulus_pa * ee * ee
            + 0.5 * hardening * dg * dg;
        near(r.stored_energy_j_m3, stored, 1e-5, 3e-14);
        near(
            r.plastic_dissipation_j_m3,
            p.yield_stress_pa * dg,
            1e-7,
            2e-14,
        );
        near(
            r.incremental_potential_j_m3,
            stored + p.yield_stress_pa * dg,
            1e-5,
            3e-14,
        );
        near(
            material.stored_energy(&f, &r.state).unwrap(),
            stored,
            1e-5,
            5e-14,
        );
        // Cauchy equivalent stress is q/J; it is not the configured yield radius.
        assert!((q / j - q).abs() > 1e7);
        material.validate_committed(&f, &r.state, 1e-3).unwrap();
        assert_eq!(old, material.initial_state());
    }
}

#[test]
fn rigid_rotations_and_noncoaxial_retained_history_are_objective() {
    let material = law();
    let initial = material.initial_state();
    let q = rotation(0.73);
    let pure_rotation = material.response(&q, &initial).unwrap();
    assert_eq!(pure_rotation.state, initial);
    assert_eq!(pure_rotation.plastic_increment, 0.0);
    assert!(pure_rotation.first_piola_pa.iter().all(|v| v.abs() < 1e-3));

    let first = [1.0, 0.22, 0.1, 0.0, 1.0, -0.12, 0.0, 0.0, 1.0];
    let second = [1.12, -0.1, 0.04, 0.03, 0.91, 0.16, -0.03, 0.02, 1.01];
    let retained = material.response(&first, &initial).unwrap().state;
    assert!(
        retained.inverse_plastic_cauchy_green[3..]
            .iter()
            .any(|x| x.abs() > 0.01)
    );
    // A hold or superposed rotation of an accepted plastic endpoint must not
    // accumulate fictitious flow from an eigensolver-sized yield residual.
    for held_f in [first, multiply(&q, &first)] {
        let held = material.response(&held_f, &retained).unwrap();
        assert_eq!(held.state, retained);
        assert_eq!(held.plastic_increment, 0.0);
        assert_eq!(held.plastic_dissipation_j_m3, 0.0);
    }
    for old in [&initial, &retained] {
        let a = material.response(&second, old).unwrap();
        let b = material.response(&multiply(&q, &second), old).unwrap();
        let qp = multiply(&q, &a.first_piola_pa);
        let qtauqt = multiply(&multiply(&q, &a.kirchhoff_pa), &transpose(&q));
        for i in 0..9 {
            near(b.first_piola_pa[i], qp[i], 2e-3, 2e-12);
            near(b.kirchhoff_pa[i], qtauqt[i], 2e-3, 2e-12);
        }
        for i in 0..6 {
            near(
                b.state.inverse_plastic_cauchy_green[i],
                a.state.inverse_plastic_cauchy_green[i],
                2e-14,
                2e-14,
            );
        }
        near(b.state.alpha, a.state.alpha, 2e-15, 2e-14);
        near(
            b.incremental_potential_j_m3,
            a.incremental_potential_j_m3,
            1e-5,
            2e-13,
        );
        near(
            b.plastic_dissipation_j_m3,
            a.plastic_dissipation_j_m3,
            1e-5,
            2e-13,
        );
    }
}

fn check_derivatives(material: &FiniteJ2, f: [f64; 9], old: &FiniteJ2State) {
    let original = old.clone();
    let r = material.response(&f, old).unwrap();
    assert_eq!(
        material.stress_and_potential(&f, old).unwrap(),
        (r.first_piola_pa, r.incremental_potential_j_m3)
    );
    let mut squared_error = 0.0;
    let mut squared_reference = 0.0;
    for column in 0..9 {
        let h = 1e-6;
        let mut plus = f;
        let mut minus = f;
        plus[column] += h;
        minus[column] -= h;
        let a = material.stress_and_potential(&plus, old).unwrap();
        let b = material.stress_and_potential(&minus, old).unwrap();
        let potential_derivative = (a.1 - b.1) / (2.0 * h);
        near(potential_derivative, r.first_piola_pa[column], 25.0, 4e-7);
        for row in 0..9 {
            let numerical = (a.0[row] - b.0[row]) / (2.0 * h);
            let analytic = r.tangent_pa[row][column];
            near(analytic, numerical, 5e3, 2e-6);
            near(analytic, r.tangent_pa[column][row], 1e-3, 3e-12);
            squared_error += (analytic - numerical).powi(2);
            squared_reference += numerical * numerical;
        }
    }
    assert!((squared_error / squared_reference).sqrt() < 2e-7);
    assert_eq!(old, &original);
    assert_eq!(material.response(&f, old).unwrap(), r);
}

#[test]
fn analytic_full_tangent_and_potential_gradient_cover_repeated_roots_and_plastic_history() {
    let material = law();
    let initial = material.initial_state();
    check_derivatives(&material, IDENTITY, &initial);
    check_derivatives(&material, coaxial(0.0004, 0.0002), &initial);
    // Two exactly repeated trial eigenvalues, active plastic return.
    check_derivatives(&material, coaxial(0.2, 0.015), &initial);
    let mut nearly_repeated = coaxial(0.2, 0.015);
    nearly_repeated[4] *= 1.0 + 2e-11;
    nearly_repeated[8] *= 1.0 - 2e-11;
    check_derivatives(&material, nearly_repeated, &initial);
    let old = material
        .response(&coaxial(0.16, 0.0), &initial)
        .unwrap()
        .state;
    check_derivatives(
        &material,
        [1.18, 0.16, 0.08, -0.04, 0.9, 0.04, 0.02, -0.02, 0.98],
        &old,
    );
    let perfect = FiniteJ2::new(
        FiniteJ2Parameters {
            hardening_modulus_pa: 0.0,
            ..parameters()
        },
        limits(),
    )
    .unwrap();
    check_derivatives(&perfect, coaxial(0.2, 0.0), &perfect.initial_state());
}

#[test]
fn unloading_and_reversal_preserve_complete_history_and_isochoric_plastic_flow() {
    let material = law();
    let initial = material.initial_state();
    let loaded = material.response(&coaxial(0.25, 0.0), &initial).unwrap();
    let relaxed_f = coaxial(loaded.state.alpha, 0.0);
    let relaxed = material.response(&relaxed_f, &loaded.state).unwrap();
    assert_eq!(relaxed.state, loaded.state);
    assert_eq!(relaxed.plastic_dissipation_j_m3, 0.0);
    assert!(relaxed.first_piola_pa.iter().all(|x| x.abs() < 1e-3));
    near(
        relaxed.stored_energy_j_m3,
        0.5 * parameters().hardening_modulus_pa * relaxed.state.alpha.powi(2),
        1e-6,
        1e-13,
    );

    let reverse_f = coaxial(-0.2, 0.0);
    let reverse = material.response(&reverse_f, &loaded.state).unwrap();
    let p = parameters();
    let q_trial = 3.0 * p.shear_modulus_pa * (0.2 + loaded.state.alpha);
    let dg = (q_trial - p.yield_stress_pa - p.hardening_modulus_pa * loaded.state.alpha)
        / (3.0 * p.shear_modulus_pa + p.hardening_modulus_pa);
    near(reverse.plastic_increment, dg, 2e-15, 1e-14);
    near(reverse.state.alpha, loaded.state.alpha + dg, 2e-15, 1e-14);
    assert!(reverse.kirchhoff_pa[0] < 0.0);
    assert!(reverse.state.alpha > loaded.state.alpha);
    near(
        reverse.state.inverse_plastic_cauchy_green[0],
        (-2.0 * (loaded.state.alpha - dg)).exp(),
        2e-14,
        2e-14,
    );

    let mut state = reverse.state;
    for f in [
        [1.0, 0.14, -0.07, 0.0, 1.0, 0.08, 0.0, 0.0, 1.0],
        coaxial(0.13, 0.02),
        [1.0, -0.2, 0.04, 0.0, 1.0, -0.11, 0.0, 0.0, 1.0],
    ] {
        let r = material.response(&f, &state).unwrap();
        assert!(r.state.alpha >= state.alpha);
        assert!(r.plastic_dissipation_j_m3 >= 0.0);
        near(determinant(&metric(&r.state)), 1.0, 3e-13, 0.0);
        material.validate_committed(&f, &r.state, 1e-3).unwrap();
        near(
            material.stored_energy(&f, &r.state).unwrap(),
            r.stored_energy_j_m3,
            1e-4,
            3e-12,
        );
        state = r.state;
    }
}

#[test]
fn hydrostatic_stress_has_no_plastic_flow_and_elastic_domain_applies_after_return() {
    let material = FiniteJ2::new(
        parameters(),
        FiniteJ2Limits {
            maximum_elastic_log_strain: 0.01,
            ..limits()
        },
    )
    .unwrap();
    let old = material.initial_state();
    // The large trial elastic strain is outside the domain, but the returned
    // elastic strain is small and the finite plastic stretch remains admitted.
    let plastic = material.response(&coaxial(0.3, 0.0), &old).unwrap();
    assert!(plastic.plastic_increment > 0.29);
    let hydro = material.response(&coaxial(0.0, 0.003), &old).unwrap();
    assert_eq!(hydro.state, old);
    assert_eq!(hydro.plastic_increment, 0.0);
    for i in [0, 4, 8] {
        near(
            hydro.kirchhoff_pa[i],
            parameters().bulk_modulus_pa * 0.009,
            1e-4,
            2e-13,
        );
    }
    assert!(matches!(
        material.response(&coaxial(0.0, 0.03), &old),
        Err(FiniteJ2Error::Domain { .. })
    ));
}

#[test]
fn malformed_histories_bounds_and_unaccepted_endpoints_are_rejected_without_mutation() {
    let material = law();
    let initial = material.initial_state();
    let mut invalid = initial.clone();
    invalid.inverse_plastic_cauchy_green[0] = 2.0;
    assert!(material.validate_history(&invalid).is_err());
    invalid.alpha = 0.1;
    assert!(material.validate_history(&invalid).is_err()); // det Cp != 1
    invalid.inverse_plastic_cauchy_green = [
        (-0.4_f64).exp(),
        0.2_f64.exp(),
        0.2_f64.exp(),
        0.0,
        0.0,
        0.0,
    ];
    assert!(material.validate_history(&invalid).is_err()); // alpha insufficient
    invalid.alpha = 0.3;
    material.validate_history(&invalid).unwrap();
    invalid.inverse_plastic_cauchy_green[0] = -1.0;
    assert!(material.validate_history(&invalid).is_err());
    invalid = initial.clone();
    invalid.alpha = f64::NAN;
    assert!(material.validate_history(&invalid).is_err());
    let mut inverted = IDENTITY;
    inverted[0] = -1.0;
    assert!(material.response(&inverted, &initial).is_err());
    assert!(material.stress_and_potential(&inverted, &initial).is_err());
    assert!(material.response(&[f64::NAN; 9], &initial).is_err());
    assert!(
        material
            .validate_committed(&coaxial(0.1, 0.0), &initial, 1e-3)
            .is_err()
    );
    assert_eq!(initial, material.initial_state());

    let accumulated_bound = FiniteJ2::new(
        parameters(),
        FiniteJ2Limits {
            maximum_accumulated_plastic_strain: 0.4,
            ..limits()
        },
    )
    .unwrap();
    let a = accumulated_bound
        .response(&coaxial(0.1, 0.0), &initial)
        .unwrap()
        .state;
    let b = accumulated_bound
        .response(&coaxial(-0.1, 0.0), &a)
        .unwrap()
        .state;
    let checkpoint = b.clone();
    assert!(matches!(
        accumulated_bound.response(&coaxial(0.1, 0.0), &b),
        Err(FiniteJ2Error::Domain {
            quantity: "accumulated plastic strain",
            ..
        })
    ));
    assert_eq!(b, checkpoint);
    let metric_bound = FiniteJ2::new(
        parameters(),
        FiniteJ2Limits {
            maximum_plastic_log_strain: 0.05,
            ..limits()
        },
    )
    .unwrap();
    assert!(metric_bound.response(&coaxial(0.1, 0.0), &initial).is_err());
    assert!(
        metric_bound
            .stress_and_potential(&coaxial(0.1, 0.0), &initial)
            .is_err()
    );
    for bad in [f64::NAN, f64::INFINITY, -1.0, 0.0] {
        assert!(
            FiniteJ2::new(
                FiniteJ2Parameters {
                    shear_modulus_pa: bad,
                    ..parameters()
                },
                limits()
            )
            .is_err()
        );
    }
    assert!(
        FiniteJ2::new(
            parameters(),
            FiniteJ2Limits {
                minimum_jacobian: 1.0,
                ..limits()
            }
        )
        .is_err()
    );
}
