//! Physical force residues of original full quadratic pencils, without a
//! diagonal modal projection. References are analytic scalar/2x2 resolvents.
use fs_math::c64::C64;
use fs_modal::quadratic_modes::{
    PoleOptions, PoleResponseOptions, PreparedQuadraticModes, PreparedQuadraticResponse,
    QuadraticForceResidue, QuadraticMode, QuadraticObservation, QuadraticResponseError,
    QuadraticResponseOptions,
};

const I: C64 = C64::new(0.0, 1.0);
fn mode_options(rate: f64) -> PoleOptions {
    PoleOptions {
        spectral_scale: rate,
        step_tolerance: 2e-13,
        residual_tolerance: 2e-12,
        max_iterations: 40,
        max_dimension: 16,
        max_dense_bytes: 16 * 1024 * 1024,
        max_work: 100_000_000,
    }
}
fn response_options(rate: f64) -> QuadraticResponseOptions {
    QuadraticResponseOptions {
        numerical: PoleResponseOptions {
            spectral_scale: rate,
            pole_step_tolerance: 3e-12,
            pole_residual_tolerance: 3e-12,
            linear_residual_tolerance: 3e-12,
            max_dimension: 16,
            max_vectors: 4,
            max_dense_bytes: 16 * 1024 * 1024,
            max_work: 100_000_000,
        },
        max_observations: 4,
        max_physical_bytes: 16 * 1024 * 1024,
        max_physical_work: 100_000_000,
    }
}
fn close(actual: C64, expected: C64, tolerance: f64) {
    let bound = if expected == C64::ZERO {
        tolerance
    } else {
        tolerance * expected.abs()
    };
    assert!(
        (actual - expected).abs() <= bound,
        "{actual:?} != {expected:?}, bound {bound}"
    );
}
fn prepare(
    m: &[f64],
    c: &[f64],
    k: &[f64],
    s: C64,
    q: &[C64],
    rate: f64,
) -> PreparedQuadraticResponse {
    let owner =
        PreparedQuadraticModes::new(m, c, k, q.len(), mode_options(rate), || false).unwrap();
    let mode = owner.refine(s, q, || false).unwrap();
    owner
        .prepare_response(&mode, response_options(rate), || false)
        .unwrap()
}
fn check_mapping(p: &PreparedQuadraticResponse, r: &QuadraticForceResidue) {
    assert!(r.linear_residual < 3e-12);
    assert!(p.physical_residual() < 3e-12);
    for i in 0..r.displacement.len() {
        close(
            r.displacement[i],
            p.right_displacement()[i] * r.participation,
            2e-12,
        );
        close(
            r.velocity[i],
            p.laplace_pole_per_s() * r.displacement[i],
            2e-12,
        );
    }
}
fn pencil(m: &[f64], c: &[f64], k: &[f64], s: C64) -> [C64; 4] {
    core::array::from_fn(|i| s * s.scale(m[i]) + s.scale(c[i]) + C64::from_re(k[i]))
}
fn analytic_residue(m: &[f64], c: &[f64], k: &[f64], s: C64, f: &[C64]) -> [C64; 2] {
    let t = pencil(m, c, k, s);
    let d: [C64; 4] = core::array::from_fn(|i| s.scale(2.0 * m[i]) + C64::from_re(c[i]));
    let determinant_derivative = d[0] * t[3] + t[0] * d[3] - d[1] * t[2] - t[1] * d[2];
    [
        (t[3] * f[0] - t[1] * f[1]) / determinant_derivative,
        (-t[2] * f[0] + t[0] * f[1]) / determinant_derivative,
    ]
}
fn direct_response(m: &[f64], c: &[f64], k: &[f64], s: C64, f: &[C64]) -> [C64; 2] {
    let t = pencil(m, c, k, s);
    let determinant = t[0] * t[3] - t[1] * t[2];
    [
        (t[3] * f[0] - t[1] * f[1]) / determinant,
        (-t[2] * f[0] + t[0] * f[1]) / determinant,
    ]
}

#[test]
fn scalar_residues_have_the_original_force_mass_and_laplace_units() {
    let force = [C64::new(2.0, -0.75)];
    for c in [2.0, -2.0] {
        let s = C64::new(-c / 8.0, (9.0_f64 - (c / 8.0) * (c / 8.0)).sqrt());
        let prepared = prepare(&[4.0], &[c], &[36.0], s, &[C64::ONE], 3.0);
        let displacement = [C64::new(1.0, 2.0)];
        let velocity = [C64::new(0.5, -0.25)];
        let observation = QuadraticObservation {
            displacement: &displacement,
            velocity: &velocity,
        };
        let results = prepared
            .forcing_residues(&[&force], &[observation], || false)
            .unwrap();
        let r = &results[0];
        let expected = force[0] / (s.scale(8.0) + C64::from_re(c));
        close(r.displacement[0], expected, 5e-13);
        close(r.participation, expected.scale(2.0), 5e-13);
        close(
            r.observations[0],
            displacement[0] * expected + velocity[0] * s * expected,
            8e-13,
        );
        assert_eq!(r.laplace_pole_per_s.re > 0.0, c < 0.0);
        check_mapping(&prepared, r);
    }
    // Both simple overdamped roots are kept, including the opposite residue sign.
    for s in [-1.0, -2.0] {
        let p = prepare(&[1.0], &[3.0], &[2.0], C64::from_re(s), &[C64::ONE], 2.0);
        let r = p.forcing_residues(&[&force], &[], || false).unwrap();
        close(
            r[0].displacement[0],
            force[0].scale(1.0 / (2.0 * s + 3.0)),
            3e-13,
        );
    }
}

#[test]
fn nonnormal_force_with_zero_right_projection_still_excites_the_mode() {
    // Triangular K gives right mode [1,0], while its left vector sees force 2.
    // A phi^H*f shortcut would return a completely wrong zero response.
    let m = [1.0, 0.0, 0.0, 1.0];
    let c = [0.4, 0.0, 0.0, 0.8];
    let k = [4.0, 7.0, 0.0, 9.0];
    let s = C64::new(-0.2, 3.96_f64.sqrt());
    let p = prepare(&m, &c, &k, s, &[C64::ONE, C64::ZERO], 3.0);
    let forces = [
        [C64::ONE, C64::ZERO],
        [C64::ZERO, C64::ONE],
        [C64::ZERO, C64::ZERO],
    ];
    let zeros = [C64::ZERO; 2];
    let observation = QuadraticObservation {
        displacement: &forces[1],
        velocity: &zeros,
    };
    let r = p
        .forcing_residues(
            &[&forces[0], &forces[1], &forces[2]],
            &[observation],
            || false,
        )
        .unwrap();
    assert_eq!(p.right_displacement()[1], C64::ZERO);
    for i in 0..2 {
        let expected = analytic_residue(&m, &c, &k, s, &forces[i]);
        for j in 0..2 {
            close(r[i].displacement[j], expected[j], 5e-13);
        }
        assert_eq!(r[i].observations[0], C64::ZERO);
        check_mapping(&p, &r[i]);
    }
    assert!(r[1].displacement[0].abs() > 0.1);
    assert_eq!(r[2].displacement, zeros);
    assert_eq!(r[2].velocity, zeros);
    assert_eq!(r[2].participation, C64::ZERO);
}

#[test]
fn nonproportional_phase_matches_resolvent_limit_and_complex_sensor_weights() {
    let m = [1.0, 0.0, 0.0, 1.0];
    let c = [1.0, 0.5, 0.5, 3.0];
    let k = [5.0, 2.5, 2.5, 5.0];
    let s = C64::new(-1.0, 2.0);
    let p = prepare(&m, &c, &k, s, &[C64::ONE, I], 3.0);
    let f = [C64::ONE, C64::ZERO];
    let phased = f.map(|z| z * C64::new(-0.25, 0.75));
    let q_weights = [C64::new(1.0, 2.0), C64::from_re(0.5)];
    let v_weights = [C64::ZERO, -I];
    let o = QuadraticObservation {
        displacement: &q_weights,
        velocity: &v_weights,
    };
    let results = p.forcing_residues(&[&f, &phased], &[o], || false).unwrap();
    let expected = analytic_residue(&m, &c, &k, s, &f);
    for i in 0..2 {
        close(results[0].displacement[i], expected[i], 2e-12);
        close(
            results[1].displacement[i],
            expected[i] * C64::new(-0.25, 0.75),
            2e-12,
        );
    }
    close(
        results[0].observations[0],
        q_weights[0] * expected[0] + q_weights[1] * expected[1] + v_weights[1] * s * expected[1],
        2e-12,
    );
    // Hermitian projection is not the complex-symmetric left-mode contraction.
    let wrong = C64::from_re(0.5_f64.sqrt()) / (s.scale(2.0) + C64::from_re(2.0));
    assert!((results[0].participation - wrong).abs() > 0.2);
    let h = C64::new(3e-6, 4e-6);
    let direct = direct_response(&m, &c, &k, s + h, &f);
    for i in 0..2 {
        close(direct[i] * h, expected[i], 2e-5);
    }
    for r in &results {
        check_mapping(&p, r);
    }
    let again = p.forcing_residues(&[&f, &phased], &[o], || false).unwrap();
    assert_eq!(results[0].displacement, again[0].displacement);
    assert_eq!(results[1].participation, again[1].participation);
}

#[test]
fn gyroscopic_handedness_changes_participation_with_full_signed_coupling() {
    let m = [2.0, 0.0, 0.0, 2.0];
    let k = [18.0, 0.0, 0.0, 18.0];
    let force = [C64::new(1.0, 0.3), C64::new(-0.4, 0.2)];
    for g in [-4.0_f64, 4.0] {
        let c = [0.0, -g, g, 0.0];
        let s = I.scale((g + (g * g + 144.0).sqrt()) / 4.0);
        let p = prepare(&m, &c, &k, s, &[C64::ONE, I], 4.0);
        let r = p.forcing_residues(&[&force], &[], || false).unwrap();
        let expected = analytic_residue(&m, &c, &k, s, &force);
        for i in 0..2 {
            close(r[0].displacement[i], expected[i], 2e-12);
        }
        check_mapping(&p, &r[0]);
    }
}

#[test]
fn original_coordinate_and_extreme_time_scales_preserve_physical_residues() {
    let m = [1.0, 0.0, 0.0, 1.0];
    let c = [1.0, 0.5, 0.5, 3.0];
    let k = [5.0, 2.5, 2.5, 5.0];
    let s = C64::new(-1.0, 2.0);
    let f = [C64::new(1.0, 0.25), C64::new(-0.5, 0.1)];
    let expected = analytic_residue(&m, &c, &k, s, &f);
    let d = [2.0_f64.powi(300), -2.0_f64.powi(-300)];
    let map = |a: &[f64; 4]| -> [f64; 4] { core::array::from_fn(|i| a[i] * d[i / 2] * d[i % 2]) };
    let transformed_force = [f[0].scale(d[0]), f[1].scale(d[1])];
    let transformed_seed = [C64::from_re(1.0 / d[0]), I.scale(1.0 / d[1])];
    let p = prepare(&map(&m), &map(&c), &map(&k), s, &transformed_seed, 3.0);
    let zeros = [C64::ZERO; 2];
    let weights = [C64::from_re(d[0]), C64::from_re(2.0 * d[1])];
    let r = p
        .forcing_residues(
            &[&transformed_force],
            &[QuadraticObservation {
                displacement: &weights,
                velocity: &zeros,
            }],
            || false,
        )
        .unwrap();
    for i in 0..2 {
        close(r[0].displacement[i].scale(d[i]), expected[i], 2e-12);
    }
    close(
        r[0].observations[0],
        expected[0] + expected[1].scale(2.0),
        2e-12,
    );
    check_mapping(&p, &r[0]);
    // A full non-diagonal coordinate change adds mutual mass, not just units.
    // q=S*z, S=[[1,0.5],[0,2]], so T_z=S^T*T_q*S and F_z=S^T*F_q.
    let transform = |a: &[f64; 4]| -> [f64; 4] {
        [
            a[0],
            0.5 * a[0] + 2.0 * a[1],
            0.5 * a[0] + 2.0 * a[2],
            0.25 * a[0] + a[1] + a[2] + 4.0 * a[3],
        ]
    };
    let mf = transform(&m);
    let cf = transform(&c);
    let kf = transform(&k);
    let ff = [f[0], f[0].scale(0.5) + f[1].scale(2.0)];
    let p = prepare(
        &mf,
        &cf,
        &kf,
        s,
        &[C64::ONE - I.scale(0.25), I.scale(0.5)],
        3.0,
    );
    let r = p.forcing_residues(&[&ff], &[], || false).unwrap();
    close(
        r[0].displacement[0] + r[0].displacement[1].scale(0.5),
        expected[0],
        3e-12,
    );
    close(r[0].displacement[1].scale(2.0), expected[1], 3e-12);
    check_mapping(&p, &r[0]);
    for e in [-500, 500] {
        let a = 2.0_f64.powi(e);
        let b = 2.0_f64.powi(-e);
        // M'=b M, C'=b*a C=C, K'=b*a^2 K=a K: Rq'=Rq/(b*a)=Rq.
        let mt = m.map(|x| x * b);
        let kt = k.map(|x| x * a);
        let p = prepare(&mt, &c, &kt, s.scale(a), &[C64::ONE, I], 3.0 * a);
        let r = p.forcing_residues(&[&f], &[], || false).unwrap();
        for i in 0..2 {
            close(r[0].displacement[i], expected[i], 3e-12);
            close(r[0].velocity[i].scale(1.0 / a), s * expected[i], 3e-12);
        }
        check_mapping(&p, &r[0]);
    }
}

#[test]
fn finite_large_complex_force_and_zero_pole_are_not_erased_or_overflowed() {
    let mass = 2.0_f64.powi(1000);
    let p = prepare(
        &[mass],
        &[0.0],
        &[4.0 * mass],
        I.scale(2.0),
        &[C64::ONE],
        2.0,
    );
    let force = [C64::new(1e308, 1e308)]; // |F| overflows; each physical component is finite.
    let r = p.forcing_residues(&[&force], &[], || false).unwrap();
    close(
        r[0].displacement[0],
        force[0].scale(1.0 / mass) / I.scale(4.0),
        3e-13,
    );
    check_mapping(&p, &r[0]);
    let free = prepare(&[1.0], &[2.0], &[0.0], C64::ZERO, &[C64::ONE], 2.0);
    let r = free
        .forcing_residues(&[&[C64::from_re(3.0)]], &[], || false)
        .unwrap();
    close(r[0].displacement[0], C64::from_re(1.5), 3e-13);
    assert_eq!(r[0].velocity[0], C64::ZERO);
}

#[test]
fn stale_modes_repeated_roots_and_resource_or_cancellation_refusals_are_atomic() {
    let owner = PreparedQuadraticModes::new(&[1.0], &[1.0], &[4.0], 1, mode_options(2.0), || false)
        .unwrap();
    let mode = owner
        .refine(C64::new(-0.5, 3.75_f64.sqrt()), &[C64::ONE], || false)
        .unwrap();
    let p = owner
        .prepare_response(&mode, response_options(2.0), || false)
        .unwrap();
    let force = [C64::new(1.0, 0.2)];
    let before = p.forcing_residues(&[&force], &[], || false).unwrap();
    let mut calls = 0;
    assert!(matches!(
        p.forcing_residues(&[&force, &force], &[], || {
            calls += 1;
            calls > 8
        }),
        Err(QuadraticResponseError::Cancelled)
    ));
    let after = p.forcing_residues(&[&force], &[], || false).unwrap();
    assert_eq!(before[0].displacement, after[0].displacement);
    assert!(p.forcing_residues(&[], &[], || false).is_err());
    assert!(
        p.forcing_residues(
            &[&force],
            &[QuadraticObservation {
                displacement: &[],
                velocity: &force
            }],
            || false
        )
        .is_err()
    );
    assert!(
        p.forcing_residues(&[&[C64::new(f64::NAN, 0.0)]], &[], || false)
            .is_err()
    );
    let mut regauged = mode.clone();
    for q in &mut regauged.displacement {
        *q = *q * C64::new(3.0, 2.0);
    }
    regauged.velocity.fill(C64::ZERO);
    let rechecked = owner
        .prepare_response(&regauged, response_options(2.0), || false)
        .unwrap();
    let r = rechecked
        .forcing_residues(&[&force], &[], || false)
        .unwrap();
    close(r[0].participation, before[0].participation, 2e-13);
    close(
        rechecked.right_displacement()[0],
        p.right_displacement()[0],
        2e-13,
    );
    let mut stale = mode.clone();
    stale.laplace_frequency_per_s = stale.laplace_frequency_per_s + C64::from_re(0.1);
    stale.relative_residual = 0.0; // A supplied receipt does not override the equation.
    assert!(
        owner
            .prepare_response(&stale, response_options(2.0), || false)
            .is_err()
    );
    for kind in 0..4 {
        let mut o = response_options(2.0);
        match kind {
            0 => o.max_physical_bytes = 1,
            1 => o.max_physical_work = 1,
            2 => o.numerical.max_dense_bytes = 1,
            _ => o.numerical.max_vectors = 0,
        }
        assert!(owner.prepare_response(&mode, o, || false).is_err());
    }
    assert!(matches!(
        owner.prepare_response(&mode, response_options(2.0), || true),
        Err(QuadraticResponseError::Cancelled)
    ));
    let critical =
        PreparedQuadraticModes::new(&[1.0], &[4.0], &[4.0], 1, mode_options(2.0), || false)
            .unwrap();
    let repeated = QuadraticMode {
        laplace_frequency_per_s: C64::from_re(-2.0),
        displacement: vec![C64::ONE],
        velocity: vec![C64::from_re(-2.0)],
        relative_residual: 0.0,
        solver_relative_residual: 0.0,
        relative_correction: 0.0,
        iterations: 0,
    };
    assert!(
        critical
            .prepare_response(&repeated, response_options(2.0), || false)
            .is_err()
    );
}
