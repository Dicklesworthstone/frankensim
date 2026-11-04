//! G1/G3 physical complex mode shapes: original polynomials, phase and units.

use fs_math::c64::C64;
use fs_modal::quadratic_modes::{
    PoleOptions, PreparedQuadraticModes, QuadraticMode, QuadraticModeError,
};

const I: C64 = C64::new(0.0, 1.0);

fn options(rate: f64) -> PoleOptions {
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

fn close(actual: C64, expected: C64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * expected.abs().max(1e-280),
        "{actual:?} != {expected:?}"
    );
}

fn mass_norm(mode: &QuadraticMode, mass: &[f64]) -> C64 {
    let n = mode.displacement.len();
    let mut result = C64::ZERO;
    for i in 0..n {
        for j in 0..n {
            result =
                result + mode.displacement[i].conj().scale(mass[i * n + j]) * mode.displacement[j];
        }
    }
    result
}

fn check_original(mode: &QuadraticMode, m: &[f64], c: &[f64], k: &[f64]) {
    let n = mode.displacement.len();
    let s = mode.laplace_frequency_per_s;
    for i in 0..n {
        let mut row = C64::ZERO;
        let mut bound = 0.0;
        for j in 0..n {
            let q = mode.displacement[j];
            let term =
                s * s.scale(m[i * n + j]) + s.scale(c[i * n + j]) + C64::from_re(k[i * n + j]);
            row = row + term * q;
            bound += (s.abs().powi(2) * m[i * n + j].abs()
                + s.abs() * c[i * n + j].abs()
                + k[i * n + j].abs())
                * q.abs();
        }
        assert!(
            row.abs() <= 2e-11 * bound.max(1e-280),
            "row {i}: {row:?}, bound {bound}"
        );
        close(mode.velocity[i], s * mode.displacement[i], 2e-15);
    }
    close(mass_norm(mode, m), C64::ONE, 8e-14);
    assert!(mode.relative_residual <= 2e-12);
}

#[test]
fn scalar_mode_has_physical_mass_and_velocity_including_active_growth() {
    for c in [2.0, -2.0] {
        let decay = -c / 8.0;
        let root = C64::new(decay, (9.0_f64 - decay * decay).sqrt());
        let damping = [c];
        let model =
            PreparedQuadraticModes::new(&[4.0], &damping, &[36.0], 1, options(3.0), || false)
                .unwrap();
        let mode = model
            .refine(root + C64::new(0.02, 0.03), &[C64::new(7.0, -2.0)], || {
                false
            })
            .unwrap();
        close(mode.laplace_frequency_per_s, root, 3e-13);
        close(mode.displacement[0], C64::from_re(0.5), 2e-15);
        check_original(&mode, &[4.0], &[c], &[36.0]);
        assert_eq!(mode.laplace_frequency_per_s.re.is_sign_positive(), c < 0.0);
    }
}

#[test]
fn simple_rigid_pole_and_unstable_stiffness_keep_their_physical_real_shapes() {
    let free =
        PreparedQuadraticModes::new(&[1.0], &[2.0], &[0.0], 1, options(2.0), || false).unwrap();
    let rigid = free
        .refine(C64::from_re(0.05), &[C64::ONE], || false)
        .unwrap();
    assert!(rigid.laplace_frequency_per_s.abs() < 1e-12);
    assert!(rigid.velocity[0].abs() < 1e-12);
    close(rigid.displacement[0], C64::ONE, 2e-15);
    assert!(rigid.relative_residual < 1e-12);

    let unstable =
        PreparedQuadraticModes::new(&[1.0], &[0.4], &[-4.0], 1, options(2.0), || false).unwrap();
    let expected = (-0.4 + 16.16_f64.sqrt()) / 2.0;
    let growth = unstable
        .refine(C64::from_re(expected * 1.03), &[C64::from_re(-5.0)], || {
            false
        })
        .unwrap();
    close(
        growth.laplace_frequency_per_s,
        C64::from_re(expected),
        3e-12,
    );
    close(growth.displacement[0], C64::ONE, 2e-15);
    check_original(&growth, &[1.0], &[0.4], &[-4.0]);
}

#[test]
fn nonproportional_damping_retains_relative_complex_motion() {
    // The original determinant is (s^2+2s+5)(s^2+2s+15/4).
    // At s=-1+2i, the ORIGINAL equations have phi proportional to [1,i].
    // C and K do not commute; no real modal basis diagonalizes both.
    let m = [1.0, 0.0, 0.0, 1.0];
    let c = [1.0, 0.5, 0.5, 3.0];
    let k = [5.0, 2.5, 2.5, 5.0];
    let model = PreparedQuadraticModes::new(&m, &c, &k, 2, options(3.0), || false).unwrap();
    let seed = [C64::ONE, C64::new(0.02, 0.9)];
    let mode = model.refine(C64::new(-0.9, 2.08), &seed, || false).unwrap();
    close(mode.laplace_frequency_per_s, C64::new(-1.0, 2.0), 3e-12);
    close(mode.displacement[1] / mode.displacement[0], I, 5e-12);
    check_original(&mode, &m, &c, &k);
    let replay = model.refine(C64::new(-0.9, 2.08), &seed, || false).unwrap();
    assert_eq!(mode.laplace_frequency_per_s, replay.laplace_frequency_per_s);
    assert_eq!(mode.displacement, replay.displacement);
    assert_eq!(mode.velocity, replay.velocity);
}

#[test]
fn gyro_sign_selects_the_physical_circular_polarization_branch() {
    let m = [2.0, 0.0, 0.0, 2.0];
    let k = [18.0, 0.0, 0.0, 18.0];
    for g in [-4.0_f64, 4.0] {
        // G*[1,i] = -i*g*[1,i], hence -2*w^2 + g*w + 18 = 0.
        let frequency = (g + (g * g + 144.0).sqrt()) / 4.0;
        let c = [0.0, -g, g, 0.0];
        let model = PreparedQuadraticModes::new(&m, &c, &k, 2, options(5.0), || false).unwrap();
        let mode = model
            .refine(C64::new(0.01, frequency * 1.01), &[C64::ONE, I], || false)
            .unwrap();
        close(
            mode.laplace_frequency_per_s,
            C64::new(0.0, frequency),
            2e-12,
        );
        close(mode.displacement[1] / mode.displacement[0], I, 3e-12);
        check_original(&mode, &m, &c, &k);
        let mut power = C64::ZERO;
        for i in 0..2 {
            for j in 0..2 {
                power = power + mode.velocity[i].conj() * mode.velocity[j].scale(c[i * 2 + j]);
            }
        }
        assert!(power.re.abs() < 1e-13); // a gyroscope adds no loss
    }
}

#[test]
fn congruent_mixed_units_and_extreme_rate_scales_preserve_normalized_shapes() {
    let mass = [1.0, 0.0, 0.0, 1.0];
    let damping = [1.0, 0.5, 0.5, 3.0];
    let stiffness = [5.0, 2.5, 2.5, 5.0];
    let two = 2.0_f64;
    for (coordinate, coefficient, rate) in [
        ([-two.powi(300), two.powi(-300)], 1.0, 1.0),
        ([1.0, 1.0], two.powi(-500), two.powi(500)),
        ([1.0, 1.0], two.powi(500), two.powi(-500)),
    ] {
        // x=D*z. Congruence scales every ORIGINAL coefficient, while a
        // time-unit change multiplies C by rate and K by rate squared.
        let mut transformed = [[0.0; 4]; 3];
        for (order, original) in [mass, damping, stiffness].iter().enumerate() {
            for i in 0..2 {
                for j in 0..2 {
                    let magnitude = match order {
                        0 => coefficient,
                        1 => coefficient * rate,
                        _ => (coefficient * rate) * rate,
                    };
                    transformed[order][i * 2 + j] =
                        original[i * 2 + j] * coordinate[i] * coordinate[j] * magnitude;
                }
            }
        }
        let model = PreparedQuadraticModes::new(
            &transformed[0],
            &transformed[1],
            &transformed[2],
            2,
            options(3.0 * rate),
            || false,
        )
        .unwrap();
        let seed = [
            C64::from_re(1.0 / coordinate[0]),
            I.scale(1.0 / coordinate[1]),
        ];
        let mode = model
            .refine(C64::new(-0.93, 2.07).scale(rate), &seed, || false)
            .unwrap();
        close(
            mode.laplace_frequency_per_s.scale(1.0 / rate),
            C64::new(-1.0, 2.0),
            3e-12,
        );
        let mapped = [
            mode.displacement[0].scale(coordinate[0] * coefficient.sqrt()),
            mode.displacement[1].scale(coordinate[1] * coefficient.sqrt()),
        ];
        close(mapped[1] / mapped[0], I, 8e-12);
        assert!((mapped[0].norm_sq() + mapped[1].norm_sq() - 1.0).abs() < 8e-14);
        close(mass_norm(&mode, &transformed[0]), C64::ONE, 8e-14);
        assert!(mode.relative_residual < 2e-12);
        for i in 0..2 {
            close(
                mode.velocity[i],
                mode.laplace_frequency_per_s * mode.displacement[i],
                3e-15,
            );
        }
    }
}

#[test]
fn mutual_inertia_is_not_replaced_by_diagonal_modal_mass() {
    // Congruence x=S*z, S=[[1,1/2],[0,2]], of the analytic nonproportional case.
    let s = [[1.0, 0.5], [0.0, 2.0]];
    let mut matrices = [[0.0; 4]; 3];
    for (order, a) in [
        [1.0, 0.0, 0.0, 1.0],
        [1.0, 0.5, 0.5, 3.0],
        [5.0, 2.5, 2.5, 5.0],
    ]
    .iter()
    .enumerate()
    {
        for i in 0..2 {
            for j in 0..2 {
                for r in 0..2 {
                    for col in 0..2 {
                        matrices[order][i * 2 + j] += s[r][i] * a[r * 2 + col] * s[col][j];
                    }
                }
            }
        }
    }
    let model = PreparedQuadraticModes::new(
        &matrices[0],
        &matrices[1],
        &matrices[2],
        2,
        options(3.0),
        || false,
    )
    .unwrap();
    let mode = model
        .refine(
            C64::new(-0.95, 2.02),
            &[C64::new(1.0, -0.25), C64::new(0.0, 0.5)],
            || false,
        )
        .unwrap();
    close(mode.laplace_frequency_per_s, C64::new(-1.0, 2.0), 3e-12);
    let mapped = [
        mode.displacement[0] + mode.displacement[1].scale(0.5),
        mode.displacement[1].scale(2.0),
    ];
    close(mapped[1] / mapped[0], I, 4e-12);
    check_original(&mode, &matrices[0], &matrices[1], &matrices[2]);
}

#[test]
fn repeated_defective_invalid_budget_and_cancelled_pencils_refuse() {
    let m = [1.0, 0.0, 0.0, 1.0];
    let c = [0.0; 4];
    let k = [4.0, 0.0, 0.0, 4.0];
    let model = PreparedQuadraticModes::new(&m, &c, &k, 2, options(2.0), || false).unwrap();
    assert!(matches!(
        model.refine(C64::new(0.0, 2.0), &[C64::ONE, C64::ZERO], || false),
        Err(QuadraticModeError::Solve(_))
    ));
    let defective =
        PreparedQuadraticModes::new(&[1.0], &[4.0], &[4.0], 1, options(2.0), || false).unwrap();
    assert!(matches!(
        defective.refine(C64::from_re(-2.0), &[C64::ONE], || false),
        Err(QuadraticModeError::Solve(_))
    ));
    let mut small = options(2.0);
    small.max_dense_bytes = 1;
    assert!(matches!(
        PreparedQuadraticModes::new(&m, &c, &k, 2, small, || false),
        Err(QuadraticModeError::Budget { .. })
    ));
    small = options(2.0);
    small.max_work = 1;
    assert!(matches!(
        PreparedQuadraticModes::new(&m, &c, &k, 2, small, || false),
        Err(QuadraticModeError::Budget { .. })
    ));
    assert!(matches!(
        PreparedQuadraticModes::new(&m, &c, &k, 2, options(2.0), || true),
        Err(QuadraticModeError::Cancelled)
    ));
    let mut polls = 0;
    assert!(matches!(
        model.refine(C64::new(0.01, 2.01), &[C64::ONE, I], || {
            polls += 1;
            polls > 5
        }),
        Err(QuadraticModeError::Cancelled)
    ));
    assert!(matches!(
        PreparedQuadraticModes::new(&[1.0, 0.0, 1.0, 1.0], &c, &k, 2, options(2.0), || false),
        Err(QuadraticModeError::MassNotPositiveDefinite)
    ));
    assert!(matches!(
        model.refine(I, &[C64::ZERO; 2], || false),
        Err(QuadraticModeError::Invalid(_))
    ));
}
