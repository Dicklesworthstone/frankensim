use super::*;

fn options() -> PoleResponseOptions {
    PoleResponseOptions {
        spectral_scale: 7.0,
        pole_step_tolerance: 1e-12,
        pole_residual_tolerance: 1e-11,
        linear_residual_tolerance: 1e-11,
        max_dimension: 16,
        max_vectors: 3,
        max_dense_bytes: 1_000_000,
        max_work: 1_000_000,
    }
}

fn prepare(generator: &[f64], scale: f64, s: C64, seed: &[C64]) -> PreparedPole {
    let n = seed.len();
    let mut matrix: Vec<_> = generator.iter().map(|a| C64::from_re(-a / scale)).collect();
    let mut derivative = vec![C64::ZERO; n * n];
    for i in 0..n {
        matrix[i * n + i] = matrix[i * n + i] + C64::new(s.re / scale, s.im / scale);
        derivative[i * n + i] = C64::from_re(1.0 / scale);
    }
    PreparedPole::new(
        s,
        seed,
        &PoleEvaluation {
            matrix,
            derivative,
            residual_scale: 32.0 / scale,
        },
        options(),
        || false,
    )
    .unwrap()
}

fn near(actual: C64, expected: C64) {
    assert!(
        (actual - expected).abs() < 2e-11,
        "{actual:?} != {expected:?}"
    );
}

// G1: an independent rational pressure response checks the native adapter's
// forcing scale, bilinear pressure observation and nonnormal excitation factor.
#[test]
fn original_ports_retain_complex_damped_participation() {
    let (omega, decay) = (7.0_f64, 2.0_f64);
    let generator = [0.0, omega, -omega, -2.0 * decay];
    let s = C64::new(-decay, -det::sqrt(omega * omega - decay * decay));
    let seed = [C64::ONE, s.scale(1.0 / omega)];
    let prepared = prepare(&generator, 100.0, s, &seed);
    let ports = [(PortId(4), vec![2.0, 0.0]), (PortId(8), vec![-5.0, 0.0])];
    let (residues, residuals) = observe_ports(&prepared, &ports, 100.0, &mut || false).unwrap();
    // B_i B_j * (s+2*decay)/(s²+2*decay*s+omega²).
    let base = (s + C64::from_re(2.0 * decay)) / (s + C64::from_re(decay)).scale(2.0);
    for (&got, weight) in residues.iter().zip([4.0, -10.0, -10.0, 25.0]) {
        near(got, base.scale(weight));
    }
    assert!(base.im.abs() > 0.1); // A squared mode amplitude cannot supply this phase.
    assert!(residuals.into_iter().all(|r| r < 1e-13));
}

#[test]
fn port_residue_preserves_pencil_scale_mode_gauge_and_dark_relaxation() {
    let generator = [0.0, 7.0, 0.0, -7.0, -4.0, 0.0, 0.0, 0.0, -3.0];
    let s = C64::new(-2.0, -det::sqrt(45.0));
    let seed = [C64::ONE, s.scale(1.0 / 7.0), C64::ZERO];
    let ports = [
        (PortId(4), vec![2.0, 0.0, 0.0]),
        (PortId(8), vec![-5.0, 0.0, 0.0]),
    ];
    let reference = observe_ports(
        &prepare(&generator, 7.0, s, &seed),
        &ports,
        7.0,
        &mut || false,
    )
    .unwrap()
    .0;
    for scale in [1e-100, 7.0, 1e100] {
        for gain in [C64::new(1e200, -3e200), C64::new(-1e-200, 2e-200)] {
            let changed = seed.map(|x| x * gain);
            let prepared = prepare(&generator, scale, s, &changed);
            let (residues, _) = observe_ports(&prepared, &ports, scale, &mut || false).unwrap();
            for (&got, &expected) in residues.iter().zip(&reference) {
                near(got, expected);
            }
        }
    }
    // Real internal relaxation remains a valid simple mode even though neither
    // port can excite or observe it; preserve exact zero without filtering it.
    let dark = prepare(
        &generator,
        7.0,
        C64::from_re(-3.0),
        &[C64::ZERO, C64::ZERO, C64::ONE],
    );
    let (residues, _) = observe_ports(&dark, &ports, 7.0, &mut || false).unwrap();
    assert_eq!(residues, vec![C64::ZERO; 4]);
    assert!(observe_ports(&dark, &ports, 7.0, &mut || true).is_err());
}

#[test]
fn complete_mapping_and_numerical_factor_are_admitted_together() {
    let limits = options();
    assert!(admit(3, 2, &limits).is_ok());
    assert!(admit(3, 4, &limits).is_err());
    assert!(admit(17, 2, &limits).is_err());
    assert!(admit(3, 0, &limits).is_err());
    assert!(
        admit(
            3,
            2,
            &PoleResponseOptions {
                max_work: 1,
                ..limits
            }
        )
        .is_err()
    );
    assert!(
        admit(
            3,
            2,
            &PoleResponseOptions {
                max_dense_bytes: 1,
                ..limits
            }
        )
        .is_err()
    );
    // Physical maps fit, but leave too little for the shared bordered factor.
    assert!(
        admit(
            3,
            2,
            &PoleResponseOptions {
                max_dense_bytes: 2000,
                ..limits
            }
        )
        .is_err()
    );
}

#[test]
fn forcing_underflow_cannot_turn_a_small_physical_port_into_a_dark_mode() {
    let s = C64::new(-2.0, -det::sqrt(45.0));
    let prepared = prepare(
        &[0.0, 7.0, -7.0, -4.0],
        1e308,
        s,
        &[C64::ONE, s.scale(1.0 / 7.0)],
    );
    let ports = [(PortId(4), vec![1e-20, 0.0])];
    // This pressure/volume residue is O(1e-40), representable in f64. The
    // intermediate scaled force is not; refuse instead of publishing zero.
    assert!(observe_ports(&prepared, &ports, 1e308, &mut || false).is_err());
}
