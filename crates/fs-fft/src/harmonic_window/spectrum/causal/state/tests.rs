use super::*;

fn limits(n: usize) -> HarmonicDecayIntervalLimits {
    HarmonicDecayIntervalLimits {
        maximum_tones: n,
        maximum_workspace_bytes: PreparedHarmonicDecayInterval::required_workspace_bytes(n).unwrap(),
        maximum_work: PreparedHarmonicDecayInterval::required_work(n).unwrap(),
    }
}

fn prepare(omega: &[f64], start: Dd, h: f64, tau: f64) -> PreparedHarmonicDecayInterval {
    PreparedHarmonicDecayInterval::new(omega, start, h, tau, limits(omega.len())).unwrap()
}

fn pair(re: f64, im: f64) -> Pair { [Dd::from_f64(re), Dd::from_f64(im)] }

fn relative(actual: Dd, expected: f64, tolerance: f64) {
    let error = (actual.to_f64() - expected).abs();
    assert!(error <= tolerance * expected.abs().max(1e-300),
        "actual={actual:?}, expected={expected:e}, error={error:e}");
}

// Independent integration of the original ODE and its two observables.
fn rk4(omega: &[f64], drive: &[Pair], work: &[Pair], t0: f64, h: f64, tau: f64, q0: f64) -> [f64; 3] {
    let waveform = |field: &[Pair], t: f64| -> f64 {
        omega.iter().zip(field).map(|(&w, a)| {
            a[0].to_f64() * (w * t).cos() + a[1].to_f64() * (w * t).sin()
        }).sum()
    };
    let rhs = |t: f64, y: [f64; 3]| {
        [waveform(drive, t) - y[0] / tau, y[0] * y[0], y[0] * waveform(work, t)]
    };
    let n = 8192;
    let dt = h / n as f64;
    let mut y = [q0, 0.0, 0.0];
    for step in 0..n {
        let t = t0 + step as f64 * dt;
        let a = rhs(t, y);
        let b = rhs(t + dt / 2.0, core::array::from_fn(|k| y[k] + dt * a[k] / 2.0));
        let c = rhs(t + dt / 2.0, core::array::from_fn(|k| y[k] + dt * b[k] / 2.0));
        let d = rhs(t + dt, core::array::from_fn(|k| y[k] + dt * c[k]));
        for k in 0..3 { y[k] += dt * (a[k] + 2.0 * b[k] + 2.0 * c[k] + d[k]) / 6.0; }
    }
    y
}

#[test]
fn complete_causal_state_and_independent_work_match_original_ode() {
    let omega = [0.0, 1.7, 2.3, -0.8];
    let drive = [pair(0.2, 0.0), pair(0.7, -0.4), pair(-0.2, 0.3), pair(0.1, 0.5)];
    let work = [pair(-0.1, 0.0), pair(0.8, 0.2), pair(0.6, -0.1), pair(-0.3, 0.4)];
    for (h, tau) in [(0.1, 0.7), (1.4, 0.7), (0.2, 0.003)] {
        let kernel = prepare(&omega, Dd::from_f64(0.37), h, tau);
        let m = kernel.evaluate(Dd::from_f64(-0.6), &drive, &work).unwrap();
        let reference = rk4(&omega, &drive, &work, 0.37, h, tau, -0.6);
        relative(m.end, reference[0], 3e-10);
        relative(m.square_integral().unwrap(), reference[1], 3e-10);
        relative(m.mixed_integral().unwrap(), reference[2], 3e-10);
        assert!(m.mean_square.hi > 0.0);
        assert!(kernel.numerical_payload_bytes()
            <= PreparedHarmonicDecayInterval::required_workspace_bytes(omega.len()).unwrap());
    }
}

#[test]
fn cosine_peak_startup_retains_quadratic_state_and_quintic_heat() {
    let h = 1e-10;
    let omega = 2.0;
    let strain = 1.3;
    let tau = 0.7;
    let drive = [pair(0.0, -omega * strain)];
    let m = prepare(&[omega], Dd::ZERO, h, tau).evaluate(Dd::ZERO, &drive, &drive).unwrap();
    relative(m.end, -0.5 * strain * omega * omega * h * h, 2e-9);
    relative(m.increment, m.end.to_f64(), 1e-28);
    relative(m.mean_square, strain * strain * omega.powi(4) * h.powi(4) / 20.0, 2e-9);
    relative(m.square_integral().unwrap() / Dd::from_f64(tau),
        strain * strain * omega.powi(4) * h.powi(5) / (20.0 * tau), 2e-9);
    relative(m.mean_mixed, strain * strain * omega.powi(4) * h.powi(3) / 8.0, 2e-9);
}

#[test]
fn retained_initial_memory_releases_energy_below_endpoint_resolution() {
    let h = 1e-20;
    let tau = 0.4;
    let m = prepare(&[], Dd::ZERO, h, tau).evaluate(Dd::ONE, &[], &[]).unwrap();
    let alpha = Dd::from_f64(h) / Dd::from_f64(tau);
    assert_eq!(m.end.hi, 1.0);
    assert!(m.end.lo < 0.0);
    assert!((m.increment + alpha).abs().hi < 1e-38);
    assert!(m.square_integral().unwrap().hi > 0.0);
    let storage_change = m.increment + Dd::from_f64(0.5) * m.increment * m.increment;
    let heat = m.mean_square * alpha;
    assert!((storage_change + heat).abs().hi < 1e-48);
    assert!(is_zero(m.mean_mixed));
}

#[test]
fn subdivision_and_original_energy_identity_retain_the_source_clock() {
    let omega = [0.73, 1.29];
    let drive = [pair(0.6, -0.2), pair(-0.3, 0.7)];
    let start = Dd::from_pair(2.0_f64.powi(40), 0.125);
    let h = 0.0625;
    let tau = 0.4;
    let q0 = Dd::from_f64(-0.2);
    let full = prepare(&omega, start, 2.0 * h, tau).evaluate(q0, &drive, &drive).unwrap();
    let a = prepare(&omega, start, h, tau).evaluate(q0, &drive, &drive).unwrap();
    let b = prepare(&omega, start + Dd::from_f64(h), h, tau).evaluate(a.end, &drive, &drive).unwrap();
    assert!((full.end - b.end).abs().hi < 1e-15);
    assert!((full.square_integral().unwrap() - a.square_integral().unwrap()
        - b.square_integral().unwrap()).abs().hi < 1e-16);
    assert!((full.mixed_integral().unwrap() - a.mixed_integral().unwrap()
        - b.mixed_integral().unwrap()).abs().hi < 1e-16);
    let storage = q0 * full.increment + Dd::from_f64(0.5) * full.increment * full.increment;
    let heat = full.square_integral().unwrap() / Dd::from_f64(tau);
    assert!((full.mixed_integral().unwrap() - storage - heat).abs().hi < 1e-15);
}

#[test]
fn endpoint_phase_change_is_not_a_rounded_absolute_clock_subtraction() {
    let omega = 0.7;
    let h = 2.0_f64.powi(-30);
    let start = Dd::from_pair(2.0_f64.powi(45), 0.0003);
    assert_eq!(start.hi + h, start.hi);
    let phase = phase_angular_interval_from_clock_dd(omega, start, h).unwrap();
    let local = multiply(phase.increment, [phase.start[0], -phase.start[1]]).unwrap();
    let beta = Dd::from_f64(omega) * Dd::from_f64(h);
    relative(local[0], (-Dd::from_f64(0.5) * beta * beta).to_f64(), 2e-14);
    relative(local[1], (-beta).to_f64(), 2e-14);
    assert!(!is_zero(phase.increment[0]) || !is_zero(phase.increment[1]));
}

#[test]
fn coherent_equal_tones_are_combined_before_quadratic_moments() {
    let large = Dd::from_pair(1e16, 1.0);
    let drive = [[large, Dd::ZERO], pair(-1e16, 0.0)];
    let work = [pair(0.2, 0.3), pair(0.1, -0.5)];
    let actual = prepare(&[2.0, 2.0], Dd::from_f64(0.1), 0.8, 0.7)
        .evaluate(Dd::from_f64(0.2), &drive, &work).unwrap();
    let reference = prepare(&[2.0], Dd::from_f64(0.1), 0.8, 0.7)
        .evaluate(Dd::from_f64(0.2), &[pair(1.0, 0.0)], &[pair(0.3, -0.2)]).unwrap();
    assert!((actual.mean_square - reference.mean_square).abs().hi < 1e-28);
    assert!((actual.mean_mixed - reference.mean_mixed).abs().hi < 1e-16);
    assert!((actual.end - reference.end).abs().hi < 1e-28);
}

#[test]
fn finite_admission_cancellation_and_unresolved_energy_refuse_explicitly() {
    let mut small = limits(1);
    small.maximum_work = 1;
    assert!(matches!(PreparedHarmonicDecayInterval::new(&[1.0], Dd::ZERO, 1.0, 1.0, small),
        Err(ProductSpectrumError::Budget)));
    assert!(PreparedHarmonicDecayInterval::required_workspace_bytes(65).is_err());
    assert!(matches!(PreparedHarmonicDecayInterval::new_cancellable(&[1.0], Dd::ZERO, 1.0, 1.0,
        limits(1), &mut || true), Err(ProductSpectrumError::Cancelled)));
    let k = prepare(&[1.0], Dd::ZERO, 0.1, 1.0);
    assert!(k.evaluate(Dd::ZERO, &[], &[]).is_err());
    assert!(k.evaluate(Dd::from_f64(f64::NAN), &[pair(0.0, 0.0)], &[pair(0.0, 0.0)]).is_err());
    assert!(matches!(k.evaluate_cancellable(Dd::ONE, &[pair(0.0, 0.0)], &[pair(0.0, 0.0)],
        &mut || true), Err(ProductSpectrumError::Cancelled)));
    assert!(prepare(&[], Dd::ZERO, 0.1, 1.0)
        .evaluate(Dd::from_f64(1e-200), &[], &[]).is_err());
    let nearly_null = prepare(&[1.0, 1.0 + 1e-15], Dd::ZERO, 10.0, 1.0);
    assert!(nearly_null.evaluate(Dd::ZERO, &[pair(1.0, 0.0), pair(-1.0, 0.0)],
        &[pair(0.0, 0.0); 2]).is_err());
    let zero = prepare(&[1.0], Dd::from_f64(1e20), 0.0, 1.0)
        .evaluate(Dd::ONE, &[pair(1.0, 0.0)], &[pair(1.0, 0.0)]).unwrap();
    assert_eq!(zero.end, Dd::ONE);
    assert!(is_zero(zero.increment) && is_zero(zero.mean_square) && is_zero(zero.mean_mixed));
}

#[test]
fn stiff_decay_does_not_round_the_factor_before_a_large_initial_state() {
    let q0: f64 = 1e150;
    let m = prepare(&[], Dd::ZERO, 740.0, 1.0)
        .evaluate(Dd::from_f64(q0), &[], &[]).unwrap();
    relative(m.end, (q0.ln() - 740.0).exp(), 2e-13);
    relative(m.mean_square, q0 * (q0 / 1480.0), 2e-15);
    assert!(m.end.hi > 0.0 && m.mean_square.hi > 0.0);
}

#[test]
fn conjugate_opposite_frequencies_share_one_real_carrier() {
    let fields = [pair(0.7, -0.2), pair(-0.7, -0.2)];
    for h in [0.1, 2.0] {
        let moment = prepare(&[1.2, -1.2], Dd::from_f64(0.3), h, 0.7)
            .evaluate(Dd::ZERO, &fields, &fields).unwrap();
        assert!(is_zero(moment.end) && is_zero(moment.increment));
        assert!(is_zero(moment.mean_square) && is_zero(moment.mean_mixed));
    }
}
