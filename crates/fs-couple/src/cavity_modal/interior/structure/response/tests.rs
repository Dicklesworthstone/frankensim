use super::*;

fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() < 3e-14 * expected.abs().max(1.0),
        "{actual} != {expected}"
    );
}

fn real_energy(a: &[f64], x: &[f64; 2]) -> f64 {
    0.5 * (a[0] * x[0] * x[0] + (a[1] + a[2]) * x[0] * x[1] + a[3] * x[1] * x[1])
}

fn block_pair(a: [f64; 4], b: [f64; 4]) -> [f64; 16] {
    [
        a[0], a[1], 0.0, 0.0, a[2], a[3], 0.0, 0.0, 0.0, 0.0, b[0], b[1], 0.0, 0.0, b[2], b[3],
    ]
}

#[test]
fn original_mass_and_compliance_storage_matches_physical_cycle_average() {
    // G1: original non-diagonal structural mass and fluid compliance are
    // different metrics. No normalization or pressure-potential energy swap.
    let ms = [2.0, 0.3, 0.3, 1.5];
    let ks = [7.0, -1.0, -1.0, 4.0];
    let ma = [0.7, 0.15, 0.15, 1.1];
    let ka = [3.0, -3.0, -3.0, 3.0];
    let u = [C64::new(0.4, 0.3), C64::new(-0.2, 0.8)];
    let psi = [C64::new(0.2, -0.1), C64::new(-0.4, 0.5)];
    let v = u.map(|z| C64::new(2.3 * z.im, -2.3 * z.re));
    let p = psi.map(|z| C64::new(2.3 * z.im, -2.3 * z.re));
    let result = energy(
        &block_pair(ms, ma),
        &block_pair(ks, ka),
        &u,
        &v,
        &psi,
        &p,
        &mut || false,
    )
    .unwrap();
    let mut expected = [0.0; 4];
    for (cos, sin) in [(1.0, 0.0), (0.0, 1.0), (-1.0, 0.0), (0.0, -1.0)] {
        // Four exact phase samples integrate these real quadratic energies.
        let real = |q: [C64; 2]| q.map(|z| z.re * cos + z.im * sin);
        expected[0] += 0.25 * real_energy(&ms, &real(v));
        expected[1] += 0.25 * real_energy(&ks, &real(u));
        expected[2] += 0.25 * real_energy(&ka, &real(psi));
        expected[3] += 0.25 * real_energy(&ma, &real(p));
    }
    for (actual, expected) in [
        result.structural_kinetic_j,
        result.structural_strain_j,
        result.fluid_kinetic_j,
        result.fluid_compression_j,
    ]
    .into_iter()
    .zip(expected)
    {
        close(actual, expected);
    }
    close(result.total_j(), expected.iter().sum());
    assert!(result.fluid_kinetic_j != result.fluid_compression_j);
    assert!(result.structural_kinetic_j != result.structural_strain_j);
}

#[test]
fn reciprocal_original_interface_retains_weak_coherent_flow_and_work_slopes() {
    // G3: a weak surviving trace is retained beside exact large cancellation.
    let (force, flow) = interface(
        &[1e16, 1.0, -1e16],
        &[C64::from_re(1.0); 3],
        &[C64::from_re(1.0)],
        &mut || false,
    )
    .unwrap();
    assert_eq!(flow, [C64::from_re(1.0)]);
    close(work(&force, &[C64::from_re(1.0); 3]), 0.5);

    let b = [0.7, -0.2, 1.1, 0.8];
    let v = [C64::new(0.2, -0.7), C64::new(-0.4, 0.3)];
    let p = [C64::new(1.2, 0.6), C64::new(-0.8, 0.5)];
    let dv = [C64::new(-0.3, 0.2), C64::new(0.7, -0.1)];
    let dp = [C64::new(0.5, -0.2), C64::new(-0.6, -0.4)];
    let (force, flow) = interface(&b, &v, &p, &mut || false).unwrap();
    let (df, dq) = interface(&b, &dv, &dp, &mut || false).unwrap();
    close(work(&force, &v), work(&p, &flow));
    let mechanical = work_derivative(&force, Some(&df), &v, &dv);
    let acoustic = work_derivative(&p, Some(&dp), &flow, &dq);
    close(mechanical, acoustic);
    let h = 0.125;
    let varied = |sign: f64| {
        let pressure = [p[0] + dp[0].scale(sign * h), p[1] + dp[1].scale(sign * h)];
        let velocity = [v[0] + dv[0].scale(sign * h), v[1] + dv[1].scale(sign * h)];
        let (_, flow) = interface(&b, &velocity, &pressure, &mut || false).unwrap();
        work(&pressure, &flow)
    };
    // The linear phasor family gives an exact quadratic work polynomial.
    close(mechanical, (varied(1.0) - varied(-1.0)) / (2.0 * h));
}

#[test]
fn independent_source_and_heat_derivatives_exclude_gyroscopic_work() {
    let c = [2.0, 4.25, -3.75, 1.0]; // D=[[2,.25],[.25,1]], gyro=4.
    let v = [C64::new(0.6, -0.9), C64::new(0.4, 0.7)];
    let dv = [C64::new(0.2, 0.3), C64::new(-0.5, 0.6)];
    let f = [C64::new(0.7, -0.8), C64::new(-0.2, 0.4)];
    let df = [C64::new(-0.1, 0.2), C64::new(0.8, -0.3)];
    let heat = |q: [C64; 2]| {
        q[0].norm_sq() + 0.5 * q[1].norm_sq() + 0.25 * (q[0].re * q[1].re + q[0].im * q[1].im)
    };
    close(dissipation(&c, 2, 0, &v, &mut || false).unwrap(), heat(v));
    let h = 0.25;
    let varied = |q: [C64; 2], dq: [C64; 2], sign| {
        [q[0] + dq[0].scale(sign * h), q[1] + dq[1].scale(sign * h)]
    };
    close(
        dissipation_derivative(&c, 2, 0, &v, &dv, &mut || false).unwrap(),
        (heat(varied(v, dv, 1.0)) - heat(varied(v, dv, -1.0))) / (2.0 * h),
    );
    close(
        work_derivative(&f, Some(&df), &v, &dv),
        (work(&varied(f, df, 1.0), &varied(v, dv, 1.0))
            - work(&varied(f, df, -1.0), &varied(v, dv, -1.0)))
            / (2.0 * h),
    );
    close(work_derivative(&f, None, &v, &dv), work(&f, &dv));
    assert!(work_derivative(&f, Some(&df), &v, &dv) != work(&f, &dv));
}

#[test]
fn physical_recovery_refuses_invalid_data_and_cancellation_without_mutation() {
    let m = [2.0, 0.0, 0.0, 0.7];
    let k = [5.0, 0.0, 0.0, 3.0];
    let u = [C64::new(0.3, -0.1)];
    let v = [C64::new(0.7, 0.2)];
    let before = energy(&m, &k, &u, &v, &v, &u, &mut || false).unwrap();
    assert_eq!(
        interface(&[0.3], &v, &u, &mut || true),
        Err(ResponseError::Cancelled)
    );
    let mut polls = 0;
    let cancelled = energy(&m, &k, &u, &v, &v, &u, &mut || {
        polls += 1;
        polls == 3
    });
    assert_eq!(cancelled, Err(ResponseError::Cancelled));
    assert!(interface(&[], &v, &u, &mut || false).is_err());
    assert!(interface(&[f64::NAN], &v, &u, &mut || false).is_err());
    assert!(energy(&m, &[5.0, 0.0, 0.0, -3.0], &u, &v, &v, &u, &mut || false).is_err());
    assert!(energy(&m[..3], &k, &u, &v, &v, &u, &mut || false).is_err());
    assert_eq!(
        before,
        energy(&m, &k, &u, &v, &v, &u, &mut || false).unwrap()
    );
}
