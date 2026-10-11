use super::*;
use crate::kirchhoff_carrier_string;

fn parameters() -> KcStringParams {
    KcStringParams {
        length: 0.83,
        tension: 650.0,
        lin_density: 0.009,
        ea: 150_000.0,
    }
}

fn relative(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * actual.abs().max(expected.abs()).max(1e-18),
        "actual={actual:.17e}, expected={expected:.17e}"
    );
}

fn spatial_slope(plane: &[f64], fraction: f64) -> f64 {
    plane[0]
        + plane[1..]
            .iter()
            .enumerate()
            .map(|(n, a)| a * det::cos((n + 1) as f64 * core::f64::consts::PI * fraction))
            .sum::<f64>()
}

#[test]
fn energy_matches_independent_span_integral_with_explicit_static_condensation() {
    let p = parameters();
    let slopes = [
        0.003, -0.009, 0.014, -0.005, 0.007, -0.002, 0.005, 0.008, -0.011, 0.004,
    ];
    let n = 4;
    let samples = 2048;
    // Sample the continuum strain directly. The independent quadrature also
    // computes the omitted axial displacements that minimize that strain.
    let geometry: Vec<f64> = (0..samples)
        .map(|k| {
            let x = (k as f64 + 0.5) / samples as f64;
            slopes
                .chunks_exact(n + 1)
                .map(|plane| 0.5 * spatial_slope(plane, x).powi(2))
                .sum()
        })
        .collect();
    let measured_harmonics: Vec<f64> = (1..=2 * n)
        .map(|mode| {
            2.0 * geometry
                .iter()
                .enumerate()
                .map(|(k, s)| {
                    let x = (k as f64 + 0.5) / samples as f64;
                    s * det::cos(mode as f64 * core::f64::consts::PI * x)
                })
                .sum::<f64>()
                / samples as f64
        })
        .collect();
    for modes in [vec![], vec![2, 5], (1..=8).collect()] {
        let model = LongitudinalString::new(p, n, 2, &modes).unwrap();
        let axial: Vec<f64> = modes
            .iter()
            .map(|m| (1.0 - 0.19 * *m as f64) * 2e-7)
            .collect();
        let integral = geometry
            .iter()
            .enumerate()
            .map(|(k, geometric)| {
                let x = (k as f64 + 0.5) / samples as f64;
                let mut strain = *geometric;
                for mode in 1..=2 * n {
                    let axial_slope = if let Some(j) = modes.iter().position(|m| *m == mode) {
                        model.axial_slope_scales()[j] * axial[j]
                    } else {
                        -measured_harmonics[mode - 1]
                    };
                    strain += axial_slope * det::cos(mode as f64 * core::f64::consts::PI * x);
                }
                0.5 * p.ea * p.length * strain * strain
            })
            .sum::<f64>()
            / samples as f64;
        relative(model.energy(&slopes, &axial).unwrap(), integral, 3e-13);
    }
}

#[test]
fn force_matches_energy_derivatives_in_both_planes_and_every_axial_coordinate() {
    let model = LongitudinalString::new(parameters(), 4, 2, &[1, 2, 5, 8]).unwrap();
    let slopes = [
        0.003, -0.009, 0.014, -0.005, 0.007, -0.002, 0.005, 0.008, -0.011, 0.004,
    ];
    let axial = [2e-7, -3e-7, 1e-7, 4e-7];
    let mut gs = [0.0; 10];
    let mut gu = [0.0; 4];
    model
        .discrete_gradient(&slopes, &slopes, &axial, &axial, &mut gs, &mut gu)
        .unwrap();
    for i in 0..slopes.len() {
        let mut minus = slopes;
        let mut plus = slopes;
        let h = 1e-7;
        minus[i] -= h;
        plus[i] += h;
        let derivative = (model.energy(&plus, &axial).unwrap()
            - model.energy(&minus, &axial).unwrap())
            / (2.0 * h);
        relative(gs[i], derivative, 2e-8);
    }
    for i in 0..axial.len() {
        let mut minus = axial;
        let mut plus = axial;
        let h = 1e-11;
        minus[i] -= h;
        plus[i] += h;
        let derivative = (model.energy(&slopes, &plus).unwrap()
            - model.energy(&slopes, &minus).unwrap())
            / (2.0 * h);
        relative(gu[i], derivative, 2e-8);
    }
}

#[test]
fn finite_step_work_is_exact_and_linear_axial_storage_is_removed_only_from_residual() {
    let model = LongitudinalString::new(parameters(), 3, 2, &[1, 2, 4, 6]).unwrap();
    let a = [0.004, -0.02, 0.014, 0.007, -0.003, 0.006, -0.013, 0.009];
    let b = [-0.002, -0.013, 0.019, -0.004, 0.007, -0.008, -0.011, 0.003];
    let u = [2e-7, -4e-7, 1e-7, -2e-7];
    let v = [-5e-7, 2e-7, 3e-7, 1e-7];
    let mut gs = [0.0; 8];
    let mut gu = [0.0; 4];
    model
        .discrete_gradient(&a, &b, &u, &v, &mut gs, &mut gu)
        .unwrap();
    let work = gs
        .iter()
        .zip(a.iter().zip(b))
        .map(|(g, (a, b))| g * (b - a))
        .sum::<f64>()
        + gu.iter()
            .zip(u.iter().zip(v))
            .map(|(g, (a, b))| g * (b - a))
            .sum::<f64>();
    let delta = model.energy(&b, &v).unwrap() - model.energy(&a, &u).unwrap();
    relative(work, delta, 2e-13);
    let mut residual_s = [0.0; 8];
    let mut residual_u = [0.0; 4];
    model
        .coupling_discrete_gradient(&a, &b, &u, &v, &mut residual_s, &mut residual_u)
        .unwrap();
    assert_eq!(gs, residual_s);
    for j in 0..gu.len() {
        let linear = model.axial_omegas()[j].powi(2) * f64::midpoint(u[j], v[j]);
        relative(gu[j], residual_u[j] + linear, 3e-15);
    }
    let residual_work = residual_s
        .iter()
        .zip(a.iter().zip(b))
        .map(|(g, (a, b))| g * (b - a))
        .sum::<f64>()
        + residual_u
            .iter()
            .zip(u.iter().zip(v))
            .map(|(g, (a, b))| g * (b - a))
            .sum::<f64>();
    let linear_delta = model
        .axial_omegas()
        .iter()
        .zip(u.iter().zip(v))
        .map(|(omega, (a, b))| 0.5 * omega.powi(2) * (b * b - a * a))
        .sum::<f64>();
    relative(residual_work + linear_delta, delta, 3e-13);
}

#[test]
fn condensed_limit_is_existing_kirchhoff_carrier_and_single_sine_drives_twice_its_index() {
    let p = parameters();
    let kc = kirchhoff_carrier_string(&p, 3).unwrap();
    let coordinates = [0.00003, -0.00002, 0.00001];
    let normalization = det::sqrt(2.0 / (p.lin_density * p.length));
    let mut slopes = vec![0.0];
    slopes.extend(
        coordinates
            .iter()
            .enumerate()
            .map(|(j, q)| (j + 1) as f64 * core::f64::consts::PI / p.length * normalization * q),
    );
    let condensed = LongitudinalString::new(p, 3, 1, &[]).unwrap();
    // Read the established owner's stress channel directly instead of
    // subtracting a much larger linear energy to extract its quartic part.
    let channel = &kc.channels[0];
    let old_strain = (0..3)
        .map(|i| {
            (0..3)
                .map(|j| coordinates[i] * channel.coupling[3 * i + j] * coordinates[j])
                .sum::<f64>()
        })
        .sum::<f64>();
    relative(
        condensed.energy(&slopes, &[]).unwrap(),
        0.25 * channel.coefficient * old_strain * old_strain,
        1e-12,
    );

    let model = LongitudinalString::new(p, 3, 1, &[1, 2, 3, 4, 5, 6]).unwrap();
    let equilibrated: Vec<f64> = model
        .axial_modes()
        .iter()
        .enumerate()
        .map(|(j, mode)| -model.geometric_harmonic(&slopes, *mode) / model.axial_slope_scales()[j])
        .collect();
    relative(
        model.energy(&slopes, &equilibrated).unwrap(),
        condensed.energy(&slopes, &[]).unwrap(),
        2e-15,
    );
    let mut slope_gradient = [0.0; 4];
    let mut axial_gradient = [0.0; 6];
    model
        .discrete_gradient(
            &slopes,
            &slopes,
            &equilibrated,
            &equilibrated,
            &mut slope_gradient,
            &mut axial_gradient,
        )
        .unwrap();
    assert!(axial_gradient.iter().all(|f| f.abs() < 1e-10));

    let single = [0.0, 0.0, 0.02, 0.0]; // transverse partial 2
    let axial = [0.0; 6];
    model
        .coupling_discrete_gradient(
            &single,
            &single,
            &axial,
            &axial,
            &mut slope_gradient,
            &mut axial_gradient,
        )
        .unwrap();
    for (j, force) in axial_gradient.iter().enumerate() {
        if j == 3 {
            assert!(*force > 0.0);
        } else {
            assert_eq!(*force, 0.0);
        }
    }
    let d4 = 4.0 * core::f64::consts::PI / p.length * normalization;
    relative(
        axial_gradient[3],
        p.ea * p.length / 2.0 * d4 * 0.25 * 0.02_f64.powi(2),
        2e-15,
    );
}

#[test]
fn chord_pullback_is_reciprocal_and_strain_bounds_hold_on_the_whole_span() {
    let p = parameters();
    let model = LongitudinalString::new(p, 2, 1, &[1, 2, 3, 4]).unwrap();
    let scale = det::sqrt(2.0 / (p.lin_density * p.length)) * core::f64::consts::PI / p.length;
    let beta = [0.013, -0.008];
    let q = [0.00004, -0.00002];
    let bridge = 0.0007;
    let pull = |b: f64| {
        [
            b / p.length,
            scale * (q[0] - beta[0] * b),
            2.0 * scale * (q[1] - beta[1] * b),
        ]
    };
    let slopes = pull(bridge);
    let axial = [2e-7, -3e-7, 1e-7, 4e-8];
    let mut gs = [0.0; 3];
    let mut gu = [0.0; 4];
    model
        .discrete_gradient(&slopes, &slopes, &axial, &axial, &mut gs, &mut gu)
        .unwrap();
    let derivative = gs[0] / p.length - gs[1] * scale * beta[0] - gs[2] * 2.0 * scale * beta[1];
    let h = 1e-8;
    let finite_difference = (model.energy(&pull(bridge + h), &axial).unwrap()
        - model.energy(&pull(bridge - h), &axial).unwrap())
        / (2.0 * h);
    relative(derivative, finite_difference, 2e-8);
    let (mean, radius) = model.strain_bounds(&slopes, &axial).unwrap();
    assert!(mean > 0.0 && radius > 0.0);
    for j in 0..=4096 {
        let strain = model.strain_at(&slopes, &axial, j as f64 / 4096.0).unwrap();
        assert!(strain >= mean - radius - 1e-18 && strain <= mean + radius + 1e-18);
    }
}

#[test]
fn bounded_admission_keeps_requested_modes_and_refuses_invalid_states() {
    let p = parameters();
    let largest = LongitudinalString::new(
        p,
        MAX_TRANSVERSE_MODES,
        2,
        &(993..=MAX_LONGITUDINAL_INDEX).collect::<Vec<_>>(),
    )
    .unwrap();
    assert_eq!(largest.slope_count(), 1026);
    assert_eq!(largest.axial_modes().len(), MAX_LONGITUDINAL_MODES);
    let mut gs = vec![0.0; largest.slope_count()];
    let mut gu = vec![0.0; MAX_LONGITUDINAL_MODES];
    let slopes = vec![0.001; largest.slope_count()];
    let axial = vec![1e-12; MAX_LONGITUDINAL_MODES];
    largest
        .discrete_gradient(&slopes, &slopes, &axial, &axial, &mut gs, &mut gu)
        .unwrap();
    assert!(largest.energy(&slopes, &axial).unwrap() > 0.0);
    assert!(gs.iter().chain(&gu).all(|x| x.is_finite()));
    let chord = LongitudinalString::new(p, 0, 1, &[1]).unwrap();
    relative(
        chord.energy(&[0.1], &[0.0]).unwrap(),
        p.ea * p.length * 0.1_f64.powi(4) / 8.0,
        2e-15,
    );
    for (n, planes, modes) in [
        (513, 1, vec![1]),
        (1, 0, vec![1]),
        (1, 3, vec![1]),
        (1, 1, vec![0]),
        (1, 1, vec![2, 1]),
        (1, 1, vec![1, 1]),
        (1, 1, vec![1025]),
        (1, 1, (1..=33).collect()),
    ] {
        assert!(LongitudinalString::new(p, n, planes, &modes).is_err());
    }
    assert!(
        LongitudinalString::new(
            KcStringParams {
                ea: f64::INFINITY,
                ..p
            },
            1,
            1,
            &[1]
        )
        .is_err()
    );
    assert!(chord.energy(&[f64::NAN], &[0.0]).is_err());
    assert!(chord.energy(&[f64::MAX], &[0.0]).is_err());
    assert!(chord.energy(&[], &[0.0]).is_err());
    assert!(chord.strain_at(&[0.0], &[0.0], -0.01).is_err());
    assert!(chord.strain_bounds(&[0.0], &[f64::MAX]).is_err());
    assert!(
        chord
            .discrete_gradient(&[0.0], &[0.0], &[0.0], &[0.0], &mut [], &mut [0.0])
            .is_err()
    );
}
