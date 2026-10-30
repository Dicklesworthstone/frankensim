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
fn g1_material_strain_rate_matches_independent_continuum_projection() {
    let p = parameters();
    let n = 3;
    let modes: Vec<_> = (1..=2*n).collect();
    let law = LongitudinalString::new(p, n, 2, &modes).unwrap();
    let q = [0.004, -0.02, 0.014, 0.007, -0.003, 0.006, -0.013, 0.009];
    let v = [0.8, -1.2, 0.7, 0.4, -0.3, 0.1, 1.4, -0.6];
    let axial_v = [1e-4, -2e-4, 4e-4, 3e-4, -1e-4, 2e-4];
    let mut rates = [0.0; 6];
    let mean = law.strain_rate(&q, &v, &axial_v, &mut rates).unwrap();
    let samples = 1024;
    let normalization = (2.0 / (p.lin_density*p.length)).sqrt();
    let continuum: Vec<_> = (0..samples).map(|i| {
        let x = (i as f64 + 0.5) / samples as f64;
        let transverse: f64 = q.chunks_exact(n+1).zip(v.chunks_exact(n+1))
            .map(|(q,v)| spatial_slope(q,x)*spatial_slope(v,x)).sum();
        transverse + axial_v.iter().enumerate().map(|(j,v)| {
            let k = (j+1) as f64*core::f64::consts::PI/p.length;
            normalization*k*v*det::cos(k*p.length*x)
        }).sum::<f64>()
    }).collect();
    relative(mean, continuum.iter().sum::<f64>() / samples as f64, 5e-13);
    for (j, rate) in rates.iter().enumerate() {
        let projected = 2.0 * continuum.iter().enumerate().map(|(i,r)|
            r*det::cos((j+1) as f64*core::f64::consts::PI*(i as f64+0.5)/samples as f64))
            .sum::<f64>() / samples as f64;
        relative(*rate, projected, 5e-13);
    }
    let physical_loss = continuum.iter().map(|r|r*r).sum::<f64>()/samples as f64;
    relative(mean*mean+0.5*rates.iter().map(|r|r*r).sum::<f64>(),physical_loss,5e-13);
    assert!(law.strain_rate(&q, &v[..7], &axial_v, &mut rates).is_err());
    assert!(law.strain_rate(&q, &v, &[f64::NAN;6], &mut rates).is_err());
}

#[test]
fn g1_midpoint_material_rate_has_exact_finite_step_viscous_work() {
    let p = parameters();
    let law = LongitudinalString::new(p, 3, 2, &[1,2,4,6]).unwrap();
    let a = [0.004, -0.02, 0.014, 0.007, -0.003, 0.006, -0.013, 0.009];
    let b = [-0.002, -0.013, 0.019, -0.004, 0.007, -0.008, -0.011, 0.003];
    let u = [2e-7, -4e-7, 1e-7, -2e-7];
    let w = [-5e-7, 2e-7, 3e-7, 1e-7];
    let dt = 1e-4;
    let midpoint: Vec<_> = a.iter().zip(b).map(|(&a,b)|f64::midpoint(a,b)).collect();
    let rate: Vec<_> = a.iter().zip(b).map(|(&a,b)|(b-a)/dt).collect();
    let axial_rate: Vec<_> = u.iter().zip(w).map(|(&a,b)|(b-a)/dt).collect();
    let mut strain_rates = [0.0;4];
    let mean = law.strain_rate(&midpoint,&rate,&axial_rate,&mut strain_rates).unwrap();
    // Independently integrate actual endpoint strain increments. No force,
    // damping-ratio or per-mode loss formula supplies this continuum reference.
    let ca = 0.37; // Manufactured sum(eta_i A_i) [N s], not measured material data.
    let mut integral = 0.0;
    for i in 0..1024 {
        let x = (i as f64+0.5)/1024.0;
        let actual = (law.strain_at(&b,&w,x).unwrap()-law.strain_at(&a,&u,x).unwrap())/dt;
        let projected = mean+strain_rates.iter().zip(law.axial_modes()).map(|(r,n)|
            r*det::cos(*n as f64*core::f64::consts::PI*x)).sum::<f64>();
        assert!((actual-projected).abs()<5e-13);
        integral += ca*p.length*dt*actual*actual/1024.0;
    }
    let c = ca*p.length;
    let weights: Vec<_> = strain_rates.iter().map(|r|0.5*c*r).collect();
    let mut gradient = [0.0;8];
    law.geometric_strain_gradient(&midpoint,c*mean,&weights,&mut gradient).unwrap();
    let work = gradient.iter().zip(a.iter().zip(b)).map(|(g,(&a,b))|g*(b-a)).sum::<f64>()
        + weights.iter().zip(law.axial_slope_scales()).zip(u.iter().zip(w))
            .map(|((g,d),(&a,b))|g*d*(b-a)).sum::<f64>();
    assert!(integral>0.0);
    relative(work,integral,1e-12);
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
fn one_state_strain_storage_and_gradient_are_the_equal_endpoint_methods() {
    // Two planes; axial indices below, at and above the transverse count (both
    // product branches) on a count long enough for the four-way sums; and an
    // empty axial set (the Kirchhoff--Carrier image). Same quantities, summed
    // in another order: equal to rounding.
    let n = 11;
    let slopes: Vec<f64> = (0..2 * (n + 1)).map(|i| 0.01 * (0.7 * i as f64 + 0.3).sin()).collect();
    for axial_modes in [&[1usize, 2, 5, 11, 15][..], &[][..]] {
        let model = LongitudinalString::new(parameters(), n, 2, axial_modes).unwrap();
        let count = axial_modes.len();
        let axial: Vec<f64> = (0..count).map(|j| 1e-7 * (j as f64 - 1.5)).collect();
        let (mut harmonics, mut strain) = (vec![0.0; count], vec![0.0; count]);
        let mean = model.geometric_strain(&slopes, &mut harmonics).unwrap();
        let (energy, nonuniform) = model.storage_from_strain(mean, &harmonics, &axial, &mut strain).unwrap();
        relative(energy, model.energy(&slopes, &axial).unwrap(), 1e-14);
        let (bound_mean, bound_nonuniform) = model.strain_bounds(&slopes, &axial).unwrap();
        assert_eq!(mean.to_bits(), bound_mean.to_bits());
        relative(nonuniform, bound_nonuniform, 1e-14);
        let c = model.strain_energy_coefficient();
        let weights: Vec<f64> = strain.iter().map(|s| c * s).collect();
        let mut one = vec![0.0; slopes.len()];
        model.geometric_strain_gradient(&slopes, 2.0 * c * mean, &weights, &mut one).unwrap();
        let (mut two, mut axial_gradient) = (vec![0.0; slopes.len()], vec![0.0; count]);
        model.discrete_gradient(&slopes, &slopes, &axial, &axial, &mut two, &mut axial_gradient).unwrap();
        let scale = two.iter().fold(0.0_f64, |m, g| m.max(g.abs()));
        for (a, b) in one.iter().zip(&two) {
            assert!((a - b).abs() <= 1e-14 * scale, "{a:e} vs {b:e}");
        }
        for ((&g, &s), &d) in axial_gradient.iter().zip(&strain).zip(model.axial_slope_scales()) {
            relative(g, c * d * s, 1e-14);
        }
    }
    let model = LongitudinalString::new(parameters(), 4, 1, &[1, 3]).unwrap();
    assert!(model.geometric_strain(&[0.0; 5], &mut [0.0; 1]).is_err());
    assert!(model.geometric_strain(&[f64::NAN, 0.0, 0.0, 0.0, 0.0], &mut [0.0; 2]).is_err());
    assert!(model.storage_from_strain(0.0, &[0.0; 2], &[0.0; 1], &mut [0.0; 2]).is_err());
    assert!(model.geometric_strain_gradient(&[0.0; 5], 0.0, &[0.0; 2], &mut [0.0; 4]).is_err());
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
