use super::*;
use fs_math::{dd::Dd, special::bessel_j1};

fn air() -> Medium {
    Medium::air()
}

/// J0 for the oracle only: double-double series to 18, then the optimally
/// truncated Hankel expansion with angle addition, as fs-math's J1 owner.
fn bessel_j0(x: f64) -> f64 {
    let x = x.abs();
    if x <= 18.0 {
        let z = -(Dd::from_f64(x) * Dd::from_f64(x)) / Dd::from_f64(4.0);
        let (mut term, mut sum) = (Dd::ONE, Dd::ONE);
        for k in 1..=200u32 {
            let k = f64::from(k);
            term = term * z / Dd::from_f64(k * k);
            sum = sum + term;
            if term.abs().hi <= 1e-34 * sum.abs().hi.max(1.0) {
                break;
            }
        }
        return sum.to_f64();
    }
    let (mut p, mut q, mut term, mut last) = (0.0, 0.0, 1.0_f64, f64::INFINITY);
    for m in 0..128_u32 {
        if term.abs() > last {
            break;
        }
        last = term.abs();
        match m % 4 {
            0 => p += term,
            1 => q += term,
            2 => p -= term,
            _ => q -= term,
        }
        if term.abs() < 1e-18 {
            break;
        }
        let mf = f64::from(m);
        term *= -(2.0 * mf + 1.0) * (2.0 * mf + 1.0) / (8.0 * (mf + 1.0)) / x;
    }
    let (s, c) = det::sin_cos(x);
    let amp = det::sqrt(core::f64::consts::FRAC_2_PI / x) * core::f64::consts::FRAC_1_SQRT_2;
    amp * ((p + q) * c + (p - q) * s)
}

/// Exact coplanar baffled-piston mutual impedance in its Hankel-transform
/// (King) representation, independently of the product's disc quadrature:
/// `Z_ab / (rho c 2 pi a b) = INT_0^inf J1(as) J1(bs) J0(ds) / (s sqrt(1-s^2)) ds`
/// with wavenumber-scaled radii and separation, and the outgoing branch
/// `sqrt(1-s^2) = +i sqrt(s^2-1)` past s = 1 (pinned below by the d = 0 limit).
/// The tail past `s_max` is bounded by `INT |f| <= c s^-3.5`; this returns
/// the integral and that bound (infinite at delta = 0, where it does not apply).
fn king_mutual(alpha: f64, beta: f64, delta: f64, s_max: f64) -> (C64, f64) {
    let f = |s: f64| bessel_j1(alpha * s) * bessel_j1(beta * s) * bessel_j0(delta * s) / s;
    let (x, w) = gauss_legendre(24);
    let panel = |lo: f64, hi: f64, g: &dyn Fn(f64) -> f64| -> f64 {
        x.iter()
            .zip(&w)
            .map(|(&x, &w)| 0.5 * (hi - lo) * w * g(0.5 * (hi - lo) * x + 0.5 * (hi + lo)))
            .sum::<f64>()
    };
    let width = (core::f64::consts::PI / (alpha + beta + delta)).min(0.25);
    // s = sin(theta) on [0, 1] and s = cosh(t) on [1, 2] remove the endpoint
    // singularities; [2, s_max] is smooth and oscillatory.
    let mut real = 0.0;
    let mut lo = 0.0;
    while lo < core::f64::consts::FRAC_PI_2 {
        let hi = (lo + width).min(core::f64::consts::FRAC_PI_2);
        real += panel(lo, hi, &|t| f(det::sin(t)));
        lo = hi;
    }
    let mut imag = 0.0;
    let t_end = det::ln(2.0 + det::sqrt(3.0));
    lo = 0.0;
    while lo < t_end {
        let hi = (lo + width).min(t_end);
        imag += panel(lo, hi, &|t| f(0.5 * (det::exp(t) + det::exp(-t))));
        lo = hi;
    }
    lo = 2.0;
    while lo < s_max {
        let hi = (lo + width).min(s_max);
        imag += panel(lo, hi, &|s| f(s) / det::sqrt(s * s - 1.0));
        lo = hi;
    }
    // Bessel moduli: |J0(y)| <= sqrt(2/(pi y)) and |J1(y)| <= sqrt(2/(pi y))
    // (1 + 3/(16 y^2) + ...), within 1.001 for y >= 25 (alpha, beta >= 0.1 and
    // s_max >= 256 below). 1/sqrt(s^2-1) <= 1.16/s for s >= 2.
    let two_over_pi = 2.0 / core::f64::consts::PI;
    let bound = two_over_pi * det::sqrt(two_over_pi) * 1.001 * 1.16
        / det::sqrt(alpha * beta * delta)
        / (2.5 * s_max * s_max * det::sqrt(s_max));
    (C64::new(real, -imag), bound)
}

/// Double `s_max` from 256 until the tail bound is below 1e-7 of the value.
fn king_mutual_resolved(alpha: f64, beta: f64, delta: f64) -> (C64, f64) {
    let mut s_max = 256.0;
    loop {
        let (value, bound) = king_mutual(alpha, beta, delta, s_max);
        if bound <= 1e-7 * value.abs() || s_max > 1.0e5 {
            return (value, bound);
        }
        s_max *= 2.0;
    }
}

#[test]
fn two_piston_mutual_vs_pritchard() {
    // The oracle is the exact King (Hankel-transform) representation that
    // Pritchard (1960) expands into his series. The series coefficients were
    // not available here, so the same exact integral is evaluated directly,
    // with an explicit tail bound, independently of the product's quadrature.
    let medium = air();
    let rho_c = medium.density * medium.sound_speed;
    // Branch and normalization: at zero separation the representation is the
    // self impedance 1 - J1(2ka)/ka - i H1(2ka)/ka of the closed-form owner.
    for ka in [0.1, 0.7, 2.0] {
        // The s^-3 tail here is not covered by the bound; 1e4 keeps it below
        // 1e-7 of the reactance for these ka.
        let (value, _) = king_mutual(ka, ka, 0.0, 1.0e4);
        let k = ka / 0.01;
        let closed = baffled_piston_face_impedance(0.01, k, medium)
            .unwrap()
            .scale(1.0 / rho_c);
        let rel = (value.scale(2.0) - closed).abs() / closed.abs();
        assert!(rel < 1e-6, "King d=0 limit at ka={ka}: {rel:.3e}");
    }
    let a = 0.01;
    let mut rows = String::new();
    for &(ka, kb, kd) in &[
        (0.1, 0.1, 0.5),
        (0.2, 0.2, 0.5),
        (0.5, 0.5, 1.5),
        (0.5, 0.5, 5.0),
        (1.0, 1.0, 2.5),
        (1.0, 0.5, 2.0),
        (1.0, 1.0, 10.0),
        (1.5, 1.5, 4.0),
        (2.0, 2.0, 5.0),
        (2.0, 2.0, 20.0),
    ] {
        let k = ka / a;
        let b = kb / k;
        let d = kd / k;
        let pistons = [
            BaffledPiston {
                center_m: [0.0, 0.0],
                radius_m: a,
            },
            BaffledPiston {
                center_m: [d * 0.6, d * 0.8],
                radius_m: b,
            },
        ];
        let z = baffled_piston_array_impedance(&pistons, k, medium).unwrap();
        let (king, bound) = king_mutual_resolved(ka, kb, kd);
        let oracle = king.scale(rho_c * core::f64::consts::TAU * a * b);
        assert!(
            bound <= 1e-6 * king.abs(),
            "oracle tail bound {bound:.2e} vs {:.2e}",
            king.abs()
        );
        for (index, value) in [(1, z[1]), (2, z[2])] {
            let rel = (value - oracle).abs() / oracle.abs();
            assert!(
                rel <= 1e-4,
                "ka={ka} kb={kb} kd={kd} entry {index}: {rel:.3e}"
            );
            if index == 1 {
                rows.push_str(&format!(
                    "{{\"ka\":{ka},\"kb\":{kb},\"kd\":{kd},\"rel\":{rel:.3e}}},"
                ));
            }
        }
        for (index, radius) in [(0, a), (3, b)] {
            let closed = baffled_piston_face_impedance(radius, k, medium)
                .unwrap()
                .scale(core::f64::consts::PI * radius * radius);
            assert!((z[index] - closed).abs() <= 1e-10 * closed.abs());
        }
    }
    rows.pop();
    println!("{{\"suite\":\"fs-bem-aperture\",\"case\":\"two-piston-mutual\",\"rows\":[{rows}]}}");
}

#[test]
fn baffled_piston_array_refuses_overlap_and_nearly_touching_discs() {
    let k = 100.0;
    let touch = [
        BaffledPiston {
            center_m: [0.0, 0.0],
            radius_m: 0.01,
        },
        BaffledPiston {
            center_m: [0.02, 0.0],
            radius_m: 0.01,
        },
    ];
    assert!(baffled_piston_array_impedance(&touch, k, air()).is_err());
    // A 1e-6 gap is admitted geometrically but the ladder cannot resolve it.
    let close = [
        BaffledPiston {
            center_m: [0.0, 0.0],
            radius_m: 0.01,
        },
        BaffledPiston {
            center_m: [0.020_000_01, 0.0],
            radius_m: 0.01,
        },
    ];
    assert!(baffled_piston_array_impedance(&close, k, air()).is_err());
    assert!(baffled_piston_array_impedance(&[], k, air()).is_err());
}

fn cap(surface: &SpherePanels, axis: [f64; 3], cos_half_angle: f64) -> Aperture {
    let shape = surface
        .centroids()
        .iter()
        .map(|c| {
            let r = det::sqrt(c[0] * c[0] + c[1] * c[1] + c[2] * c[2]);
            let along = (c[0] * axis[0] + c[1] * axis[1] + c[2] * axis[2]) / r;
            if along > cos_half_angle { 1.0 } else { 0.0 }
        })
        .collect();
    Aperture {
        shape,
        reference_point: axis,
        normal: axis,
    }
}

#[test]
fn multi_aperture_reciprocity_passivity() {
    // Pistons and Galerkin ports are reciprocal and passive by construction:
    // the gate is roundoff. Collocation is measured and reported.
    let pistons = [
        BaffledPiston {
            center_m: [0.0, 0.0],
            radius_m: 0.01,
        },
        BaffledPiston {
            center_m: [0.035, 0.0],
            radius_m: 0.008,
        },
        BaffledPiston {
            center_m: [0.0, 0.05],
            radius_m: 0.012,
        },
    ];
    for k in [10.0, 60.0, 150.0] {
        let z = baffled_piston_array_impedance(&pistons, k, air()).unwrap();
        let (defect, margin) = reciprocity_and_passivity(&z, 3);
        assert!(defect <= 1e-12, "piston reciprocity {defect:.3e} at k={k}");
        assert!(margin >= -1e-12, "piston passivity {margin:.3e} at k={k}");
    }
    let surface = SpherePanels::icosphere(1.0, 2).unwrap();
    let apertures = partition(&surface);
    let medium = Medium {
        density: 1.0,
        sound_speed: 1.0,
    };
    let gate = CancelGate::new_clock_free();
    let omegas = [0.3, 0.6, 0.9];
    let galerkin = aperture_impedance(
        &surface,
        &apertures,
        &omegas,
        medium,
        &ApertureRoute::Galerkin { budget: budget() },
        &gate,
    )
    .unwrap();
    let collocation = aperture_impedance(
        &surface,
        &apertures,
        &omegas,
        medium,
        &ApertureRoute::Collocation { iterative: None },
        &gate,
    )
    .unwrap();
    assert_ne!(galerkin.declaration, collocation.declaration);
    let area: f64 = surface.areas().iter().sum();
    for (g, c) in galerkin.rows.iter().zip(&collocation.rows) {
        assert_eq!((g.route, c.route), ("galerkin", "dense"));
        assert!(
            g.reciprocity_defect <= 1e-12,
            "galerkin reciprocity {:.3e}",
            g.reciprocity_defect
        );
        assert!(
            g.passivity_margin >= -1e-12,
            "galerkin passivity {:.3e}",
            g.passivity_margin
        );
        // The apertures partition the sphere, so unit velocity on all of them
        // is the pulsating sphere: the sum of every entry is exact for the
        // discrete shapes, whatever their staircase boundaries.
        let oracle = pulsating(g.omega).scale(area);
        let error = |row: &ApertureRow| {
            (row.impedance.iter().fold(C64::ZERO, |s, &z| s + z) - oracle).abs() / oracle.abs()
        };
        let (eg, ec) = (error(g), error(c));
        println!(
            "{{\"suite\":\"fs-bem-aperture\",\"case\":\"reciprocity-passivity\",\"omega\":{},\"galerkin_defect\":{:.3e},\"galerkin_margin\":{:.3e},\"collocation_defect\":{:.3e},\"collocation_margin\":{:.3e},\"galerkin_pulsating_rel\":{eg:.3e},\"collocation_pulsating_rel\":{ec:.3e}}}",
            g.omega,
            g.reciprocity_defect,
            g.passivity_margin,
            c.reciprocity_defect,
            c.passivity_margin
        );
        // Measured on trj (b179): 7.4e-3 / 6.8e-3 / 6.1e-3 Galerkin and
        // 7.4e-3 / 6.7e-3 / 5.8e-3 collocation at omega 0.3 / 0.6 / 0.9.
        assert!(
            eg < 0.03 && ec < 0.03,
            "pulsating partition: {eg:.3e} {ec:.3e}"
        );
    }
}

/// Mean pressure per unit velocity of a pulsating unit sphere (rho = c = 1).
fn pulsating(ka: f64) -> C64 {
    C64::new(ka * ka, -ka).scale(1.0 / (1.0 + ka * ka))
}

fn budget() -> Budget {
    Budget {
        maximum_unknowns: 4096,
        maximum_cache_bytes: 1 << 30,
        maximum_work: 200_000_000_000,
        angular_refinement: 0,
    }
}

/// Two disjoint caps and the rest of the sphere.
fn partition(surface: &SpherePanels) -> Vec<Aperture> {
    let mut apertures = vec![
        cap(surface, [0.0, 0.0, 1.0], 0.85),
        cap(surface, [-1.0, 0.0, 0.0], 0.9),
    ];
    let rest = (0..surface.areas().len())
        .map(|i| 1.0 - apertures[0].shape[i] - apertures[1].shape[i])
        .collect();
    apertures.push(Aperture {
        shape: rest,
        reference_point: [0.0, 1.0, 0.0],
        normal: [0.0, 1.0, 0.0],
    });
    apertures
}

#[test]
fn aperture_inputs_refuse_before_solving() {
    let surface = SpherePanels::icosphere(1.0, 1).unwrap();
    let n = surface.areas().len();
    let good = Aperture {
        shape: vec![1.0; n],
        reference_point: [0.0, 0.0, 1.0],
        normal: [0.0, 0.0, 1.0],
    };
    let medium = Medium {
        density: 1.0,
        sound_speed: 1.0,
    };
    let gate = CancelGate::new_clock_free();
    let route = ApertureRoute::Collocation { iterative: None };
    let run = |apertures: &[Aperture], omegas: &[f64]| {
        aperture_impedance(&surface, apertures, omegas, medium, &route, &gate)
    };
    assert!(run(&[], &[0.5]).is_err());
    assert!(
        run(
            &[Aperture {
                shape: vec![1.0; n - 1],
                ..good.clone()
            }],
            &[0.5]
        )
        .is_err()
    );
    assert!(
        run(
            &[Aperture {
                shape: vec![0.0; n],
                ..good.clone()
            }],
            &[0.5]
        )
        .is_err()
    );
    let mut nan = good.clone();
    nan.shape[2] = f64::NAN;
    assert!(run(&[nan], &[0.5]).is_err());
    assert!(
        run(
            &[Aperture {
                normal: [0.0, 0.0, 2.0],
                ..good.clone()
            }],
            &[0.5]
        )
        .is_err()
    );
    assert!(run(&[good.clone()], &[]).is_err());
    assert!(run(&[good.clone()], &[0.5, 0.5]).is_err());
    assert!(run(&[good.clone()], &[-0.5]).is_err());
    assert!(run(&[good.clone()], &[0.5]).is_ok());
    let stopped = CancelGate::new_clock_free();
    stopped.request();
    assert!(matches!(
        aperture_impedance(&surface, &[good], &[0.5], medium, &route, &stopped),
        Err(ApertureError::Helmholtz(HelmholtzError::Cancelled))
    ));
}
