use super::*;

fn close(a: C64, b: C64, tolerance: f64) {
    assert!(
        (a - b).abs() < tolerance * (1.0 + a.abs() + b.abs()),
        "{a:?} != {b:?}"
    );
}

// Independent exponential pair, including a wave launched at the far face.
fn reference(kz: C64, y: C64, f: C64, b: C64, d: f64, z: f64) -> [C64; 2] {
    let exponential = |x: f64| {
        let (s, c) = det::sin_cos(kz.re * x);
        C64::new(c, -s).scale(det::exp(kz.im * x))
    };
    let forward = f * exponential(z);
    let backward = b * exponential(d - z);
    [forward + backward, y * (forward - backward)]
}

#[test]
fn internal_fields_match_independent_bidirectional_waves_through_critical_neighborhood() {
    for kz in [
        C64::from_re(2.0),
        C64::new(1.1, -0.4),
        C64::new(0.0, -1.2),
        C64::new(1e-4, -1e-5),
        C64::from_re(0.499),
        C64::from_re(0.501),
    ] {
        for tm in [false, true] {
            let (a, b, y) = if tm {
                let b = C64::new(2.2, -0.5);
                ((kz * kz) / b, b, b / kz)
            } else {
                let a = C64::from_re(1.7);
                (a, (kz * kz) / a, kz / a)
            };
            let f = C64::new(0.7, 0.2);
            let back = C64::new(0.3, -0.1);
            let d = 1.0;
            let left = reference(kz, y, f, back, d, 0.0);
            let right = reference(kz, y, f, back, d, d);
            let kernel = LayerFieldKernel::new(d, kz, a, b, left, right).unwrap();
            for z in [0.0, 1e-8, 0.1, 0.5, 0.9, 1.0] {
                let expected = reference(kz, y, f, back, d, z);
                let got = kernel.field(z).unwrap();
                for c in 0..2 {
                    close(got[c], expected[c], 3e-13);
                }
            }
        }
    }
}

#[test]
fn exactly_critical_te_and_tm_keep_the_linear_maxwell_limit() {
    let left = [C64::new(0.8, -0.3), C64::new(0.4, 0.7)];
    let d = 3.2;
    for (a, b) in [
        (C64::from_re(1.3), C64::ZERO),
        (C64::ZERO, C64::from_re(2.0)),
    ] {
        let expected = |z: f64| {
            [
                left[0] - C64::new(0.0, z) * a * left[1],
                left[1] - C64::new(0.0, z) * b * left[0],
            ]
        };
        let kernel = LayerFieldKernel::new(d, C64::ZERO, a, b, left, expected(d)).unwrap();
        for z in [0.0, 0.01, 0.7, d] {
            let got = kernel.field(z).unwrap();
            let reference = expected(z);
            close(got[0], reference[0], 2e-15);
            close(got[1], reference[1], 2e-15);
        }
        for z in [-0.1, d + 0.1, f64::NAN, f64::INFINITY] {
            assert!(kernel.field(z).is_none());
        }
    }
}

#[test]
fn integrated_vector_dielectric_heat_matches_signed_poynting_loss_for_both_polarizations() {
    let kz = C64::new(0.9, -0.2);
    let omega_mu = 1.3;
    let kt = 0.7;
    let omega_epsilon = (kz * kz + C64::from_re(kt * kt)).scale(1.0 / omega_mu);
    let thickness = 1.6;
    for tm in [false, true] {
        let (a, b, y) = if tm {
            ((kz * kz) / omega_epsilon, omega_epsilon, omega_epsilon / kz)
        } else {
            (
                C64::from_re(omega_mu),
                (kz * kz).scale(1.0 / omega_mu),
                kz.scale(1.0 / omega_mu),
            )
        };
        let left = reference(
            kz,
            y,
            C64::new(0.9, 0.2),
            C64::new(-0.3, 0.4),
            thickness,
            0.0,
        );
        let right = reference(
            kz,
            y,
            C64::new(0.9, 0.2),
            C64::new(-0.3, 0.4),
            thickness,
            thickness,
        );
        let kernel = LayerFieldKernel::new(thickness, kz, a, b, left, right).unwrap();
        let flux = |field: [C64; 2]| 0.5 * (field[0] * field[1].conj()).re;
        let heat = |z: f64| {
            let field = kernel.field(z).unwrap();
            let ez = if tm {
                -field[1].scale(kt) / omega_epsilon
            } else {
                C64::ZERO
            };
            -0.5 * omega_epsilon.im * (field[0].norm_sq() + ez.norm_sq())
        };
        // Composite Simpson is an independent volume integration, not a
        // rearrangement of the boundary flux or positive slab loss factors.
        let count = 512;
        let dz = thickness / count as f64;
        let mut integral = heat(0.0) + heat(thickness);
        for i in 1..count {
            integral += if i % 2 == 0 { 2.0 } else { 4.0 } * heat(i as f64 * dz);
        }
        integral *= dz / 3.0;
        let expected = flux(left) - flux(right);
        assert!(expected > 0.0);
        assert!(
            (integral - expected).abs() < 2e-11 * expected,
            "{integral} != {expected}"
        );
    }
}

#[test]
fn negative_real_permittivity_retains_metal_fields_and_positive_volume_heat() {
    let omega_epsilon = C64::new(-2.3, -0.45);
    let omega_mu = 1.3;
    let kt = 0.7;
    let kz = (omega_epsilon.scale(omega_mu) - C64::from_re(kt * kt)).sqrt();
    assert!(kz.re > 0.0 && kz.im < 0.0);
    // Exercise both the entire small-depth propagator and face-decaying waves.
    for thickness in [0.05, 1.6] {
        for tm in [false, true] {
            let (a, b, y) = if tm {
                ((kz * kz) / omega_epsilon, omega_epsilon, omega_epsilon / kz)
            } else {
                (
                    C64::from_re(omega_mu),
                    (kz * kz).scale(1.0 / omega_mu),
                    kz.scale(1.0 / omega_mu),
                )
            };
            let independent =
                |z| reference(kz, y, C64::new(0.9, 0.2), C64::new(-0.3, 0.4), thickness, z);
            let left = independent(0.0);
            let right = independent(thickness);
            let kernel = LayerFieldKernel::new(thickness, kz, a, b, left, right).unwrap();
            let heat = |z| {
                let field = kernel.field(z).unwrap();
                let expected = independent(z);
                close(field[0], expected[0], 3e-13);
                close(field[1], expected[1], 3e-13);
                let ez = if tm {
                    -field[1].scale(kt) / omega_epsilon
                } else {
                    C64::ZERO
                };
                -0.5 * omega_epsilon.im * (field[0].norm_sq() + ez.norm_sq())
            };
            let intervals = 512;
            let dz = thickness / intervals as f64;
            let mut integrated = heat(0.0) + heat(thickness);
            for i in 1..intervals {
                integrated += if i % 2 == 0 { 2.0 } else { 4.0 } * heat(i as f64 * dz);
            }
            integrated *= dz / 3.0;
            let expected = 0.5 * ((left[0] * left[1].conj()).re - (right[0] * right[1].conj()).re);
            assert!(expected > 0.0 && integrated > 0.0);
            assert!(
                (integrated - expected).abs() < 2e-10 * expected,
                "metal volume heat {integrated} != boundary work {expected}"
            );
        }
    }
}

#[test]
fn opaque_film_never_needs_a_growing_wave_or_an_unobservable_underflowed_phase() {
    let kz = C64::new(2.0, -800.0);
    let a = C64::ONE;
    let b = kz * kz;
    let kernel = LayerFieldKernel::new(1e308, kz, a, b, [C64::ONE, kz], [C64::ZERO; 2]).unwrap();
    let near = kernel.field(0.001).unwrap();
    assert!(near[0].abs() > 0.4 && near[0].abs() < 0.5);
    assert_eq!(kernel.field(0.999e308), Some([C64::ZERO; 2]));
    assert_eq!(kernel.field(1e308), Some([C64::ZERO; 2]));
    assert!(
        LayerFieldKernel::new(1.0, C64::new(1.0, 0.1), a, b, [C64::ONE; 2], [C64::ONE; 2])
            .is_none()
    );
}
