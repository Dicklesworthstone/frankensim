use super::*;

const MU: f64 = 1.256_637_062_12e-6;
fn input(ni: f64, no: f64, ne: f64, nt: f64, depth: f64, axis: [f64; 3]) -> Input {
    Input {
        crystal: ResolvedFrequency {
            omega: C0,
            ni,
            no,
            ne,
            mui: MU,
            mut_: MU,
        },
        exit_index: nt,
        exit_mu_h_m: MU,
        thickness_m: depth,
        axis,
    }
}
fn near(a: f64, b: f64, t: f64) {
    assert!((a - b).abs() <= t, "{a:.17e} != {b:.17e}; tolerance={t}");
}
fn complex_near(a: C64, b: C64, t: f64) {
    assert!((a - b).abs() <= t, "{a:?} != {b:?}; tolerance={t}");
}
fn phase(v: f64) -> C64 {
    C64::new(v.cos(), -v.sin())
}
fn cross(a: Vector, b: Vector) -> Vector {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

// Independent scalar Airy sum, using tangential-E admittances. The complex
// propagation factor is supplied from its analytic homogeneous dispersion.
fn airy(yi: C64, yc: C64, yt: C64, p: C64) -> (C64, C64) {
    let r0 = (yi - yc) / (yi + yc);
    let r1 = (yc - yt) / (yc + yt);
    let d = C64::ONE + r0 * r1 * p * p;
    (
        (r0 + r1 * p * p) / d,
        yi.scale(2.0) / (yi + yc) * yc.scale(2.0) / (yc + yt) * p / d,
    )
}

#[test]
fn isotropic_oblique_slab_matches_independent_te_tm_airy_and_internal_flux() {
    for mu_ratio in [0.7, 1.0, 1.3] {
        let mut d = input(1.2, 1.7, 1.7, 1.4, 0.83, [1.0, 2.0, 3.0]);
        d.crystal.mut_ = MU * mu_ratio;
        let p = 0.65;
        let qi = (d.crystal.ni * d.crystal.ni - p * p).sqrt();
        let qc = (d.crystal.no * d.crystal.no - p * p).sqrt();
        let qt = (d.exit_index * d.exit_index - p * p).sqrt();
        for tm in [false, true] {
            let e = if tm {
                [C64::new(0.7, -0.4), C64::ZERO]
            } else {
                [C64::ZERO, C64::new(0.7, -0.4)]
            };
            let index = if tm { 0 } else { 1 };
            let (yi, yc, yt) = if tm {
                (
                    d.crystal.ni * d.crystal.ni / qi,
                    d.crystal.no * d.crystal.no / (qc * mu_ratio),
                    d.exit_index * d.exit_index / qt,
                )
            } else {
                (qi, qc / mu_ratio, qt)
            };
            let (r, t) = airy(
                C64::from_re(yi),
                C64::from_re(yc),
                C64::from_re(yt),
                phase(qc * d.thickness_m),
            );
            let actual = solve(d, Incidence::TangentialIndex(p), e).unwrap();
            complex_near(actual.reflected_electric_v_m[index], e[index] * r, 3e-15);
            complex_near(actual.transmitted_electric_v_m[index], e[index] * t, 3e-15);
            near(actual.reflected_power_fraction, r.norm_sq(), 3e-15);
            near(
                actual.transmitted_power_fraction,
                t.norm_sq() * yt / yi,
                4e-15,
            );
            for z in [0.0, 0.19, 0.57, d.thickness_m] {
                let field = actual
                    .field_at(SlabRegion::Crystal, [0.4, -0.2, z])
                    .unwrap();
                near(
                    field.poynting_w_m2[2],
                    actual.incident_normal_power_w_m2 * actual.transmitted_power_fraction,
                    3e-17,
                );
            }
        }
    }
}

#[test]
fn birefringent_etalon_retains_both_complex_jones_channels_and_resonances() {
    let d = input(1.0, 1.5, 2.0, 1.0, core::f64::consts::PI, [1.0, 0.0, 0.0]);
    let e = [C64::new(0.7, 0.2), C64::new(-0.3, 0.5)];
    let result = solve(d, Incidence::Cosine(1.0), e).unwrap();
    for (j, n) in [(0, 2.0), (1, 1.5)] {
        let (r, t) = airy(
            C64::ONE,
            C64::from_re(n),
            C64::ONE,
            phase(n * d.thickness_m),
        );
        complex_near(result.reflected_electric_v_m[j], e[j] * r, 4e-15);
        complex_near(result.transmitted_electric_v_m[j], e[j] * t, 4e-15);
    }
    // Extraordinary phase 2*pi is an exact lossless etalon resonance; the
    // ordinary channel is antiresonant and retains a finite reflected field.
    assert!(result.reflected_electric_v_m[0].abs() < 2e-15);
    assert!(result.reflected_electric_v_m[1].abs() > 0.1);
    near(
        result.reflected_power_fraction + result.transmitted_power_fraction,
        1.0,
        3e-15,
    );
}

#[test]
fn azimuthal_crystal_slab_has_analytic_cross_polarized_reflection() {
    let d = input(1.0, 1.4, 1.9, 1.2, 0.7, [1.0, 1.0, 0.0]);
    let result = solve(d, Incidence::Cosine(1.0), [C64::ONE, C64::ZERO]).unwrap();
    let (ro, to) = airy(
        C64::ONE,
        C64::from_re(1.4),
        C64::from_re(1.2),
        phase(1.4 * 0.7),
    );
    let (re, te) = airy(
        C64::ONE,
        C64::from_re(1.9),
        C64::from_re(1.2),
        phase(1.9 * 0.7),
    );
    complex_near(
        result.reflected_electric_v_m[0],
        (re + ro).scale(0.5),
        4e-15,
    );
    complex_near(
        result.reflected_electric_v_m[1],
        (re - ro).scale(0.5),
        4e-15,
    );
    complex_near(
        result.transmitted_electric_v_m[0],
        (te + to).scale(0.5),
        4e-15,
    );
    complex_near(
        result.transmitted_electric_v_m[1],
        (te - to).scale(0.5),
        4e-15,
    );
    assert!(result.reflected_electric_v_m[1].abs() > 0.05);
}

#[test]
fn backward_tilted_modes_obey_axial_h_reversal_and_both_maxwell_curls() {
    // Extraordinary positive energy has a negative phase normal here. The
    // backward normal is the other quadratic root, not minus the forward one.
    for axis in [[1.0, 0.0, 1.0], [0.6, 0.4, 0.8]] {
        let d = input(2.0, 1.0, 2.0, 1.8, 0.53, axis);
        let p = 1.5;
        let result = solve(
            d,
            Incidence::TangentialIndex(p),
            [C64::new(0.7, 0.2), C64::new(-0.3, 0.5)],
        )
        .unwrap();
        let length = axis[0].hypot(axis[1]).hypot(axis[2]);
        let a = axis.map(|v| v / length);
        for (forward, waves) in [(true, result.forward), (false, result.backward)] {
            for wave in waves {
                let k = [C64::from_re(p), C64::ZERO, wave.normal_wavenumber_rad_m];
                let projection =
                    (0..3).fold(C64::ZERO, |v, j| v + wave.electric_v_m[j].scale(a[j]));
                let displacement: Vector = core::array::from_fn(|j| {
                    (wave.electric_v_m[j].scale(d.crystal.no * d.crystal.no)
                        + projection.scale(
                            (d.crystal.ne * d.crystal.ne - d.crystal.no * d.crystal.no) * a[j],
                        ))
                    .scale(1.0 / (MU * C0 * C0))
                });
                let faraday = cross(k, wave.electric_v_m);
                let ampere = cross(k, wave.magnetic_a_m);
                for j in 0..3 {
                    complex_near(faraday[j], wave.magnetic_a_m[j].scale(C0 * MU), 8e-15);
                    complex_near(ampere[j], displacement[j].scale(-C0), 8e-17);
                }
                if wave.wave_kind == UniaxialWaveKind::Propagating {
                    assert!(if forward {
                        flux(wave.electric_v_m, wave.magnetic_a_m) >= 0.0
                    } else {
                        flux(wave.electric_v_m, wave.magnetic_a_m) <= 0.0
                    });
                }
            }
        }
        if axis[1] == 0.0 {
            assert!(result.forward[1].normal_wavenumber_rad_m.re < 0.0);
            assert!(
                (result.backward[1].normal_wavenumber_rad_m
                    + result.forward[1].normal_wavenumber_rad_m)
                    .abs()
                    > 0.1
            );
        }
        for z in [0.0, 0.11, 0.37, d.thickness_m] {
            let field = result.field_at(SlabRegion::Crystal, [0.3, 0.0, z]).unwrap();
            near(
                field.poynting_w_m2[2],
                result.incident_normal_power_w_m2 * result.transmitted_power_fraction,
                5e-17,
            );
        }
    }
}

#[test]
fn finite_evanescent_crystal_tunnels_without_growing_exponentials() {
    let p = 1.2;
    let qi = (1.5_f64 * 1.5 - p * p).sqrt();
    let gamma = (p * p - 1.0).sqrt();
    for depth in [0.05, 0.7, 8.0, 200.0] {
        let d = input(1.5, 1.0, 1.0, 1.5, depth, [1.0, 2.0, 3.0]);
        let (r, t) = airy(
            C64::from_re(qi),
            C64::new(0.0, -gamma),
            C64::from_re(qi),
            C64::from_re((-gamma * depth).exp()),
        );
        let actual = solve(d, Incidence::TangentialIndex(p), [C64::ZERO, C64::ONE]).unwrap();
        complex_near(actual.reflected_electric_v_m[1], r, 4e-15);
        assert!((actual.transmitted_electric_v_m[1] / t - C64::ONE).abs() < 2e-13);
        assert!((actual.transmitted_power_fraction / t.norm_sq() - 1.0).abs() < 4e-13);
        assert!(actual.transmitted_power_fraction > 0.0);
        assert!(
            actual
                .forward
                .iter()
                .all(|w| w.normal_wavenumber_rad_m.im < 0.0)
        );
        assert!(
            actual
                .backward
                .iter()
                .all(|w| w.normal_wavenumber_rad_m.im > 0.0)
        );
        if depth < 10.0 {
            let center = actual
                .field_at(SlabRegion::Crystal, [0.0, 0.0, 0.5 * depth])
                .unwrap();
            near(
                center.poynting_w_m2[2],
                actual.incident_normal_power_w_m2 * actual.transmitted_power_fraction,
                4e-17,
            );
        }
    }
    // A finite barrier has positive tunneling power. Once that positive flux
    // underflows, refuse instead of reporting a fabricated exact extinction.
    assert!(
        solve(
            input(1.5, 1.0, 1.0, 1.5, 1000.0, [1.0, 2.0, 3.0]),
            Incidence::TangentialIndex(p),
            [C64::ZERO, C64::ONE]
        )
        .is_err()
    );
}

#[test]
fn critical_or_evanescent_exit_retains_fields_with_zero_normal_power() {
    for nt in [0.9, 1.2] {
        let d = input(1.5, 1.7, 2.0, nt, 0.6, [1.0, 2.0, 3.0]);
        let result = solve(
            d,
            Incidence::TangentialIndex(1.2),
            [C64::ONE, C64::new(0.2, 0.7)],
        )
        .unwrap();
        near(result.reflected_power_fraction, 1.0, 6e-15);
        assert_eq!(result.transmitted_power_fraction, 0.0);
        assert_eq!(
            result.exit_wave_kind,
            if nt == 1.2 {
                UniaxialWaveKind::Critical
            } else {
                UniaxialWaveKind::Evanescent
            }
        );
        let front = result
            .field_at(SlabRegion::Transmitted, [0.0, 0.0, 0.6])
            .unwrap();
        let further = result
            .field_at(SlabRegion::Transmitted, [0.0, 0.0, 1.1])
            .unwrap();
        let factor = (result.exit_normal_wavenumber_rad_m.im * 0.5).exp();
        for j in 0..3 {
            complex_near(
                further.electric_v_m[j],
                front.electric_v_m[j].scale(factor),
                4e-15,
            );
        }
        assert!(further.poynting_w_m2[2].abs() < 1e-17);
    }
}

#[test]
fn matched_and_thin_slab_limits_keep_phase_and_physical_domains() {
    let d = input(1.3, 1.3, 1.3, 1.3, 0.8, [1.0, 2.0, 3.0]);
    let e = [C64::new(0.4, 0.7), C64::new(-0.2, 0.1)];
    let result = solve(d, Incidence::TangentialIndex(0.5), e).unwrap();
    assert!(result.reflected_power_fraction < 1e-29);
    near(result.transmitted_power_fraction, 1.0, 3e-15);
    for j in 0..2 {
        complex_near(
            result.transmitted_electric_v_m[j],
            e[j] * phase(1.2 * 0.8),
            3e-15,
        );
    }
    let thin = solve(
        input(1.0, 1.5, 1.9, 1.4, 1e-10, [1.0, 2.0, 3.0]),
        Incidence::Cosine(1.0),
        e,
    )
    .unwrap();
    for j in 0..2 {
        complex_near(
            thin.reflected_electric_v_m[j],
            e[j].scale((1.0 - 1.4) / (1.0 + 1.4)),
            5e-10,
        );
        complex_near(
            thin.transmitted_electric_v_m[j],
            e[j].scale(2.0 / (1.0 + 1.4)),
            5e-10,
        );
    }
    assert!(
        result
            .field_at(SlabRegion::Crystal, [0.0, 0.0, -0.1])
            .is_err()
    );
    assert!(
        result
            .field_at(SlabRegion::Transmitted, [0.0, 0.0, 0.7])
            .is_err()
    );
    assert!(
        result
            .field_at(SlabRegion::Incident, [0.0, 0.0, 0.1])
            .is_err()
    );
    assert!(
        solve(
            input(2.0, 1.0, 1.7, 2.0, 0.5, [0.0, 0.0, 1.0]),
            Incidence::TangentialIndex(1.0),
            e
        )
        .is_err()
    );
    assert!(
        solve(
            Input {
                thickness_m: 0.0,
                ..d
            },
            Incidence::Cosine(1.0),
            e
        )
        .is_err()
    );
}

#[test]
fn optical_exit_fields_retain_the_authored_face_distance() {
    let mut d = input(1.0, 1.0, 1.0, 1.0, 1.000_000_000_000_000_2, [1.0, 2.0, 3.0]);
    d.crystal.omega = 1e7 * C0;
    let r = solve(d, Incidence::Cosine(1.0), [C64::ZERO, C64::ONE]).unwrap();
    let field = r
        .field_at(SlabRegion::Transmitted, [0.0, 0.0, 1000.000_000_000_000_1])
        .unwrap();
    // Independent 110-digit reduction for the exact binary64 k*z. A matched
    // slab must preserve that complete physical path across its exit anchor.
    let expected = C64::new(0.873_120_176_906_475_4, 0.487_505_032_465_107_1);
    complex_near(field.electric_v_m[1], expected, 5e-16);
    complex_near(
        field.magnetic_a_m[0],
        expected.scale(-1.0 / (MU * C0)),
        3e-18,
    );
    near(field.poynting_w_m2[2], 0.5 / (MU * C0), 2e-18);
}

#[test]
fn original_slab_exit_retains_a_late_evanescent_field_without_fictitious_heat() {
    let d = input(1.5, 1.7, 2.0, 0.9, 0.6, [1.0, 2.0, 3.0]);
    let r = solve(
        d,
        Incidence::TangentialIndex(1.2),
        [C64::ZERO, C64::from_re(1e100)],
    )
    .unwrap();
    assert_eq!(r.exit_wave_kind, UniaxialWaveKind::Evanescent);
    let z = r.thickness_m - 800.0 / r.exit_normal_wavenumber_rad_m.im;
    let tail = r.field_at(SlabRegion::Transmitted, [0.0, 0.0, z]).unwrap();
    let exposure = -r.exit_normal_wavenumber_rad_m.im * (z - r.thickness_m);
    let expected = (r.transmitted_electric_v_m[1].abs().ln() - exposure).exp();
    assert!(tail.electric_v_m[1].abs() > 0.0);
    near(tail.electric_v_m[1].abs() / expected, 1.0, 4e-13);
    assert_eq!(r.transmitted_power_fraction, 0.0);
    assert_eq!(tail.poynting_w_m2[2], 0.0);
}
