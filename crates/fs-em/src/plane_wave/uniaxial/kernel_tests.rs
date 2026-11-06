use super::*;

const MU: f64 = 1.256_637_062_12e-6;

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:.17e} != {expected:.17e}; tolerance={tolerance:.3e}"
    );
}

fn complex_near(actual: C64, expected: C64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:?} != {expected:?}; tolerance={tolerance:.3e}"
    );
}

fn input(ni: f64, no: f64, ne: f64, p: f64, axis: [f64; 2], mode: UniaxialMode) -> Input {
    Input {
        omega_rad_s: core::f64::consts::TAU * LIGHT_SPEED_M_S / 550e-9,
        incident_index: ni,
        incident_mu_h_m: MU,
        ordinary_index: no,
        extraordinary_index: ne,
        crystal_mu_h_m: MU,
        optic_axis: axis,
        incidence: Incidence::TangentialIndex(p),
        mode,
    }
}

fn norm(v: [f64; 2]) -> [f64; 2] {
    let length = v[0].hypot(v[1]);
    [v[0] / length, v[1] / length]
}

fn compare_direction(actual: Option<[f64; 2]>, expected: [f64; 2], tolerance: f64) {
    let actual = actual.unwrap();
    near(actual[0], expected[0], tolerance);
    near(actual[1], expected[1], tolerance);
}

// Independent field reconstruction from the actual permittivity tensor.
// This does not reuse the production dispersion metric, root, or admittance.
fn maxwell_check(data: Input, result: UniaxialInterfaceResponse) {
    let k0 = data.omega_rad_s / LIGHT_SPEED_M_S;
    let p = result.tangential_wavenumber_rad_m / k0;
    let q = result.transmitted_normal_wavenumber_rad_m.scale(1.0 / k0);
    let qi = result.incident_normal_wavenumber_rad_m / k0;
    let ni2 = data.incident_index * data.incident_index;
    let no2 = data.ordinary_index * data.ordinary_index;
    let ne2 = data.extraordinary_index * data.extraordinary_index;
    let axis = norm(data.optic_axis);
    let contrast = ne2 - no2;
    let exx = no2 + contrast * axis[0] * axis[0];
    let ezz = no2 + contrast * axis[1] * axis[1];
    let exz = contrast * axis[0] * axis[1];
    let determinant = exx * ezz - exz * exz;
    let ixx = ezz / determinant;
    let izz = exx / determinant;
    let ixz = -exz / determinant;
    let et = result.transmitted_tangential_e;
    let yi = match data.mode {
        UniaxialMode::Ordinary => qi / (data.incident_mu_h_m * LIGHT_SPEED_M_S),
        UniaxialMode::Extraordinary => ni2 / (qi * data.incident_mu_h_m * LIGHT_SPEED_M_S),
    };
    let ht = result.transmitted_conjugate_h.scale(yi);
    complex_near(et, C64::ONE + result.reflected_tangential_e, 3e-14);
    complex_near(
        result.transmitted_conjugate_h,
        C64::ONE - result.reflected_tangential_e,
        3e-14,
    );
    let (sx, sz) = match data.mode {
        UniaxialMode::Ordinary => {
            complex_near(q * q + C64::from_re(p * p), C64::from_re(no2), 3e-13);
            let hx = -(q * et).scale(1.0 / (data.crystal_mu_h_m * LIGHT_SPEED_M_S));
            let hz = et.scale(p / (data.crystal_mu_h_m * LIGHT_SPEED_M_S));
            complex_near(-hx, ht, 3e-14 * yi.max(1.0));
            (0.5 * (et * hz.conj()).re, -0.5 * (et * hx.conj()).re)
        }
        UniaxialMode::Extraordinary => {
            let material_scale = data.crystal_mu_h_m * LIGHT_SPEED_M_S;
            let ex = (q.scale(ixx) - C64::from_re(ixz * p)) * ht;
            let ez = (q.scale(ixz) - C64::from_re(izz * p)) * ht;
            let ex = ex.scale(material_scale);
            let ez = ez.scale(material_scale);
            complex_near(ex, et, 5e-14);
            // Faraday's law verifies the anisotropic dispersion independently.
            let reconstructed_h = (q * ex - ez.scale(p)).scale(1.0 / material_scale);
            complex_near(reconstructed_h, ht, 5e-13 * yi.max(1.0));
            // D is transverse even though E generally is not.
            let dx = ex.scale(exx) + ez.scale(exz);
            let dz = ex.scale(exz) + ez.scale(ezz);
            complex_near(dx.scale(p) + q * dz, C64::ZERO, 8e-13);
            (-0.5 * (ez * ht.conj()).re, 0.5 * (ex * ht.conj()).re)
        }
    };
    near(sz / (0.5 * yi), result.transmitted_power_fraction, 2e-13);
    near(
        result.reflected_power_fraction + result.transmitted_power_fraction,
        1.0,
        4e-14,
    );
    if result.wave_kind == UniaxialWaveKind::Propagating {
        assert!(sz > 0.0);
        compare_direction(result.phase_direction, norm([p, q.re]), 8e-14);
        compare_direction(result.energy_direction, norm([sx, sz]), 8e-13);
    } else if result.wave_kind == UniaxialWaveKind::Evanescent {
        assert!(q.im < 0.0);
        assert!(result.phase_direction.is_none());
        assert!(result.energy_direction.is_none());
        assert_eq!(result.transmitted_power_fraction, 0.0);
    }
}

#[test]
fn isotropic_limit_matches_independent_fresnel_for_both_modes() {
    for ni in [1.0, 1.5, 2.0] {
        for nt in [1.0, 1.7, 2.3] {
            for cosine in [0.05, 0.3, 0.8, 1.0] {
                let p = ni * (1.0_f64 - cosine * cosine).sqrt();
                for mode in [UniaxialMode::Ordinary, UniaxialMode::Extraordinary] {
                    let mut data = input(ni, nt, nt, p, [0.37, -0.91], mode);
                    data.incidence = Incidence::Cosine(cosine);
                    data.crystal_mu_h_m = 1.9 * MU;
                    let result = solve(data).unwrap();
                    let qi = ni * cosine;
                    let squared = nt * nt - p * p;
                    let qt = if squared < 0.0 {
                        C64::new(0.0, -(-squared).sqrt())
                    } else {
                        C64::from_re(squared.sqrt())
                    };
                    let yi = match mode {
                        UniaxialMode::Ordinary => qi / MU,
                        UniaxialMode::Extraordinary => ni * ni / (qi * MU),
                    };
                    let yt = match mode {
                        UniaxialMode::Ordinary => qt.scale(1.0 / data.crystal_mu_h_m),
                        UniaxialMode::Extraordinary => {
                            C64::from_re(nt * nt / data.crystal_mu_h_m) / qt
                        }
                    };
                    let reflected = (C64::from_re(yi) - yt) / (C64::from_re(yi) + yt);
                    complex_near(result.reflected_tangential_e, reflected, 8e-14);
                    maxwell_check(data, result);
                }
            }
        }
    }
}

#[test]
fn tilted_axis_retains_positive_energy_with_negative_normal_phase() {
    let data = input(2.0, 1.0, 2.0, 1.5, [1.0, 1.0], UniaxialMode::Extraordinary);
    let result = solve(data).unwrap();
    let k0 = data.omega_rad_s / LIGHT_SPEED_M_S;
    near(
        result.transmitted_normal_wavenumber_rad_m.re / k0,
        -0.5,
        2e-15,
    );
    compare_direction(result.phase_direction, norm([1.5, -0.5]), 3e-15);
    compare_direction(result.energy_direction, norm([1.5, 0.5]), 3e-15);
    assert!(result.transmitted_power_fraction > 0.0);
    maxwell_check(data, result);
}

#[test]
fn normal_incidence_can_still_have_extraordinary_walkoff() {
    let data = input(1.0, 1.5, 2.0, 0.0, [1.0, 1.0], UniaxialMode::Extraordinary);
    let result = solve(data).unwrap();
    compare_direction(result.phase_direction, [0.0, 1.0], 0.0);
    // At p=0, Sx/Sz=G/F=(ne²-no²)/(ne²+no²).
    compare_direction(result.energy_direction, norm([1.75 / 6.25, 1.0]), 4e-15);
    maxwell_check(data, result);
}

#[test]
fn aligned_axes_match_analytic_wave_normals_and_walkoff() {
    let (ni, no, ne, p) = (1.5, 1.7, 2.1, 0.8);
    for axis in [[0.0, 1.0], [1.0, 0.0]] {
        let data = input(ni, no, ne, p, axis, UniaxialMode::Extraordinary);
        let result = solve(data).unwrap();
        let (q, energy) = if axis[0] == 0.0 {
            let q = no * (1.0_f64 - p * p / (ne * ne)).sqrt();
            (q, norm([p / (ne * ne), q / (no * no)]))
        } else {
            let q = ne * (1.0_f64 - p * p / (no * no)).sqrt();
            (q, norm([p / (no * no), q / (ne * ne)]))
        };
        let k0 = data.omega_rad_s / LIGHT_SPEED_M_S;
        near(result.transmitted_normal_wavenumber_rad_m.re / k0, q, 3e-15);
        compare_direction(result.energy_direction, energy, 3e-15);
        maxwell_check(data, result);
    }
}

#[test]
fn critical_te_and_tm_keep_distinct_finite_fields_without_infinite_admittance() {
    let ordinary = input(2.0, 1.0, 1.4, 1.0, [0.0, 1.0], UniaxialMode::Ordinary);
    let extraordinary = input(2.0, 1.4, 1.0, 1.0, [0.0, 1.0], UniaxialMode::Extraordinary);
    for (data, r, te, th) in [(ordinary, 1.0, 2.0, 0.0), (extraordinary, -1.0, 0.0, 2.0)] {
        let result = solve(data).unwrap();
        assert_eq!(result.wave_kind, UniaxialWaveKind::Critical);
        assert_eq!(result.transmitted_power_fraction, 0.0);
        assert_eq!(result.reflected_tangential_e, C64::from_re(r));
        assert_eq!(result.transmitted_tangential_e, C64::from_re(te));
        assert_eq!(result.transmitted_conjugate_h, C64::from_re(th));
        compare_direction(result.phase_direction, [1.0, 0.0], 0.0);
        compare_direction(result.energy_direction, [1.0, 0.0], 0.0);
        maxwell_check(data, result);
    }
}

#[test]
fn adjacent_tangential_inputs_preserve_both_sides_of_each_critical_threshold() {
    for mode in [UniaxialMode::Ordinary, UniaxialMode::Extraordinary] {
        for (p, kind) in [
            (1.0_f64.next_down(), UniaxialWaveKind::Propagating),
            (1.0, UniaxialWaveKind::Critical),
            (1.0_f64.next_up(), UniaxialWaveKind::Evanescent),
        ] {
            let mut data = input(2.0, 1.0, 1.0, p, [0.6, 0.8], mode);
            if mode == UniaxialMode::Extraordinary {
                data.ordinary_index = 1.3;
                data.optic_axis = [0.0, 1.0];
            }
            let result = solve(data).unwrap();
            assert_eq!(result.wave_kind, kind);
            if kind == UniaxialWaveKind::Propagating {
                assert!(result.transmitted_power_fraction > 0.0);
            }
            maxwell_check(data, result);
        }
    }
}

#[test]
fn evanescent_extraordinary_field_keeps_tilted_phase_and_reflection_phase() {
    let data = input(2.0, 1.0, 2.0, 1.75, [1.0, 1.0], UniaxialMode::Extraordinary);
    let result = solve(data).unwrap();
    assert_eq!(result.wave_kind, UniaxialWaveKind::Evanescent);
    let k0 = data.omega_rad_s / LIGHT_SPEED_M_S;
    complex_near(
        result.transmitted_normal_wavenumber_rad_m.scale(1.0 / k0),
        C64::new(-1.05, -0.6),
        4e-15,
    );
    let yi = 4.0 / (4.0_f64 - 1.75 * 1.75).sqrt();
    let yt = C64::new(0.0, 2.0 / 0.75);
    let expected = (C64::from_re(yi) - yt) / (C64::from_re(yi) + yt);
    complex_near(result.reflected_tangential_e, expected, 5e-15);
    assert!(result.reflected_tangential_e.im.abs() > 0.1);
    maxwell_check(data, result);
}

#[test]
fn ordinary_and_extraordinary_can_have_different_total_reflection_thresholds() {
    let ordinary = input(2.0, 1.2, 1.8, 1.5, [0.0, 1.0], UniaxialMode::Ordinary);
    let mut extraordinary = ordinary;
    extraordinary.mode = UniaxialMode::Extraordinary;
    let o = solve(ordinary).unwrap();
    let e = solve(extraordinary).unwrap();
    assert_eq!(o.wave_kind, UniaxialWaveKind::Evanescent);
    assert_eq!(e.wave_kind, UniaxialWaveKind::Propagating);
    assert_eq!(o.transmitted_power_fraction, 0.0);
    assert!(e.transmitted_power_fraction > 0.5);
}

#[test]
fn maxwell_tensor_boundary_and_poynting_checks_cover_signed_incidence_and_axis_tilt() {
    for angle in [-1.4_f64, -0.9, -0.35, 0.0, 0.2, 0.7, 1.3] {
        let (s, c) = angle.sin_cos();
        for p in [-1.68, -1.3, -0.7, 0.0, 0.7, 1.3, 1.68] {
            for mode in [UniaxialMode::Ordinary, UniaxialMode::Extraordinary] {
                let mut data = input(1.7, 1.35, 2.1, p, [s, c], mode);
                data.crystal_mu_h_m = 2.3 * MU;
                let result = solve(data).unwrap();
                maxwell_check(data, result);
            }
        }
    }
}

#[test]
fn axis_sign_and_scale_preserve_the_same_physical_crystal() {
    for mode in [UniaxialMode::Ordinary, UniaxialMode::Extraordinary] {
        let mut data = input(1.5, 1.7, 2.0, 1.1, [3.0, -4.0], mode);
        let reference = solve(data).unwrap();
        for axis in [[-3.0, 4.0], [3e150, -4e150], [3e-150, -4e-150]] {
            data.optic_axis = axis;
            let result = solve(data).unwrap();
            complex_near(
                result.reflected_tangential_e,
                reference.reflected_tangential_e,
                5e-15,
            );
            compare_direction(
                result.energy_direction,
                reference.energy_direction.unwrap(),
                5e-15,
            );
        }
    }
}

#[test]
fn matched_media_keep_exact_transmission_at_extreme_grazing() {
    for mode in [UniaxialMode::Ordinary, UniaxialMode::Extraordinary] {
        let mut data = input(1.7, 1.7, 1.7, 0.0, [0.3, 0.7], mode);
        for cosine in [1.0, 0.2, 1e-12, 1e-100, 1e-200] {
            data.incidence = Incidence::Cosine(cosine);
            let result = solve(data).unwrap();
            assert_eq!(result.wave_kind, UniaxialWaveKind::Propagating);
            assert_eq!(result.reflected_tangential_e, C64::ZERO);
            assert_eq!(result.transmitted_tangential_e, C64::ONE);
            assert_eq!(result.transmitted_conjugate_h, C64::ONE);
            assert_eq!(result.reflected_power_fraction, 0.0);
            assert_eq!(result.transmitted_power_fraction, 1.0);
            assert!(result.energy_direction.unwrap()[1] > 0.0);
        }
    }
}

#[test]
fn magnetic_contrast_is_retained_even_when_optical_indices_match() {
    for mode in [UniaxialMode::Ordinary, UniaxialMode::Extraordinary] {
        let mut data = input(1.5, 1.5, 1.5, 0.9, [0.6, 0.8], mode);
        data.crystal_mu_h_m = 2.0 * MU;
        let result = solve(data).unwrap();
        near(result.reflected_tangential_e.re, 1.0 / 3.0, 1e-15);
        near(result.transmitted_power_fraction, 8.0 / 9.0, 1e-15);
        maxwell_check(data, result);
    }
}

#[test]
fn numerical_range_scaling_preserves_refraction_and_weak_power() {
    let data = input(1.5, 1.7, 2.1, 0.9, [0.6, 0.8], UniaxialMode::Extraordinary);
    let reference = solve(data).unwrap();
    for scale in [1e-100, 1e100] {
        let mut scaled = data;
        scaled.incident_index *= scale;
        scaled.ordinary_index *= scale;
        scaled.extraordinary_index *= scale;
        scaled.incidence = Incidence::TangentialIndex(0.9 * scale);
        let result = solve(scaled).unwrap();
        complex_near(
            result.reflected_tangential_e,
            reference.reflected_tangential_e,
            4e-15,
        );
        compare_direction(
            result.energy_direction,
            reference.energy_direction.unwrap(),
            4e-15,
        );
    }
    let mut weak = input(1.0, 1.0, 1.0, 0.5, [1.0, 0.0], UniaxialMode::Extraordinary);
    weak.crystal_mu_h_m = MU * 1e-100;
    let result = solve(weak).unwrap();
    assert_eq!(result.reflected_power_fraction, 1.0);
    assert!(result.transmitted_power_fraction > 0.0);
    near(result.transmitted_power_fraction / 4e-100, 1.0, 4e-15);
}

#[test]
fn malformed_or_unrepresentable_ports_refuse_instead_of_inventing_a_mode() {
    let base = input(1.5, 1.7, 2.1, 0.9, [0.6, 0.8], UniaxialMode::Extraordinary);
    for axis in [[0.0, 0.0], [f64::NAN, 1.0], [1.0, f64::INFINITY]] {
        let mut data = base;
        data.optic_axis = axis;
        assert!(solve(data).is_err());
    }
    for c in [0.0, -0.1, 1.0_f64.next_up(), f64::NAN] {
        let mut data = base;
        data.incidence = Incidence::Cosine(c);
        assert!(solve(data).is_err());
    }
    for p in [-1.5, 1.5, f64::INFINITY, f64::NAN] {
        let mut data = base;
        data.incidence = Incidence::TangentialIndex(p);
        assert!(solve(data).is_err());
    }
    for value in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        let mut data = base;
        data.ordinary_index = value;
        assert!(solve(data).is_err());
        data = base;
        data.crystal_mu_h_m = value;
        assert!(solve(data).is_err());
        data = base;
        data.omega_rad_s = value;
        assert!(solve(data).is_err());
    }
}
