use super::super::{C0, Incidence, Input, UniaxialWaveKind, solve, unit};
use super::*;

const MU: f64 = 1.256_637_062_12e-6;

fn input(axis: [f64; 3], cosine: f64) -> Input {
    Input {
        omega: core::f64::consts::TAU * C0 / 550e-9,
        ni: 1.0,
        no: 1.5,
        ne: 1.9,
        mui: MU,
        mut_: MU,
        axis,
        incidence: Incidence::Cosine(cosine),
        electric: [C64::new(0.7, 0.2), C64::new(-0.3, 0.5)],
    }
}

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:e} != {expected:e}"
    );
}

fn complex_near(actual: C64, expected: C64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:?} != {expected:?}"
    );
}

fn direct_phase(kx: f64, kz: C64, position: [f64; 3]) -> C64 {
    // Separate platform-libm reference, used only for modest phase products.
    let (sine, cosine) = (kx * position[0] + kz.re * position[2]).sin_cos();
    C64::new(cosine, -sine).scale((kz.im * position[2]).exp())
}

#[test]
fn authored_coordinates_propagate_every_component_with_its_actual_wavevector() {
    use OrientedCrystalFieldSide::*;
    let mut data = input([1.0, 0.2, 1.0], 0.8);
    data.ni = 2.0;
    data.no = 1.0;
    data.ne = 2.0;
    data.incidence = Incidence::TangentialIndex(1.5);
    let r = solve(data).unwrap();
    assert!(r.transmitted[1].normal_wavenumber_rad_m.re < 0.0);
    for (side, normal, electric, magnetic, z) in [
        (
            Incident,
            C64::from_re(r.incident_normal_wavenumber_rad_m),
            r.incident_electric_v_m,
            r.incident_magnetic_a_m,
            -0.21e-6,
        ),
        (
            Reflected,
            C64::from_re(-r.incident_normal_wavenumber_rad_m),
            r.reflected_electric_v_m,
            r.reflected_magnetic_a_m,
            -0.21e-6,
        ),
        (
            Ordinary,
            r.transmitted[0].normal_wavenumber_rad_m,
            r.transmitted[0].electric_v_m,
            r.transmitted[0].magnetic_a_m,
            0.13e-6,
        ),
        (
            Extraordinary,
            r.transmitted[1].normal_wavenumber_rad_m,
            r.transmitted[1].electric_v_m,
            r.transmitted[1].magnetic_a_m,
            0.13e-6,
        ),
    ] {
        let point = [0.17e-6, 3.0, z];
        let actual = field(&r, side, point).unwrap();
        let rotation = direct_phase(r.tangential_wavenumber_rad_m, normal, point);
        for j in 0..3 {
            complex_near(actual.electric_v_m[j], electric[j] * rotation, 2e-15);
            complex_near(actual.magnetic_a_m[j], magnetic[j] * rotation, 2e-17);
        }
        assert_eq!(actual, field(&r, side, [point[0], -1e200, z]).unwrap());
        if side == Extraordinary {
            // Negative normal phase advances with z, while energy enters z>0.
            assert!(actual.poynting_w_m2[2] > 0.0);
        }
    }
}

#[test]
fn incident_side_standing_wave_changes_intensity_but_keeps_normal_transport() {
    use OrientedCrystalFieldSide::*;
    let mut data = input([1.0, 2.0, 3.0], 1.0);
    data.no = 2.0;
    data.ne = 2.0;
    data.electric = [C64::ZERO, C64::ONE];
    let r = solve(data).unwrap();
    let at_face = field(&r, TotalIncidentSide, [0.0; 3]).unwrap();
    let quarter = core::f64::consts::FRAC_PI_2 / r.incident_normal_wavenumber_rad_m;
    let away = field(&r, TotalIncidentSide, [0.0, 0.0, -quarter]).unwrap();
    complex_near(at_face.electric_v_m[1], C64::from_re(2.0 / 3.0), 2e-16);
    complex_near(away.electric_v_m[1], C64::new(0.0, 4.0 / 3.0), 5e-16);
    near(
        away.electric_v_m[1].norm_sq() / at_face.electric_v_m[1].norm_sq(),
        4.0,
        4e-15,
    );
    let expected = 0.5 / (MU * C0) * (1.0 - 1.0 / 9.0);
    near(at_face.poynting_w_m2[2], expected, 2e-18);
    near(away.poynting_w_m2[2], expected, 2e-18);
    let reflected = field(&r, Reflected, [0.0, 0.0, -quarter]).unwrap();
    near(reflected.poynting_w_m2[2], -0.5 / (MU * C0 * 9.0), 1e-18);
    let inside = field(&r, Transmitted, [0.0; 3]).unwrap();
    for j in 0..2 {
        complex_near(at_face.electric_v_m[j], inside.electric_v_m[j], 3e-16);
        complex_near(at_face.magnetic_a_m[j], inside.magnetic_a_m[j], 3e-18);
    }
}

#[test]
fn birefringent_retardation_rotates_the_coherent_transmitted_polarization() {
    let mut data = input([1.0, 1.0, 0.0], 1.0);
    data.ne = 2.0;
    data.electric = [C64::ONE, C64::ZERO];
    let r = solve(data).unwrap();
    let ko = r.transmitted[0].normal_wavenumber_rad_m.re;
    let ke = r.transmitted[1].normal_wavenumber_rad_m.re;
    let to = 2.0 * data.ni / (data.ni + data.no);
    let te = 2.0 * data.ni / (data.ni + data.ne);
    let half_beat = core::f64::consts::PI / (ke - ko);
    let expected_power = (data.no * to * to + data.ne * te * te) / (4.0 * MU * C0);
    for z in [0.0, half_beat / 2.0, half_beat] {
        let point = [0.0, 0.0, z];
        let actual = field(&r, OrientedCrystalFieldSide::Transmitted, point).unwrap();
        let ordinary = direct_phase(0.0, C64::from_re(ko), point).scale(0.5 * to);
        let extraordinary = direct_phase(0.0, C64::from_re(ke), point).scale(0.5 * te);
        complex_near(actual.electric_v_m[0], ordinary + extraordinary, 2e-15);
        complex_near(actual.electric_v_m[1], extraordinary - ordinary, 2e-15);
        near(actual.poynting_w_m2[2], expected_power, 3e-18);
    }
    let face = field(&r, OrientedCrystalFieldSide::Transmitted, [0.0; 3]).unwrap();
    let rotated = field(
        &r,
        OrientedCrystalFieldSide::Transmitted,
        [0.0, 0.0, half_beat],
    )
    .unwrap();
    assert!(rotated.electric_v_m[1].norm_sq() > 50.0 * face.electric_v_m[1].norm_sq());
    assert!(rotated.electric_v_m[0].norm_sq() < 0.02 * face.electric_v_m[0].norm_sq());
}

#[test]
fn vector_poynting_retains_cross_mode_interference() {
    use OrientedCrystalFieldSide::*;
    let mut data = input([1.0, 1.0, 1.0], 1.0);
    data.electric = [C64::ONE, C64::ZERO];
    let r = solve(data).unwrap();
    let face = field(&r, Transmitted, [0.0; 3]).unwrap();
    let ordinary = field(&r, Ordinary, [0.0; 3]).unwrap();
    let extraordinary = field(&r, Extraordinary, [0.0; 3]).unwrap();
    near(
        face.poynting_w_m2[2],
        ordinary.poynting_w_m2[2] + extraordinary.poynting_w_m2[2],
        3e-18,
    );
    let cross = (0..2)
        .map(|j| {
            (face.poynting_w_m2[j] - ordinary.poynting_w_m2[j] - extraordinary.poynting_w_m2[j])
                .abs()
        })
        .fold(0.0_f64, f64::max);
    assert!(
        cross > 1e-5,
        "cross-mode transverse transport must survive: {cross:e}"
    );
    let beat = (r.transmitted[1].normal_wavenumber_rad_m.re
        - r.transmitted[0].normal_wavenumber_rad_m.re)
        .abs();
    let shifted = field(
        &r,
        Transmitted,
        [0.0, 0.0, core::f64::consts::FRAC_PI_2 / beat],
    )
    .unwrap();
    near(shifted.poynting_w_m2[2], face.poynting_w_m2[2], 3e-18);
    assert!((0..2).any(|j| (shifted.poynting_w_m2[j] - face.poynting_w_m2[j]).abs() > 1e-5));
}

#[test]
fn spatial_maxwell_curls_hold_in_both_halfspaces_and_with_an_evanescent_mode() {
    use OrientedCrystalFieldSide::*;
    for (ni, no, ne) in [(1.0, 1.5, 1.9), (2.0, 1.0, 1.7)] {
        let data = Input {
            ni,
            no,
            ne,
            ..input([1.0, 2.0, 0.7], 0.7)
        };
        let r = solve(data).unwrap();
        let axis = unit(data.axis).unwrap();
        let h = 550e-9 / 4000.0;
        for (side, z) in [(TotalIncidentSide, -0.11e-6), (Transmitted, 0.11e-6)] {
            let point = [0.17e-6, -0.3e-6, z];
            let at = field(&r, side, point).unwrap();
            let mut de = [[C64::ZERO; 3]; 3];
            let mut dh = de;
            // Independent fourth-order Cartesian finite differences of the
            // physical query, not production modal derivatives or phase code.
            for d in 0..3 {
                for (offset, weight) in [(-2.0, 1.0), (-1.0, -8.0), (1.0, 8.0), (2.0, -1.0)] {
                    let mut p = point;
                    p[d] += offset * h;
                    let f = field(&r, side, p).unwrap();
                    for j in 0..3 {
                        de[d][j] = de[d][j] + f.electric_v_m[j].scale(weight / (12.0 * h));
                        dh[d][j] = dh[d][j] + f.magnetic_a_m[j].scale(weight / (12.0 * h));
                    }
                }
            }
            let mu = if side == Transmitted {
                data.mut_
            } else {
                data.mui
            };
            let epsilon_o = if side == Transmitted {
                no * no
            } else {
                ni * ni
            } / (mu * C0 * C0);
            let epsilon_e = if side == Transmitted {
                ne * ne / (mu * C0 * C0)
            } else {
                epsilon_o
            };
            let projection = (0..3).fold(C64::ZERO, |s, j| s + at.electric_v_m[j].scale(axis[j]));
            for j in 0..3 {
                let a = (j + 1) % 3;
                let b = (j + 2) % 3;
                let displacement = at.electric_v_m[j].scale(epsilon_o)
                    + projection.scale((epsilon_e - epsilon_o) * axis[j]);
                complex_near(
                    de[a][b] - de[b][a],
                    at.magnetic_a_m[j] * C64::new(0.0, -data.omega * mu),
                    0.02,
                );
                complex_near(
                    dh[a][b] - dh[b][a],
                    displacement * C64::new(0.0, data.omega),
                    1e-4,
                );
            }
        }
    }
}

#[test]
fn evanescent_fields_decay_with_real_phase_and_critical_fields_remain_finite() {
    use OrientedCrystalFieldSide::*;
    let data = Input {
        ni: 2.0,
        no: 1.0,
        ne: 1.7,
        ..input([1.0, 2.0, 0.7], 0.7)
    };
    let r = solve(data).unwrap();
    assert_eq!(r.transmitted[0].wave_kind, UniaxialWaveKind::Evanescent);
    assert_eq!(r.transmitted[1].wave_kind, UniaxialWaveKind::Propagating);
    let depth = -1.0 / r.transmitted[0].normal_wavenumber_rad_m.im;
    let face = field(&r, Ordinary, [0.0; 3]).unwrap();
    for exposures in [1.0, 2.0, 10.0] {
        let at = field(&r, Ordinary, [0.0, 0.0, depth * exposures]).unwrap();
        for j in 0..3 {
            complex_near(
                at.electric_v_m[j],
                face.electric_v_m[j].scale((-exposures).exp()),
                3e-16,
            );
        }
        near(at.poynting_w_m2[2], 0.0, 1e-18);
        let total = field(&r, Transmitted, [0.0, 0.0, depth * exposures]).unwrap();
        near(
            total.poynting_w_m2[2],
            r.incident_normal_power_w_m2 * r.transmitted_power_fraction,
            2e-17,
        );
    }
    let critical = solve(Input {
        axis: [0.0, 1.0, 0.0],
        incidence: Incidence::TangentialIndex(data.ne),
        ..data
    })
    .unwrap();
    assert_eq!(
        critical.transmitted[1].wave_kind,
        UniaxialWaveKind::Critical
    );
    let face = field(&critical, Extraordinary, [0.0; 3]).unwrap();
    assert_eq!(
        face,
        field(&critical, Extraordinary, [0.0, 0.0, 3.0]).unwrap()
    );
    near(face.poynting_w_m2[2], 0.0, 1e-18);
}

#[test]
fn deep_evanescence_preserves_scaled_tails_and_bypasses_extinguished_phase() {
    use OrientedCrystalFieldSide::*;
    let data = Input {
        ni: 2.0,
        no: 1.0,
        ne: 1.7,
        electric: [C64::new(1e100, 0.0), C64::new(0.0, 1e100)],
        ..input([1.0, 2.0, 0.7], 0.1)
    };
    let r = solve(data).unwrap();
    let z = -800.0 / r.transmitted[0].normal_wavenumber_rad_m.im;
    let tail = field(&r, Ordinary, [0.0, 0.0, z]).unwrap();
    let original = r.transmitted[0].electric_v_m[0].abs();
    let expected = (original.ln() + r.transmitted[0].normal_wavenumber_rad_m.im * z).exp();
    assert!(tail.electric_v_m[0].abs() > 0.0);
    near(tail.electric_v_m[0].abs() / expected, 1.0, 3e-13);
    let gone = field(&r, Transmitted, [f64::MAX; 3]).unwrap();
    assert_eq!(gone.electric_v_m, [C64::ZERO; 3]);
    assert_eq!(gone.magnetic_a_m, [C64::ZERO; 3]);
    assert_eq!(gone.poynting_w_m2, [0.0; 3]);
}

#[test]
fn optical_product_residuals_and_shared_carrier_keep_coherent_phase() {
    use OrientedCrystalFieldSide::*;
    let matched = solve(Input {
        omega: 1e7 * C0,
        ni: 1.0,
        no: 1.0,
        ne: 1.0,
        electric: [C64::ONE, C64::ZERO],
        ..input([1.0, 2.0, 3.0], 1.0)
    })
    .unwrap();
    assert_eq!(matched.transmitted[0].normal_wavenumber_rad_m.re, 1e7);
    // Independent 110-digit decimal values for exact binary64 k*z. These
    // distinguish product compensation from a prematurely rounded phase.
    for (z, expected) in [
        (
            1.000_000_000_000_000_2,
            C64::new(-0.907_270_387_115_543_2, -0.420_547_791_176_237_53),
        ),
        (
            1000.000_000_000_000_1,
            C64::new(0.873_120_176_906_475_4, 0.487_505_032_465_107_1),
        ),
    ] {
        complex_near(
            field(&matched, Transmitted, [0.0, 0.0, z])
                .unwrap()
                .electric_v_m[0],
            expected,
            2e-16,
        );
        complex_near(
            field(&matched, Incident, [0.0, 0.0, -z])
                .unwrap()
                .electric_v_m[0],
            expected.conj(),
            2e-16,
        );
    }
    let r = solve(input([1.0, 2.0, 0.7], 0.7)).unwrap();
    let close = field(&r, Transmitted, [0.0, 0.0, 0.4e-6]).unwrap();
    let distant = field(&r, Transmitted, [1e200, 0.0, 0.4e-6]).unwrap();
    assert_eq!(close.poynting_w_m2, distant.poynting_w_m2);
    for a in 0..3 {
        for b in 0..3 {
            complex_near(
                close.electric_v_m[a] * close.electric_v_m[b].conj(),
                distant.electric_v_m[a] * distant.electric_v_m[b].conj(),
                5e-16,
            );
        }
    }
}

#[test]
fn wrong_halfspaces_and_nonfinite_coordinates_refuse_without_changing_preparation() {
    use OrientedCrystalFieldSide::*;
    let r = solve(input([1.0, 2.0, 0.7], 0.7)).unwrap();
    for (side, bad_z) in [
        (Incident, 1e-6),
        (Reflected, 1e-6),
        (TotalIncidentSide, 1e-6),
        (Ordinary, -1e-6),
        (Extraordinary, -1e-6),
        (Transmitted, -1e-6),
    ] {
        assert!(field(&r, side, [0.0, 0.0, bad_z]).is_err());
    }
    let before = field(&r, Transmitted, [0.0, 0.0, 0.4e-6]).unwrap();
    for dimension in 0..3 {
        for value in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            let mut point = [0.0; 3];
            point[dimension] = value;
            assert!(field(&r, Transmitted, point).is_err());
        }
    }
    assert!(field(&r, Incident, [0.0, 0.0, -f64::MAX]).is_err());
    assert_eq!(before, field(&r, Transmitted, [0.0, 0.0, 0.4e-6]).unwrap());
}

#[test]
fn shared_accumulator_keeps_distinct_direction_dark_and_longitudinal_fields() {
    let data = Input {
        omega: C0,
        ni: 1.0, no: 1.0, ne: 1.0,
        electric: [C64::ZERO, C64::ONE],
        ..input([1.0, 2.0, 3.0], 0.8)
    };
    let a = solve(Input { incidence: Incidence::TangentialIndex(0.6), ..data }).unwrap();
    let b = solve(Input { incidence: Incidence::TangentialIndex(-0.6), ..data }).unwrap();
    let sample = |x| {
        let mut sum = FieldSum::new(OrientedCrystalFieldSide::Transmitted, [x, 0.0, 0.3]).unwrap();
        sum.add_response(&a, [1.0, 0.0]).unwrap();
        sum.add_response(&b, [1.0, 0.0]).unwrap();
        sum.finish().unwrap()
    };
    let bright = sample(0.0);
    near(bright.poynting_w_m2[2], 4.0 * 0.8 / (2.0 * MU * C0), 4e-18);
    let dark = sample(core::f64::consts::PI / 1.2);
    assert!(dark.electric_v_m[1].abs() < 1e-15);
    assert!(dark.poynting_w_m2[2].abs() < 1e-30);
    near(dark.magnetic_a_m[2].abs(), 1.2 / (MU * C0), 4e-18);
}

#[test]
fn shared_accumulator_retains_sub_ulp_relative_direction_phase() {
    let p = 0.6_f64;
    let other = p.next_up();
    let data = Input {
        omega: C0,
        ni: 1.0, no: 1.0, ne: 1.0,
        electric: [C64::ZERO, C64::ONE],
        incidence: Incidence::TangentialIndex(p),
        ..input([1.0, 2.0, 3.0], 0.8)
    };
    let a = solve(data).unwrap();
    let b = solve(Input {
        incidence: Incidence::TangentialIndex(other),
        electric: [C64::ZERO, -C64::ONE],
        ..data
    }).unwrap();
    let x = 1e-4;
    let mut sum = FieldSum::new(OrientedCrystalFieldSide::Incident, [x, 0.0, 0.0]).unwrap();
    sum.add_response(&a, [1.0, 0.0]).unwrap();
    sum.add_response(&b, [1.0, 0.0]).unwrap();
    let actual = sum.finish().unwrap();
    // exp(-i*p*x)*(1-exp(-i*delta*x)); its omitted quadratic term is <1e-40.
    // Separately rounded C64 phasors lose much of this ~1e-20 physical field.
    let delta = (other - p) * x;
    let expected = direct_phase(p, C64::ZERO, [x, 0.0, 0.0]) * C64::new(0.0, delta);
    assert!(actual.electric_v_m[1].abs() > 0.0);
    complex_near(actual.electric_v_m[1], expected, 3e-35);
}

#[test]
fn shared_accumulator_uses_the_global_incidence_azimuth_and_all_three_coordinates() {
    let plane = [0.6, 0.8];
    let authored = [C64::new(0.8, 0.3), C64::new(-0.2, 0.7)];
    let p = -0.7;
    let data = Input {
        omega: C0,
        ni: 1.5, no: 1.5, ne: 1.5,
        electric: [authored[0].scale(plane[0]) + authored[1].scale(plane[1]),
            authored[1].scale(plane[0]) - authored[0].scale(plane[1])],
        incidence: Incidence::TangentialIndex(p),
        ..input([1.0, 2.0, 3.0], 0.8)
    };
    let mut r = solve(data).unwrap();
    let rotate = |v: Vector| [v[0].scale(plane[0]) - v[1].scale(plane[1]),
        v[0].scale(plane[1]) + v[1].scale(plane[0]), v[2]];
    for wave in &mut r.transmitted {
        wave.electric_v_m = rotate(wave.electric_v_m);
        wave.magnetic_a_m = rotate(wave.magnetic_a_m);
    }
    let point = [0.3, -0.5, 0.7];
    let mut sum = FieldSum::new(OrientedCrystalFieldSide::Transmitted, point).unwrap();
    sum.add_response(&r, plane).unwrap();
    let actual = sum.finish().unwrap();
    // Independent homogeneous-medium transversality and Faraday law.
    let q = (1.5_f64 * 1.5 - p * p).sqrt();
    let k = [p * plane[0], p * plane[1], q];
    let e = [authored[0], authored[1],
        -(authored[0].scale(k[0]) + authored[1].scale(k[1])).scale(1.0 / q)];
    let h = super::super::cross(k.map(C64::from_re), e).map(|v| v.scale(1.0 / (C0 * MU)));
    let (sine, cosine) = (k[0] * point[0] + k[1] * point[1] + k[2] * point[2]).sin_cos();
    let phase = C64::new(cosine, -sine);
    for j in 0..3 {
        complex_near(actual.electric_v_m[j], e[j] * phase, 4e-15);
        complex_near(actual.magnetic_a_m[j], h[j] * phase, 4e-17);
    }
}

#[test]
fn already_phased_fields_share_the_coherent_sum_without_halfspace_admission() {
    let electric = [C64::ZERO, C64::ONE, C64::ZERO];
    let magnetic = [-C64::ONE, C64::ZERO, C64::ZERO];
    let mut sum = FieldSum::from_phasors();
    sum.add_phasor(electric, magnetic, C64::ONE);
    sum.add_phasor(electric, magnetic, C64::new(0.0, 1.0));
    let field = sum.finish().unwrap();
    assert_eq!(field.electric_v_m[1], C64::new(1.0, 1.0));
    assert_eq!(field.magnetic_a_m[0], C64::new(-1.0, -1.0));
    assert_eq!(field.poynting_w_m2, [0.0, 0.0, 1.0]);

    let mut invalid = FieldSum::from_phasors();
    invalid.add_phasor(electric, magnetic, C64::new(f64::NAN, 0.0));
    assert!(invalid.finish().is_err());
}
