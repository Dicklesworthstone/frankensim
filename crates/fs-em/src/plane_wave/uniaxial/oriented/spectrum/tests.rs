use super::*;

const MU: f64 = 1.256_637_062_12e-6;
const C0: f64 = 299_792_458.0;
fn frequency(ni: f64, no: f64, ne: f64) -> ResolvedFrequency {
    ResolvedFrequency {
        omega: C0,
        ni,
        no,
        ne,
        mui: MU,
        mut_: MU,
    }
}
fn mode(p: f64, axis: [f64; 2], e: [C64; 2]) -> UniaxialIncidentMode {
    UniaxialIncidentMode {
        tangential_phase_index: p,
        incidence_plane_axis: axis,
        electric_v_m: e,
    }
}
fn near(x: f64, y: f64, t: f64) {
    assert!((x - y).abs() <= t, "{x:.17e} != {y:.17e}; tolerance={t}");
}
fn complex_near(x: C64, y: C64, t: f64) {
    assert!((x - y).abs() <= t, "{x:?} != {y:?}; tolerance={t}");
}
fn cross(a: [C64; 3], b: [C64; 3]) -> [C64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
fn phase(theta: f64) -> C64 {
    C64::new(theta.cos(), -theta.sin())
}

#[test]
fn single_mode_agrees_with_original_boundary_fields_and_signed_phase() {
    let data = frequency(1.0, 1.5, 1.9);
    let axis = [1.0, 2.0, 3.0];
    let electric = [C64::new(0.7, 0.2), C64::new(-0.3, 0.5)];
    let direct = data
        .solve(axis, Incidence::TangentialIndex(-0.6), electric)
        .unwrap();
    let modes = prepare_modes(data, axis, &[mode(-0.6, [1.0, 0.0], electric)], &mut || {
        false
    })
    .unwrap();
    for (side, e, h, kz) in [
        (
            AngularSpectrumSide::Incident,
            direct.incident_electric_v_m,
            direct.incident_magnetic_a_m,
            direct.incident_normal_wavenumber_rad_m,
        ),
        (
            AngularSpectrumSide::Reflected,
            direct.reflected_electric_v_m,
            direct.reflected_magnetic_a_m,
            -direct.incident_normal_wavenumber_rad_m,
        ),
    ] {
        let point = [0.4, -0.2, -0.7];
        let actual = field_at(&modes, side, point, &mut || false).unwrap();
        let expected_phase = phase(-0.6 * point[0] + kz * point[2]);
        for j in 0..3 {
            complex_near(actual.electric_v_m[j], e[j] * expected_phase, 4e-15);
            complex_near(actual.magnetic_a_m[j], h[j] * expected_phase, 4e-17);
        }
    }
    let actual = field_at(
        &modes,
        AngularSpectrumSide::Transmitted,
        [0.0; 3],
        &mut || false,
    )
    .unwrap();
    for j in 0..3 {
        complex_near(
            actual.electric_v_m[j],
            direct.transmitted[0].electric_v_m[j] + direct.transmitted[1].electric_v_m[j],
            2e-16,
        );
        complex_near(
            actual.magnetic_a_m[j],
            direct.transmitted[0].magnetic_a_m[j] + direct.transmitted[1].magnetic_a_m[j],
            2e-18,
        );
    }
}

#[test]
fn arbitrary_azimuth_reduces_to_exact_isotropic_maxwell_fields() {
    // Independent homogeneous-medium transversality and k cross E/(omega mu).
    let data = frequency(1.5, 1.5, 1.5);
    let input = mode(-0.7, [3.0, 4.0], [C64::new(0.8, 0.3), C64::new(-0.2, 0.7)]);
    let modes = prepare_modes(data, [1.0, 2.0, 3.0], &[input], &mut || false).unwrap();
    let q = (1.5_f64 * 1.5 - 0.7 * 0.7).sqrt();
    let k = [
        C64::from_re(-0.7 * 0.6),
        C64::from_re(-0.7 * 0.8),
        C64::from_re(q),
    ];
    let e = [
        input.electric_v_m[0],
        input.electric_v_m[1],
        -(input.electric_v_m[0] * k[0] + input.electric_v_m[1] * k[1]).scale(1.0 / q),
    ];
    let h = cross(k, e).map(|z| z.scale(1.0 / (C0 * MU)));
    for point in [[0.0; 3], [0.3, -0.5, 0.7]] {
        let value = field_at(&modes, AngularSpectrumSide::Transmitted, point, &mut || {
            false
        })
        .unwrap();
        let factor = phase(k[0].re * point[0] + k[1].re * point[1] + q * point[2]);
        for j in 0..3 {
            complex_near(value.electric_v_m[j], e[j] * factor, 4e-15);
            complex_near(value.magnetic_a_m[j], h[j] * factor, 4e-17);
        }
    }
    let reflected = field_at(
        &modes,
        AngularSpectrumSide::Reflected,
        [0.2, -0.4, -0.7],
        &mut || false,
    )
    .unwrap();
    assert!(reflected.electric_v_m.iter().all(|z| z.abs() < 1e-15));
}

#[test]
fn fixed_global_axis_preserves_cross_polarization_maxwell_and_interface_continuity() {
    let data = frequency(1.0, 1.5, 1.9);
    let axis = [1.0, 1.0, 0.0];
    let a = axis.map(|v| v / 2.0_f64.sqrt());
    let modes = prepare_modes(
        data,
        axis,
        &[
            mode(0.0, [1.0, 0.0], [C64::ONE, C64::ZERO]),
            mode(-0.4, [3.0, -4.0], [C64::new(0.2, 0.7), C64::new(-0.3, 0.1)]),
        ],
        &mut || false,
    )
    .unwrap();
    // At normal incidence the 45-degree crystal axes give the independently
    // known off-diagonal reflection (r_extraordinary-r_ordinary)/2.
    complex_near(
        modes[0].response.reflected_electric_v_m[1],
        C64::from_re(0.5 * ((1.0 - data.ne) / (1.0 + data.ne) - (1.0 - data.no) / (1.0 + data.no))),
        3e-15,
    );
    for source in &modes {
        let r = &source.response;
        for wave in &r.transmitted {
            let k = [
                C64::from_re(r.tangential_wavenumber_rad_m * source.axis[0]),
                C64::from_re(r.tangential_wavenumber_rad_m * source.axis[1]),
                wave.normal_wavenumber_rad_m,
            ];
            let dot = (0..3).fold(C64::ZERO, |v, j| v + wave.electric_v_m[j].scale(a[j]));
            let d: [C64; 3] = core::array::from_fn(|j| {
                (wave.electric_v_m[j].scale(data.no * data.no)
                    + dot.scale((data.ne * data.ne - data.no * data.no) * a[j]))
                .scale(1.0 / (MU * C0 * C0))
            });
            let faraday = cross(k, wave.electric_v_m);
            let ampere = cross(k, wave.magnetic_a_m);
            for j in 0..3 {
                complex_near(faraday[j], wave.magnetic_a_m[j].scale(C0 * MU), 3e-15);
                complex_near(ampere[j], d[j].scale(-C0), 3e-17);
            }
        }
    }
    for point in [[0.0; 3], [0.3, -0.5, 0.0], [-1.1, 0.2, 0.0]] {
        let left = field_at(
            &modes,
            AngularSpectrumSide::TotalIncidentSide,
            point,
            &mut || false,
        )
        .unwrap();
        let right = field_at(&modes, AngularSpectrumSide::Transmitted, point, &mut || {
            false
        })
        .unwrap();
        for j in 0..2 {
            complex_near(left.electric_v_m[j], right.electric_v_m[j], 3e-15);
            complex_near(left.magnetic_a_m[j], right.magnetic_a_m[j], 3e-17);
        }
        near(left.poynting_w_m2[2], right.poynting_w_m2[2], 3e-17);
    }
}

#[test]
fn crystal_depth_produces_birefringent_quarter_wave_beating_without_power_loss() {
    let data = frequency(1.0, 1.5, 2.0);
    let modes = prepare_modes(
        data,
        [1.0, 0.0, 0.0],
        &[mode(0.0, [1.0, 0.0], [C64::ONE, C64::ONE])],
        &mut || false,
    )
    .unwrap();
    // Normal-incidence Fresnel amplitudes and independent retardance (ne-no)z.
    let ex = 2.0 / 3.0;
    let ey = 2.0 / 2.5;
    let power = (2.0 * ex * ex + 1.5 * ey * ey) / (2.0 * MU * C0);
    for z in [
        0.0,
        core::f64::consts::PI,
        2.0 * core::f64::consts::PI,
        4.0 * core::f64::consts::PI,
    ] {
        let f = field_at(
            &modes,
            AngularSpectrumSide::Transmitted,
            [0.0, 0.0, z],
            &mut || false,
        )
        .unwrap();
        complex_near(f.electric_v_m[0], phase(2.0 * z).scale(ex), 2e-15);
        complex_near(f.electric_v_m[1], phase(1.5 * z).scale(ey), 2e-15);
        complex_near(
            f.electric_v_m[0] * f.electric_v_m[1].conj(),
            phase(0.5 * z).scale(ex * ey),
            2e-15,
        );
        near(f.poynting_w_m2[2], power, 4e-18);
    }
}

#[test]
fn distinct_incident_directions_keep_bright_dark_and_longitudinal_interference() {
    let data = frequency(1.0, 1.0, 1.0);
    let p = 0.6;
    let modes = prepare_modes(
        data,
        [1.0, 2.0, 3.0],
        &[
            mode(p, [1.0, 0.0], [C64::ZERO, C64::ONE]),
            mode(-p, [1.0, 0.0], [C64::ZERO, C64::ONE]),
        ],
        &mut || false,
    )
    .unwrap();
    let bright = field_at(
        &modes,
        AngularSpectrumSide::Transmitted,
        [0.0, 0.0, 0.3],
        &mut || false,
    )
    .unwrap();
    let one_mode_normal_power = 0.8 / (2.0 * MU * C0);
    near(bright.poynting_w_m2[2], 4.0 * one_mode_normal_power, 8e-18);
    let dark = field_at(
        &modes,
        AngularSpectrumSide::Transmitted,
        [core::f64::consts::PI / (2.0 * p), 0.0, 0.3],
        &mut || false,
    )
    .unwrap();
    assert!(dark.electric_v_m[1].abs() < 1e-15);
    assert!(dark.poynting_w_m2[2].abs() < 1e-30);
    near(dark.magnetic_a_m[2].abs(), 2.0 * p / (MU * C0), 4e-18);
}

#[test]
fn evanescent_modes_decay_and_unphysical_half_space_queries_refuse() {
    let data = frequency(2.0, 1.0, 1.4);
    let p = 1.7;
    let modes = prepare_modes(
        data,
        [0.0, 0.0, 1.0],
        &[mode(p, [3.0, 4.0], [C64::from_re(-0.8), C64::from_re(0.6)])],
        &mut || false,
    )
    .unwrap();
    let boundary = field_at(
        &modes,
        AngularSpectrumSide::Transmitted,
        [0.0; 3],
        &mut || false,
    )
    .unwrap();
    let depth = 0.7;
    let inside = field_at(
        &modes,
        AngularSpectrumSide::Transmitted,
        [0.0, 0.0, depth],
        &mut || false,
    )
    .unwrap();
    let attenuation = (-((p * p - 1.0).sqrt()) * depth).exp();
    for j in 0..3 {
        complex_near(
            inside.electric_v_m[j],
            boundary.electric_v_m[j].scale(attenuation),
            4e-15,
        );
        complex_near(
            inside.magnetic_a_m[j],
            boundary.magnetic_a_m[j].scale(attenuation),
            4e-17,
        );
    }
    assert!(inside.poynting_w_m2[2].abs() < 1e-18);
    assert!(
        field_at(
            &modes,
            AngularSpectrumSide::Transmitted,
            [0.0, 0.0, -0.1],
            &mut || false
        )
        .is_err()
    );
    assert!(
        field_at(
            &modes,
            AngularSpectrumSide::Incident,
            [0.0, 0.0, 0.1],
            &mut || false
        )
        .is_err()
    );
    let far = field_at(
        &modes,
        AngularSpectrumSide::Transmitted,
        [0.0, 0.0, 1000.0],
        &mut || false,
    )
    .unwrap();
    assert_eq!(far.electric_v_m, [C64::ZERO; 3]);
}

#[test]
fn resource_source_support_and_cancellation_have_no_partial_field() {
    let limits = UniaxialSpectrumLimits {
        maximum_modes: 2,
        maximum_workspace_bytes: usize::MAX,
    };
    assert!(admit(3, limits).is_err());
    assert!(
        admit(
            2,
            UniaxialSpectrumLimits {
                maximum_workspace_bytes: 1,
                ..limits
            }
        )
        .is_err()
    );
    assert!(
        admit(
            usize::MAX,
            UniaxialSpectrumLimits {
                maximum_modes: usize::MAX,
                ..limits
            }
        )
        .is_err()
    );
    let data = frequency(1.0, 1.5, 1.9);
    assert!(
        prepare_modes(
            data,
            [1.0, 2.0, 3.0],
            &[mode(1.0, [1.0, 0.0], [C64::ONE, C64::ZERO])],
            &mut || false
        )
        .is_err()
    );
    let empty = prepare_modes(
        data,
        [1.0, 2.0, 3.0],
        &[mode(10.0, [1.0, 0.0], [C64::ZERO; 2])],
        &mut || false,
    )
    .unwrap();
    assert!(empty.is_empty());
    assert!(matches!(
        field_at(
            &empty,
            AngularSpectrumSide::Transmitted,
            [0.0; 3],
            &mut || true
        ),
        Err(AngularSpectrumError::Cancelled)
    ));
    let modes = prepare_modes(
        data,
        [1.0, 2.0, 3.0],
        &[mode(0.6, [1.0, 0.0], [C64::ONE, C64::ZERO])],
        &mut || false,
    )
    .unwrap();
    let mut polls = 0;
    assert!(matches!(
        field_at(
            &modes,
            AngularSpectrumSide::Transmitted,
            [0.0; 3],
            &mut || {
                polls += 1;
                polls == 2
            }
        ),
        Err(AngularSpectrumError::Cancelled)
    ));
    let recovered = field_at(
        &modes,
        AngularSpectrumSide::Transmitted,
        [0.0; 3],
        &mut || false,
    )
    .unwrap();
    assert!(recovered.poynting_w_m2[2] > 0.0);
}
