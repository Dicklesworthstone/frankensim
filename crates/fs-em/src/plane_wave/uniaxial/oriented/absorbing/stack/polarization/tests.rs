use super::super::{C0, Frequency, ResolvedLayer};
use super::*;
use crate::plane_wave::uniaxial::oriented::kernel::unit;

const MU: f64 = 1.256_637_062_12e-6;

fn frequency() -> Frequency {
    let mut entrance = 0.0;
    let layers = [
        (
            [0.3, 0.7, 0.5],
            [C64::new(1.5, -0.03), C64::new(1.9, -0.2)],
            MU,
            0.6,
        ),
        (
            [-0.6, 0.4, 0.7],
            [C64::new(1.7, -0.1), C64::new(1.45, -0.025)],
            1.2 * MU,
            0.9,
        ),
    ]
    .into_iter()
    .map(|(axis, indices, mu, thickness)| {
        let exit = entrance + thickness;
        let layer = ResolvedLayer {
            axis: unit(axis).unwrap(),
            indices,
            mu,
            thickness,
            entrance,
            exit,
            loss_f_m: indices.map(|n| -2.0 * n.re * n.im / (mu * C0 * C0)),
        };
        entrance = exit;
        layer
    })
    .collect();
    Frequency {
        omega: C0,
        ni: 1.2,
        nt: 1.4,
        mui: MU,
        mut_: 1.3 * MU,
        layers,
    }
}

fn near(a: f64, b: f64, relative: f64) {
    assert!(
        (a - b).abs() <= relative * b.abs().max(1e-30),
        "{a:e} != {b:e}"
    );
}

#[test]
fn cached_elliptical_layer_heat_matches_direct_stack_solves() {
    let data = frequency();
    for plane in [[1.0, 0.0], [0.6, 0.8]] {
        for cosine in [1.0, 0.7, 1e-6] {
            let prepare = |electric| {
                kernel::solve(
                    &data,
                    Incidence::Cosine(cosine),
                    electric,
                    plane,
                    &mut || false,
                )
                .unwrap_or_else(|error| {
                    panic!("plane={plane:?}, cosine={cosine:e}, electric={electric:?}: {error:?}")
                })
            };
            let basis = [
                prepare([C64::ONE, C64::ZERO]),
                prepare([C64::ZERO, C64::ONE]),
            ];
            for electric in [
                [C64::new(0.3, 0.7), C64::new(-0.5, 0.2)],
                [C64::ONE, C64::new(0.0, 1.0)],
                [C64::from_re(1e-4), C64::from_re(-2.0)],
            ] {
                let expected = prepare(electric);
                near(
                    incident_power([&basis[0], &basis[1]], electric).unwrap(),
                    expected.incident_normal_power_w_m2,
                    3e-12,
                );
                for (j, layer) in expected.layers.iter().enumerate() {
                    let start = 0.125 * layer.thickness_m;
                    let end = 0.875 * layer.thickness_m;
                    near(
                        layer_power(
                            [&basis[0].layers[j], &basis[1].layers[j]],
                            electric,
                            1.7,
                            start,
                            end,
                        )
                        .unwrap(),
                        layer.absorbed_power_w(1.7, start, end).unwrap(),
                        4e-9,
                    );
                }
            }
        }
    }
}

#[test]
fn cached_stack_preserves_weak_absorption_beside_lossless_polarization() {
    let mut data = frequency();
    for layer in &mut data.layers {
        layer.axis = [0.0, 1.0, 0.0];
        layer.indices[0].im = 0.0;
        layer.loss_f_m[0] = 0.0;
    }
    let prepare = |electric| {
        kernel::solve(
            &data,
            Incidence::Cosine(1.0),
            electric,
            [1.0, 0.0],
            &mut || false,
        )
        .unwrap()
    };
    let basis = [
        prepare([C64::ONE, C64::ZERO]),
        prepare([C64::ZERO, C64::ONE]),
    ];
    let electric = [C64::from_re(1e150), C64::from_re(1e-150)];
    assert!(
        incident_power([&basis[0], &basis[1]], electric)
            .unwrap()
            .is_finite()
    );
    for (j, layer) in basis[1].layers.iter().enumerate() {
        let expected = layer.absorbed_power_w(1.0, 0.0, layer.thickness_m).unwrap();
        assert!(expected > 0.0);
        near(
            layer_power(
                [&basis[0].layers[j], layer],
                electric,
                1e300,
                0.0,
                layer.thickness_m,
            )
            .unwrap(),
            expected,
            5e-13,
        );
    }
}

#[test]
fn cached_stack_refuses_invalid_coefficients_and_layer_volume() {
    let data = frequency();
    let prepare = |electric| {
        kernel::solve(
            &data,
            Incidence::Cosine(0.8),
            electric,
            [1.0, 0.0],
            &mut || false,
        )
        .unwrap()
    };
    let basis = [
        prepare([C64::ONE, C64::ZERO]),
        prepare([C64::ZERO, C64::ONE]),
    ];
    let layers = [&basis[0].layers[0], &basis[1].layers[0]];
    for electric in [
        [C64::ZERO; 2],
        [C64::new(f64::NAN, 0.0), C64::ONE],
        [C64::from_re(f64::INFINITY), C64::ONE],
    ] {
        assert!(incident_power([&basis[0], &basis[1]], electric).is_err());
        assert!(layer_power(layers, electric, 1.0, 0.0, 0.6).is_err());
    }
    assert!(incident_power([&basis[0], &basis[1]], [C64::from_re(f64::MAX), C64::ONE]).is_err());
    for (area, start, end) in [(0.0, 0.0, 0.6), (1.0, -0.1, 0.6), (1.0, 0.0, 0.7)] {
        assert!(layer_power(layers, [C64::ONE; 2], area, start, end).is_err());
    }
    let mut wrong = *layers[1];
    wrong.thickness_m = 0.7;
    assert!(layer_power([layers[0], &wrong], [C64::ONE; 2], 1.0, 0.0, 0.6).is_err());
}

fn same_field(actual: AbsorbingCrystalField, expected: AbsorbingCrystalField) {
    for (a, b) in [
        (actual.field.electric_v_m, expected.field.electric_v_m),
        (actual.field.magnetic_a_m, expected.field.magnetic_a_m),
    ] {
        let scale = b
            .iter()
            .fold(1e-30_f64, |s, v| s.max(v.re.abs()).max(v.im.abs()));
        for (a, b) in a.into_iter().zip(b) {
            assert!((a.re - b.re).abs().max((a.im - b.im).abs()) < 2e-8 * scale);
        }
    }
    let scale = expected
        .field
        .poynting_w_m2
        .iter()
        .fold(1e-30_f64, |s, v| s.max(v.abs()));
    for (a, b) in actual
        .field
        .poynting_w_m2
        .into_iter()
        .zip(expected.field.poynting_w_m2)
    {
        assert!((a - b).abs() < 2e-8 * scale);
    }
    near(
        actual.absorbed_power_density_w_m3,
        expected.absorbed_power_density_w_m3,
        2e-8,
    );
}

#[test]
fn g1_cached_stack_fields_keep_elliptical_interference_and_layer_heat() {
    let data = frequency();
    for plane in [[1.0, 0.0], [0.6, 0.8]] {
        for cosine in [1.0, 0.7, 1e-6] {
            let prepare = |electric| {
                kernel::solve(
                    &data,
                    Incidence::Cosine(cosine),
                    electric,
                    plane,
                    &mut || false,
                )
                .unwrap()
            };
            let basis = [
                prepare([C64::ONE, C64::ZERO]),
                prepare([C64::ZERO, C64::ONE]),
            ];
            let electric = [C64::new(0.3, 0.7), C64::new(-0.5, 0.2)];
            let direct = prepare(electric);
            for (region, depth) in [
                (StackRegion::Incident, -0.2),
                (StackRegion::Reflected, -0.2),
                (StackRegion::TotalIncidentSide, -0.2),
                (StackRegion::Layer(0), 0.3),
                (StackRegion::Layer(1), 0.9),
                (StackRegion::Transmitted, 1.8),
            ] {
                let point = [0.17, -0.21, depth];
                same_field(
                    field_at([&basis[0], &basis[1]], electric, region, point).unwrap(),
                    direct.field_at(region, point).unwrap(),
                );
            }
            // Maxwell tangential E/H continuity is an independent physical
            // constraint on the recovered polarization at a shared layer face.
            let left = field_at(
                [&basis[0], &basis[1]],
                electric,
                StackRegion::Layer(0),
                [0.0, 0.0, 0.6],
            )
            .unwrap()
            .field;
            let right = field_at(
                [&basis[0], &basis[1]],
                electric,
                StackRegion::Layer(1),
                [0.0, 0.0, 0.6],
            )
            .unwrap()
            .field;
            for (a, b) in [
                (left.electric_v_m, right.electric_v_m),
                (left.magnetic_a_m, right.magnetic_a_m),
            ] {
                for j in 0..2 {
                    near(a[j].re, b[j].re, 2e-8);
                    near(a[j].im, b[j].im, 2e-8);
                }
            }
        }
    }
}

#[test]
fn g1_isotropic_matched_stack_transports_circular_and_unpolarized_light() {
    let mut data = frequency();
    data.nt = data.ni;
    data.mut_ = data.mui;
    for layer in &mut data.layers {
        layer.indices = [C64::from_re(data.ni); 2];
        layer.mu = data.mui;
        layer.loss_f_m = [0.0; 2];
    }
    let prepare = |electric| {
        kernel::solve(
            &data,
            Incidence::Cosine(1.0),
            electric,
            [1.0, 0.0],
            &mut || false,
        )
        .unwrap()
    };
    let basis = [
        prepare([C64::ONE, C64::ZERO]),
        prepare([C64::ZERO, C64::ONE]),
    ];
    let admittance = data.ni / (MU * C0);
    for (region, z) in [
        (StackRegion::Layer(0), 0.3),
        (StackRegion::Layer(1), 0.9),
        (StackRegion::Transmitted, 1.8),
    ] {
        let circular = field_at(
            [&basis[0], &basis[1]],
            [C64::ONE, C64::new(0.0, 1.0)],
            region,
            [0.0, 0.0, z],
        )
        .unwrap();
        let phase = C64::new((data.ni * z).cos(), -(data.ni * z).sin());
        let expected_y = C64::new(-phase.im, phase.re);
        near(circular.field.electric_v_m[0].re, phase.re, 2e-12);
        near(circular.field.electric_v_m[0].im, phase.im, 2e-12);
        near(circular.field.electric_v_m[1].re, expected_y.re, 2e-12);
        near(circular.field.electric_v_m[1].im, expected_y.im, 2e-12);
        near(circular.field.poynting_w_m2[2], admittance, 2e-12);
        assert_eq!(circular.absorbed_power_density_w_m3, 0.0);
        let unpolarized = field_statistics(
            [&basis[0], &basis[1]],
            &[[C64::ONE, C64::ZERO], [C64::ZERO, C64::ONE]],
            2,
            region,
            [0.0, 0.0, z],
            &mut || false,
        )
        .unwrap();
        near(unpolarized.electric_coherency_v2_m2[0][0].re, 0.5, 2e-12);
        near(unpolarized.electric_coherency_v2_m2[1][1].re, 0.5, 2e-12);
        assert!(unpolarized.electric_coherency_v2_m2[0][1].norm_sq() < 1e-26);
        near(unpolarized.poynting_w_m2[2], admittance, 2e-12);
        assert_eq!(unpolarized.absorbed_power_density_w_m3, 0.0);
    }
}

#[test]
fn g3_stack_partial_polarization_is_factorization_invariant_and_bounded() {
    let data = frequency();
    let prepare = |electric| {
        kernel::solve(
            &data,
            Incidence::Cosine(0.7),
            electric,
            [0.6, 0.8],
            &mut || false,
        )
        .unwrap()
    };
    let basis = [
        prepare([C64::ONE, C64::ZERO]),
        prepare([C64::ZERO, C64::ONE]),
    ];
    let a = [C64::new(0.3, 0.7), C64::new(-0.5, 0.2)];
    let b = [C64::new(0.2, -0.4), C64::new(0.1, 0.6)];
    let imaginary = C64::new(0.0, 1.0);
    let root_half = core::f64::consts::FRAC_1_SQRT_2;
    let mixed = [
        core::array::from_fn(|j| (a[j] + imaginary * b[j]).scale(root_half)),
        core::array::from_fn(|j| (imaginary * a[j] + b[j]).scale(root_half)),
    ];
    let region = StackRegion::Layer(1);
    let point = [0.2, 0.3, 0.9];
    let query = |factors: &[[C64; 2]]| {
        field_statistics(
            [&basis[0], &basis[1]],
            factors,
            2,
            region,
            point,
            &mut || false,
        )
        .unwrap()
    };
    let expected = query(&[a, b]);
    let actual = query(&mixed);
    for (a, b) in actual
        .electric_coherency_v2_m2
        .into_iter()
        .flatten()
        .zip(expected.electric_coherency_v2_m2.into_iter().flatten())
    {
        assert!((a - b).norm_sq() < 1e-26);
    }
    near(
        actual.absorbed_power_density_w_m3,
        expected.absorbed_power_density_w_m3,
        2e-12,
    );
    for (a, b) in actual.poynting_w_m2.into_iter().zip(expected.poynting_w_m2) {
        near(a, b, 2e-12);
    }
    let opposite = query(&[a, a.map(|v| -v)]);
    let single = query(&[a]);
    near(
        opposite.absorbed_power_density_w_m3,
        2.0 * single.absorbed_power_density_w_m3,
        2e-12,
    );
    assert!(opposite.absorbed_power_density_w_m3 > 0.0);
    assert_eq!(query(&[]).poynting_w_m2, [0.0; 3]);
    assert_eq!(query(&[[C64::ZERO; 2]]).absorbed_power_density_w_m3, 0.0);
    assert!(
        field_statistics(
            [&basis[0], &basis[1]],
            &[a, b],
            1,
            region,
            point,
            &mut || false
        )
        .is_err()
    );
    assert!(
        field_statistics(
            [&basis[0], &basis[1]],
            &[],
            0,
            region,
            [0.0; 3],
            &mut || false
        )
        .is_err()
    );
    assert!(
        field_statistics(
            [&basis[0], &basis[1]],
            &[[C64::new(f64::NAN, 0.0); 2]],
            1,
            region,
            point,
            &mut || false
        )
        .is_err()
    );
    assert!(matches!(
        field_statistics(
            [&basis[0], &basis[1]],
            &[a, b],
            2,
            region,
            point,
            &mut || true
        ),
        Err(super::super::AngularSpectrumError::Cancelled)
    ));
    let mut polls = 0;
    assert!(matches!(
        field_statistics(
            [&basis[0], &basis[1]],
            &[a, b],
            2,
            region,
            point,
            &mut || {
                polls += 1;
                polls == 3
            }
        ),
        Err(super::super::AngularSpectrumError::Cancelled)
    ));
}
