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
