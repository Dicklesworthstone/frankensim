use super::super::{absorption::integrated_heat, kernel::solve};
use super::*;
use crate::plane_wave::uniaxial::oriented::{
    absorbing::kernel::Input, kernel::unit, slab::SlabRegion,
};

const MU: f64 = 1.256_637_062_12e-6;

fn input(axis: [f64; 3], indices: [C64; 2], cosine: f64) -> Input {
    Input {
        omega: C0,
        ni: 1.2,
        mui: MU,
        mut_: MU,
        axis,
        indices,
        incidence: Incidence::Cosine(cosine),
        electric: [C64::ONE, C64::ZERO],
        loss_f_m: indices.map(|n| -2.0 * n.re * n.im / (MU * C0 * C0)),
    }
}

fn near(a: f64, b: f64, relative: f64) {
    assert!(
        (a - b).abs() <= relative * b.abs().max(1e-30),
        "{a:e} != {b:e}"
    );
}

#[test]
fn cached_elliptical_response_matches_direct_maxwell_fields_and_heat() {
    let global_axis = unit([0.3, 0.7, 0.5]).unwrap();
    for plane in [[1.0, 0.0], [0.6, 0.8]] {
        let axis = [
            global_axis[0] * plane[0] + global_axis[1] * plane[1],
            global_axis[1] * plane[0] - global_axis[0] * plane[1],
            global_axis[2],
        ];
        for cosine in [1.0, 0.7, 1e-6] {
            let data = input(axis, [C64::new(1.5, -0.03), C64::new(1.9, -0.2)], cosine);
            let local = |v: [C64; 2]| {
                [
                    v[0].scale(plane[0]) + v[1].scale(plane[1]),
                    v[1].scale(plane[0]) - v[0].scale(plane[1]),
                ]
            };
            let prepare = |electric| {
                let mut response = solve(
                    Input {
                        electric: local(electric),
                        ..data
                    },
                    1.4,
                    MU,
                    0.9,
                )
                .unwrap();
                response.rotate(plane, global_axis).unwrap();
                response
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
                let actual = combine_responses([&basis[0], &basis[1]], electric).unwrap();
                let expected = prepare(electric);
                near(
                    actual.incident_normal_power_w_m2,
                    expected.incident_normal_power_w_m2,
                    3e-12,
                );
                near(
                    actual.reflected_power_fraction,
                    expected.reflected_power_fraction,
                    3e-12,
                );
                near(
                    actual.transmitted_power_fraction,
                    expected.transmitted_power_fraction,
                    3e-9,
                );
                for (region, z) in [
                    (SlabRegion::TotalIncidentSide, -0.2),
                    (SlabRegion::Crystal, 0.3),
                    (SlabRegion::Crystal, 0.8),
                    (SlabRegion::Transmitted, 1.2),
                ] {
                    let a = actual.field_at(region, [0.2, -0.1, z]).unwrap();
                    let b = expected.field_at(region, [0.2, -0.1, z]).unwrap();
                    let scale = b
                        .field
                        .electric_v_m
                        .iter()
                        .map(|e| e.abs())
                        .fold(0.0_f64, f64::max);
                    for j in 0..3 {
                        assert!(
                            (a.field.electric_v_m[j] - b.field.electric_v_m[j]).abs()
                                <= 2e-9 * scale
                        );
                    }
                    if region == SlabRegion::Crystal {
                        near(
                            a.absorbed_power_density_w_m3,
                            b.absorbed_power_density_w_m3,
                            3e-9,
                        );
                    }
                }
                near(
                    integrated_heat(&actual, 1.7, 0.1, 0.8).unwrap(),
                    integrated_heat(&expected, 1.7, 0.1, 0.8).unwrap(),
                    3e-9,
                );
            }
        }
    }
}

#[test]
fn cached_polarization_preserves_weak_orthogonal_tensor_heat() {
    let data = input(
        [0.0, 1.0, 0.0],
        [C64::from_re(1.4), C64::new(1.8, -0.1)],
        1.0,
    );
    let basis = [
        solve(data, 1.0, MU, 0.7).unwrap(),
        solve(
            Input {
                electric: [C64::ZERO, C64::ONE],
                ..data
            },
            1.0,
            MU,
            0.7,
        )
        .unwrap(),
    ];
    let result = combine_responses(
        [&basis[0], &basis[1]],
        [C64::from_re(1e150), C64::from_re(1e-150)],
    )
    .unwrap();
    let expected = integrated_heat(&basis[1], 1.0, 0.0, 0.7).unwrap();
    assert!(expected > 0.0);
    near(
        integrated_heat(&result, 1e300, 0.0, 0.7).unwrap(),
        expected,
        4e-13,
    );
}

#[test]
fn cached_polarization_refuses_invalid_inputs_and_incompatible_modes() {
    let data = input(
        [0.3, 0.7, 0.5],
        [C64::new(1.5, -0.03), C64::new(1.9, -0.2)],
        0.8,
    );
    let a = solve(data, 1.0, MU, 0.7).unwrap();
    let b = solve(
        Input {
            electric: [C64::ZERO, C64::ONE],
            ..data
        },
        1.0,
        MU,
        0.7,
    )
    .unwrap();
    for electric in [
        [C64::ZERO; 2],
        [C64::new(f64::NAN, 0.0), C64::ONE],
        [C64::from_re(f64::MAX), C64::ONE],
    ] {
        assert!(combine_responses([&a, &b], electric).is_err());
    }
    let mut wrong = b;
    wrong.thickness_m = 0.8;
    assert!(combine_responses([&a, &wrong], [C64::ONE; 2]).is_err());
}
