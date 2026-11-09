use super::super::super::{Incidence, ResolvedFrequency};
use super::super::kernel::{Input, solve, solve_thickness};
use super::*;

const C0: f64 = 299_792_458.0;
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

fn close(actual: C64, expected: C64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual={actual:?}, expected={expected:?}, tolerance={tolerance:.4e}"
    );
}

fn all_fields(r: UniaxialSlabResponse) -> [[C64; 6]; 6] {
    let join = |e: Vector, h: Vector| [e[0], e[1], e[2], h[0], h[1], h[2]];
    [
        join(r.reflected_electric_v_m, r.reflected_magnetic_a_m),
        join(r.transmitted_electric_v_m, r.transmitted_magnetic_a_m),
        join(r.forward[0].electric_v_m, r.forward[0].magnetic_a_m),
        join(r.forward[1].electric_v_m, r.forward[1].magnetic_a_m),
        join(r.backward[0].electric_v_m, r.backward[0].magnetic_a_m),
        join(r.backward[1].electric_v_m, r.backward[1].magnetic_a_m),
    ]
}

fn jets(r: UniaxialSlabThicknessDerivatives) -> [[[C64; 6]; 6]; 2] {
    let values = [
        r.reflected,
        r.transmitted,
        r.forward[0],
        r.forward[1],
        r.backward[0],
        r.backward[1],
    ];
    core::array::from_fn(|order| {
        values.map(|j| {
            let (e, h) = if order == 0 {
                (j.electric_first_v_m2, j.magnetic_first_a_m2)
            } else {
                (j.electric_second_v_m3, j.magnetic_second_a_m3)
            };
            [e[0], e[1], e[2], h[0], h[1], h[2]]
        })
    })
}

#[test]
fn homogeneous_thickness_jet_has_native_moving_exit_phase_and_constant_power() {
    let electric = [C64::new(0.7, -0.4), C64::new(-0.2, 0.6)];
    for depth in [1e-8, 0.83, 10_000.0] {
        let data = input(1.5, 1.5, 1.5, 1.5, depth, [0.6, -0.4, 0.8]);
        let actual = solve_thickness(data, Incidence::Cosine(0.8), electric).unwrap();
        assert_eq!(
            actual.response,
            solve(data, Incidence::Cosine(0.8), electric).unwrap()
        );
        let normal = actual.response.exit_normal_wavenumber_rad_m.re;
        let generator = C64::new(0.0, -normal);
        for (field, first, second) in [
            (
                actual.response.transmitted_electric_v_m,
                actual.transmitted.electric_first_v_m2,
                actual.transmitted.electric_second_v_m3,
            ),
            (
                actual.response.transmitted_magnetic_a_m,
                actual.transmitted.magnetic_first_a_m2,
                actual.transmitted.magnetic_second_a_m3,
            ),
        ] {
            for j in 0..3 {
                close(first[j], generator * field[j], 3e-14);
                close(second[j], generator * generator * field[j], 8e-14);
            }
        }
        for value in [
            actual.reflected_power_first_per_m,
            actual.reflected_power_second_per_m2,
            actual.transmitted_power_first_per_m,
            actual.transmitted_power_second_per_m2,
        ] {
            assert!(
                value.abs() < 2e-14,
                "constant homogeneous power derivative={value}"
            );
        }
    }
}

#[test]
fn coupled_and_tunneling_thickness_jets_match_independent_five_point_controls() {
    let electric = [C64::new(0.7, 0.2), C64::new(-0.3, 0.5)];
    let step = 0.001;
    let cases = [
        (input(1.0, 1.4, 1.9, 1.2, 0.7, [1.0, 2.0, 3.0]), 0.6),
        (input(2.0, 1.0, 2.0, 1.8, 0.53, [0.6, 0.4, 0.8]), -1.5),
        (input(1.5, 1.0, 1.0, 1.5, 0.83, [0.0, 0.0, 1.0]), 1.2),
        (input(1.5, 1.7, 1.9, 0.8, 0.71, [1.0, 1.0, 1.0]), 1.2),
    ];
    for (data, tangent) in cases {
        let incidence = Incidence::TangentialIndex(tangent);
        let actual = solve_thickness(data, incidence, electric).unwrap();
        let samples: [_; 5] = core::array::from_fn(|i| {
            solve(
                Input {
                    thickness_m: data.thickness_m + (i as f64 - 2.0) * step,
                    ..data
                },
                incidence,
                electric,
            )
            .unwrap()
        });
        assert_eq!(actual.response, samples[2]);
        let fields = samples.map(all_fields);
        let actual_fields = jets(actual);
        for wave in 0..6 {
            for component in 0..6 {
                let [mm, m, center, p, pp] = fields.map(|v| v[wave][component]);
                let expected = [
                    (mm - m.scale(8.0) + p.scale(8.0) - pp).scale(1.0 / (12.0 * step)),
                    (-pp + p.scale(16.0) - center.scale(30.0) + m.scale(16.0) - mm)
                        .scale(1.0 / (12.0 * step * step)),
                ];
                for order in 0..2 {
                    let tolerance = if component < 3 {
                        [3e-10, 8e-8][order]
                    } else {
                        [3e-12, 8e-10][order]
                    };
                    close(
                        actual_fields[order][wave][component],
                        expected[order],
                        tolerance,
                    );
                }
            }
        }
        for (values, first, second) in [
            (
                samples.map(|s| s.reflected_power_fraction),
                actual.reflected_power_first_per_m,
                actual.reflected_power_second_per_m2,
            ),
            (
                samples.map(|s| s.transmitted_power_fraction),
                actual.transmitted_power_first_per_m,
                actual.transmitted_power_second_per_m2,
            ),
        ] {
            let [mm, m, center, p, pp] = values;
            let fd1 = (mm - 8.0 * m + 8.0 * p - pp) / (12.0 * step);
            let fd2 = (-pp + 16.0 * p - 30.0 * center + 16.0 * m - mm) / (12.0 * step * step);
            assert!((first - fd1).abs() < 5e-10, "{first} != {fd1}");
            assert!((second - fd2).abs() < 2e-8, "{second} != {fd2}");
        }
        assert!(
            (actual.reflected_power_first_per_m + actual.transmitted_power_first_per_m).abs()
                < 2e-13
        );
        assert!(
            (actual.reflected_power_second_per_m2 + actual.transmitted_power_second_per_m2).abs()
                < 2e-12
        );
    }
}

#[test]
fn etalon_reflection_zero_retains_nonzero_field_slope_and_positive_power_curvature() {
    let data = input(1.0, 1.5, 2.0, 1.0, core::f64::consts::PI, [1.0, 0.0, 0.0]);
    let actual = solve_thickness(data, Incidence::Cosine(1.0), [C64::ONE, C64::ZERO]).unwrap();
    assert!(actual.response.reflected_electric_v_m[0].abs() < 2e-15);
    // Scalar Airy r=(r01+r12 exp(-2 i n d))/(1+r01*r12 exp(-2 i n d));
    // with n=2 and symmetric vacuum ports, r'(pi)=-1.5 i.
    close(
        actual.reflected.electric_first_v_m2[0],
        C64::new(0.0, -1.5),
        8e-15,
    );
    assert!(actual.reflected_power_first_per_m.abs() < 2e-14);
    assert!((actual.reflected_power_second_per_m2 - 4.5).abs() < 3e-13);
    assert!((actual.transmitted_power_second_per_m2 + 4.5).abs() < 3e-13);
}

#[test]
fn thickness_units_amplitude_scaling_and_opaque_tunneling_are_retained() {
    let data = input(1.5, 1.0, 1.0, 1.5, 280.0, [0.0, 0.0, 1.0]);
    let electric = [C64::ZERO, C64::new(0.7, 0.2)];
    let base = solve_thickness(data, Incidence::TangentialIndex(1.2), electric).unwrap();
    assert!(base.transmitted.electric_first_v_m2[1].abs() > 0.0);
    assert!(base.transmitted_power_first_per_m < 0.0);
    assert!(base.transmitted_power_second_per_m2 > 0.0);
    for amplitude in [1e-60, 1e60] {
        let scaled = solve_thickness(
            data,
            Incidence::TangentialIndex(1.2),
            electric.map(|z| z.scale(amplitude)),
        )
        .unwrap();
        close(
            scaled.transmitted.electric_first_v_m2[1].scale(1.0 / amplitude),
            base.transmitted.electric_first_v_m2[1],
            base.transmitted.electric_first_v_m2[1].abs() * 2e-13,
        );
        assert!(
            (scaled.transmitted_power_first_per_m / base.transmitted_power_first_per_m - 1.0).abs()
                < 2e-13
        );
    }
    let mut scaled = data;
    scaled.crystal.omega *= 1e6;
    scaled.thickness_m /= 1e6;
    let optical = solve_thickness(scaled, Incidence::TangentialIndex(1.2), electric).unwrap();
    close(
        optical.transmitted.electric_first_v_m2[1].scale(1e-6),
        base.transmitted.electric_first_v_m2[1],
        base.transmitted.electric_first_v_m2[1].abs() * 3e-13,
    );
    close(
        optical.transmitted.electric_second_v_m3[1].scale(1e-12),
        base.transmitted.electric_second_v_m3[1],
        base.transmitted.electric_second_v_m3[1].abs() * 3e-13,
    );
    for depth in [0.0, -1.0, f64::INFINITY] {
        assert!(
            solve_thickness(
                Input {
                    thickness_m: depth,
                    ..data
                },
                Incidence::TangentialIndex(1.2),
                electric
            )
            .is_err()
        );
    }
    let extinguished = Input {
        thickness_m: 2000.0,
        exit_index: 0.8,
        ..data
    };
    assert!(solve(extinguished, Incidence::TangentialIndex(1.2), electric).is_ok());
    assert!(solve_thickness(extinguished, Incidence::TangentialIndex(1.2), electric).is_err());
}
