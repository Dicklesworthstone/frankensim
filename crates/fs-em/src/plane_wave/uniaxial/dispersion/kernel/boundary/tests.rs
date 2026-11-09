//! Independent controls against original carrier fields and closed-form ports.
use super::*;
use crate::plane_wave::uniaxial::kernel::{self as carrier, Input as CarrierInput};

const MU: f64 = 1.256_637_062_12e-6;
const OMEGA: f64 = core::f64::consts::TAU * 299_792_458.0 / 550e-9;

#[derive(Clone, Copy)]
struct Law {
    n: f64,
    slope: f64,
    curvature: f64,
}
impl Law {
    fn constant(n: f64) -> Self {
        Self {
            n,
            slope: 0.0,
            curvature: 0.0,
        }
    }
    fn at(self, relative_frequency: f64) -> f64 {
        let x = relative_frequency - 1.0;
        self.n + self.slope * x + 0.5 * self.curvature * x * x
    }
    fn relative_epsilon(self) -> [f64; 2] {
        [
            2.0 * self.slope / self.n,
            2.0 * (self.slope * self.slope + self.n * self.curvature) / (self.n * self.n),
        ]
    }
}
fn native(
    laws: [Law; 3],
    incidence: Incidence,
    axis: [f64; 2],
    mode: UniaxialMode,
    mu_ratio: f64,
    relative_frequency: f64,
) -> UniaxialInterfaceResponse {
    carrier::solve(CarrierInput {
        omega_rad_s: OMEGA * relative_frequency,
        incident_index: laws[0].at(relative_frequency),
        ordinary_index: laws[1].at(relative_frequency),
        extraordinary_index: laws[2].at(relative_frequency),
        incident_mu_h_m: MU,
        crystal_mu_h_m: MU * mu_ratio,
        optic_axis: axis,
        incidence,
        mode,
    })
    .unwrap()
}
fn input(
    laws: [Law; 3],
    incidence: Incidence,
    axis: [f64; 2],
    mode: UniaxialMode,
    constraint: IncidenceFrequencyConstraint,
    mu_ratio: f64,
) -> Input {
    Input {
        indices: laws.map(|law| law.n),
        relative_permittivity_derivatives: laws.map(Law::relative_epsilon),
        incidence,
        axis,
        constraint,
        carrier: native(laws, incidence, axis, mode, mu_ratio, 1.0),
    }
}
fn near(actual: C64, expected: C64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:?} != {expected:?}; tolerance={tolerance:.3e}"
    );
}
fn root(value: f64) -> C64 {
    if value >= 0.0 {
        C64::from_re(value.sqrt())
    } else {
        C64::new(0.0, -(-value).sqrt())
    }
}

#[test]
fn fixed_angle_nondispersive_fields_are_stationary_on_both_outgoing_branches() {
    for (indices, incidence, axis) in [
        ([1.0, 1.5, 1.8], Incidence::Cosine(0.8), [0.6, 0.8]),
        ([2.0, 1.0, 2.0], Incidence::TangentialIndex(1.5), [1.0, 1.0]),
        (
            [2.2, 1.1, 1.6],
            Incidence::TangentialIndex(-1.9),
            [-0.8, 0.6],
        ),
        (
            [8.0, 5.0, 35.0],
            Incidence::TangentialIndex(7.0),
            [1.0, 1.0],
        ),
    ] {
        for mode in [UniaxialMode::Ordinary, UniaxialMode::Extraordinary] {
            let data = input(
                indices.map(Law::constant),
                incidence,
                axis,
                mode,
                IncidenceFrequencyConstraint::FixedAngle,
                1.6,
            );
            if indices[0] == 2.0 && mode == UniaxialMode::Extraordinary {
                assert!(data.carrier.transmitted_normal_wavenumber_rad_m.re < 0.0);
                assert_eq!(data.carrier.wave_kind, UniaxialWaveKind::Propagating);
            }
            if indices[0] == 8.0 && mode == UniaxialMode::Extraordinary {
                assert_eq!(data.carrier.transmitted_normal_wavenumber_rad_m, C64::ZERO);
                assert_eq!(data.carrier.wave_kind, UniaxialWaveKind::Propagating);
            }
            assert_eq!(amplitude_derivatives(data).unwrap(), [C64::ZERO; 2]);
        }
    }
}

#[test]
fn fixed_transverse_isotropic_fields_match_differentiated_admittances() {
    for (ni, nt, p) in [(1.5, 1.4, 0.8), (2.0, 1.0, 1.5)] {
        for mu_ratio in [1.0, 1.6] {
            let qi = root(ni * ni - p * p);
            let qt = root(nt * nt - p * p);
            let qi1 = C64::from_re(ni * ni) / qi;
            let qt1 = C64::from_re(nt * nt) / qt;
            let qi2 = C64::from_re(-ni * ni * p * p) / (qi * qi * qi);
            let qt2 = C64::from_re(-nt * nt * p * p) / (qt * qt * qt);
            for mode in [UniaxialMode::Ordinary, UniaxialMode::Extraordinary] {
                // These q are physical normal wavenumbers in fixed k0 units:
                // q(s)^2=n^2*s^2-p^2, with s=omega/omega0.
                let (a, b, a1, b1, a2, b2) = match mode {
                    UniaxialMode::Ordinary => (
                        qi.scale(mu_ratio),
                        qt,
                        qi1.scale(mu_ratio),
                        qt1,
                        qi2.scale(mu_ratio),
                        qt2,
                    ),
                    UniaxialMode::Extraordinary => (
                        qt.scale(mu_ratio * ni * ni),
                        qi.scale(nt * nt),
                        qt1.scale(mu_ratio * ni * ni),
                        qi1.scale(nt * nt),
                        qt2.scale(mu_ratio * ni * ni),
                        qi2.scale(nt * nt),
                    ),
                };
                let sum = a + b;
                let cross = a1 * b - a * b1;
                let expected_first = cross.scale(2.0) / (sum * sum);
                let expected_second = (a2 * b - a * b2).scale(2.0) / (sum * sum)
                    - (cross * (a1 + b1)).scale(4.0) / (sum * sum * sum);
                let data = input(
                    [ni, nt, nt].map(Law::constant),
                    Incidence::TangentialIndex(p),
                    [0.6, -0.8],
                    mode,
                    IncidenceFrequencyConstraint::FixedTransverseWavenumber,
                    mu_ratio,
                );
                let [first, second] = amplitude_derivatives(data).unwrap();
                near(first, expected_first, 8e-14 * expected_first.abs().max(1.0));
                near(
                    second,
                    expected_second,
                    2e-13 * expected_second.abs().max(1.0),
                );
            }
        }
    }
}

#[test]
fn distinct_anisotropic_source_jets_follow_original_frequency_resolved_fields() {
    let laws = [
        Law {
            n: 2.0,
            slope: -0.13,
            curvature: 0.06,
        },
        Law {
            n: 1.2,
            slope: 0.21,
            curvature: -0.17,
        },
        Law {
            n: 1.6,
            slope: -0.09,
            curvature: 0.11,
        },
    ];
    for cosine in [0.25_f64, 0.9] {
        let p = laws[0].n * ((1.0 - cosine) * (1.0 + cosine)).sqrt();
        for mode in [UniaxialMode::Ordinary, UniaxialMode::Extraordinary] {
            for constraint in [
                IncidenceFrequencyConstraint::FixedAngle,
                IncidenceFrequencyConstraint::FixedTransverseWavenumber,
            ] {
                let data = input(
                    laws,
                    Incidence::Cosine(cosine),
                    [0.6, 0.8],
                    mode,
                    constraint,
                    1.6,
                );
                let [first, second] = amplitude_derivatives(data).unwrap();
                // Test-only frequency differences replay actual carrier solves.
                // Fixed physical kx means the carrier coordinate kx/k0 is p/s.
                let h = 2e-4;
                let values: [C64; 5] = core::array::from_fn(|j| {
                    let relative = 1.0 + (j as f64 - 2.0) * h;
                    let incidence = match constraint {
                        IncidenceFrequencyConstraint::FixedAngle => Incidence::Cosine(cosine),
                        IncidenceFrequencyConstraint::FixedTransverseWavenumber => {
                            Incidence::TangentialIndex(p / relative)
                        }
                    };
                    native(laws, incidence, [0.6, 0.8], mode, 1.6, relative).reflected_tangential_e
                });
                let [a, b, c, d, e] = values;
                let expected_first = (a - b.scale(8.0) + d.scale(8.0) - e).scale(1.0 / (12.0 * h));
                let expected_second = (-a + b.scale(16.0) - c.scale(30.0) + d.scale(16.0) - e)
                    .scale(1.0 / (12.0 * h * h));
                near(first, expected_first, 5e-10);
                near(second, expected_second, 2e-7);
            }
        }
    }
}

#[test]
fn exact_zero_reflection_keeps_nonzero_raw_derivatives() {
    let n = 1.5;
    let law = Law {
        n,
        slope: 0.2,
        curvature: -0.1,
    };
    for mode in [UniaxialMode::Ordinary, UniaxialMode::Extraordinary] {
        for constraint in [
            IncidenceFrequencyConstraint::FixedAngle,
            IncidenceFrequencyConstraint::FixedTransverseWavenumber,
        ] {
            let data = input(
                [Law::constant(n), law, law],
                Incidence::Cosine(1.0),
                [0.0, 1.0],
                mode,
                constraint,
                1.0,
            );
            assert_eq!(data.carrier.reflected_tangential_e, C64::ZERO);
            let [first, second] = amplitude_derivatives(data).unwrap();
            near(first, C64::from_re(-law.slope / (2.0 * n)), 3e-15);
            near(
                second,
                C64::from_re(-law.curvature / (2.0 * n) + law.slope * law.slope / (2.0 * n * n)),
                3e-15,
            );
        }
    }
}

#[test]
fn exact_local_matches_keep_stationary_ports_at_authored_extreme_grazing() {
    for cosine in [1e-20, 1e-100, 1e-200] {
        for (indices, axis, mode) in [
            ([1.5, 1.5, 1.5], [0.37, -0.81], UniaxialMode::Ordinary),
            ([1.5, 1.5, 1.5], [0.37, -0.81], UniaxialMode::Extraordinary),
            ([5.0, 1.0, 7.0], [1.0, 1.0], UniaxialMode::Extraordinary),
        ] {
            let laws = indices.map(|n| Law {
                n,
                slope: n * 0.125,
                curvature: n * 0.25,
            });
            for constraint in [
                IncidenceFrequencyConstraint::FixedAngle,
                IncidenceFrequencyConstraint::FixedTransverseWavenumber,
            ] {
                let data = input(laws, Incidence::Cosine(cosine), axis, mode, constraint, 1.6);
                assert_eq!(data.carrier.wave_kind, UniaxialWaveKind::Propagating);
                assert_eq!(amplitude_derivatives(data).unwrap(), [C64::ZERO; 2]);
            }
        }
    }
}

#[test]
fn unmatched_fixed_transverse_grazing_keeps_finite_second_field_derivatives() {
    let cosine = 1e-100;
    let transmitted_index: f64 = 1.5;
    let transmitted_normal = (transmitted_index * transmitted_index - 1.0).sqrt();
    for mode in [UniaxialMode::Ordinary, UniaxialMode::Extraordinary] {
        let data = input(
            [1.0, transmitted_index, transmitted_index].map(Law::constant),
            Incidence::Cosine(cosine),
            [0.6, 0.8],
            mode,
            IncidenceFrequencyConstraint::FixedTransverseWavenumber,
            1.0,
        );
        let [first, second] = amplitude_derivatives(data).unwrap();
        // The exact admittance expansion has relative corrections O(cosine).
        // Log curvature itself is O(cosine^-4), outside f64, while these
        // physical field derivatives are O(cosine^-1) and O(cosine^-3).
        let coefficient = match mode {
            UniaxialMode::Ordinary => 2.0 / transmitted_normal,
            UniaxialMode::Extraordinary => {
                -2.0 * transmitted_index * transmitted_index / transmitted_normal
            }
        };
        let expected_first = coefficient / cosine;
        let expected_second = -coefficient / (cosine * cosine * cosine);
        near(
            first,
            C64::from_re(expected_first),
            2e-12 * expected_first.abs(),
        );
        near(
            second,
            C64::from_re(expected_second),
            3e-12 * expected_second.abs(),
        );
    }
}

#[test]
fn original_critical_kind_refuses_even_when_tilted_normal_phase_is_nonzero() {
    for (indices, p, axis, mode) in [
        ([2.0, 1.0, 1.0], 1.0, [0.6, 0.8], UniaxialMode::Ordinary),
        (
            [2.0, 1.0, 1.0],
            1.0,
            [0.6, 0.8],
            UniaxialMode::Extraordinary,
        ),
        (
            [6.0, 1.0, 7.0],
            5.0,
            [1.0, 1.0],
            UniaxialMode::Extraordinary,
        ),
    ] {
        let data = input(
            indices.map(Law::constant),
            Incidence::TangentialIndex(p),
            axis,
            mode,
            IncidenceFrequencyConstraint::FixedTransverseWavenumber,
            1.0,
        );
        assert_eq!(data.carrier.wave_kind, UniaxialWaveKind::Critical);
        if indices[0] == 6.0 {
            assert_ne!(data.carrier.transmitted_normal_wavenumber_rad_m, C64::ZERO);
        }
        assert!(amplitude_derivatives(data).is_err());
    }
}
