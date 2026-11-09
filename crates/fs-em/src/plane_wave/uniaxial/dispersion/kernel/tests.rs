use super::*;
use crate::plane_wave::uniaxial::kernel::{
    self as carrier, Input as CarrierInput, UniaxialInterfaceResponse,
};

const MU: f64 = 1.256_637_062_12e-6;
const OMEGA: f64 = core::f64::consts::TAU * C0 / 550e-9;

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

    // n(omega) = n0 + slope*(omega/omega0-1)
    //            + curvature*(omega/omega0-1)^2/2.
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

    fn bulk(self) -> [f64; 2] {
        [
            (self.n + self.slope) / C0,
            (2.0 * self.slope + self.curvature) / C0,
        ]
    }
}

fn native(
    laws: [Law; 3],
    incidence: Incidence,
    axis: [f64; 2],
    mode: UniaxialMode,
    relative_frequency: f64,
) -> UniaxialInterfaceResponse {
    carrier::solve(CarrierInput {
        omega_rad_s: OMEGA * relative_frequency,
        incident_index: laws[0].at(relative_frequency),
        ordinary_index: laws[1].at(relative_frequency),
        extraordinary_index: laws[2].at(relative_frequency),
        incident_mu_h_m: MU,
        crystal_mu_h_m: MU,
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
) -> Input {
    let response = native(laws, incidence, axis, mode, 1.0);
    let axis = carrier::unit_axis(axis).unwrap();
    Input {
        omega: OMEGA,
        indices: laws.map(|law| law.n),
        relative_permittivity_derivatives: laws.map(Law::relative_epsilon),
        principal_wavenumber_derivatives: [laws[1].bulk(), laws[2].bulk()],
        incidence,
        axis: [axis[0], 0.0, axis[1]],
        mode,
        kind: response.wave_kind,
        normal_wavenumber: response.transmitted_normal_wavenumber_rad_m,
        constraint,
    }
}

fn near(actual: C64, expected: C64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:?} != {expected:?}; tolerance={tolerance:.3e}"
    );
}

#[test]
fn nondispersive_fixed_angle_keeps_signed_and_evanescent_carriers_linear() {
    for (indices, incidence, axis) in [
        ([1.0, 1.5, 1.8], Incidence::Cosine(0.8), [0.6, 0.8]),
        ([2.0, 1.0, 2.0], Incidence::TangentialIndex(1.5), [1.0, 1.0]),
        (
            [2.2, 1.1, 1.6],
            Incidence::TangentialIndex(1.9),
            [-0.8, 0.6],
        ),
    ] {
        for mode in [UniaxialMode::Ordinary, UniaxialMode::Extraordinary] {
            let data = input(
                indices.map(Law::constant),
                incidence,
                axis,
                mode,
                IncidenceFrequencyConstraint::FixedAngle,
            );
            let expected = data.normal_wavenumber.scale(1.0 / OMEGA);
            let [first, second] = normal_derivatives(data).unwrap();
            near(first, expected, 8.0 * f64::EPSILON * expected.abs());
            assert_eq!(second, C64::ZERO);
        }
    }
}

#[test]
fn fixed_transverse_isotropic_jets_match_closed_form_on_both_root_branches() {
    for (ni, n, p) in [(1.5, 1.4, 0.8), (2.0, 1.0, 1.5)] {
        for mode in [UniaxialMode::Ordinary, UniaxialMode::Extraordinary] {
            let data = input(
                [ni, n, n].map(Law::constant),
                Incidence::TangentialIndex(p),
                [0.6, 0.8],
                mode,
                IncidenceFrequencyConstraint::FixedTransverseWavenumber,
            );
            let q = data.normal_wavenumber.scale(C0 / OMEGA);
            let expected_first = C64::from_re(n * n) / q;
            let expected_second = C64::from_re(-n * n * p * p) / (q * q * q);
            let [first, second] = normal_derivatives(data).unwrap();
            near(
                first.scale(C0),
                expected_first,
                3e-14 * expected_first.abs(),
            );
            near(
                second.scale(C0 * OMEGA),
                expected_second,
                5e-14 * expected_second.abs(),
            );
            if q.im < 0.0 {
                assert!(first.im > 0.0 && second.im > 0.0);
            } else {
                assert!(first.re > 0.0 && second.re < 0.0);
            }
        }
    }
}

#[test]
fn tilted_fixed_transverse_uses_outgoing_energy_even_for_negative_phase() {
    let no: f64 = 1.0;
    let ne: f64 = 2.0;
    for p in [-1.5, 1.5, 1.8] {
        let data = input(
            [2.0, no, ne].map(Law::constant),
            Incidence::TangentialIndex(p),
            [1.0, 1.0],
            UniaxialMode::Extraordinary,
            IncidenceFrequencyConstraint::FixedTransverseWavenumber,
        );
        if p == 1.5 {
            assert!(data.normal_wavenumber.re < 0.0);
            assert_eq!(data.kind, UniaxialWaveKind::Propagating);
        }
        // Independent inverse-index metric M, with A=Mzz and
        // H=(Mxx*Mzz-Mxz^2)*p^2. The shifted outgoing root obeys
        // r=sqrt(A-H), C0*kz'=1/r, C0*omega*kz''=-H/r^3.
        let a = 0.5 / (no * no) + 0.5 / (ne * ne);
        let cross = 0.5 / (no * no) - 0.5 / (ne * ne);
        let h = (a * a - cross * cross) * p * p;
        let r = if a >= h {
            C64::from_re((a - h).sqrt())
        } else {
            C64::new(0.0, -(h - a).sqrt())
        };
        let [first, second] = normal_derivatives(data).unwrap();
        near(first.scale(C0), C64::ONE / r, 3e-13);
        near(
            second.scale(C0 * OMEGA),
            C64::from_re(-h) / (r * r * r),
            5e-12,
        );
    }
}

#[test]
fn distinct_original_index_jets_follow_native_frequency_laws() {
    let laws = [
        Law {
            n: 1.7,
            slope: -0.13,
            curvature: 0.06,
        },
        Law {
            n: 1.4,
            slope: 0.21,
            curvature: -0.17,
        },
        Law {
            n: 1.9,
            slope: -0.09,
            curvature: 0.11,
        },
    ];
    let cosine: f64 = 0.72;
    let p = laws[0].n * (1.0 - cosine * cosine).sqrt();
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
            );
            let [first, second] = normal_derivatives(data).unwrap();
            // Five-point differences are a test-only independent oracle against
            // original carrier solves, with physical kx fixed only when authored.
            let h = 2e-4;
            let values: [C64; 5] = core::array::from_fn(|j| {
                let relative = 1.0 + (j as f64 - 2.0) * h;
                let incidence = match constraint {
                    IncidenceFrequencyConstraint::FixedAngle => Incidence::Cosine(cosine),
                    IncidenceFrequencyConstraint::FixedTransverseWavenumber => {
                        Incidence::TangentialIndex(p / relative)
                    }
                };
                native(laws, incidence, [0.6, 0.8], mode, relative)
                    .transmitted_normal_wavenumber_rad_m
                    .scale(C0 / OMEGA)
            });
            let [a, b, c, d, e] = values;
            let expected_first = (a - b.scale(8.0) + d.scale(8.0) - e).scale(1.0 / (12.0 * h));
            let expected_second = (-a + b.scale(16.0) - c.scale(30.0) + d.scale(16.0) - e)
                .scale(1.0 / (12.0 * h * h));
            near(first.scale(C0), expected_first, 3e-11);
            near(second.scale(C0 * OMEGA), expected_second, 2e-7);
        }
    }
}

#[test]
fn matched_axial_grazing_keeps_bulk_curvature_and_does_not_match_distinct_slopes() {
    let incident = Law {
        n: 2.0,
        slope: -0.25,
        curvature: 0.4,
    };
    let ordinary = Law {
        n: 1.3,
        slope: 0.15,
        curvature: -0.3,
    };
    for cosine in [0.4, 1e-20, 1e-100] {
        let data = input(
            [incident, ordinary, incident],
            Incidence::Cosine(cosine),
            [0.0, 1.0],
            UniaxialMode::Extraordinary,
            IncidenceFrequencyConstraint::FixedAngle,
        );
        let [first, second] = normal_derivatives(data).unwrap();
        near(
            first,
            C64::from_re(ordinary.bulk()[0] * cosine),
            2e-15 * first.abs(),
        );
        // This authored local law has exactly zero K'', despite nonzero n'.
        assert_eq!(second, C64::ZERO);
    }
    let equal_phase_different_slope = Law {
        n: 2.0,
        slope: 0.3,
        curvature: -0.07,
    };
    let data = input(
        [incident, ordinary, equal_phase_different_slope],
        Incidence::Cosine(0.4),
        [0.0, 1.0],
        UniaxialMode::Extraordinary,
        IncidenceFrequencyConstraint::FixedAngle,
    );
    let [first, second] = normal_derivatives(data).unwrap();
    assert!((first.re - ordinary.bulk()[0] * 0.4).abs() > 1e-10);
    assert!(second.re.abs() > 1e-26);
}

#[test]
fn tilted_matched_cutoff_keeps_the_authored_cosine_in_fixed_transverse_jets() {
    // f=(no^2+ne^2)/2=25 exactly matches ni^2. At grazing, the
    // forward normal energy survives in the authored cosine even though the
    // normal phase approaches the nonzero shifted critical root.
    let indices = [5.0, 1.0, 7.0];
    for cosine in [1e-20, 1e-12, 1e-4] {
        let data = input(
            indices.map(Law::constant),
            Incidence::Cosine(cosine),
            [1.0, 1.0],
            UniaxialMode::Extraordinary,
            IncidenceFrequencyConstraint::FixedTransverseWavenumber,
        );
        assert_eq!(data.kind, UniaxialWaveKind::Propagating);
        let [first, second] = normal_derivatives(data).unwrap();
        // The inverse metric has A=25/49, det(M)=1/49. Therefore
        // r=sqrt(A)*cosine=(5/7)*cosine and H=A*sin(theta)^2.
        let r = (5.0 / 7.0) * cosine;
        let h = (25.0 / 49.0) * (1.0 - cosine) * (1.0 + cosine);
        near(first.scale(C0), C64::from_re(1.0 / r), 2e-12 / r);
        near(
            second.scale(C0 * OMEGA),
            C64::from_re(-h / (r * r * r)),
            5e-12 * h / (r * r * r),
        );
    }
}

#[test]
fn exact_native_critical_kind_refuses_a_frequency_jet() {
    for mode in [UniaxialMode::Ordinary, UniaxialMode::Extraordinary] {
        let data = input(
            [2.0, 1.0, 1.0].map(Law::constant),
            Incidence::TangentialIndex(1.0),
            [0.6, 0.8],
            mode,
            IncidenceFrequencyConstraint::FixedTransverseWavenumber,
        );
        assert_eq!(data.kind, UniaxialWaveKind::Critical);
        assert!(normal_derivatives(data).is_err());
    }
    let tilted = input(
        [6.0, 1.0, 7.0].map(Law::constant),
        Incidence::TangentialIndex(5.0),
        [1.0, 1.0],
        UniaxialMode::Extraordinary,
        IncidenceFrequencyConstraint::FixedTransverseWavenumber,
    );
    assert_eq!(tilted.kind, UniaxialWaveKind::Critical);
    assert!(tilted.normal_wavenumber.re < 0.0);
    assert!(normal_derivatives(tilted).is_err());
}

#[test]
fn zero_normal_phase_is_not_a_critical_energy_mode() {
    // h*p^2=no^2*ne^2 exactly: p=7, no=5, ne=35, a=[1,1].
    // The native stable root is zero, but its energy crosses the interface.
    let data = input(
        [8.0, 5.0, 35.0].map(Law::constant),
        Incidence::TangentialIndex(7.0),
        [1.0, 1.0],
        UniaxialMode::Extraordinary,
        IncidenceFrequencyConstraint::FixedTransverseWavenumber,
    );
    assert_eq!(data.normal_wavenumber, C64::ZERO);
    assert_eq!(data.kind, UniaxialWaveKind::Propagating);
    let [first, second] = normal_derivatives(data).unwrap();
    near(first.scale(C0), C64::from_re(175.0 / 24.0), 3e-14);
    near(
        second.scale(C0 * OMEGA),
        C64::from_re(-8575.0 / 13824.0),
        3e-13,
    );
}

#[test]
fn original_bulk_cancellations_and_negative_slopes_survive_fixed_transverse_projection() {
    for law in [
        Law {
            n: 1.4,
            slope: -1.4,
            curvature: 2.8,
        },
        Law {
            n: 1.5,
            slope: -2.0,
            curvature: 4.0,
        },
    ] {
        let p = 0.8;
        let data = input(
            [Law::constant(1.7), law, Law::constant(1.9)],
            Incidence::TangentialIndex(p),
            [0.6, 0.8],
            UniaxialMode::Ordinary,
            IncidenceFrequencyConstraint::FixedTransverseWavenumber,
        );
        let q = data.normal_wavenumber.re * C0 / OMEGA;
        let group = law.n + law.slope;
        let [first, second] = normal_derivatives(data).unwrap();
        if group == 0.0 {
            // This local n proportional to 1/omega has constant original K.
            // Projecting a fixed kx cannot invent group delay or GVD.
            assert_eq!(first, C64::ZERO);
            assert_eq!(second, C64::ZERO);
        } else {
            assert!(first.re < 0.0);
            near(first.scale(C0), C64::from_re(law.n * group / q), 3e-14);
            near(
                second.scale(C0 * OMEGA),
                C64::from_re(-group * group * p * p / (q * q * q)),
                3e-14,
            );
        }
    }
}
