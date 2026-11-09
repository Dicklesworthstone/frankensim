use super::super::kernel as carrier;
use crate::plane_wave::dispersion::IncidenceFrequencyConstraint;
use crate::plane_wave::uniaxial::dispersion::kernel as normal;
use crate::plane_wave::uniaxial::kernel::{Incidence, UniaxialMode, UniaxialWaveKind};
use fs_math::c64::C64;

const C0: f64 = 299_792_458.0;
const OMEGA: f64 = core::f64::consts::TAU * C0 / 550e-9;
const MU: f64 = 1.256_637_062_12e-6;

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
    fn at(self, relative: f64) -> f64 {
        let x = relative - 1.0;
        self.n + self.slope * x + 0.5 * self.curvature * x * x
    }
    fn epsilon(self) -> [f64; 2] {
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
    axis: [f64; 3],
    incidence: Incidence,
    relative: f64,
) -> carrier::OrientedCrystalResponse {
    let [ni, no, ne] = laws.map(|law| law.at(relative));
    let cosine = match incidence {
        Incidence::Cosine(c) => c,
        Incidence::TangentialIndex(p) => ((ni - p) * (ni + p)).sqrt() / ni,
    };
    carrier::solve(carrier::Input {
        omega: OMEGA * relative,
        ni,
        no,
        ne,
        mui: MU,
        mut_: MU,
        axis,
        incidence,
        electric: [C64::new(0.7 * cosine, 0.1 * cosine), C64::new(0.3, -0.4)],
    })
    .unwrap()
}

fn derivatives(
    laws: [Law; 3],
    axis: [f64; 3],
    incidence: Incidence,
    mode: usize,
    constraint: IncidenceFrequencyConstraint,
) -> [C64; 2] {
    let original = native(laws, axis, incidence, 1.0).transmitted[mode];
    normal::normal_derivatives(normal::Input {
        omega: OMEGA,
        indices: laws.map(|law| law.n),
        relative_permittivity_derivatives: laws.map(Law::epsilon),
        principal_wavenumber_derivatives: [laws[1].bulk(), laws[2].bulk()],
        axis,
        incidence,
        mode: if mode == 0 {
            UniaxialMode::Ordinary
        } else {
            UniaxialMode::Extraordinary
        },
        kind: original.wave_kind,
        normal_wavenumber: original.normal_wavenumber_rad_m,
        constraint,
    })
    .unwrap()
}

fn near(actual: C64, expected: C64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:?} != {expected:?}; tolerance={tolerance:.3e}"
    );
}

#[test]
fn arbitrary_axis_jets_match_inverse_metric_and_preserve_matched_grazing() {
    let laws = [2.4, 1.1, 1.8].map(Law::constant);
    let axis = [1.0, 2.0, 3.0];
    let [x, _, z] = carrier::unit(axis).unwrap();
    for p in [-1.0, 1.0, 2.2] {
        for mode in 0..2 {
            let (a, h) = if mode == 0 {
                let a = 1.0 / (laws[1].n * laws[1].n);
                (a, a * a * p * p)
            } else {
                // Independent inverse-index metric M=(I-aa^T)/ne^2+aa^T/no^2.
                let ordinary = 1.0 / (laws[1].n * laws[1].n);
                let extraordinary = 1.0 / (laws[2].n * laws[2].n);
                let a = extraordinary + z * z * (ordinary - extraordinary);
                let b = x * z * (ordinary - extraordinary);
                let d = extraordinary + x * x * (ordinary - extraordinary);
                (a, (a * d - b * b) * p * p)
            };
            let root = if a > h {
                C64::from_re((a - h).sqrt())
            } else {
                C64::new(0.0, -(h - a).sqrt())
            };
            let result = derivatives(
                laws,
                axis,
                Incidence::TangentialIndex(p),
                mode,
                IncidenceFrequencyConstraint::FixedTransverseWavenumber,
            );
            near(result[0].scale(C0), C64::ONE / root, 2e-13);
            near(
                result[1].scale(C0 * OMEGA),
                C64::from_re(-h) / (root * root * root),
                2e-12,
            );
            assert_eq!(
                result,
                derivatives(
                    laws,
                    axis.map(|v| -v),
                    Incidence::TangentialIndex(p),
                    mode,
                    IncidenceFrequencyConstraint::FixedTransverseWavenumber
                )
            );
        }
    }

    // Exact native out-of-plane cutoff: a=[1,1,0], ni=7,no=35,ne=5.
    // Its extraordinary metric has Mzz=1/25 and Mxx=1/49.
    let laws = [7.0, 35.0, 5.0].map(Law::constant);
    for cosine in [1e-10, 1e-16, 1e-20] {
        let incidence = Incidence::Cosine(cosine);
        let original = native(laws, [1.0, 1.0, 0.0], incidence, 1.0).transmitted[1];
        assert_eq!(original.wave_kind, UniaxialWaveKind::Propagating);
        let fixed = derivatives(
            laws,
            [1.0, 1.0, 0.0],
            incidence,
            1,
            IncidenceFrequencyConstraint::FixedTransverseWavenumber,
        );
        near(
            fixed[0].scale(C0),
            C64::from_re(5.0 / cosine),
            2e-13 * 5.0 / cosine,
        );
        let curvature = -5.0 * (1.0 - cosine) * (1.0 + cosine) / cosine.powi(3);
        near(
            fixed[1].scale(C0 * OMEGA),
            C64::from_re(curvature),
            5e-13 * curvature.abs(),
        );
        let angled = derivatives(
            laws,
            [1.0, 1.0, 0.0],
            incidence,
            1,
            IncidenceFrequencyConstraint::FixedAngle,
        );
        near(
            angled[0].scale(C0),
            C64::from_re(5.0 * cosine),
            2e-13 * 5.0 * cosine,
        );
        assert_eq!(angled[1], C64::ZERO);
    }

    // The native exactly isotropic carrier skips anisotropic coefficient
    // contractions. Its original tiny cosine must survive that same limit.
    let n = 1.5;
    let laws = [n; 3].map(Law::constant);
    let cosine = 1e-20;
    let fixed = derivatives(
        laws,
        [1.0, 2.0, 3.0],
        Incidence::Cosine(cosine),
        1,
        IncidenceFrequencyConstraint::FixedTransverseWavenumber,
    );
    near(
        fixed[0].scale(C0),
        C64::from_re(n / cosine),
        2e-13 * n / cosine,
    );
    let curvature = -n * (1.0 - cosine) * (1.0 + cosine) / cosine.powi(3);
    near(
        fixed[1].scale(C0 * OMEGA),
        C64::from_re(curvature),
        5e-13 * curvature.abs(),
    );
}

#[test]
fn arbitrary_axis_distinct_material_jets_follow_both_native_frequency_paths() {
    let laws = [
        Law {
            n: 1.6,
            slope: -0.2,
            curvature: 0.08,
        },
        Law {
            n: 1.8,
            slope: 0.1,
            curvature: -0.13,
        },
        Law {
            n: 1.3,
            slope: -0.12,
            curvature: 0.17,
        },
    ];
    let axis = [0.3, 0.4, 0.5];
    let cosine: f64 = 0.8;
    let p = laws[0].n * (1.0 - cosine * cosine).sqrt();
    for mode in 0..2 {
        for constraint in [
            IncidenceFrequencyConstraint::FixedAngle,
            IncidenceFrequencyConstraint::FixedTransverseWavenumber,
        ] {
            let result = derivatives(laws, axis, Incidence::Cosine(cosine), mode, constraint);
            let h = 2e-4;
            // Test-only independent differentiation of actual vector Maxwell
            // carriers. The source laws and physical fixed-kx path are authored.
            let values: [C64; 5] = core::array::from_fn(|j| {
                let relative = 1.0 + (j as f64 - 2.0) * h;
                let incidence = match constraint {
                    IncidenceFrequencyConstraint::FixedAngle => Incidence::Cosine(cosine),
                    IncidenceFrequencyConstraint::FixedTransverseWavenumber => {
                        Incidence::TangentialIndex(p / relative)
                    }
                };
                native(laws, axis, incidence, relative).transmitted[mode]
                    .normal_wavenumber_rad_m
                    .scale(C0 / OMEGA)
            });
            let [a, b, c, d, e] = values;
            let expected_first = (a - b.scale(8.0) + d.scale(8.0) - e).scale(1.0 / (12.0 * h));
            let expected_second = (-a + b.scale(16.0) - c.scale(30.0) + d.scale(16.0) - e)
                .scale(1.0 / (12.0 * h * h));
            near(result[0].scale(C0), expected_first, 3e-11);
            near(result[1].scale(C0 * OMEGA), expected_second, 2e-7);
        }
    }
}
