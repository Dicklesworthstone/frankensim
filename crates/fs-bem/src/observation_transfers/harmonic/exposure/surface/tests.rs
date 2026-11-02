use super::*;
use fs_fft::harmonic_window::HarmonicWindowProfile;

type Law = HarmonicExposureFactorStatistics;

fn exposure_limits() -> HarmonicExposureLimits {
    HarmonicExposureLimits {
        maximum_lines: 32,
        maximum_factors: 8,
        maximum_transfer_products: 8192,
        maximum_receiver_frequency_pairs: 8192,
        maximum_workspace_bytes: 1_000_000,
    }
}

fn limits() -> HarmonicSurfacePowerLimits {
    HarmonicSurfacePowerLimits {
        maximum_quadrature_points: 32,
        maximum_quadrature_frequency_pairs: 8192,
        maximum_generated_terms: 128,
        maximum_workspace_bytes: 1_000_000,
    }
}

fn owner(omega: f64, fields: &[[C64; 4]]) -> PreparedHarmonicObservation {
    PreparedHarmonicObservation {
        omega,
        medium: crate::helmholtz::Medium {
            density: 1.2,
            sound_speed: 343.0,
        },
        sources: vec![SourceId::StructuralMode("original-mode".into())],
        receivers: fields
            .iter()
            .enumerate()
            .map(|(i, _)| ReceiverSpec {
                id: format!("physical-{i}"),
                position_m: [0.2 * i as f64, -0.3, 0.1],
            })
            .collect(),
        producer_hash: ContentHash([43; 32]),
        transfers: fields.to_vec(),
    }
}

struct Fixture {
    owners: [PreparedHarmonicObservation; 2],
    amplitudes: [[C64; 2]; 3],
}

impl Fixture {
    fn new() -> Self {
        Self {
            owners: [
                owner(
                    5.0,
                    &[
                        [
                            C64::new(1.1, 0.2),
                            C64::new(-0.3, 0.7),
                            C64::new(0.2, -0.4),
                            C64::new(0.5, 0.1),
                        ],
                        [
                            C64::new(-0.2, 0.9),
                            C64::new(0.6, -0.1),
                            C64::new(-0.4, 0.3),
                            C64::new(0.1, -0.8),
                        ],
                    ],
                ),
                owner(
                    9.0,
                    &[
                        [
                            C64::new(0.4, -0.3),
                            C64::new(0.2, 0.1),
                            C64::new(-0.6, -0.2),
                            C64::new(0.3, -0.4),
                        ],
                        [
                            C64::new(0.8, 0.2),
                            C64::new(-0.4, -0.5),
                            C64::new(0.2, 0.7),
                            C64::new(-0.1, 0.2),
                        ],
                    ],
                ),
            ],
            amplitudes: [
                [C64::new(0.7, 0.2), C64::new(-0.4, 0.5)],
                [C64::new(0.3, -0.1), C64::new(0.6, 0.8)],
                [C64::new(-0.2, 0.9), C64::new(0.4, -0.3)],
            ],
        }
    }

    fn exposure(&self) -> PreparedHarmonicExposure<'_> {
        let banks: [[HarmonicExposureLine<'_>; 2]; 3] = core::array::from_fn(|b| {
            core::array::from_fn(|f| HarmonicExposureLine {
                observation: &self.owners[f],
                coefficients: &self.amplitudes[b][f..f + 1],
            })
        });
        PreparedHarmonicExposure::new_with_statistics(
            &self.owners[0],
            "physical-clock",
            &banks[0],
            &[&banks[1], &banks[2]],
            &[Law::Real, Law::ProperComplex],
            exposure_limits(),
        )
        .unwrap()
    }

    // Exact ±1 real coordinate and four proper phases; all original frequencies
    // and receivers share the same coordinate in each finite realization.
    fn signal(&self, state: usize, r: usize, c: usize, time: f64) -> f64 {
        let real = if state % 2 == 0 { C64::ONE } else { -C64::ONE };
        let proper = [C64::ONE, C64::new(0.0, 1.0), -C64::ONE, C64::new(0.0, -1.0)][state / 2];
        let mut value = 0.0;
        for f in 0..2 {
            let a = self.owners[f].transfers[r][c]
                * (self.amplitudes[0][f]
                    + real * self.amplitudes[1][f]
                    + proper * self.amplitudes[2][f]);
            let (sin, cos) = (self.owners[f].omega * time).sin_cos();
            value += a.re * cos + a.im * sin;
        }
        value
    }

    fn direct_power(&self, quadrature: &[HarmonicSurfaceQuadrature], time: f64) -> f64 {
        let mut power = 0.0;
        for state in 0..8 {
            for sample in quadrature {
                let p = self.signal(state, sample.receiver, 0, time);
                for c in 0..3 {
                    power += sample.oriented_area_m2[c]
                        * p
                        * self.signal(state, sample.receiver, c + 1, time)
                        / 8.0;
                }
            }
        }
        power
    }

    fn direct_mean(&self, quadrature: &[HarmonicSurfaceQuadrature], window: HarmonicWindow) -> f64 {
        let steps = 2048;
        let mut mean = 0.0;
        for i in 0..=steps {
            let u = i as f64 / steps as f64;
            let simpson = if i == 0 || i == steps {
                1.0
            } else if i % 2 == 0 {
                2.0
            } else {
                4.0
            };
            let sensitivity = match window.profile() {
                HarmonicWindowProfile::Uniform => 1.0,
                HarmonicWindowProfile::Triangular => 2.0 * (1.0 - (2.0 * u - 1.0).abs()),
            };
            mean += simpson
                * sensitivity
                * self.direct_power(quadrature, window.start_s() + window.duration_s() * u)
                / (3.0 * steps as f64);
        }
        mean
    }
}

fn quadrature() -> Vec<HarmonicSurfaceQuadrature> {
    vec![
        HarmonicSurfaceQuadrature {
            receiver: 0,
            oriented_area_m2: [0.4, -0.3, 0.2],
        },
        HarmonicSurfaceQuadrature {
            receiver: 1,
            oriented_area_m2: [0.1, 0.7, -0.2],
        },
        HarmonicSurfaceQuadrature {
            receiver: 0,
            oriented_area_m2: [-0.15, 0.1, 0.05],
        },
    ]
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:e} != {expected:e}, tolerance {tolerance:e}"
    );
}

#[test]
fn g1_mixed_surface_spectrum_matches_direct_real_states_and_original_window_flux() {
    let fixture = Fixture::new();
    let exposure = fixture.exposure();
    let quadrature = quadrature();
    let surface = PreparedHarmonicSurfacePower::new(&exposure, &quadrature, limits()).unwrap();
    assert_eq!(surface.factor_statistics(), [Law::Real, Law::ProperComplex]);
    assert_eq!(surface.generated_terms(), 20);
    assert_eq!(surface.quadrature_frequency_pairs(), 36);
    let frequencies: Vec<_> = surface
        .spectrum()
        .tones()
        .iter()
        .map(|t| t.angular_frequency_dd().to_f64())
        .collect();
    assert_eq!(frequencies, [4.0, 10.0, 14.0, 18.0]);
    for time in [-0.2, 0.0, 0.13, 0.7, 1.2] {
        close(
            surface.evaluate(time).unwrap(),
            fixture.direct_power(&quadrature, time),
            3e-14,
        );
    }
    for profile in [
        HarmonicWindowProfile::Uniform,
        HarmonicWindowProfile::Triangular,
    ] {
        let window = HarmonicWindow::with_profile(0.13, 0.37, profile).unwrap();
        let power = surface.mean_power_w(window).unwrap();
        close(power, fixture.direct_mean(&quadrature, window), 2e-10);
        let points = exposure.observe_window(window).unwrap();
        let original: f64 = quadrature
            .iter()
            .map(|q| {
                (0..3)
                    .map(|c| {
                        q.oriented_area_m2[c] * points.points()[q.receiver].mean_intensity_w_m2[c]
                    })
                    .sum::<f64>()
            })
            .sum();
        close(power, original, 3e-14);
    }
    let uniform = HarmonicWindow::new(0.13, 0.37).unwrap();
    close(
        surface.energy_j_dd(0.13, 0.37).unwrap().to_f64(),
        fixture.direct_mean(&quadrature, uniform) * 0.37,
        1e-10,
    );
    let first = surface.energy_j_dd(0.125, 0.125).unwrap();
    let second = surface
        .energy_j_dd_from_clock(Dd::from_f64(0.125) + Dd::from_f64(0.125), 0.25)
        .unwrap();
    close(
        (first + second).to_f64(),
        surface.energy_j_dd(0.125, 0.375).unwrap().to_f64(),
        3e-15,
    );
    assert_eq!(surface.energy_j_dd(0.13, 0.0).unwrap(), Dd::ZERO);
}

#[test]
fn g3_orientation_subdivision_permutation_and_original_scope_are_preserved() {
    let fixture = Fixture::new();
    // The surface owns its quadrature and statistics; only original physical
    // identity/clock are borrowed after the exposure is dropped.
    let surface = {
        let temporary = fixture.exposure();
        PreparedHarmonicSurfacePower::new(&temporary, &quadrature(), limits()).unwrap()
    };
    assert_eq!(surface.physical_receivers(), &fixture.owners[0].receivers);
    assert_eq!(surface.sources(), &fixture.owners[0].sources);
    assert_eq!(surface.producer_hash(), fixture.owners[0].producer_hash);
    assert_eq!(surface.source_clock(), "physical-clock");
    assert_eq!(surface.quadrature(), quadrature());
    assert_eq!(
        surface.producer_uncertainty(),
        ProducerUncertainty::Unresolved {
            source_hash: fixture.owners[0].producer_hash
        }
    );
    let exposure = fixture.exposure();
    let reversed: Vec<_> = quadrature()
        .iter()
        .map(|q| HarmonicSurfaceQuadrature {
            receiver: q.receiver,
            oriented_area_m2: q.oriented_area_m2.map(|a| -a),
        })
        .collect();
    let reversed = PreparedHarmonicSurfacePower::new(&exposure, &reversed, limits()).unwrap();
    let mut divided = Vec::new();
    for q in quadrature() {
        for fraction in [0.25, 0.75] {
            divided.push(HarmonicSurfaceQuadrature {
                receiver: q.receiver,
                oriented_area_m2: q.oriented_area_m2.map(|a| a * fraction),
            });
        }
    }
    divided.reverse();
    divided.push(HarmonicSurfaceQuadrature {
        receiver: 1,
        oriented_area_m2: [0.0; 3],
    });
    let divided = PreparedHarmonicSurfacePower::new(&exposure, &divided, limits()).unwrap();
    assert_eq!(divided.quadrature().len(), 7);
    assert_eq!(divided.generated_terms(), surface.generated_terms());
    for time in [0.0, 0.11, 0.37, 0.7] {
        let expected = surface.evaluate(time).unwrap();
        close(reversed.evaluate(time).unwrap(), -expected, 0.0);
        close(divided.evaluate(time).unwrap(), expected, 3e-14);
    }
    let window = HarmonicWindow::with_profile(0.2, 0.3, HarmonicWindowProfile::Triangular).unwrap();
    close(
        reversed.mean_power_w(window).unwrap(),
        -surface.mean_power_w(window).unwrap(),
        0.0,
    );
    close(
        divided.mean_power_w(window).unwrap(),
        surface.mean_power_w(window).unwrap(),
        3e-14,
    );
    close(
        reversed.energy_j_dd(0.2, 0.3).unwrap().to_f64(),
        -surface.energy_j_dd(0.2, 0.3).unwrap().to_f64(),
        0.0,
    );
}

#[test]
fn g1_standing_wave_surface_retains_signed_reactive_exchange_and_real_statistics() {
    let reference = owner(3.0, &[[C64::ONE, C64::new(0.0, 1.0), C64::ZERO, C64::ZERO]]);
    let coefficients = [C64::ONE];
    let lines = [HarmonicExposureLine {
        observation: &reference,
        coefficients: &coefficients,
    }];
    let mean =
        PreparedHarmonicExposure::new(&reference, "physical-clock", &lines, &[], exposure_limits())
            .unwrap();
    let real = PreparedHarmonicExposure::new_with_statistics(
        &reference,
        "physical-clock",
        &[],
        &[&lines],
        &[Law::Real],
        exposure_limits(),
    )
    .unwrap();
    let proper = PreparedHarmonicExposure::new(
        &reference,
        "physical-clock",
        &[],
        &[&lines],
        exposure_limits(),
    )
    .unwrap();
    let q = [HarmonicSurfaceQuadrature {
        receiver: 0,
        oriented_area_m2: [2.0, 0.0, 0.0],
    }];
    let mean = PreparedHarmonicSurfacePower::new(&mean, &q, limits()).unwrap();
    let real = PreparedHarmonicSurfacePower::new(&real, &q, limits()).unwrap();
    let proper = PreparedHarmonicSurfacePower::new(&proper, &q, limits()).unwrap();
    assert_eq!(mean.spectrum().dc_dd(), Dd::ZERO);
    assert_eq!(mean.spectrum().tones().len(), 1);
    assert_eq!(
        mean.spectrum().tones()[0].angular_frequency_dd().to_f64(),
        6.0
    );
    close(
        mean.evaluate(core::f64::consts::PI / 12.0).unwrap(),
        1.0,
        3e-15,
    );
    close(
        mean.evaluate(core::f64::consts::PI / 4.0).unwrap(),
        -1.0,
        3e-15,
    );
    for time in [0.1, 0.4, 0.8] {
        close(
            real.evaluate(time).unwrap(),
            mean.evaluate(time).unwrap(),
            0.0,
        );
        close(proper.evaluate(time).unwrap(), 0.0, 0.0);
    }
    let h = core::f64::consts::PI / 12.0;
    close(
        mean.mean_power_w(HarmonicWindow::new(0.0, h).unwrap())
            .unwrap(),
        2.0 / core::f64::consts::PI,
        3e-15,
    );
    close(mean.energy_j_dd(0.0, h).unwrap().to_f64(), 1.0 / 6.0, 3e-15);
    close(
        mean.mean_power_w(HarmonicWindow::new(0.0, core::f64::consts::PI / 3.0).unwrap())
            .unwrap(),
        0.0,
        3e-15,
    );
}

#[test]
fn g1_spatial_contraction_keeps_weak_flux_and_only_frequency_pair_storage() {
    let fields = [[C64::ONE, C64::ONE, C64::ZERO, C64::ZERO]];
    let a = owner(4_503_599_627_370_496.0, &fields);
    let b = owner(4_503_599_627_370_497.0, &fields);
    let positive = [C64::ONE];
    let negative = [-C64::ONE];
    let factor = [
        HarmonicExposureLine {
            observation: &a,
            coefficients: &positive,
        },
        HarmonicExposureLine {
            observation: &b,
            coefficients: &negative,
        },
    ];
    let exposure =
        PreparedHarmonicExposure::new(&a, "physical-clock", &[], &[&factor], exposure_limits())
            .unwrap();
    let q: Vec<_> = [1e16, 1.0, -1e16]
        .into_iter()
        .map(|area| HarmonicSurfaceQuadrature {
            receiver: 0,
            oriented_area_m2: [area, 0.0, 0.0],
        })
        .collect();
    let mut bounded = limits();
    bounded.maximum_generated_terms = 4; // Four pairs, not twelve quadrature products.
    let surface = PreparedHarmonicSurfacePower::new(&exposure, &q, bounded).unwrap();
    assert_eq!(surface.generated_terms(), 4);
    assert_eq!(surface.quadrature_frequency_pairs(), 12);
    assert_eq!(surface.spectrum().tones().len(), 1);
    assert_eq!(
        surface.spectrum().tones()[0]
            .angular_frequency_dd()
            .to_f64(),
        1.0
    );
    let h = 1e-8;
    close(surface.evaluate(h).unwrap(), h * h / 2.0, 3e-31);
    close(
        surface.energy_j_dd(0.0, h).unwrap().to_f64(),
        h * h * h / 6.0,
        3e-40,
    );
}

#[test]
fn g0_surface_quadrature_admission_cancellation_and_darkness_are_atomic() {
    let fixture = Fixture::new();
    let exposure = fixture.exposure();
    let q = quadrature();
    let surface = PreparedHarmonicSurfacePower::new(&exposure, &q, limits()).unwrap();
    for (index, limit) in [
        (0, q.len() - 1),
        (1, surface.quadrature_frequency_pairs() - 1),
        (2, surface.generated_terms() - 1),
        (3, surface.workspace_bytes() - 1),
    ] {
        let mut bounded = limits();
        match index {
            0 => bounded.maximum_quadrature_points = limit,
            1 => bounded.maximum_quadrature_frequency_pairs = limit,
            2 => bounded.maximum_generated_terms = limit,
            _ => bounded.maximum_workspace_bytes = limit,
        }
        assert!(matches!(
            PreparedHarmonicSurfacePower::new(&exposure, &q, bounded),
            Err(TransferError::Cap { .. })
        ));
    }
    for bad_sample in [
        HarmonicSurfaceQuadrature {
            receiver: 2,
            oriented_area_m2: [1.0, 0.0, 0.0],
        },
        HarmonicSurfaceQuadrature {
            receiver: 0,
            oriented_area_m2: [f64::NAN, 0.0, 0.0],
        },
        HarmonicSurfaceQuadrature {
            receiver: 0,
            oriented_area_m2: [0.0, f64::INFINITY, 0.0],
        },
    ] {
        assert!(matches!(
            PreparedHarmonicSurfacePower::new(&exposure, &[bad_sample], limits()),
            Err(TransferError::Invalid(_))
        ));
    }
    let mut polls = 0;
    let _ = PreparedHarmonicSurfacePower::new_cancellable(&exposure, &q, limits(), || {
        polls += 1;
        false
    })
    .unwrap();
    for stop in [1, 12, polls - 1] {
        let mut count = 0;
        assert!(matches!(
            PreparedHarmonicSurfacePower::new_cancellable(&exposure, &q, limits(), || {
                count += 1;
                count >= stop
            },),
            Err(TransferError::Cancelled)
        ));
    }
    let window = HarmonicWindow::new(0.1, 0.2).unwrap();
    assert!(matches!(
        surface.evaluate_cancellable(0.1, || true),
        Err(TransferError::Cancelled)
    ));
    assert!(matches!(
        surface.mean_power_w_cancellable(window, || true),
        Err(TransferError::Cancelled)
    ));
    assert!(matches!(
        surface.energy_j_dd_cancellable(0.1, 0.2, || true),
        Err(TransferError::Cancelled)
    ));
    close(
        surface.evaluate(0.1).unwrap(),
        fixture.direct_power(&q, 0.1),
        3e-14,
    );
    let mut empty_limits = limits();
    empty_limits.maximum_quadrature_frequency_pairs = 0;
    empty_limits.maximum_generated_terms = 0;
    let dark = PreparedHarmonicSurfacePower::new(&exposure, &[], empty_limits).unwrap();
    assert_eq!(dark.generated_terms(), 0);
    assert_eq!(dark.quadrature_frequency_pairs(), 0);
    assert_eq!(dark.factor_statistics(), [Law::Real, Law::ProperComplex]);
    assert_eq!(dark.physical_receivers(), surface.physical_receivers());
    assert_eq!(dark.evaluate(0.1).unwrap(), 0.0);
    assert_eq!(dark.mean_power_w(window).unwrap(), 0.0);
    assert_eq!(dark.energy_j_dd(0.1, 0.2).unwrap(), Dd::ZERO);
    let zero = [HarmonicSurfaceQuadrature {
        receiver: 0,
        oriented_area_m2: [0.0; 3],
    }];
    let zero = PreparedHarmonicSurfacePower::new(&exposure, &zero, limits()).unwrap();
    assert_eq!(zero.quadrature().len(), 1);
    assert_eq!(zero.generated_terms(), surface.generated_terms());
    assert_eq!(zero.evaluate(0.1).unwrap(), 0.0);
    assert!(matches!(
        surface.energy_j_dd(0.1, -0.1),
        Err(TransferError::Invalid(_))
    ));
}
