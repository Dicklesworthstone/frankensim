//! G1/G3/G4 source-bound finite-window Maxwell stress and traction.
use super::*;
use fs_em::plane_wave::angular_spectrum::{
    AngularSpectrumExposureStressLimits, AngularSpectrumIlluminationLimits,
    AngularSpectrumIlluminationSide as Face, PreparedAngularSpectrumIllumination,
};

fn stress_limits() -> AngularSpectrumExposureStressLimits {
    AngularSpectrumExposureStressLimits {
        maximum_workspace_bytes: 8 << 20,
    }
}

fn illumination<'b, 's, 'm>(
    beam: &'b PlanarAngularSpectrum<'s, 'm>,
    right: &[[C64; 2]],
) -> PreparedAngularSpectrumIllumination<'b, 's, 'm> {
    beam.prepare_illumination(
        right,
        AngularSpectrumIlluminationLimits {
            maximum_layer_mode_pairs: 4096,
            maximum_workspace_bytes: 8 << 20,
            maximum_fft_work: 1 << 20,
        },
        || false,
    )
    .unwrap()
}

fn phasor(angle: f64) -> C64 {
    let (s, c) = angle.sin_cos();
    C64::new(c, s)
}

fn at(field: &fs_em::plane_wave::angular_spectrum::AngularSpectrumField, i: usize) -> [C64; 6] {
    let e = field.electric_v_per_m()[i];
    let h = field.magnetic_a_per_m()[i];
    [e[0], e[1], e[2], h[0], h[1], h[2]]
}

// Independent integration of real instantaneous stress. A common enumerated
// proper coefficient multiplies every frequency of its authored factor.
fn reference(
    lines: &[([C64; 6], [C64; 6], f64)],
    window: HarmonicWindow,
    epsilon: f64,
    mu: f64,
    phases: &[C64],
) -> [[f64; 3]; 3] {
    let intervals = 2048;
    let mut stress = [[0.0; 3]; 3];
    for &phase in phases {
        for step in 0..=intervals {
            let t = window.start_s() + window.duration_s() * step as f64 / intervals as f64;
            let mut real = [0.0; 6];
            for (mean, factor, omega) in lines {
                let carrier = phasor(omega * t);
                for c in 0..6 {
                    real[c] += ((mean[c] + phase * factor[c]) * carrier).re;
                }
            }
            let weight = if step == 0 || step == intervals {
                1.0
            } else if step % 2 == 0 {
                2.0
            } else {
                4.0
            };
            let trace = (0..3)
                .map(|c| epsilon * real[c] * real[c] + mu * real[c + 3] * real[c + 3])
                .sum::<f64>();
            for i in 0..3 {
                for j in 0..3 {
                    stress[i][j] += weight
                        * (epsilon * real[i] * real[j] + mu * real[i + 3] * real[j + 3]
                            - if i == j { 0.5 * trace } else { 0.0 });
                }
            }
        }
    }
    stress.map(|row| row.map(|v| v / (3.0 * intervals as f64 * phases.len() as f64)))
}

fn close_tensor(actual: [[f64; 3]; 3], expected: [[f64; 3]; 3], tolerance: f64) {
    for i in 0..3 {
        for j in 0..3 {
            assert!(
                (actual[i][j] - expected[i][j]).abs() < tolerance,
                "stress[{i},{j}] {} != {}",
                actual[i][j],
                expected[i][j]
            );
            assert_eq!(actual[i][j], actual[j][i]);
        }
    }
}

#[test]
fn g1_two_face_vector_stress_retains_shear_and_standing_wave_interference() {
    let medium = material(1.0, 1.0);
    let film = material(2.4, 1.1);
    let frequencies = [2.0, 2.3];
    let stacks = frequencies.map(|omega| {
        PlanarDielectricLayers::try_new(
            &medium,
            DielectricResponse::Lossless(&medium),
            omega,
            &[DielectricLayer {
                response: DielectricResponse::Lossless(&film),
                thickness_m: 0.37,
            }],
        )
        .unwrap()
    });
    let first = (0..4)
        .map(|i| {
            let p = phasor(-(i as f64) * core::f64::consts::FRAC_PI_2);
            [C64::new(0.8, 0.2) * p, C64::new(-0.3, 0.6)]
        })
        .collect::<Vec<_>>();
    let second = (0..4)
        .map(|i| {
            let p = phasor((i as f64) * core::f64::consts::FRAC_PI_2);
            [C64::new(-0.2, 0.4), C64::new(0.5, -0.1) * p]
        })
        .collect::<Vec<_>>();
    let left = [prepare(&stacks[0], &first), prepare(&stacks[1], &second)];
    let lines = [
        illumination(&left[0], &[[C64::new(0.4, -0.2), C64::new(0.1, 0.3)]; 4]),
        illumination(&left[1], &[[C64::new(-0.3, 0.1), C64::new(0.2, 0.5)]; 4]),
    ];
    let bank = AngularSpectrumExposure::try_new_illumination(&left[0], &lines, limits(), || false)
        .unwrap();
    let window = HarmonicWindow::new(0.17, 0.31).unwrap();
    let distance = 0.23;
    let observed = bank
        .illumination_stress_field(Face::TotalLeft, distance, window, stress_limits(), || false)
        .unwrap();
    let fields = lines
        .iter()
        .map(|line| line.field(Face::TotalLeft, distance, || false).unwrap())
        .collect::<Vec<_>>();
    let mut integral = [[0.0; 3]; 3];
    for (i, point) in observed.points().iter().enumerate() {
        let phasors = fields
            .iter()
            .zip(frequencies)
            .map(|(f, w)| (at(f, i), [C64::ZERO; 6], w))
            .collect::<Vec<_>>();
        let expected = reference(&phasors, window, 1.0, 1.0, &[C64::ONE]);
        close_tensor(point.maxwell_stress_pa, expected, 3e-12);
        for a in 0..3 {
            for b in 0..3 {
                integral[a][b] += expected[a][b] * 24.0;
            }
        }
        let normal = [0.6, 0.0, 0.8];
        let traction = point.traction_pa(normal).unwrap();
        for a in 0..3 {
            near(
                traction[a],
                expected[a][0] * 0.6 + expected[a][2] * 0.8,
                4e-12,
            );
        }
    }
    close_tensor(observed.integrated_stress_n(), integral, 3e-10);
    assert!(
        observed
            .points()
            .iter()
            .any(|p| p.maxwell_stress_pa[0][2].abs() > 1e-3)
    );
    let incoming = bank
        .illumination_stress_field(
            Face::IncomingLeft,
            distance,
            window,
            stress_limits(),
            || false,
        )
        .unwrap();
    let outgoing = bank
        .illumination_stress_field(
            Face::OutgoingLeft,
            distance,
            window,
            stress_limits(),
            || false,
        )
        .unwrap();
    assert!(
        observed
            .points()
            .iter()
            .zip(incoming.points())
            .zip(outgoing.points())
            .any(|((total, l), r)| {
                (total.maxwell_stress_pa[0][0]
                    - l.maxwell_stress_pa[0][0]
                    - r.maxwell_stress_pa[0][0])
                    .abs()
                    > 1e-3
            })
    );
}

#[test]
fn g1_stress_of_mean_and_joint_face_covariance_matches_authored_ensemble() {
    let medium = material(1.0, 1.0);
    let film = material(2.4, 1.1);
    let frequencies = [2.0, 2.3];
    let stacks = frequencies.map(|omega| {
        PlanarDielectricLayers::try_new(
            &medium,
            DielectricResponse::Lossless(&medium),
            omega,
            &[DielectricLayer {
                response: DielectricResponse::Lossless(&film),
                thickness_m: 0.37,
            }],
        )
        .unwrap()
    });
    let mean_lines = [
        constant(&stacks[0], [C64::new(0.2, -0.1), C64::new(0.1, 0.2)]),
        constant(&stacks[1], [C64::new(-0.1, 0.3), C64::new(0.2, 0.1)]),
    ];
    let left = [
        constant(&stacks[0], [C64::new(0.8, 0.2), C64::new(0.3, -0.1)]),
        constant(&stacks[1], [C64::new(-0.3, 0.6), C64::new(0.1, 0.4)]),
    ];
    let lines = [
        illumination(&left[0], &[[C64::new(-0.2, 0.4), C64::new(0.4, 0.2)]; 4]),
        illumination(&left[1], &[[C64::new(0.7, -0.1), C64::new(0.3, -0.5)]; 4]),
    ];
    let mean =
        AngularSpectrumExposure::try_new(&mean_lines[0], &mean_lines, limits(), || false).unwrap();
    let factors =
        [
            AngularSpectrumExposure::try_new_illumination(&left[0], &lines, limits(), || false)
                .unwrap(),
        ];
    let ensemble = AngularSpectrumExposureEnsemble::try_new(
        &left[0],
        Some(&mean),
        &factors,
        ensemble_limits(),
        || false,
    )
    .unwrap();
    let window = HarmonicWindow::new(0.17, 0.12).unwrap();
    let distance = 0.23;
    let observed = ensemble
        .illumination_stress_field(Face::TotalRight, distance, window, stress_limits(), || {
            false
        })
        .unwrap();
    let mean_fields = mean_lines
        .iter()
        .map(|line| line.field(Side::Transmitted, distance, || false).unwrap())
        .collect::<Vec<_>>();
    let factor_fields = lines
        .iter()
        .map(|line| line.field(Face::TotalRight, distance, || false).unwrap())
        .collect::<Vec<_>>();
    let phases = [C64::ONE, -C64::ONE, C64::new(0.0, 1.0), C64::new(0.0, -1.0)];
    for (i, point) in observed.points().iter().enumerate() {
        let phasors = (0..2)
            .map(|f| {
                (
                    at(&mean_fields[f], i),
                    at(&factor_fields[f], i),
                    frequencies[f],
                )
            })
            .collect::<Vec<_>>();
        close_tensor(
            point.maxwell_stress_pa,
            reference(&phasors, window, 1.0, 1.0, &phases),
            3e-12,
        );
    }
    let left_stress = ensemble
        .illumination_stress_field(Face::TotalLeft, 0.1, window, stress_limits(), || false)
        .unwrap();
    let resultant = ensemble
        .boundary_stress_resultant_n([0.1, distance], window, stress_limits(), || false)
        .unwrap();
    let l = left_stress.integrated_traction_n([0.0, 0.0, -1.0]).unwrap();
    let r = observed.integrated_traction_n([0.0, 0.0, 1.0]).unwrap();
    for i in 0..3 {
        near(resultant[i], l[i] + r[i], 2e-12);
    }
}

#[test]
fn g1_empty_region_short_window_stress_is_field_momentum_storage_not_body_force() {
    // Vacuum-value fixture, prescribed traveling wave, no films and no body.
    // Surface stress over two separated planes is nonzero until its stored
    // field-momentum change is included. A cycle-average hides this distinction.
    let epsilon = 8.8541878128e-12;
    let mu = 1.25663706212e-6;
    let medium = material(epsilon, mu);
    let omega = 1.0e9;
    let stack =
        PlanarDielectricLayers::try_new(&medium, DielectricResponse::Lossless(&medium), omega, &[])
            .unwrap();
    let lines = [constant(&stack, [C64::ONE, C64::from_re(0.3)])];
    let bank = AngularSpectrumExposure::try_new(&lines[0], &lines, limits(), || false).unwrap();
    let window = HarmonicWindow::new(0.17e-9, 0.11e-9).unwrap();
    let distances = [0.23, 0.31];
    let result = bank
        .boundary_stress_resultant_n(distances, window, stress_limits(), || false)
        .unwrap();
    let speed = 1.0 / (epsilon * mu).sqrt();
    let k = omega / speed;
    let area = 4.0 * 24.0;
    let momentum = |t: f64| {
        let z_left = -distances[0];
        let z_right = distances[1];
        let integral = (z_right - z_left) * 0.5
            + ((2.0 * (omega * t - k * z_left)).sin() - (2.0 * (omega * t - k * z_right)).sin())
                / (4.0 * k);
        area * epsilon / speed * 1.09 * integral
    };
    let rate = (momentum(window.start_s() + window.duration_s()) - momentum(window.start_s()))
        / window.duration_s();
    assert!(result[2].abs() > 1e-13);
    near(result[2], rate, 3e-10 * rate.abs());
    near(result[0], 0.0, 1e-20);
    near(result[1], 0.0, 1e-20);
    let cycle = bank
        .boundary_stress_resultant_n(
            distances,
            HarmonicWindow::new(0.0, core::f64::consts::TAU / omega).unwrap(),
            stress_limits(),
            || false,
        )
        .unwrap();
    assert!(cycle[2].abs() < 1e-21);
}

#[test]
fn g4_stress_admission_cancellation_darkness_and_replay() {
    let medium = material(1.0, 1.0);
    let other = material(2.0, 1.0);
    let stack =
        PlanarDielectricLayers::try_new(&medium, DielectricResponse::Lossless(&medium), 2.0, &[])
            .unwrap();
    let lines = [constant(&stack, [C64::ONE, C64::ZERO])];
    let bank = AngularSpectrumExposure::try_new(&lines[0], &lines, limits(), || false).unwrap();
    let window = HarmonicWindow::new(0.17, 0.31).unwrap();
    let exact = AngularSpectrumExposureStressLimits {
        maximum_workspace_bytes: bank.stress_workspace_bytes().unwrap(),
    };
    let first = bank
        .stress_field(Side::Incident, 0.0, window, exact, || false)
        .unwrap();
    assert!(
        bank.stress_field(
            Side::Incident,
            0.0,
            window,
            AngularSpectrumExposureStressLimits {
                maximum_workspace_bytes: exact.maximum_workspace_bytes - 1,
            },
            || false
        )
        .is_err()
    );
    for distance in [-1.0, f64::NAN] {
        assert!(
            bank.stress_field(Side::Incident, distance, window, exact, || false)
                .is_err()
        );
    }
    for normal in [[0.0; 3], [0.0, 0.0, 2.0], [f64::NAN, 0.0, 1.0]] {
        assert!(first.points()[0].traction_pa(normal).is_err());
        assert!(first.integrated_traction_n(normal).is_err());
    }
    let mut polls = 0;
    assert!(matches!(
        bank.stress_field(Side::Incident, 0.0, window, exact, || {
            polls += 1;
            polls > 5
        }),
        Err(AngularSpectrumError::Cancelled)
    ));
    assert_eq!(
        first.points(),
        bank.stress_field(Side::Incident, 0.0, window, exact, || false)
            .unwrap()
            .points()
    );
    let none = [];
    let dark =
        AngularSpectrumExposureEnsemble::try_new(&lines[0], None, &none, ensemble_limits(), || {
            false
        })
        .unwrap();
    let zero = dark
        .stress_field(Side::Incident, 0.0, window, stress_limits(), || false)
        .unwrap();
    assert!(
        zero.points()
            .iter()
            .all(|p| p.maxwell_stress_pa == [[0.0; 3]; 3])
    );
    assert_eq!(zero.integrated_stress_n(), [[0.0; 3]; 3]);
    let mismatched =
        PlanarDielectricLayers::try_new(&medium, DielectricResponse::Lossless(&other), 2.0, &[])
            .unwrap();
    let wrong_lines = [constant(&mismatched, [C64::ONE, C64::ZERO])];
    let wrong = AngularSpectrumExposure::try_new(&wrong_lines[0], &wrong_lines, limits(), || false)
        .unwrap();
    assert!(
        wrong
            .stress_field(Side::Incident, 0.0, window, stress_limits(), || false)
            .is_err()
    );
}
