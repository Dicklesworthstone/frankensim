//! G1/G3/G4 joint stress and momentum through original Maxwell source owners.
use super::*;
use fs_em::plane_wave::angular_spectrum::{
    AngularSpectrumExposureMomentumField, AngularSpectrumExposureMomentumLimits,
    AngularSpectrumExposureStressField, AngularSpectrumExposureStressLimits,
    AngularSpectrumExposureStressMomentumField, AngularSpectrumExposureStressMomentumLimits,
    AngularSpectrumFactorStatistics, AngularSpectrumIlluminationLimits,
    AngularSpectrumIlluminationSide as Face, PreparedAngularSpectrumIllumination,
};
use fs_fft::harmonic_window::HarmonicWindowProfile::{Triangular, Uniform};

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

fn same_observables(
    joint: &AngularSpectrumExposureStressMomentumField,
    stress: &AngularSpectrumExposureStressField,
    momentum: &AngularSpectrumExposureMomentumField,
) {
    assert_eq!(joint.grid(), stress.grid());
    assert_eq!(joint.grid(), momentum.grid());
    assert_eq!(joint.stress().points(), stress.points());
    assert_eq!(
        joint.stress().integrated_stress_n(),
        stress.integrated_stress_n()
    );
    assert_eq!(
        joint
            .stress()
            .integrated_traction_n([0.6, 0.0, 0.8])
            .unwrap(),
        stress.integrated_traction_n([0.6, 0.0, 0.8]).unwrap()
    );
    assert_eq!(joint.momentum().points(), momentum.points());
    assert_eq!(
        joint.momentum().momentum_per_depth_n_s_m(),
        momentum.momentum_per_depth_n_s_m()
    );
    assert_eq!(
        joint.momentum().momentum_rate_per_depth_n_m(),
        momentum.momentum_rate_per_depth_n_m()
    );
}

#[test]
fn g3_joint_query_preserves_mean_real_proper_and_both_face_interference() {
    // Different exterior owners also exercise the observed-side constitutive
    // selection; the joint query must inherit each original field's SI law.
    let left_medium = material(1.3, 0.7);
    let right_medium = material(2.1, 1.2);
    let film = material(1.8, 0.9);
    let stacks = [2.0, 2.3].map(|omega| {
        PlanarDielectricLayers::try_new(
            &left_medium,
            DielectricResponse::Lossless(&right_medium),
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
    let factors = [
        AngularSpectrumExposure::try_new_illumination(&left[0], &lines, limits(), || false)
            .unwrap(),
        AngularSpectrumExposure::try_new_illumination(&left[0], &lines, limits(), || false)
            .unwrap(),
    ];
    let statistics = [
        AngularSpectrumFactorStatistics::Real,
        AngularSpectrumFactorStatistics::ProperComplex,
    ];
    let ensemble = AngularSpectrumExposureEnsemble::try_new_with_statistics(
        &left[0],
        Some(&mean),
        &factors,
        &statistics,
        ensemble_limits(),
        || false,
    )
    .unwrap();
    let joint_limits = AngularSpectrumExposureStressMomentumLimits {
        maximum_workspace_bytes: ensemble.stress_momentum_workspace_bytes().unwrap(),
    };
    let stress_limits = AngularSpectrumExposureStressLimits {
        maximum_workspace_bytes: ensemble.stress_workspace_bytes().unwrap(),
    };
    let momentum_limits = AngularSpectrumExposureMomentumLimits {
        maximum_workspace_bytes: ensemble.momentum_workspace_bytes().unwrap(),
    };
    assert_eq!(
        joint_limits.maximum_workspace_bytes,
        stress_limits.maximum_workspace_bytes + momentum_limits.maximum_workspace_bytes
            - ensemble.workspace_bytes()
    );
    for profile in [Uniform, Triangular] {
        let window = HarmonicWindow::with_profile(0.17, 0.12, profile).unwrap();
        for face in [
            Face::IncomingLeft,
            Face::OutgoingLeft,
            Face::TotalLeft,
            Face::IncomingRight,
            Face::OutgoingRight,
            Face::TotalRight,
        ] {
            let joint = ensemble
                .illumination_stress_momentum_field(face, 0.23, window, joint_limits, || false)
                .unwrap();
            let stress = ensemble
                .illumination_stress_field(face, 0.23, window, stress_limits, || false)
                .unwrap();
            let momentum = ensemble
                .illumination_momentum_field(face, 0.23, window, momentum_limits, || false)
                .unwrap();
            same_observables(&joint, &stress, &momentum);
            // The deterministic two-sided owner takes the same shared route.
            let deterministic = factors[0]
                .illumination_stress_momentum_field(
                    face,
                    0.23,
                    window,
                    AngularSpectrumExposureStressMomentumLimits {
                        maximum_workspace_bytes: factors[0]
                            .stress_momentum_workspace_bytes()
                            .unwrap(),
                    },
                    || false,
                )
                .unwrap();
            same_observables(
                &deterministic,
                &factors[0]
                    .illumination_stress_field(face, 0.23, window, stress_limits, || false)
                    .unwrap(),
                &factors[0]
                    .illumination_momentum_field(face, 0.23, window, momentum_limits, || false)
                    .unwrap(),
            );
        }
        for side in [
            Side::Incident,
            Side::Reflected,
            Side::TotalIncidentSide,
            Side::Transmitted,
        ] {
            let joint = ensemble
                .stress_momentum_field(side, 0.23, window, joint_limits, || false)
                .unwrap();
            same_observables(
                &joint,
                &ensemble
                    .stress_field(side, 0.23, window, stress_limits, || false)
                    .unwrap(),
                &ensemble
                    .momentum_field(side, 0.23, window, momentum_limits, || false)
                    .unwrap(),
            );
        }
    }
}

#[test]
fn g1_joint_nyquist_wave_keeps_exact_cell_stress_and_momentum_storage() {
    let medium = material(1.0, 1.0);
    let omega = 2.0;
    let stack =
        PlanarDielectricLayers::try_new(&medium, DielectricResponse::Lossless(&medium), omega, &[])
            .unwrap();
    let samples = (0..4)
        .map(|i| [C64::ZERO, C64::from_re(if i % 2 == 0 { 1.0 } else { -1.0 })])
        .collect::<Vec<_>>();
    let lines = [prepare(&stack, &samples)];
    let bank = AngularSpectrumExposure::try_new(&lines[0], &lines, limits(), || false).unwrap();
    let exact = AngularSpectrumExposureStressMomentumLimits {
        maximum_workspace_bytes: bank.stress_momentum_workspace_bytes().unwrap(),
    };
    let cell = grid();
    let kx = -core::f64::consts::PI / cell.spacing_m[0];
    let direction = [kx / omega, 0.0, (omega * omega - kx * kx).sqrt() / omega];
    let area = (cell.width * cell.height) as f64 * cell.spacing_m[0] * cell.spacing_m[1];
    for profile in [Uniform, Triangular] {
        let window = HarmonicWindow::with_profile(0.0, 0.2 / omega, profile).unwrap();
        let joint = bank
            .stress_momentum_field(Side::Incident, 0.0, window, exact, || false)
            .unwrap();
        for i in 0..3 {
            near(
                joint.momentum().momentum_per_depth_n_s_m()[i],
                0.5 * area * direction[i],
                3e-13 * area,
            );
            near(
                joint.momentum().momentum_rate_per_depth_n_m()[i],
                0.0,
                3e-13 * area,
            );
            for j in 0..3 {
                near(
                    joint.stress().integrated_stress_n()[i][j],
                    -0.5 * area * direction[i] * direction[j],
                    3e-13 * area,
                );
            }
        }
        let sampled_density = joint
            .momentum()
            .points()
            .iter()
            .map(|p| p.minkowski_density_n_s_m3[2])
            .sum::<f64>()
            * cell.spacing_m[0]
            * cell.spacing_m[1];
        assert!(sampled_density > 1.9 * joint.momentum().momentum_per_depth_n_s_m()[2]);
        assert!(
            joint
                .momentum()
                .points()
                .iter()
                .any(|p| p.minkowski_rate_n_m3[2].abs() > 0.1)
        );
        let (stress, momentum) = joint.into_fields();
        assert_eq!(stress.grid(), momentum.grid());
        near(momentum.momentum_rate_per_depth_n_m()[2], 0.0, 3e-13 * area);
    }
}

#[test]
fn g4_joint_budget_cancellation_darkness_and_retry_publish_complete_pairs() {
    let medium = material(1.0, 1.0);
    let stack =
        PlanarDielectricLayers::try_new(&medium, DielectricResponse::Lossless(&medium), 2.0, &[])
            .unwrap();
    let lines = [constant(&stack, [C64::ONE, C64::new(0.2, 0.3)])];
    let bank = AngularSpectrumExposure::try_new(&lines[0], &lines, limits(), || false).unwrap();
    let exact = AngularSpectrumExposureStressMomentumLimits {
        maximum_workspace_bytes: bank.stress_momentum_workspace_bytes().unwrap(),
    };
    let separate =
        bank.stress_workspace_bytes().unwrap() + bank.momentum_workspace_bytes().unwrap();
    assert_eq!(
        exact.maximum_workspace_bytes,
        separate - bank.workspace_bytes()
    );
    assert!(exact.maximum_workspace_bytes < separate);
    let window = HarmonicWindow::new(0.17, 0.31).unwrap();
    let mut complete_polls = 0;
    let baseline = bank
        .stress_momentum_field(Side::Incident, 0.0, window, exact, || {
            complete_polls += 1;
            false
        })
        .unwrap();
    assert!(
        bank.stress_momentum_field(
            Side::Incident,
            0.0,
            window,
            AngularSpectrumExposureStressMomentumLimits {
                maximum_workspace_bytes: exact.maximum_workspace_bytes - 1,
            },
            || false,
        )
        .is_err()
    );
    for distance in [-1.0, f64::NAN, f64::INFINITY] {
        assert!(
            bank.stress_momentum_field(Side::Incident, distance, window, exact, || false,)
                .is_err()
        );
    }
    for stop in [1, complete_polls / 2, complete_polls] {
        let mut polls = 0;
        assert!(matches!(
            bank.stress_momentum_field(Side::Incident, 0.0, window, exact, || {
                polls += 1;
                polls >= stop
            },),
            Err(AngularSpectrumError::Cancelled)
        ));
    }
    let replay = bank
        .stress_momentum_field(Side::Incident, 0.0, window, exact, || false)
        .unwrap();
    same_observables(&replay, baseline.stress(), baseline.momentum());
    let none = [];
    let dark =
        AngularSpectrumExposureEnsemble::try_new(&lines[0], None, &none, ensemble_limits(), || {
            false
        })
        .unwrap();
    let zero = dark
        .illumination_stress_momentum_field(
            Face::TotalRight,
            0.0,
            window,
            AngularSpectrumExposureStressMomentumLimits {
                maximum_workspace_bytes: dark.stress_momentum_workspace_bytes().unwrap(),
            },
            || false,
        )
        .unwrap();
    assert_eq!(zero.stress().integrated_stress_n(), [[0.0; 3]; 3]);
    assert_eq!(zero.momentum().momentum_per_depth_n_s_m(), [0.0; 3]);
    assert_eq!(zero.momentum().momentum_rate_per_depth_n_m(), [0.0; 3]);
    assert!(
        zero.stress()
            .points()
            .iter()
            .all(|p| p.maxwell_stress_pa == [[0.0; 3]; 3])
    );
    assert!(
        zero.momentum()
            .points()
            .iter()
            .all(|p| p.minkowski_density_n_s_m3 == [0.0; 3] && p.minkowski_rate_n_m3 == [0.0; 3])
    );
}
