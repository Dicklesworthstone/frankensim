//! G1/G3: finite-window stress and storage in each original exterior material.
use super::*;
use fs_em::plane_wave::angular_spectrum::{
    AngularSpectrumExposureMomentumField, AngularSpectrumExposureMomentumLimits,
    AngularSpectrumExposureStressField, AngularSpectrumExposureStressLimits,
    AngularSpectrumFactorStatistics, AngularSpectrumIlluminationLimits,
    AngularSpectrumIlluminationSide as Face,
};
use fs_fft::harmonic_window::HarmonicWindowProfile;

// Independent averages of cos²(omega*t) and its actual time derivative.
// The triangular shutter has a squared sinc, not an endpoint-difference rate.
fn moments(omega: f64, window: HarmonicWindow, profile: HarmonicWindowProfile) -> (f64, f64) {
    let sinc = |x: f64| if x == 0.0 { 1.0 } else { x.sin() / x };
    let response = match profile {
        HarmonicWindowProfile::Uniform => sinc(omega * window.duration_s()),
        HarmonicWindowProfile::Triangular => sinc(0.5 * omega * window.duration_s()).powi(2),
    };
    let phase = 2.0 * omega * (window.start_s() + 0.5 * window.duration_s());
    (
        0.5 * (1.0 + response * phase.cos()),
        -omega * response * phase.sin(),
    )
}

fn check(
    stress: &AngularSpectrumExposureStressField,
    momentum: &AngularSpectrumExposureMomentumField,
    epsilon: f64,
    mu: f64,
    ex: f64,
    hy: f64,
    square: f64,
    rate: f64,
) {
    // Only Ex and Hy are present in this independent normally incident field.
    let electric = epsilon * ex * ex;
    let magnetic = mu * hy * hy;
    let diagonal = [
        0.5 * (electric - magnetic) * square,
        0.5 * (magnetic - electric) * square,
        -0.5 * (electric + magnetic) * square,
    ];
    let density = epsilon * mu * ex * hy * square;
    let storage = epsilon * mu * ex * hy * rate;
    let close =
        |actual: f64, expected: f64| near(actual, expected, 3e-12 * expected.abs().max(1.0));
    for point in stress.points() {
        for i in 0..3 {
            for j in 0..3 {
                close(
                    point.maxwell_stress_pa[i][j],
                    if i == j { diagonal[i] } else { 0.0 },
                );
            }
        }
    }
    for point in momentum.points() {
        for c in 0..3 {
            close(
                point.minkowski_density_n_s_m3[c],
                if c == 2 { density } else { 0.0 },
            );
            close(
                point.minkowski_rate_n_m3[c],
                if c == 2 { storage } else { 0.0 },
            );
        }
    }
    let grid = stress.grid();
    let area = (grid.width * grid.height) as f64 * grid.spacing_m[0] * grid.spacing_m[1];
    let tensor = stress.integrated_stress_n();
    for i in 0..3 {
        for j in 0..3 {
            close(tensor[i][j], if i == j { area * diagonal[i] } else { 0.0 });
        }
    }
    for c in 0..3 {
        close(
            momentum.momentum_per_depth_n_s_m()[c],
            if c == 2 { area * density } else { 0.0 },
        );
        close(
            momentum.momentum_rate_per_depth_n_m()[c],
            if c == 2 { area * storage } else { 0.0 },
        );
    }
}

#[test]
fn g1_refracting_interfaces_use_actual_electric_and_magnetic_exteriors() {
    let left = material(1.0, 1.0);
    let omega = 2.0;
    for (epsilon, mu) in [(4.0_f64, 1.0_f64), (1.0, 4.0), (4.0, 4.0)] {
        let right = material(epsilon, mu);
        let stack = PlanarDielectricLayers::try_new(
            &left,
            DielectricResponse::Lossless(&right),
            omega,
            &[],
        )
        .unwrap();
        let lines = [constant(&stack, [C64::ONE, C64::ZERO])];
        let bank = AngularSpectrumExposure::try_new(&lines[0], &lines, limits(), || false).unwrap();
        let stress_limits = AngularSpectrumExposureStressLimits {
            maximum_workspace_bytes: bank.stress_workspace_bytes().unwrap(),
        };
        let momentum_limits = AngularSpectrumExposureMomentumLimits {
            maximum_workspace_bytes: bank.momentum_workspace_bytes().unwrap(),
        };
        // Independent normal-incidence continuity of tangential E and H.
        let y = (epsilon / mu).sqrt();
        let r = (1.0 - y) / (1.0 + y);
        let t = 2.0 / (1.0 + y);
        for profile in [
            HarmonicWindowProfile::Uniform,
            HarmonicWindowProfile::Triangular,
        ] {
            let window = HarmonicWindow::with_profile(0.17, 0.31, profile).unwrap();
            let (square, rate) = moments(omega, window, profile);
            for (side, e, m, ex, hy) in [
                (Side::Incident, 1.0, 1.0, 1.0, 1.0),
                (Side::Reflected, 1.0, 1.0, r, -r),
                (Side::TotalIncidentSide, 1.0, 1.0, 1.0 + r, 1.0 - r),
                (Side::Transmitted, epsilon, mu, t, y * t),
            ] {
                let stress = bank
                    .stress_field(side, 0.0, window, stress_limits, || false)
                    .unwrap();
                let momentum = bank
                    .momentum_field(side, 0.0, window, momentum_limits, || false)
                    .unwrap();
                check(&stress, &momentum, e, m, ex, hy, square, rate);
            }
            let area = 4.0 * grid().spacing_m[0] * grid().spacing_m[1];
            let surface = bank
                .boundary_stress_resultant_n([0.0; 2], window, stress_limits, || false)
                .unwrap();
            let expected = area * (1.0 + r * r - epsilon * t * t) * square;
            near(surface[2], expected, 3e-12 * expected.abs().max(1.0));
            near(surface[0], 0.0, 1e-12);
            near(surface[1], 0.0, 1e-12);
        }
    }
}

#[test]
fn g3_joint_faces_and_real_proper_factors_keep_each_medium_and_direction() {
    let left = material(1.0, 1.0);
    let right = material(4.0, 1.0);
    let omega = 2.0;
    let stack =
        PlanarDielectricLayers::try_new(&left, DielectricResponse::Lossless(&right), omega, &[])
            .unwrap();
    let a = 0.6;
    let b = -0.4;
    let beam = constant(&stack, [C64::from_re(a), C64::ZERO]);
    let lines = [beam
        .prepare_illumination(
            &[[C64::from_re(b), C64::ZERO]; 4],
            AngularSpectrumIlluminationLimits {
                maximum_layer_mode_pairs: 4096,
                maximum_workspace_bytes: 8 << 20,
                maximum_fft_work: 1 << 20,
            },
            || false,
        )
        .unwrap()];
    let mean =
        AngularSpectrumExposure::try_new_illumination(&beam, &lines, limits(), || false).unwrap();
    let factors = [
        AngularSpectrumExposure::try_new_illumination(&beam, &lines, limits(), || false).unwrap(),
        AngularSpectrumExposure::try_new_illumination(&beam, &lines, limits(), || false).unwrap(),
    ];
    let statistics = [
        AngularSpectrumFactorStatistics::Real,
        AngularSpectrumFactorStatistics::ProperComplex,
    ];
    let ensemble = AngularSpectrumExposureEnsemble::try_new_with_statistics(
        &beam,
        Some(&mean),
        &factors,
        &statistics,
        ensemble_limits(),
        || false,
    )
    .unwrap();
    let stress_limits = AngularSpectrumExposureStressLimits {
        maximum_workspace_bytes: ensemble.stress_workspace_bytes().unwrap(),
    };
    let momentum_limits = AngularSpectrumExposureMomentumLimits {
        maximum_workspace_bytes: ensemble.momentum_workspace_bytes().unwrap(),
    };
    let y = 2.0;
    let r = -1.0 / 3.0;
    let outgoing_left = r * a + (1.0 - r) * b;
    let outgoing_right = (1.0 + r) * a - r * b;
    for profile in [
        HarmonicWindowProfile::Uniform,
        HarmonicWindowProfile::Triangular,
    ] {
        let window = HarmonicWindow::with_profile(0.17, 0.31, profile).unwrap();
        let (square, rate) = moments(omega, window, profile);
        // Mean plus an independent unit-variance real factor contribute two
        // real copies. The proper factor contributes mean cos²=1/2 and rate=0.
        for (side, epsilon, ex, hy) in [
            (Face::IncomingLeft, 1.0, a, a),
            (Face::OutgoingLeft, 1.0, outgoing_left, -outgoing_left),
            (Face::TotalLeft, 1.0, a + outgoing_left, a - outgoing_left),
            (Face::IncomingRight, 4.0, b, -y * b),
            (Face::OutgoingRight, 4.0, outgoing_right, y * outgoing_right),
            (
                Face::TotalRight,
                4.0,
                b + outgoing_right,
                y * (outgoing_right - b),
            ),
        ] {
            let stress = ensemble
                .illumination_stress_field(side, 0.0, window, stress_limits, || false)
                .unwrap();
            let momentum = ensemble
                .illumination_momentum_field(side, 0.0, window, momentum_limits, || false)
                .unwrap();
            check(
                &stress,
                &momentum,
                epsilon,
                1.0,
                ex,
                hy,
                2.0 * square + 0.5,
                2.0 * rate,
            );
        }
    }
}

#[test]
fn g4_lossy_or_dispersive_exit_refuses_only_unsupported_observation_sides() {
    let left = material(1.0, 1.0);
    let (dielectric, conductor) = heat_material();
    for response in [
        DielectricResponse::Debye(&dielectric),
        DielectricResponse::Ohmic(&conductor),
    ] {
        let stack = PlanarDielectricLayers::try_new(&left, response, 2.0, &[]).unwrap();
        let lines = [constant(&stack, [C64::ONE, C64::ZERO])];
        let bank = AngularSpectrumExposure::try_new(&lines[0], &lines, limits(), || false).unwrap();
        let window = HarmonicWindow::new(0.17, 0.31).unwrap();
        let stress_limits = AngularSpectrumExposureStressLimits {
            maximum_workspace_bytes: bank.stress_workspace_bytes().unwrap(),
        };
        let momentum_limits = AngularSpectrumExposureMomentumLimits {
            maximum_workspace_bytes: bank.momentum_workspace_bytes().unwrap(),
        };
        let before = bank
            .stress_field(Side::TotalIncidentSide, 0.0, window, stress_limits, || {
                false
            })
            .unwrap();
        let momentum_before = bank
            .momentum_field(
                Side::TotalIncidentSide,
                0.0,
                window,
                momentum_limits,
                || false,
            )
            .unwrap();
        assert!(
            bank.stress_field(Side::Transmitted, 0.0, window, stress_limits, || false,)
                .is_err()
        );
        assert!(
            bank.momentum_field(Side::Transmitted, 0.0, window, momentum_limits, || false,)
                .is_err()
        );
        for side in [Face::IncomingRight, Face::OutgoingRight, Face::TotalRight] {
            assert!(
                bank.illumination_stress_field(side, 0.0, window, stress_limits, || false,)
                    .is_err()
            );
            assert!(
                bank.illumination_momentum_field(side, 0.0, window, momentum_limits, || false,)
                    .is_err()
            );
        }
        assert!(
            bank.boundary_stress_resultant_n([0.0; 2], window, stress_limits, || false,)
                .is_err()
        );
        assert_eq!(
            before.points(),
            bank.stress_field(Side::TotalIncidentSide, 0.0, window, stress_limits, || {
                false
            },)
                .unwrap()
                .points()
        );
        assert_eq!(
            momentum_before.points(),
            bank.momentum_field(
                Side::TotalIncidentSide,
                0.0,
                window,
                momentum_limits,
                || false,
            )
            .unwrap()
            .points()
        );
    }
}
