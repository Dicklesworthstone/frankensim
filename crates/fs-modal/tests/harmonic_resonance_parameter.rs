//! Native physical jets connected to actual driven peak and notch sensitivities.

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_math::c64::C64;
use fs_modal::harmonic::resonance::{
    ResonanceControls, TransferExtremumKind, refine_transfer_extremum_with_parameter_derivative,
};
use fs_modal::harmonic::{
    HarmonicFrequencyControls, HarmonicFrequencyParameterControls, HarmonicResponseControls,
    PreparedHarmonicResponse,
};
use fs_modal::quadratic_modes::QuadraticParameterDirection;

fn with_cx(f: impl FnOnce(&Cx<'_>)) {
    let gate = CancelGate::new();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        f(&Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 37,
                kernel_id: 25,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        ))
    })
}

fn response_controls() -> HarmonicResponseControls {
    HarmonicResponseControls {
        maximum_coordinates: 4,
        maximum_payload_bytes: 1 << 20,
        maximum_relative_residual: 2e-11,
    }
}
fn first_controls() -> HarmonicFrequencyControls {
    HarmonicFrequencyControls {
        maximum_additional_payload_bytes: 1 << 20,
        maximum_work: 1_000_000,
    }
}
fn jet_controls() -> HarmonicFrequencyParameterControls {
    HarmonicFrequencyParameterControls {
        maximum_additional_payload_bytes: 1 << 20,
        maximum_work: 1_000_000,
    }
}
fn query_controls() -> ResonanceControls {
    ResonanceControls {
        absolute_angular_tolerance_rad_s: 1e-12,
        relative_angular_tolerance: 1e-12,
        maximum_evaluations: 200,
    }
}

#[test]
fn actual_displacement_and_velocity_peak_motion_follow_mass_stiffness_and_damping() {
    with_cx(|cx| {
        let (m, c, k) = (2.0_f64, 0.3_f64, 8.0_f64);
        let displacement_peak = (k / m - c * c / (2.0 * m * m)).sqrt();
        let velocity_peak = (k / m).sqrt();
        let force = [C64::new(1.3, -0.7)];
        for (dm, dc, dk) in [(1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0)] {
            let direction = QuadraticParameterDirection {
                mass: &[dm],
                damping: &[dc],
                stiffness: &[dk],
            };
            for velocity in [false, true] {
                let prepare = |w| {
                    PreparedHarmonicResponse::try_new(
                        &[m],
                        &[c],
                        &[k],
                        1,
                        w,
                        response_controls(),
                        cx,
                    )
                };
                let result = refine_transfer_extremum_with_parameter_derivative(
                    (1.7, 2.3),
                    query_controls(),
                    cx,
                    |w| {
                        let first = prepare(w)?.solve_with_frequency_derivative(
                            &force,
                            None,
                            None,
                            first_controls(),
                            cx,
                        )?;
                        if velocity {
                            first.velocity_transfer(&[1.0])
                        } else {
                            first.displacement_transfer(&[1.0])
                        }
                    },
                    |w| {
                        let jet = prepare(w)?.solve_with_frequency_parameter_derivative(
                            &force,
                            direction,
                            None,
                            None,
                            None,
                            None,
                            jet_controls(),
                            cx,
                        )?;
                        if velocity {
                            jet.velocity_transfer_parameter(&[1.0])
                        } else {
                            jet.displacement_transfer_parameter(&[1.0])
                        }
                    },
                )
                .unwrap();
                let (expected_w, expected_slope) = if velocity {
                    (
                        velocity_peak,
                        (-dm * k / (m * m) + dk / m) / (2.0 * velocity_peak),
                    )
                } else {
                    (
                        displacement_peak,
                        (dm * (-k / (m * m) + c * c / (m * m * m)) - dc * c / (m * m) + dk / m)
                            / (2.0 * displacement_peak),
                    )
                };
                assert_eq!(result.extremum.kind, TransferExtremumKind::Maximum);
                assert!((result.extremum.angular_frequency_rad_s - expected_w).abs() < 4e-12);
                assert!((result.angular_frequency_derivative_rad_s - expected_slope).abs() < 4e-11);
                assert!(result.normalized_curvature < 0.0);
            }
        }
    });
}

#[test]
fn a_physical_cancellation_notch_tracks_its_local_spring_and_mass() {
    with_cx(|cx| {
        let mass = [1.0, 0.0, 0.0, 1.0];
        let damping = [0.2, 0.0, 0.0, 0.0];
        let stiffness = [2.0, -1.0, -1.0, 4.0];
        let zero = [0.0; 4];
        let local = [0.0, 0.0, 0.0, 1.0];
        let force = [C64::ONE, C64::ZERO];
        let prepare = |w| {
            PreparedHarmonicResponse::try_new(
                &mass,
                &damping,
                &stiffness,
                2,
                w,
                response_controls(),
                cx,
            )
        };
        for (dm, dk, expected_slope) in [(&zero, &local, 0.25), (&local, &zero, -1.0)] {
            let direction = QuadraticParameterDirection {
                mass: dm,
                damping: &zero,
                stiffness: dk,
            };
            let result = refine_transfer_extremum_with_parameter_derivative(
                (1.95, 2.05),
                query_controls(),
                cx,
                |w| {
                    prepare(w)?
                        .solve_with_frequency_derivative(&force, None, None, first_controls(), cx)?
                        .displacement_transfer(&[1.0, 0.0])
                },
                |w| {
                    prepare(w)?
                        .solve_with_frequency_parameter_derivative(
                            &force,
                            direction,
                            None,
                            None,
                            None,
                            None,
                            jet_controls(),
                            cx,
                        )?
                        .displacement_transfer_parameter(&[1.0, 0.0])
                },
            )
            .unwrap();
            assert_eq!(result.extremum.kind, TransferExtremumKind::Minimum);
            assert!((result.extremum.angular_frequency_rad_s - 2.0).abs() < 3e-12);
            assert!((result.angular_frequency_derivative_rad_s - expected_slope).abs() < 3e-11);
            assert!(result.extremum.sample.amplitude.abs() < 2e-11);
            assert!(result.amplitude_derivative.abs() < 2e-9);
            assert!(result.normalized_curvature > 0.0);
        }
        // The vanishing receiver is a true internal cancellation: the other
        // physical coordinate still moves and the full system is nonsingular.
        let at_notch = prepare(2.0).unwrap().solve(&force, cx).unwrap();
        assert!(at_notch.displacement[0].abs() < 1e-14);
        assert!((at_notch.displacement[1] + C64::ONE).abs() < 1e-14);
    });
}
