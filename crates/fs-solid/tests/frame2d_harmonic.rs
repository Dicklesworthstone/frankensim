//! G1/G3 native frame resonance: analytic amplitudes, phase, power and rotation.

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_math::c64::C64;
use fs_modal::harmonic::{HarmonicResponseControls, HarmonicResponseError};
use fs_solid::frame2d::{FrameMassFormulation, harmonic::FrameHarmonicError};
use fs_solid::{Frame2d, Geometry, Rayleigh, SectionModel};

fn with_cx<R>(gate: &CancelGate, f: impl FnOnce(&Cx<'_>) -> R) -> R {
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        f(&Cx::new(
            gate,
            arena,
            StreamKey {
                seed: 91,
                kernel_id: 21,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        ))
    })
}

fn controls() -> HarmonicResponseControls {
    HarmonicResponseControls {
        maximum_coordinates: 12,
        maximum_payload_bytes: 1 << 20,
        maximum_relative_residual: 1e-12,
    }
}

fn frame(c: f64, s: f64, guided: bool, formulation: FrameMassFormulation) -> Frame2d {
    let mut frame = Frame2d::with_mass_formulation(Geometry::Linear, formulation);
    frame.add_node(0.0, 0.0);
    frame.add_node(c, s);
    frame.fix(0, [true; 3]);
    if guided {
        frame.fix(1, [false, true, true]);
    }
    frame
        .add_element(0, 1, 5, 6.0, &|| SectionModel::Elastic { ea: 8.0, ei: 2.0 })
        .unwrap();
    frame
}

fn close(actual: C64, expected: C64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * expected.abs().max(1e-30),
        "{actual:?} != {expected:?}"
    );
}

#[test]
fn g1_native_axial_resonance_has_the_original_mass_phase_and_power() {
    with_cx(&CancelGate::new(), |cx| {
        let frame = frame(1.0, 0.0, true, FrameMassFormulation::ConsistentReference);
        let accepted = frame.displacements().to_vec();
        // Unit-length axial interpolation gives tip mass rho*A*L/3 = 2,
        // stiffness EA/L = 8, and authored C = a0*M+a1*K = 1.
        let rayleigh = Rayleigh {
            a0: 0.25,
            a1: 0.0625,
        };
        for omega in [0.0, 0.3, 1.8, 2.0, 2.2, 5.0] {
            let prepared = frame
                .prepare_harmonic_response(omega, rayleigh, controls(), cx)
                .unwrap();
            assert!(std::ptr::eq(prepared.source(), &frame));
            assert_eq!(prepared.free_dofs(), &[3]);
            for force in [C64::new(3.0, -4.0), C64::new(-0.7, 0.2)] {
                let response = prepared.solve(&[force], cx).unwrap();
                let q = force / C64::new(8.0 - 2.0 * omega * omega, -omega);
                let v = C64::new(q.im, -q.re).scale(omega);
                close(response.displacement[0], q, 2e-13);
                close(response.velocity[0], v, 2e-13);
                let heat = 0.5 * v.norm_sq();
                assert!((response.average_viscous_power_w - heat).abs() < 2e-12 * heat.max(1e-30));
                assert!(response.power_balance_residual_w.abs() < 2e-12 * heat.max(1e-30));
                if omega == 2.0 {
                    close(response.velocity[0], force, 2e-13);
                }
            }
        }
        assert_eq!(frame.displacements(), accepted);
    });
}

#[test]
fn g1_complete_consistent_bending_pencil_retains_translation_rotation_coupling() {
    with_cx(&CancelGate::new(), |cx| {
        let frame = frame(1.0, 0.0, false, FrameMassFormulation::ConsistentReference);
        let rayleigh = Rayleigh { a0: 0.1, a1: 0.02 };
        for omega in [0.0, 1.0, 3.0, 9.0] {
            let prepared = frame
                .prepare_harmonic_response(omega, rayleigh, controls(), cx)
                .unwrap();
            assert_eq!(prepared.free_dofs(), &[3, 4, 5]);
            let force = [C64::ZERO, C64::new(1.0, -0.3), C64::new(-0.2, 0.4)];
            let response = prepared.solve(&force, cx).unwrap();
            // Independent exact Hermite integrals for a unit-length beam.
            let entry = |k: f64, m: f64| {
                C64::new(
                    k - omega * omega * m,
                    -omega * (rayleigh.a0 * m + rayleigh.a1 * k),
                )
            };
            let a = entry(24.0, 78.0 / 35.0);
            let b = entry(-12.0, -11.0 / 35.0);
            let d = entry(8.0, 2.0 / 35.0);
            let determinant = a * d - b * b;
            close(
                response.displacement[1],
                (d * force[1] - b * force[2]) / determinant,
                3e-12,
            );
            close(
                response.displacement[2],
                (a * force[2] - b * force[1]) / determinant,
                3e-12,
            );
            assert!(response.relative_residual < 1e-12);
        }
    });
}

#[test]
fn g3_oblique_frame_keeps_reciprocity_and_physical_rotation() {
    with_cx(&CancelGate::new(), |cx| {
        let reference = frame(1.0, 0.0, false, FrameMassFormulation::ConsistentReference);
        let rotated = frame(0.6, 0.8, false, FrameMassFormulation::ConsistentReference);
        let rayleigh = Rayleigh { a0: 0.2, a1: 0.03 };
        let base = reference
            .prepare_harmonic_response(3.0, rayleigh, controls(), cx)
            .unwrap();
        let other = rotated
            .prepare_harmonic_response(3.0, rayleigh, controls(), cx)
            .unwrap();
        let force = [C64::new(0.2, -0.1), C64::new(1.0, 0.3), C64::new(-0.2, 0.4)];
        let turn = |x: &[C64]| {
            [
                x[0].scale(0.6) - x[1].scale(0.8),
                x[0].scale(0.8) + x[1].scale(0.6),
                x[2],
            ]
        };
        let expected = base.solve(&force, cx).unwrap();
        let actual = other.solve(&turn(&force), cx).unwrap();
        for (got, want) in actual.displacement.iter().zip(turn(&expected.displacement)) {
            close(*got, want, 5e-12);
        }
        assert!(
            (actual.average_input_power_w / expected.average_input_power_w - 1.0).abs() < 5e-12
        );
        let x = other.solve(&[C64::ONE, C64::ZERO, C64::ZERO], cx).unwrap();
        let y = other.solve(&[C64::ZERO, C64::ONE, C64::ZERO], cx).unwrap();
        close(x.displacement[1], y.displacement[0], 3e-13);
    });
}

#[test]
fn g0_refusals_preserve_original_state_and_do_not_invent_mass_or_loss() {
    let gate = CancelGate::new();
    with_cx(&gate, |cx| {
        let frame = frame(1.0, 0.0, true, FrameMassFormulation::ConsistentReference);
        let lossless = Rayleigh { a0: 0.0, a1: 0.0 };
        let accepted = frame.displacements().to_vec();
        assert!(matches!(
            frame.prepare_harmonic_response(2.0, lossless, controls(), cx),
            Err(FrameHarmonicError::Harmonic(
                HarmonicResponseError::Singular { .. }
            ))
        ));
        for bad in [
            HarmonicResponseControls {
                maximum_coordinates: 5,
                ..controls()
            },
            HarmonicResponseControls {
                maximum_payload_bytes: 1,
                ..controls()
            },
        ] {
            assert!(matches!(
                frame.prepare_harmonic_response(1.0, lossless, bad, cx),
                Err(FrameHarmonicError::Harmonic(
                    HarmonicResponseError::Budget { .. }
                ))
            ));
        }
        let prepared = frame
            .prepare_harmonic_response(1.0, lossless, controls(), cx)
            .unwrap();
        assert!(prepared.solve(&[C64::new(f64::NAN, 0.0)], cx).is_err());
        let a = prepared.solve(&[C64::ONE], cx).unwrap();
        let b = prepared.solve(&[C64::ONE], cx).unwrap();
        assert_eq!(a.displacement, b.displacement);
        assert_eq!(frame.displacements(), accepted);
        let massless = self::frame(1.0, 0.0, false, FrameMassFormulation::Lumped);
        assert!(
            massless
                .prepare_harmonic_response(1.0, lossless, controls(), cx)
                .is_err()
        );
        gate.request();
        assert!(matches!(
            prepared.solve(&[C64::ONE], cx),
            Err(FrameHarmonicError::Harmonic(
                HarmonicResponseError::Cancelled
            ))
        ));
        assert!(matches!(
            frame.prepare_harmonic_response(1.0, lossless, controls(), cx),
            Err(FrameHarmonicError::Harmonic(
                HarmonicResponseError::Cancelled
            ))
        ));
    });
}
