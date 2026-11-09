//! Admission reserves bounded refinement work only for coupled factors.

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_math::c64::C64;
use fs_modal::harmonic::{
    HarmonicFrequencyControls, HarmonicResponseControls, HarmonicResponseError,
    PreparedHarmonicResponse,
};

#[test]
fn coupled_frequency_admission_reserves_refinement_without_inflating_scalar_work() {
    let gate = CancelGate::new();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 91,
                kernel_id: 23,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        let controls = HarmonicResponseControls {
            maximum_coordinates: 2,
            maximum_payload_bytes: 1 << 20,
            maximum_relative_residual: 2e-12,
        };
        let frequency_controls = |maximum_work| HarmonicFrequencyControls {
            maximum_additional_payload_bytes: 1 << 20,
            maximum_work,
        };
        let force = [C64::ONE, C64::ZERO];
        let coupled = PreparedHarmonicResponse::try_new(
            &[1.0, 0.0, 0.0, 1.0],
            &[0.2, 0.0, 0.0, 0.0],
            &[2.0, -1.0, -1.0, 4.0],
            2,
            1.7,
            controls,
            &cx,
        )
        .unwrap();
        // Two coupled coordinates cost 64*n^2 + 256*n = 768 units per
        // pass. Admit all three bounded residual-refinement passes up front.
        assert!(matches!(
            coupled.solve_with_frequency_derivative(
                &force,
                None,
                None,
                frequency_controls(2303),
                &cx,
            ),
            Err(HarmonicResponseError::WorkBudget {
                requested: 2304,
                cap: 2303,
            })
        ));
        let response = coupled
            .solve_with_frequency_derivative(&force, None, None, frequency_controls(2304), &cx)
            .unwrap();
        assert!(response.relative_derivative_residual <= controls.maximum_relative_residual);

        // Scalar factors have no coupled refinement loop and retain 320*n.
        let diagonal = PreparedHarmonicResponse::try_diagonal_modal(
            &[1.0, 2.0],
            &[0.2, 0.4],
            1.7,
            controls,
            &cx,
        )
        .unwrap();
        assert!(matches!(
            diagonal.solve_with_frequency_derivative(
                &force,
                None,
                None,
                frequency_controls(639),
                &cx,
            ),
            Err(HarmonicResponseError::WorkBudget {
                requested: 640,
                cap: 639,
            })
        ));
        let response = diagonal
            .solve_with_frequency_derivative(&force, None, None, frequency_controls(640), &cx)
            .unwrap();
        assert!(response.relative_derivative_residual <= controls.maximum_relative_residual);
    });
}
