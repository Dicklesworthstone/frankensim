//! Spatial refinement of a computed latent-heat front at fixed time resolution.

mod support;

#[allow(dead_code)]
#[path = "../examples/enthalpy_stefan.rs"]
mod stefan;

#[test]
fn computed_stefan_front_refines_spatially_and_conserves_applied_heat() {
    support::with_cx(|cx| {
        // Keep time resolution fixed: this comparison makes no temporal or
        // formal convergence-order claim. Only the axial mesh is refined.
        let coarse = stefan::run(cx, 40, 480).unwrap();
        let fine = stefan::run(cx, 80, 480).unwrap();

        // Independent similarity data: alpha=1, lambda=0.5, t_end=0.16
        // give s=2*lambda*sqrt(alpha*t_end)=0.4 m. Integrating the imposed
        // q_in(t)=0.5*exp(0.25)/sqrt(t) from 0.04 to 0.16 over area 0.02^2
        // gives this heat input. Neither uses the computed phase field.
        let reference_front = 0.4;
        let applied_heat = 0.00008 * 0.25_f64.exp();
        for (cells, result) in [(40, &coarse), (80, &fine)] {
            assert!((result.reference_front_m - reference_front).abs() < 1e-14);
            assert!((result.applied_heat_j - applied_heat).abs() < 1e-16);
            let relative_energy_defect =
                (result.stored_energy_change_j - applied_heat).abs() / applied_heat;
            assert!(
                relative_energy_defect < 5e-8,
                "{cells} cells: relative energy defect {relative_energy_defect}"
            );
            assert!(result.krylov_iterations > 0);
            assert!(result.temperature_rms_error_k > 0.0);
            println!(
                "cells={cells}: front={:.9} m; molten length={:.9} m; RMS T error={:.9e} K; relative energy defect={:.9e}; Krylov work={}",
                result.front_m,
                result.equivalent_molten_length_m,
                result.temperature_rms_error_k,
                relative_energy_defect,
                result.krylov_iterations,
            );
        }

        let coarse_front_error = (coarse.front_m - reference_front).abs();
        let fine_front_error = (fine.front_m - reference_front).abs();
        assert!(
            fine_front_error < 2e-4,
            "fine front error {fine_front_error}"
        );
        assert!(
            fine_front_error < 0.6 * coarse_front_error,
            "front errors: coarse={coarse_front_error}, fine={fine_front_error}"
        );
        let coarse_molten_error = (coarse.equivalent_molten_length_m - reference_front).abs();
        let fine_molten_error = (fine.equivalent_molten_length_m - reference_front).abs();
        assert!(
            fine_molten_error < coarse_molten_error,
            "molten-length errors: coarse={coarse_molten_error}, fine={fine_molten_error}"
        );
        assert!(
            fine.temperature_rms_error_k < 0.6 * coarse.temperature_rms_error_k,
            "RMS temperature errors: coarse={}, fine={}",
            coarse.temperature_rms_error_k,
            fine.temperature_rms_error_k,
        );
    });
}
