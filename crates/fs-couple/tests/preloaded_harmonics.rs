//! G1/G3: a held structural preload changes the actual forced spectrum.
//! The source is one retained Duffing coordinate with known physical units;
//! fs-nlmodal owns its law, fs-orbit owns the solve, and fs-couple observes work.

use fs_couple::phs_orbit::spectrum::{PortSpectrumControls, observe_port_spectrum};
use fs_couple::phs_orbit::{
    PeriodicPortInput, audit_forced_orbit, solve_forced_harmonic_balance_adaptive_scaled_seeded,
    solve_forced_harmonic_balance_scaled_seeded,
};
use fs_math::c64::C64;
use fs_nlmodal::{PreloadedModalStorage, SosModalStorage, StressChannel, assemble_storage};
use fs_orbit::{HbBudget, HbRefinementBudget, HbScales};
use fs_phs::PortHamiltonian;

const STIFFNESS: f64 = 4.0;
const CUBIC: f64 = 10.0;
const DAMPING: f64 = 0.3;
const OMEGA: f64 = 1.3;

fn system(bias: f64) -> PortHamiltonian {
    let source = SosModalStorage {
        omegas: vec![2.0],
        channels: vec![StressChannel {
            coefficient: CUBIC,
            coupling: vec![1.0],
        }],
    };
    let storage = PreloadedModalStorage::new(source, &[bias]).unwrap();
    close(
        storage.preload_forces()[0],
        STIFFNESS * bias + CUBIC * bias.powi(3),
        1e-14,
    );
    assemble_storage(
        1,
        &[2.0],
        &[DAMPING / 4.0],
        1,
        vec![0.0, 1.0],
        Box::new(storage),
    )
    .unwrap()
}

fn scales() -> HbScales {
    HbScales {
        states: vec![0.2, 0.4],
        equations: vec![0.4, 2.0],
        angular_frequency: OMEGA,
    }
}

fn phase(angle: f64) -> C64 {
    C64::new(angle.cos(), angle.sin())
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:e} != {expected:e}"
    );
}

#[test]
fn translated_storage_recovers_an_independently_manufactured_multiharmonic_load() {
    let (bias, amplitude, angle) = (0.4, 0.15, 0.37);
    let source = system(bias);
    // q=A cos(theta+angle) has exactly this applied force for
    // q''+C q'+(K+3 beta b^2)q+3 beta b q^2+beta q^3=F.
    // The static holding force K*b+beta*b^3 is already subtracted by the owner.
    let linear = STIFFNESS + 3.0 * CUBIC * bias * bias;
    let quadratic = 3.0 * CUBIC * bias;
    let input = PeriodicPortInput::new(
        OMEGA,
        vec![vec![
            C64::from_re(0.5 * quadratic * amplitude * amplitude),
            phase(angle)
                * C64::new(
                    0.5 * ((linear - OMEGA * OMEGA) * amplitude + 0.75 * CUBIC * amplitude.powi(3)),
                    0.5 * DAMPING * OMEGA * amplitude,
                ),
            phase(2.0 * angle).scale(0.25 * quadratic * amplitude * amplitude),
            phase(3.0 * angle).scale(0.125 * CUBIC * amplitude.powi(3)),
        ]],
    )
    .unwrap();
    let budget = HbBudget {
        harmonics: 5,
        tolerance: 1e-12,
        ..HbBudget::default()
    };
    let solved = solve_forced_harmonic_balance_scaled_seeded(
        &source,
        &input,
        vec![0.0; 2 * (2 * budget.harmonics + 1)],
        &budget,
        &scales(),
    )
    .unwrap()
    .orbit;
    let q = phase(angle).scale(amplitude / 2.0);
    let p = C64::new(0.0, OMEGA) * q;
    for (row, expected) in [(0, q), (1, p)] {
        for (order, actual) in solved.coeffs[row].iter().enumerate() {
            let wanted = if order == 1 { expected } else { C64::ZERO };
            assert!(
                (*actual - wanted).abs() < 1e-10,
                "state {row}, harmonic {order}"
            );
        }
    }
    let audit = audit_forced_orbit(&source, &input, &solved, 256).unwrap();
    let loss = 0.5 * DAMPING * OMEGA * OMEGA * amplitude * amplitude;
    close(audit.mean_supplied_power_w, loss, 1e-11);
    close(audit.mean_resistance_loss_w, loss, 1e-11);
    assert!(audit.relative_power_residual < 1e-9);

    let spectrum = observe_port_spectrum(
        &source,
        &solved,
        &PortSpectrumControls {
            harmonics: 5,
            initial_samples: 32,
            max_samples: 256,
            max_source_evaluations: 480,
            max_workspace_bytes: 1 << 20,
            rms_tolerances: vec![1e-10],
        },
        || false,
    )
    .unwrap();
    // The same velocity has a negative-time PEAK phasor at the output seam.
    let expected = C64::new(0.0, -OMEGA * amplitude) * phase(-angle);
    assert!((spectrum.peak_phasors[0][0] - expected).abs() < 1e-10);
    assert!(
        spectrum.peak_phasors[0][1..]
            .iter()
            .all(|z| z.abs() < 1e-10)
    );
}

#[test]
fn one_physical_drive_generates_even_harmonics_only_when_the_structure_is_biased() {
    let input = PeriodicPortInput::new(OMEGA, vec![vec![C64::ZERO, C64::from_re(0.2)]]).unwrap();
    let budget = HbBudget {
        harmonics: 2,
        tolerance: 1e-12,
        ..HbBudget::default()
    };
    let refinement = HbRefinementBudget {
        max_harmonics: 15,
        residual_tolerances: vec![1e-8, 1e-8],
        state_tolerances: vec![1e-8, 1e-8],
        frequency_tolerance: 1e-12,
    };
    let solve = |bias| {
        let source = system(bias);
        let result = solve_forced_harmonic_balance_adaptive_scaled_seeded(
            &source,
            &input,
            vec![0.0; 10],
            &budget,
            &refinement,
            &scales(),
        )
        .unwrap();
        assert!(result.refinement_trace.len() >= 2);
        let audit = audit_forced_orbit(&source, &input, &result.orbit, 512).unwrap();
        assert!(audit.mean_resistance_loss_w > 0.0);
        assert!(audit.relative_power_residual < 1e-8);
        result.orbit
    };
    let positive = solve(0.4);
    let negative = solve(-0.4);
    let centered = solve(0.0);
    assert!(
        positive.coeffs[0][2].abs() > 1e-4,
        "the even partial must be physical and nonzero"
    );
    assert!(centered.coeffs[0][2].abs() < 1e-10);
    // Reversing the static displacement reverses even harmonics and DC;
    // odd harmonics keep their amplitude and phase under the same drive.
    for state in 0..2 {
        for (order, (&plus, &minus)) in positive.coeffs[state]
            .iter()
            .zip(&negative.coeffs[state])
            .enumerate()
        {
            let expected = plus.scale(if order % 2 == 0 { -1.0 } else { 1.0 });
            assert!(
                (minus - expected).abs() < 1e-9,
                "state {state}, harmonic {order}"
            );
        }
    }
}
