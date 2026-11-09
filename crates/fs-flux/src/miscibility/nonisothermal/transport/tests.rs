//! Complete-interval limits, conservation and an independent linear mode.

use super::*;
use crate::miscibility::nonisothermal::advection::tests::{enthalpy, fixture, near, opposing_flow};

fn transport(
    initial: &BinaryNonisothermalField,
    flow: &BinaryMiscibilityAdvection,
    duration: f64,
) -> BinaryNonisothermalTransportCandidate {
    initial
        .transport_candidate(
            flow,
            duration,
            BinaryMiscibilityAdvectionControls::default(),
            BinaryNonisothermalSolve::default(),
            || false,
        )
        .unwrap()
}

#[test]
fn g0_zero_flow_is_original_local_heat_and_species_diffusion() {
    let initial = fixture(
        [20.0, 60.0],
        0.0,
        1e-5,
        0.1,
        &[[0.001, 0.3, 340.0], [0.001, 0.7, 520.0]],
        &[[0, 1]],
    );
    let flow = opposing_flow(&initial, 0.0);
    let expected = initial
        .step_candidate(0.01, BinaryNonisothermalSolve::default(), || false)
        .unwrap();
    let candidate = transport(&initial, &flow, 0.01);
    assert_eq!(candidate.field(), expected.field());
    assert_eq!(&candidate.step().diffusion, expected.step());
    assert_eq!(
        candidate.step().advection.transferred_component_moles,
        [[0.0; 2]]
    );
    assert_eq!(candidate.step().advection.transferred_energy_j, [0.0]);
    assert_eq!(candidate.field().elapsed_s(), 0.01);
}

#[test]
fn g0_zero_molecular_coefficients_retain_advection_and_one_physical_clock() {
    let initial = fixture(
        [20.0, 60.0],
        -100.0,
        0.0,
        0.0,
        &[[0.001, 0.3, 340.0], [0.001, 0.7, 520.0]],
        &[[0, 1]],
    );
    let flow = opposing_flow(&initial, 0.0004);
    let expected = initial
        .advect_candidate(
            &flow,
            0.5,
            BinaryMiscibilityAdvectionControls::default(),
            BinaryNonisothermalSolve::default(),
            || false,
        )
        .unwrap();
    let candidate = transport(&initial, &flow, 0.5);
    assert_eq!(candidate.field(), expected.field());
    assert_eq!(&candidate.step().advection, expected.step());
    assert_eq!(candidate.step().diffusion.transferred_a_mol, [0.0]);
    assert_eq!(candidate.step().diffusion.transferred_energy_j, [0.0]);
    assert_eq!(
        transport(candidate.field(), &flow, 0.5).field().elapsed_s(),
        1.0
    );
}

#[test]
fn g0_complete_receipts_close_actual_species_and_nonideal_parcel_energy() {
    let cp = [20.0, 60.0];
    let b = -100.0;
    let initial = fixture(
        cp,
        b,
        1e-5,
        0.03,
        &[[0.001, 0.3, 340.0], [0.002, 0.7, 520.0]],
        &[[0, 1]],
    );
    let flow = opposing_flow(&initial, 0.0003);
    let candidate = transport(&initial, &flow, 0.03);
    let report = candidate.step();
    for (cell, sign) in [(0, -1.0), (1, 1.0)] {
        for s in 0..2 {
            near(
                candidate.field().component_moles()[cell][s] - initial.component_moles()[cell][s],
                sign * report.transferred_component_moles[0][s],
                2e-13,
            );
        }
        let delta_h = enthalpy(
            candidate.field().component_moles()[cell],
            candidate.field().temperatures()[cell].value(),
            cp,
            b,
        ) - enthalpy(
            initial.component_moles()[cell],
            initial.temperatures()[cell].value(),
            cp,
            b,
        );
        near(delta_h, sign * report.transferred_energy_j[0], 1e-8);
    }
    for (s, sign) in [(0, 1.0), (1, -1.0)] {
        near(
            report.transferred_component_moles[0][s],
            report.advection.transferred_component_moles[0][s]
                + sign * report.diffusion.transferred_a_mol[0],
            2e-14,
        );
    }
    near(
        report.transferred_energy_j[0],
        report.advection.transferred_energy_j[0] + report.diffusion.transferred_energy_j[0],
        1e-10,
    );
    near(report.energy_residual_j, 0.0, 1e-8);
    assert!(report.entropy_change_j_per_k > 0.0);
    assert_eq!(candidate.field().elapsed_s(), 0.03);
}

#[test]
fn g1_macro_time_refinement_converges_to_independent_ideal_mixing_mode() {
    let amplitude = 1e-4;
    let initial = fixture(
        [40.0; 2],
        0.0,
        1e-4,
        0.0,
        &[
            [0.001, 0.5 + amplitude, 400.0],
            [0.001, 0.5 - amplitude, 400.0],
        ],
        &[[0, 1]],
    );
    let flow = opposing_flow(&initial, 0.0002);
    // N=1 mol, c=1000 mol/m³, A/d=1 m. The linear ideal mode obeys
    // a'=-2(Q/V+c*D/N)a=-0.6a. This checks the complete split operator.
    let expected = 0.5 + amplitude * (-0.6_f64).exp();
    let mut errors = Vec::new();
    for steps in [4, 16, 64] {
        let mut field = initial.clone();
        for _ in 0..steps {
            field = transport(&field, &flow, 1.0 / f64::from(steps))
                .into_parts()
                .0;
        }
        errors.push((field.mole_fractions()[0][0] - expected).abs());
        for t in field.temperatures() {
            near(t.value(), 400.0, 1e-8);
        }
        near(field.elapsed_s(), 1.0, 2e-14);
    }
    assert!(errors[1] < errors[0] / 3.0, "{errors:?}");
    assert!(errors[2] < errors[1] / 3.0, "{errors:?}");
    assert!(errors[2] < 1e-7, "{errors:?}");
}
