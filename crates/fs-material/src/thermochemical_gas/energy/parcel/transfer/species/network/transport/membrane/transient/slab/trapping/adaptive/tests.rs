use super::super::super::finite::MolarPressureCapacity;
use super::super::super::{
    MembraneDiffusion, MembraneDiffusionLayer, MembraneDiffusivity, MolarHenrySolubility,
    MolarPermeance,
};
use super::super::LangmuirMembraneTrap;
use super::*;

const AREA: f64 = 0.01;
const WIDTH: f64 = 2e-4;
const SOLUBILITY: f64 = 1e-5;
const PRESSURE: f64 = 80_000.0;

fn state(
    cells: usize,
    pressure: f64,
    loading: f64,
    sealed: bool,
) -> TrappingMembraneDiffusionState {
    let seal = sealed.then(|| MolarPermeance::try_new(0.0).unwrap());
    let model = MembraneDiffusion::new_with_surfaces(
        AREA,
        &[MembraneDiffusionLayer {
            thickness_m: WIDTH,
            cells,
            diffusivity: MembraneDiffusivity::try_new(2e-10).unwrap(),
            solubility: MolarHenrySolubility::try_new(SOLUBILITY).unwrap(),
        }],
        &[],
        [seal; 2],
        || false,
    )
    .unwrap();
    TrappingMembraneDiffusionState::new(
        model
            .initial_state(&vec![pressure; cells], || false)
            .unwrap(),
        &[Some(LangmuirMembraneTrap::try_new(2.0, 2e-5, 0.5).unwrap())],
        &vec![loading; cells],
        || false,
    )
    .unwrap()
}
fn controls() -> TrappingMembraneAdvanceControls {
    TrappingMembraneAdvanceControls {
        temporal: MembraneDiffusionAdvanceControls {
            step: MembraneDiffusionSolve {
                absolute_amount_tolerance_mol: 1e-14,
                relative_amount_tolerance: 1e-10,
            },
            initial_step_s: 1.0,
            minimum_step_s: 1e-9,
            maximum_step_s: 10.0,
            maximum_evaluations: 100_000,
            absolute_pressure_tolerance_pa: 5.0,
            absolute_amount_tolerance_mol: 1e-10,
            relative_tolerance: 0.0,
            amount_balance_tolerance_mol: 1e-13,
        },
        maximum_iterations: 64,
        maximum_backtracks: 32,
    }
}
fn near(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a:.16e} != {b:.16e} (tolerance {tolerance:.3e})"
    );
}

#[test]
fn sealed_kinetics_refines_separate_stores_even_when_total_loading_is_constant() {
    let initial = state(1, PRESSURE, 0.2e-6, true);
    // Independent continuous Riccati solution of finite mobile/site exchange.
    let cm = AREA * WIDTH * SOLUBILITY;
    let ct = AREA * WIDTH * 2.0;
    let total = cm * PRESSURE + 0.2e-6;
    let a = 0.5 * 2e-5 / cm;
    let b = a * (total + ct) + 0.5;
    let c = a * total * ct;
    let discriminant = (b * b - 4.0 * a * c).sqrt();
    let low = 2.0 * c / (b + discriminant);
    let high = (b + discriminant) / (2.0 * a);
    let z = ((0.2e-6 - low) / (0.2e-6 - high)) * (-discriminant).exp();
    let exact = (low - z * high) / (1.0 - z);
    let mut errors = Vec::new();
    let mut counts = Vec::new();
    for amount in [1e-8, 1e-11] {
        let mut c = controls();
        c.temporal.absolute_pressure_tolerance_pa = 1e12;
        c.temporal.absolute_amount_tolerance_mol = amount;
        let mut work = initial.adaptive_workspace().unwrap();
        let candidate = initial
            .advance_candidate([0.0; 2], 1.0, c, &mut work, || false)
            .unwrap();
        let report = candidate.report();
        errors.push((candidate.state().trapped_moles()[0] - exact).abs());
        counts.push(report.transport.accepted_intervals);
        assert!(report.transport.rejected_intervals > 0);
        assert_eq!(report.transport.face_amounts_mol, vec![0.0; 2]);
        assert!(report.trapped_change_mol > 0.0);
        assert!(report.mobile_change_mol < 0.0);
        near(
            report.mobile_change_mol + report.trapped_change_mol,
            0.0,
            1e-18,
        );
        near(report.transport.stored_before_mol, total, 1e-20);
        near(report.transport.stored_after_mol, total, 1e-18);
        near(candidate.state().stored_amount_mol().unwrap(), total, 1e-18);
        assert!(report.transport.maximum_accepted_error_ratio <= 1.0);
    }
    assert!(counts[1] > counts[0]);
    assert!(
        errors[1] < 0.1 * errors[0],
        "kinetic refinement: {errors:?}"
    );
    assert!(errors[1] < 2e-9, "kinetic error: {errors:?}");
}

#[test]
fn finite_reservoir_adaptation_preserves_traps_and_accumulated_transport() {
    let capacity = MolarPressureCapacity::try_new(AREA * WIDTH * SOLUBILITY).unwrap();
    let initial = FiniteTrappingMembraneDiffusionState::new(
        state(4, 0.0, 0.0, false),
        [capacity; 2],
        [PRESSURE, 0.0],
    )
    .unwrap();
    let total = initial.total_amount_mol().unwrap();
    let mut reference = initial.clone();
    let mut reference_work = reference.workspace().unwrap();
    let mut transmitted = 0.0;
    for _ in 0..2048 {
        let next = reference
            .step_candidate(
                40.0 / 2048.0,
                TrappingMembraneSolve {
                    balance: MembraneDiffusionSolve {
                        absolute_amount_tolerance_mol: 1e-19,
                        relative_amount_tolerance: 1e-11,
                    },
                    maximum_iterations: 64,
                    maximum_backtracks: 32,
                },
                &mut reference_work,
                || false,
            )
            .unwrap();
        transmitted += next.step().boundary_transfers_mol[1];
        reference = next.into_parts().0;
    }
    let mut work = initial.adaptive_workspace().unwrap();
    let candidate = initial
        .advance_candidate(40.0, controls(), &mut work, || false)
        .unwrap();
    let report = candidate.report();
    assert_eq!(candidate.state().membrane().mobile().elapsed_s(), 40.0);
    assert!(report.trapped_after_mol > 0.0);
    near(candidate.state().total_amount_mol().unwrap(), total, 1e-16);
    let finite = report.transport.finite_reservoirs.as_ref().unwrap();
    assert!(
        finite.closed_amount_residual_mol.abs() <= controls().temporal.amount_balance_tolerance_mol
    );
    near(
        finite.amount_changes_mol[0]
            + report.mobile_change_mol
            + report.trapped_change_mol
            + finite.amount_changes_mol[1],
        0.0,
        1e-16,
    );
    for (&a, &b) in candidate
        .state()
        .membrane()
        .mobile()
        .partial_pressures_pa()
        .iter()
        .zip(reference.membrane().mobile().partial_pressures_pa())
    {
        near(a, b, 100.0);
    }
    for (&a, &b) in candidate
        .state()
        .membrane()
        .trapped_moles()
        .iter()
        .zip(reference.membrane().trapped_moles())
    {
        near(a, b, 2e-9);
    }
    near(report.transport.face_amounts_mol[4], transmitted, 1e-9);
    for (&signed, &gross) in report
        .transport
        .face_amounts_mol
        .iter()
        .zip(&report.transport.gross_face_amounts_mol)
    {
        assert!(gross >= signed.abs());
    }
}

#[test]
fn work_native_refusal_and_late_cancellation_preserve_both_inventories_and_clock() {
    let initial = state(4, 0.0, 0.0, false);
    let saved = initial.clone();
    let mut work = initial.adaptive_workspace().unwrap();
    let mut limited = controls();
    limited.temporal.initial_step_s = 0.25;
    limited.temporal.maximum_step_s = 0.25;
    limited.temporal.maximum_evaluations = 3;
    limited.temporal.absolute_pressure_tolerance_pa = 1e12;
    limited.temporal.absolute_amount_tolerance_mol = 1.0;
    assert!(
        matches!(initial.advance_candidate([PRESSURE, 0.0], 1.0, limited, &mut work, || false),
        Err(MembraneDiffusionAdvanceError::WorkLimit { time_s, .. }) if time_s > 0.0)
    );
    assert_eq!(initial, saved);
    let mut native_limit = controls();
    native_limit.maximum_iterations = 1;
    assert!(matches!(
        initial.advance_candidate([PRESSURE, 0.0], 1.0, native_limit, &mut work, || false),
        Err(MembraneDiffusionAdvanceError::Native(
            TrappingMembraneError::NotConverged { .. }
        ))
    ));
    assert_eq!(initial, saved);
    let mut polls = 0;
    let expected = initial
        .advance_candidate([PRESSURE, 0.0], 1.0, controls(), &mut work, || {
            polls += 1;
            false
        })
        .unwrap();
    let mut calls = 0;
    assert!(matches!(
        initial.advance_candidate([PRESSURE, 0.0], 1.0, controls(), &mut work, || {
            calls += 1;
            calls == polls
        }),
        Err(MembraneDiffusionAdvanceError::Cancelled)
    ));
    assert_eq!(calls, polls);
    assert_eq!(initial, saved);
    assert_eq!(
        expected,
        initial
            .advance_candidate([PRESSURE, 0.0], 1.0, controls(), &mut work, || false)
            .unwrap()
    );
}

#[test]
fn equilibrium_traps_reach_the_exact_rounded_target_with_both_stores_reported() {
    let mut initial = state(1, PRESSURE, 0.0, true);
    let bq = 2e-5 * PRESSURE;
    let capacity = AREA * WIDTH * 2.0;
    initial.trapped_moles[0] = capacity * bq / (1.0 + bq);
    let mut c = controls();
    c.temporal.initial_step_s = 0.01;
    c.temporal.maximum_step_s = 0.01;
    let mut work = initial.adaptive_workspace().unwrap();
    let next = initial
        .advance_candidate([0.0; 2], 0.1, c, &mut work, || false)
        .unwrap();
    assert_eq!(
        next.state().mobile().elapsed_s().to_bits(),
        0.1_f64.to_bits()
    );
    assert_eq!(next.report().transport.accepted_intervals, 10);
    assert_eq!(next.report().transport.rejected_intervals, 0);
    near(
        next.report().trapped_after_mol,
        initial.trapped_moles()[0],
        1e-20,
    );
    near(
        next.report().transport.stored_before_mol,
        initial.stored_amount_mol().unwrap(),
        1e-20,
    );
    near(
        next.report().transport.stored_after_mol,
        initial.stored_amount_mol().unwrap(),
        1e-19,
    );
}
