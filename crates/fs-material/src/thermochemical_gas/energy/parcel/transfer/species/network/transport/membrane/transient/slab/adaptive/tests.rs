//! G1 Fourier/eigenmode limits, G3 storage sensitivity, and G4 atomic refusal.
use super::super::finite::MolarPressureCapacity;
use super::super::{
    MembraneDiffusion, MembraneDiffusionLayer, MembraneDiffusivity, MolarHenrySolubility,
    MolarPermeance,
};
use super::*;
use fs_math::det;

const AREA: f64 = 0.01;
const L: f64 = 2e-4;
const D: f64 = 2e-10;
const S: f64 = 1e-5;
const P: f64 = 80000.0;
const TAU: f64 = L * L / D;

fn controls() -> MembraneDiffusionAdvanceControls {
    MembraneDiffusionAdvanceControls {
        step: MembraneDiffusionSolve {
            absolute_amount_tolerance_mol: 1e-14,
            relative_amount_tolerance: 1e-10,
        },
        initial_step_s: 0.5 * TAU,
        minimum_step_s: 1e-8,
        maximum_step_s: 0.5 * TAU,
        maximum_evaluations: 100_000,
        absolute_pressure_tolerance_pa: 20.0,
        absolute_amount_tolerance_mol: 1e-9,
        relative_tolerance: 0.0,
        amount_balance_tolerance_mol: 1e-12,
    }
}
fn layer(cells: usize, solubility: f64) -> MembraneDiffusionLayer {
    MembraneDiffusionLayer {
        thickness_m: L,
        cells,
        diffusivity: MembraneDiffusivity::try_new(D).unwrap(),
        solubility: MolarHenrySolubility::try_new(solubility).unwrap(),
    }
}
fn slab(cells: usize) -> MembraneDiffusionState {
    MembraneDiffusion::new(AREA, &[layer(cells, S)], &[], || false)
        .unwrap()
        .initial_state(&vec![0.0; cells], || false)
        .unwrap()
}
fn symmetric() -> FiniteMembraneDiffusionState {
    let capacity = MolarPressureCapacity::try_new(AREA * S * L).unwrap();
    FiniteMembraneDiffusionState::new(slab(1), [capacity; 2], [P, 0.0], || false).unwrap()
}
fn outlet_fraction(time: f64) -> f64 {
    let mut result = 1.0;
    for n in 1..=100 {
        let sign = if n % 2 == 0 { 1.0 } else { -1.0 };
        result +=
            2.0 * sign * det::exp(-(n * n) as f64 * core::f64::consts::PI.powi(2) * time / TAU);
    }
    result
}
fn outlet_amount(time: f64) -> f64 {
    let mut integral = time / TAU - 1.0 / 6.0;
    for n in 1..=100 {
        let sign = if n % 2 == 0 { 1.0 } else { -1.0 };
        integral -= 2.0 * sign / ((n * n) as f64 * core::f64::consts::PI.powi(2))
            * det::exp(-(n * n) as f64 * core::f64::consts::PI.powi(2) * time / TAU);
    }
    AREA * S * P * L * integral
}

#[test]
fn controlled_breakthrough_and_cumulative_permeation_converge_to_fourier_solution() {
    let state = slab(96);
    let time = 0.2 * TAU;
    let mut flux_errors = Vec::new();
    let mut amount_errors = Vec::new();
    for pressure in [200.0, 5.0] {
        let c = MembraneDiffusionAdvanceControls {
            absolute_pressure_tolerance_pa: pressure,
            absolute_amount_tolerance_mol: pressure * 1e-11,
            ..controls()
        };
        let mut work = state.adaptive_workspace().unwrap();
        let candidate = state
            .advance_candidate([P, 0.0], time, c, &mut work, || false)
            .unwrap();
        let report = candidate.report();
        assert_eq!(candidate.state().elapsed_s(), time);
        assert!(report.rejected_intervals > 0);
        assert_eq!(
            report.evaluations,
            3 * (report.accepted_intervals + report.rejected_intervals)
        );
        assert_eq!(report.accepted_substeps, 2 * report.accepted_intervals);
        assert!(report.maximum_accepted_error_ratio <= 1.0);
        assert!(report.amount_residual_mol.abs() <= c.amount_balance_tolerance_mol);
        assert!(report.maximum_cell_amount_residual_mol <= c.amount_balance_tolerance_mol);
        for (&signed, &gross) in report
            .face_amounts_mol
            .iter()
            .zip(&report.gross_face_amounts_mol)
        {
            assert!(gross >= signed.abs());
        }
        let normalized = report.last_substep.face_molar_fluxes_mol_m2_s[96] / (D * S * P / L);
        flux_errors.push((normalized - outlet_fraction(time)).abs());
        amount_errors.push((report.face_amounts_mol[96] - outlet_amount(time)).abs());
    }
    assert!(flux_errors[1] < 0.4 * flux_errors[0], "{flux_errors:?}");
    assert!(
        amount_errors[1] < 0.4 * amount_errors[0],
        "{amount_errors:?}"
    );
    assert!(flux_errors[1] < 0.003, "{flux_errors:?}");
}

#[test]
fn finite_reservoir_accuracy_uses_the_continuous_three_inventory_eigenmodes() {
    let state = symmetric();
    let time = 0.4 * TAU;
    let lambda = 2.0 / TAU;
    let slow = det::exp(-lambda * time);
    let fast = det::exp(-3.0 * lambda * time);
    let exact = [
        P * (1.0 / 3.0 + 0.5 * slow + fast / 6.0),
        P * (1.0 - fast) / 3.0,
        P * (1.0 / 3.0 - 0.5 * slow + fast / 6.0),
    ];
    let mut errors = Vec::new();
    for pressure in [100.0, 1.0] {
        let c = MembraneDiffusionAdvanceControls {
            absolute_pressure_tolerance_pa: pressure,
            absolute_amount_tolerance_mol: 1.0,
            ..controls()
        };
        let mut work = state.adaptive_workspace().unwrap();
        let candidate = state
            .advance_candidate(time, c, &mut work, || false)
            .unwrap();
        let p = candidate.state().reservoir_pressures_pa();
        let actual = [
            p[0],
            candidate.state().membrane().partial_pressures_pa()[0],
            p[1],
        ];
        errors.push(
            actual
                .iter()
                .zip(exact)
                .map(|(a, b)| (a - b).abs())
                .fold(0.0_f64, f64::max),
        );
        assert!(
            candidate
                .report()
                .finite_reservoirs
                .as_ref()
                .unwrap()
                .closed_amount_residual_mol
                .abs()
                <= c.amount_balance_tolerance_mol
        );
        assert!(
            (candidate.state().total_amount_mol().unwrap() - state.total_amount_mol().unwrap())
                .abs()
                < 1e-12
        );
    }
    assert!(errors[1] < 0.2 * errors[0], "{errors:?}");
    assert!(errors[1] < 100.0, "{errors:?}");
}

#[test]
fn amount_accuracy_remains_active_when_pressure_accuracy_is_loose() {
    let state = symmetric();
    let mut loose = controls();
    loose.absolute_pressure_tolerance_pa = 1e12;
    loose.absolute_amount_tolerance_mol = 1.0;
    let mut work = state.adaptive_workspace().unwrap();
    let coarse = state
        .advance_candidate(0.25 * TAU, loose, &mut work, || false)
        .unwrap();
    assert_eq!(coarse.report().evaluations, 3);
    let tight = MembraneDiffusionAdvanceControls {
        absolute_amount_tolerance_mol: 1e-10,
        ..loose
    };
    let fine = state
        .advance_candidate(0.25 * TAU, tight, &mut work, || false)
        .unwrap();
    assert!(fine.report().rejected_intervals > 0);
    assert!(fine.report().accepted_intervals > 10);
    assert!(fine.report().maximum_accepted_error_ratio <= 1.0);
}

#[test]
fn algebraic_accuracy_cannot_be_loosened_and_native_refusals_are_counted_once() {
    let state = slab(8);
    let saved = state.clone();
    let mut work = state.adaptive_workspace().unwrap();
    let reference = state
        .advance_candidate([P, 0.0], 0.05 * TAU, controls(), &mut work, || false)
        .unwrap();
    let mut c = controls();
    c.step = MembraneDiffusionSolve {
        absolute_amount_tolerance_mol: 1.0,
        relative_amount_tolerance: 0.5,
    };
    assert_eq!(
        reference,
        state
            .advance_candidate([P, 0.0], 0.05 * TAU, c, &mut work, || false)
            .unwrap()
    );

    // An impossible temporal request fails the additional algebraic audit on
    // its first probe even though the authored native residual controls are lax.
    c.absolute_amount_tolerance_mol = 1e-40;
    let mut calls = 0;
    assert!(matches!(
        probe(
            &state,
            [P, 0.0],
            1.0,
            1.0,
            c,
            &mut work.full,
            &mut calls,
            &mut || false
        ),
        Err(Error::AmountBalance { .. })
    ));
    assert_eq!(calls, 1);

    // Original native admission is preserved before the stricter adaptive audit.
    c = controls();
    c.step.absolute_amount_tolerance_mol = 1e-40;
    c.step.relative_amount_tolerance = 0.0;
    calls = 0;
    assert!(matches!(
        probe(
            &state,
            [P, 0.0],
            1.0,
            1.0,
            c,
            &mut work.full,
            &mut calls,
            &mut || false
        ),
        Err(Error::Native(MembraneDiffusionError::AmountBalance { .. }))
    ));
    assert_eq!(calls, 1);
    assert_eq!(state, saved);
}

#[test]
fn sealed_layer_equilibria_preserve_exact_storage_and_complete_the_requested_clock() {
    let seal = Some(MolarPermeance::try_new(0.0).unwrap());
    let state = MembraneDiffusion::new_with_surfaces(
        AREA,
        &[layer(3, S), layer(5, 1e3 * S)],
        &[seal],
        [seal; 2],
        || false,
    )
    .unwrap()
    .initial_state(&[1234.0, 1234.0, 1234.0, P, P, P, P, P], || false)
    .unwrap();
    let mut work = state.adaptive_workspace().unwrap();
    let candidate = state
        .advance_candidate([1e6, 0.0], 0.1, controls(), &mut work, || false)
        .unwrap();
    assert_eq!(
        candidate.state().partial_pressures_pa(),
        state.partial_pressures_pa()
    );
    assert_eq!(candidate.report().face_amounts_mol, vec![0.0; 9]);
    assert_eq!(candidate.report().stored_change_mol, 0.0);
    assert_eq!(candidate.report().interval_s, [0.0, 0.1]);
    assert_eq!(candidate.report().evaluations, 3);
}

#[test]
fn boundary_history_is_split_explicitly_and_unloading_retains_both_signed_streams() {
    let state = slab(16);
    let mut work = state.adaptive_workspace().unwrap();
    let c = MembraneDiffusionAdvanceControls {
        absolute_pressure_tolerance_pa: 100.0,
        ..controls()
    };
    let loaded = state
        .advance_candidate([P; 2], 0.1 * TAU, c, &mut work, || false)
        .unwrap();
    let unloaded = loaded
        .state()
        .advance_candidate([0.0; 2], 0.1 * TAU, c, &mut work, || false)
        .unwrap();
    assert_eq!(unloaded.report().interval_s, [0.1 * TAU, 0.2 * TAU]);
    assert!(unloaded.report().face_amounts_mol[0] < 0.0);
    assert!(unloaded.report().face_amounts_mol[16] > 0.0);
    assert!(unloaded.report().stored_change_mol < 0.0);
    assert!(unloaded.report().amount_residual_mol.abs() < c.amount_balance_tolerance_mol);
}

#[test]
fn work_accuracy_and_late_cancellation_publish_no_partial_state_and_caches_reuse() {
    let state = slab(8);
    let saved = state.clone();
    let mut work = state.adaptive_workspace().unwrap();
    let mut c = controls();
    c.maximum_evaluations = 3;
    assert!(matches!(
        state.advance_candidate([P, 0.0], TAU, c, &mut work, || false),
        Err(Error::WorkLimit { .. })
    ));
    c = controls();
    c.minimum_step_s = c.initial_step_s;
    assert!(matches!(
        state.advance_candidate([P, 0.0], TAU, c, &mut work, || false),
        Err(Error::Accuracy { .. })
    ));
    assert_eq!(state, saved);
    let mut calls = 0;
    let result = state
        .advance_candidate([P, 0.0], 0.05 * TAU, controls(), &mut work, || {
            calls += 1;
            false
        })
        .unwrap();
    for stop in [1, calls / 2, calls] {
        let mut seen = 0;
        assert!(matches!(
            state.advance_candidate([P, 0.0], 0.05 * TAU, controls(), &mut work, || {
                seen += 1;
                seen == stop
            }),
            Err(Error::Cancelled)
        ));
        assert_eq!(state, saved);
    }
    let mut fresh = state.adaptive_workspace().unwrap();
    assert_eq!(
        result,
        state
            .advance_candidate([P, 0.0], 0.05 * TAU, controls(), &mut fresh, || false)
            .unwrap()
    );
    assert_eq!(
        result,
        state
            .advance_candidate([P, 0.0], 0.05 * TAU, controls(), &mut work, || false)
            .unwrap()
    );
    c = controls();
    c.amount_balance_tolerance_mol = 1e-40;
    assert!(matches!(
        state.advance_candidate([P, 0.0], 0.05 * TAU, c, &mut work, || false),
        Err(Error::AmountBalance { .. })
    ));
    assert_eq!(state, saved);
}
