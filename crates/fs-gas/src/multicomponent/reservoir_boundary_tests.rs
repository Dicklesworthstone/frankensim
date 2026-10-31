//! G0/G1/G3/G4 prescribed-mixture pressure boundaries on the actual RK owner.

use super::super::characteristic_tests::caloric_gas;
use super::*;

fn flow(gas: MixtureGas, boundaries: Boundaries) -> MixtureEuler1d {
    let initial = gas.conserved(&[1.0, 0.0], 0.0, 1.0).unwrap();
    MixtureEuler1d::new(
        gas,
        Length::new(0.05),
        0.7,
        vec![initial; 12],
        boundaries,
        || false,
    )
    .unwrap()
}

fn walls() -> Boundaries {
    Boundaries::Ends {
        left: Boundary::Reflecting,
        right: Boundary::Reflecting,
    }
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * expected.abs().max(1.0),
        "actual={actual:.16e}, expected={expected:.16e}"
    );
}

fn balance(report: &StepReport) {
    for value in report.conservation_residual.species_masses[..report.before.species_count]
        .iter()
        .chain([
            &report.conservation_residual.momentum,
            &report.conservation_residual.energy,
        ])
    {
        close(*value, 0.0, 3e-12);
    }
}

#[test]
fn absent_reservoirs_preserve_original_field_report_clock_and_cancellation_sequence() {
    for boundaries in [walls(), Boundaries::Periodic] {
        let gas = caloric_gas(true, 0.0);
        let states: Vec<_> = (0..12)
            .map(|i| {
                let density = 1.0 + 0.08 * (i as f64 * 0.7).sin();
                gas.conserved(&[0.4 * density, 0.6 * density], 0.13, 1.0)
                    .unwrap()
            })
            .collect();
        let mut original = MixtureEuler1d::new(
            gas.clone(),
            Length::new(0.05),
            0.7,
            states.clone(),
            boundaries,
            || false,
        )
        .unwrap();
        let mut explicit =
            MixtureEuler1d::new(gas, Length::new(0.05), 0.7, states, boundaries, || false).unwrap();
        for _ in 0..5 {
            let (mut original_polls, mut explicit_polls) = (0, 0);
            let a = original
                .step(Time::new(0.001), StepControl::default(), || {
                    original_polls += 1;
                    false
                })
                .unwrap();
            let b = explicit
                .step_with_reservoirs(Time::new(0.001), StepControl::default(), [None; 2], || {
                    explicit_polls += 1;
                    false
                })
                .unwrap();
            assert_eq!(original_polls, explicit_polls);
            assert_eq!(original.states, explicit.states);
            assert_eq!(original.time, explicit.time);
            assert_eq!(a, b);
        }
    }
}

#[test]
fn pressure_reservoir_refreshes_both_rk_stages_and_transports_actual_donor_enthalpy() {
    let gas = caloric_gas(true, 0.0);
    let initial = gas.conserved(&[1.0, 0.0], 0.0, 1.0).unwrap();
    let mut injected = Vec::new();
    for densities in [[0.0, 2.0], [0.0, 1.0], [1.0, 0.0]] {
        let reservoir = gas.conserved(&densities, 0.0, 2.0).unwrap();
        let mut actual = flow(gas.clone(), walls());
        let report = actual
            .step_with_reservoirs(
                Time::new(0.001),
                StepControl::default(),
                [Some(reservoir.as_cell()), None],
                || false,
            )
            .unwrap();
        balance(&report);
        assert_eq!(report.halvings, 0);
        assert!(report.left_boundary_change.mass() > 0.0);
        assert_eq!(report.right_boundary_change.mass(), 0.0);
        assert_eq!(report.right_boundary_change.energy, 0.0);
        // The first stage's only nonuniform incoming face is the reservoir.
        // Independently form its predictor average, then evaluate the original
        // characteristic owner on that NEW state. Reusing the initial flux for
        // the corrector would fail this check.
        let mut first = [0.0; 4];
        gas.reservoir_flux(
            initial.as_cell(),
            reservoir.as_cell(),
            -1.0,
            &mut first,
            || false,
        )
        .unwrap();
        let factor = report.elapsed.value() / 0.05;
        let predictor = Conserved {
            partial_densities: vec![
                initial.partial_densities[0] + factor * first[0],
                initial.partial_densities[1] + factor * first[1],
            ],
            momentum: factor * (first[2] - 1.0),
            energy: initial.energy + factor * first[3],
        };
        let mut second = [0.0; 4];
        gas.reservoir_flux(
            predictor.as_cell(),
            reservoir.as_cell(),
            -1.0,
            &mut second,
            || false,
        )
        .unwrap();
        assert!((second[0] + second[1] - first[0] - first[1]).abs() > 1e-7);
        let scale = 0.5 * report.elapsed.value() * 0.7;
        for k in 0..2 {
            close(
                report.left_boundary_change.species_masses[k],
                scale * (first[k] + second[k]),
                3e-12,
            );
        }
        close(
            report.left_boundary_change.energy,
            scale * (first[3] + second[3]),
            3e-12,
        );
        let pure_b = densities[0] == 0.0;
        let (r, cp_constant, cp_linear) = if pure_b {
            (0.5, 2.5, 0.35)
        } else {
            (1.0, 3.5, 0.25)
        };
        let temperature = 2.0 / (r * (densities[0] + densities[1]));
        // Independent NASA fixture h(T)=a*T+b*T²/2, with each species'
        // retained e(Tmin) reference removed once from the advected enthalpy.
        let expected_enthalpy = cp_constant * temperature
            + 0.5 * cp_linear * temperature * temperature
            - (cp_constant - r) * 0.001
            - 0.5 * cp_linear * 0.001_f64.powi(2);
        close(
            report.left_boundary_change.energy / report.left_boundary_change.mass(),
            expected_enthalpy,
            2e-12,
        );
        let absent = usize::from(!pure_b);
        assert_eq!(report.left_boundary_change.species_masses[absent], 0.0);
        if pure_b {
            assert!(actual.cell(0).unwrap().partial_densities[1] > 0.0);
            injected.push((
                report.left_boundary_change.mass(),
                report.left_boundary_change.energy,
            ));
        } else {
            assert_eq!(actual.totals().species_masses[1], 0.0);
        }
    }
    assert!((injected[0].0 - injected[1].0).abs() > 1e-6);
    assert!((injected[0].1 - injected[1].1).abs() > 1e-6);
}

#[test]
fn pressure_inflow_reverses_end_orientation_without_changing_species_or_heat() {
    let gas = caloric_gas(true, 0.0);
    let reservoir = gas.conserved(&[0.0, 2.0], 0.0, 2.0).unwrap();
    let mut left = flow(gas.clone(), walls());
    let mut right = flow(gas, walls());
    for _ in 0..6 {
        let a = left
            .step_with_reservoirs(
                Time::new(0.001),
                StepControl::default(),
                [Some(reservoir.as_cell()), None],
                || false,
            )
            .unwrap();
        let b = right
            .step_with_reservoirs(
                Time::new(0.001),
                StepControl::default(),
                [None, Some(reservoir.as_cell())],
                || false,
            )
            .unwrap();
        balance(&a);
        balance(&b);
        close(
            a.left_boundary_change.mass(),
            b.right_boundary_change.mass(),
            3e-12,
        );
        close(
            a.left_boundary_change.energy,
            b.right_boundary_change.energy,
            3e-12,
        );
        assert!(a.left_boundary_change.momentum > 0.0);
        assert!(b.right_boundary_change.momentum < 0.0);
        for i in 0..left.cell_count() {
            let a = left.cell(i).unwrap();
            let b = right.cell(right.cell_count() - 1 - i).unwrap();
            for (x, y) in a.partial_densities.iter().zip(b.partial_densities) {
                close(*x, *y, 3e-12);
            }
            close(a.energy, b.energy, 3e-12);
            close(a.momentum, -b.momentum, 3e-12);
        }
    }
}

#[test]
fn prescribed_mixture_pressure_is_invariant_to_each_species_formation_reference() {
    let mut expected: Option<(Vec<f64>, Vec<StepReport>)> = None;
    for formation in [0.0, -1e18, 1e18] {
        let gas = caloric_gas(true, formation);
        let supply = gas.conserved(&[0.2, 1.8], 0.0, 2.0).unwrap();
        let sink = gas.conserved(&[0.7, 0.3], 0.0, 0.5).unwrap();
        let mut actual = flow(gas, walls());
        let mut reports = Vec::new();
        for _ in 0..4 {
            let report = actual
                .step_with_reservoirs(
                    Time::new(0.001),
                    StepControl::default(),
                    [Some(supply.as_cell()), Some(sink.as_cell())],
                    || false,
                )
                .unwrap();
            balance(&report);
            assert!(report.left_boundary_change.mass() > 0.0);
            assert!(report.right_boundary_change.mass() < 0.0);
            reports.push(report);
        }
        if let Some((states, previous)) = &expected {
            assert_eq!(&actual.states, states);
            assert_eq!(&reports, previous);
        } else {
            expected = Some((actual.states, reports));
        }
    }
}

#[test]
fn reservoir_face_cfl_retries_recompute_fluxes_from_the_accepted_start_state() {
    let gas = caloric_gas(false, 0.0);
    let reservoir = gas.conserved(&[0.0, 100.0], 0.0, 100.0).unwrap();
    let mut actual = flow(gas.clone(), walls());
    let initial = actual.states.clone();
    assert!(matches!(
        actual.step_with_reservoirs(
            Time::new(1.0),
            StepControl {
                maximum_halvings: 0,
                ..StepControl::default()
            },
            [Some(reservoir.as_cell()), None],
            || false
        ),
        Err(Error::StepRetries { .. })
    ));
    assert_eq!(actual.states, initial);
    assert_eq!(actual.time, 0.0);
    let accepted = actual
        .step_with_reservoirs(
            Time::new(1.0),
            StepControl::default(),
            [Some(reservoir.as_cell()), None],
            || false,
        )
        .unwrap();
    assert!(accepted.halvings > 0);
    assert!(accepted.maximum_cfl <= 0.4 * (1.0 + 64.0 * f64::EPSILON));
    balance(&accepted);
    let mut direct = flow(gas, walls());
    let replay = direct
        .step_with_reservoirs(
            accepted.elapsed,
            StepControl {
                reconstruction: if accepted.first_order_retry {
                    Reconstruction::FirstOrder
                } else {
                    Reconstruction::Muscl
                },
                ..StepControl::default()
            },
            [Some(reservoir.as_cell()), None],
            || false,
        )
        .unwrap();
    assert_eq!(replay.halvings, 0);
    // Temperature caches may retain a rejected-stage Newton hint. A retry
    // must reproduce the accepted conservative evolution within the original
    // caloric inversion tolerance, without demanding identical hints.
    for (retried, replayed) in actual.states.iter().zip(&direct.states) {
        close(*retried, *replayed, 3e-12);
    }
    assert_eq!(actual.time, direct.time);
    for k in 0..4 {
        close(
            accepted.boundary_change.get(k),
            replay.boundary_change.get(k),
            3e-12,
        );
    }
}

#[test]
fn invalid_reservoir_and_early_late_cancellation_publish_no_field_or_clock() {
    let gas = caloric_gas(true, 0.0);
    let reservoir = gas.conserved(&[0.0, 2.0], 0.0, 2.0).unwrap();
    let mut actual = flow(gas.clone(), walls());
    let before = actual.states.clone();
    let mut wrong_dimension = reservoir.clone();
    wrong_dimension.partial_densities.pop();
    let mut moving = reservoir.clone();
    moving.momentum = 0.01;
    let mut negative = reservoir.clone();
    negative.partial_densities[0] = -1.0;
    let mut nonfinite = reservoir.clone();
    nonfinite.energy = f64::NAN;
    let mut outside_source = reservoir.clone();
    outside_source.energy = 1e20;
    for invalid in [wrong_dimension, moving, negative, nonfinite, outside_source] {
        assert!(
            actual
                .step_with_reservoirs(
                    Time::new(0.001),
                    StepControl::default(),
                    [Some(invalid.as_cell()), None],
                    || false
                )
                .is_err()
        );
        assert_eq!(actual.states, before);
        assert_eq!(actual.time, 0.0);
    }
    let mut periodic = flow(gas.clone(), Boundaries::Periodic);
    assert!(
        periodic
            .step_with_reservoirs(
                Time::new(0.001),
                StepControl::default(),
                [Some(reservoir.as_cell()), None],
                || false
            )
            .is_err()
    );
    assert_eq!(periodic.states, before);
    let mut admitted = flow(gas.clone(), walls());
    let mut polls = 0;
    admitted
        .step_with_reservoirs(
            Time::new(0.001),
            StepControl::default(),
            [Some(reservoir.as_cell()), None],
            || {
                polls += 1;
                false
            },
        )
        .unwrap();
    for stop in [1, polls / 2, polls] {
        let mut cancelled_flow = flow(gas.clone(), walls());
        let mut count = 0;
        assert_eq!(
            cancelled_flow
                .step_with_reservoirs(
                    Time::new(0.001),
                    StepControl::default(),
                    [Some(reservoir.as_cell()), None],
                    || {
                        count += 1;
                        count == stop
                    }
                )
                .err(),
            Some(Error::Cancelled)
        );
        assert_eq!(cancelled_flow.states, before);
        assert_eq!(cancelled_flow.time, 0.0);
        cancelled_flow
            .step_with_reservoirs(
                Time::new(0.001),
                StepControl::default(),
                [Some(reservoir.as_cell()), None],
                || false,
            )
            .unwrap();
        for (replayed, expected) in cancelled_flow.states.iter().zip(&admitted.states) {
            close(*replayed, *expected, 3e-12);
        }
        assert_eq!(cancelled_flow.time, admitted.time);
    }
    actual
        .step_with_reservoirs(
            Time::new(0.001),
            StepControl::default(),
            [Some(reservoir.as_cell()), None],
            || false,
        )
        .unwrap();
    assert_eq!(actual.states, admitted.states);
    assert_eq!(actual.time, admitted.time);
}
