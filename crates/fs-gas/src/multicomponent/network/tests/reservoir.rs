//! External pressure drives coexist with the actual finite-plenum RK owner.
use super::*;
use crate::multicomponent::Reconstruction;

fn drives<'a>(
    supply: &'a Conserved,
    sink: &'a Conserved,
    reversed: bool,
) -> [[Option<Cell<'a>>; 2]; 3] {
    let mut ends = [[None; 2]; 3];
    ends[0][usize::from(reversed)] = Some(supply.as_cell());
    ends[1][usize::from(!reversed)] = Some(sink.as_cell());
    ends
}

fn near(a: f64, b: f64) {
    close(a, b, a.abs().max(b.abs()).max(1.0), 3e-12);
}

fn fields_near(a: &MixtureNetwork1d, b: &MixtureNetwork1d, reversed: bool) {
    assert_eq!(a.time(), b.time());
    for (left, right) in a.branches().iter().zip(b.branches()) {
        assert_eq!(left.time(), a.time());
        assert_eq!(right.time(), b.time());
        for i in 0..left.cell_count() {
            let p = left.cell(i).unwrap();
            let q = right
                .cell(if reversed {
                    right.cell_count() - 1 - i
                } else {
                    i
                })
                .unwrap();
            for k in 0..2 {
                near(p.partial_densities[k], q.partial_densities[k]);
            }
            near(p.momentum, if reversed { -q.momentum } else { q.momentum });
            near(p.energy, q.energy);
        }
    }
    for (left, right) in a.plenums().iter().zip(b.plenums()) {
        near(left.volume_m3(), right.volume_m3());
        for k in 0..2 {
            near(
                left.state().partial_densities[k],
                right.state().partial_densities[k],
            );
        }
        near(left.state().energy, right.state().energy);
        assert_eq!(left.state().momentum, 0.0);
        assert_eq!(right.state().momentum, 0.0);
    }
}

fn check_balance(network: &MixtureNetwork1d, before: ([f64; 2], f64), report: NetworkStepReport) {
    let after = integrated(network);
    for k in 0..2 {
        near(
            after.0[k] - before.0[k],
            report.boundary_change.species_masses[k],
        );
        near(report.conservation_residual.species_masses[k], 0.0);
    }
    near(
        after.1 - before.1,
        report
            .boundary_change
            .absolute_energy(network.gas())
            .unwrap()
            + report.boundary_work_on_gas_j,
    );
    near(report.conservation_residual.energy, 0.0);
    for (connection, port) in network
        .connections()
        .iter()
        .zip(network.accepted_port_exchanges().unwrap())
    {
        let branch = &network.accepted_branch_reports().unwrap()[connection.branch];
        let into_branch = match connection.end {
            BranchEnd::Left => branch.left_boundary_change,
            BranchEnd::Right => branch.right_boundary_change,
        };
        for k in 0..2 {
            assert_eq!(
                port.species_masses_to_plenum[k],
                -into_branch.species_masses[k]
            );
        }
        assert_eq!(port.energy_to_plenum, -into_branch.energy);
        assert_eq!(port.momentum_to_plenum, -into_branch.momentum);
    }
    for branch in network.branches() {
        assert_eq!(branch.time(), network.time());
    }
}

#[test]
fn network_absent_reservoirs_preserve_reports_polls_and_clear_previous_drives() {
    let gas = caloric_gas(true, 0.0);
    let supply = gas.conserved(&[0.0, 2.0], 0.0, 2.0).unwrap();
    let sink = gas.conserved(&[1.0, 0.0], 0.0, 0.7).unwrap();
    let sources = drives(&supply, &sink, false);
    let mut original = connected(&gas, 4, false);
    let mut explicit = connected(&gas, 4, false);
    // The third call must not retain either endpoint flux from the second.
    for active in [false, true, false, false] {
        let (mut a_polls, mut b_polls) = (0, 0);
        let mut a_cancel = || {
            a_polls += 1;
            false
        };
        let a = if active {
            original.step_with_reservoirs(
                Time::new(0.002),
                StepControl::default(),
                &sources,
                &mut a_cancel,
            )
        } else {
            original.step(Time::new(0.002), StepControl::default(), &mut a_cancel)
        }
        .unwrap();
        let b = explicit
            .step_with_reservoirs(
                Time::new(0.002),
                StepControl::default(),
                if active { &sources } else { &[[None; 2]; 3] },
                || {
                    b_polls += 1;
                    false
                },
            )
            .unwrap();
        assert_eq!(a_polls, b_polls);
        assert_eq!(a, b);
        assert_eq!(snapshot(&original), snapshot(&explicit));
    }
}

#[test]
fn network_pressure_drive_and_finite_plenum_refresh_both_stages_and_close_actual_storage() {
    pressure_drive_and_finite_plenum_oracle(0.0);
}

fn pressure_drive_and_finite_plenum_oracle(volume_rate: f64) {
    let gas = caloric_gas(true, 0.0);
    let initial = gas.conserved(&[1.0, 0.0], 0.0, 1.0).unwrap();
    let supply = gas.conserved(&[0.0, 2.0], 0.0, 2.0).unwrap();
    let finite = gas.conserved(&[1.0, 0.0], 0.0, 0.5).unwrap();
    let duct = MixtureEuler1d::new(
        gas.clone(),
        Length::new(0.05),
        0.7,
        vec![initial.clone(); 2],
        walls(),
        || false,
    )
    .unwrap();
    let mut network = MixtureNetwork1d::new(
        vec![duct],
        vec![Plenum::new(0.1, finite.clone()).unwrap()],
        vec![PlenumConnection {
            branch: 0,
            end: BranchEnd::Right,
            plenum: 0,
        }],
        || false,
    )
    .unwrap();
    let flux = |inside: &Conserved, outside: &Conserved, normal| {
        let mut result = [0.0; 4];
        gas.reservoir_flux(
            inside.as_cell(),
            outside.as_cell(),
            normal,
            &mut result,
            || false,
        )
        .unwrap();
        result
    };
    let first_left = flux(&initial, &supply, -1.0);
    let first_right = flux(&initial, &finite, 1.0);
    let before = integrated(&network);
    let report = network
        .step_with_reservoirs_and_volume_rates(
            Time::new(0.001),
            StepControl::default(),
            &[[Some(supply.as_cell()), None]],
            &[volume_rate],
            || false,
        )
        .unwrap();
    assert_eq!(report.halvings, 0);
    check_balance(&network, before, report);
    let factor = report.elapsed.value() / 0.05;
    let final_volume = 0.1 + report.elapsed.value() * volume_rate;
    let mut left = initial.clone();
    let mut right = initial.clone();
    let mut chamber = finite.clone();
    for k in 0..2 {
        left.partial_densities[k] += factor * first_left[k];
        right.partial_densities[k] -= factor * first_right[k];
        chamber.partial_densities[k] = (0.1 * finite.partial_densities[k]
            + report.elapsed.value() * 0.7 * first_right[k])
            / final_volume;
    }
    left.momentum = factor * (first_left[2] - 1.0);
    right.momentum = factor * (1.0 - first_right[2]);
    left.energy += factor * first_left[3];
    right.energy -= factor * first_right[3];
    chamber.energy = (0.1 * finite.energy
        + report.elapsed.value() * (0.7 * first_right[3] - 0.5 * volume_rate))
        / final_volume;
    let second_left = flux(&left, &supply, -1.0);
    let second_right = flux(&right, &chamber, 1.0);
    let second_pressure = gas.primitive(chamber.as_cell(), || false).unwrap().pressure;
    let work = -0.5 * report.elapsed.value() * volume_rate * (0.5 + second_pressure);
    near(report.boundary_work_on_gas_j, work);
    let mechanical = network.accepted_plenum_boundary_work().unwrap()[0];
    assert_eq!(mechanical.initial_volume_m3, 0.1);
    assert_eq!(mechanical.final_volume_m3, final_volume);
    assert_eq!(network.plenums()[0].volume_m3(), final_volume);
    near(mechanical.work_on_gas_j, work);
    assert!((second_left[1] - first_left[1]).abs() > 1e-7);
    assert!((second_right[0] - first_right[0]).abs() > 1e-7);
    let scale = 0.5 * report.elapsed.value() * 0.7;
    let branch = &network.accepted_branch_reports().unwrap()[0];
    for k in 0..4 {
        near(
            branch.left_boundary_change.get(k),
            scale * (first_left[k] + second_left[k]),
        );
        near(
            branch.right_boundary_change.get(k),
            -scale * (first_right[k] + second_right[k]),
        );
    }
    for k in 0..2 {
        near(
            report.boundary_change.species_masses[k],
            branch.left_boundary_change.species_masses[k],
        );
        near(
            final_volume * network.plenums()[0].state().partial_densities[k]
                - 0.1 * finite.partial_densities[k],
            network.accepted_port_exchanges().unwrap()[0].species_masses_to_plenum[k],
        );
    }
    near(
        report.boundary_change.energy,
        branch.left_boundary_change.energy,
    );
    near(
        final_volume * network.plenums()[0].state().energy - 0.1 * finite.energy,
        scale * (first_right[3] + second_right[3]) + work,
    );
    assert_eq!(report.boundary_change.species_masses[0], 0.0);
    assert!(report.boundary_change.species_masses[1] > 0.0);
    assert!(network.branches()[0].cell(0).unwrap().partial_densities[1] > 0.0);
    assert!(
        final_volume * network.plenums()[0].state().density() > 0.1 * finite.as_cell().density()
    );
    // Pure B: R=0.5, Cp=2.5+0.35T, T0=2, with e(Tmin=0.001)
    // subtracted once. The external energy is the true donor stagnation enthalpy.
    let h0 = 2.5 * 2.0 + 0.175 * 4.0 - 2.0 * 0.001 - 0.175 * 0.001_f64.powi(2);
    near(
        report.boundary_change.energy / report.boundary_change.mass(),
        h0,
    );
}

#[test]
fn network_reservoir_and_volume_drives_share_stages_work_and_atomic_publication() {
    for rate in [-0.02, 0.03] {
        pressure_drive_and_finite_plenum_oracle(rate);
    }
    let gas = caloric_gas(true, 0.0);
    let supply = gas.conserved(&[0.0, 2.0], 0.0, 2.0).unwrap();
    let sink = gas.conserved(&[1.0, 0.0], 0.0, 0.7).unwrap();
    let sources = drives(&supply, &sink, false);
    // Each degenerate combined path preserves its original owner's arithmetic.
    for moving in [false, true] {
        let mut original = connected(&gas, 4, false);
        let mut combined = connected(&gas, 4, false);
        let rates = [if moving { -0.02 } else { 0.0 }];
        let ends = if moving { [[None; 2]; 3] } else { sources };
        let a = if moving {
            original.step_with_volume_rates(
                Time::new(0.002),
                StepControl::default(),
                &rates,
                || false,
            )
        } else {
            original.step_with_reservoirs(Time::new(0.002), StepControl::default(), &ends, || false)
        }
        .unwrap();
        let b = combined
            .step_with_reservoirs_and_volume_rates(
                Time::new(0.002),
                StepControl::default(),
                &ends,
                &rates,
                || false,
            )
            .unwrap();
        assert_eq!(a, b);
        assert_eq!(snapshot(&original), snapshot(&combined));
        assert_eq!(original.plenums(), combined.plenums());
        assert_eq!(
            original.accepted_plenum_boundary_work(),
            combined.accepted_plenum_boundary_work()
        );
    }
    let prepare = || {
        let mut network = connected(&gas, 4, false);
        network
            .step_with_reservoirs_and_volume_rates(
                Time::new(0.002),
                StepControl::default(),
                &sources,
                &[-0.02],
                || false,
            )
            .unwrap();
        network
    };
    let mut expected = prepare();
    let mut polls = 0;
    expected
        .step_with_reservoirs_and_volume_rates(
            Time::new(0.002),
            StepControl::default(),
            &sources,
            &[-0.02],
            || {
                polls += 1;
                false
            },
        )
        .unwrap();
    for stop in [1, polls / 2, polls] {
        let mut actual = prepare();
        let before = snapshot(&actual);
        let volumes = actual.plenums().to_vec();
        let work = actual.accepted_plenum_boundary_work().unwrap().to_vec();
        let mut visits = 0;
        assert_eq!(
            actual.step_with_reservoirs_and_volume_rates(
                Time::new(0.002),
                StepControl::default(),
                &sources,
                &[-0.02],
                || {
                    visits += 1;
                    visits == stop
                }
            ),
            Err(Error::Cancelled)
        );
        assert_eq!(snapshot(&actual), before);
        assert_eq!(actual.plenums(), volumes);
        assert_eq!(actual.accepted_plenum_boundary_work().unwrap(), work);
        actual
            .step_with_reservoirs_and_volume_rates(
                Time::new(0.002),
                StepControl::default(),
                &sources,
                &[-0.02],
                || false,
            )
            .unwrap();
        fields_near(&actual, &expected, false);
    }
}

#[test]
fn network_reservoir_end_reversal_and_species_formation_gauges_preserve_physics() {
    let run = |formation, reversed| {
        let gas = caloric_gas(true, formation);
        let supply = gas.conserved(&[0.0, 2.0], 0.0, 2.0).unwrap();
        let sink = gas.conserved(&[1.0, 0.0], 0.0, 0.7).unwrap();
        let mut network = connected(&gas, 4, reversed);
        for _ in 0..4 {
            let before = integrated(&network);
            let report = network
                .step_with_reservoirs(
                    Time::new(0.002),
                    StepControl::default(),
                    &drives(&supply, &sink, reversed),
                    || false,
                )
                .unwrap();
            if formation == 0.0 {
                check_balance(&network, before, report);
            }
            let branch = &network.accepted_branch_reports().unwrap()[1];
            let outward = if reversed {
                branch.left_boundary_change
            } else {
                branch.right_boundary_change
            };
            assert!(outward.mass() < 0.0);
            // The low-pressure source is pure A, but outgoing daughter gas is
            // initially pure B and must continue to carry its interior donor.
            assert!(outward.species_masses[1] < 0.0);
        }
        network
    };
    let reference = run(0.0, false);
    for formation in [-1e18, 1e18] {
        assert_eq!(snapshot(&reference), snapshot(&run(formation, false)));
    }
    let reversed = run(0.0, true);
    fields_near(&reference, &reversed, true);
    for (a, b) in reference
        .accepted_port_exchanges()
        .unwrap()
        .iter()
        .zip(reversed.accepted_port_exchanges().unwrap())
    {
        for k in 0..2 {
            near(a.species_masses_to_plenum[k], b.species_masses_to_plenum[k]);
        }
        near(a.energy_to_plenum, b.energy_to_plenum);
        near(a.momentum_to_plenum, -b.momentum_to_plenum);
    }
}

#[test]
fn network_reservoir_cfl_retry_replays_all_owners_from_accepted_storage() {
    let gas = caloric_gas(false, 0.0);
    let supply = gas.conserved(&[0.0, 100.0], 0.0, 100.0).unwrap();
    let sink = gas.conserved(&[1.0, 0.0], 0.0, 0.7).unwrap();
    let sources = drives(&supply, &sink, false);
    let mut actual = connected(&gas, 4, false);
    actual
        .step(Time::new(0.002), StepControl::default(), || false)
        .unwrap();
    let before = snapshot(&actual);
    assert!(matches!(
        actual.step_with_reservoirs(
            Time::new(1.0),
            StepControl {
                maximum_halvings: 0,
                ..StepControl::default()
            },
            &sources,
            || false
        ),
        Err(Error::StepRetries { .. })
    ));
    assert_eq!(snapshot(&actual), before);
    let stored = integrated(&actual);
    let accepted = actual
        .step_with_reservoirs(Time::new(1.0), StepControl::default(), &sources, || false)
        .unwrap();
    assert!(accepted.halvings > 0);
    assert!(accepted.maximum_cfl <= 0.4 * (1.0 + 64.0 * f64::EPSILON));
    check_balance(&actual, stored, accepted);
    let mut replay = connected(&gas, 4, false);
    replay
        .step(Time::new(0.002), StepControl::default(), || false)
        .unwrap();
    let repeated = replay
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
            &sources,
            || false,
        )
        .unwrap();
    assert_eq!(repeated.halvings, 0);
    fields_near(&actual, &replay, false);
    for k in 0..2 {
        near(
            accepted.boundary_change.species_masses[k],
            repeated.boundary_change.species_masses[k],
        );
    }
    near(
        accepted.boundary_change.energy,
        repeated.boundary_change.energy,
    );
}

#[test]
fn network_reservoir_admission_refuses_shape_connected_periodic_and_invalid_source_ends() {
    let gas = caloric_gas(true, 0.0);
    let supply = gas.conserved(&[0.0, 2.0], 0.0, 2.0).unwrap();
    let mut actual = connected(&gas, 4, false);
    actual
        .step(Time::new(0.002), StepControl::default(), || false)
        .unwrap();
    let before = snapshot(&actual);
    for count in [0, 2, 4] {
        assert!(
            actual
                .step_with_reservoirs(
                    Time::new(0.002),
                    StepControl::default(),
                    &vec![[None; 2]; count],
                    || false
                )
                .is_err()
        );
        assert_eq!(snapshot(&actual), before);
    }
    // The first source is valid; the next targets an internally owned end.
    let conflicting = [
        [Some(supply.as_cell()), None],
        [Some(supply.as_cell()), None],
        [None; 2],
    ];
    assert!(
        actual
            .step_with_reservoirs(
                Time::new(0.002),
                StepControl::default(),
                &conflicting,
                || false
            )
            .is_err()
    );
    let mut moving = supply.clone();
    moving.momentum = 0.01;
    let mut foreign_axis = supply.clone();
    foreign_axis.partial_densities.pop();
    let mut negative = supply.clone();
    negative.partial_densities[0] = -1.0;
    let mut nonfinite = supply.clone();
    nonfinite.energy = f64::NAN;
    let mut outside_domain = supply.clone();
    outside_domain.energy = 1e20;
    for invalid in [moving, foreign_axis, negative, nonfinite, outside_domain] {
        assert!(
            actual
                .step_with_reservoirs(
                    Time::new(0.002),
                    StepControl::default(),
                    &[[Some(invalid.as_cell()), None], [None; 2], [None; 2]],
                    || false
                )
                .is_err()
        );
        assert_eq!(snapshot(&actual), before);
    }
    let state = gas.conserved(&[1.0, 0.0], 0.0, 1.0).unwrap();
    let periodic = MixtureEuler1d::new(
        gas,
        Length::new(0.1),
        0.7,
        vec![state; 2],
        Boundaries::Periodic,
        || false,
    )
    .unwrap();
    let mut periodic = MixtureNetwork1d::new(vec![periodic], vec![], vec![], || false).unwrap();
    periodic
        .step_with_reservoirs(
            Time::new(0.002),
            StepControl::default(),
            &[[None; 2]],
            || false,
        )
        .unwrap();
    let saved = snapshot(&periodic);
    assert!(
        periodic
            .step_with_reservoirs(
                Time::new(0.002),
                StepControl::default(),
                &[[Some(supply.as_cell()), None]],
                || false
            )
            .is_err()
    );
    assert_eq!(snapshot(&periodic), saved);
}

#[test]
fn network_reservoir_cancellation_preserves_every_field_clock_and_accepted_report() {
    let gas = caloric_gas(true, 0.0);
    let supply = gas.conserved(&[0.0, 2.0], 0.0, 2.0).unwrap();
    let sink = gas.conserved(&[1.0, 0.0], 0.0, 0.7).unwrap();
    let sources = drives(&supply, &sink, false);
    let mut expected = connected(&gas, 4, false);
    expected
        .step(Time::new(0.002), StepControl::default(), || false)
        .unwrap();
    let mut polls = 0;
    expected
        .step_with_reservoirs(Time::new(0.002), StepControl::default(), &sources, || {
            polls += 1;
            false
        })
        .unwrap();
    for stop in [1, polls / 2, polls] {
        let mut actual = connected(&gas, 4, false);
        actual
            .step(Time::new(0.002), StepControl::default(), || false)
            .unwrap();
        let before = snapshot(&actual);
        let mut visits = 0;
        assert_eq!(
            actual.step_with_reservoirs(Time::new(0.002), StepControl::default(), &sources, || {
                visits += 1;
                visits == stop
            }),
            Err(Error::Cancelled)
        );
        assert_eq!(snapshot(&actual), before);
        for branch in actual.branches() {
            assert_eq!(branch.time().value(), before.2);
        }
        actual
            .step_with_reservoirs(Time::new(0.002), StepControl::default(), &sources, || false)
            .unwrap();
        fields_near(&actual, &expected, false);
    }
}
