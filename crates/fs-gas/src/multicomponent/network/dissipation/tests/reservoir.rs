//! Exterior baths, finite storage and complete wave/diffusion transactions.
use super::*;
use crate::multicomponent::network::dissipating::{DissipatingNetwork1d, DissipatingStepControl};
use crate::multicomponent::network::{BranchEnd, Plenum, PlenumConnection};

pub(in crate::multicomponent::network) fn bath_network(
    gas: &MixtureGas,
    reversed: bool,
) -> MixtureNetwork1d {
    let branch = MixtureEuler1d::new(
        gas.clone(),
        Length::new(0.4),
        0.2,
        vec![gas.conserved(&[0.4, 0.6], 0.0, 2.0).unwrap(); 2],
        Boundaries::Ends {
            left: Boundary::Reflecting,
            right: Boundary::Reflecting,
        },
        || false,
    )
    .unwrap();
    MixtureNetwork1d::new(
        vec![branch],
        vec![Plenum::new(0.07, gas.conserved(&[0.1, 0.9], 0.0, 2.0).unwrap()).unwrap()],
        vec![PlenumConnection {
            branch: 0,
            end: if reversed {
                BranchEnd::Left
            } else {
                BranchEnd::Right
            },
            plenum: 0,
        }],
        || false,
    )
    .unwrap()
}

// Independent three-storage finite-volume ODE, with rho=1, equal molecular
// masses/Cp, D=.03, A=.2, dx=.4, V_plenum=.07, port mixing length=.15.
// x[0] is the exterior cell, x[1] the port cell and x[2] the finite bath.
// The fixed exterior bath has x_A=.9. Its molecular transfer has no bulk flow.
pub(in crate::multicomponent::network) fn bath_step(x: [f64; 3], dt: f64) -> ([f64; 3], f64) {
    let rhs = |q: [f64; 3]| {
        let left = 0.03 * (0.9 - q[0]) / 0.2;
        let middle = 0.03 * (q[0] - q[1]) / 0.4;
        let right = 0.03 * (q[1] - q[2]) / 0.35;
        [
            (left - middle) / 0.4,
            (middle - right) / 0.4,
            0.2 * right / 0.07,
        ]
    };
    let a = rhs(x);
    let predictor = core::array::from_fn(|i| x[i] + dt * a[i]);
    let b = rhs(predictor);
    let next = core::array::from_fn(|i| x[i] + 0.5 * dt * (a[i] + b[i]));
    let exterior = 0.5 * dt * 0.2 * 0.03 * ((0.9 - x[0]) + (0.9 - predictor[0])) / 0.2;
    (next, exterior)
}

fn baths(cell: Cell<'_>, reversed: bool) -> [[Option<Cell<'_>>; 2]; 1] {
    if reversed {
        [[None, Some(cell)]]
    } else {
        [[Some(cell), None]]
    }
}

#[test]
fn external_bath_and_finite_plenum_match_independent_conservative_three_storage_step() {
    let dt = 0.04;
    let (expected, exterior) = bath_step([0.4, 0.4, 0.1], dt);
    for reversed in [false, true] {
        let mut baseline = None;
        for reference in [0.0, 1e18] {
            let gas = gas(2, false, reference);
            let bath = gas.conserved(&[0.9, 0.1], 0.0, 2.0).unwrap();
            let mut network = bath_network(&gas, reversed);
            let mut op = NetworkDissipationOperator::new(
                &network,
                law(&gas, 0.0, 0.0, 0.0, 0.03),
                vec![Length::new(0.15)],
                || false,
            )
            .unwrap();
            let report = op
                .advance_with_reservoirs_for(
                    &mut network,
                    Time::new(dt),
                    DissipationControl::default(),
                    &baths(bath.as_cell(), reversed),
                    || false,
                )
                .unwrap();
            assert_eq!(report.accepted_steps, 1);
            assert_eq!(network.time().value(), 0.0);
            for (i, &value) in expected[..2].iter().enumerate() {
                let q = network.branches()[0]
                    .cell(if reversed { 1 - i } else { i })
                    .unwrap();
                close(q.partial_densities[0], value, 3e-13);
                close(q.partial_densities[1], 1.0 - value, 3e-13);
                close(gas.primitive(q, || false).unwrap().temperature, 2.0, 3e-12);
            }
            close(
                network.plenums()[0].state().partial_densities[0],
                expected[2],
                3e-13,
            );
            close(report.boundary_change.species_masses[0], exterior, 3e-14);
            close(report.boundary_change.species_masses[1], -exterior, 3e-14);
            close(report.boundary_change.mass(), 0.0, 3e-14);
            close(report.boundary_change.energy, 0.0, 3e-13);
            close(report.conservation_residual.species_masses[0], 0.0, 3e-14);
            close(report.conservation_residual.energy, 0.0, 3e-13);
            let port = op.accepted_port_exchanges().unwrap()[0];
            close(
                port.species_masses_to_plenum[0],
                0.07 * (expected[2] - 0.1),
                3e-14,
            );
            let branch = op.accepted_branch_reports().unwrap()[0];
            let (external, internal) = if reversed {
                (branch.right_boundary_change, branch.left_boundary_change)
            } else {
                (branch.left_boundary_change, branch.right_boundary_change)
            };
            close(external.species_masses[0], exterior, 3e-14);
            close(
                internal.species_masses[0],
                -port.species_masses_to_plenum[0],
                3e-14,
            );
            if let Some((original, boundary)) = &baseline {
                assert_eq!(&fields(&network), original);
                assert_eq!(&report.boundary_change, boundary);
            } else {
                baseline = Some((fields(&network), report.boundary_change));
            }
        }
    }
}

#[test]
fn thermal_bath_and_finite_plenum_close_the_independent_fourier_storage_balance() {
    let gas = gas(2, false, 0.0);
    let mut network = bath_network(&gas, false);
    let branch = gas.conserved(&[0.4, 0.6], 0.0, 1.4).unwrap();
    for q in network.branches[0].states.chunks_exact_mut(4) {
        q.copy_from_slice(&[0.4, 0.6, 0.0, branch.energy]);
    }
    network.plenums[0].state = gas.conserved(&[0.4, 0.6], 0.0, 1.1).unwrap();
    let bath = gas.conserved(&[0.4, 0.6], 0.0, 1.9).unwrap();
    // rho*Cv=2.5, lambda=.075: alpha=.03, the same independently specified
    // three-storage ODE for T-1. The pressure gradient creates no bulk flow in D.
    let mut op = NetworkDissipationOperator::new(
        &network,
        law(&gas, 0.0, 0.0, 0.075, 0.03),
        vec![Length::new(0.15)],
        || false,
    )
    .unwrap();
    let report = op
        .advance_with_reservoirs_for(
            &mut network,
            Time::new(0.04),
            DissipationControl::default(),
            &baths(bath.as_cell(), false),
            || false,
        )
        .unwrap();
    let (expected, transfer) = bath_step([0.4, 0.4, 0.1], 0.04);
    for (i, &temperature) in expected[..2].iter().enumerate() {
        let q = network.branches()[0].cell(i).unwrap();
        close(
            gas.primitive(q, || false).unwrap().temperature,
            1.0 + temperature,
            3e-12,
        );
        close(q.partial_densities[0], 0.4, 3e-14);
        close(q.partial_densities[1], 0.6, 3e-14);
    }
    close(
        gas.primitive(network.plenums()[0].state(), || false)
            .unwrap()
            .temperature,
        1.0 + expected[2],
        3e-12,
    );
    close(report.boundary_change.energy, 2.5 * transfer, 3e-13);
    close(
        report.after.energy - report.before.energy,
        2.5 * transfer,
        3e-13,
    );
    close(
        op.accepted_port_exchanges().unwrap()[0].energy_to_plenum,
        0.07 * 2.5 * (expected[2] - 0.1),
        3e-13,
    );
    close(report.absolute_energy_residual_j, 0.0, 3e-13);
    close(report.boundary_change.mass(), 0.0, 3e-14);
}

#[test]
fn existing_wave_reservoir_api_advances_both_molecular_halves_without_pressure_flow() {
    let gas = gas(2, false, 0.0);
    let bath = gas.conserved(&[0.9, 0.1], 0.0, 2.0).unwrap();
    let mut flow = DissipatingNetwork1d::new(
        bath_network(&gas, false),
        law(&gas, 0.0, 0.0, 0.0, 0.03),
        vec![Length::new(0.15)],
        || false,
    )
    .unwrap();
    let report = flow
        .step_with_reservoirs(
            Time::new(0.04),
            DissipatingStepControl::default(),
            &baths(bath.as_cell(), false),
            || false,
        )
        .unwrap();
    let dt = report.elapsed.value();
    let (half, first) = bath_step([0.4, 0.4, 0.1], 0.5 * dt);
    let (expected, second) = bath_step(half, 0.5 * dt);
    for (i, &value) in expected[..2].iter().enumerate() {
        close(
            flow.branches()[0].cell(i).unwrap().partial_densities[0],
            value,
            3e-12,
        );
    }
    close(
        flow.plenums()[0].state().partial_densities[0],
        expected[2],
        3e-12,
    );
    close(
        report.molecular.first.boundary_change.species_masses[0],
        first,
        3e-13,
    );
    close(
        report.molecular.second.boundary_change.species_masses[0],
        second,
        3e-13,
    );
    assert!(first > 0.0 && second > 0.0 && first != second);
    close(
        report.transport.boundary_change.species_masses[0],
        0.0,
        3e-13,
    );
    close(
        report.boundary_change.species_masses[0],
        first + second,
        3e-13,
    );
    close(report.conservation_residual.species_masses[0], 0.0, 3e-13);
    close(report.absolute_energy_residual_j, 0.0, 3e-12);
    close(flow.time().value(), dt, 1e-16);
    let d = flow.accepted_molecular_port_exchanges().unwrap()[0];
    let a = flow.accepted_port_exchanges().unwrap()[0];
    close(
        d.species_masses_to_plenum[0] + a.species_masses_to_plenum[0],
        0.07 * (expected[2] - 0.1),
        3e-13,
    );

    // A final-poll refusal after both bath-driven D halves must preserve the
    // accepted network and its composed molecular/advective reports together.
    let mut oracle = DissipatingNetwork1d::new(
        bath_network(&gas, false),
        law(&gas, 0.0, 0.0, 0.0, 0.03),
        vec![Length::new(0.15)],
        || false,
    )
    .unwrap();
    oracle
        .step_with_reservoirs(
            Time::new(0.04),
            DissipatingStepControl::default(),
            &baths(bath.as_cell(), false),
            || false,
        )
        .unwrap();
    let mut polls = 0;
    oracle
        .step_with_reservoirs(
            Time::new(0.04),
            DissipatingStepControl::default(),
            &baths(bath.as_cell(), false),
            || {
                polls += 1;
                false
            },
        )
        .unwrap();
    let snapshot = |owner: &DissipatingNetwork1d<ConstantTransport>| {
        let branch = (0..2)
            .map(|i| owned(owner.branches()[0].cell(i).unwrap()))
            .collect::<Vec<_>>();
        (
            branch,
            owned(owner.plenums()[0].state()),
            owner.time(),
            owner.last_report(),
            owner.accepted_molecular_port_exchanges().unwrap().to_vec(),
            owner.accepted_port_exchanges().unwrap().to_vec(),
        )
    };
    let before = snapshot(&flow);
    let mut seen = 0;
    assert!(matches!(
        flow.step_with_reservoirs(
            Time::new(0.04),
            DissipatingStepControl::default(),
            &baths(bath.as_cell(), false),
            || {
                seen += 1;
                seen == polls
            }
        ),
        Err(crate::multicomponent::network::dissipating::DissipatingNetworkError::Cancelled)
    ));
    assert_eq!(snapshot(&flow), before);
    flow.step_with_reservoirs(
        Time::new(0.04),
        DissipatingStepControl::default(),
        &baths(bath.as_cell(), false),
        || false,
    )
    .unwrap();
    assert_eq!(snapshot(&flow), snapshot(&oracle));
}

#[test]
fn network_bath_refusal_and_late_cancellation_preserve_all_storage_and_accepted_reports() {
    let gas = gas(2, false, 0.0);
    let bath = gas.conserved(&[0.9, 0.1], 0.0, 2.0).unwrap();
    let controls = DissipationControl {
        maximum_step_s: 0.01,
        ..Default::default()
    };
    let make = || {
        let network = bath_network(&gas, false);
        let op = NetworkDissipationOperator::new(
            &network,
            law(&gas, 0.0, 0.0, 0.0, 0.03),
            vec![Length::new(0.15)],
            || false,
        )
        .unwrap();
        (network, op)
    };
    let (mut network, mut op) = make();
    let (mut expected, mut oracle) = make();
    for (n, o) in [(&mut network, &mut op), (&mut expected, &mut oracle)] {
        o.advance_with_reservoirs_for(
            n,
            Time::new(0.02),
            controls,
            &baths(bath.as_cell(), false),
            || false,
        )
        .unwrap();
    }
    let mut polls = 0;
    let expected_report = oracle
        .advance_with_reservoirs_for(
            &mut expected,
            Time::new(0.04),
            controls,
            &baths(bath.as_cell(), false),
            || {
                polls += 1;
                false
            },
        )
        .unwrap();
    let initial = fields(&network);
    let last = op.last_report();
    let ports = op.accepted_port_exchanges().unwrap().to_vec();
    let branches = op.accepted_branch_reports().unwrap().to_vec();
    let mut seen = 0;
    assert_eq!(
        op.advance_with_reservoirs_for(
            &mut network,
            Time::new(0.04),
            controls,
            &baths(bath.as_cell(), false),
            || {
                seen += 1;
                seen == polls
            }
        ),
        Err(DissipationError::Cancelled)
    );
    assert_eq!(fields(&network), initial);
    assert_eq!(op.last_report(), last);
    assert_eq!(op.accepted_port_exchanges().unwrap(), ports);
    assert_eq!(op.accepted_branch_reports().unwrap(), branches);
    // A prescribed source must never silently replace a finite-capacity port.
    for invalid in [
        vec![],
        vec![[None, Some(bath.as_cell())]],
        vec![[
            Some(Cell {
                momentum: 0.1,
                ..bath.as_cell()
            }),
            None,
        ]],
    ] {
        assert!(
            op.advance_with_reservoirs_for(
                &mut network,
                Time::new(0.04),
                controls,
                &invalid,
                || false
            )
            .is_err()
        );
        assert_eq!(fields(&network), initial);
        assert_eq!(op.last_report(), last);
    }
    assert!(matches!(
        op.advance_with_reservoirs_for(
            &mut network,
            Time::new(0.04),
            DissipationControl {
                maximum_steps: 1,
                ..controls
            },
            &baths(bath.as_cell(), false),
            || false
        ),
        Err(DissipationError::WorkLimit { .. })
    ));
    assert_eq!(fields(&network), initial);
    assert_eq!(op.last_report(), last);
    let report = op
        .advance_with_reservoirs_for(
            &mut network,
            Time::new(0.04),
            controls,
            &baths(bath.as_cell(), false),
            || false,
        )
        .unwrap();
    assert_eq!(fields(&network), fields(&expected));
    assert_eq!(report, expected_report);
}

#[test]
fn absent_network_baths_preserve_existing_molecular_fields_reports_and_polls() {
    let gas = gas(2, false, 0.0);
    let mut a = bath_network(&gas, false);
    let mut b = bath_network(&gas, false);
    let build = |n: &MixtureNetwork1d| {
        NetworkDissipationOperator::new(
            n,
            law(&gas, 0.0, 0.0, 0.0, 0.03),
            vec![Length::new(0.15)],
            || false,
        )
        .unwrap()
    };
    let (mut oa, mut ob) = (build(&a), build(&b));
    let (mut pa, mut pb) = (0, 0);
    let ra = oa
        .advance_for(
            &mut a,
            Time::new(0.02),
            DissipationControl::default(),
            || {
                pa += 1;
                false
            },
        )
        .unwrap();
    let rb = ob
        .advance_with_reservoirs_for(
            &mut b,
            Time::new(0.02),
            DissipationControl::default(),
            &[[None; 2]],
            || {
                pb += 1;
                false
            },
        )
        .unwrap();
    assert_eq!(ra, rb);
    assert_eq!(pa, pb);
    assert_eq!(fields(&a), fields(&b));
}
