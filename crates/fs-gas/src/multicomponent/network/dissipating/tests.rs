//! Complete physical exchange and atomic D/A/D network composition.
use super::super::dissipation::tests::{close, fields, gas, law, owned};
use super::*;
use crate::multicomponent::{Boundaries, Boundary, Conserved};

fn walls() -> Boundaries {
    Boundaries::Ends {
        left: Boundary::Reflecting,
        right: Boundary::Reflecting,
    }
}

fn closed(g: &MixtureGas) -> MixtureNetwork1d {
    let states = (0..8)
        .map(|i| {
            g.conserved(&[0.9, 0.1], 0.0, if i < 4 { 1.6 } else { 0.8 })
                .unwrap()
        })
        .collect();
    let branch =
        MixtureEuler1d::new(g.clone(), Length::new(0.1), 0.08, states, walls(), || false).unwrap();
    MixtureNetwork1d::new(
        vec![branch],
        vec![Plenum::new(0.005, g.conserved(&[0.2, 0.8], 0.0, 1.0).unwrap()).unwrap()],
        vec![PlenumConnection {
            branch: 0,
            end: super::super::BranchEnd::Right,
            plenum: 0,
        }],
        || false,
    )
    .unwrap()
}

fn owner(
    g: &MixtureGas,
) -> DissipatingNetwork1d<crate::multicomponent::dissipation::ConstantTransport> {
    DissipatingNetwork1d::new(
        closed(g),
        law(g, 0.003, 0.001, 0.02, 0.01),
        vec![Length::new(0.05)],
        || false,
    )
    .unwrap()
}

// Independent integration of the accepted fields, including physical formation
// references. No whole-step residual or network-total helper enters this sum.
fn inventories(
    flow: &DissipatingNetwork1d<crate::multicomponent::dissipation::ConstantTransport>,
) -> ([f64; 2], f64) {
    let mut mass = [0.0; 2];
    let mut energy = 0.0;
    for branch in flow.branches() {
        for (i, volume) in branch.mesh().volumes().iter().enumerate() {
            let q = branch.cell(i).unwrap();
            for k in 0..2 {
                mass[k] += volume * q.partial_densities[k];
            }
            energy += volume * flow.gas().absolute_energy_density(q).unwrap();
        }
    }
    for plenum in flow.plenums() {
        let q = plenum.state();
        for k in 0..2 {
            mass[k] += plenum.volume_m3() * q.partial_densities[k];
        }
        energy += plenum.volume_m3() * flow.gas().absolute_energy_density(q).unwrap();
    }
    (mass, energy)
}

#[test]
fn zero_molecular_law_retains_the_original_network_at_the_actual_accepted_duration() {
    let g = gas(1, false, 0.0);
    let make = || {
        let states = (0..12)
            .map(|i| {
                g.conserved(&[1.0], 0.0, if i < 6 { 1.3 } else { 1.0 })
                    .unwrap()
            })
            .collect();
        let branch = MixtureEuler1d::new(
            g.clone(),
            Length::new(1.0 / 12.0),
            0.2,
            states,
            walls(),
            || false,
        )
        .unwrap();
        MixtureNetwork1d::new(vec![branch], Vec::new(), Vec::new(), || false).unwrap()
    };
    let mut reference = make();
    let mut flow =
        DissipatingNetwork1d::new(make(), law(&g, 0.0, 0.0, 0.0, 0.0), Vec::new(), || false)
            .unwrap();
    for _ in 0..12 {
        let report = flow
            .step(Time::new(0.01), DissipatingStepControl::default(), || false)
            .unwrap();
        let plain = reference
            .step(report.elapsed, StepControl::default(), || false)
            .unwrap();
        assert_eq!(plain.elapsed, report.elapsed);
        // The composing owner swaps prepared workspaces, so the independent
        // advection owner can retain different temperature-inversion guesses.
        // Require roundoff agreement of A and exact zero-change D budgets below.
        for (actual, expected) in flow.branches().iter().zip(reference.branches()) {
            for i in 0..actual.cell_count() {
                let a = actual.cell(i).unwrap();
                let b = expected.cell(i).unwrap();
                close(
                    a.partial_densities[0],
                    b.partial_densities[0],
                    64.0 * f64::EPSILON,
                );
                close(a.momentum, b.momentum, 64.0 * f64::EPSILON);
                close(a.energy, b.energy, 64.0 * f64::EPSILON);
            }
        }
        assert_eq!(flow.time(), reference.time());
        assert_eq!(
            report.molecular.first.elapsed.value(),
            0.5 * report.elapsed.value()
        );
        assert_eq!(
            report.molecular.second.elapsed.value(),
            0.5 * report.elapsed.value()
        );
        assert_eq!(report.molecular.first.before, report.molecular.first.after);
        assert_eq!(
            report.molecular.second.before,
            report.molecular.second.after
        );
        close(report.mass_residual_kg, 0.0, 3e-13);
        close(report.absolute_energy_residual_j, 0.0, 3e-12);
    }
}

#[test]
fn shock_and_molecular_port_transfers_close_independent_plenum_and_global_budgets() {
    let g = gas(2, true, 5.0);
    let mut flow = owner(&g);
    let initial = inventories(&flow);
    let mut maximum_advective_species = 0.0_f64;
    let mut maximum_molecular_species = 0.0_f64;
    for _ in 0..20 {
        let before = owned(flow.plenums()[0].state());
        let time = flow.time().value();
        let report = flow
            .step(Time::new(0.004), DissipatingStepControl::default(), || {
                false
            })
            .unwrap();
        assert_eq!(flow.last_report(), Some(report));
        assert_eq!(report.elapsed, report.transport.elapsed);
        close(flow.time().value() - time, report.elapsed.value(), 5e-17);
        assert_eq!(
            report.molecular.first.elapsed.value(),
            0.5 * report.elapsed.value()
        );
        assert_eq!(
            report.molecular.second.elapsed.value(),
            0.5 * report.elapsed.value()
        );
        assert!(report.molecular.first.maximum_parabolic_number <= 0.40000000000001);
        assert!(report.molecular.second.maximum_parabolic_number <= 0.40000000000001);
        let a = flow.accepted_port_exchanges().unwrap()[0];
        let d = flow.accepted_molecular_port_exchanges().unwrap()[0];
        let after = flow.plenums()[0].state();
        for k in 0..2 {
            let actual = 0.005 * (after.partial_densities[k] - before.partial_densities[k]);
            close(
                actual,
                a.species_masses_to_plenum[k] + d.species_masses_to_plenum[k],
                3e-14,
            );
            assert_eq!(report.boundary_change.species_masses[k], 0.0);
            close(report.conservation_residual.species_masses[k], 0.0, 3e-13);
        }
        let actual_energy = 0.005
            * (g.absolute_energy_density(after).unwrap()
                - g.absolute_energy_density(before.as_cell()).unwrap());
        close(
            actual_energy,
            a.absolute_energy_to_plenum(&g).unwrap() + d.absolute_energy_to_plenum(&g).unwrap(),
            3e-12,
        );
        assert_eq!(report.boundary_change.energy, 0.0);
        close(report.absolute_energy_residual_j, 0.0, 3e-12);
        maximum_advective_species =
            maximum_advective_species.max(a.species_masses_to_plenum[0].abs());
        maximum_molecular_species =
            maximum_molecular_species.max(d.species_masses_to_plenum[0].abs());
        for branch in flow.branches() {
            assert_eq!(branch.time(), flow.time());
        }
    }
    assert!(
        maximum_advective_species > 1e-6,
        "no pressure-driven species transfer"
    );
    assert!(
        maximum_molecular_species > 1e-6,
        "no molecular transfer across the network port"
    );
    let final_storage = inventories(&flow);
    for k in 0..2 {
        close(final_storage.0[k], initial.0[k], 5e-13);
    }
    close(final_storage.1, initial.1, 1e-11);
    let before_release = fields(&flow.accepted);
    let clock = flow.time();
    let released = flow.into_network();
    assert_eq!(fields(&released), before_release);
    assert_eq!(released.time(), clock);
    assert!(released.last_report().is_none());
    assert!(released.accepted_port_exchanges().is_none());
}

type Snapshot = (
    (Vec<Vec<Conserved>>, Vec<Conserved>),
    Time,
    Option<DissipatingNetworkReport>,
    Option<Vec<PortExchange>>,
    Option<Vec<PortExchange>>,
    Option<Vec<StepReport>>,
);
fn snapshot(
    flow: &DissipatingNetwork1d<crate::multicomponent::dissipation::ConstantTransport>,
) -> Snapshot {
    (
        fields(&flow.accepted),
        flow.time(),
        flow.last_report(),
        flow.accepted_port_exchanges().map(<[_]>::to_vec),
        flow.accepted_molecular_port_exchanges().map(<[_]>::to_vec),
        flow.accepted_branch_transport_reports().map(<[_]>::to_vec),
    )
}

#[test]
fn refused_or_late_cancelled_compositions_preserve_both_port_reports_and_reuse_workspace() {
    let g = gas(2, true, 5.0);
    let make = || {
        let mut flow = owner(&g);
        flow.step(Time::new(0.004), DissipatingStepControl::default(), || {
            false
        })
        .unwrap();
        flow
    };
    let mut reference = make();
    let mut polls = 0;
    reference
        .step(Time::new(0.004), DissipatingStepControl::default(), || {
            polls += 1;
            false
        })
        .unwrap();
    assert!(polls > 100);
    let mut flow = make();
    let before = snapshot(&flow);
    let operator_report = flow.molecular.last_report();
    let operator_ports = flow.molecular.accepted_port_exchanges().unwrap().to_vec();
    let operator_branches = flow.molecular.accepted_branch_reports().unwrap().to_vec();
    let mut visits = 0;
    assert_eq!(
        flow.step(Time::new(0.004), DissipatingStepControl::default(), || {
            visits += 1;
            visits >= polls
        }),
        Err(DissipatingNetworkError::Cancelled)
    );
    assert_eq!(snapshot(&flow), before);
    assert_eq!(flow.molecular.last_report(), operator_report);
    assert_eq!(
        flow.molecular.accepted_port_exchanges().unwrap(),
        operator_ports
    );
    assert_eq!(
        flow.molecular.accepted_branch_reports().unwrap(),
        operator_branches
    );
    // An incomplete molecular half cannot count as an accepted whole step.
    let refused = flow.step(
        Time::new(0.004),
        DissipatingStepControl {
            molecular: DissipationControl {
                maximum_steps: 1,
                maximum_step_s: 1e-5,
                ..DissipationControl::default()
            },
            maximum_coupling_retries: 0,
            ..DissipatingStepControl::default()
        },
        || false,
    );
    assert!(matches!(
        refused,
        Err(DissipatingNetworkError::RetryLimit { .. })
    ));
    assert_eq!(snapshot(&flow), before);
    assert_eq!(flow.molecular.last_report(), operator_report);
    flow.step(Time::new(0.004), DissipatingStepControl::default(), || {
        false
    })
    .unwrap();
    let a = fields(&flow.accepted);
    let b = fields(&reference.accepted);
    assert_eq!(flow.time(), reference.time());
    for (q, expected) in
        a.0.iter()
            .flatten()
            .chain(a.1.iter())
            .zip(b.0.iter().flatten().chain(b.1.iter()))
    {
        for k in 0..2 {
            close(q.partial_densities[k], expected.partial_densities[k], 3e-12);
        }
        close(q.momentum, expected.momentum, 3e-12);
        close(q.energy, expected.energy, 3e-11);
    }
}


#[test]
fn prescribed_pressure_reservoirs_participate_in_molecular_network_retries_and_balances() {
    let g = gas(2, true, 5.0);
    let mut flow = owner(&g);
    let supply = g.conserved(&[0.0, 2.0], 0.0, 2.0).unwrap();
    let reservoirs = [[Some(supply.as_cell()), None]];
    let before = inventories(&flow);
    let report = flow
        .step_with_reservoirs(
            Time::new(0.004),
            DissipatingStepControl::default(),
            &reservoirs,
            || false,
        )
        .unwrap();

    assert_eq!(report.elapsed, report.transport.elapsed);
    assert!(report.boundary_change.species_masses[1] > 0.0);
    assert!(flow.branches()[0].cell(0).unwrap().partial_densities[1] > 0.1);
    assert_eq!(report.molecular.first.elapsed.value(), 0.5 * report.elapsed.value());
    assert_eq!(report.molecular.second.elapsed.value(), 0.5 * report.elapsed.value());
    for k in 0..2 {
        close(report.conservation_residual.species_masses[k], 0.0, 6e-12);
    }
    close(report.mass_residual_kg, 0.0, 6e-12);
    close(report.absolute_energy_residual_j, 0.0, 6e-11);

    let after = inventories(&flow);
    for k in 0..2 {
        close(after.0[k] - before.0[k], report.boundary_change.species_masses[k], 6e-12);
    }
    close(
        after.1 - before.1,
        report.boundary_change.absolute_energy(&g).unwrap(),
        6e-11,
    );
}
