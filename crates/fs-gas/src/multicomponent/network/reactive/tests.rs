//! Independent coupled chemistry, physical-energy and network-owner checks.
mod molecular;
use super::{ReactingNetwork1d, ReactiveNetworkError, ReactiveNetworkReport};
use crate::chemistry::{IsochoricKineticsSolve, advance_adiabatic};
use crate::multicomponent::network::{
    BranchEnd, MixtureNetwork1d, Plenum, PlenumConnection, PortExchange,
};
use crate::multicomponent::reactive::ReactiveStepControl;
use crate::multicomponent::tests::{bank, source};
use crate::multicomponent::{
    Boundaries, Boundary, CaloricDomain, Cell, Conserved, MixtureEuler1d, MixtureGas, StepReport,
};
use fs_qty::{Amount, Length, Temperature, Time, Volume};
use fs_thermochem::kinetics::{
    ArrheniusRateV1, ElementaryReactionV1, IdealGasKineticsV1, ReactionProgressRateQuantityV1,
};
use fs_thermochem::{
    MolarEnergyQuantityV1, ReactionId, SpeciesId, UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R,
};

fn mechanism(heat: f64, rate: f64, reference: f64) -> IdealGasKineticsV1 {
    let thermodynamics = bank(vec![
        source("A", R, 3.5, 0.0, reference),
        source("B", R, 3.5, 0.0, reference - heat),
    ]);
    let reaction = ElementaryReactionV1::new(
        ReactionId::new("network-isomerization").unwrap(),
        vec![(SpeciesId::new("A").unwrap(), 1)],
        vec![(SpeciesId::new("B").unwrap(), 1)],
        ArrheniusRateV1::new(
            Temperature::new(1.0),
            ReactionProgressRateQuantityV1::new(rate),
            0.0,
            MolarEnergyQuantityV1::new(0.0),
        )
        .unwrap(),
    )
    .unwrap();
    IdealGasKineticsV1::new(
        thermodynamics,
        "reacting-network-isomers",
        1,
        vec![reaction],
        1e-12,
    )
    .unwrap()
}
fn chemistry() -> IsochoricKineticsSolve {
    IsochoricKineticsSolve {
        minimum_temperature_k: 0.001,
        maximum_temperature_k: 50.0,
        energy_tolerance_j: 1e-10,
        relative_tolerance: 2e-7,
        absolute_amount_tolerance_mol: 1e-14,
        entropy_tolerance_j_k: 1e-11,
        initial_step_s: 1e-4,
        maximum_step_s: 1.0,
        maximum_steps: 30_000,
        ..IsochoricKineticsSolve::default()
    }
}
fn network(
    mechanism: &IdealGasKineticsV1,
    n: usize,
    nonuniform: bool,
    plenum_temperature: f64,
) -> MixtureNetwork1d {
    let gas = MixtureGas::new(
        mechanism.thermodynamics().clone(),
        CaloricDomain::new(Temperature::new(0.001), Temperature::new(50.0)),
        || false,
    )
    .unwrap();
    let walls = Boundaries::Ends {
        left: Boundary::Reflecting,
        right: Boundary::Reflecting,
    };
    let branches = (0..2)
        .map(|branch| {
            let states = (0..n)
                .map(|i| {
                    let (densities, pressure) = if nonuniform {
                        if branch == 0 {
                            ([0.9, 0.1], if i * 2 >= n { 2.0 } else { 1.0 })
                        } else {
                            ([0.2, 0.8], 1.0)
                        }
                    } else {
                        ([1.0, 0.0], 1.0)
                    };
                    gas.conserved(&densities, 0.0, pressure).unwrap()
                })
                .collect();
            MixtureEuler1d::new(
                gas.clone(),
                Length::new(1.0 / n as f64),
                0.1,
                states,
                walls,
                || false,
            )
            .unwrap()
        })
        .collect();
    let densities = if nonuniform { [0.65, 0.35] } else { [1.0, 0.0] };
    MixtureNetwork1d::new(
        branches,
        vec![
            Plenum::new(
                0.03,
                gas.conserved(&densities, 0.0, plenum_temperature).unwrap(),
            )
            .unwrap(),
        ],
        vec![
            PlenumConnection {
                branch: 0,
                end: BranchEnd::Right,
                plenum: 0,
            },
            PlenumConnection {
                branch: 1,
                end: BranchEnd::Left,
                plenum: 0,
            },
        ],
        || false,
    )
    .unwrap()
}
fn reacting(mechanism: &IdealGasKineticsV1, n: usize, nonuniform: bool) -> ReactingNetwork1d {
    ReactingNetwork1d::new(
        network(mechanism, n, nonuniform, if nonuniform { 1.4 } else { 1.0 }),
        mechanism.clone(),
        chemistry(),
        || false,
    )
    .unwrap()
}
fn owned(q: Cell<'_>) -> Conserved {
    Conserved {
        partial_densities: q.partial_densities.to_vec(),
        momentum: q.momentum,
        energy: q.energy,
    }
}
fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual={actual:.15e}, expected={expected:.15e}, tolerance={tolerance:.3e}"
    );
}

#[test]
fn connected_homogeneous_reaction_matches_actual_duration_oracles_including_hot_restart() {
    for (heat, rate) in [(0.0, 2.0), (40.0, 500.0)] {
        let mechanism = mechanism(heat, rate, 0.0);
        let mut flow = reacting(&mechanism, 4, false);
        let initial = flow.totals();
        let report = flow
            .step(Time::new(0.2), ReactiveStepControl::default(), || false)
            .unwrap();
        assert_eq!(report.elapsed, report.transport.elapsed);
        assert_eq!(flow.time(), report.elapsed);
        assert_eq!(
            report.chemistry.cell_solves, 18,
            "both halves must react all eight branch cells and the finite plenum"
        );
        assert!(report.chemistry.entropy_change_j_k > 0.0);
        close(report.mass_residual_kg, 0.0, 3e-13);
        close(report.absolute_energy_residual_j, 0.0, 4e-12);
        close(
            report.after.absolute_energy(flow.gas()).unwrap(),
            initial.absolute_energy(flow.gas()).unwrap(),
            4e-12,
        );
        if heat > 0.0 {
            assert!(
                report.coupling_retries > 0,
                "heat release must force a complete acoustic-duration restart"
            );
            assert!(report.chemistry.referenced_energy_change_j > 0.0);
        }
        let check = |q: Cell<'_>, volume: f64| {
            let p = flow.gas().primitive(q, || false).unwrap();
            close(p.velocity, 0.0, 2e-8);
            close(q.density(), 1.0, 2e-8);
            if heat == 0.0 {
                // A <=> B has equal forward/reverse rate k for equal source
                // free energies: Y_A=1/2+1/2 exp(-2kt), independent of volume.
                let exact_a = 0.5 + 0.5 * (-2.0 * rate * report.elapsed.value()).exp();
                close(q.partial_densities[0], exact_a, 3e-5);
                close(q.partial_densities[1], 1.0 - exact_a, 3e-5);
                close(p.temperature, 1.0, 2e-8);
            } else {
                // An independently owned UV reactor runs the whole accepted
                // duration, not either rejected trial duration or a half step.
                let reference = advance_adiabatic(
                    &mechanism,
                    &[(SpeciesId::new("A").unwrap(), Amount::new(volume / R))],
                    Volume::new(volume),
                    Temperature::new(1.0),
                    report.elapsed,
                    chemistry(),
                    || false,
                )
                .unwrap();
                close(p.temperature, reference.temperature().value(), 3e-4);
                close(
                    q.partial_densities[1],
                    reference.amounts()[1].value() * R / volume,
                    3e-5,
                );
                assert!(p.temperature > 5.0);
            }
        };
        for branch in flow.branches() {
            assert_eq!(branch.time(), flow.time());
            for (i, &volume) in branch.mesh().volumes().iter().enumerate() {
                check(branch.cell(i).unwrap(), volume);
            }
        }
        for plenum in flow.plenums() {
            check(plenum.state(), plenum.volume_m3());
        }
    }
}

#[test]
fn reacting_plenum_and_branches_close_separate_chemical_and_actual_port_energy_budgets() {
    let mechanism = mechanism(8.0, 20.0, 3.0);
    let mut flow = reacting(&mechanism, 6, true);
    let initial = flow.totals();
    let mut plenum_production = 0.0;
    let mut branch_production = [0.0; 2];
    let (mut inward, mut outward) = (false, false);
    for _ in 0..3 {
        let old_plenum = owned(flow.plenums()[0].state());
        let old_branches: Vec<_> = flow.branches().iter().map(MixtureEuler1d::totals).collect();
        let report = flow
            .step(Time::new(0.01), ReactiveStepControl::default(), || false)
            .unwrap();
        assert_eq!(report.chemistry.cell_solves, 26);
        assert!(report.chemistry.species_change_kg[0] < 0.0);
        assert!(report.chemistry.species_change_kg[1] > 0.0);
        close(
            report.chemistry.species_change_kg[0] + report.chemistry.species_change_kg[1],
            0.0,
            3e-13,
        );
        for k in 0..2 {
            close(
                report.after.species_masses[k] - report.before.species_masses[k],
                report.chemistry.species_change_kg[k],
                4e-12,
            );
            close(report.conservation_residual.species_masses[k], 0.0, 4e-12);
            assert_eq!(report.transport.boundary_change.species_masses[k], 0.0);
        }
        close(report.mass_residual_kg, 0.0, 4e-12);
        close(report.absolute_energy_residual_j, 0.0, 4e-11);
        let ports = flow.accepted_port_exchanges().unwrap();
        let volume = flow.plenums()[0].volume_m3();
        let q = flow.plenums()[0].state();
        let received_b: f64 = ports.iter().map(|p| p.species_masses_to_plenum[1]).sum();
        plenum_production +=
            volume * (q.partial_densities[1] - old_plenum.partial_densities[1]) - received_b;
        let received_energy: f64 = ports
            .iter()
            .map(|p| p.absolute_energy_to_plenum(flow.gas()).unwrap())
            .sum();
        let actual_energy_change = volume
            * (flow.gas().absolute_energy_density(q).unwrap()
                - flow
                    .gas()
                    .absolute_energy_density(old_plenum.as_cell())
                    .unwrap());
        close(actual_energy_change, received_energy, 4e-11);
        assert_eq!(q.momentum, 0.0);
        for (i, port) in ports.iter().enumerate() {
            inward |= port.mass_to_plenum() > 1e-6;
            outward |= port.mass_to_plenum() < -1e-6;
            let current = flow.branches()[i].totals();
            let actual = current.absolute_energy(flow.gas()).unwrap()
                - old_branches[i].absolute_energy(flow.gas()).unwrap();
            close(
                actual,
                -port.absolute_energy_to_plenum(flow.gas()).unwrap(),
                4e-11,
            );
            branch_production[i] += current.species_masses[1] - old_branches[i].species_masses[1]
                + port.species_masses_to_plenum[1];
        }
    }
    assert!(
        inward && outward,
        "the fixture must exchange real gas both into and out of its reacting plenum"
    );
    assert!(
        plenum_production > 1e-5,
        "the finite plenum did not react: {plenum_production}"
    );
    assert!(
        branch_production.iter().all(|v| *v > 1e-5),
        "one branch did not react: {branch_production:?}"
    );
    close(flow.totals().mass(), initial.mass(), 4e-12);
    close(
        flow.totals().absolute_energy(flow.gas()).unwrap(),
        initial.absolute_energy(flow.gas()).unwrap(),
        8e-11,
    );
    assert!(
        flow.gas()
            .primitive(flow.plenums()[0].state(), || false)
            .unwrap()
            .temperature
            > 1.4
    );
}

type Snapshot = (
    Vec<Vec<Conserved>>,
    Vec<Conserved>,
    Vec<f64>,
    f64,
    Option<ReactiveNetworkReport>,
    Option<Vec<PortExchange>>,
    Option<Vec<StepReport>>,
);
fn snapshot(flow: &ReactingNetwork1d) -> Snapshot {
    (
        flow.branches()
            .iter()
            .map(|b| {
                (0..b.cell_count())
                    .map(|i| owned(b.cell(i).unwrap()))
                    .collect()
            })
            .collect(),
        flow.plenums().iter().map(|p| owned(p.state())).collect(),
        flow.branches().iter().map(|b| b.time().value()).collect(),
        flow.time().value(),
        flow.last_report(),
        flow.accepted_port_exchanges().map(<[_]>::to_vec),
        flow.accepted_branch_transport_reports().map(<[_]>::to_vec),
    )
}

#[test]
fn late_network_chemistry_cancellation_preserves_all_owners_reports_and_workspace_reuse() {
    let mechanism = mechanism(8.0, 20.0, 0.0);
    let mut probe = reacting(&mechanism, 3, true);
    probe
        .step(Time::new(0.002), ReactiveStepControl::default(), || false)
        .unwrap();
    let mut polls = 0;
    let expected = probe
        .step(Time::new(0.002), ReactiveStepControl::default(), || {
            polls += 1;
            false
        })
        .unwrap();
    assert!(polls > 100);
    for stop in [polls / 2, polls - 2, polls] {
        let mut flow = reacting(&mechanism, 3, true);
        flow.step(Time::new(0.002), ReactiveStepControl::default(), || false)
            .unwrap();
        let original = snapshot(&flow);
        let mut visits = 0;
        let result = flow.step(Time::new(0.002), ReactiveStepControl::default(), || {
            visits += 1;
            visits >= stop
        });
        assert!(
            matches!(result, Err(ReactiveNetworkError::Cancelled)),
            "stop={stop}, result={result:?}"
        );
        assert_eq!(snapshot(&flow), original);
        let recovered = flow
            .step(Time::new(0.002), ReactiveStepControl::default(), || false)
            .unwrap();
        assert_eq!(recovered.elapsed, expected.elapsed);
        assert_eq!(flow.time(), probe.time());
        for (a, b) in flow.branches().iter().zip(probe.branches()) {
            for i in 0..a.cell_count() {
                let qa = a.cell(i).unwrap();
                let qb = b.cell(i).unwrap();
                for k in 0..2 {
                    close(qa.partial_densities[k], qb.partial_densities[k], 2e-11);
                }
                close(qa.momentum, qb.momentum, 2e-11);
                close(qa.energy, qb.energy, 2e-11);
            }
        }
        let a = flow.plenums()[0].state();
        let b = probe.plenums()[0].state();
        for k in 0..2 {
            close(a.partial_densities[k], b.partial_densities[k], 2e-11);
        }
        close(a.energy, b.energy, 2e-11);
    }
}

#[test]
fn reacting_network_admission_checks_original_source_and_every_plenum_temperature() {
    let original = mechanism(8.0, 20.0, 0.0);
    let other = mechanism(8.0, 20.0, 1.0);
    assert!(matches!(
        ReactingNetwork1d::new(
            network(&original, 3, false, 1.0),
            other,
            chemistry(),
            || false
        ),
        Err(ReactiveNetworkError::Invalid(_))
    ));
    let mut wider = chemistry();
    wider.maximum_temperature_k = 60.0;
    assert!(
        ReactingNetwork1d::new(
            network(&original, 3, false, 1.0),
            original.clone(),
            wider,
            || false
        )
        .is_err()
    );
    let mut narrow = chemistry();
    narrow.maximum_temperature_k = 2.0;
    // Every branch is at T=1, inside the declared chemistry interval; only
    // the finite plenum is at T=3. Admission must inspect that owner too.
    assert!(
        ReactingNetwork1d::new(
            network(&original, 3, false, 3.0),
            original.clone(),
            narrow,
            || false
        )
        .is_err()
    );
    assert!(matches!(
        ReactingNetwork1d::new(
            network(&original, 3, false, 1.0),
            original,
            chemistry(),
            || true
        ),
        Err(ReactiveNetworkError::Cancelled)
    ));
}
