//! Independent chemical/source entropy and whole-network publication checks.
mod exterior;
use super::*;
use crate::multicomponent::dissipation::{ConstantTransport, DissipationControl};
use crate::multicomponent::network::dissipation::NetworkDissipationOperator;
use crate::multicomponent::network::dissipation::tests::fields;
use crate::multicomponent::{MixtureEntropyControl, Reconstruction, StepControl};
use fs_thermochem::{Nasa9RegionV1, Nasa9StandardStateModelV1};

mod stabilized;

fn controls() -> ReactiveStepControl {
    ReactiveStepControl {
        transport: StepControl {
            reconstruction: Reconstruction::FirstOrder,
            ..Default::default()
        },
        ..Default::default()
    }
}

// All these sources have M=R, p0=1 and Cp/R=a+b*T, s0/R=a*ln(T)+b*T+c.
// The c_i dm_i contribution is real chemical entropy when species react.
fn source_entropy(flow: &ReactingNetwork1d, calorics: &[[f64; 3]]) -> f64 {
    let density = |q: Cell<'_>| {
        let t = flow.gas().primitive(q, || false).unwrap().temperature;
        q.partial_densities
            .iter()
            .zip(calorics)
            .filter(|(rho, _)| **rho > 0.0)
            .map(|(&rho, &[a, b, c])| rho * ((a - 1.0) * t.ln() + b * t + c - rho.ln()))
            .sum::<f64>()
    };
    flow.branches()
        .iter()
        .map(|branch| {
            branch
                .mesh()
                .volumes()
                .iter()
                .enumerate()
                .map(|(i, &v)| v * density(branch.cell(i).unwrap()))
                .sum::<f64>()
        })
        .sum::<f64>()
        + flow
            .plenums()
            .iter()
            .map(|plenum| plenum.volume_m3() * density(plenum.state()))
            .sum::<f64>()
}

fn entropy_datum(model: Nasa9StandardStateModelV1, value: f64) -> Nasa9StandardStateModelV1 {
    let mut card = model.regions()[0].card().clone();
    card.parameters.get_mut("a8").unwrap().value = value;
    Nasa9StandardStateModelV1::new(
        model.species().clone(),
        model.molar_mass(),
        model.convention().clone(),
        vec![Nasa9RegionV1::from_card(card).unwrap()],
    )
    .unwrap()
}

fn variable_caloric() -> (
    ReactingNetwork1d,
    NetworkDissipationOperator<ConstantTransport>,
) {
    let thermodynamics = bank(vec![
        source("A", R, 3.5, 0.2, 0.0),
        entropy_datum(source("B", R, 4.0, 0.1, -4.0), 1.5),
    ]);
    let reaction = ElementaryReactionV1::new(
        ReactionId::new("network-entropy-source-isomers").unwrap(),
        vec![(SpeciesId::new("A").unwrap(), 1)],
        vec![(SpeciesId::new("B").unwrap(), 1)],
        ArrheniusRateV1::new(
            Temperature::new(1.0),
            ReactionProgressRateQuantityV1::new(20.0),
            0.0,
            MolarEnergyQuantityV1::new(0.0),
        )
        .unwrap(),
    )
    .unwrap();
    let mechanism = IdealGasKineticsV1::new(
        thermodynamics,
        "variable-caloric-network-entropy",
        1,
        vec![reaction],
        1e-12,
    )
    .unwrap();
    let network = network(&mechanism, 4, true, 1.4);
    let molecular = super::molecular::operator(&network, 0.04, 0.01, 0.04);
    let flow = ReactingNetwork1d::new(network, mechanism, chemistry(), || false).unwrap();
    (flow, molecular)
}

#[test]
fn g1_reaction_created_species_and_heat_release_keep_the_original_complete_duration() {
    let mechanism = mechanism(40.0, 500.0, 0.0);
    let mut flow = reacting(&mechanism, 4, false);
    let mut native = reacting(&mechanism, 4, false);
    let calorics = [[3.5, 0.0, 0.0]; 2];
    let initial_entropy = source_entropy(&flow, &calorics);
    let initial = flow.totals();
    assert_eq!(initial.species_masses[1], 0.0);
    let report = flow
        .step_with_entropy(
            Time::new(0.2),
            controls(),
            MixtureEntropyControl::default(),
            || false,
        )
        .unwrap();
    let step = report.coupled;
    let original = native.step(Time::new(0.2), controls(), || false).unwrap();
    assert_eq!(step, original);
    assert_eq!(fields(&flow.accepted), fields(&native.accepted));
    assert!(step.coupling_retries > 0);
    assert_eq!(step.elapsed, step.transport.elapsed);
    assert_eq!(flow.time(), step.elapsed);
    assert_eq!(step.chemistry.cell_solves, 18);
    assert!(step.after.species_masses[1] > 0.01);
    assert!(step.dissipation.is_none() && report.dissipation_entropy.is_none());
    for k in 0..2 {
        close(
            step.after.species_masses[k] - initial.species_masses[k],
            step.chemistry.species_change_kg[k],
            3e-12,
        );
    }
    close(
        report.entropy.source_entropy_change_j_k,
        source_entropy(&flow, &calorics) - initial_entropy,
        3e-9,
    );
    close(
        report.entropy.entropy_change_j_k,
        report.entropy.source_entropy_change_j_k,
        3e-10,
    );
    assert!(report.entropy.entropy_change_j_k > 0.01);
    assert!(report.entropy.maximum_relative_species_residual < 1e-11);
    close(report.transport_entropy.entropy_change_j_k, 0.0, 3e-11);
    close(step.mass_residual_kg, 0.0, 3e-13);
    close(step.absolute_energy_residual_j, 0.0, 4e-11);
}

#[test]
fn g1_variable_caloric_network_composition_retains_chemical_species_entropy_datums() {
    let (mut flow, mut molecular) = variable_caloric();
    let calorics = [[3.5, 0.2, 0.0], [4.0, 0.1, 1.5]];
    let stripped = [[3.5, 0.2, 0.0], [4.0, 0.1, 0.0]];
    let initial_entropy = source_entropy(&flow, &calorics);
    let initial_stripped = source_entropy(&flow, &stripped);
    let before = flow.totals();
    let report = flow
        .step_with_dissipation_and_entropy(
            &mut molecular,
            Time::new(0.002),
            controls(),
            DissipationControl::default(),
            MixtureEntropyControl::default(),
            || false,
        )
        .unwrap();
    let step = report.coupled;
    let halves = step.dissipation.unwrap();
    let entropy_halves = report.dissipation_entropy.unwrap();
    assert_eq!(step.elapsed.value(), 0.002);
    assert_eq!(flow.time(), step.elapsed);
    assert_eq!(halves.first.elapsed.value(), 0.5 * step.elapsed.value());
    assert_eq!(halves.second.elapsed.value(), 0.5 * step.elapsed.value());
    assert!(entropy_halves.first.entropy_change_j_k > 0.0);
    assert!(entropy_halves.second.entropy_change_j_k > 0.0);
    assert!(report.transport_entropy.entropy_change_j_k >= -1e-10);
    assert!(step.chemistry.species_change_kg[1] > 1e-4);
    let actual = source_entropy(&flow, &calorics) - initial_entropy;
    let without_datum = source_entropy(&flow, &stripped) - initial_stripped;
    close(
        actual - without_datum,
        1.5 * (step.after.species_masses[1] - before.species_masses[1]),
        2e-13,
    );
    close(report.entropy.source_entropy_change_j_k, actual, 3e-10);
    close(
        report.entropy.entropy_change_j_k,
        actual - report.entropy.reference_entropy_residual_j_k,
        3e-10,
    );
    assert!(
        (report.entropy.entropy_change_j_k - without_datum).abs() > 1e-4,
        "physical chemical entropy datum must not be cancelled as numerical species drift"
    );
    assert!(report.entropy.maximum_relative_species_residual < 1e-11);
    for k in 0..2 {
        close(
            step.after.species_masses[k] - before.species_masses[k],
            step.chemistry.species_change_kg[k],
            3e-12,
        );
    }
    close(step.mass_residual_kg, 0.0, 3e-12);
    close(step.absolute_energy_residual_j, 0.0, 5e-11);
    assert_eq!(step.boundary_change.energy, 0.0);
    assert_eq!(molecular.last_report(), Some(halves.second));
    assert!(
        flow.accepted_molecular_port_exchanges()
            .unwrap()
            .iter()
            .any(|port| port.species_masses_to_plenum[0] != 0.0)
    );
}

#[test]
fn g4_final_cancellation_preserves_reacting_network_and_both_operator_receipts() {
    let entropy = MixtureEntropyControl::default();
    let (mut probe, mut probe_op) = variable_caloric();
    probe
        .step_with_dissipation_and_entropy(
            &mut probe_op,
            Time::new(0.001),
            controls(),
            DissipationControl::default(),
            entropy,
            || false,
        )
        .unwrap();
    let mut polls = 0;
    let expected = probe
        .step_with_dissipation_and_entropy(
            &mut probe_op,
            Time::new(0.002),
            controls(),
            DissipationControl::default(),
            entropy,
            || {
                polls += 1;
                false
            },
        )
        .unwrap();
    let (mut flow, mut molecular) = variable_caloric();
    flow.step_with_dissipation_and_entropy(
        &mut molecular,
        Time::new(0.001),
        controls(),
        DissipationControl::default(),
        entropy,
        || false,
    )
    .unwrap();
    let before = fields(&flow.accepted);
    let clock = flow.time();
    let report = flow.last_report();
    let ports = flow.accepted_port_exchanges().unwrap().to_vec();
    let molecular_ports = flow.accepted_molecular_port_exchanges().unwrap().to_vec();
    let branches = flow.accepted_branch_transport_reports().unwrap().to_vec();
    let operator_report = molecular.last_report();
    let operator_ports = molecular.accepted_port_exchanges().unwrap().to_vec();
    let mut visits = 0;
    let result = flow.step_with_dissipation_and_entropy(
        &mut molecular,
        Time::new(0.002),
        controls(),
        DissipationControl::default(),
        entropy,
        || {
            visits += 1;
            visits == polls
        },
    );
    assert_eq!(result, Err(ReactiveNetworkError::Cancelled));
    assert_eq!(fields(&flow.accepted), before);
    assert_eq!(flow.time(), clock);
    assert_eq!(flow.last_report(), report);
    assert_eq!(flow.accepted_port_exchanges().unwrap(), ports);
    assert_eq!(
        flow.accepted_molecular_port_exchanges().unwrap(),
        molecular_ports
    );
    assert_eq!(flow.accepted_branch_transport_reports().unwrap(), branches);
    assert_eq!(molecular.last_report(), operator_report);
    assert_eq!(molecular.accepted_port_exchanges().unwrap(), operator_ports);
    let recovered = flow
        .step_with_dissipation_and_entropy(
            &mut molecular,
            Time::new(0.002),
            controls(),
            DissipationControl::default(),
            entropy,
            || false,
        )
        .unwrap();
    assert_eq!(flow.time(), probe.time());
    let actual = fields(&flow.accepted);
    let reference = fields(&probe.accepted);
    for (a, b) in actual
        .0
        .iter()
        .flatten()
        .chain(&actual.1)
        .zip(reference.0.iter().flatten().chain(&reference.1))
    {
        for k in 0..2 {
            close(a.partial_densities[k], b.partial_densities[k], 3e-10);
        }
        close(a.momentum, b.momentum, 3e-10);
        close(a.energy, b.energy, 3e-9);
    }
    close(
        recovered.entropy.source_entropy_change_j_k,
        expected.entropy.source_entropy_change_j_k,
        3e-9,
    );
    close(recovered.coupled.absolute_energy_residual_j, 0.0, 5e-11);
}
