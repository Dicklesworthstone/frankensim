//! Original polynomial-Cp state entropy and stationary reservoir-effort oracles.
use super::*;
use crate::multicomponent::Totals;

const CALORICS: [[f64; 3]; 2] = [[3.5, 0.2, 0.0], [4.0, 0.1, 1.5]];
const MINIMUM_TEMPERATURE: f64 = 0.001;

// Each source has M=R, p0=1 and Cp/R=a+b*T. This quadratic inverse uses only
// the authored coefficients and original referenced conservative fields.
fn temperature(q: Cell<'_>, calorics: &[[f64; 3]]) -> f64 {
    let density: f64 = q.partial_densities.iter().sum();
    let a: f64 = q
        .partial_densities
        .iter()
        .zip(calorics)
        .map(|(&rho, &[a, _, _])| rho * (a - 1.0))
        .sum();
    let b: f64 = q
        .partial_densities
        .iter()
        .zip(calorics)
        .map(|(&rho, &[_, b, _])| 0.5 * rho * b)
        .sum();
    let target = q.energy - 0.5 * q.momentum * q.momentum / density
        + a * MINIMUM_TEMPERATURE
        + b * MINIMUM_TEMPERATURE.powi(2);
    2.0 * target / (a + (a * a + 4.0 * b * target).sqrt())
}

fn independent_entropy(flow: &ReactingNetwork1d, calorics: &[[f64; 3]]) -> f64 {
    let density = |q: Cell<'_>| {
        let t = temperature(q, calorics);
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
                .map(|(i, &volume)| volume * density(branch.cell(i).unwrap()))
                .sum::<f64>()
        })
        .sum::<f64>()
        + flow
            .plenums()
            .iter()
            .map(|plenum| plenum.volume_m3() * density(plenum.state()))
            .sum::<f64>()
}

fn reservoir_entropy(source: Cell<'_>, transfer: &Totals, calorics: &[[f64; 3]]) -> f64 {
    assert_eq!(source.momentum, 0.0);
    let t = temperature(source, calorics);
    let mut input = transfer.energy / t;
    for (k, &[a, b, c]) in calorics.iter().enumerate() {
        if transfer.species_masses[k] == 0.0 {
            continue;
        }
        let rho = source.partial_densities[k];
        assert!(rho > 0.0);
        let u =
            (a - 1.0) * (t - MINIMUM_TEMPERATURE) + 0.5 * b * (t * t - MINIMUM_TEMPERATURE.powi(2));
        let s = (a - 1.0) * t.ln() + b * t + c - rho.ln();
        input -= (u / t + 1.0 - s) * transfer.species_masses[k];
    }
    input
}

fn boundary_entropy(
    flow: &ReactingNetwork1d,
    reservoirs: &[[Option<Cell<'_>>; 2]],
    calorics: &[[f64; 3]],
) -> [f64; 2] {
    let mut input = [0.0; 2];
    for (branch, ends) in flow
        .accepted_branch_transport_reports()
        .unwrap()
        .iter()
        .zip(reservoirs)
    {
        for (end, transfer) in [branch.left_boundary_change, branch.right_boundary_change]
            .iter()
            .enumerate()
        {
            if let Some(source) = ends[end] {
                input[end] += reservoir_entropy(source, transfer, calorics);
            }
        }
    }
    input
}

fn reservoirs(flow: &ReactingNetwork1d) -> [Conserved; 2] {
    [
        flow.gas().conserved(&[0.6, 1.2], 0.0, 4.5).unwrap(),
        flow.gas().conserved(&[0.15, 0.15], 0.0, 0.24).unwrap(),
    ]
}

fn ends(states: &[Conserved; 2]) -> [[Option<Cell<'_>>; 2]; 2] {
    [
        [Some(states[0].as_cell()), None],
        [None, Some(states[1].as_cell())],
    ]
}

#[test]
fn g1_reacting_reservoir_network_retains_original_chemical_and_boundary_entropy() {
    let (mut flow, _) = variable_caloric();
    let (mut native, _) = variable_caloric();
    let states = reservoirs(&flow);
    let sources = ends(&states);
    let before_s = independent_entropy(&flow, &CALORICS);
    let stripped = [[3.5, 0.2, 0.0], [4.0, 0.1, 0.0]];
    let before_stripped = independent_entropy(&flow, &stripped);
    let report = flow
        .step_with_reservoirs_and_entropy(
            Time::new(0.002),
            controls(),
            &sources,
            MixtureEntropyControl::default(),
            || false,
        )
        .unwrap();
    let ordinary = native
        .step_with_reservoirs(Time::new(0.002), controls(), &sources, || false)
        .unwrap();
    assert_eq!(report.coupled, ordinary);
    assert_eq!(fields(&flow.accepted), fields(&native.accepted));
    let input = boundary_entropy(&flow, &sources, &CALORICS);
    let actual = independent_entropy(&flow, &CALORICS) - before_s;
    for end in 0..2 {
        close(
            report.entropy.boundary_entropy_input_j_k[end],
            input[end],
            5e-13,
        );
    }
    close(report.entropy.source_entropy_change_j_k, actual, 5e-12);
    close(
        report.entropy.source_entropy_generation_j_k,
        actual - input.iter().sum::<f64>(),
        5e-12,
    );
    assert!(report.entropy.entropy_generation_j_k > 0.0);
    assert!(report.transport_entropy.entropy_generation_j_k > 0.0);
    assert!(report.entropy.maximum_relative_species_residual < 2e-11);
    assert!(report.coupled.chemistry.species_change_kg[1] > 1e-4);
    // Real chemical entropy datums survive boundary cancellation. Only the
    // externally transported B datum cancels against its reservoir supply.
    let stripped_input = boundary_entropy(&flow, &sources, &stripped);
    let stripped_generation = independent_entropy(&flow, &stripped)
        - before_stripped
        - stripped_input.iter().sum::<f64>();
    close(
        actual - input.iter().sum::<f64>() - stripped_generation,
        1.5 * report.coupled.chemistry.species_change_kg[1],
        7e-13,
    );
    assert!((report.entropy.entropy_generation_j_k - stripped_generation).abs() > 1e-4);
    let branches = flow.accepted_branch_transport_reports().unwrap();
    assert!(branches[0].left_boundary_change.mass() > 0.0);
    assert!(branches[1].right_boundary_change.mass() < 0.0);
    assert!(
        flow.accepted_port_exchanges()
            .unwrap()
            .iter()
            .any(|port| port.energy_to_plenum != 0.0)
    );
    assert_eq!(report.coupled.boundary_work_on_gas_j, 0.0);
    close(report.coupled.mass_residual_kg, 0.0, 5e-12);
    close(report.coupled.absolute_energy_residual_j, 0.0, 5e-11);
}

#[test]
fn g1_reacting_outflow_can_lower_state_entropy_with_positive_generation() {
    let mechanism = mechanism(0.0, 1e-8, 0.0);
    let gas = MixtureGas::new(
        mechanism.thermodynamics().clone(),
        CaloricDomain::new(Temperature::new(0.001), Temperature::new(50.0)),
        || false,
    )
    .unwrap();
    // At rho_A=1 and T=2, s_A=2.5*ln(2)>0. Removing that material
    // lowers the network's entropy while the receiving reservoir can generate
    // entropy. The T=1 source datum has s_A=0 and cannot prove this distinction.
    let hot = gas.conserved(&[1.0, 0.0], 0.0, 2.0).unwrap();
    let branches = (0..2)
        .map(|_| {
            MixtureEuler1d::new(
                gas.clone(),
                Length::new(0.25),
                0.1,
                vec![hot.clone(); 4],
                Boundaries::Ends {
                    left: Boundary::Reflecting,
                    right: Boundary::Reflecting,
                },
                || false,
            )
            .unwrap()
        })
        .collect();
    let network = MixtureNetwork1d::new(
        branches,
        vec![Plenum::new(0.03, hot).unwrap()],
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
    .unwrap();
    let mut flow = ReactingNetwork1d::new(network, mechanism, chemistry(), || false).unwrap();
    let sink = flow.gas().conserved(&[0.3, 0.01], 0.0, 0.31).unwrap();
    let sources = [[Some(sink.as_cell()), None], [None, Some(sink.as_cell())]];
    let calorics = [[3.5, 0.0, 0.0]; 2];
    let before_s = independent_entropy(&flow, &calorics);
    let report = flow
        .step_with_reservoirs_and_entropy(
            Time::new(0.002),
            controls(),
            &sources,
            MixtureEntropyControl::default(),
            || false,
        )
        .unwrap();
    assert!(report.coupled.boundary_change.mass() < 0.0);
    assert!(report.entropy.source_entropy_change_j_k < 0.0);
    assert!(report.entropy.entropy_generation_j_k > 0.0);
    let input = boundary_entropy(&flow, &sources, &calorics);
    for end in 0..2 {
        close(
            report.entropy.boundary_entropy_input_j_k[end],
            input[end],
            5e-13,
        );
    }
    close(
        report.entropy.source_entropy_change_j_k,
        independent_entropy(&flow, &calorics) - before_s,
        5e-12,
    );
    assert!(report.entropy.maximum_relative_species_residual < 2e-11);
}

#[test]
fn g1_open_heat_release_retries_and_absent_reservoirs_keep_native_physics() {
    let mechanism = mechanism(40.0, 500.0, 0.0);
    for active in [false, true] {
        let mut flow = reacting(&mechanism, 4, false);
        let mut native = reacting(&mechanism, 4, false);
        let source = flow.gas().conserved(&[0.9, 0.1], 0.0, 1.0).unwrap();
        let trace = active.then_some(source.as_cell());
        let sources = [[trace, None], [None, trace]];
        let report = flow
            .step_with_reservoirs_and_entropy(
                Time::new(0.2),
                controls(),
                &sources,
                MixtureEntropyControl::default(),
                || false,
            )
            .unwrap();
        let ordinary = native
            .step_with_reservoirs(Time::new(0.2), controls(), &sources, || false)
            .unwrap();
        assert_eq!(report.coupled, ordinary);
        assert_eq!(fields(&flow.accepted), fields(&native.accepted));
        assert!(report.coupled.coupling_retries > 0);
        assert_eq!(flow.time(), report.coupled.elapsed);
        assert_eq!(report.coupled.elapsed, report.coupled.transport.elapsed);
        assert!(report.entropy.entropy_generation_j_k > 0.0);
        if !active {
            assert_eq!(report.entropy.boundary_entropy_input_j_k, [0.0; 2]);
        }
        close(report.coupled.absolute_energy_residual_j, 0.0, 5e-11);
    }
}

#[test]
fn g4_open_reacting_final_cancellation_and_zero_activity_preserve_accepted_receipts() {
    let (mut probe, _) = variable_caloric();
    let states = reservoirs(&probe);
    let sources = ends(&states);
    let advance = |flow: &mut ReactingNetwork1d| {
        flow.step_with_reservoirs_and_entropy(
            Time::new(0.001),
            controls(),
            &sources,
            MixtureEntropyControl::default(),
            || false,
        )
        .unwrap()
    };
    advance(&mut probe);
    let mut polls = 0;
    let expected = probe
        .step_with_reservoirs_and_entropy(
            Time::new(0.002),
            controls(),
            &sources,
            MixtureEntropyControl::default(),
            || {
                polls += 1;
                false
            },
        )
        .unwrap();
    let (mut flow, _) = variable_caloric();
    advance(&mut flow);
    let before = fields(&flow.accepted);
    let clock = flow.time();
    let receipt = flow.last_report();
    let ports = flow.accepted_port_exchanges().unwrap().to_vec();
    let branches = flow.accepted_branch_transport_reports().unwrap().to_vec();
    let mut visits = 0;
    let result = flow.step_with_reservoirs_and_entropy(
        Time::new(0.002),
        controls(),
        &sources,
        MixtureEntropyControl::default(),
        || {
            visits += 1;
            visits == polls
        },
    );
    assert_eq!(result, Err(ReactiveNetworkError::Cancelled));
    assert_eq!(fields(&flow.accepted), before);
    assert_eq!(flow.time(), clock);
    assert_eq!(flow.last_report(), receipt);
    assert_eq!(flow.accepted_port_exchanges().unwrap(), ports);
    assert_eq!(flow.accepted_branch_transport_reports().unwrap(), branches);
    let sink = flow.gas().conserved(&[0.1, 0.0], 0.0, 0.1).unwrap();
    let invalid = [[Some(sink.as_cell()), None], [None, Some(sink.as_cell())]];
    assert!(
        flow.step_with_reservoirs_and_entropy(
            Time::new(0.002),
            controls(),
            &invalid,
            MixtureEntropyControl::default(),
            || false,
        )
        .is_err()
    );
    assert_eq!(fields(&flow.accepted), before);
    assert_eq!(flow.time(), clock);
    assert_eq!(flow.last_report(), receipt);
    assert_eq!(flow.accepted_port_exchanges().unwrap(), ports);
    assert_eq!(flow.accepted_branch_transport_reports().unwrap(), branches);
    let recovered = flow
        .step_with_reservoirs_and_entropy(
            Time::new(0.002),
            controls(),
            &sources,
            MixtureEntropyControl::default(),
            || false,
        )
        .unwrap();
    close(
        recovered.entropy.source_entropy_change_j_k,
        expected.entropy.source_entropy_change_j_k,
        5e-12,
    );
    close(
        recovered.entropy.entropy_generation_j_k,
        expected.entropy.entropy_generation_j_k,
        5e-12,
    );
    for (actual, expected) in flow.branches().iter().zip(probe.branches()) {
        for i in 0..actual.cell_count() {
            let (a, b) = (actual.cell(i).unwrap(), expected.cell(i).unwrap());
            for k in 0..2 {
                close(a.partial_densities[k], b.partial_densities[k], 5e-12);
            }
            close(a.momentum, b.momentum, 5e-12);
            close(a.energy, b.energy, 5e-11);
        }
    }
    assert_eq!(flow.time(), probe.time());
}
