//! Finite-port entropy, physical storage and complete-transaction regressions.
use super::*;
use crate::multicomponent::dissipation::MolecularEntropyControl;

mod adaptive;
mod stabilized;

// Independently integrate Cp/R=a+b*T and ideal-mixture partial pressure.
// Each tuple is (R/M, Cv(T=0)/R, b). These sources have p0=1 and s0 datum=0.
fn source_entropy(network: &MixtureNetwork1d, calorics: &[(f64, f64, f64)]) -> f64 {
    let density = |q: Cell<'_>| {
        let t = network.gas().primitive(q, || false).unwrap().temperature;
        q.partial_densities
            .iter()
            .zip(calorics)
            .filter(|(rho, _)| **rho > 0.0)
            .map(|(&rho, &(r, cv, b))| rho * r * (cv * t.ln() + b * t - (rho * r).ln()))
            .sum::<f64>()
    };
    network
        .branches()
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
        + network
            .plenums()
            .iter()
            .map(|plenum| plenum.volume_m3() * density(plenum.state()))
            .sum::<f64>()
}

fn contrasting_pools(
    formation: f64,
    common_absence: bool,
    one_sided_zero: bool,
) -> MixtureNetwork1d {
    let mut models = vec![
        source("A", R, 3.5, 0.0, formation),
        source("B", R, 6.5, 0.0, -2.0 * formation),
    ];
    if common_absence {
        models.push(source("C", R, 4.0, 0.0, 3.0 * formation));
    }
    let gas = source_gas(models);
    let mut cold = if one_sided_zero {
        vec![1.0, 0.0]
    } else {
        vec![0.1, 0.9]
    };
    let mut hot = vec![0.9, 0.1];
    if common_absence {
        cold.push(0.0);
        hot.push(0.0);
    }
    two_pools(
        &gas,
        gas.conserved(&cold, 0.0, 1.0).unwrap(),
        gas.conserved(&hot, 0.0, 10.0).unwrap(),
        0.07,
    )
}

fn diffusion(network: &MixtureNetwork1d) -> NetworkDissipationOperator<ConstantTransport> {
    NetworkDissipationOperator::new(
        network,
        law(network.gas(), 0.0, 0.0, 0.0, 0.01),
        vec![Length::new(PLENUM_DISTANCE); 2],
        || false,
    )
    .unwrap()
}

#[test]
fn g1_source_conjugate_finite_ports_repair_conservative_negative_entropy_mixing() {
    let controls = DissipationControl {
        maximum_step_s: 1e-4,
        maximum_halvings: 3,
        ..Default::default()
    };
    let entropy = MolecularEntropyControl {
        absolute_tolerance_j_k: 1e-14,
        ..Default::default()
    };
    let calorics = [(1.0, 2.5, 0.0), (1.0, 5.5, 0.0)];
    let mut native = contrasting_pools(0.0, false, false);
    let initial_entropy = source_entropy(&native, &calorics);
    let report = diffusion(&native)
        .advance_for(&mut native, Time::new(1e-4), controls, || false)
        .unwrap();
    assert!(source_entropy(&native, &calorics) < initial_entropy - 1e-6);
    close(report.absolute_energy_residual_j, 0.0, 2e-13);

    let mut refused = contrasting_pools(0.0, false, false);
    let before = fields(&refused);
    let error = diffusion(&refused)
        .advance_for_with_entropy(&mut refused, Time::new(1e-4), controls, entropy, || false)
        .unwrap_err();
    assert!(matches!(error, DissipationError::StepFailure { last, .. }
        if matches!(*last, DissipationError::EntropyDecrease { change_j_k, .. } if change_j_k < -1e-14)));
    assert_eq!(fields(&refused), before);

    let mut reference = None;
    for formation in [0.0, 1e18] {
        let mut network = contrasting_pools(formation, false, false);
        let initial_entropy = source_entropy(&network, &calorics);
        let (report, measured) = diffusion(&network)
            .advance_entropy_consistent_for(
                &mut network,
                Time::new(1e-4),
                controls,
                entropy,
                || false,
            )
            .unwrap();
        assert_eq!(report.accepted_steps, 1);
        assert_eq!(report.rejected_steps, 0);
        assert_eq!(report.boundary_change.species_masses[..2], [0.0; 2]);
        assert_eq!(report.boundary_change.energy, 0.0);
        assert!(measured.entropy_change_j_k > 0.0);
        close(
            measured.source_entropy_change_j_k,
            source_entropy(&network, &calorics) - initial_entropy,
            3e-13,
        );
        assert!(network.branches()[0].cell(0).unwrap().partial_densities[0] > 0.1);
        assert!(network.plenums()[0].state().partial_densities[0] < 0.9);
        assert_eq!(network.time().value(), 0.0);
        if let Some((state, receipt)) = &reference {
            assert_eq!(&fields(&network), state);
            assert_eq!(&measured, receipt);
        } else {
            close(report.absolute_energy_residual_j, 0.0, 3e-13);
            reference = Some((fields(&network), measured));
        }
    }
}

#[test]
fn g1_variable_cp_molecular_network_retains_all_storage_and_stationary_support_work() {
    let gas = gas(2, true, 0.0);
    let controls = DissipationControl {
        maximum_step_s: 0.005,
        ..Default::default()
    };
    let calorics = [(1.0, 2.5, 0.1), (0.5, 3.3, 0.05)];
    let mut receipts = Vec::new();
    for reversed in [false, true] {
        let mut network = one_port(&gas, reversed);
        let before = network.totals();
        let initial_entropy = source_entropy(&network, &calorics);
        let mut op = NetworkDissipationOperator::new(
            &network,
            law(&gas, 0.02, 0.01, 0.04, 0.015),
            vec![Length::new(PLENUM_DISTANCE)],
            || false,
        )
        .unwrap();
        let (report, entropy) = op
            .advance_entropy_consistent_for(
                &mut network,
                Time::new(0.02),
                controls,
                MolecularEntropyControl::default(),
                || false,
            )
            .unwrap();
        assert_eq!(report.accepted_steps, 4);
        assert_eq!(report.rejected_steps, 0);
        close(
            entropy.source_entropy_change_j_k,
            source_entropy(&network, &calorics) - initial_entropy,
            2e-12,
        );
        assert!(entropy.entropy_change_j_k > 0.0);
        assert!(entropy.minimum_step_entropy_change_j_k > 0.0);
        assert!(entropy.maximum_relative_species_residual < 2e-14);
        let port = op.accepted_port_exchanges().unwrap()[0];
        assert_ne!(port.momentum_to_plenum, 0.0);
        assert_eq!(network.plenums()[0].state().momentum, 0.0);
        for k in 0..2 {
            close(
                report.after.species_masses[k],
                before.species_masses[k],
                3e-14,
            );
            close(
                port.species_masses_to_plenum[k],
                0.07 * (network.plenums()[0].state().partial_densities[k] - [0.2, 0.8][k]),
                3e-14,
            );
        }
        close(report.after.energy, before.energy, 3e-13);
        close(report.absolute_energy_residual_j, 0.0, 3e-13);
        assert_eq!(report.boundary_change.energy, 0.0);
        assert_eq!(network.time().value(), 0.0);
        receipts.push((port, entropy));
    }
    close(
        receipts[0].0.momentum_to_plenum,
        -receipts[1].0.momentum_to_plenum,
        3e-14,
    );
    close(
        receipts[0].0.energy_to_plenum,
        receipts[1].0.energy_to_plenum,
        3e-14,
    );
    close(
        receipts[0].1.entropy_change_j_k,
        receipts[1].1.entropy_change_j_k,
        3e-14,
    );
}

#[test]
fn g3_measured_engineering_lane_preserves_native_accepted_network_trajectory() {
    let gas = gas(2, true, 0.0);
    let mut measured = one_port(&gas, false);
    let mut native = one_port(&gas, false);
    let make = |network: &MixtureNetwork1d| {
        NetworkDissipationOperator::new(
            network,
            law(&gas, 0.02, 0.01, 0.04, 0.015),
            vec![Length::new(PLENUM_DISTANCE)],
            || false,
        )
        .unwrap()
    };
    let mut measured_op = make(&measured);
    let mut native_op = make(&native);
    let controls = DissipationControl {
        maximum_step_s: 0.005,
        ..Default::default()
    };
    let (report, entropy) = measured_op
        .advance_for_with_entropy(
            &mut measured,
            Time::new(0.02),
            controls,
            MolecularEntropyControl::default(),
            || false,
        )
        .unwrap();
    let original = native_op
        .advance_for(&mut native, Time::new(0.02), controls, || false)
        .unwrap();
    assert_eq!(report, original);
    assert_eq!(fields(&measured), fields(&native));
    assert_eq!(
        measured_op.accepted_port_exchanges(),
        native_op.accepted_port_exchanges()
    );
    assert_eq!(
        measured_op.accepted_branch_reports(),
        native_op.accepted_branch_reports()
    );
    assert!(entropy.entropy_change_j_k > 0.0);
}

#[test]
fn g4_whole_network_interval_cannot_multiply_the_substep_entropy_allowance() {
    let controls = DissipationControl {
        maximum_step_s: 1e-4,
        ..Default::default()
    };
    let allowance = MolecularEntropyControl {
        absolute_tolerance_j_k: 5e-6,
        ..Default::default()
    };
    let mut split = contrasting_pools(0.0, false, false);
    let mut split_op = diffusion(&split);
    for _ in 0..2 {
        let (report, entropy) = split_op
            .advance_for_with_entropy(&mut split, Time::new(1e-4), controls, allowance, || false)
            .unwrap();
        assert_eq!(report.accepted_steps, 1);
        assert!(entropy.entropy_change_j_k < 0.0);
        assert!(entropy.entropy_change_j_k >= -allowance.absolute_tolerance_j_k);
    }
    let mut whole = contrasting_pools(0.0, false, false);
    let before = fields(&whole);
    let mut op = diffusion(&whole);
    // Install accepted receipts before the later complete-interval refusal.
    op.advance_for(&mut whole, Time::new(0.0), controls, || false)
        .unwrap();
    let report = op.last_report();
    let ports = op.accepted_port_exchanges().unwrap().to_vec();
    let branches = op.accepted_branch_reports().unwrap().to_vec();
    let error = op
        .advance_for_with_entropy(&mut whole, Time::new(2e-4), controls, allowance, || false)
        .unwrap_err();
    assert!(
        matches!(error, DissipationError::EntropyDecrease { change_j_k, .. }
        if change_j_k < -allowance.absolute_tolerance_j_k)
    );
    assert_eq!(fields(&whole), before);
    assert_eq!(whole.time().value(), 0.0);
    assert_eq!(op.last_report(), report);
    assert_eq!(op.accepted_port_exchanges().unwrap(), ports);
    assert_eq!(op.accepted_branch_reports().unwrap(), branches);
}

#[test]
fn g0_common_absence_zero_interval_and_one_sided_port_refusal_keep_physical_inventory() {
    let controls = DissipationControl::default();
    let entropy = MolecularEntropyControl::default();
    let mut network = contrasting_pools(0.0, true, false);
    let mut op = diffusion(&network);
    let before = fields(&network);
    let (report, zero) = op
        .advance_entropy_consistent_for(&mut network, Time::new(0.0), controls, entropy, || false)
        .unwrap();
    assert_eq!(report.accepted_steps, 0);
    assert_eq!(zero.entropy_change_j_k, 0.0);
    assert_eq!(zero.source_entropy_change_j_k, 0.0);
    assert_eq!(fields(&network), before);
    let (report, entropy) = op
        .advance_entropy_consistent_for(&mut network, Time::new(1e-4), controls, entropy, || false)
        .unwrap();
    assert_eq!(report.after.species_masses[2], 0.0);
    assert!(entropy.entropy_change_j_k > 0.0);
    for branch in network.branches() {
        for i in 0..branch.cell_count() {
            assert_eq!(branch.cell(i).unwrap().partial_densities[2], 0.0);
        }
    }
    assert_eq!(network.plenums()[0].state().partial_densities[2], 0.0);
    let mut refused = contrasting_pools(0.0, false, true);
    let before = fields(&refused);
    let mut op = diffusion(&refused);
    let error = op
        .advance_entropy_consistent_for(
            &mut refused,
            Time::new(1e-4),
            controls,
            MolecularEntropyControl::default(),
            || false,
        )
        .unwrap_err();
    assert!(matches!(error, DissipationError::Invalid(_)));
    assert_eq!(fields(&refused), before);
    assert_eq!(refused.time().value(), 0.0);
    assert!(op.last_report().is_none());
}

#[test]
fn g4_last_poll_cancellation_keeps_every_network_receipt_and_reuses_the_prepared_owner() {
    let controls = DissipationControl {
        maximum_step_s: 1e-4,
        ..Default::default()
    };
    let entropy = MolecularEntropyControl::default();
    let mut probe = contrasting_pools(0.0, false, false);
    let mut probe_op = diffusion(&probe);
    let mut polls = 0;
    let expected = probe_op
        .advance_entropy_consistent_for(&mut probe, Time::new(2e-4), controls, entropy, || {
            polls += 1;
            false
        })
        .unwrap();
    let mut network = contrasting_pools(0.0, false, false);
    let mut op = diffusion(&network);
    op.advance_for(&mut network, Time::new(0.0), controls, || false)
        .unwrap();
    let before = fields(&network);
    let report = op.last_report();
    let ports = op.accepted_port_exchanges().unwrap().to_vec();
    let branches = op.accepted_branch_reports().unwrap().to_vec();
    let mut visits = 0;
    let result =
        op.advance_entropy_consistent_for(&mut network, Time::new(2e-4), controls, entropy, || {
            visits += 1;
            visits == polls
        });
    assert_eq!(result, Err(DissipationError::Cancelled));
    assert_eq!(fields(&network), before);
    assert_eq!(network.time().value(), 0.0);
    assert_eq!(op.last_report(), report);
    assert_eq!(op.accepted_port_exchanges().unwrap(), ports);
    assert_eq!(op.accepted_branch_reports().unwrap(), branches);
    let recovered = op
        .advance_entropy_consistent_for(&mut network, Time::new(2e-4), controls, entropy, || false)
        .unwrap();
    assert_eq!(fields(&network), fields(&probe));
    assert_eq!(recovered.0.elapsed, expected.0.elapsed);
    assert_eq!(recovered.0.accepted_steps, expected.0.accepted_steps);
    assert_eq!(recovered.0.rejected_steps, expected.0.rejected_steps);
    assert_eq!(
        recovered.0.maximum_gradient_residual_per_m,
        expected.0.maximum_gradient_residual_per_m
    );
    assert_eq!(
        recovered.0.maximum_mass_frame_residual_kg_m2_s,
        expected.0.maximum_mass_frame_residual_kg_m2_s
    );
    assert_eq!(recovered.0.before, expected.0.before);
    assert_eq!(recovered.0.after, expected.0.after);
    assert_eq!(recovered.0.boundary_change, expected.0.boundary_change);
    assert_eq!(
        recovered.0.conservation_residual,
        expected.0.conservation_residual
    );
    assert_eq!(
        recovered.0.absolute_energy_residual_j,
        expected.0.absolute_energy_residual_j
    );
    // Cancelled trials retain only private EOS temperature guesses. Reusing a
    // guess can choose an adjacent within-tolerance root; the fields above are
    // still bit-identical. Compare temperature-derived diagnostics at roundoff,
    // while cancellation itself must retain every accepted receipt exactly.
    close(
        recovered.0.maximum_parabolic_number,
        expected.0.maximum_parabolic_number,
        8.0 * f64::EPSILON * expected.0.maximum_parabolic_number,
    );
    let entropy_roundoff = 64.0
        * f64::EPSILON
        * (2.5 * expected.0.before.species_masses[0] + 5.5 * expected.0.before.species_masses[1]);
    close(
        recovered.1.entropy_change_j_k,
        expected.1.entropy_change_j_k,
        entropy_roundoff,
    );
    close(
        recovered.1.source_entropy_change_j_k,
        expected.1.source_entropy_change_j_k,
        entropy_roundoff,
    );
    close(
        recovered.1.minimum_step_entropy_change_j_k,
        expected.1.minimum_step_entropy_change_j_k,
        entropy_roundoff,
    );
    assert_eq!(
        recovered.1.reference_entropy_residual_j_k,
        expected.1.reference_entropy_residual_j_k
    );
    assert_eq!(
        recovered.1.maximum_relative_species_residual,
        expected.1.maximum_relative_species_residual
    );
    assert_eq!(
        recovered.1.rejected_entropy_steps,
        expected.1.rejected_entropy_steps
    );
}
