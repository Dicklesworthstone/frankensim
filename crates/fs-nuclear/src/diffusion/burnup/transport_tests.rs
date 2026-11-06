//! Synthetic event/transport equations; no fitted reactor or isotope claim.

use super::tests::{budget, close};
use super::*;
use crate::burnup::transport::{
    BurnupTransportLimits, BurnupTransportState, NuclideOutflow, NuclideTransfer, NuclideTransport,
};

fn limits() -> BurnupTransportLimits {
    BurnupTransportLimits {
        burnup: budget().burnup.limits,
        ..BurnupTransportLimits::default()
    }
}

fn movement(rate: f64) -> NuclideTransport {
    NuclideTransport {
        transfers: vec![NuclideTransfer {
            from_compartment: 0,
            to_compartment: 1,
            rates_per_s: vec![rate, 0.0],
        }],
        outflows: vec![],
    }
}

fn fixture(
    prompt_yield: f64,
    neutrons: [f64; 2],
    parents: [f64; 2],
    flows: NuclideTransport,
) -> DiffusionBurnupState {
    let base = super::tests::fixture(2, true, false);
    let mut data = base.neutrons().model().data().clone();
    data.links.clear();
    for reaction in &mut data.reactions {
        reaction.kind = DiffusionReactionKind::Fission {
            prompt_yields: vec![prompt_yield],
            delayed_yields: vec![],
        };
    }
    let model = Arc::new(DiffusionModel::new(data, budget().diffusion.limits).unwrap());
    let inventories = parents
        .into_iter()
        .map(|parent| {
            BurnupInventory::new(Arc::clone(base.inventories()[0].model()), vec![parent, 0.0])
                .unwrap()
        })
        .collect();
    let transported = BurnupTransportState::new(
        inventories,
        model.cell_volumes_m3().to_vec(),
        flows,
        limits(),
        || Ok(()),
    )
    .unwrap();
    DiffusionBurnupState::new_with_transport(
        DiffusionState::new(model, neutrons.to_vec(), vec![]).unwrap(),
        transported,
        base.bindings().to_vec(),
        limits(),
    )
    .unwrap()
}

#[test]
fn g1_moving_fuel_fissions_in_receiver_during_same_interval_and_heat_stays_local() {
    // One prompt neutron per event keeps the synthetic flux exactly constant.
    // A'= -k A, B'=k A - r B, F'=r B, with A(0)=100 and B(0)=0.
    let mut flows = movement(0.6);
    flows.transfers.push(NuclideTransfer {
        from_compartment: 1,
        to_compartment: 0,
        rates_per_s: vec![0.0, 0.5],
    });
    let initial = fixture(1.0, [0.0, 40.0], [100.0, 0.0], flows);
    let dt = 1.2_f64;
    let a = 100.0 * (-0.6 * dt).exp();
    let b = 100.0 * 0.6 * ((-0.4 * dt).exp() - (-0.6 * dt).exp()) / 0.2;
    let fissions = 100.0 - a - b;
    let (next, report) = initial.trial(dt, &[300.0; 2], budget(), || Ok(())).unwrap();
    close(next.inventories()[0].populations()[0], a, 2e-9);
    close(next.inventories()[1].populations()[0], b, 2e-9);
    close(report.burnup[1].fissions, fissions, 2e-9);
    close(report.neutrons.fissions, fissions, 2e-9);
    assert_eq!(report.burnup[0].fissions, 0.0);
    assert_eq!(next.neutrons().neutron_counts()[0], 0.0);
    close(next.neutrons().neutron_counts()[1], 40.0, 1e-10);
    assert!(report.iterations > 1);
    assert!(report.reaction_defect <= budget().reaction_tolerance);
    let movement = report.transport.as_ref().unwrap();
    close(movement.transfer_counts[0][0], 100.0 - a, 2e-9);
    assert!(movement.transfer_counts[1][1] > 0.0);
    assert!(next.inventories()[0].populations()[1] > 0.0);
    close(
        next.inventories()
            .iter()
            .map(|cell| cell.populations()[1])
            .sum(),
        2.0 * fissions,
        2e-9,
    );
    // Integrated residence includes fuel arriving continuously and stable
    // products with zero local removal; an endpoint rectangle is incorrect.
    let ia = (100.0 - a) / 0.6;
    let ib = fissions / 0.4;
    close(movement.integrated_populations_s[0][0], ia, 2e-9);
    close(movement.integrated_populations_s[1][0], ib, 2e-9);
    close(
        movement
            .integrated_populations_s
            .iter()
            .map(|row| row[1])
            .sum(),
        2.0 * (100.0 * dt - ia - ib),
        3e-9,
    );
    assert_eq!(next.inventories()[0].cumulative_released().total_j(), 0.0);
    close(
        next.inventories()[1].cumulative_released().deposited_j,
        8.0 * fissions,
        2e-8,
    );
    close(report.burnup[1].released.radiation_j, 2.0 * fissions, 3e-9);
    close(report.burnup[1].released.escaped_j, fissions, 2e-9);
    assert_eq!(
        initial.transport().unwrap().cumulative_transfer_counts(),
        vec![vec![0.0; 2]; 2]
    );
    assert_eq!(
        next.transport().unwrap().cumulative_transfer_counts(),
        movement.transfer_counts
    );
}

#[test]
fn g3_no_movement_preserves_independent_spatial_solution() {
    let joint = fixture(
        2.0,
        [10.0, 20.0],
        [100.0, 80.0],
        NuclideTransport::default(),
    );
    let independent = DiffusionBurnupState::new(
        joint.neutrons().clone(),
        joint.inventories().to_vec(),
        joint.bindings().to_vec(),
    )
    .unwrap();
    let (a, ar) = joint.trial(0.4, &[300.0; 2], budget(), || Ok(())).unwrap();
    let (b, br) = independent
        .trial(0.4, &[300.0; 2], budget(), || Ok(()))
        .unwrap();
    for cell in 0..2 {
        close(
            a.neutrons().neutron_counts()[cell],
            b.neutrons().neutron_counts()[cell],
            1e-9,
        );
        for (a, b) in a.inventories()[cell]
            .populations()
            .iter()
            .zip(b.inventories()[cell].populations())
        {
            close(*a, *b, 1e-9);
        }
        close(ar.burnup[cell].fissions, br.burnup[cell].fissions, 1e-9);
        close(
            ar.burnup[cell].released.total_j(),
            br.burnup[cell].released.total_j(),
            1e-8,
        );
    }
    assert!(b.transport().is_none());
    assert!(br.transport.is_none());
    assert!(ar.transport.unwrap().transfer_counts.is_empty());
}

#[test]
fn g1_continuous_feed_and_outflow_use_finite_donors_and_true_residence() {
    let mut flows = movement(0.3);
    flows.outflows.push(NuclideOutflow {
        compartment: 0,
        rates_per_s: vec![0.2, 0.0],
    });
    let initial = fixture(1.0, [0.0; 2], [10.0, 0.0], flows);
    let sources = DiffusionBurnupSources {
        neutron_rates_s: vec![0.0; 2],
        precursor_rates_s: vec![],
        nuclide_rates_s: vec![vec![6.0, 0.0], vec![0.0, 0.0]],
    };
    let dt = 2.0_f64;
    let a = 12.0 + (10.0 - 12.0) * (-0.5 * dt).exp();
    let residence = (10.0 + 6.0 * dt - a) / 0.5;
    let (next, report) = initial
        .trial_with_sources(dt, &[300.0; 2], &sources, budget(), || Ok(()))
        .unwrap();
    close(next.inventories()[0].populations()[0], a, 1e-9);
    close(
        next.inventories()[1].populations()[0],
        0.3 * residence,
        1e-9,
    );
    assert_eq!(report.produced_nuclei_by_cell, [12.0, 0.0]);
    assert_eq!(report.neutrons.fissions, 0.0);
    let transported = report.transport.unwrap();
    close(transported.integrated_populations_s[0][0], residence, 1e-9);
    close(transported.transfer_counts[0][0], 0.3 * residence, 1e-9);
    close(transported.outflow_counts[0][0], 0.2 * residence, 1e-9);
    for defect in transported.cohort_balance_residuals {
        close(defect, 0.0, 1e-9);
    }
    // Zero-removal receiver's residence is positive despite beginning empty.
    assert!(transported.integrated_populations_s[1][0] > 0.0);
    assert!(
        transported.integrated_populations_s[1][0] < next.inventories()[1].populations()[0] * dt
    );
}

fn nonlinear_reference(duration: f64) -> [f64; 3] {
    fn rhs(y: [f64; 3]) -> [f64; 3] {
        let transfer = 0.6 * y[0];
        let fission = 0.01 * y[1] * y[2];
        [-transfer, transfer - fission, fission]
    }
    fn plus(y: [f64; 3], k: [f64; 3], dt: f64) -> [f64; 3] {
        std::array::from_fn(|i| y[i] + dt * k[i])
    }
    let dt = duration / 10_000.0;
    let mut y = [100.0, 0.0, 20.0];
    for _ in 0..10_000 {
        let a = rhs(y);
        let b = rhs(plus(y, a, dt * 0.5));
        let c = rhs(plus(y, b, dt * 0.5));
        let d = rhs(plus(y, c, dt));
        y = std::array::from_fn(|i| y[i] + dt * (a[i] + 2.0 * b[i] + 2.0 * c[i] + d[i]) / 6.0);
    }
    y
}

#[test]
fn g1_transport_changes_neutron_feedback_and_converges_to_independent_nonlinear_ode() {
    let initial = fixture(2.0, [0.0, 20.0], [100.0, 0.0], movement(0.6));
    let reference = nonlinear_reference(1.0);
    let mut errors = Vec::new();
    for steps in [4, 8, 16, 32] {
        let mut state = initial.clone();
        for _ in 0..steps {
            state
                .advance(1.0 / f64::from(steps), &[300.0; 2], budget(), || Ok(()))
                .unwrap();
        }
        let a = state.inventories()[0].populations()[0];
        let b = state.inventories()[1].populations()[0];
        let n = state.neutrons().neutron_counts()[1];
        close(a + b + n, 120.0, 2e-8);
        close(
            state.inventories()[1].cumulative_reaction_counts()[0],
            n - 20.0,
            2e-8,
        );
        assert!(n > 20.0);
        errors.push((a - reference[0]).abs() + (b - reference[1]).abs() + (n - reference[2]).abs());
    }
    for pair in errors.windows(2) {
        assert!(pair[1] < pair[0] * 0.4, "{errors:?}");
    }
    assert!(errors[3] < 0.002, "{errors:?}");
}

#[test]
fn g4_transport_budget_nonconvergence_and_late_cancellation_are_atomic() {
    let initial = fixture(2.0, [0.0, 20.0], [100.0, 0.0], movement(0.6));
    let mut state = initial.clone();
    let mut too_small = budget();
    too_small.burnup.limits.propagation.maximum_nodes = 3;
    assert!(
        state
            .advance(0.3, &[300.0; 2], too_small, || Ok(()))
            .is_err()
    );
    let mut one_iteration = budget();
    one_iteration.max_iterations = 1;
    assert!(matches!(
        state.advance(0.3, &[300.0; 2], one_iteration, || Ok(())),
        Err(DiffusionBurnupError::Nonconvergence { .. })
    ));
    let mut calls = 0;
    state
        .trial(0.3, &[300.0; 2], budget(), || {
            calls += 1;
            Ok(())
        })
        .unwrap();
    let mut observed = 0;
    assert!(
        state
            .advance(0.3, &[300.0; 2], budget(), || {
                observed += 1;
                if observed == calls {
                    Err(FissionError::Cancelled)
                } else {
                    Ok(())
                }
            })
            .is_err()
    );
    assert_eq!(state.inventories(), initial.inventories());
    assert_eq!(state.transport(), initial.transport());
    assert_eq!(
        state.neutrons().neutron_counts(),
        initial.neutrons().neutron_counts()
    );
    assert!(
        DiffusionBurnupState::new_with_transport(
            state.neutrons().clone(),
            state.transport().unwrap().clone(),
            state.bindings().to_vec(),
            BurnupTransportLimits {
                max_transfers: 0,
                ..limits()
            },
        )
        .is_err()
    );
    let mut geometry = state.neutrons().model().data().clone();
    geometry.cell_volumes_m3[0] = 2.0;
    let geometry = Arc::new(DiffusionModel::new(geometry, budget().diffusion.limits).unwrap());
    assert!(
        DiffusionBurnupState::new_with_transport(
            DiffusionState::new(geometry, vec![0.0, 20.0], vec![]).unwrap(),
            state.transport().unwrap().clone(),
            state.bindings().to_vec(),
            limits(),
        )
        .is_err()
    );
}
