//! Moving finite fuel and products on the original spatial/thermal clock.

use super::*;
use crate::burnup::transport::{
    BurnupTransportLimits, BurnupTransportState, NuclideTransfer, NuclideTransport,
};

fn moving_fuel(prompt_yield: f64) -> Arc<DiffusionBurnupState> {
    let base = spatial(900.0);
    let fuel = Arc::new(
        BurnupModel::new(
            BurnupData {
                data_id: "synthetic-moving-thermal-fuel".into(),
                groups: Arc::clone(base.neutrons().model().groups()),
                temperature_min_k: 300.0,
                temperature_max_k: 900.0,
                nuclides: base.inventories()[0].model().nuclides().to_vec(),
                reactions: vec![BurnupReaction {
                    data_id: "synthetic-moving-fission".into(),
                    parent: 0,
                    rate: BurnupRate::Fission {
                        cross_sections_m2: vec![0.0001],
                    },
                    products: vec![BurnupProduct {
                        nuclide: 1,
                        yield_per_event: 2.0,
                    }],
                    energy_per_event: EnergyChannels {
                        deposited_j: 8.0,
                        radiation_j: 2.0,
                        escaped_j: 1.0,
                    },
                }],
            },
            spatial_budget().burnup.limits,
        )
        .unwrap(),
    );
    let limits = BurnupTransportLimits {
        burnup: spatial_budget().burnup.limits,
        ..BurnupTransportLimits::default()
    };
    let transport = BurnupTransportState::new(
        vec![
            BurnupInventory::new(Arc::clone(&fuel), vec![0.0, 0.0]).unwrap(),
            BurnupInventory::new(fuel, vec![100.0, 0.0]).unwrap(),
        ],
        vec![0.01; 2],
        NuclideTransport {
            transfers: vec![
                NuclideTransfer {
                    from_compartment: 1,
                    to_compartment: 0,
                    rates_per_s: vec![0.6, 0.0],
                },
                NuclideTransfer {
                    from_compartment: 0,
                    to_compartment: 1,
                    rates_per_s: vec![0.0, 0.3],
                },
            ],
            outflows: vec![],
        },
        limits,
        || Ok(()),
    )
    .unwrap();
    let mut neutron_data = base.neutrons().model().data().clone();
    neutron_data.links.clear();
    for reaction in &mut neutron_data.reactions {
        reaction.kind = DiffusionReactionKind::Fission {
            prompt_yields: vec![prompt_yield],
            delayed_yields: vec![],
        };
    }
    let model =
        Arc::new(DiffusionModel::new(neutron_data, spatial_budget().diffusion.limits).unwrap());
    Arc::new(
        DiffusionBurnupState::new_with_transport(
            DiffusionState::new(model, vec![20.0, 0.0], vec![]).unwrap(),
            transport,
            base.bindings().to_vec(),
            limits,
        )
        .unwrap(),
    )
}

fn sources(state: Arc<DiffusionBurnupState>) -> Vec<DiffusionCellSource> {
    DiffusionCellSource::bind(
        state,
        vec![(CellId(20), MaterialId(1)), (CellId(10), MaterialId(2))],
    )
    .unwrap()
}

#[test]
fn g1_native_heat_uses_arriving_fuel_and_keeps_energy_at_original_event_location() {
    let thermal = two_cells(0.0);
    let original = moving_fuel(1.0);
    let source =
        DiffusionHeating::new(&thermal, sources(Arc::clone(&original)), spatial_budget()).unwrap();
    let (heated, accepted, audit) = source
        .advance_to(&thermal, 1.0, controls(1.0), || Ok(()))
        .unwrap();
    let a = accepted.source(CellId(20)).unwrap();
    let b = accepted.source(CellId(10)).unwrap();
    assert!(Arc::ptr_eq(a.state(), b.state()));
    let donor = 100.0 * (-0.6_f64).exp();
    let receiver = 150.0 * ((-0.2_f64).exp() - (-0.6_f64).exp());
    let fissions = 100.0 - donor - receiver;
    close(a.state().inventories()[0].populations()[0], receiver, 2e-9);
    close(b.state().inventories()[1].populations()[0], donor, 2e-9);
    close(accepted.audit().fissions, fissions, 2e-9);
    close(
        heated.temperature_k(CellId(20)).unwrap(),
        300.0 + fissions,
        2e-8,
    );
    assert_eq!(heated.temperature_k(CellId(10)).unwrap(), 300.0);
    close(audit.external_heat_j, 8.0 * fissions, 2e-8);
    close(accepted.audit().released.radiation_j, 2.0 * fissions, 3e-9);
    close(accepted.audit().released.escaped_j, fissions, 2e-9);
    assert!(b.state().inventories()[1].populations()[1] > 0.0);
    assert_eq!(
        b.state().inventories()[1].cumulative_released().total_j(),
        0.0
    );
    assert_eq!(heated.cells()[&CellId(20)].mass_kg, 2.0);
    assert_eq!(heated.cells()[&CellId(10)].mass_kg, 3.0);
    assert_eq!(
        original.transport().unwrap().cumulative_transfer_counts(),
        vec![vec![0.0; 2]; 2]
    );
    assert!(
        accepted
            .last_step(CellId(20))
            .unwrap()
            .spatial()
            .transport
            .is_some()
    );
    // Each local model fits three nodes, but the joint graph needs four.
    let mut too_small = spatial_budget();
    too_small.burnup.limits.propagation.maximum_nodes = 3;
    assert!(DiffusionHeating::new(&thermal, sources(original), too_small).is_err());
}

#[test]
fn g3_adaptive_spatial_fission_preserves_shared_transport_history_and_path_counts() {
    let thermal = two_cells(0.0);
    let heating =
        DiffusionHeating::new(&thermal, sources(moving_fuel(2.0)), spatial_budget()).unwrap();
    let (reference_heat, reference, _) = heating
        .advance_to(&thermal, 1.0, controls(1.0 / 512.0), || Ok(()))
        .unwrap();
    let control = NuclearAdaptiveControls {
        step: StepConfig {
            max_dt_s: 1.0,
            energy_tolerance_j: 1e-6,
        },
        initial_dt_s: 1.0,
        minimum_dt_s: 1e-8,
        maximum_trials: 4096,
        maximum_substeps: 4096,
        tolerances: NuclearTimeTolerances {
            specific_energy_j_kg: 1e-5,
            temperature_k: 1e-5,
            nuclei: 1e-6,
            nuclear_energy_j: 1e-5,
            neutron_count: 1e-6,
            relative: 0.0,
            ..NuclearTimeTolerances::default()
        },
    };
    let (heated, accepted, report) = heating
        .advance_adaptive_to(&thermal, 1.0, control, || Ok(()))
        .unwrap();
    let a = accepted.source(CellId(20)).unwrap();
    let b = accepted.source(CellId(10)).unwrap();
    assert!(Arc::ptr_eq(a.state(), b.state()));
    assert!(report.accuracy_rejections > 0);
    assert!(report.maximum_accepted_normalized_error <= 1.0);
    assert_eq!(report.thermal.substeps, 2 * report.accepted_comparisons);
    close(
        heated.temperature_k(CellId(20)).unwrap(),
        reference_heat.temperature_k(CellId(20)).unwrap(),
        2e-4,
    );
    let fissions = a.state().neutrons().neutron_counts()[0] - 20.0;
    close(accepted.audit().fissions, fissions, 2e-8);
    close(accepted.audit().released.deposited_j, 8.0 * fissions, 2e-8);
    let actual = a.state().transport().unwrap().cumulative_transfer_counts();
    let expected = reference
        .source(CellId(20))
        .unwrap()
        .state()
        .transport()
        .unwrap()
        .cumulative_transfer_counts();
    for (actual, expected) in actual.iter().flatten().zip(expected.iter().flatten()) {
        close(*actual, *expected, 2e-4);
    }
    assert!(a.state().inventories()[1].populations()[1] > 0.0);
    assert_eq!(heated.temperature_k(CellId(10)).unwrap(), 300.0);
}
