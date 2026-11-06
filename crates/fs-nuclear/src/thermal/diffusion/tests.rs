use super::*;
use crate::burnup::{
    BurnupBudget, BurnupData, BurnupInventory, BurnupLimits, BurnupModel, BurnupNuclide,
    BurnupProduct, BurnupRate, BurnupReaction,
};
use crate::cross_sections::{
    CrossSectionData, CrossSectionLimits, CrossSectionTable, NeutronReaction,
};
use crate::diffusion::burnup::{
    DiffusionBurnupBudget, DiffusionBurnupError, DiffusionBurnupState, ReactionBinding,
};
use crate::diffusion::{
    DiffusionBudget, DiffusionData, DiffusionLimits, DiffusionLink, DiffusionModel,
    DiffusionReaction, DiffusionReactionKind, DiffusionState,
};
use fs_qty::Frequency;
use fs_thermochem::decay::DecayLawV1;

#[path = "forcing_tests.rs"]
mod forcing_tests;

fn spatial_budget() -> DiffusionBurnupBudget {
    DiffusionBurnupBudget {
        diffusion: DiffusionBudget {
            limits: DiffusionLimits::default(),
            population_tolerance: 1e-8,
        },
        burnup: BurnupBudget {
            limits: BurnupLimits::default(),
            nuclei_tolerance: 1e-8,
            energy_tolerance_j: 1e-8,
        },
        max_iterations: 64,
        reaction_tolerance: 1e-11,
        relaxation: 1.0,
    }
}
fn spatial(maximum_k: f64) -> Arc<DiffusionBurnupState> {
    let groups = Arc::new(
        NeutronGroups::new("manufactured-thermal-neutron-grid".into(), vec![0.0, 1.0]).unwrap(),
    );
    let parent = NuclideId::new(10, 25, 0).unwrap();
    let table = Arc::new(
        CrossSectionTable::new(
            CrossSectionData {
                data_id: "manufactured-heating-dependent-capture".into(),
                nuclide: parent,
                reaction: NeutronReaction::Capture,
                groups: Arc::clone(&groups),
                temperatures_k: vec![300.0, maximum_k],
                cross_sections_m2: vec![vec![0.0001], vec![0.0002]],
            },
            CrossSectionLimits::default(),
        )
        .unwrap(),
    );
    let model = Arc::new(
        BurnupModel::new(
            BurnupData {
                data_id: "manufactured-capture-and-product-heat".into(),
                groups: Arc::clone(&groups),
                temperature_min_k: 300.0,
                temperature_max_k: maximum_k,
                nuclides: vec![parent, NuclideId::new(10, 26, 0).unwrap()]
                    .into_iter()
                    .map(|id| BurnupNuclide {
                        id,
                        law: DecayLawV1::from_decay_constant(Frequency::new(0.0)).unwrap(),
                    })
                    .collect(),
                reactions: vec![BurnupReaction {
                    data_id: "manufactured-capture".into(),
                    parent: 0,
                    rate: BurnupRate::TabulatedCapture {
                        cross_sections: table,
                    },
                    products: vec![BurnupProduct {
                        nuclide: 1,
                        yield_per_event: 1.0,
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
    let inventories = vec![0.0, 100.0]
        .into_iter()
        .map(|fuel| BurnupInventory::new(Arc::clone(&model), vec![fuel, 0.0]).unwrap())
        .collect();
    let neutron_model = Arc::new(
        DiffusionModel::new(
            DiffusionData {
                data_id: "manufactured-two-cell-neutron-field".into(),
                groups,
                speeds_m_s: vec![1.0],
                cell_volumes_m3: vec![0.01, 0.01],
                delayed_groups: vec![],
                reactions: (0..2)
                    .map(|cell| DiffusionReaction {
                        data_id: format!("capture-{cell}"),
                        cell,
                        incident_group: 0,
                        macroscopic_cross_section_m_inv: if cell == 0 { 0.0 } else { 1.0 },
                        kind: DiffusionReactionKind::Capture,
                    })
                    .collect(),
                links: vec![DiffusionLink {
                    cell_a: 0,
                    cell_b: 1,
                    conductance_m2: vec![0.005],
                }],
                leakage: vec![],
            },
            spatial_budget().diffusion.limits,
        )
        .unwrap(),
    );
    Arc::new(
        DiffusionBurnupState::new(
            DiffusionState::new(neutron_model, vec![20.0, 0.0], vec![]).unwrap(),
            inventories,
            vec![
                ReactionBinding {
                    diffusion_reaction: 0,
                    burnup_reaction: 0,
                },
                ReactionBinding {
                    diffusion_reaction: 1,
                    burnup_reaction: 0,
                },
            ],
        )
        .unwrap(),
    )
}
fn bound(
    thermal: &DepositionSimulation<'_>,
    maximum_k: f64,
    budget: DiffusionBurnupBudget,
) -> DiffusionHeating {
    let sources = DiffusionCellSource::bind(
        spatial(maximum_k),
        vec![(CellId(20), MaterialId(1)), (CellId(10), MaterialId(2))],
    )
    .unwrap();
    DiffusionHeating::new(thermal, sources, budget).unwrap()
}
fn two_cells(boundary_source: f64) -> DepositionSimulation<'static> {
    thermal(
        vec![
            cell(20, 1, 2.0, 1200.0, boundary_source),
            cell(10, 2, 3.0, 600.0, 0.0),
        ],
        &[(MaterialId(1), &HOT), (MaterialId(2), &COOL)],
        vec![],
    )
}

#[test]
fn g1_spatial_fluence_fuel_products_and_original_material_heat_share_one_clock() {
    let initial = two_cells(1.0);
    let source = bound(&initial, 900.0, spatial_budget());
    let (heated, accepted, _) = source
        .advance_to(&initial, 1.0, controls(0.1), || Ok(()))
        .unwrap();
    let left = accepted.source(CellId(20)).unwrap();
    let right = accepted.source(CellId(10)).unwrap();
    assert!(Arc::ptr_eq(left.state(), right.state()));
    let state = right.state();
    let captures = 100.0 - state.inventories()[1].populations()[0];
    assert!(captures > 1.0);
    close(state.inventories()[1].populations()[1], captures, 1e-8);
    close(
        20.0 - state.neutrons().neutron_counts().iter().sum::<f64>(),
        captures,
        1e-8,
    );
    close(
        heated.temperature_k(CellId(10)).unwrap(),
        300.0 + 8.0 * captures / 6.0,
        1e-8,
    );
    close(
        heated.temperature_k(CellId(20)).unwrap(),
        300.0 + 1.0 / 8.0,
        1e-10,
    );
    close(accepted.audit().thermal_credit_j, 8.0 * captures, 1e-8);
    close(accepted.audit().captures, captures, 1e-8);
    close(accepted.audit().released.radiation_j, 2.0 * captures, 1e-8);
    assert_eq!(accepted.time_s(), heated.time_s());
    assert_eq!(heated.cells()[&CellId(10)].mass_kg, 3.0);
    assert_eq!(heated.cells()[&CellId(20)].mass_kg, 2.0);
    // The material's changed temperature actually rebuilt its microscopic rate.
    let sigma = state.neutrons().model().reactions()[1].macroscopic_cross_section_m_inv;
    assert!(sigma > state.inventories()[1].populations()[0] * 0.0001 / 0.01);
    let report = accepted.last_step(CellId(10)).unwrap().spatial();
    assert!(report.reaction_defect <= spatial_budget().reaction_tolerance);
    assert_eq!(
        source.source(CellId(10)).unwrap().state().inventories()[1].populations(),
        &[100.0, 0.0]
    );
}

#[test]
fn g4_material_domain_and_nonlinear_refusals_roll_back_all_spatial_and_thermal_history() {
    let initial = two_cells(0.0);
    let source = bound(&initial, 300.01, spatial_budget());
    assert!(
        source
            .advance_to(&initial, 1.0, controls(0.5), || Ok(()))
            .is_err()
    );
    assert_eq!(initial.temperature_k(CellId(10)).unwrap(), 300.0);
    assert_eq!(source.time_s(), 0.0);
    assert_eq!(
        source
            .source(CellId(10))
            .unwrap()
            .state()
            .neutrons()
            .neutron_counts(),
        &[20.0, 0.0]
    );
    assert_eq!(
        source.source(CellId(10)).unwrap().state().inventories()[1].populations(),
        &[100.0, 0.0]
    );
    let mut strict = spatial_budget();
    strict.max_iterations = 1;
    let source = bound(&initial, 900.0, strict);
    assert!(matches!(
        source.advance_to(&initial, 1.0, controls(0.5), || Ok(())),
        Err(FissionThermalError::Diffusion(
            DiffusionBurnupError::Nonconvergence { .. }
        ))
    ));
    assert!(matches!(
        source.advance_to(&initial, 0.1, controls(0.1), || Err(
            DepositionError::MaterialDomain
        )),
        Err(FissionThermalError::Thermal(
            DepositionError::MaterialDomain
        ))
    ));
    assert_eq!(source.time_s(), 0.0);
}

#[test]
fn g0_geometry_and_complete_shared_cell_correspondence_are_required() {
    let initial = two_cells(0.0);
    let mut wrong_cell = cell(20, 1, 2.0, 1200.0, 0.0);
    wrong_cell.geometry.volume_m3 = 0.02;
    let wrong = thermal(
        vec![wrong_cell, cell(10, 2, 3.0, 600.0, 0.0)],
        &[(MaterialId(1), &HOT), (MaterialId(2), &COOL)],
        vec![],
    );
    let state = spatial(900.0);
    let sources = DiffusionCellSource::bind(
        Arc::clone(&state),
        vec![(CellId(20), MaterialId(1)), (CellId(10), MaterialId(2))],
    )
    .unwrap();
    assert!(matches!(
        DiffusionHeating::new(&wrong, sources.clone(), spatial_budget()),
        Err(FissionThermalError::BindingChanged(CellId(20)))
    ));
    let incomplete =
        DiffusionHeating::new(&initial, vec![sources[0].clone()], spatial_budget()).unwrap();
    assert!(matches!(
        incomplete.advance_to(&initial, 0.1, controls(0.1), || Ok(())),
        Err(FissionThermalError::Invalid(_))
    ));
    assert!(DiffusionCellSource::bind(state, vec![(CellId(20), MaterialId(1))]).is_err());
}
