//! Coupled temperature-sensitive reaction oracles on the original caloric law.
//! All coefficients are manufactured and do not calibrate a named isotope.
use super::*;
use crate::burnup::{
    BurnupBudget, BurnupData, BurnupInventory, BurnupLimits, BurnupModel, BurnupNuclide,
    BurnupProduct, BurnupRate, BurnupReaction,
};
use crate::cross_sections::{
    CrossSectionData, CrossSectionLimits, CrossSectionTable, NeutronReaction,
};
use crate::irradiation::TabulatedIrradiationData;
use crate::kinetics::{
    CompositionFeedback, KineticsBudget, PointKineticsModel, PointKineticsState,
    TemperatureFeedback, TemperatureFeedbackLaw,
};
use fs_thermochem::decay::DecayLawV1;

fn temperature_table(reaction: NeutronReaction) -> Arc<CrossSectionTable> {
    Arc::new(
        CrossSectionTable::new(
            CrossSectionData {
                data_id: format!("manufactured-coupled-T-{reaction:?}"),
                nuclide: NuclideId::new(92, 235, 0).unwrap(),
                reaction,
                groups: Arc::new(
                    NeutronGroups::new("manufactured-coupled-T-grid".into(), vec![0.0, 1.0])
                        .unwrap(),
                ),
                temperatures_k: vec![300.0, 400.0],
                cross_sections_m2: match reaction {
                    NeutronReaction::Fission => vec![vec![0.2], vec![0.1]],
                    NeutronReaction::Capture
                    | NeutronReaction::ProtonProduction
                    | NeutronReaction::DeuteronProduction
                    | NeutronReaction::TritonProduction
                    | NeutronReaction::Helium3Production
                    | NeutronReaction::AlphaProduction => vec![vec![0.0], vec![0.0]],
                },
            },
            CrossSectionLimits {
                max_temperatures: 2,
                max_neutron_groups: 1,
                max_values: 2,
            },
        )
        .unwrap(),
    )
}
fn burnup_budget() -> BurnupHeatingBudget {
    BurnupHeatingBudget {
        source: BurnupBudget {
            limits: BurnupLimits::default(),
            nuclei_tolerance: 1e-7,
            energy_tolerance_j: 1e-7,
        },
        max_composition_feedback_terms: 2,
    }
}
fn burnup_source(maximum_k: f64) -> BurnupCellSource {
    let xs = temperature_table(NeutronReaction::Fission);
    let groups = xs.groups().clone();
    let nuclides = [(92, 235), (54, 135), (38, 98)]
        .into_iter()
        .map(|(z, a)| BurnupNuclide {
            id: NuclideId::new(z, a, 0).unwrap(),
            law: DecayLawV1::stable(),
        })
        .collect();
    let model = Arc::new(
        BurnupModel::new(
            BurnupData {
                data_id: "manufactured-coupled-T-products".into(),
                groups: groups.clone(),
                temperature_min_k: 300.0,
                temperature_max_k: maximum_k,
                nuclides,
                reactions: vec![BurnupReaction {
                    data_id: "manufactured-T-fission".into(),
                    parent: 0,
                    rate: BurnupRate::TabulatedFission { cross_sections: xs },
                    products: vec![
                        BurnupProduct {
                            nuclide: 1,
                            yield_per_event: 1.0,
                        },
                        BurnupProduct {
                            nuclide: 2,
                            yield_per_event: 1.0,
                        },
                    ],
                    energy_per_event: EnergyChannels {
                        deposited_j: 8.0,
                        ..EnergyChannels::default()
                    },
                }],
            },
            burnup_budget().source.limits,
        )
        .unwrap(),
    );
    BurnupCellSource {
        cell_id: CellId(1),
        material_id: MaterialId(1),
        inventory: BurnupInventory::new(model, vec![100.0, 0.0, 0.0]).unwrap(),
        fluxes: NeutronFluxSchedule::new(GroupFlux::new(groups, vec![1.0]).unwrap(), vec![])
            .unwrap(),
        composition_feedback: vec![],
    }
}
fn parent_source() -> IrradiationCellSource {
    let fission = temperature_table(NeutronReaction::Fission);
    let groups = fission.groups().clone();
    let model = Arc::new(
        IrradiationModel::new_tabulated(TabulatedIrradiationData {
            data_id: "manufactured-coupled-T-parent".into(),
            fission,
            capture: temperature_table(NeutronReaction::Capture),
            temperature_min_k: 300.0,
            temperature_max_k: 400.0,
            fission_model: Arc::new(
                FissionModel::new(
                    "manufactured-T-energy".into(),
                    EnergyChannels {
                        deposited_j: 8.0,
                        ..EnergyChannels::default()
                    },
                    vec![],
                )
                .unwrap(),
            ),
            capture_energy_per_event: EnergyChannels::default(),
        })
        .unwrap(),
    );
    IrradiationCellSource {
        cell_id: CellId(1),
        material_id: MaterialId(1),
        inventories: vec![IrradiationState::new(model, 100.0, vec![]).unwrap()],
        fluxes: NeutronFluxSchedule::new(GroupFlux::new(groups, vec![1.0]).unwrap(), vec![])
            .unwrap(),
    }
}

#[test]
fn original_material_temperature_changes_both_reaction_owners_and_refines_to_exact_cooling_feedback()
 {
    let initial = thermal(
        vec![cell(1, 1, 2.0, 1200.0, 0.0)],
        &[(MaterialId(1), &HOT)],
        vec![],
    );
    let products =
        BurnupHeating::new(&initial, vec![burnup_source(400.0)], burnup_budget()).unwrap();
    let parents =
        IrradiationHeating::new(&initial, vec![parent_source()], irradiation_budget()).unwrap();
    // Native heat capacity is 8 J/K and each event deposits8 J: T=300+100-N.
    // sigma(T)=0.2-0.001*(T-300), hence N'=-(0.1+0.001*N)*N.
    let exact = 100.0 / (2.0 * (0.1_f64 * 2.0).exp() - 1.0);
    let mut errors = Vec::new();
    for dt in [0.2, 0.1, 0.05] {
        let (heated, products, _) = products
            .advance_to(&initial, 2.0, controls(dt), || Ok(()))
            .unwrap();
        let (parent_heat, parents, _) = parents
            .advance_to(&initial, 2.0, controls(dt), || Ok(()))
            .unwrap();
        let remaining = products.source(CellId(1)).unwrap().inventory.populations()[0];
        close(
            remaining,
            parents.source(CellId(1)).unwrap().inventories[0].remaining_parent_nuclei(),
            2e-10,
        );
        close(
            heated.temperature_k(CellId(1)).unwrap(),
            400.0 - remaining,
            2e-10,
        );
        close(heated.energy_j(), parent_heat.energy_j(), 3e-9);
        close(products.audit().fissions, 100.0 - remaining, 2e-10);
        assert!(remaining > 100.0 * (-0.2_f64 * 2.0).exp());
        errors.push((remaining - exact).abs());
    }
    assert!(
        errors
            .windows(2)
            .all(|p| (1.8..2.2).contains(&(p[0] / p[1]))),
        "T-dependent refinement {errors:?}"
    );
    assert!(errors[2] < 0.1, "T-dependent refinement {errors:?}");
}

fn kinetic_budget() -> KineticHeatingBudget<BurnupHeatingBudget> {
    KineticHeatingBudget {
        irradiation: burnup_budget(),
        kinetics: KineticsBudget {
            max_delayed_groups: 4,
            max_squarings: 32,
            max_taylor_terms: 24,
            population_tolerance: 1e-8,
            integral_tolerance_population_s: 1e-8,
        },
        max_feedback_terms: 2,
    }
}
#[test]
fn tabulated_reactions_material_feedback_and_product_worth_share_one_neutron_history() {
    let initial = thermal(
        vec![cell(1, 1, 2.0, 1200.0, 0.0)],
        &[(MaterialId(1), &HOT)],
        vec![],
    );
    let mut irradiation = burnup_source(400.0);
    irradiation.composition_feedback.push(
        CompositionFeedback::new(
            "manufactured-T-poison-worth".into(),
            NuclideId::new(54, 135, 0).unwrap(),
            0.0,
            0.0,
            100.0,
            -0.01,
        )
        .unwrap(),
    );
    let source = KineticCellSource {
        irradiation,
        kinetics: PointKineticsState::new(
            Arc::new(
                PointKineticsModel::new("manufactured-T-kinetics".into(), 1.0, vec![]).unwrap(),
            ),
            1.0,
            vec![],
        )
        .unwrap(),
        controls: KineticsControlSchedule::new(
            KineticsControl {
                reactivity: 0.05,
                external_source_per_s: 0.0,
            },
            vec![],
        )
        .unwrap(),
        feedback: vec![
            TemperatureFeedback::new(
                "manufactured-T-reactivity".into(),
                300.0,
                300.0,
                400.0,
                TemperatureFeedbackLaw::Linear {
                    coefficient_per_kelvin: -0.01,
                },
            )
            .unwrap(),
        ],
    };
    let heating = KineticHeating::new(&initial, vec![source], kinetic_budget()).unwrap();
    let mut errors = Vec::new();
    for dt in [0.1, 0.05, 0.025] {
        let (heated, state, _) = heating
            .advance_to(&initial, 1.0, controls(dt), || Ok(()))
            .unwrap();
        let source = state.source(CellId(1)).unwrap();
        let fuel = source.irradiation.inventory.populations()[0];
        let products = source.irradiation.inventory.populations()[1];
        close(fuel + products, 100.0, 3e-10);
        close(
            heated.temperature_k(CellId(1)).unwrap(),
            300.0 + products,
            3e-10,
        );
        // Eliminate t from F'=-(0.1+0.001F)nF and n'=(-1.95+0.02F)n.
        // Both explicit feedback terms contribute; no duplicate numerical solve.
        let exact_n = 1.0 + 19.5 * (fuel / 100.0).ln() - 39.5 * ((0.1 + 0.001 * fuel) / 0.2).ln();
        errors.push((source.kinetics.neutron_population() - exact_n).abs());
        assert!(source.kinetics.neutron_population() < 1.0);
        assert!(state.last_step(CellId(1)).unwrap().applied_reactivity < 0.0);
    }
    assert!(
        errors
            .windows(2)
            .all(|p| (1.8..2.2).contains(&(p[0] / p[1]))),
        "combined feedback refinement {errors:?}"
    );
    assert!(errors[2] < 0.01, "combined feedback refinement {errors:?}");
}

#[test]
fn a_thermal_endpoint_outside_the_tabulated_model_domain_publishes_no_inventory_or_heat() {
    let initial = thermal(
        vec![cell(1, 1, 2.0, 1200.0, 0.0)],
        &[(MaterialId(1), &HOT)],
        vec![],
    );
    let heating =
        BurnupHeating::new(&initial, vec![burnup_source(300.05)], burnup_budget()).unwrap();
    assert!(matches!(
        heating.advance_to(&initial, 1.0, controls(1.0), || Ok(())),
        Err(FissionThermalError::Burnup(
            crate::burnup::BurnupError::Invalid(_)
        ))
    ));
    assert_eq!(initial.energy_j(), 2400.0);
    assert_eq!(
        heating.source(CellId(1)).unwrap().inventory.populations(),
        &[100.0, 0.0, 0.0]
    );
    assert_eq!(heating.audit().released, EnergyChannels::default());
}
