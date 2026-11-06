// Child of thermal::tests: reuse its original material and thermal receiver.
// All isotope laws, cross sections, energies and worths below are manufactured
// equation data, not evaluated data for the named nuclides.
use super::*;
use crate::burnup::production::{BurnupProductionChange, BurnupProductionSchedule};
use crate::burnup::{
    BurnupBudget, BurnupData, BurnupError, BurnupInventory, BurnupLimits, BurnupModel,
    BurnupNuclide, BurnupProduct, BurnupRate, BurnupReaction,
};
use crate::kinetics::{
    CompositionFeedback, KineticsBudget, KineticsError, PointKineticsModel, PointKineticsState,
};
use fs_la::decay::DecayLimits;
use fs_thermochem::decay::DecayLawV1;

fn production_budget() -> BurnupHeatingBudget {
    BurnupHeatingBudget {
        source: BurnupBudget {
            limits: BurnupLimits {
                max_neutron_groups: 4,
                max_reactions: 8,
                max_products: 8,
                propagation: DecayLimits::default(),
            },
            nuclei_tolerance: 1e-8,
            energy_tolerance_j: 1e-8,
        },
        max_composition_feedback_terms: 4,
    }
}

fn energy(deposited_j: f64, radiation_j: f64, escaped_j: f64) -> EnergyChannels {
    EnergyChannels {
        deposited_j,
        radiation_j,
        escaped_j,
    }
}
fn rates(parent: f64, tracer: f64) -> Vec<f64> {
    vec![parent, 0.0, 0.0, tracer]
}
fn receiver(boundary_w: f64) -> DepositionSimulation<'static> {
    thermal(
        vec![cell(1, 1, 2.0, 1200.0, boundary_w)],
        &[(MaterialId(1), &HOT)],
        vec![],
    )
}

fn fed_source(
    capture_m2: f64,
    maximum_k: f64,
    initial_rates: (f64, f64),
    production_changes: &[(f64, f64, f64)],
    flux_changes: &[(f64, f64)],
) -> BurnupProductionCellSource {
    let groups = Arc::new(
        NeutronGroups::new("manufactured-fed-burnup-grid-v1".into(), vec![0.0, 1.0]).unwrap(),
    );
    let model = Arc::new(
        BurnupModel::new(
            BurnupData {
                data_id: "manufactured-fed-burnup-v1".into(),
                groups: Arc::clone(&groups),
                temperature_min_k: 200.0,
                temperature_max_k: maximum_k,
                nuclides: [(27, 60), (27, 61), (28, 60), (26, 56)]
                    .into_iter()
                    .enumerate()
                    .map(|(index, (z, a))| BurnupNuclide {
                        id: NuclideId::new(z, a, 0).unwrap(),
                        law: if index == 0 {
                            DecayLawV1::from_half_life(fs_qty::Time::new(
                                2.0 * core::f64::consts::LN_2,
                            ))
                            .unwrap()
                        } else {
                            DecayLawV1::stable()
                        },
                    })
                    .collect(),
                reactions: vec![
                    BurnupReaction {
                        data_id: "manufactured-fed-capture-v1".into(),
                        parent: 0,
                        rate: BurnupRate::Capture {
                            cross_sections_m2: vec![capture_m2],
                        },
                        products: vec![BurnupProduct {
                            nuclide: 1,
                            yield_per_event: 1.0,
                        }],
                        energy_per_event: energy(2.0, 1.0, 1.0),
                    },
                    BurnupReaction {
                        data_id: "manufactured-fed-decay-v1".into(),
                        parent: 0,
                        rate: BurnupRate::Decay {
                            branch_fraction: 1.0,
                        },
                        products: vec![BurnupProduct {
                            nuclide: 2,
                            yield_per_event: 1.0,
                        }],
                        energy_per_event: energy(4.0, 1.0, 2.0),
                    },
                ],
            },
            production_budget().source.limits,
        )
        .unwrap(),
    );
    let production = Arc::new(
        BurnupProductionSchedule::new(
            Arc::clone(&model),
            rates(initial_rates.0, initial_rates.1),
            production_changes
                .iter()
                .map(|&(time_s, parent, tracer)| BurnupProductionChange {
                    time_s,
                    rates_nuclei_s: rates(parent, tracer),
                })
                .collect(),
            8,
        )
        .unwrap(),
    );
    BurnupProductionCellSource::new(
        BurnupCellSource {
            cell_id: CellId(1),
            material_id: MaterialId(1),
            inventory: BurnupInventory::new(model, vec![0.0, 0.0, 0.0, 10.0]).unwrap(),
            fluxes: NeutronFluxSchedule::new(
                GroupFlux::new(Arc::clone(&groups), vec![1.0]).unwrap(),
                flux_changes
                    .iter()
                    .map(|&(time_s, value)| NeutronFluxChange {
                        time_s,
                        flux: GroupFlux::new(Arc::clone(&groups), vec![value]).unwrap(),
                    })
                    .collect(),
            )
            .unwrap(),
            composition_feedback: vec![],
        },
        production,
    )
    .unwrap()
}

fn with_worth(
    source: BurnupProductionCellSource,
    maximum_nuclei: f64,
) -> BurnupProductionCellSource {
    let mut original = source.source().clone();
    original.composition_feedback = vec![
        CompositionFeedback::new(
            "manufactured-fed-tracer-worth-v1".into(),
            original.inventory.model().nuclides()[3].id,
            10.0,
            0.0,
            maximum_nuclei,
            0.01,
        )
        .unwrap(),
    ];
    BurnupProductionCellSource::new(original, Arc::clone(source.production())).unwrap()
}

// Independent scalar solution for P'=s-(lambda+k_capture)P and its event
// integrals. Both stable products are retained in separate original nodes.
fn interval(parent: f64, supply: f64, capture: f64, duration: f64) -> (f64, f64, f64) {
    let rate = 0.5 + capture;
    let loss = -(-rate * duration).exp_m1();
    let population_time = parent * loss / rate + supply * (duration - loss / rate) / rate;
    (
        parent * (-rate * duration).exp() + supply * loss / rate,
        capture * population_time,
        0.5 * population_time,
    )
}

fn check_populations(source: &BurnupProductionCellSource, expected: [f64; 4]) {
    for (&actual, wanted) in source.inventory().populations().iter().zip(expected) {
        close(actual, wanted, 2e-9);
    }
}
fn check_release(actual: EnergyChannels, captures: f64, decays: f64) {
    close(actual.deposited_j, 2.0 * captures + 4.0 * decays, 1e-8);
    close(actual.radiation_j, captures + decays, 1e-8);
    close(actual.escaped_j, captures + 2.0 * decays, 1e-8);
}

#[test]
fn supplied_nuclei_react_and_heat_the_material_during_the_same_interval() {
    let initial = receiver(3.0);
    let supplied = fed_source(0.5, 1000.0, (100.0, 0.0), &[], &[]);
    let (untimed, _) = supplied
        .trial_flux(
            supplied.reference_flux(0.0).unwrap(),
            1.0,
            300.0,
            production_budget(),
            &mut || Ok(()),
        )
        .unwrap();
    let heating =
        BurnupProductionHeating::new(&initial, vec![supplied], production_budget()).unwrap();
    let (heated, evolved, audit) = heating
        .advance_to(&initial, 1.0, controls(10.0), || Ok(()))
        .unwrap();
    let (parent, captures, decays) = interval(0.0, 100.0, 0.5, 1.0);
    let accepted = evolved.source(CellId(1)).unwrap();
    check_populations(accepted, [parent, captures, decays, 10.0]);
    assert_eq!(accepted.inventory(), untimed.inventory());
    close(evolved.audit().produced_nuclei, 100.0, 1e-9);
    close(evolved.audit().captures, captures, 1e-9);
    close(evolved.audit().decays, decays, 1e-9);
    assert_eq!(evolved.audit().fissions, 0.0);
    assert_eq!(
        accepted.inventory().cumulative_produced_nuclei(),
        rates(100.0, 0.0)
    );
    check_release(evolved.audit().released, captures, decays);
    let heat = 2.0 * captures + 4.0 * decays;
    close(audit.external_heat_j, heat + 3.0, 1e-8);
    close(
        heated.temperature_k(CellId(1)).unwrap(),
        300.0 + (heat + 3.0) / 8.0,
        1e-9,
    );
    assert_eq!(audit.substeps, 1);
    assert_eq!(heated.cells()[&CellId(1)].mass_kg, 2.0);
    assert_eq!(heated.cells()[&CellId(1)].boundary.source_w, 3.0);
}

#[test]
fn feed_and_flux_transitions_clip_together_and_shutdown_retains_decay_heat() {
    let initial = receiver(3.0);
    let supplied = fed_source(
        0.5,
        1000.0,
        (100.0, 0.0),
        &[(0.5, 200.0, 0.0), (1.0, 0.0, 0.0)],
        &[(0.75, 0.0)],
    );
    let flux = supplied.reference_flux(0.0).unwrap();
    assert!(matches!(
        supplied.trial_flux(flux, 0.25, 300.0, production_budget(), &mut || Ok(())),
        Err(FissionThermalError::Invalid(
            "scheduled production requires an absolute-time flux trial"
        ))
    ));
    assert!(matches!(
        supplied.trial_flux_at(0.0, flux, 0.75, 300.0, production_budget(), &mut || Ok(())),
        Err(FissionThermalError::Invalid(
            "production trial crosses a scheduled rate change"
        ))
    ));
    let heating =
        BurnupProductionHeating::new(&initial, vec![supplied], production_budget()).unwrap();
    let (heated, evolved, audit) = heating
        .advance_to(&initial, 2.0, controls(10.0), || Ok(()))
        .unwrap();
    let (mut parent, mut captures, mut decays) = (0.0, 0.0, 0.0);
    for (duration, supply, capture_rate) in [
        (0.5, 100.0, 0.5),
        (0.25, 200.0, 0.5),
        (0.25, 200.0, 0.0),
        (1.0, 0.0, 0.0),
    ] {
        let (next, c, d) = interval(parent, supply, capture_rate, duration);
        parent = next;
        captures += c;
        decays += d;
    }
    check_populations(
        evolved.source(CellId(1)).unwrap(),
        [parent, captures, decays, 10.0],
    );
    check_release(evolved.audit().released, captures, decays);
    close(evolved.audit().produced_nuclei, 150.0, 1e-9);
    close(
        audit.external_heat_j,
        2.0 * captures + 4.0 * decays + 6.0,
        1e-8,
    );
    assert_eq!(audit.substeps, 4);
    let (later, continued, late_audit) = evolved
        .advance_to(&heated, 3.0, controls(10.0), || Ok(()))
        .unwrap();
    let late_decays = parent * -(-0.5_f64).exp_m1();
    assert_eq!(continued.audit().produced_nuclei, 0.0);
    assert_eq!(continued.audit().captures, 0.0);
    close(continued.audit().decays, late_decays, 1e-9);
    close(late_audit.external_heat_j, 4.0 * late_decays + 3.0, 1e-8);
    assert!(later.temperature_k(CellId(1)).unwrap() > heated.temperature_k(CellId(1)).unwrap());
}

#[test]
fn neutron_driven_production_uses_absolute_time_and_original_composition_worth() {
    let initial = receiver(3.0);
    let supplied = with_worth(
        fed_source(
            0.25,
            1000.0,
            (100.0, 2.0),
            &[(0.5, 200.0, 4.0), (1.0, 0.0, 0.0)],
            &[(1.0, 0.0)],
        ),
        100.0,
    );
    let neutron_model = Arc::new(
        PointKineticsModel::new("manufactured-fed-neutron-region-v1".into(), 1.0, vec![]).unwrap(),
    );
    let heating = KineticHeating::<BurnupProductionCellSource>::new(
        &initial,
        vec![KineticCellSource {
            irradiation: supplied,
            kinetics: PointKineticsState::new(neutron_model, 2.0, vec![]).unwrap(),
            controls: KineticsControlSchedule::new(
                KineticsControl {
                    reactivity: 0.0,
                    external_source_per_s: 0.0,
                },
                vec![],
            )
            .unwrap(),
            feedback: vec![],
        }],
        KineticHeatingBudget {
            irradiation: production_budget(),
            kinetics: KineticsBudget {
                max_delayed_groups: 0,
                max_squarings: 32,
                max_taylor_terms: 24,
                population_tolerance: 1e-8,
                integral_tolerance_population_s: 1e-8,
            },
            max_feedback_terms: 0,
        },
    )
    .unwrap();
    let (heated, evolved, audit) = heating
        .advance_to(&initial, 1.5, controls(10.0), || Ok(()))
        .unwrap();
    // Stable tracer supply makes its populations 10 -> 11 -> 13 -> 13.
    // The original worth therefore gives frozen rho=0, 0.01 and 0.03.
    // Check the bridge's declared mean-flux partition analytically; no claim
    // of exact product birth timing under continuously changing flux is made.
    let mean_second = 2.0 * 0.005_f64.exp_m1() / 0.005;
    let n_before_last = 2.0 * 0.005_f64.exp();
    let (p1, c1, d1) = interval(0.0, 100.0, 0.5, 0.5);
    let (p2, c2, d2) = interval(p1, 200.0, 0.25 * mean_second, 0.5);
    let (p3, c3, d3) = interval(p2, 0.0, 0.0, 0.5);
    let (captures, decays) = (c1 + c2 + c3, d1 + d2 + d3);
    let accepted = evolved.source(CellId(1)).unwrap();
    check_populations(&accepted.irradiation, [p3, captures, decays, 13.0]);
    close(evolved.audit().produced_nuclei, 153.0, 1e-9);
    close(
        accepted.kinetics.neutron_population(),
        2.0 * 0.02_f64.exp(),
        1e-10,
    );
    let last = evolved.last_step(CellId(1)).unwrap();
    close(last.applied_reactivity, 0.03, 1e-12);
    close(
        last.kinetics.integrated_population_s,
        n_before_last * 0.015_f64.exp_m1() / 0.03,
        1e-10,
    );
    check_release(evolved.audit().released, captures, decays);
    let heat = 2.0 * captures + 4.0 * decays;
    close(audit.external_heat_j, heat + 4.5, 1e-8);
    close(
        heated.temperature_k(CellId(1)).unwrap(),
        300.0 + (heat + 4.5) / 8.0,
        1e-9,
    );
    assert_eq!(audit.substeps, 3);
}

#[test]
fn cloned_production_cache_and_restart_retain_the_original_interval_history() {
    let initial = receiver(0.0);
    let heating = BurnupProductionHeating::new(
        &initial,
        vec![fed_source(
            0.5,
            1000.0,
            (100.0, 0.0),
            &[(0.5, 200.0, 0.0)],
            &[],
        )],
        production_budget(),
    )
    .unwrap();
    let (whole_t, whole_s, _) = heating
        .advance_to(&initial, 1.0, controls(0.25), || Ok(()))
        .unwrap();
    let mut cold_polls = 0;
    let (first_t, first_s, _) = heating
        .advance_to(&initial, 0.25, controls(0.25), || {
            cold_polls += 1;
            Ok(())
        })
        .unwrap();
    let mut warm_polls = 0;
    let (middle_t, middle_s, _) = first_s
        .advance_to(&first_t, 0.5, controls(0.25), || {
            warm_polls += 1;
            Ok(())
        })
        .unwrap();
    assert!(
        warm_polls < cold_polls,
        "an unchanged operator must reuse production preparation"
    );
    let checkpoint_t = middle_t.clone();
    let checkpoint_s = middle_s.clone();
    let (split_t, split_s, _) = checkpoint_s
        .advance_to(&checkpoint_t, 1.0, controls(0.25), || Ok(()))
        .unwrap();
    assert_eq!(whole_t.cells(), split_t.cells());
    assert_eq!(
        whole_s.source(CellId(1)).unwrap().inventory(),
        split_s.source(CellId(1)).unwrap().inventory()
    );
    assert_eq!(
        whole_s.cumulative_released().unwrap(),
        split_s.cumulative_released().unwrap()
    );
    assert_eq!(whole_s.time_s(), split_s.time_s());
}

#[test]
fn model_work_temperature_worth_and_checkpoint_refusals_preserve_the_accepted_source() {
    let initial = receiver(0.0);
    let supplied = fed_source(0.5, 301.0, (100.0, 0.0), &[], &[]);
    let saved = supplied.inventory().clone();
    let other_model = Arc::new((**supplied.inventory().model()).clone());
    let wrong_schedule =
        Arc::new(BurnupProductionSchedule::new(other_model, rates(100.0, 0.0), vec![], 0).unwrap());
    assert!(matches!(
        BurnupProductionCellSource::new(supplied.source().clone(), wrong_schedule),
        Err(BurnupError::Invalid(_))
    ));
    let mut limited = production_budget();
    limited.source.limits.max_reactions = 1;
    assert!(BurnupProductionHeating::new(&initial, vec![supplied.clone()], limited).is_err());
    assert!(matches!(
        supplied.trial_source(0.0, 0.1, 300.0, production_budget(), &mut || Err(
            FissionError::Budget("original production checkpoint")
        )),
        Err(FissionThermalError::Fission(FissionError::Budget(
            "original production checkpoint"
        )))
    ));
    let heating =
        BurnupProductionHeating::new(&initial, vec![supplied], production_budget()).unwrap();
    assert!(matches!(
        heating.advance_to(&initial, 1.0, controls(10.0), || Ok(())),
        Err(FissionThermalError::Burnup(BurnupError::Invalid(_)))
    ));
    assert_eq!(heating.source(CellId(1)).unwrap().inventory(), &saved);
    assert_eq!(heating.time_s(), 0.0);
    assert!(heating.last_step(CellId(1)).is_none());

    let bounded = with_worth(fed_source(0.5, 1000.0, (100.0, 2.0), &[], &[]), 10.25);
    let original = bounded.inventory().clone();
    let heating =
        BurnupProductionHeating::new(&initial, vec![bounded], production_budget()).unwrap();
    assert!(matches!(
        heating.advance_to(&initial, 1.0, controls(10.0), || Ok(())),
        Err(FissionThermalError::Kinetics(KineticsError::InvalidInput(
            _
        )))
    ));
    assert_eq!(heating.source(CellId(1)).unwrap().inventory(), &original);
    assert_eq!(initial.time_s(), 0.0);
    assert_eq!(initial.temperature_k(CellId(1)).unwrap(), 300.0);
}
