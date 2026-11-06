// Reuse the original thermal test fixtures. Nuclear constants below are
// manufactured equation data, not an evaluated isotope validation deck.
use super::*;
use crate::burnup::production::{BurnupProductionChange, BurnupProductionSchedule};
use crate::burnup::{
    BurnupBudget, BurnupData, BurnupInventory, BurnupLimits, BurnupModel, BurnupNuclide,
    BurnupProduct, BurnupRate, BurnupReaction,
};
use crate::kinetics::{
    KineticsBudget, PointKineticsModel, PointKineticsState, TemperatureFeedback,
    TemperatureFeedbackLaw,
};
use fs_thermochem::decay::DecayLawV1;

fn adaptive(maximum: f64, temperature: f64) -> NuclearAdaptiveControls {
    NuclearAdaptiveControls {
        step: StepConfig {
            max_dt_s: maximum,
            energy_tolerance_j: 1e-6,
        },
        initial_dt_s: maximum,
        minimum_dt_s: 1e-8,
        maximum_trials: 20_000,
        maximum_substeps: 20_000,
        tolerances: NuclearTimeTolerances {
            specific_energy_j_kg: 4.0 * temperature,
            temperature_k: temperature,
            nuclei: 1e-3,
            nuclear_energy_j: 1e-3,
            neutron_population: 1e-4,
            neutron_count: 1e-4,
            relative: 0.0,
        },
    }
}

fn cooling() -> (DepositionSimulation<'static>, FissionHeating) {
    let mut initial_cell = cell(1, 1, 2.0, 1600.0, 3.0);
    initial_cell.boundary.htc_w_m2_k = 2.0;
    let initial = thermal(vec![initial_cell], &[(MaterialId(1), &HOT)], vec![]);
    let mut law = source(1, 1, 0.0, &[], true);
    law.state = FissionState::new(model(true), 1000.0, vec![40.0]).unwrap();
    let heating = FissionHeating::new(&initial, vec![law], budget()).unwrap();
    (initial, heating)
}

#[test]
fn g1_adaptive_decay_cooling_refines_native_cfl_and_keeps_full_accounts() {
    let (initial, heating) = cooling();
    // C=8 J/K, G=2 W/K, delayed deposited energy=20 J with lambda=.5/s.
    // The independently integrated ODE includes the authored 3 W boundary.
    let time = 4.0;
    let exact = 300.0
        + 100.0 * (-0.25_f64 * time).exp()
        + 1.5 * -(-0.25_f64 * time).exp_m1()
        + 1.25 * ((-0.5_f64 * time).exp() - (-0.25_f64 * time).exp()) / (0.25 - 0.5);
    let mut errors = Vec::new();
    for tolerance in [0.05, 0.0005] {
        let (heated, evolved, report) = heating
            .advance_adaptive_to(&initial, time, adaptive(4.0, tolerance), || Ok(()))
            .unwrap();
        errors.push((heated.temperature_k(CellId(1)).unwrap() - exact).abs());
        assert!(
            report.native_interval_cuts > 0,
            "native CFL must cut the initial 4 s span"
        );
        assert!(report.accuracy_rejections > 0);
        assert!(report.maximum_accepted_normalized_error <= 1.0);
        assert_eq!(report.thermal.substeps, 2 * report.accepted_comparisons);
        assert_eq!(evolved.audit().substeps, report.thermal.substeps);
        close(
            evolved.audit().released.deposited_j,
            20.0 * -(-0.5 * time).exp_m1(),
            2e-12,
        );
        close(
            evolved.source(CellId(1)).unwrap().state.delayed_energy_j()[0],
            40.0 * (-0.5 * time).exp(),
            2e-12,
        );
        close(
            report.thermal.stored_energy_change_j,
            heated.energy_j() - initial.energy_j(),
            2e-10,
        );
        assert!(report.thermal.residual_j.abs() < 1e-6);
        assert_eq!(heated.cells()[&CellId(1)].boundary.source_w, 3.0);
        assert_eq!(evolved.audit().fissions, 0.0);
    }
    assert!(errors[1] < errors[0] / 5.0, "actual ODE errors {errors:?}");
    assert!(errors[1] < 0.1, "actual ODE errors {errors:?}");
}

fn feed_budget() -> BurnupHeatingBudget {
    BurnupHeatingBudget {
        source: BurnupBudget {
            limits: BurnupLimits::default(),
            nuclei_tolerance: 1e-8,
            energy_tolerance_j: 1e-8,
        },
        max_composition_feedback_terms: 0,
    }
}

fn fed_kinetics() -> KineticCellSource<BurnupProductionCellSource> {
    let groups =
        Arc::new(NeutronGroups::new("manufactured-adaptive-grid".into(), vec![0.0, 1.0]).unwrap());
    let model = Arc::new(
        BurnupModel::new(
            BurnupData {
                data_id: "manufactured-adaptive-fed-capture-decay".into(),
                groups: Arc::clone(&groups),
                temperature_min_k: 200.0,
                temperature_max_k: 1e6,
                nuclides: vec![
                    BurnupNuclide {
                        id: NuclideId::new(1, 3, 0).unwrap(),
                        law: DecayLawV1::from_half_life(fs_qty::Time::new(1.0)).unwrap(),
                    },
                    BurnupNuclide {
                        id: NuclideId::new(2, 3, 0).unwrap(),
                        law: DecayLawV1::stable(),
                    },
                ],
                reactions: vec![
                    BurnupReaction {
                        data_id: "manufactured-capture".into(),
                        parent: 0,
                        rate: BurnupRate::Capture {
                            cross_sections_m2: vec![1.0],
                        },
                        products: vec![BurnupProduct {
                            nuclide: 1,
                            yield_per_event: 1.0,
                        }],
                        energy_per_event: EnergyChannels {
                            deposited_j: 2.0,
                            radiation_j: 1.0,
                            escaped_j: 0.0,
                        },
                    },
                    BurnupReaction {
                        data_id: "manufactured-decay".into(),
                        parent: 0,
                        rate: BurnupRate::Decay {
                            branch_fraction: 1.0,
                        },
                        products: vec![BurnupProduct {
                            nuclide: 1,
                            yield_per_event: 1.0,
                        }],
                        energy_per_event: EnergyChannels {
                            deposited_j: 2.0,
                            radiation_j: 1.0,
                            escaped_j: 0.0,
                        },
                    },
                ],
            },
            feed_budget().source.limits,
        )
        .unwrap(),
    );
    let production = Arc::new(
        BurnupProductionSchedule::new(
            Arc::clone(&model),
            vec![10.0, 0.0],
            vec![BurnupProductionChange {
                time_s: 1.0,
                rates_nuclei_s: vec![0.0, 0.0],
            }],
            2,
        )
        .unwrap(),
    );
    let source = BurnupProductionCellSource::new(
        BurnupCellSource {
            cell_id: CellId(1),
            material_id: MaterialId(1),
            inventory: BurnupInventory::new(model, vec![0.0, 0.0]).unwrap(),
            fluxes: NeutronFluxSchedule::new(
                GroupFlux::new(Arc::clone(&groups), vec![1.0]).unwrap(),
                vec![NeutronFluxChange {
                    time_s: 1.5,
                    flux: GroupFlux::new(groups, vec![0.0]).unwrap(),
                }],
            )
            .unwrap(),
            composition_feedback: vec![],
        },
        production,
    )
    .unwrap();
    let neutrons = Arc::new(
        PointKineticsModel::new("manufactured-linear-neutron-source".into(), 1.0, vec![]).unwrap(),
    );
    KineticCellSource {
        irradiation: source,
        kinetics: PointKineticsState::new(neutrons, 1.0, vec![]).unwrap(),
        controls: KineticsControlSchedule::new(
            KineticsControl {
                reactivity: 0.0,
                external_source_per_s: 1.0,
            },
            vec![KineticsControlChange {
                time_s: 1.0,
                control: KineticsControl {
                    reactivity: 0.0,
                    external_source_per_s: 0.0,
                },
            }],
        )
        .unwrap(),
        feedback: vec![],
    }
}

fn fed_budget() -> KineticHeatingBudget<BurnupHeatingBudget> {
    KineticHeatingBudget {
        irradiation: feed_budget(),
        kinetics: neutron_budget(),
        max_feedback_terms: 2,
    }
}

fn neutron_budget() -> KineticsBudget {
    KineticsBudget {
        max_delayed_groups: 8,
        max_squarings: 32,
        max_taylor_terms: 24,
        population_tolerance: 1e-8,
        integral_tolerance_population_s: 1e-8,
    }
}

fn fed_exact_at_one() -> f64 {
    // Integrating factor for A'=10-[ln(2)+1+t] A, A(0)=0.
    // Independent composite Simpson quadrature; it does not call a propagator.
    let a = core::f64::consts::LN_2 + 1.0;
    let n = 4096;
    let h = 1.0 / n as f64;
    let mut integral = 1.0 + (a + 0.5).exp();
    for i in 1..n {
        let t = i as f64 * h;
        integral += if i % 2 == 0 { 2.0 } else { 4.0 } * (a * t + 0.5 * t * t).exp();
    }
    10.0 * (-a - 0.5).exp() * integral * h / 3.0
}

#[test]
fn g1_neutron_driven_supply_refines_birth_times_shutdown_and_restart() {
    let initial = thermal(
        vec![cell(1, 1, 2.0, 1200.0, 0.0)],
        &[(MaterialId(1), &HOT)],
        vec![],
    );
    let heating = KineticHeating::new(&initial, vec![fed_kinetics()], fed_budget()).unwrap();
    let mut control = adaptive(2.0, 1e-4);
    control.tolerances.nuclei = 1e-5;
    control.tolerances.nuclear_energy_j = 2e-5;
    let (heated, evolved, report) = heating
        .advance_adaptive_to(&initial, 2.0, control, || Ok(()))
        .unwrap();
    let source = evolved.source(CellId(1)).unwrap();
    let expected_parent = fed_exact_at_one()
        * (-(core::f64::consts::LN_2 + 2.0) * 0.5).exp()
        * (-core::f64::consts::LN_2 * 0.5).exp();
    close(
        source.irradiation.inventory().populations()[0],
        expected_parent,
        1e-4,
    );
    close(
        source.irradiation.inventory().populations().iter().sum(),
        10.0,
        1e-11,
    );
    close(evolved.audit().produced_nuclei, 10.0, 1e-12);
    close(
        evolved.audit().released.deposited_j,
        2.0 * (10.0 - expected_parent),
        2e-4,
    );
    close(
        heated.temperature_k(CellId(1)).unwrap(),
        300.0 + (10.0 - expected_parent) / 4.0,
        3e-5,
    );
    close(source.kinetics.neutron_population(), 2.0, 1e-13);
    assert!(
        evolved
            .last_step(CellId(1))
            .unwrap()
            .irradiation
            .burnup
            .decays
            > 0.0
    );
    assert_eq!(
        evolved
            .last_step(CellId(1))
            .unwrap()
            .irradiation
            .produced_nuclei,
        [0.0, 0.0]
    );
    assert!(report.native_interval_cuts >= 2);
    assert!(report.accuracy_rejections > 0);
    assert_eq!(evolved.audit().substeps, report.thermal.substeps);

    let (middle, history, first) = heating
        .advance_adaptive_to(&initial, 1.0, control, || Ok(()))
        .unwrap();
    control.initial_dt_s = first.suggested_dt_s;
    let (restarted, after, second) = history
        .advance_adaptive_to(&middle, 2.0, control, || Ok(()))
        .unwrap();
    assert_eq!(
        after
            .source(CellId(1))
            .unwrap()
            .irradiation
            .inventory()
            .populations(),
        source.irradiation.inventory().populations()
    );
    assert_eq!(restarted.cells(), heated.cells());
    close(
        history.audit().released.deposited_j + after.audit().released.deposited_j,
        evolved.audit().released.deposited_j,
        1e-12,
    );
    assert_eq!(
        first.thermal.substeps + second.thermal.substeps,
        report.thermal.substeps
    );
}

#[test]
fn g3_temperature_feedback_accuracy_improves_an_independent_invariant() {
    let initial = thermal(
        vec![cell(1, 1, 2.0, 1200.0, 0.0)],
        &[(MaterialId(1), &HOT)],
        vec![],
    );
    let source = KineticCellSource {
        irradiation: irradiated_source(1000.0, false),
        kinetics: PointKineticsState::new(
            Arc::new(PointKineticsModel::new("manufactured-feedback".into(), 1.0, vec![]).unwrap()),
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
                "manufactured-temperature-worth".into(),
                300.0,
                290.0,
                1000.0,
                TemperatureFeedbackLaw::Linear {
                    coefficient_per_kelvin: -0.02,
                },
            )
            .unwrap(),
        ],
    };
    let heating = KineticHeating::new(
        &initial,
        vec![source],
        KineticHeatingBudget {
            irradiation: irradiation_budget(),
            kinetics: neutron_budget(),
            max_feedback_terms: 2,
        },
    )
    .unwrap();
    let mut errors = Vec::new();
    for amplitude in [2e-3, 2e-5] {
        let mut control = adaptive(1.0, 1e4);
        control.tolerances.nuclei = 1e4;
        control.tolerances.nuclear_energy_j = 1e4;
        control.tolerances.neutron_population = amplitude;
        let (heated, evolved, report) = heating
            .advance_adaptive_to(&initial, 1.0, control, || Ok(()))
            .unwrap();
        let source = evolved.source(CellId(1)).unwrap();
        let i = -(source.irradiation.inventories[0].remaining_parent_nuclei() / 100.0).ln() / 0.2;
        let invariant = 1.0 + 0.05 * i - 0.02 * (1015.0 / 8.0) * (i + (-0.2 * i).exp_m1() / 0.2);
        errors.push((source.kinetics.neutron_population() - invariant).abs());
        assert!(report.accuracy_rejections > 0);
        assert!(report.maximum_accepted_normalized_error <= 1.0);
        assert!(source.kinetics.neutron_population() < 1.0);
        close(
            heated.energy_j() - initial.energy_j(),
            evolved.audit().released.deposited_j,
            1e-8,
        );
    }
    assert!(
        errors[1] < errors[0] / 5.0,
        "feedback invariant errors {errors:?}"
    );
    assert!(errors[1] < 0.001, "feedback invariant errors {errors:?}");
}

struct LatentPlateau;
impl ThermalLaw for LatentPlateau {
    fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
        if h.is_finite() {
            Ok(300.0)
        } else {
            Err(DepositionError::MaterialDomain)
        }
    }
    fn maximum_temperature_slope(&self) -> f64 {
        0.0
    }
}
static PLATEAU: LatentPlateau = LatentPlateau;

#[test]
fn g1_existing_decay_cells_use_the_shared_adaptive_owner_without_production() {
    use crate::decay::thermal::DecayCellSource;
    use crate::decay::{DecayBranch, DecayBudget, DecayInventory, DecayModel, DecayNuclide};
    let mut initial_cell = cell(1, 1, 2.0, 1200.0, 0.0);
    initial_cell.boundary.htc_w_m2_k = 2.0;
    let initial = thermal(vec![initial_cell], &[(MaterialId(1), &HOT)], vec![]);
    let limits = fs_la::decay::DecayLimits::default();
    let model = Arc::new(
        DecayModel::new(
            "manufactured-pure-decay-refinement".into(),
            vec![
                DecayNuclide {
                    id: "parent".into(),
                    law: DecayLawV1::from_half_life(fs_qty::Time::new(
                        2.0 * core::f64::consts::LN_2,
                    ))
                    .unwrap(),
                },
                DecayNuclide {
                    id: "daughter".into(),
                    law: DecayLawV1::stable(),
                },
            ],
            vec![DecayBranch {
                parent: 0,
                daughter: 1,
                probability: 1.0,
                energy_per_decay: EnergyChannels {
                    deposited_j: 2.0,
                    radiation_j: 1.0,
                    escaped_j: 0.0,
                },
            }],
            limits,
        )
        .unwrap(),
    );
    let source = DecayCellSource {
        cell_id: CellId(1),
        material_id: MaterialId(1),
        inventory: DecayInventory::new(model, vec![10.0, 0.0]).unwrap(),
    };
    let heating = NuclearHeating::new(
        &initial,
        vec![source],
        DecayBudget {
            propagation: limits,
            nuclei_tolerance: 1e-8,
            energy_tolerance_j: 1e-8,
        },
    )
    .unwrap();
    let (heated, accepted, report) = heating
        .advance_adaptive_to(&initial, 2.0, adaptive(2.0, 1e-5), || Ok(()))
        .unwrap();
    let remaining = 10.0 * (-1.0_f64).exp();
    let expected_temperature = 300.0 + 1.25 * ((-1.0_f64).exp() - (-0.5_f64).exp()) / (0.25 - 0.5);
    close(
        accepted.source(CellId(1)).unwrap().inventory.populations()[0],
        remaining,
        2e-12,
    );
    close(accepted.audit().decays, 10.0 - remaining, 2e-12);
    close(
        accepted.audit().released.deposited_j,
        2.0 * (10.0 - remaining),
        2e-12,
    );
    close(
        heated.temperature_k(CellId(1)).unwrap(),
        expected_temperature,
        0.002,
    );
    assert!(report.accuracy_rejections > 0);
    assert_eq!(accepted.audit().substeps, report.thermal.substeps);
}

#[test]
fn g3_latent_storage_and_reference_shift_are_visible_without_temperature_error() {
    let mut answers = Vec::new();
    let mut trials = Vec::new();
    for datum in [0.0, 1_048_576.0] {
        let initial = thermal(
            vec![cell(1, 1, 2.0, datum, 0.0)],
            &[(MaterialId(1), &PLATEAU)],
            vec![],
        );
        let heating = KineticHeating::new(&initial, vec![fed_kinetics()], fed_budget()).unwrap();
        let mut control = adaptive(1.0, 1e4);
        control.step.energy_tolerance_j = 1e-4;
        control.tolerances.specific_energy_j_kg = 1e-4;
        control.tolerances.nuclei = 1e4;
        control.tolerances.nuclear_energy_j = 1e4;
        control.tolerances.neutron_population = 1e4;
        let (heated, _, report) = heating
            .advance_adaptive_to(&initial, 1.0, control, || Ok(()))
            .unwrap();
        assert_eq!(heated.temperature_k(CellId(1)).unwrap(), 300.0);
        assert_eq!(report.last_error.unwrap().temperature_difference_k, 0.0);
        assert!(
            report.accuracy_rejections > 0,
            "latent h must refine even at identical temperatures"
        );
        answers.push(heated.cells()[&CellId(1)].specific_h_j_kg - datum);
        trials.push(report.comparisons);
    }
    close(answers[0], answers[1], 1e-7);
    close(answers[0], 10.0 - fed_exact_at_one(), 0.002);
    assert_eq!(
        trials[0], trials[1],
        "caloric datum must not alter the time-error scale"
    );
}

#[test]
fn g4_adaptive_limits_and_late_callback_preserve_original_history() {
    let (initial, heating) = cooling();
    let mut control = adaptive(4.0, 1e-10);
    control.maximum_trials = 1;
    assert!(matches!(
        heating.advance_adaptive_to(&initial, 4.0, control, || Ok(())),
        Err(NuclearAdaptiveError::TrialBudget)
    ));
    control.maximum_trials = 1000;
    control.minimum_dt_s = 4.0;
    assert!(matches!(
        heating.advance_adaptive_to(&initial, 4.0, control, || Ok(())),
        Err(NuclearAdaptiveError::TimeAccuracy { .. })
    ));
    control = adaptive(4.0, 0.001);
    control.maximum_substeps = 1;
    assert!(matches!(
        heating.advance_adaptive_to(&initial, 4.0, control, || Ok(())),
        Err(NuclearAdaptiveError::SubstepBudget)
    ));
    control.maximum_substeps = 1000;
    let mut calls = 0;
    let result = heating.advance_adaptive_to(&initial, 4.0, control, || {
        calls += 1;
        if calls == 500 {
            Err(DepositionError::InvalidInput("caller stop"))
        } else {
            Ok(())
        }
    });
    assert!(matches!(
        result,
        Err(NuclearAdaptiveError::Source(FissionThermalError::Thermal(
            DepositionError::InvalidInput("caller stop")
        )))
    ));
    assert_eq!(calls, 500);
    assert_eq!(heating.time_s(), 0.0);
    assert_eq!(heating.audit(), FissionHeatingAudit::default());
    assert_eq!(
        heating.source(CellId(1)).unwrap().state.delayed_energy_j(),
        [40.0]
    );
    assert_eq!(initial.cells()[&CellId(1)].specific_h_j_kg, 1600.0);

    let limited = DepositionSimulation::new(
        vec![TimedEvent {
            time_s: 0.0,
            event: Event::Birth(cell(1, 1, 2.0, 1200.0, 0.0)),
        }],
        &[(MaterialId(1), &HOT as &dyn ThermalLaw)],
        DepositionBudget {
            max_cells: 1,
            max_links: 0,
            max_events: 1,
            max_substeps: 1,
        },
    )
    .unwrap()
    .advance_to(0.0, controls(1.0), || Ok(()))
    .unwrap()
    .0;
    let limited_source =
        FissionHeating::new(&limited, vec![source(1, 1, 0.0, &[], false)], budget()).unwrap();
    assert!(
        matches!(
            limited_source.advance_adaptive_to(&limited, 1.0, adaptive(1.0, 1.0), || Ok(())),
            Err(NuclearAdaptiveError::SubstepBudget)
        ),
        "original thermal budget must remain a whole-call cap"
    );
}
