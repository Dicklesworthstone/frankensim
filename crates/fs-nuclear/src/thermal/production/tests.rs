use super::*;
use crate::decay::production::DecayProductionChange;
use crate::decay::{DecayBranch, DecayModel, DecayNuclide};
use fs_qty::Time;
use fs_thermal::enthalpy::{
    Boundary, CellGeometry, DepositedCell, DepositionBudget, DepositionError, DepositionObserver,
    DepositionSimulation, Event, StepConfig, ThermalLaw, TimedEvent,
};
use fs_thermochem::decay::DecayLawV1;

// Manufactured calorics and production laws isolate the source/thermal
// equations. No isotope dataset, fission-yield calibration or source history
// is inferred from this receiver material.
struct LinearLaw {
    identity: u8,
    maximum_k: f64,
}

impl ThermalLaw for LinearLaw {
    fn caloric_identity(&self) -> Option<[u8; 32]> {
        Some([self.identity; 32])
    }
    fn temperature_k(&self, energy: f64) -> Result<f64, DepositionError> {
        let temperature = energy / 4.0;
        if temperature.is_finite() && temperature > 0.0 && temperature <= self.maximum_k {
            Ok(temperature)
        } else {
            Err(DepositionError::MaterialDomain)
        }
    }
    fn maximum_temperature_slope(&self) -> f64 {
        0.25
    }
}

static NORMAL: LinearLaw = LinearLaw {
    identity: 1,
    maximum_k: 1e6,
};
static LIMITED: LinearLaw = LinearLaw {
    identity: 2,
    maximum_k: 301.0,
};

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:?} != {expected:?}, tolerance {tolerance}"
    );
}

fn controls(max_dt_s: f64) -> StepConfig {
    StepConfig {
        max_dt_s,
        energy_tolerance_j: 1e-8,
    }
}

fn budget() -> DecayBudget {
    DecayBudget {
        propagation: DecayLimits::default(),
        nuclei_tolerance: 1e-8,
        energy_tolerance_j: 1e-8,
    }
}

fn energy(deposited_j: f64, radiation_j: f64, escaped_j: f64) -> EnergyChannels {
    EnergyChannels {
        deposited_j,
        radiation_j,
        escaped_j,
    }
}

fn thermal(law: &'static dyn ThermalLaw, authored_w: f64) -> DepositionSimulation<'static> {
    let cell = DepositedCell {
        id: CellId(1),
        segment_id: 11,
        body_id: 101,
        material_id: MaterialId(1),
        mass_kg: 2.0,
        specific_h_j_kg: 1200.0,
        geometry: CellGeometry {
            start_m: [0.0, 0.0, 0.0],
            end_m: [0.0, 1.0, 0.0],
            volume_m3: 0.01,
            surface_area_m2: 1.0,
        },
        boundary: Boundary {
            ambient_k: 300.0,
            htc_w_m2_k: 0.0,
            reservoir: None,
            source_w: authored_w,
        },
    };
    DepositionSimulation::new(
        vec![TimedEvent {
            time_s: 0.0,
            event: Event::Birth(cell),
        }],
        &[(MaterialId(1), law)],
        DepositionBudget {
            max_cells: 2,
            max_links: 2,
            max_events: 8,
            max_substeps: 10_000,
        },
    )
    .unwrap()
    .advance_to(0.0, controls(1.0), || Ok(()))
    .unwrap()
    .0
}

fn model(half_life_s: f64, yield_j: EnergyChannels) -> Arc<DecayModel> {
    Arc::new(
        DecayModel::new(
            "manufactured-radioactive-production-v1".into(),
            vec![
                DecayNuclide {
                    id: "parent".into(),
                    law: DecayLawV1::from_half_life(Time::new(half_life_s)).unwrap(),
                },
                DecayNuclide {
                    id: "stable".into(),
                    law: DecayLawV1::stable(),
                },
            ],
            vec![DecayBranch {
                parent: 0,
                daughter: 1,
                probability: 1.0,
                energy_per_decay: yield_j,
            }],
            DecayLimits::default(),
        )
        .unwrap(),
    )
}

fn source(
    inventory: DecayInventory,
    rates: Vec<f64>,
    changes: Vec<DecayProductionChange>,
) -> DecayProductionCellSource {
    let schedule = Arc::new(
        DecayProductionSchedule::new(Arc::clone(inventory.model()), rates, changes, 16).unwrap(),
    );
    DecayProductionCellSource::new(CellId(1), MaterialId(1), inventory, schedule).unwrap()
}

#[test]
fn g1_production_shutdown_clips_the_native_step_and_keeps_residual_decay_heat() {
    let initial = thermal(&NORMAL, 3.0);
    let law = model(core::f64::consts::LN_2, energy(2.0, 1.0, 3.0));
    let supplied = source(
        DecayInventory::new(law, vec![0.0, 0.0]).unwrap(),
        vec![10.0, 0.0],
        vec![DecayProductionChange {
            time_s: 1.0,
            rates_nuclei_s: vec![0.0, 0.0],
        }],
    );
    assert_eq!(supplied.schedule().rates_at(1.0).unwrap(), [0.0, 0.0]);
    let heating = DecayProductionHeating::new(&initial, vec![supplied], budget()).unwrap();
    // The maximum step spans shutdown; the shared transition hook must cut it.
    let (heated, evolved, audit) = heating
        .advance_to(&initial, 3.0, controls(10.0), || Ok(()))
        .unwrap();
    let surviving = 10.0 * -(-1.0_f64).exp_m1() * (-2.0_f64).exp();
    let transformations = 10.0 - surviving;
    let inventory = evolved.source(CellId(1)).unwrap().inventory();
    close(inventory.populations()[0], surviving, 1e-9);
    close(inventory.populations()[1], transformations, 1e-9);
    assert_eq!(inventory.cumulative_produced_nuclei(), [10.0, 0.0]);
    close(evolved.audit().produced_nuclei, 10.0, 1e-12);
    close(evolved.audit().decays, transformations, 1e-9);
    close(
        evolved.audit().released.deposited_j,
        2.0 * transformations,
        1e-8,
    );
    close(evolved.audit().released.radiation_j, transformations, 1e-9);
    close(
        evolved.audit().released.escaped_j,
        3.0 * transformations,
        1e-8,
    );
    close(audit.external_heat_j, 2.0 * transformations + 9.0, 1e-8);
    close(
        heated.temperature_k(CellId(1)).unwrap(),
        300.0 + (2.0 * transformations + 9.0) / 8.0,
        1e-9,
    );
    assert_eq!(audit.substeps, 2);
    assert_eq!(evolved.audit().substeps, 2);
    assert_eq!(evolved.last_step(CellId(1)).unwrap().decay.duration_s, 2.0);
    assert_eq!(heated.cells()[&CellId(1)].boundary.source_w, 3.0);
    assert_eq!(heated.cells()[&CellId(1)].mass_kg, 2.0);

    let (later, continued, later_audit) = evolved
        .advance_to(&heated, 4.0, controls(10.0), || Ok(()))
        .unwrap();
    let late_decays = surviving * -(-1.0_f64).exp_m1();
    assert_eq!(continued.audit().produced_nuclei, 0.0);
    assert!(continued.audit().released.deposited_j > 0.0);
    close(
        continued.audit().released.deposited_j,
        2.0 * late_decays,
        1e-8,
    );
    close(later_audit.external_heat_j, 2.0 * late_decays + 3.0, 1e-8);
    assert!(later.temperature_k(CellId(1)).unwrap() > heated.temperature_k(CellId(1)).unwrap());
    assert_eq!(initial.time_s(), 0.0);
    assert_eq!(
        heating.source(CellId(1)).unwrap().inventory().populations(),
        [0.0, 0.0]
    );
}

#[test]
fn g1_continuous_production_includes_same_interval_daughter_ingrowth_and_heat() {
    let rate = DecayLawV1::from_half_life(Time::new(core::f64::consts::LN_2)).unwrap();
    let law = Arc::new(
        DecayModel::new(
            "manufactured-fed-equal-rate-chain-v1".into(),
            vec![
                DecayNuclide {
                    id: "parent".into(),
                    law: rate,
                },
                DecayNuclide {
                    id: "daughter".into(),
                    law: rate,
                },
                DecayNuclide {
                    id: "stable".into(),
                    law: DecayLawV1::stable(),
                },
            ],
            vec![
                DecayBranch {
                    parent: 0,
                    daughter: 1,
                    probability: 1.0,
                    energy_per_decay: energy(2.0, 1.0, 0.0),
                },
                DecayBranch {
                    parent: 1,
                    daughter: 2,
                    probability: 1.0,
                    energy_per_decay: energy(4.0, 0.0, 2.0),
                },
            ],
            DecayLimits::default(),
        )
        .unwrap(),
    );
    let initial = thermal(&NORMAL, 0.0);
    let heating = DecayProductionHeating::new(
        &initial,
        vec![source(
            DecayInventory::new(law, vec![0.0; 3]).unwrap(),
            vec![100.0, 0.0, 0.0],
            vec![],
        )],
        budget(),
    )
    .unwrap();
    let (heated, evolved, _) = heating
        .advance_to(&initial, 1.0, controls(5.0), || Ok(()))
        .unwrap();
    // Convolution of a unit-rate equal-lambda chain with a uniform 100/s input.
    let exponential = (-1.0_f64).exp();
    let parent = 100.0 * (1.0 - exponential);
    let daughter = 100.0 * (1.0 - 2.0 * exponential);
    let terminal = 100.0 * (-1.0 + 3.0 * exponential);
    let parent_decays = 100.0 - parent;
    let deposited = 2.0 * parent_decays + 4.0 * terminal;
    let populations = evolved.source(CellId(1)).unwrap().inventory().populations();
    for (&actual, expected) in populations.iter().zip([parent, daughter, terminal]) {
        close(actual, expected, 3e-9);
    }
    let report = evolved.last_step(CellId(1)).unwrap();
    close(report.decay.decays[1], terminal, 3e-9);
    close(report.decay.released.deposited_j, deposited, 1e-8);
    close(report.decay.released.radiation_j, parent_decays, 3e-9);
    close(report.decay.released.escaped_j, 2.0 * terminal, 3e-9);
    close(heated.energy_j() - initial.energy_j(), deposited, 1e-8);
    assert_eq!(report.produced_nuclei, [100.0, 0.0, 0.0]);

    let (refined, refined_source, _) = heating
        .advance_to(&initial, 1.0, controls(0.125), || Ok(()))
        .unwrap();
    close(refined.energy_j(), heated.energy_j(), 1e-8);
    for (&actual, &expected) in refined_source
        .source(CellId(1))
        .unwrap()
        .inventory()
        .populations()
        .iter()
        .zip(populations)
    {
        close(actual, expected, 3e-9);
    }
}

#[test]
fn g3_rate_changes_reuse_preparation_and_restart_with_an_independent_fission_owner() {
    use crate::thermal::{FissionCellSource, FissionHeating, FissionRateSchedule};
    use crate::{FissionBudget, FissionModel, FissionState};
    use fs_thermal::observers::CoupledObservers;

    let initial = thermal(&NORMAL, 3.0);
    let decay = DecayProductionHeating::new(
        &initial,
        vec![source(
            DecayInventory::new(
                model(core::f64::consts::LN_2, energy(2.0, 1.0, 3.0)),
                vec![0.0, 0.0],
            )
            .unwrap(),
            vec![10.0, 0.0],
            vec![DecayProductionChange {
                time_s: 1.0,
                rates_nuclei_s: vec![20.0, 0.0],
            }],
        )],
        budget(),
    )
    .unwrap();
    // Independent supplied production and prescribed fission are deliberately
    // different inventories; this composition does not infer fission yields.
    let fission_model = Arc::new(
        FissionModel::new(
            "manufactured-independent-prompt-fission-v1".into(),
            energy(2.0, 1.0, 0.0),
            vec![],
        )
        .unwrap(),
    );
    let fission = FissionHeating::new(
        &initial,
        vec![FissionCellSource {
            cell_id: CellId(1),
            material_id: MaterialId(1),
            state: FissionState::new(fission_model, 100.0, vec![]).unwrap(),
            rates: FissionRateSchedule::new(5.0, vec![]).unwrap(),
        }],
        FissionBudget {
            max_decay_groups: 1,
            energy_tolerance_j: 1e-8,
            fuel_tolerance_nuclei: 1e-8,
        },
    )
    .unwrap();
    let owners = CoupledObservers::new(decay, fission);
    let (whole_t, whole_s, whole_audit) = initial
        .advance_observed(&owners, 2.0, controls(0.5), || Ok(()))
        .unwrap();
    let (part_t, part_s, _) = initial
        .advance_observed(&owners, 1.0, controls(0.5), || Ok(()))
        .unwrap();
    let cached = Arc::clone(
        &part_s
            .first
            .source(CellId(1))
            .unwrap()
            .prepared
            .as_ref()
            .unwrap()
            .plan,
    );
    let (split_t, split_s, _) = part_t
        .advance_observed(&part_s, 2.0, controls(0.5), || Ok(()))
        .unwrap();
    assert!(
        Arc::ptr_eq(
            &cached,
            &split_s
                .first
                .source(CellId(1))
                .unwrap()
                .prepared
                .as_ref()
                .unwrap()
                .plan,
        ),
        "a production-rate change does not change the prepared propagation map"
    );
    assert_eq!(whole_t.cells(), split_t.cells());
    let whole = whole_s.first.source(CellId(1)).unwrap().inventory();
    let split = split_s.first.source(CellId(1)).unwrap().inventory();
    assert_eq!(whole.populations(), split.populations());
    assert_eq!(
        whole.cumulative_produced_nuclei(),
        split.cumulative_produced_nuclei()
    );
    assert_eq!(whole.cumulative_released(), split.cumulative_released());
    assert_eq!(whole.cumulative_produced_nuclei(), [30.0, 0.0]);
    close(whole_s.first.audit().produced_nuclei, 30.0, 1e-12);
    assert_eq!(
        whole_s
            .second
            .source(CellId(1))
            .unwrap()
            .state
            .remaining_fuel_nuclei(),
        90.0
    );
    let surviving = 10.0 * -(-1.0_f64).exp_m1() * (-1.0_f64).exp() + 20.0 * -(-1.0_f64).exp_m1();
    close(
        whole_audit.external_heat_j,
        2.0 * (30.0 - surviving) + 20.0 + 6.0,
        1e-8,
    );
    assert_eq!(whole_audit.substeps, 4);
}

#[test]
fn g4_production_thermal_refusal_and_checkpoint_errors_publish_neither_history() {
    let initial = thermal(&LIMITED, 0.0);
    let heating = DecayProductionHeating::new(
        &initial,
        vec![source(
            DecayInventory::new(
                model(core::f64::consts::LN_2, energy(2.0, 1.0, 3.0)),
                vec![0.0, 0.0],
            )
            .unwrap(),
            vec![100.0, 0.0],
            vec![],
        )],
        budget(),
    )
    .unwrap();
    assert!(matches!(
        heating.advance_to(&initial, 1.0, controls(10.0), || Ok(())),
        Err(FissionThermalError::Thermal(
            DepositionError::MaterialDomain
        ))
    ));
    assert_eq!(initial.time_s(), 0.0);
    assert_eq!(
        heating.source(CellId(1)).unwrap().inventory().populations(),
        [0.0, 0.0]
    );
    assert_eq!(
        heating
            .source(CellId(1))
            .unwrap()
            .inventory()
            .cumulative_produced_nuclei(),
        [0.0, 0.0]
    );

    let mut polls = 0;
    let expected = heating
        .advance_to(&initial, 0.1, controls(0.025), || {
            polls += 1;
            Ok(())
        })
        .unwrap();
    for cutoff in [polls / 2, polls] {
        let mut calls = 0;
        let result = heating.advance_to(&initial, 0.1, controls(0.025), || {
            calls += 1;
            if calls == cutoff {
                Err(DepositionError::Cancelled)
            } else {
                Ok(())
            }
        });
        assert!(matches!(
            result,
            Err(FissionThermalError::Thermal(DepositionError::Cancelled))
        ));
        assert_eq!(
            heating.source(CellId(1)).unwrap().inventory().populations(),
            [0.0, 0.0]
        );
        assert_eq!(
            heating.cumulative_released().unwrap(),
            EnergyChannels::default()
        );
    }
    let retry = heating
        .advance_to(&initial, 0.1, controls(0.025), || Ok(()))
        .unwrap();
    assert_eq!(retry.0.cells(), expected.0.cells());
    assert_eq!(
        retry.1.source(CellId(1)).unwrap().inventory(),
        expected.1.source(CellId(1)).unwrap().inventory()
    );
    assert_eq!(retry.2, expected.2);

    let mut candidate = heating.clone();
    let mut calls = 0;
    let refusal = candidate.step_sources(&initial, 0.0, 0.01, 0.01, &mut || {
        calls += 1;
        if calls > 3 {
            Err(DepositionError::Budget("caller checkpoint"))
        } else {
            Ok(())
        }
    });
    assert!(matches!(
        refusal,
        Err(FissionThermalError::Thermal(DepositionError::Budget(
            "caller checkpoint"
        )))
    ));
    assert!(candidate.pending.is_none());
    assert_eq!(
        candidate
            .source(CellId(1))
            .unwrap()
            .inventory()
            .populations(),
        [0.0, 0.0]
    );
}

#[test]
fn g0_unresolved_production_shares_one_population_allowance_for_the_entire_advance() {
    let law = Arc::new(
        DecayModel::new(
            "manufactured-stable-production-v1".into(),
            vec![DecayNuclide {
                id: "stable".into(),
                law: DecayLawV1::stable(),
            }],
            vec![],
            DecayLimits::default(),
        )
        .unwrap(),
    );
    let supplied = source(
        DecayInventory::new(law, vec![1e16]).unwrap(),
        vec![0.01],
        vec![],
    );
    let initial = thermal(&NORMAL, 0.0);
    // A 0.01 birth is resolved in the input/account but below one ULP of the
    // stored 1e16 population. Measure one interval's noncancelling owner account;
    // the numerical owner may include both aggregate and per-node defects.
    let mut loose = budget();
    loose.nuclei_tolerance = 1.0;
    let probe = DecayProductionHeating::new(&initial, vec![supplied.clone()], loose).unwrap();
    let (_, once, _) = probe
        .advance_to(&initial, 1.0, controls(1.0), || Ok(()))
        .unwrap();
    let one_defect = once.audit().population_roundoff_nuclei;
    assert!(one_defect >= 0.01 && one_defect < 0.1);
    close(once.audit().produced_nuclei, 0.01, 1e-12);
    assert_eq!(
        once.source(CellId(1)).unwrap().inventory().populations(),
        [1e16]
    );

    let mut tight = budget();
    tight.nuclei_tolerance = 1.5 * one_defect;
    let heating = DecayProductionHeating::new(&initial, vec![supplied], tight).unwrap();
    heating
        .advance_to(&initial, 1.0, controls(1.0), || Ok(()))
        .unwrap();
    assert!(matches!(
        heating.advance_to(&initial, 2.0, controls(1.0), || Ok(())),
        Err(FissionThermalError::NucleiBudget { .. })
    ));
    assert_eq!(
        heating.source(CellId(1)).unwrap().inventory().populations(),
        [1e16]
    );
    assert_eq!(
        heating
            .source(CellId(1))
            .unwrap()
            .inventory()
            .cumulative_produced_nuclei(),
        [0.0]
    );
    assert_eq!(initial.time_s(), 0.0);
}

#[test]
fn g0_zero_production_preserves_closed_decay_and_refusals_retain_their_types() {
    use crate::decay::thermal::{DecayCellSource, DecayHeating};

    let initial = thermal(&NORMAL, 0.0);
    let law = model(1.0, energy(8.0, 1.0, 3.0));
    let inventory = DecayInventory::new(Arc::clone(&law), vec![100.0, 0.0]).unwrap();
    let supplied = source(inventory.clone(), vec![0.0, 0.0], vec![]);
    let production =
        DecayProductionHeating::new(&initial, vec![supplied.clone()], budget()).unwrap();
    let closed = DecayHeating::new(
        &initial,
        vec![DecayCellSource {
            cell_id: CellId(1),
            material_id: MaterialId(1),
            inventory: inventory.clone(),
        }],
        budget(),
    )
    .unwrap();
    let (production_t, production_s, production_audit) = production
        .advance_to(&initial, 2.0, controls(0.25), || Ok(()))
        .unwrap();
    let (closed_t, closed_s, closed_audit) = closed
        .advance_to(&initial, 2.0, controls(0.25), || Ok(()))
        .unwrap();
    close(production_t.energy_j(), closed_t.energy_j(), 1e-8);
    close(
        production_audit.external_heat_j,
        closed_audit.external_heat_j,
        1e-8,
    );
    for (&actual, &expected) in production_s
        .source(CellId(1))
        .unwrap()
        .inventory()
        .populations()
        .iter()
        .zip(closed_s.source(CellId(1)).unwrap().inventory.populations())
    {
        close(actual, expected, 1e-9);
    }
    assert_eq!(production_s.audit().produced_nuclei, 0.0);
    assert_eq!(
        production_s
            .source(CellId(1))
            .unwrap()
            .inventory()
            .cumulative_produced_nuclei(),
        [0.0, 0.0]
    );
    assert!(matches!(
        initial.power_balance_observed(&production, || Ok(())),
        Err(FissionThermalError::Thermal(DepositionError::InvalidInput(
            _
        )))
    ));

    let mut invalid = budget();
    invalid.propagation.maximum_work = 0;
    assert!(matches!(
        DecayProductionHeating::new(&initial, vec![supplied], invalid),
        Err(FissionThermalError::Decay(_))
    ));
    let different = model(1.0, energy(8.0, 1.0, 3.0));
    let schedule =
        Arc::new(DecayProductionSchedule::new(different, vec![0.0, 0.0], vec![], 16).unwrap());
    assert!(matches!(
        DecayProductionCellSource::new(CellId(1), MaterialId(1), inventory, schedule),
        Err(DecayError::Invalid(_))
    ));
}
