use super::super::{DecayBranch, DecayModel, DecayNuclide};
use super::*;
use fs_la::decay::DecayLimits;
use fs_qty::Time;
use fs_thermal::enthalpy::{
    Boundary, DepositedCell, DepositionBudget, ThermalLaw, ThermalLink, TimedEvent,
};
use fs_thermochem::decay::DecayLawV1;

// Manufactured caloric and nuclear data verify the equations and coupling.
// These coefficients are not an isotope catalog or empirical validation.
struct LinearLaw {
    cp: f64,
    identity: u8,
    maximum_k: f64,
    basis: ThermalStorageBasis,
}

impl ThermalLaw for LinearLaw {
    fn storage_basis(&self) -> ThermalStorageBasis {
        self.basis
    }
    fn caloric_identity(&self) -> Option<[u8; 32]> {
        Some([self.identity; 32])
    }
    fn temperature_k(&self, energy: f64) -> Result<f64, DepositionError> {
        let temperature = energy / self.cp;
        if temperature.is_finite() && temperature > 0.0 && temperature <= self.maximum_k {
            Ok(temperature)
        } else {
            Err(DepositionError::MaterialDomain)
        }
    }
    fn maximum_temperature_slope(&self) -> f64 {
        self.cp.recip()
    }
}

static HOT: LinearLaw = LinearLaw {
    cp: 4.0,
    identity: 1,
    maximum_k: 1e6,
    basis: ThermalStorageBasis::Enthalpy,
};
static COOL: LinearLaw = LinearLaw {
    cp: 2.0,
    identity: 2,
    maximum_k: 1e6,
    basis: ThermalStorageBasis::Enthalpy,
};
static LIMITED: LinearLaw = LinearLaw {
    cp: 4.0,
    identity: 3,
    maximum_k: 301.0,
    basis: ThermalStorageBasis::Enthalpy,
};
static INTERNAL: LinearLaw = LinearLaw {
    cp: 4.0,
    identity: 1,
    maximum_k: 1e6,
    basis: ThermalStorageBasis::InternalEnergy,
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

fn cell(id: u64, material: u64, mass: f64, energy: f64, source: f64) -> DepositedCell {
    DepositedCell {
        id: CellId(id),
        segment_id: id + 10,
        body_id: id + 100,
        material_id: MaterialId(material),
        mass_kg: mass,
        specific_h_j_kg: energy,
        geometry: CellGeometry {
            start_m: [id as f64, 0.0, 0.0],
            end_m: [id as f64, 1.0, 0.0],
            volume_m3: 0.01,
            surface_area_m2: 1.0,
        },
        boundary: Boundary {
            ambient_k: 300.0,
            htc_w_m2_k: 0.0,
            reservoir: None,
            source_w: source,
        },
    }
}

fn thermal(
    cells: Vec<DepositedCell>,
    materials: &[(MaterialId, &'static dyn ThermalLaw)],
    extra_events: Vec<TimedEvent>,
) -> DepositionSimulation<'static> {
    let mut events: Vec<_> = cells
        .into_iter()
        .map(|cell| TimedEvent {
            time_s: 0.0,
            event: Event::Birth(cell),
        })
        .collect();
    events.extend(extra_events);
    DepositionSimulation::new(
        events,
        materials,
        DepositionBudget {
            max_cells: 8,
            max_links: 16,
            max_events: 64,
            max_substeps: 10_000,
        },
    )
    .unwrap()
    .advance_to(0.0, controls(1.0), || Ok(()))
    .unwrap()
    .0
}

fn parent_inventory(population: f64, half_life_s: f64, yield_j: EnergyChannels) -> DecayInventory {
    let model = Arc::new(
        DecayModel::new(
            "manufactured-parent-calorimetry-v1".into(),
            vec![
                DecayNuclide {
                    id: "parent".into(),
                    law: DecayLawV1::from_half_life(Time::new(half_life_s)).unwrap(),
                },
                DecayNuclide {
                    id: "stable-daughter".into(),
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
    );
    DecayInventory::new(model, vec![population, 0.0]).unwrap()
}

fn source(id: u64, material: u64, inventory: DecayInventory) -> DecayCellSource {
    DecayCellSource {
        cell_id: CellId(id),
        material_id: MaterialId(material),
        inventory,
    }
}

#[test]
fn g1_half_life_calorimetry_deposits_only_its_channel_and_adds_authored_heat_once() {
    let initial = thermal(
        vec![cell(1, 1, 2.0, 1200.0, 3.0)],
        &[(MaterialId(1), &HOT)],
        vec![],
    );
    let heating = DecayHeating::new(
        &initial,
        vec![source(
            1,
            1,
            parent_inventory(100.0, 1.0, energy(8.0, 1.0, 3.0)),
        )],
        budget(),
    )
    .unwrap();
    // The final thermal endpoint clips the much larger proposed maximum step.
    let (heated, decayed, audit) = heating
        .advance_to(&initial, 1.0, controls(10.0), || Ok(()))
        .unwrap();
    let inventory = &decayed.source(CellId(1)).unwrap().inventory;
    close(inventory.populations()[0], 50.0, 1e-9);
    close(inventory.populations()[1], 50.0, 1e-9);
    close(decayed.audit().released.deposited_j, 400.0, 1e-8);
    close(decayed.audit().released.radiation_j, 50.0, 1e-9);
    close(decayed.audit().released.escaped_j, 150.0, 1e-9);
    close(decayed.audit().decays, 50.0, 1e-9);
    close(audit.external_heat_j, 403.0, 1e-8);
    close(
        heated.temperature_k(CellId(1)).unwrap(),
        300.0 + 403.0 / 8.0,
        1e-9,
    );
    assert_eq!(decayed.last_step(CellId(1)).unwrap().duration_s, 1.0);
    assert_eq!(audit.substeps, 1);
    assert_eq!(heated.cells()[&CellId(1)].boundary.source_w, 3.0);
    assert_eq!(heated.cells()[&CellId(1)].mass_kg, 2.0);
    assert_eq!(initial.time_s(), 0.0);
    assert_eq!(
        heating.source(CellId(1)).unwrap().inventory.populations(),
        [100.0, 0.0]
    );

    let (later, continued, later_audit) = decayed
        .advance_to(&heated, 2.0, controls(10.0), || Ok(()))
        .unwrap();
    close(
        continued.source(CellId(1)).unwrap().inventory.populations()[0],
        25.0,
        1e-9,
    );
    close(continued.audit().released.deposited_j, 200.0, 1e-8);
    close(
        continued.cumulative_released().unwrap().deposited_j,
        600.0,
        1e-8,
    );
    close(later_audit.external_heat_j, 203.0, 1e-8);
    close(later.energy_j() - initial.energy_j(), 606.0, 1e-8);
}

#[test]
fn g1_newly_produced_daughters_decay_and_heat_during_the_same_interval() {
    let law = DecayLawV1::from_half_life(Time::new(core::f64::consts::LN_2)).unwrap();
    let model = Arc::new(
        DecayModel::new(
            "manufactured-equal-rate-chain-v1".into(),
            vec![
                DecayNuclide {
                    id: "parent".into(),
                    law,
                },
                DecayNuclide {
                    id: "daughter".into(),
                    law,
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
    let initial = thermal(
        vec![cell(1, 1, 2.0, 1200.0, 0.0)],
        &[(MaterialId(1), &HOT)],
        vec![],
    );
    let heating = DecayHeating::new(
        &initial,
        vec![source(
            1,
            1,
            DecayInventory::new(model, vec![100.0, 0.0, 0.0]).unwrap(),
        )],
        budget(),
    )
    .unwrap();
    let (heated, evolved, _) = heating
        .advance_to(&initial, 1.0, controls(5.0), || Ok(()))
        .unwrap();
    let parent = 100.0 * (-1.0_f64).exp();
    let parent_decays = 100.0 - parent;
    let daughter_decays = 100.0 - 2.0 * parent;
    let deposited = 2.0 * parent_decays + 4.0 * daughter_decays;
    let populations = evolved.source(CellId(1)).unwrap().inventory.populations();
    close(populations[0], parent, 3e-9);
    close(populations[1], parent, 3e-9);
    close(populations[2], daughter_decays, 3e-9);
    let report = evolved.last_step(CellId(1)).unwrap();
    close(report.decays[1], daughter_decays, 3e-9);
    close(report.released.deposited_j, deposited, 1e-8);
    close(report.released.radiation_j, parent_decays, 3e-9);
    close(report.released.escaped_j, 2.0 * daughter_decays, 3e-9);
    close(heated.energy_j() - initial.energy_j(), deposited, 1e-8);

    let (refined, refined_source, _) = heating
        .advance_to(&initial, 1.0, controls(0.125), || Ok(()))
        .unwrap();
    close(refined.energy_j(), heated.energy_j(), 1e-8);
    for (&actual, &expected) in refined_source
        .source(CellId(1))
        .unwrap()
        .inventory
        .populations()
        .iter()
        .zip(populations)
    {
        close(actual, expected, 3e-9);
    }
}

#[test]
fn g1_decay_heat_uses_original_contact_exchange_with_first_order_time_refinement() {
    let initial = thermal(
        vec![cell(1, 1, 2.0, 1200.0, 0.0), cell(2, 2, 3.0, 560.0, 0.0)],
        &[(MaterialId(1), &HOT), (MaterialId(2), &COOL)],
        vec![TimedEvent {
            time_s: 0.0,
            event: Event::Contact(ThermalLink {
                id: 7,
                a: CellId(1),
                b: CellId(2),
                area_m2: 0.1,
                conductance_w_k: 1.0,
            }),
        }],
    );
    let heating = DecayHeating::new(
        &initial,
        vec![source(
            1,
            1,
            parent_inventory(100.0, core::f64::consts::LN_2, energy(2.0, 1.0, 1.0)),
        )],
        budget(),
    )
    .unwrap();
    let end = 4.0_f64;
    let released = 200.0 * -(-end).exp_m1();
    let average = (8.0 * 300.0 + 6.0 * 280.0 + released) / 14.0;
    let rate = 7.0 / 24.0;
    let difference =
        20.0 * (-rate * end).exp() + 25.0 * ((-end).exp() - (-rate * end).exp()) / (rate - 1.0);
    let exact_hot = average + 6.0 / 14.0 * difference;
    let exact_cool = average - 8.0 / 14.0 * difference;
    let mut errors = Vec::new();
    for dt in [0.125, 0.0625, 0.03125] {
        let (result, source, audit) = heating
            .advance_to(&initial, end, controls(dt), || Ok(()))
            .unwrap();
        close(audit.external_heat_j, released, 1e-8);
        close(result.energy_j() - initial.energy_j(), released, 1e-8);
        close(source.audit().thermal_credit_j, released, 1e-8);
        close(
            source.audit().released.radiation_j + source.audit().released.escaped_j,
            released,
            1e-8,
        );
        let hot = result.temperature_k(CellId(1)).unwrap();
        let cool = result.temperature_k(CellId(2)).unwrap();
        assert!(cool > 280.0);
        errors.push((hot - exact_hot).abs().max((cool - exact_cool).abs()));
        assert_eq!(result.links(), initial.links());
    }
    assert!(errors[2] < 0.02);
    for pair in errors.windows(2) {
        assert!((1.9..2.1).contains(&(pair[0] / pair[1])), "{errors:?}");
    }
}

#[test]
fn g3_multiple_inventories_restart_on_identical_intervals_and_keep_original_storage_basis() {
    let initial = thermal(
        vec![cell(1, 1, 2.0, 1200.0, 0.0), cell(2, 2, 3.0, 600.0, 0.0)],
        &[(MaterialId(1), &INTERNAL), (MaterialId(2), &COOL)],
        vec![],
    );
    let heating = DecayHeating::new(
        &initial,
        vec![
            source(2, 2, parent_inventory(50.0, 2.0, energy(4.0, 1.0, 1.0))),
            source(1, 1, parent_inventory(100.0, 1.0, energy(8.0, 1.0, 3.0))),
        ],
        budget(),
    )
    .unwrap();
    assert_eq!(
        heating
            .sources()
            .map(|value| value.cell_id)
            .collect::<Vec<_>>(),
        [CellId(1), CellId(2)]
    );
    let (whole_t, whole_s, whole_audit) = heating
        .advance_to(&initial, 2.0, controls(0.5), || Ok(()))
        .unwrap();
    let (part_t, part_s, _) = heating
        .advance_to(&initial, 1.0, controls(0.5), || Ok(()))
        .unwrap();
    let checkpoint_t = part_t.clone();
    let checkpoint_s = part_s.clone();
    let (split_t, split_s, _) = checkpoint_s
        .advance_to(&checkpoint_t, 2.0, controls(0.5), || Ok(()))
        .unwrap();
    assert_eq!(whole_t.cells(), split_t.cells());
    for cell_id in [CellId(1), CellId(2)] {
        let whole = &whole_s.source(cell_id).unwrap().inventory;
        let split = &split_s.source(cell_id).unwrap().inventory;
        assert_eq!(whole.populations(), split.populations());
        assert_eq!(whole.cumulative_released(), split.cumulative_released());
    }
    close(whole_audit.external_heat_j, 700.0, 1e-8);
    assert_eq!(
        whole_t.material_storage_basis(MaterialId(1)).unwrap(),
        ThermalStorageBasis::InternalEnergy
    );
    assert_eq!(whole_s.audit().substeps, 4);

    let (idle_t, idle_s, idle_audit) = whole_s
        .advance_to(&whole_t, 2.0, controls(1.0), || Ok(()))
        .unwrap();
    assert_eq!(idle_t.cells(), whole_t.cells());
    assert_eq!(idle_s.audit(), DecayHeatingAudit::default());
    assert_eq!(idle_audit.substeps, 0);
    assert!(matches!(
        heating.advance_to(&whole_t, 2.0, controls(1.0), || Ok(())),
        Err(DecayThermalError::ClockMismatch { .. })
    ));
}

#[test]
fn g4_thermal_refusal_and_late_checkpoint_failure_roll_back_both_histories() {
    let initial = thermal(
        vec![cell(1, 1, 2.0, 1200.0, 0.0)],
        &[(MaterialId(1), &LIMITED)],
        vec![],
    );
    let heating = DecayHeating::new(
        &initial,
        vec![source(
            1,
            1,
            parent_inventory(100.0, 1.0, energy(8.0, 1.0, 3.0)),
        )],
        budget(),
    )
    .unwrap();
    assert!(matches!(
        heating.advance_to(&initial, 1.0, controls(10.0), || Ok(())),
        Err(DecayThermalError::Thermal(DepositionError::MaterialDomain))
    ));
    assert_eq!(initial.time_s(), 0.0);
    assert_eq!(
        heating.source(CellId(1)).unwrap().inventory.populations(),
        [100.0, 0.0]
    );
    assert_eq!(
        heating.cumulative_released().unwrap(),
        EnergyChannels::default()
    );

    let mut polls = 0;
    let expected = heating
        .advance_to(&initial, 0.01, controls(0.0025), || {
            polls += 1;
            Ok(())
        })
        .unwrap();
    for cutoff in [polls / 2, polls] {
        let mut calls = 0;
        let result = heating.advance_to(&initial, 0.01, controls(0.0025), || {
            calls += 1;
            if calls == cutoff {
                Err(DepositionError::Cancelled)
            } else {
                Ok(())
            }
        });
        assert!(matches!(
            result,
            Err(DecayThermalError::Thermal(DepositionError::Cancelled))
        ));
        assert_eq!(
            heating.source(CellId(1)).unwrap().inventory.populations(),
            [100.0, 0.0]
        );
        assert_eq!(
            heating.cumulative_released().unwrap(),
            EnergyChannels::default()
        );
    }
    let retry = heating
        .advance_to(&initial, 0.01, controls(0.0025), || Ok(()))
        .unwrap();
    assert_eq!(retry.0.cells(), expected.0.cells());
    assert_eq!(
        retry.1.source(CellId(1)).unwrap().inventory.populations(),
        expected
            .1
            .source(CellId(1))
            .unwrap()
            .inventory
            .populations()
    );
    assert_eq!(retry.2, expected.2);

    // The preparation/application checkpoint must retain the caller's actual error.
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
        Err(DecayThermalError::Thermal(DepositionError::Budget(
            "caller checkpoint"
        )))
    ));
    assert!(candidate.pending.is_none());
    assert_eq!(
        candidate.source(CellId(1)).unwrap().inventory.populations(),
        [100.0, 0.0]
    );
}

#[test]
fn g0_unresolved_population_debits_share_one_allowance_across_thermal_substeps() {
    let initial = thermal(
        vec![cell(1, 1, 2.0, 1200.0, 0.0)],
        &[(MaterialId(1), &HOT)],
        vec![],
    );
    // Each second resolves 0.01 expected decays and their daughter/heat credits,
    // but that debit is below one ULP of the 1e16-parent inventory. The explicit
    // allowance admits one interval; renewing it at each substep would hide loss.
    let mut allowance = budget();
    allowance.nuclei_tolerance = 0.015;
    let heating = DecayHeating::new(
        &initial,
        vec![source(
            1,
            1,
            parent_inventory(1e16, core::f64::consts::LN_2 * 1e18, energy(1.0, 0.0, 0.0)),
        )],
        allowance,
    )
    .unwrap();
    let (_, accepted, _) = heating
        .advance_to(&initial, 1.0, controls(1.0), || Ok(()))
        .unwrap();
    close(accepted.audit().population_roundoff_nuclei, 0.01, 1e-12);
    close(accepted.audit().thermal_credit_j, 0.01, 1e-12);
    assert!(matches!(
        heating.advance_to(&initial, 2.0, controls(1.0), || Ok(())),
        Err(DecayThermalError::NucleiBudget { .. })
    ));
    assert_eq!(
        heating.source(CellId(1)).unwrap().inventory.populations(),
        [1e16, 0.0]
    );
    assert_eq!(
        heating.cumulative_released().unwrap(),
        EnergyChannels::default()
    );
    assert_eq!(initial.time_s(), 0.0);
}

#[test]
fn g0_bad_bindings_and_incomplete_instantaneous_balance_are_refused() {
    let initial = thermal(
        vec![cell(1, 1, 2.0, 1200.0, 0.0), cell(2, 2, 3.0, 600.0, 0.0)],
        &[(MaterialId(1), &HOT), (MaterialId(2), &COOL)],
        vec![],
    );
    let inventory = parent_inventory(100.0, 1.0, energy(8.0, 1.0, 3.0));
    assert!(matches!(
        DecayHeating::new(&initial, vec![source(1, 2, inventory.clone())], budget()),
        Err(DecayThermalError::MaterialMismatch { .. })
    ));
    assert!(matches!(
        DecayHeating::new(&initial, vec![source(99, 1, inventory.clone())], budget()),
        Err(DecayThermalError::InactiveCell(CellId(99)))
    ));
    assert!(
        DecayHeating::new(
            &initial,
            vec![
                source(1, 1, inventory.clone()),
                source(1, 1, inventory.clone()),
            ],
            budget()
        )
        .is_err()
    );
    let heating = DecayHeating::new(&initial, vec![source(1, 1, inventory)], budget()).unwrap();
    assert!(matches!(
        initial.power_balance_observed(&heating, || Ok(())),
        Err(DecayThermalError::Thermal(DepositionError::InvalidInput(_)))
    ));
    let changed_mass = thermal(
        vec![cell(1, 1, 3.0, 1200.0, 0.0)],
        &[(MaterialId(1), &HOT)],
        vec![],
    );
    assert!(matches!(
        heating.advance_to(&changed_mass, 0.0, controls(1.0), || Ok(())),
        Err(DecayThermalError::BindingChanged(CellId(1)))
    ));
    let changed_calorics = thermal(
        vec![cell(1, 1, 2.0, 1200.0, 0.0)],
        &[(MaterialId(1), &LIMITED)],
        vec![],
    );
    assert!(matches!(
        heating.advance_to(&changed_calorics, 0.0, controls(1.0), || Ok(())),
        Err(DecayThermalError::BindingChanged(CellId(1)))
    ));
    let changed_basis = thermal(
        vec![cell(1, 1, 2.0, 1200.0, 0.0)],
        &[(MaterialId(1), &INTERNAL)],
        vec![],
    );
    assert!(matches!(
        heating.advance_to(&changed_basis, 0.0, controls(1.0), || Ok(())),
        Err(DecayThermalError::BindingChanged(CellId(1)))
    ));
}

#[test]
fn g3_independent_decay_and_fission_sources_share_one_native_thermal_advance() {
    use crate::thermal::{FissionCellSource, FissionHeating, FissionRateSchedule};
    use crate::{FissionBudget, FissionModel, FissionState};
    use fs_thermal::observers::CoupledObservers;

    let initial = thermal(
        vec![cell(1, 1, 2.0, 1200.0, 3.0)],
        &[(MaterialId(1), &HOT)],
        vec![],
    );
    let decay = DecayHeating::new(
        &initial,
        vec![source(
            1,
            1,
            parent_inventory(100.0, 1.0, energy(8.0, 1.0, 3.0)),
        )],
        budget(),
    )
    .unwrap();
    let model = Arc::new(
        FissionModel::new(
            "manufactured-independent-fission-v1".into(),
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
            state: FissionState::new(model, 100.0, vec![]).unwrap(),
            rates: FissionRateSchedule::new(3.0, vec![]).unwrap(),
        }],
        FissionBudget {
            max_decay_groups: 1,
            energy_tolerance_j: 1e-8,
            fuel_tolerance_nuclei: 1e-8,
        },
    )
    .unwrap();
    let owners = CoupledObservers::new(decay, fission);
    let (heated, accepted, audit) = initial
        .advance_observed(&owners, 1.0, controls(10.0), || Ok(()))
        .unwrap();
    close(audit.external_heat_j, 400.0 + 6.0 + 3.0, 1e-8);
    close(
        heated.energy_j() - initial.energy_j(),
        audit.external_heat_j,
        1e-8,
    );
    close(accepted.first.audit().thermal_credit_j, 400.0, 1e-8);
    close(accepted.second.audit().thermal_credit_j, 6.0, 1e-12);
    close(
        accepted
            .first
            .source(CellId(1))
            .unwrap()
            .inventory
            .populations()[0],
        50.0,
        1e-9,
    );
    assert_eq!(
        accepted
            .second
            .source(CellId(1))
            .unwrap()
            .state
            .remaining_fuel_nuclei(),
        97.0
    );
    assert_eq!(accepted.first.time_s(), 1.0);
    assert_eq!(accepted.second.time_s(), 1.0);
    assert_eq!(audit.substeps, 1);
}
