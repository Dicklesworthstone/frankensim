use super::*;
use crate::decay::activation::CaptureProduct;
use crate::decay::{DecayBranch, DecayBudget, DecayInventory, DecayModel, DecayNuclide};
use crate::irradiation::{
    GroupCrossSections, GroupFlux, IrradiationBudget, IrradiationData, IrradiationModel,
    IrradiationState, NeutronGroups, NuclideId,
};
use crate::kinetics::{KineticsBudget, PointKineticsModel, PointKineticsState};
use crate::thermal::{
    IrradiationCellSource, IrradiationHeating, IrradiationHeatingBudget, KineticCellSource,
    KineticHeating, KineticHeatingBudget, KineticsControl, KineticsControlSchedule,
    NeutronFluxChange,
};
use crate::{FissionBudget, FissionModel};
use fs_la::decay::DecayLimits;
use fs_qty::Time;
use fs_thermal::enthalpy::{
    Boundary, CellGeometry, DepositedCell, DepositionBudget, DepositionError, DepositionSimulation,
    Event, StepConfig, ThermalLaw, TimedEvent,
};
use fs_thermochem::decay::DecayLawV1;
use std::sync::Arc;

// These cross sections, half-lives and energies are manufactured equation data.
// The typed nuclide association checks capture bookkeeping, not isotope calibration.
struct LinearLaw;
impl ThermalLaw for LinearLaw {
    fn caloric_identity(&self) -> Option<[u8; 32]> {
        Some([71; 32])
    }
    fn temperature_k(&self, energy: f64) -> Result<f64, DepositionError> {
        let temperature = energy / 4.0;
        if temperature.is_finite() && temperature > 0.0 && temperature <= 1e6 {
            Ok(temperature)
        } else {
            Err(DepositionError::MaterialDomain)
        }
    }
    fn maximum_temperature_slope(&self) -> f64 {
        0.25
    }
}
static MATERIAL: LinearLaw = LinearLaw;

fn energy(deposited_j: f64, radiation_j: f64, escaped_j: f64) -> EnergyChannels {
    EnergyChannels {
        deposited_j,
        radiation_j,
        escaped_j,
    }
}
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
fn budget() -> CaptureActivationBudget {
    CaptureActivationBudget {
        irradiation: IrradiationBudget {
            max_neutron_groups: 4,
            source: FissionBudget {
                max_decay_groups: 4,
                energy_tolerance_j: 1e-8,
                fuel_tolerance_nuclei: 1e-8,
            },
        },
        decay: DecayBudget {
            propagation: DecayLimits::default(),
            nuclei_tolerance: 1e-8,
            energy_tolerance_j: 1e-8,
        },
    }
}
fn thermal(authored_w: f64) -> DepositionSimulation<'static> {
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
        &[(MaterialId(1), &MATERIAL as &dyn ThermalLaw)],
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
fn source(
    fission_rate: f64,
    capture_rate: f64,
    maximum_k: f64,
    initial_products: f64,
    changes: &[(f64, f64)],
) -> CaptureActivationCellSource {
    // Unit incident flux makes the manufactured cross sections equal to the
    // desired inverse-second rates, with no spectrum-collapse assumption.
    let groups = Arc::new(
        NeutronGroups::new("manufactured-activation-grid-v1".into(), vec![1.0, 2.0]).unwrap(),
    );
    let flux = GroupFlux::new(Arc::clone(&groups), vec![1.0]).unwrap();
    let fission = Arc::new(
        FissionModel::new(
            "manufactured-prompt-fission-v1".into(),
            energy(5.0, 2.0, 1.0),
            vec![],
        )
        .unwrap(),
    );
    let irradiation = Arc::new(
        IrradiationModel::new(IrradiationData {
            data_id: "manufactured-capture-parent-v1".into(),
            nuclide: NuclideId::new(27, 59, 0).unwrap(),
            groups: Arc::clone(&groups),
            cross_sections: vec![GroupCrossSections {
                fission_m2: fission_rate,
                capture_m2: capture_rate,
            }],
            temperature_min_k: 200.0,
            temperature_max_k: maximum_k,
            fission_model: fission,
            capture_energy_per_event: energy(2.0, 1.0, 1.0),
        })
        .unwrap(),
    );
    let parent = IrradiationState::new(irradiation, 100.0, vec![]).unwrap();
    let decay = Arc::new(
        DecayModel::new(
            "manufactured-capture-product-decay-v1".into(),
            vec![
                DecayNuclide {
                    id: "activation-product".into(),
                    law: DecayLawV1::from_half_life(Time::new(core::f64::consts::LN_2)).unwrap(),
                },
                DecayNuclide {
                    id: "stable-product".into(),
                    law: DecayLawV1::stable(),
                },
            ],
            vec![DecayBranch {
                parent: 0,
                daughter: 1,
                probability: 1.0,
                energy_per_decay: energy(4.0, 1.0, 2.0),
            }],
            DecayLimits::default(),
        )
        .unwrap(),
    );
    let products = DecayInventory::new(decay, vec![initial_products, 0.0]).unwrap();
    let state = CaptureActivationState::new(
        parent,
        products,
        CaptureProduct {
            node: 0,
            nuclide: NuclideId::new(27, 60, 0).unwrap(),
            data_id: "manufactured-radiative-capture-binding-v1".into(),
        },
    )
    .unwrap();
    let fluxes = NeutronFluxSchedule::new(
        flux,
        changes
            .iter()
            .map(|&(time_s, value)| NeutronFluxChange {
                time_s,
                flux: GroupFlux::new(Arc::clone(&groups), vec![value]).unwrap(),
            })
            .collect(),
    )
    .unwrap();
    CaptureActivationCellSource {
        cell_id: CellId(1),
        material_id: MaterialId(1),
        state,
        fluxes,
    }
}

#[test]
fn g1_depleted_capture_history_and_shutdown_supply_the_actual_product_heat() {
    let initial = thermal(3.0);
    let heating = CaptureActivationHeating::new(
        &initial,
        vec![source(0.0, 1.0, 1000.0, 0.0, &[(1.0, 0.0)])],
        budget(),
    )
    .unwrap();
    // The equal parent-removal/product-decay rates give X(1)=100 exp(-1).
    // Irradiation then stops, while X cools through two more decay lifetimes.
    let (heated, evolved, audit) = heating
        .advance_to(&initial, 3.0, controls(10.0), || Ok(()))
        .unwrap();
    let captures = 100.0 * -(-1.0_f64).exp_m1();
    let parent = 100.0 * (-1.0_f64).exp();
    let product = 100.0 * (-3.0_f64).exp();
    let product_decays = captures - product;
    let deposited = 2.0 * captures + 4.0 * product_decays;
    let accepted = &evolved.source(CellId(1)).unwrap().state;
    close(accepted.parent().remaining_parent_nuclei(), parent, 1e-9);
    close(accepted.products().populations()[0], product, 1e-9);
    close(accepted.products().populations()[1], product_decays, 1e-9);
    close(
        accepted.products().cumulative_produced_nuclei()[0],
        captures,
        1e-9,
    );
    close(evolved.audit().captures, captures, 1e-9);
    close(evolved.audit().produced_nuclei, captures, 1e-9);
    close(evolved.audit().decays, product_decays, 1e-9);
    assert_eq!(evolved.audit().fissions, 0.0);
    close(evolved.audit().released.deposited_j, deposited, 1e-8);
    close(
        evolved.audit().released.radiation_j,
        captures + product_decays,
        1e-8,
    );
    close(
        evolved.audit().released.escaped_j,
        captures + 2.0 * product_decays,
        1e-8,
    );
    close(audit.external_heat_j, deposited + 9.0, 1e-8);
    close(
        heated.temperature_k(CellId(1)).unwrap(),
        300.0 + (deposited + 9.0) / 8.0,
        1e-9,
    );
    assert_eq!(audit.substeps, 2);
    assert_eq!(heated.cells()[&CellId(1)].boundary.source_w, 3.0);
    assert_eq!(heated.cells()[&CellId(1)].mass_kg, 2.0);

    let (later, continued, later_audit) = evolved
        .advance_to(&heated, 4.0, controls(10.0), || Ok(()))
        .unwrap();
    assert_eq!(continued.audit().captures, 0.0);
    assert_eq!(continued.audit().produced_nuclei, 0.0);
    assert!(continued.audit().released.deposited_j > 0.0);
    let late_heat = 4.0 * product * -(-1.0_f64).exp_m1();
    close(later_audit.external_heat_j, late_heat + 3.0, 1e-8);
    assert_eq!(
        continued
            .source(CellId(1))
            .unwrap()
            .state
            .parent()
            .remaining_parent_nuclei(),
        accepted.parent().remaining_parent_nuclei()
    );
    assert!(later.temperature_k(CellId(1)).unwrap() > heated.temperature_k(CellId(1)).unwrap());
    assert_eq!(initial.time_s(), 0.0);
}

#[test]
fn g1_competing_fission_keeps_the_original_irradiation_partition_and_adds_product_heat_once() {
    let initial = thermal(0.0);
    let supplied = source(1.0, 1.0, 1000.0, 0.0, &[]);
    let irradiation = IrradiationHeating::new(
        &initial,
        vec![IrradiationCellSource {
            cell_id: CellId(1),
            material_id: MaterialId(1),
            inventories: vec![supplied.state.parent().clone()],
            fluxes: supplied.fluxes.clone(),
        }],
        IrradiationHeatingBudget {
            max_isotopes_per_cell: 1,
            source: budget().irradiation,
        },
    )
    .unwrap();
    let activation = CaptureActivationHeating::new(&initial, vec![supplied], budget()).unwrap();
    let (heated, evolved, audit) = activation
        .advance_to(&initial, 1.0, controls(10.0), || Ok(()))
        .unwrap();
    let (baseline_t, baseline_s, baseline_audit) = irradiation
        .advance_to(&initial, 1.0, controls(10.0), || Ok(()))
        .unwrap();
    let reactions = 50.0 * -(-2.0_f64).exp_m1();
    let product = 100.0 * ((-1.0_f64).exp() - (-2.0_f64).exp());
    let product_decays = reactions - product;
    close(evolved.audit().fissions, reactions, 1e-9);
    close(evolved.audit().captures, reactions, 1e-9);
    close(evolved.audit().produced_nuclei, reactions, 1e-9);
    close(
        evolved
            .source(CellId(1))
            .unwrap()
            .state
            .products()
            .populations()[0],
        product,
        1e-9,
    );
    assert_eq!(
        evolved.source(CellId(1)).unwrap().state.parent(),
        &baseline_s.source(CellId(1)).unwrap().inventories[0]
    );
    close(
        audit.external_heat_j - baseline_audit.external_heat_j,
        4.0 * product_decays,
        1e-8,
    );
    close(
        heated.energy_j() - baseline_t.energy_j(),
        4.0 * product_decays,
        1e-8,
    );
    close(
        evolved.audit().released.deposited_j,
        7.0 * reactions + 4.0 * product_decays,
        1e-8,
    );
    close(
        evolved.audit().released.radiation_j,
        3.0 * reactions + product_decays,
        1e-8,
    );
    close(
        evolved.audit().released.escaped_j,
        2.0 * reactions + 2.0 * product_decays,
        1e-8,
    );
}

#[test]
fn g4_temperature_grid_and_checkpoint_refusals_preserve_both_nuclear_histories() {
    let initial = thermal(0.0);
    let supplied = source(0.0, 1.0, 301.0, 0.0, &[]);
    let original_parent = supplied.state.parent().clone();
    let original_products = supplied.state.products().clone();
    let heating =
        CaptureActivationHeating::new(&initial, vec![supplied.clone()], budget()).unwrap();
    // The material chart accepts the temperature, but the nuclear data do not.
    // Endpoint applicability must be checked before publishing either inventory.
    assert!(matches!(
        heating.advance_to(&initial, 1.0, controls(10.0), || Ok(())),
        Err(FissionThermalError::Activation(
            CaptureActivationError::Irradiation(FissionError::InvalidInput(_),)
        ))
    ));
    assert_eq!(
        heating.source(CellId(1)).unwrap().state.parent(),
        &original_parent
    );
    assert_eq!(
        heating.source(CellId(1)).unwrap().state.products(),
        &original_products
    );
    assert_eq!(initial.time_s(), 0.0);

    let mut wrong_grid = supplied.clone();
    let other_groups =
        Arc::new(NeutronGroups::new("different-grid-identity-v1".into(), vec![1.0, 2.0]).unwrap());
    wrong_grid.fluxes =
        NeutronFluxSchedule::new(GroupFlux::new(other_groups, vec![0.0]).unwrap(), vec![]).unwrap();
    assert!(matches!(
        CaptureActivationHeating::new(&initial, vec![wrong_grid], budget()),
        Err(FissionThermalError::Activation(
            CaptureActivationError::Irradiation(FissionError::InvalidInput(_),)
        ))
    ));

    let refusal = supplied.trial_source(0.0, 0.01, 300.0, budget(), &mut || {
        Err(FissionError::Budget("caller checkpoint"))
    });
    assert!(matches!(
        refusal,
        Err(FissionThermalError::Fission(FissionError::Budget(
            "caller checkpoint"
        )))
    ));
    assert_eq!(supplied.state.parent(), &original_parent);
    assert_eq!(supplied.state.products(), &original_products);
}

#[test]
fn g3_capture_cache_survives_thermal_cloning_and_identical_partition_restart() {
    let initial = thermal(0.0);
    let heating =
        CaptureActivationHeating::new(&initial, vec![source(0.2, 0.8, 1000.0, 2.0, &[])], budget())
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
        "unchanged duration/rates reuse the original prepared map: warm {warm_polls}, cold {cold_polls}"
    );
    let checkpoint_t = middle_t.clone();
    let checkpoint_s = middle_s.clone();
    let (split_t, split_s, _) = checkpoint_s
        .advance_to(&checkpoint_t, 1.0, controls(0.25), || Ok(()))
        .unwrap();
    assert_eq!(whole_t.cells(), split_t.cells());
    let whole = &whole_s.source(CellId(1)).unwrap().state;
    let split = &split_s.source(CellId(1)).unwrap().state;
    assert_eq!(whole.parent(), split.parent());
    assert_eq!(whole.products(), split.products());
    assert_eq!(whole.capture_product(), split.capture_product());
    assert_eq!(
        whole_s.cumulative_released().unwrap(),
        split_s.cumulative_released().unwrap()
    );
    assert_eq!(whole_s.time_s(), split_s.time_s());
}

#[test]
fn g1_neutron_amplitude_drives_capture_births_and_heat_on_the_same_transaction() {
    let initial = thermal(3.0);
    let model = Arc::new(
        PointKineticsModel::new(
            "manufactured-constant-neutron-amplitude-v1".into(),
            1.0,
            vec![],
        )
        .unwrap(),
    );
    let supplied = KineticCellSource {
        irradiation: source(0.0, 0.5, 1000.0, 0.0, &[]),
        kinetics: PointKineticsState::new(model, 2.0, vec![]).unwrap(),
        controls: KineticsControlSchedule::new(
            KineticsControl {
                reactivity: 0.0,
                external_source_per_s: 0.0,
            },
            vec![],
        )
        .unwrap(),
        feedback: vec![],
    };
    let original = supplied.clone();
    let allowance = KineticHeatingBudget {
        irradiation: budget(),
        kinetics: KineticsBudget {
            max_delayed_groups: 0,
            max_squarings: 32,
            max_taylor_terms: 24,
            population_tolerance: 1e-8,
            integral_tolerance_population_s: 1e-8,
        },
        max_feedback_terms: 0,
    };
    let heating =
        KineticHeating::<CaptureActivationCellSource>::new(&initial, vec![supplied], allowance)
            .unwrap();
    let mut accepted_polls = 0;
    let (heated, evolved, audit) = heating
        .advance_to(&initial, 1.0, controls(10.0), || {
            accepted_polls += 1;
            Ok(())
        })
        .unwrap();

    // rho=beta=0 keeps n(t)=2. The reference capture rate 0.5 therefore
    // becomes 1/s, equal to the product's decay rate. With P(0)=100,
    // the independent Bateman limit is P(1)=X(1)=100 exp(-1).
    let survivor = 100.0 * (-1.0_f64).exp();
    let captures = 100.0 * -(-1.0_f64).exp_m1();
    let product_decays = captures - survivor;
    let deposited = 2.0 * captures + 4.0 * product_decays;
    let accepted = evolved.source(CellId(1)).unwrap();
    let report = evolved.last_step(CellId(1)).unwrap();
    close(accepted.kinetics.neutron_population(), 2.0, 1e-12);
    close(report.kinetics.integrated_population_s, 2.0, 1e-12);
    close(report.mean_neutron_population, 2.0, 1e-12);
    assert_eq!(report.applied_reactivity, 0.0);
    let inventory = &accepted.irradiation.state;
    close(inventory.parent().remaining_parent_nuclei(), survivor, 1e-9);
    close(inventory.products().populations()[0], survivor, 1e-9);
    close(inventory.products().populations()[1], product_decays, 1e-9);
    close(evolved.audit().captures, captures, 1e-9);
    close(evolved.audit().produced_nuclei, captures, 1e-9);
    close(evolved.audit().decays, product_decays, 1e-9);
    assert_eq!(evolved.audit().fissions, 0.0);
    close(evolved.audit().released.deposited_j, deposited, 1e-8);
    close(
        evolved.audit().released.radiation_j,
        captures + product_decays,
        1e-8,
    );
    close(
        evolved.audit().released.escaped_j,
        captures + 2.0 * product_decays,
        1e-8,
    );
    close(audit.external_heat_j, deposited + 3.0, 1e-8);
    close(
        heated.temperature_k(CellId(1)).unwrap(),
        300.0 + (deposited + 3.0) / 8.0,
        1e-9,
    );
    assert_eq!(audit.substeps, 1);

    // Refuse at the final checkpoint after prospective neutron, activation and
    // thermal work. All three original histories must still be available.
    assert!(accepted_polls > 1);
    let mut refused_polls = 0;
    assert!(
        heating
            .advance_to(&initial, 1.0, controls(10.0), || {
                refused_polls += 1;
                if refused_polls == accepted_polls {
                    Err(DepositionError::Cancelled)
                } else {
                    Ok(())
                }
            })
            .is_err()
    );
    assert_eq!(refused_polls, accepted_polls);
    let retained = heating.source(CellId(1)).unwrap();
    assert_eq!(retained.kinetics, original.kinetics);
    assert_eq!(
        retained.irradiation.state.parent(),
        original.irradiation.state.parent()
    );
    assert_eq!(
        retained.irradiation.state.products(),
        original.irradiation.state.products()
    );
    assert_eq!(heating.time_s(), 0.0);
    assert!(heating.last_step(CellId(1)).is_none());
    assert_eq!(initial.time_s(), 0.0);
    assert_eq!(initial.temperature_k(CellId(1)).unwrap(), 300.0);
}
