//! Manufactured event energies and calorics, not calibrated isotope data.

use std::sync::Arc;

use crate::EnergyChannels;
use crate::burnup::stochastic::{
    StochasticActivationBudget, StochasticActivationChange, StochasticActivationSchedule,
    StochasticActivationState,
};
use crate::burnup::{
    BurnupData, BurnupLimits, BurnupModel, BurnupNuclide, BurnupProduct, BurnupRate, BurnupReaction,
};
use crate::decay::stochastic::StreamKey;
use crate::irradiation::{GroupFlux, NeutronGroups, NuclideId};
use crate::thermal::{
    FissionThermalError, StochasticActivationCellSource, StochasticActivationHeating,
};
use fs_qty::Time;
use fs_thermal::enthalpy::{
    Boundary, CellGeometry, CellId, DepositedCell, DepositionBudget, DepositionError,
    DepositionSimulation, Event, MaterialId, StepConfig, ThermalLaw, TimedEvent,
};
use fs_thermochem::decay::DecayLawV1;

struct Caloric(f64);
impl ThermalLaw for Caloric {
    fn temperature_k(&self, energy: f64) -> Result<f64, DepositionError> {
        let temperature = energy / 4.0;
        if temperature.is_finite() && temperature > 0.0 && temperature <= self.0 {
            Ok(temperature)
        } else {
            Err(DepositionError::MaterialDomain)
        }
    }
    fn maximum_temperature_slope(&self) -> f64 {
        0.25
    }
}
static MATERIAL: Caloric = Caloric(1e6);
static LIMITED: Caloric = Caloric(300.25);

fn energy(deposited_j: f64, radiation_j: f64, escaped_j: f64) -> EnergyChannels {
    EnergyChannels {
        deposited_j,
        radiation_j,
        escaped_j,
    }
}
fn controls(max_dt_s: f64) -> StepConfig {
    StepConfig {
        max_dt_s,
        energy_tolerance_j: 1e-8,
    }
}
fn budget() -> StochasticActivationBudget {
    StochasticActivationBudget {
        maximum_events: 1024,
        maximum_clock_checks: 100_000,
        energy_tolerance_j: 1e-8,
    }
}
fn thermal(
    law: &'static dyn ThermalLaw,
    cells: u64,
    authored_w: f64,
) -> DepositionSimulation<'static> {
    let events = (1..=cells)
        .map(|id| TimedEvent {
            time_s: 0.0,
            event: Event::Birth(DepositedCell {
                id: CellId(id),
                segment_id: id,
                body_id: 1,
                material_id: MaterialId(1),
                mass_kg: 2.0,
                specific_h_j_kg: 1200.0,
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
                    source_w: authored_w,
                },
            }),
        })
        .collect();
    DepositionSimulation::new(
        events,
        &[(MaterialId(1), law)],
        DepositionBudget {
            max_cells: 4,
            max_links: 4,
            max_events: 8,
            max_substeps: 10_000,
        },
    )
    .unwrap()
    .advance_to(0.0, controls(1.0), || Ok(()))
    .unwrap()
    .0
}

fn source(nuclei: u64, shutdown: Option<f64>) -> StochasticActivationCellSource {
    let groups = Arc::new(
        NeutronGroups::new("manufactured activation grid".into(), vec![0.0, 1.0]).unwrap(),
    );
    let stable = DecayLawV1::stable();
    let model = Arc::new(
        BurnupModel::new(
            BurnupData {
                data_id: "manufactured integer activation calorimetry".into(),
                groups: Arc::clone(&groups),
                temperature_min_k: 200.0,
                temperature_max_k: 1000.0,
                nuclides: vec![
                    BurnupNuclide {
                        id: NuclideId::new(27, 59, 0).unwrap(),
                        law: stable,
                    },
                    BurnupNuclide {
                        id: NuclideId::new(27, 60, 0).unwrap(),
                        law: DecayLawV1::from_half_life(Time::new(1.0)).unwrap(),
                    },
                    BurnupNuclide {
                        id: NuclideId::new(28, 60, 0).unwrap(),
                        law: stable,
                    },
                ],
                reactions: vec![
                    BurnupReaction {
                        data_id: "manufactured target capture".into(),
                        parent: 0,
                        rate: BurnupRate::Capture {
                            cross_sections_m2: vec![1.0],
                        },
                        products: vec![BurnupProduct {
                            nuclide: 1,
                            yield_per_event: 1.0,
                        }],
                        energy_per_event: energy(4.0, 3.0, 0.5),
                    },
                    BurnupReaction {
                        data_id: "manufactured activated daughter decay".into(),
                        parent: 1,
                        rate: BurnupRate::Decay {
                            branch_fraction: 1.0,
                        },
                        products: vec![BurnupProduct {
                            nuclide: 2,
                            yield_per_event: 1.0,
                        }],
                        energy_per_event: energy(8.0, 2.0, 1.0),
                    },
                ],
            },
            BurnupLimits::default(),
        )
        .unwrap(),
    );
    let changes = shutdown
        .into_iter()
        .map(|time_s| StochasticActivationChange {
            time_s,
            flux: GroupFlux::new(Arc::clone(&groups), vec![0.0]).unwrap(),
            temperature_k: 300.0,
        })
        .collect();
    let schedule = Arc::new(
        StochasticActivationSchedule::new(
            model,
            GroupFlux::new(groups, vec![1.0]).unwrap(),
            300.0,
            changes,
            BurnupLimits::default(),
            1,
            || Ok(()),
        )
        .unwrap(),
    );
    StochasticActivationCellSource {
        cell_id: CellId(1),
        material_id: MaterialId(1),
        state: StochasticActivationState::new(
            schedule,
            vec![nuclei, 0, 0],
            StreamKey {
                seed: 7,
                kernel: 0x6163_7469,
                tile: 1,
            },
        )
        .unwrap(),
    }
}

#[test]
fn capture_heat_precedes_actual_daughter_decay_and_each_event_is_credited_once() {
    let initial = thermal(&MATERIAL, 1, 0.0);
    let source = source(1, None);
    let capture_time = source.state.next_transition_time_s().unwrap();
    let heating = StochasticActivationHeating::new(&initial, vec![source], budget()).unwrap();
    let (before, heating, audit) = heating
        .advance_to(&initial, capture_time / 2.0, controls(10.0), || Ok(()))
        .unwrap();
    assert_eq!(before.temperature_k(CellId(1)).unwrap(), 300.0);
    assert_eq!(audit.external_heat_j, 0.0);
    let (captured, heating, audit) = heating
        .advance_to(&before, capture_time, controls(10.0), || Ok(()))
        .unwrap();
    let state = &heating.source(CellId(1)).unwrap().state;
    assert_eq!(state.populations(), &[0, 1, 0]);
    assert_eq!(heating.audit().captures, 1.0);
    assert_eq!(heating.audit().decays, 0.0);
    assert_eq!(heating.audit().released, energy(4.0, 3.0, 0.5));
    assert_eq!(captured.temperature_k(CellId(1)).unwrap(), 300.5);
    assert_eq!(audit.external_heat_j, 4.0);
    let decay_time = state.next_transition_time_s().unwrap();
    assert!(decay_time > capture_time);
    let (decayed, heating, audit) = heating
        .advance_to(&captured, decay_time, controls(10.0), || Ok(()))
        .unwrap();
    assert_eq!(
        heating.source(CellId(1)).unwrap().state.populations(),
        &[0, 0, 1]
    );
    assert_eq!(heating.audit().captures, 0.0);
    assert_eq!(heating.audit().decays, 1.0);
    assert_eq!(audit.external_heat_j, 8.0);
    assert_eq!(decayed.temperature_k(CellId(1)).unwrap(), 301.5);
    let (late, heating, _) = heating
        .advance_to(&decayed, decay_time + 10.0, controls(10.0), || Ok(()))
        .unwrap();
    assert_eq!(late.temperature_k(CellId(1)).unwrap(), 301.5);
    assert_eq!(
        heating.cumulative_released().unwrap(),
        energy(12.0, 5.0, 1.5)
    );
    let state = &heating.source(CellId(1)).unwrap().state;
    assert_eq!(state.cumulative_capture_released(), energy(4.0, 3.0, 0.5));
    assert_eq!(state.cumulative_decay_released(), energy(8.0, 2.0, 1.0));
    assert_eq!(late.cells()[&CellId(1)].boundary.source_w, 0.0);
}

#[test]
fn shutdown_and_thermal_subdivision_preserve_integer_history_and_decay_tail() {
    let initial = thermal(&MATERIAL, 1, 2.0);
    let heating =
        StochasticActivationHeating::new(&initial, vec![source(12, Some(1.3))], budget()).unwrap();
    let (whole_thermal, whole, _) = heating
        .advance_to(&initial, 100.0, controls(100.0), || Ok(()))
        .unwrap();
    let (mut material, mut split) = (initial, heating);
    let mut captures_at_shutdown = 0;
    for endpoint in [0.1, 0.7, 1.3, 2.8, 100.0] {
        let (next_material, next, _) = split
            .advance_to(&material, endpoint, controls(0.17), || Ok(()))
            .unwrap();
        material = next_material;
        split = next;
        if endpoint == 1.3 {
            captures_at_shutdown = split
                .source(CellId(1))
                .unwrap()
                .state
                .cumulative_reaction_counts()[0];
        }
    }
    let accepted = &whole.source(CellId(1)).unwrap().state;
    assert_eq!(accepted, &split.source(CellId(1)).unwrap().state);
    assert_eq!(
        accepted.cumulative_reaction_counts(),
        &[captures_at_shutdown, captures_at_shutdown]
    );
    assert!(captures_at_shutdown > 0 && captures_at_shutdown < 12);
    assert_eq!(
        accepted.populations(),
        &[12 - captures_at_shutdown, 0, captures_at_shutdown]
    );
    let expected_k = 300.0 + (200.0 + 12.0 * captures_at_shutdown as f64) / 8.0;
    assert!((whole_thermal.temperature_k(CellId(1)).unwrap() - expected_k).abs() < 1e-9);
    assert!((material.temperature_k(CellId(1)).unwrap() - expected_k).abs() < 1e-9);
    assert_eq!(material.cells()[&CellId(1)].boundary.source_w, 2.0);
}

#[test]
fn coupled_refusals_preserve_material_rng_and_pending_capture() {
    let initial = thermal(&LIMITED, 2, 0.0);
    let source = source(1, None);
    let original = source.state.clone();
    let capture_time = source.state.next_transition_time_s().unwrap();
    let heating =
        StochasticActivationHeating::new(&initial, vec![source.clone()], budget()).unwrap();
    assert!(
        heating
            .advance_to(&initial, capture_time, controls(10.0), || Ok(()))
            .is_err()
    );
    assert_eq!(heating.source(CellId(1)).unwrap().state, original);
    assert_eq!(initial.time_s(), 0.0);
    assert_eq!(initial.temperature_k(CellId(1)).unwrap(), 300.0);
    assert!(matches!(
        heating.advance_to(&initial, capture_time, controls(10.0), || Err(
            DepositionError::Cancelled
        )),
        Err(FissionThermalError::Thermal(DepositionError::Cancelled))
    ));
    assert_eq!(heating.source(CellId(1)).unwrap().state, original);
    let mut duplicate = source.clone();
    duplicate.cell_id = CellId(2);
    assert!(StochasticActivationHeating::new(&initial, vec![source, duplicate], budget()).is_err());
    let mut advanced = original.clone();
    advanced
        .advance_to(capture_time / 2.0, budget(), || Ok(()))
        .unwrap();
    assert!(matches!(
        StochasticActivationHeating::new(
            &initial,
            vec![StochasticActivationCellSource {
                cell_id: CellId(1),
                material_id: MaterialId(1),
                state: advanced,
            }],
            budget()
        ),
        Err(FissionThermalError::ClockMismatch { .. })
    ));
}

#[test]
fn grouped_radiation_uses_actual_decay_counts_and_excludes_prompt_capture() {
    use crate::endf::decay::model::emissions::{DecayEmissionLimits, PreparedDecayEmission};
    use crate::endf::decay::model::{DecayModelEntry, EvaluatedDecayAssembly};
    use crate::endf::decay::spectra::DecayRadiation;
    use crate::endf::decay::{DecayEvaluation, EndfDecayReadLimits, EndfDecaySelection};
    use crate::radiation::RadiationDeposit;
    use crate::radiation::decay::{
        DecayRadiationClosure, PreparedDecayRadiation, tests::gamma_chain,
    };
    use crate::radiation::grouped::{
        GroupedRadiationResponse, GroupedRadiationResponseLimits, GroupedRadiationSourceResponse,
        RadiationGroupResponse,
    };
    use crate::thermal::{DecayRadiationEmitter, RadiationReceiver};

    // Deliberately large manufactured eV values make integer-event calorimetry
    // representable in this macroscopic thermal fixture. The parsed values,
    // never the pre-rounded text inputs, own the expected energy account.
    let ev_j = 1.602_176_634e-19;
    let chain = gamma_chain(1e20, 2e20, 0.5, 0.25);
    let rows = [
        (["1.002900+4", "2.900000+1", "0", "0", "1", "0"], 457),
        (["0", "0", "0", "0", "6", "0"], 457),
        (["0"; 6], 457),
        (["0", "1", "0", "0", "6", "0"], 457),
        (["0"; 6], 457),
        (["0"; 6], 0),
    ];
    let precursor_text = rows
        .into_iter()
        .map(|(fields, mt)| {
            let values = fields.map(|value| format!("{value:>11}")).concat();
            format!("{values}{:4}{:2}{mt:3}{:5}", 778, 8, 1)
        })
        .collect::<Vec<_>>()
        .join("\n");
    let precursor = Arc::new(
        DecayEvaluation::parse(
            "manufactured explicitly stable activation precursor".into(),
            EndfDecaySelection {
                material_number: 778,
                nuclide: NuclideId::new(10, 29, 0).unwrap(),
            },
            &precursor_text,
            EndfDecayReadLimits::default(),
            || Ok(()),
        )
        .unwrap(),
    );
    let mut entries: Vec<_> = chain
        .evaluations()
        .iter()
        .zip(chain.routing())
        .map(|(evaluation, &routing)| DecayModelEntry {
            evaluation: Arc::clone(evaluation),
            routing,
        })
        .collect();
    entries.push(DecayModelEntry {
        evaluation: precursor,
        routing: chain.routing()[2],
    });
    let assembly = Arc::new(
        EvaluatedDecayAssembly::new(
            "manufactured evaluated activation network".into(),
            entries,
            fs_la::decay::DecayLimits::default(),
            || Ok(()),
        )
        .unwrap(),
    );
    let limits = DecayEmissionLimits::default();
    let partition = Arc::new(
        PreparedDecayRadiation::new(
            Arc::new(
                PreparedDecayEmission::new(
                    Arc::clone(&assembly),
                    DecayRadiation::Gamma,
                    vec![0.0, 1.5e20 * ev_j, 3e20 * ev_j],
                    limits,
                    || Ok(()),
                )
                .unwrap(),
            ),
            DecayRadiationClosure {
                absolute_j_per_decay: 0.0,
                relative: 1e-12,
            },
            limits,
            || Ok(()),
        )
        .unwrap(),
    );
    let groups = Arc::new(
        NeutronGroups::new(
            "manufactured spectral activation grid".into(),
            vec![0.0, 1.0],
        )
        .unwrap(),
    );
    let mut data = assembly
        .burnup_data(Arc::clone(&groups), 200.0, 1000.0, || Ok(()))
        .unwrap();
    // All four identities and their half-lives remain bound to the evaluations.
    data.reactions.insert(
        0,
        BurnupReaction {
            data_id: "manufactured prompt capture radiation".into(),
            parent: 3,
            rate: BurnupRate::Capture {
                cross_sections_m2: vec![1.0],
            },
            products: vec![BurnupProduct {
                nuclide: 0,
                yield_per_event: 1.0,
            }],
            energy_per_event: energy(4.0, 3.0, 0.5),
        },
    );
    let model = Arc::new(BurnupModel::new(data, BurnupLimits::default()).unwrap());
    let schedule = Arc::new(
        StochasticActivationSchedule::new(
            Arc::clone(&model),
            GroupFlux::new(groups, vec![1.0]).unwrap(),
            300.0,
            vec![],
            BurnupLimits::default(),
            0,
            || Ok(()),
        )
        .unwrap(),
    );
    let state = StochasticActivationState::new(
        schedule,
        vec![0, 0, 0, 1],
        StreamKey {
            seed: 19,
            kernel: 0x6163_7469,
            tile: 1,
        },
    )
    .unwrap();
    let original = state.clone();
    let (expected, report) = state.trial_to(100.0, budget(), || Ok(())).unwrap();
    assert_eq!(report.reaction_counts, vec![1, 1, 1]);
    let response = Arc::new(
        GroupedRadiationResponse::new(
            "manufactured energy-dependent capture-product response".into(),
            DecayRadiation::Gamma,
            partition.group_edges_j().to_vec(),
            vec![GroupedRadiationSourceResponse {
                source_id: 1,
                groups: [0.75, 0.25]
                    .map(|fraction| RadiationGroupResponse {
                        deposits: vec![RadiationDeposit {
                            receiver_id: 2,
                            fraction,
                        }],
                        escape_fraction: 1.0 - fraction,
                    })
                    .to_vec(),
            }],
            GroupedRadiationResponseLimits {
                max_sources: 1,
                max_groups: 2,
                max_receivers: 1,
                max_source_groups: 2,
                max_transfers: 2,
            },
        )
        .unwrap(),
    );
    let binding = partition
        .bind_burnup(Arc::clone(&model), || Ok(()))
        .unwrap();
    let projected = binding
        .from_stochastic_activation_step(&model, &report, || Ok(()))
        .unwrap();
    assert!((projected.original_radiation_j - report.decay_released.radiation_j).abs() < 1e-12);
    let mut tampered = report.clone();
    tampered.captures += 1;
    assert!(
        binding
            .from_stochastic_activation_step(&model, &tampered, || Ok(()))
            .is_err()
    );
    let initial = thermal(&MATERIAL, 2, 0.0);
    let heating = StochasticActivationHeating::new(
        &initial,
        vec![StochasticActivationCellSource {
            cell_id: CellId(1),
            material_id: MaterialId(1),
            state,
        }],
        budget(),
    )
    .unwrap()
    .with_grouped_decay_radiation(
        &initial,
        response,
        vec![DecayRadiationEmitter {
            cell_id: CellId(1),
            partition,
        }],
        vec![RadiationReceiver {
            cell_id: CellId(2),
            material_id: MaterialId(1),
        }],
        || Ok(()),
    )
    .unwrap();
    let (heated, accepted, thermal_audit) = heating
        .advance_to(&initial, 100.0, controls(100.0), || Ok(()))
        .unwrap();
    assert_eq!(accepted.source(CellId(1)).unwrap().state, expected);
    assert_eq!(heating.source(CellId(1)).unwrap().state, original);
    let received = 0.75 * model.reactions()[1].energy_per_event.radiation_j
        + 0.25 * model.reactions()[2].energy_per_event.radiation_j;
    let audit = accepted.audit();
    assert_eq!(audit.captures, 1.0);
    assert_eq!(audit.decays, 2.0);
    assert_eq!(audit.radiation_ungrouped_j, 3.0);
    assert!((audit.radiation_deposited_j - received).abs() < 1e-12);
    assert!((audit.radiation_escaped_j - (report.released.radiation_j - received)).abs() < 1e-12);
    assert!(
        (thermal_audit.external_heat_j - (report.released.deposited_j + received)).abs() < 1e-10
    );
    assert!(
        (heated.temperature_k(CellId(1)).unwrap() - (300.0 + report.released.deposited_j / 8.0))
            .abs()
            < 1e-10
    );
    assert!((heated.temperature_k(CellId(2)).unwrap() - (300.0 + received / 8.0)).abs() < 1e-10);
}
