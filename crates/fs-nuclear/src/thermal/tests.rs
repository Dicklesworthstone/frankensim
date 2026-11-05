use super::*;
use crate::{DecayHeatGroup, FissionBudget, FissionModel};
use fs_thermal::enthalpy::{
    Boundary, DepositedCell, DepositionBudget, ThermalLaw, ThermalLink, TimedEvent,
};

use crate::irradiation::{
    GroupCrossSections, GroupFlux, IrradiationBudget, IrradiationData, IrradiationModel,
    IrradiationState, NeutronGroups, NuclideId,
};

// Manufactured caloric/source data isolate the coupling equations. These
// coefficients are not nuclear isotope data or an empirical validation deck.
struct LinearLaw {
    cp: f64,
    identity: u8,
    maximum_k: f64,
}
impl ThermalLaw for LinearLaw {
    fn caloric_identity(&self) -> Option<[u8; 32]> {
        Some([self.identity; 32])
    }
    fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
        let temperature = h / self.cp;
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
};
static COOL: LinearLaw = LinearLaw {
    cp: 2.0,
    identity: 2,
    maximum_k: 1e6,
};
static LIMITED: LinearLaw = LinearLaw {
    cp: 4.0,
    identity: 3,
    maximum_k: 301.0,
};

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:?} != {expected:?}, tol {tolerance}"
    );
}
fn controls(max_dt_s: f64) -> StepConfig {
    StepConfig {
        max_dt_s,
        energy_tolerance_j: 1e-8,
    }
}
fn budget() -> FissionBudget {
    FissionBudget {
        max_decay_groups: 8,
        energy_tolerance_j: 1e-8,
        fuel_tolerance_nuclei: 1e-9,
    }
}
fn cell(id: u64, material: u64, mass: f64, h: f64, source: f64) -> DepositedCell {
    DepositedCell {
        id: CellId(id),
        segment_id: id + 10,
        body_id: id + 100,
        material_id: MaterialId(material),
        mass_kg: mass,
        specific_h_j_kg: h,
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
fn model(delayed: bool) -> Arc<FissionModel> {
    Arc::new(
        FissionModel::new(
            "manufactured-fission-calorimetry-v1".into(),
            EnergyChannels {
                deposited_j: 8.0,
                radiation_j: 1.0,
                escaped_j: 1.0,
            },
            if delayed {
                vec![DecayHeatGroup {
                    decay_constant_s_inv: 0.5,
                    energy_per_fission: EnergyChannels {
                        deposited_j: 2.0,
                        radiation_j: 1.0,
                        escaped_j: 1.0,
                    },
                }]
            } else {
                vec![]
            },
        )
        .unwrap(),
    )
}
fn source(
    id: u64,
    material: u64,
    initial_rate: f64,
    changes: &[(f64, f64)],
    delayed: bool,
) -> FissionCellSource {
    FissionCellSource {
        cell_id: CellId(id),
        material_id: MaterialId(material),
        state: FissionState::new(
            model(delayed),
            1000.0,
            if delayed { vec![0.0] } else { vec![] },
        )
        .unwrap(),
        rates: FissionRateSchedule::new(
            initial_rate,
            changes
                .iter()
                .map(|&(time_s, fissions_per_s)| FissionRateChange {
                    time_s,
                    fissions_per_s,
                })
                .collect(),
        )
        .unwrap(),
    }
}

#[test]
fn g1_shutdown_integrates_delayed_heat_and_adds_authored_heat_exactly_once() {
    let initial = thermal(
        vec![cell(1, 1, 2.0, 1200.0, 3.0)],
        &[(MaterialId(1), &HOT)],
        vec![],
    );
    let heating = FissionHeating::new(
        &initial,
        vec![source(1, 1, 5.0, &[(2.0, 0.0)], true)],
        budget(),
    )
    .unwrap();
    // max_dt exceeds the shutdown time: the source's transition hook MUST cut it.
    let (final_thermal, final_source, audit) = heating
        .advance_to(&initial, 5.0, controls(10.0), || Ok(()))
        .unwrap();
    let stored = 40.0 * (1.0 - (-1.0_f64).exp()) * (-1.5_f64).exp();
    let deposited = 80.0 + 0.5 * (40.0 - stored);
    close(
        final_source
            .source(CellId(1))
            .unwrap()
            .state
            .delayed_energy_j()[0],
        stored,
        1e-12,
    );
    close(final_source.audit().released.deposited_j, deposited, 1e-12);
    close(
        final_source.audit().released.radiation_j,
        10.0 + 0.25 * (40.0 - stored),
        1e-12,
    );
    close(
        final_source.audit().released.escaped_j,
        10.0 + 0.25 * (40.0 - stored),
        1e-12,
    );
    close(audit.external_heat_j, deposited + 15.0, 1e-12);
    close(
        final_thermal.temperature_k(CellId(1)).unwrap(),
        300.0 + (deposited + 15.0) / 8.0,
        1e-12,
    );
    assert_eq!(
        final_source
            .source(CellId(1))
            .unwrap()
            .state
            .cumulative_fissions(),
        10.0
    );
    assert_eq!(final_source.time_s(), 5.0);
    assert_eq!(audit.substeps, 2);
    assert_eq!(final_source.audit().substeps, 2);
    assert_eq!(final_thermal.cells()[&CellId(1)].boundary.source_w, 3.0);
    assert_eq!(
        heating
            .source(CellId(1))
            .unwrap()
            .state
            .remaining_fuel_nuclei(),
        1000.0
    );
    assert_eq!(initial.time_s(), 0.0);

    let (heated_later, source_later, later_audit) = final_source
        .advance_to(&final_thermal, 9.0, controls(10.0), || Ok(()))
        .unwrap();
    assert!(
        source_later.audit().released.deposited_j > 0.0,
        "shutdown does not erase delayed heat"
    );
    assert_eq!(source_later.audit().fissions, 0.0);
    close(
        later_audit.external_heat_j,
        source_later.audit().thermal_credit_j + 12.0,
        1e-12,
    );
    assert!(
        heated_later.temperature_k(CellId(1)).unwrap()
            > final_thermal.temperature_k(CellId(1)).unwrap()
    );
}

#[test]
fn g1_finite_coolant_uses_original_contact_heat_and_converges_to_two_body_solution() {
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
    let heating =
        FissionHeating::new(&initial, vec![source(1, 1, 1.0, &[], false)], budget()).unwrap();
    let end = 4.0;
    let average = (8.0 * 300.0 + 6.0 * 280.0 + 8.0 * end) / 14.0;
    let difference = 24.0 / 7.0 + (20.0 - 24.0 / 7.0) * (-7.0_f64 * end / 24.0).exp();
    let exact_hot = average + 6.0 / 14.0 * difference;
    let exact_cool = average - 8.0 / 14.0 * difference;
    let mut errors = Vec::new();
    for max_dt in [0.125, 0.0625, 0.03125] {
        let (result, nuclear, audit) = heating
            .advance_to(&initial, end, controls(max_dt), || Ok(()))
            .unwrap();
        close(audit.external_heat_j, 32.0, 1e-11);
        close(result.energy_j() - initial.energy_j(), 32.0, 1e-9);
        close(nuclear.audit().thermal_credit_j, 32.0, 1e-12);
        // Radiation/escape are exports; depositing them too would add another 8 J.
        close(
            nuclear.audit().released.radiation_j + nuclear.audit().released.escaped_j,
            8.0,
            1e-12,
        );
        let hot = result.temperature_k(CellId(1)).unwrap();
        let cool = result.temperature_k(CellId(2)).unwrap();
        assert!(cool > 280.0);
        errors.push((hot - exact_hot).abs().max((cool - exact_cool).abs()));
        assert_eq!(result.links(), initial.links());
    }
    assert!(errors[2] < 0.02);
    for pair in errors.windows(2) {
        assert!((1.9..2.1).contains(&(pair[0] / pair[1])));
    }
}

#[test]
fn g3_multiple_material_sources_preserve_rate_transitions_and_split_restart() {
    let initial = thermal(
        vec![cell(1, 1, 2.0, 1200.0, 0.0), cell(2, 2, 3.0, 600.0, 0.0)],
        &[(MaterialId(1), &HOT), (MaterialId(2), &COOL)],
        vec![],
    );
    let heating = FissionHeating::new(
        &initial,
        vec![
            source(2, 2, 0.0, &[(1.0, 2.0), (3.0, 0.0)], true),
            source(1, 1, 5.0, &[(2.0, 0.0)], true),
        ],
        budget(),
    )
    .unwrap();
    assert_eq!(
        heating.sources().map(|s| s.cell_id).collect::<Vec<_>>(),
        [CellId(1), CellId(2)]
    );
    let (whole_t, whole_s, whole_audit) = heating
        .advance_to(&initial, 5.0, controls(20.0), || Ok(()))
        .unwrap();
    assert_eq!(whole_audit.substeps, 4);
    assert_eq!(
        whole_s
            .source(CellId(1))
            .unwrap()
            .state
            .cumulative_fissions(),
        10.0
    );
    assert_eq!(
        whole_s
            .source(CellId(2))
            .unwrap()
            .state
            .cumulative_fissions(),
        4.0
    );
    let (part_t, part_s, _) = heating
        .advance_to(&initial, 2.0, controls(20.0), || Ok(()))
        .unwrap();
    let checkpoint_s = part_s.clone();
    let checkpoint_t = part_t.clone();
    let (split_t, split_s, _) = checkpoint_s
        .advance_to(&checkpoint_t, 5.0, controls(20.0), || Ok(()))
        .unwrap();
    assert_eq!(
        whole_t.cells(),
        split_t.cells(),
        "same physical intervals replay bitwise"
    );
    for cell in [CellId(1), CellId(2)] {
        assert_eq!(
            whole_s.source(cell).unwrap().state,
            split_s.source(cell).unwrap().state
        );
    }
    close(
        whole_audit.external_heat_j,
        whole_s.cumulative_released().unwrap().deposited_j,
        1e-11,
    );
    let (repeat_t, repeat_s, repeat_audit) = heating
        .advance_to(&initial, 5.0, controls(20.0), || Ok(()))
        .unwrap();
    assert_eq!(repeat_t.cells(), whole_t.cells());
    assert_eq!(repeat_s.audit(), whole_s.audit());
    assert_eq!(repeat_audit, whole_audit);

    let (idle_t, idle_s, idle_audit) = whole_s
        .advance_to(&whole_t, 5.0, controls(1.0), || Ok(()))
        .unwrap();
    assert_eq!(idle_t.cells(), whole_t.cells());
    assert_eq!(idle_s.audit(), FissionHeatingAudit::default());
    assert_eq!(idle_audit.substeps, 0);
    assert!(matches!(
        heating.advance_to(&whole_t, 5.0, controls(1.0), || Ok(())),
        Err(FissionThermalError::ClockMismatch { .. })
    ));
}

#[test]
fn g0_fuel_exhaustion_inside_large_thermal_interval_still_releases_decay_heat() {
    let initial = thermal(
        vec![cell(1, 1, 2.0, 1200.0, 0.0)],
        &[(MaterialId(1), &HOT)],
        vec![],
    );
    let mut fuel = source(1, 1, 5.0, &[], true);
    fuel.state = FissionState::new(model(true), 10.0, vec![0.0]).unwrap();
    let heating = FissionHeating::new(&initial, vec![fuel], budget()).unwrap();
    let (result, nuclear, audit) = heating
        .advance_to(&initial, 5.0, controls(10.0), || Ok(()))
        .unwrap();
    assert_eq!(audit.substeps, 1);
    let step = nuclear.last_step(CellId(1)).unwrap();
    assert_eq!(step.fission_duration_s, 2.0);
    assert_eq!(
        nuclear
            .source(CellId(1))
            .unwrap()
            .state
            .remaining_fuel_nuclei(),
        0.0
    );
    let residual = 40.0 * (1.0 - (-1.0_f64).exp()) * (-1.5_f64).exp();
    close(
        result.energy_j() - initial.energy_j(),
        100.0 - 0.5 * residual,
        1e-10,
    );
}

#[test]
fn g4_thermal_refusal_and_late_cancellation_publish_neither_history() {
    let initial = thermal(
        vec![cell(1, 1, 2.0, 1200.0, 0.0)],
        &[(MaterialId(1), &LIMITED)],
        vec![],
    );
    let heating = FissionHeating::new(
        &initial,
        vec![source(1, 1, 5.0, &[(2.0, 0.0)], true)],
        budget(),
    )
    .unwrap();
    let source_before = heating.source(CellId(1)).unwrap().state.clone();
    assert!(matches!(
        heating.advance_to(&initial, 5.0, controls(10.0), || Ok(())),
        Err(FissionThermalError::Thermal(
            DepositionError::MaterialDomain
        ))
    ));
    assert_eq!(heating.source(CellId(1)).unwrap().state, source_before);
    assert_eq!(initial.time_s(), 0.0);

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
        assert_eq!(heating.source(CellId(1)).unwrap().state, source_before);
    }
    let retry = heating
        .advance_to(&initial, 0.1, controls(0.025), || Ok(()))
        .unwrap();
    assert_eq!(retry.0.cells(), expected.0.cells());
    assert_eq!(
        retry.1.source(CellId(1)).unwrap().state,
        expected.1.source(CellId(1)).unwrap().state
    );
    assert_eq!(retry.2, expected.2);

    // The callback inside the nuclear trial preserves the original thermal error,
    // rather than erasing a caller budget refusal into generic cancellation.
    let mut candidate = heating.clone();
    let mut calls = 0;
    let result = candidate.step_sources(&initial, 0.0, 0.01, 0.01, &mut || {
        calls += 1;
        if calls == 4 {
            Err(DepositionError::Budget("caller checkpoint"))
        } else {
            Ok(())
        }
    });
    assert!(matches!(
        result,
        Err(FissionThermalError::Thermal(DepositionError::Budget(
            "caller checkpoint"
        )))
    ));
    assert_eq!(candidate.source(CellId(1)).unwrap().state, source_before);
    assert!(candidate.pending.is_none());
}

#[test]
fn g0_bad_schedules_and_stale_material_bindings_are_refused() {
    for initial in [-1.0, f64::NAN, f64::INFINITY] {
        assert!(FissionRateSchedule::new(initial, vec![]).is_err());
    }
    for changes in [
        vec![(0.0, 1.0)],
        vec![(2.0, 1.0), (2.0, 2.0)],
        vec![(f64::INFINITY, 1.0)],
        vec![(1.0, -1.0)],
        vec![(1.0, f64::NAN)],
        vec![(2.0, 1.0), (1.0, 2.0)],
    ] {
        assert!(
            FissionRateSchedule::new(
                0.0,
                changes
                    .into_iter()
                    .map(|(time_s, fissions_per_s)| FissionRateChange {
                        time_s,
                        fissions_per_s
                    })
                    .collect()
            )
            .is_err()
        );
    }
    let rates = FissionRateSchedule::new(
        3.0,
        vec![FissionRateChange {
            time_s: 2.0,
            fissions_per_s: 0.0,
        }],
    )
    .unwrap();
    assert_eq!(rates.rate_at(0.0).unwrap(), 3.0);
    assert_eq!(rates.rate_at(2.0).unwrap(), 0.0);
    assert!(rates.rate_at(f64::NAN).is_err());

    let initial = thermal(
        vec![cell(1, 1, 2.0, 1200.0, 0.0), cell(2, 2, 3.0, 600.0, 0.0)],
        &[(MaterialId(1), &HOT), (MaterialId(2), &COOL)],
        vec![],
    );
    assert!(matches!(
        FissionHeating::new(&initial, vec![source(1, 2, 0.0, &[], true)], budget()),
        Err(FissionThermalError::MaterialMismatch { .. })
    ));
    assert!(matches!(
        FissionHeating::new(&initial, vec![source(99, 1, 0.0, &[], true)], budget()),
        Err(FissionThermalError::InactiveCell(CellId(99)))
    ));
    assert!(
        FissionHeating::new(
            &initial,
            vec![source(1, 1, 0.0, &[], true), source(1, 1, 0.0, &[], true)],
            budget()
        )
        .is_err()
    );
    let heating =
        FissionHeating::new(&initial, vec![source(1, 1, 0.0, &[], true)], budget()).unwrap();
    let different_calorics = thermal(
        vec![cell(1, 1, 2.0, 1200.0, 0.0)],
        &[(MaterialId(1), &LIMITED)],
        vec![],
    );
    assert!(matches!(
        heating.advance_to(&different_calorics, 0.0, controls(1.0), || Ok(())),
        Err(FissionThermalError::BindingChanged(CellId(1)))
    ));
    let different_mass = thermal(
        vec![cell(1, 1, 3.0, 1200.0, 0.0)],
        &[(MaterialId(1), &HOT)],
        vec![],
    );
    assert!(matches!(
        heating.advance_to(&different_mass, 0.0, controls(1.0), || Ok(())),
        Err(FissionThermalError::BindingChanged(CellId(1)))
    ));
}

fn irradiation_budget() -> IrradiationHeatingBudget {
    IrradiationHeatingBudget {
        max_isotopes_per_cell: 4,
        source: IrradiationBudget {
            max_neutron_groups: 2,
            source: budget(),
        },
    }
}

fn irradiated_source(maximum_k: f64, delayed: bool) -> IrradiationCellSource {
    let groups = Arc::new(
        NeutronGroups::new(
            "manufactured-two-group-grid".into(),
            vec![0.0, 1e-15, 1e-12],
        )
        .unwrap(),
    );
    let flux = GroupFlux::new(groups.clone(), vec![2e27, 1e27]).unwrap();
    let zero = GroupFlux::new(groups.clone(), vec![0.0, 0.0]).unwrap();
    let mut inventories = Vec::new();
    // Manufactured coefficients, not measured U-235/U-238 nuclear data.
    for (a, atoms, cross_sections) in [
        (
            235,
            100.0,
            vec![
                GroupCrossSections {
                    fission_m2: 5e-29,
                    capture_m2: 0.0,
                },
                GroupCrossSections {
                    fission_m2: 5e-29,
                    capture_m2: 5e-29,
                },
            ],
        ),
        (
            238,
            80.0,
            vec![
                GroupCrossSections {
                    fission_m2: 0.0,
                    capture_m2: 5e-29,
                },
                GroupCrossSections {
                    fission_m2: 5e-29,
                    capture_m2: 5e-29,
                },
            ],
        ),
    ] {
        let data = Arc::new(
            IrradiationModel::new(IrradiationData {
                data_id: format!("manufactured-parent-{a}"),
                nuclide: NuclideId::new(92, a, 0).unwrap(),
                groups: groups.clone(),
                cross_sections,
                temperature_min_k: 290.0,
                temperature_max_k: maximum_k,
                fission_model: model(delayed),
                capture_energy_per_event: EnergyChannels {
                    deposited_j: 3.0,
                    radiation_j: 2.0,
                    escaped_j: 0.0,
                },
            })
            .unwrap(),
        );
        inventories.push(
            IrradiationState::new(data, atoms, if delayed { vec![0.0] } else { vec![] }).unwrap(),
        );
    }
    IrradiationCellSource {
        cell_id: CellId(1),
        material_id: MaterialId(1),
        inventories,
        fluxes: NeutronFluxSchedule::new(
            flux,
            vec![NeutronFluxChange {
                time_s: 2.0,
                flux: zero,
            }],
        )
        .unwrap(),
    }
}

#[test]
fn g1_multigroup_isotope_mixture_heats_by_actual_fissions_and_captures() {
    let initial = thermal(
        vec![cell(1, 1, 2.0, 1200.0, 3.0)],
        &[(MaterialId(1), &HOT)],
        vec![],
    );
    let heating = IrradiationHeating::new(
        &initial,
        vec![irradiated_source(1000.0, false)],
        irradiation_budget(),
    )
    .unwrap();
    let (result, heated, thermal_audit) = heating
        .advance_to(&initial, 5.0, controls(10.0), || Ok(()))
        .unwrap();
    // Both parents have kf+kc=0.2/s; their competing fractions differ.
    let loss = -(-0.4_f64).exp_m1();
    let fissions = 95.0 * loss;
    let captures = 85.0 * loss;
    let deposited = 8.0 * fissions + 3.0 * captures;
    close(heated.audit().fissions, fissions, 2e-13);
    close(heated.audit().captures, captures, 2e-13);
    close(heated.audit().released.deposited_j, deposited, 2e-12);
    close(
        heated.audit().released.radiation_j,
        fissions + 2.0 * captures,
        2e-12,
    );
    close(heated.audit().released.escaped_j, fissions, 2e-12);
    close(
        result.temperature_k(CellId(1)).unwrap(),
        300.0 + (deposited + 15.0) / 8.0,
        1e-12,
    );
    close(thermal_audit.external_heat_j, deposited + 15.0, 2e-12);
    assert_eq!(
        heated.audit().substeps,
        2,
        "clip the original step exactly at flux shutdown"
    );
    let source = heated.source(CellId(1)).unwrap();
    for (inventory, initial_atoms) in source.inventories.iter().zip([100.0, 80.0]) {
        close(
            inventory.remaining_parent_nuclei(),
            initial_atoms * (-0.4_f64).exp(),
            1e-13,
        );
    }
    assert_eq!(
        result.cells()[&CellId(1)].mass_kg,
        2.0,
        "consumed parents remain material products"
    );
    assert_eq!(result.cells()[&CellId(1)].boundary.source_w, 3.0);
    let (_, later, _) = heated
        .advance_to(&result, 8.0, controls(10.0), || Ok(()))
        .unwrap();
    assert_eq!(later.audit().fissions, 0.0);
    assert_eq!(later.audit().captures, 0.0);
    assert_eq!(later.audit().released, EnergyChannels::default());
}

#[test]
fn g3_irradiated_mixture_keeps_decay_heat_and_restart_at_flux_shutdown() {
    let initial = thermal(
        vec![cell(1, 1, 2.0, 1200.0, 0.0)],
        &[(MaterialId(1), &HOT)],
        vec![],
    );
    let heating = IrradiationHeating::new(
        &initial,
        vec![irradiated_source(1000.0, true)],
        irradiation_budget(),
    )
    .unwrap();
    let (whole, whole_source, _) = heating
        .advance_to(&initial, 5.0, controls(10.0), || Ok(()))
        .unwrap();
    let (at_shutdown, shutdown_source, _) = heating
        .advance_to(&initial, 2.0, controls(10.0), || Ok(()))
        .unwrap();
    let (split, split_source, after) = shutdown_source
        .advance_to(&at_shutdown, 5.0, controls(10.0), || Ok(()))
        .unwrap();
    assert_eq!(whole.energy_j(), split.energy_j());
    assert_eq!(
        whole_source.source(CellId(1)).unwrap().inventories,
        split_source.source(CellId(1)).unwrap().inventories
    );
    assert_eq!(split_source.audit().fissions, 0.0);
    assert_eq!(split_source.audit().captures, 0.0);
    assert!(after.external_heat_j > 0.0);
    close(
        whole.energy_j() - initial.energy_j(),
        whole_source.cumulative_released().unwrap().deposited_j,
        1e-9,
    );
    let report = whole_source.last_step(CellId(1)).unwrap();
    assert_eq!(report.isotopes.len(), 2);
}

#[test]
fn g4_irradiation_rejects_a_thermal_endpoint_outside_nuclear_data() {
    // The thermal caloric law admits this temperature; the nuclear data do not.
    let initial = thermal(
        vec![cell(1, 1, 2.0, 1200.0, 0.0)],
        &[(MaterialId(1), &HOT)],
        vec![],
    );
    let heating = IrradiationHeating::new(
        &initial,
        vec![irradiated_source(300.1, false)],
        irradiation_budget(),
    )
    .unwrap();
    assert!(matches!(
        heating.advance_to(&initial, 0.5, controls(1.0), || Ok(())),
        Err(FissionThermalError::Fission(FissionError::InvalidInput(_)))
    ));
    assert_eq!(initial.time_s(), 0.0);
    assert_eq!(initial.temperature_k(CellId(1)).unwrap(), 300.0);
    assert_eq!(heating.time_s(), 0.0);
    for inventory in &heating.source(CellId(1)).unwrap().inventories {
        assert_eq!(inventory.cumulative_fissions(), 0.0);
        assert_eq!(inventory.cumulative_captures(), 0.0);
    }
}

#[test]
fn g0_irradiation_bindings_refuse_duplicate_isotopes_and_changed_group_identity() {
    let initial = thermal(
        vec![cell(1, 1, 2.0, 1200.0, 0.0)],
        &[(MaterialId(1), &HOT)],
        vec![],
    );
    let mut source = irradiated_source(1000.0, false);
    source.inventories[1] = source.inventories[0].clone();
    assert!(IrradiationHeating::new(&initial, vec![source], irradiation_budget()).is_err());
    let mut too_small = irradiation_budget();
    too_small.max_isotopes_per_cell = 1;
    assert!(
        IrradiationHeating::new(&initial, vec![irradiated_source(1000.0, false)], too_small)
            .is_err()
    );
    let source = irradiated_source(1000.0, false);
    let initial_flux = source.fluxes.flux_at(0.0).unwrap().clone();
    let changed = Arc::new(
        NeutronGroups::new("different-grid-identity".into(), vec![0.0, 1e-15, 1e-12]).unwrap(),
    );
    let wrong = GroupFlux::new(changed, vec![0.0, 0.0]).unwrap();
    assert!(
        NeutronFluxSchedule::new(
            initial_flux,
            vec![NeutronFluxChange {
                time_s: 2.0,
                flux: wrong
            },]
        )
        .is_err()
    );
}

#[test]
fn g0_mixed_isotope_heat_reduction_cannot_drop_a_joule_outside_the_budget() {
    let initial = thermal(
        vec![cell(1, 1, 1e16, 1200.0, 0.0)],
        &[(MaterialId(1), &HOT)],
        vec![],
    );
    let groups =
        Arc::new(NeutronGroups::new("manufactured-roundoff-grid".into(), vec![0.0, 1.0]).unwrap());
    let flux = GroupFlux::new(groups.clone(), vec![core::f64::consts::LN_2]).unwrap();
    let mut inventories = Vec::new();
    for (a, deposited_j) in [(235, 1e16), (238, 1.0)] {
        let fission_model = Arc::new(
            FissionModel::new(
                format!("manufactured-roundoff-yield-{a}"),
                EnergyChannels {
                    deposited_j,
                    radiation_j: 0.0,
                    escaped_j: 0.0,
                },
                vec![],
            )
            .unwrap(),
        );
        let data = Arc::new(
            IrradiationModel::new(IrradiationData {
                data_id: format!("manufactured-roundoff-parent-{a}"),
                nuclide: NuclideId::new(92, a, 0).unwrap(),
                groups: groups.clone(),
                cross_sections: vec![GroupCrossSections {
                    fission_m2: 1.0,
                    capture_m2: 0.0,
                }],
                temperature_min_k: 290.0,
                temperature_max_k: 1000.0,
                fission_model,
                capture_energy_per_event: EnergyChannels::default(),
            })
            .unwrap(),
        );
        inventories.push(IrradiationState::new(data, 2.0, vec![]).unwrap());
    }
    let source = IrradiationCellSource {
        cell_id: CellId(1),
        material_id: MaterialId(1),
        inventories,
        fluxes: NeutronFluxSchedule::new(flux, vec![]).unwrap(),
    };
    let mut exact = irradiation_budget();
    exact.source.source.energy_tolerance_j = 0.0;
    let (_, report) = source
        .trial_source(0.0, 1.0, 300.0, exact, &mut || Ok(()))
        .unwrap();
    assert_eq!(report.fissions, 2.0);
    assert!(
        report.energy_roundoff_j >= 1.0,
        "the cell reduction must charge its lost joule"
    );
    let heating = IrradiationHeating::new(&initial, vec![source], exact).unwrap();
    assert!(matches!(
        heating.advance_to(&initial, 1.0, controls(1.0), || Ok(())),
        Err(FissionThermalError::EnergyBudget { .. })
    ));
    assert_eq!(
        heating.source(CellId(1)).unwrap().inventories[0].cumulative_fissions(),
        0.0
    );
}
