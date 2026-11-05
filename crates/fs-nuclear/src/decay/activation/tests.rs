use super::*;
use crate::decay::{DecayBranch, DecayModel, DecayNuclide};
use crate::irradiation::{GroupCrossSections, IrradiationData, IrradiationModel, NeutronGroups};
use crate::{DecayHeatGroup, FissionBudget, FissionModel};
use fs_qty::Frequency;
use fs_thermochem::decay::DecayLawV1;

fn energy(deposited_j: f64, radiation_j: f64, escaped_j: f64) -> EnergyChannels {
    EnergyChannels {
        deposited_j,
        radiation_j,
        escaped_j,
    }
}

fn budget() -> CaptureActivationBudget {
    CaptureActivationBudget {
        irradiation: IrradiationBudget {
            max_neutron_groups: 2,
            source: FissionBudget {
                max_decay_groups: 2,
                energy_tolerance_j: 1e-8,
                fuel_tolerance_nuclei: 1e-8,
            },
        },
        decay: DecayBudget {
            propagation: DecayLimits::default(),
            nuclei_tolerance: 1e-7,
            energy_tolerance_j: 1e-7,
        },
    }
}

fn fixture(k_f: f64, k_c: f64, product_rate: f64) -> (CaptureActivationState, GroupFlux) {
    let groups =
        Arc::new(NeutronGroups::new("synthetic-group-grid/v1".into(), vec![0.0, 1.0]).unwrap());
    let fission = Arc::new(
        FissionModel::new(
            "synthetic-fission-energy/v1".into(),
            energy(3.0, 2.0, 1.0),
            vec![DecayHeatGroup {
                decay_constant_s_inv: 0.4,
                energy_per_fission: energy(2.0, 1.0, 1.0),
            }],
        )
        .unwrap(),
    );
    let parent = IrradiationState::new(
        Arc::new(
            IrradiationModel::new(IrradiationData {
                data_id: "synthetic-reaction-coefficients/v1".into(),
                nuclide: NuclideId::new(10, 20, 0).unwrap(),
                groups: Arc::clone(&groups),
                cross_sections: vec![GroupCrossSections {
                    fission_m2: k_f,
                    capture_m2: k_c,
                }],
                temperature_min_k: 250.0,
                temperature_max_k: 1000.0,
                fission_model: fission,
                capture_energy_per_event: energy(1.0, 2.0, 0.0),
            })
            .unwrap(),
        ),
        100.0,
        vec![0.0],
    )
    .unwrap();
    let branches = if product_rate > 0.0 {
        vec![DecayBranch {
            parent: 0,
            daughter: 1,
            probability: 1.0,
            energy_per_decay: energy(5.0, 7.0, 11.0),
        }]
    } else {
        vec![]
    };
    let products = DecayInventory::new(
        Arc::new(
            DecayModel::new(
                "synthetic-product-decays/v1".into(),
                vec![
                    DecayNuclide {
                        id: "declared-capture-product".into(),
                        law: DecayLawV1::from_decay_constant(Frequency::new(product_rate)).unwrap(),
                    },
                    DecayNuclide {
                        id: "stable-descendant".into(),
                        law: DecayLawV1::stable(),
                    },
                ],
                branches,
                DecayLimits::default(),
            )
            .unwrap(),
        ),
        vec![0.0, 0.0],
    )
    .unwrap();
    let binding = CaptureProduct {
        node: 0,
        nuclide: NuclideId::new(10, 21, 1).unwrap(),
        data_id: "synthetic-radiative-capture-binding/v1".into(),
    };
    (
        CaptureActivationState::new(parent, products, binding).unwrap(),
        GroupFlux::new(groups, vec![1.0]).unwrap(),
    )
}

fn close(actual: f64, expected: f64) {
    if expected == 0.0 {
        assert_eq!(actual, 0.0);
    } else {
        assert!(
            (actual - expected).abs() < 5e-12 * expected.abs(),
            "actual {actual:.17e}, expected {expected:.17e}"
        );
    }
}

fn close_energy(actual: EnergyChannels, expected: EnergyChannels) {
    for (a, b) in channel_array(actual)
        .into_iter()
        .zip(channel_array(expected))
    {
        close(a, b);
    }
}

fn unchanged(actual: &CaptureActivationState, saved: &CaptureActivationState) {
    assert_eq!(actual.parent(), saved.parent());
    assert_eq!(actual.products(), saved.products());
    assert_eq!(actual.capture_product(), saved.capture_product());
    match (&actual.prepared, &saved.prepared) {
        (Some(a), Some(b)) => assert!(Arc::ptr_eq(a, b)),
        (None, None) => (),
        _ => panic!("a refused trial changed the accepted preparation cache"),
    }
}

#[test]
fn g1_equal_removal_and_product_rates_use_actual_capture_ages() {
    let (mut state, flux) = fixture(0.0, 0.5, 0.5);
    let report = state
        .advance(2.0, &flux, 300.0, budget(), || Ok(()))
        .unwrap();
    let survivor = 100.0 * (-1.0_f64).exp();
    let captures = 100.0 - survivor;
    let product_decays = captures - survivor;
    close(state.parent().remaining_parent_nuclei(), survivor);
    close(state.products().populations()[0], survivor);
    close(state.products().populations()[1], product_decays);
    close(report.irradiation.captures, captures);
    assert_eq!(report.irradiation.fissions, 0.0);
    close(report.products.produced_nuclei[0], captures);
    close(report.products.decay.decays[0], product_decays);
    close(state.products().cumulative_produced_nuclei()[0], captures);
    close_energy(
        report.released,
        energy(
            captures + 5.0 * product_decays,
            2.0 * captures + 7.0 * product_decays,
            11.0 * product_decays,
        ),
    );
    assert_eq!(state.capture_product().nuclide.isomer(), 1);
    // A uniform average of all captures gives a different product population.
    let averaged_source = captures * (1.0 - (-1.0_f64).exp());
    assert!((state.products().populations()[0] - averaged_source).abs() > 1.0);
}

#[test]
fn g1_distinct_capture_and_product_rates_match_the_independent_bateman_solution() {
    let (mut state, flux) = fixture(0.0, 0.4, 1.7);
    let time: f64 = 1.3;
    let parent = 100.0 * (-0.4 * time).exp();
    let daughter = 100.0 * 0.4 * ((-0.4 * time).exp() - (-1.7 * time).exp()) / (1.7 - 0.4);
    let captures = 100.0 - parent;
    let report = state
        .advance(time, &flux, 300.0, budget(), || Ok(()))
        .unwrap();
    close(state.parent().remaining_parent_nuclei(), parent);
    close(state.products().populations()[0], daughter);
    close(state.products().populations()[1], captures - daughter);
    close(report.products.decay.decays[0], captures - daughter);
    close(
        report.products.decay.released.deposited_j,
        5.0 * (captures - daughter),
    );
    close(report.irradiation.released.deposited_j, captures);
    assert!(report.population_roundoff_nuclei <= budget().decay.nuclei_tolerance);
    assert!(report.coupling_roundoff_nuclei <= report.population_roundoff_nuclei);
}

#[test]
fn g1_competing_fission_capture_and_delayed_energy_keep_disjoint_owners() {
    let (mut state, flux) = fixture(0.3, 0.7, 0.5);
    let time: f64 = 1.2;
    let parent = 100.0 * (-time).exp();
    let fissions = 0.3 * (100.0 - parent);
    let captures = 0.7 * (100.0 - parent);
    let product = 100.0 * 0.7 * ((-time).exp() - (-0.5 * time).exp()) / (0.5 - 1.0);
    let product_decays = captures - product;
    let retained_fission_events = 100.0 * 0.3 * ((-time).exp() - (-0.4 * time).exp()) / (0.4 - 1.0);
    let delayed_fission_events = fissions - retained_fission_events;
    let report = state
        .advance(time, &flux, 300.0, budget(), || Ok(()))
        .unwrap();
    close(report.irradiation.fissions, fissions);
    close(report.irradiation.captures, captures);
    close(state.products().populations()[0], product);
    close(state.products().populations()[1], product_decays);
    close(
        state.parent().stored_delayed_energy_j(),
        4.0 * retained_fission_events,
    );
    close(
        state.parent().remaining_parent_nuclei()
            + state.products().populations().iter().sum::<f64>()
            + state.parent().cumulative_fissions(),
        100.0,
    );
    close_energy(
        report.released,
        energy(
            3.0 * fissions + captures + 2.0 * delayed_fission_events + 5.0 * product_decays,
            2.0 * fissions + 2.0 * captures + delayed_fission_events + 7.0 * product_decays,
            fissions + delayed_fission_events + 11.0 * product_decays,
        ),
    );
    close_energy(
        report.products.decay.released,
        energy(
            5.0 * product_decays,
            7.0 * product_decays,
            11.0 * product_decays,
        ),
    );
    close(report.products.produced_nuclei.iter().sum(), captures);
}

#[test]
fn g1_zero_flux_stops_production_but_preserves_both_decay_histories() {
    let (mut state, flux) = fixture(0.2, 0.6, 1.0);
    let zero_flux = GroupFlux::new(Arc::clone(flux.groups()), vec![0.0]).unwrap();
    let idle = state
        .advance(0.5, &zero_flux, 300.0, budget(), || Ok(()))
        .unwrap();
    assert_eq!(idle.released, EnergyChannels::default());
    assert_eq!(state.parent().remaining_parent_nuclei(), 100.0);
    state
        .advance(1.0, &flux, 300.0, budget(), || Ok(()))
        .unwrap();
    let before = state.clone();
    let report = state
        .advance(2.0, &zero_flux, 300.0, budget(), || Ok(()))
        .unwrap();
    assert_eq!(report.irradiation.captures, 0.0);
    assert_eq!(report.irradiation.fissions, 0.0);
    assert_eq!(report.products.produced_nuclei, [0.0, 0.0]);
    assert_eq!(
        state.parent().cumulative_captures(),
        before.parent().cumulative_captures()
    );
    assert_eq!(
        state.parent().cumulative_fissions(),
        before.parent().cumulative_fissions()
    );
    close(
        state.products().populations()[0],
        before.products().populations()[0] * (-2.0_f64).exp(),
    );
    assert!(report.irradiation.released.deposited_j > 0.0);
    assert!(report.products.decay.released.deposited_j > 0.0);
    close(
        report.released.deposited_j,
        report.irradiation.released.deposited_j + report.products.decay.released.deposited_j,
    );
}

#[test]
fn g3_split_intervals_restart_and_exact_cache_keys_preserve_histories() {
    let (initial, flux) = fixture(0.3, 0.7, 0.8);
    let (whole, _) = initial
        .trial(1.0, &flux, 300.0, budget(), || Ok(()))
        .unwrap();
    let (mut split, _) = initial
        .trial(0.5, &flux, 300.0, budget(), || Ok(()))
        .unwrap();
    assert!(initial.prepared.is_none());
    let cached = Arc::clone(split.prepared.as_ref().unwrap());
    let mut restarted = split.clone();
    let report = split
        .advance(0.5, &flux, 400.0, budget(), || Ok(()))
        .unwrap();
    let replay = restarted
        .advance(0.5, &flux, 400.0, budget(), || Ok(()))
        .unwrap();
    assert_eq!(report, replay);
    unchanged(&split, &restarted);
    assert!(Arc::ptr_eq(split.prepared.as_ref().unwrap(), &cached));
    close(
        split.parent().remaining_parent_nuclei(),
        whole.parent().remaining_parent_nuclei(),
    );
    close(
        split.parent().cumulative_captures(),
        whole.parent().cumulative_captures(),
    );
    for (&a, &b) in split
        .products()
        .populations()
        .iter()
        .zip(whole.products().populations())
    {
        close(a, b);
    }
    close_energy(
        split.parent().cumulative_released(),
        whole.parent().cumulative_released(),
    );
    close_energy(
        split.products().cumulative_released(),
        whole.products().cumulative_released(),
    );
    split
        .advance(0.25, &flux, 400.0, budget(), || Ok(()))
        .unwrap();
    assert!(!Arc::ptr_eq(split.prepared.as_ref().unwrap(), &cached));
    let shorter = Arc::clone(split.prepared.as_ref().unwrap());
    let doubled_flux = GroupFlux::new(Arc::clone(flux.groups()), vec![2.0]).unwrap();
    split
        .advance(0.25, &doubled_flux, 400.0, budget(), || Ok(()))
        .unwrap();
    assert!(!Arc::ptr_eq(split.prepared.as_ref().unwrap(), &shorter));
    let changed_rate = Arc::clone(split.prepared.as_ref().unwrap());
    let mut new_limits = budget();
    new_limits.decay.propagation.maximum_work -= 1;
    split
        .advance(0.25, &doubled_flux, 400.0, new_limits, || Ok(()))
        .unwrap();
    assert!(!Arc::ptr_eq(
        split.prepared.as_ref().unwrap(),
        &changed_rate
    ));
}

#[test]
fn g4_binding_range_work_and_checkpoint_refusals_preserve_both_histories_and_cache() {
    let (mut state, flux) = fixture(0.3, 0.7, 0.8);
    for binding in [
        CaptureProduct {
            data_id: "".into(),
            ..state.capture_product().clone()
        },
        CaptureProduct {
            node: 2,
            ..state.capture_product().clone()
        },
        CaptureProduct {
            nuclide: NuclideId::new(11, 21, 0).unwrap(),
            ..state.capture_product().clone()
        },
        CaptureProduct {
            nuclide: NuclideId::new(10, 22, 0).unwrap(),
            ..state.capture_product().clone()
        },
    ] {
        assert!(
            CaptureActivationState::new(state.parent().clone(), state.products().clone(), binding)
                .is_err()
        );
    }
    state
        .advance(0.5, &flux, 300.0, budget(), || Ok(()))
        .unwrap();
    let saved = state.clone();
    let mut calls = 0;
    let (_, report) = state
        .trial(0.6, &flux, 300.0, budget(), || {
            calls += 1;
            Ok(())
        })
        .unwrap();
    assert!(report.population_roundoff_nuclei > 0.0);
    for stop in [1, calls / 2, calls] {
        let mut count = 0;
        let result = state.advance(0.6, &flux, 300.0, budget(), || {
            count += 1;
            if count == stop {
                Err(CaptureActivationError::Cancelled)
            } else {
                Ok(())
            }
        });
        assert_eq!(result, Err(CaptureActivationError::Cancelled));
        unchanged(&state, &saved);
    }
    let mut count = 0;
    assert_eq!(
        state.advance(0.6, &flux, 300.0, budget(), || {
            count += 1;
            if count == 2 {
                Err(CaptureActivationError::Invalid("original checkpoint error"))
            } else {
                Ok(())
            }
        }),
        Err(CaptureActivationError::Invalid("original checkpoint error"))
    );
    let mut limited = budget();
    limited.decay.propagation.maximum_work = 1;
    assert!(matches!(
        state.advance(0.6, &flux, 300.0, limited, || Ok(())),
        Err(CaptureActivationError::Propagation(_))
    ));
    unchanged(&state, &saved);
    let mut too_few_nodes = budget();
    too_few_nodes.decay.propagation.maximum_nodes = 3;
    assert!(matches!(
        state.advance(0.6, &flux, 300.0, too_few_nodes, || Ok(())),
        Err(CaptureActivationError::Propagation(_))
    ));
    let mut strict = budget();
    strict.decay.nuclei_tolerance = 0.0;
    assert!(state.advance(0.6, &flux, 300.0, strict, || Ok(())).is_err());
    strict = budget();
    strict.decay.energy_tolerance_j = 0.0;
    assert!(state.advance(0.6, &flux, 300.0, strict, || Ok(())).is_err());
    assert!(
        state
            .advance(0.6, &flux, 1200.0, budget(), || Ok(()))
            .is_err()
    );
    unchanged(&state, &saved);
    let (stable, stable_flux) = fixture(0.0, 1.0, 0.0);
    let huge_parent =
        IrradiationState::new(Arc::new(stable.parent().model().clone()), 1e200, vec![0.0]).unwrap();
    let mut huge = CaptureActivationState::new(
        huge_parent,
        stable.products().clone(),
        stable.capture_product().clone(),
    )
    .unwrap();
    let huge_before = huge.clone();
    let mut generous = budget();
    generous.irradiation.source.energy_tolerance_j = f64::MAX;
    generous.irradiation.source.fuel_tolerance_nuclei = f64::MAX;
    generous.decay.energy_tolerance_j = f64::MAX;
    generous.decay.nuclei_tolerance = f64::MAX;
    assert!(matches!(
        huge.advance(800.0, &stable_flux, 300.0, generous, || Ok(())),
        Err(CaptureActivationError::Propagation(
            fs_la::decay::DecayError::Unrepresentable(_)
        ))
    ));
    unchanged(&huge, &huge_before);
}

#[test]
fn g0_channel_export_defects_cannot_cancel_or_round_below_a_budget_boundary() {
    let large = 9_007_199_254_740_992.0; // 2^53: consecutive values are two apart.
    let mut defects = Accumulator::new(0.0);
    let released = join_released(
        energy(large, large, 0.0),
        energy(1.0, 3.0, 0.0),
        &mut defects,
    )
    .unwrap();
    assert_eq!(released, energy(large, large + 4.0, 0.0));
    // Total energy is exact while deposition loses one joule and radiation
    // gains one. Both destination errors must consume the allowance.
    assert_eq!(released.total_j(), 2.0 * large + 4.0);
    assert_eq!(defect_bound(defects).unwrap(), 2.0);
    let mut boundary = Accumulator::new(large);
    add(&mut boundary, 1.0).unwrap();
    assert_eq!(boundary.value(), large);
    assert_eq!(defect_bound(boundary).unwrap(), large.next_up());
    assert_eq!(defect_bound(Accumulator::new(0.0)).unwrap(), 0.0);
}

#[test]
fn g0_shared_decay_ledger_charges_each_emission_channel_before_accepting_heat() {
    let large = 9_007_199_254_740_992.0;
    let law = DecayLawV1::from_decay_constant(Frequency::new(core::f64::consts::LN_2)).unwrap();
    let model = Arc::new(
        DecayModel::new(
            "manufactured-channel-export/v1".into(),
            vec![
                DecayNuclide {
                    id: "A".into(),
                    law,
                },
                DecayNuclide {
                    id: "B".into(),
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
                    daughter: 2,
                    probability: 1.0,
                    energy_per_decay: energy(large, large, 0.0),
                },
                DecayBranch {
                    parent: 1,
                    daughter: 2,
                    probability: 1.0,
                    energy_per_decay: energy(1.0, 3.0, 0.0),
                },
            ],
            DecayLimits::default(),
        )
        .unwrap(),
    );
    let inventory = DecayInventory::new(model, vec![2.0, 2.0, 0.0]).unwrap();
    // Exact analytic one-half-life event integrals isolate the shared physical
    // ledger from propagator approximation: each parent has one transformation.
    let evolved = DecayAdvance {
        populations: vec![1.0, 1.0, 2.0],
        decays: vec![1.0, 1.0, 0.0],
        population_balance_residual: 0.0,
    };
    let mut allowance = budget().decay;
    allowance.energy_tolerance_j = 2.0;
    let (_, report) = inventory
        .account_evolution(evolved.clone(), 1.0, None, 0.0, allowance, &mut || Ok(()))
        .unwrap();
    assert_eq!(report.released, energy(large, large + 4.0, 0.0));
    assert_eq!(report.energy_roundoff_j, 2.0);
    allowance.energy_tolerance_j = 1.0;
    assert!(matches!(
        inventory.account_evolution(evolved, 1.0, None, 0.0, allowance, &mut || Ok(())),
        Err(DecayError::Budget("emission energy accounting roundoff"))
    ));
    assert_eq!(inventory.cumulative_released(), EnergyChannels::default());
}
