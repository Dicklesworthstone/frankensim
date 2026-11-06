//! Manufactured continuous-source laws and independent nonlinear solutions.

use super::*;
use crate::EnergyChannels;
use crate::burnup::{
    BurnupData, BurnupLimits, BurnupModel, BurnupNuclide, BurnupProduct, BurnupReaction,
};
use crate::diffusion::{
    DelayedGroup, DiffusionData, DiffusionLimits, DiffusionReaction, PrecursorOutflow,
    PrecursorTransport,
};
use crate::irradiation::{NeutronGroups, NuclideId};
use fs_qty::Frequency;
use fs_thermochem::decay::DecayLawV1;

fn budget() -> DiffusionBurnupBudget {
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
        reaction_tolerance: 1e-10,
        relaxation: 1.0,
    }
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:.17e} != {expected:.17e}, tolerance {tolerance}"
    );
}

fn nuclide(index: u16, lambda: f64) -> BurnupNuclide {
    BurnupNuclide {
        id: NuclideId::new(10, 25 + index, 0).unwrap(),
        law: DecayLawV1::from_decay_constant(Frequency::new(lambda)).unwrap(),
    }
}

fn state(
    nuclides: Vec<BurnupNuclide>,
    reactions: Vec<BurnupReaction>,
    populations: Vec<f64>,
    neutron_reactions: Vec<DiffusionReaction>,
    neutrons: f64,
    precursors: Option<f64>,
    transport: PrecursorTransport,
) -> DiffusionBurnupState {
    let groups =
        Arc::new(NeutronGroups::new("manufactured-source-group".into(), vec![0.0, 1.0]).unwrap());
    let nuclear = Arc::new(
        BurnupModel::new(
            BurnupData {
                data_id: "manufactured-continuous-spatial-feed".into(),
                groups: Arc::clone(&groups),
                temperature_min_k: 300.0,
                temperature_max_k: 600.0,
                nuclides,
                reactions,
            },
            budget().burnup.limits,
        )
        .unwrap(),
    );
    let bindings = (0..neutron_reactions.len())
        .map(|index| ReactionBinding {
            diffusion_reaction: index,
            burnup_reaction: index,
        })
        .collect();
    let model = Arc::new(
        DiffusionModel::new_with_transport(
            DiffusionData {
                data_id: "manufactured-source-spatial-neutrons".into(),
                groups,
                speeds_m_s: vec![1.0],
                cell_volumes_m3: vec![1.0],
                delayed_groups: if precursors.is_some() {
                    vec![DelayedGroup {
                        data_id: "manufactured-source-precursor".into(),
                        decay_constant_s_inv: 0.4,
                        spectrum: vec![1.0],
                    }]
                } else {
                    vec![]
                },
                reactions: neutron_reactions,
                links: vec![],
                leakage: vec![],
            },
            transport,
            budget().diffusion.limits,
        )
        .unwrap(),
    );
    DiffusionBurnupState::new(
        DiffusionState::new(model, vec![neutrons], precursors.into_iter().collect()).unwrap(),
        vec![BurnupInventory::new(nuclear, populations).unwrap()],
        bindings,
    )
    .unwrap()
}

fn capture(
    neutrons: f64,
    fuel: f64,
    precursors: Option<f64>,
    transport: PrecursorTransport,
) -> DiffusionBurnupState {
    state(
        vec![nuclide(0, 0.0), nuclide(1, 0.0)],
        vec![BurnupReaction {
            data_id: "manufactured-fed-capture".into(),
            parent: 0,
            rate: BurnupRate::Capture {
                cross_sections_m2: vec![0.01],
            },
            products: vec![BurnupProduct {
                nuclide: 1,
                yield_per_event: 1.0,
            }],
            energy_per_event: EnergyChannels {
                deposited_j: 1.0,
                radiation_j: 0.0,
                escaped_j: 0.0,
            },
        }],
        vec![fuel, 0.0],
        vec![DiffusionReaction {
            data_id: "manufactured-fed-capture".into(),
            cell: 0,
            incident_group: 0,
            macroscopic_cross_section_m_inv: fuel * 0.01,
            kind: DiffusionReactionKind::Capture,
        }],
        neutrons,
        precursors,
        transport,
    )
}

fn equilibrium() -> (DiffusionBurnupState, DiffusionBurnupSources) {
    let state = state(
        vec![
            nuclide(0, 0.0),
            nuclide(1, 0.4),
            nuclide(2, 0.0),
            nuclide(3, 0.0),
        ],
        vec![
            BurnupReaction {
                data_id: "manufactured-fed-fission".into(),
                parent: 0,
                rate: BurnupRate::Fission {
                    cross_sections_m2: vec![0.006],
                },
                products: vec![BurnupProduct {
                    nuclide: 1,
                    yield_per_event: 2.0,
                }],
                energy_per_event: EnergyChannels {
                    deposited_j: 8.0,
                    radiation_j: 2.0,
                    escaped_j: 1.0,
                },
            },
            BurnupReaction {
                data_id: "manufactured-fed-capture".into(),
                parent: 0,
                rate: BurnupRate::Capture {
                    cross_sections_m2: vec![0.004],
                },
                products: vec![BurnupProduct {
                    nuclide: 2,
                    yield_per_event: 1.0,
                }],
                energy_per_event: EnergyChannels {
                    deposited_j: 1.0,
                    radiation_j: 0.0,
                    escaped_j: 0.0,
                },
            },
            BurnupReaction {
                data_id: "manufactured-fed-fragment-decay".into(),
                parent: 1,
                rate: BurnupRate::Decay {
                    branch_fraction: 1.0,
                },
                products: vec![BurnupProduct {
                    nuclide: 3,
                    yield_per_event: 1.0,
                }],
                energy_per_event: EnergyChannels {
                    deposited_j: 2.0,
                    radiation_j: 1.0,
                    escaped_j: 0.5,
                },
            },
        ],
        vec![100.0, 0.0, 0.0, 0.0],
        vec![
            DiffusionReaction {
                data_id: "manufactured-fed-fission".into(),
                cell: 0,
                incident_group: 0,
                macroscopic_cross_section_m_inv: 0.6,
                kind: DiffusionReactionKind::Fission {
                    prompt_yields: vec![1.0],
                    delayed_yields: vec![],
                },
            },
            DiffusionReaction {
                data_id: "manufactured-fed-capture".into(),
                cell: 0,
                incident_group: 0,
                macroscopic_cross_section_m_inv: 0.4,
                kind: DiffusionReactionKind::Capture,
            },
        ],
        20.0,
        None,
        PrecursorTransport::default(),
    );
    (
        state,
        DiffusionBurnupSources {
            neutron_rates_s: vec![8.0],
            precursor_rates_s: vec![],
            nuclide_rates_s: vec![vec![20.0, 3.0, 0.0, 0.0]],
        },
    )
}

fn unchanged(actual: &DiffusionBurnupState, saved: &DiffusionBurnupState) {
    assert_eq!(
        actual.neutrons().neutron_counts(),
        saved.neutrons().neutron_counts()
    );
    assert_eq!(
        actual.neutrons().precursor_counts(),
        saved.neutrons().precursor_counts()
    );
    assert_eq!(actual.neutrons().model(), saved.neutrons().model());
    assert_eq!(actual.inventories(), saved.inventories());
    assert_eq!(actual.bindings(), saved.bindings());
}

#[test]
fn g1_continuous_sources_hold_equilibrium_and_include_same_interval_product_ages() {
    let (initial, sources) = equilibrium();
    let saved = initial.clone();
    let duration = 1.25;
    let (next, report) = initial
        .trial_with_sources(duration, &[300.0], &sources, budget(), || Ok(()))
        .unwrap();
    // Fuel: 20/s feed exactly balances 12 fissions/s + 8 captures/s.
    // One prompt neutron/fission returns its incident neutron; 8/s injection
    // replaces capture. Thus F=100 and n=20 are a coupled exact equilibrium.
    let fissions = 12.0 * duration;
    let captures = 8.0 * duration;
    let fragments = 27.0 / 0.4 * -(-0.4 * duration).exp_m1();
    let decays = 27.0 * duration - fragments;
    assert!(report.iterations > 1);
    assert!(report.reaction_defect <= budget().reaction_tolerance);
    close(next.neutrons().neutron_counts()[0], 20.0, 2e-9);
    close(next.inventories()[0].populations()[0], 100.0, 2e-9);
    close(next.inventories()[0].populations()[1], fragments, 3e-9);
    close(next.inventories()[0].populations()[2], captures, 2e-9);
    close(next.inventories()[0].populations()[3], decays, 3e-9);
    close(report.neutrons.fissions, fissions, 2e-9);
    close(report.neutrons.captures, captures, 2e-9);
    close(report.burnup[0].fissions, fissions, 2e-9);
    close(report.burnup[0].captures, captures, 2e-9);
    close(report.burnup[0].decays, decays, 3e-9);
    close(
        report.burnup[0].released.deposited_j,
        8.0 * fissions + captures + 2.0 * decays,
        3e-8,
    );
    close(
        report.burnup[0].released.radiation_j,
        2.0 * fissions + decays,
        1e-8,
    );
    close(
        report.burnup[0].released.escaped_j,
        fissions + 0.5 * decays,
        1e-8,
    );
    assert_eq!(report.neutrons.injected_neutron_counts, [10.0]);
    assert_eq!(report.produced_nuclei, [vec![25.0, 3.75, 0.0, 0.0]]);
    assert_eq!(report.produced_nuclei_by_cell, [28.75]);
    assert_eq!(
        next.inventories()[0].cumulative_produced_nuclei(),
        [25.0, 3.75, 0.0, 0.0]
    );
    unchanged(&initial, &saved);
}

#[test]
fn g1_feed_changes_an_empty_parents_mean_even_before_neutrons_arrive() {
    let initial = capture(0.0, 0.0, None, PrecursorTransport::default());
    let sources = DiffusionBurnupSources {
        neutron_rates_s: vec![0.0],
        precursor_rates_s: vec![],
        nuclide_rates_s: vec![vec![6.0, 0.0]],
    };
    let (fed, report) = initial
        .trial_with_sources(2.0, &[300.0], &sources, budget(), || Ok(()))
        .unwrap();
    assert_eq!(fed.inventories()[0].populations(), [12.0, 0.0]);
    assert_eq!(report.produced_nuclei_by_cell, [12.0]);
    assert_eq!(report.burnup[0].captures, 0.0);
    // The trial's mean is s*dt/2=6, rather than the zero initial value or
    // endpoint inventory 12. Its coefficient is still correct at zero flux.
    close(
        fed.neutrons().model().reactions()[0].macroscopic_cross_section_m_inv,
        0.06,
        1e-15,
    );
    let (cooled, after_feed) = fed.trial(0.5, &[300.0], budget(), || Ok(())).unwrap();
    assert_eq!(cooled.inventories()[0].populations(), [12.0, 0.0]);
    assert_eq!(
        cooled.inventories()[0].cumulative_produced_nuclei(),
        [12.0, 0.0]
    );
    assert_eq!(after_feed.produced_nuclei, [vec![0.0, 0.0]]);
    assert_eq!(after_feed.produced_nuclei_by_cell, [0.0]);
}

#[test]
fn g3_source_startup_refines_to_the_independent_nonlinear_tanh_solution() {
    // Equal neutron and parent sources from vacuum, capture only:
    // n=N, n'=s-sigma*n², n=sqrt(s/sigma)*tanh(sqrt(s*sigma)*t).
    let duration = 2.0;
    let feed = 10.0_f64;
    let sigma = 0.01_f64;
    let exponential = (-2.0 * (feed * sigma).sqrt() * duration).exp();
    let exact = (feed / sigma).sqrt() * (1.0 - exponential) / (1.0 + exponential);
    let sources = DiffusionBurnupSources {
        neutron_rates_s: vec![feed],
        precursor_rates_s: vec![],
        nuclide_rates_s: vec![vec![feed, 0.0]],
    };
    let mut errors = Vec::new();
    for steps in [2, 4, 8] {
        let mut evolving = capture(0.0, 0.0, None, PrecursorTransport::default());
        let mut captures = 0.0;
        for _ in 0..steps {
            let report = evolving
                .advance_with_sources(
                    duration / steps as f64,
                    &[300.0],
                    &sources,
                    budget(),
                    || Ok(()),
                )
                .unwrap();
            assert!(report.reaction_defect <= budget().reaction_tolerance);
            captures += report.burnup[0].captures;
        }
        let neutron_count = evolving.neutrons().neutron_counts()[0];
        let parent = evolving.inventories()[0].populations()[0];
        close(neutron_count, parent, 2e-8);
        close(parent + captures, feed * duration, 2e-8);
        close(evolving.inventories()[0].populations()[1], captures, 2e-8);
        assert_eq!(
            evolving.inventories()[0].cumulative_produced_nuclei(),
            [20.0, 0.0]
        );
        errors.push((parent - exact).abs());
        close(
            (captures - (feed * duration - exact)).abs(),
            errors[errors.len() - 1],
            2e-8,
        );
    }
    assert!(errors[1] < errors[0] / 1.8, "{errors:?}");
    assert!(errors[2] < errors[1] / 1.8, "{errors:?}");
}

#[test]
fn g1_precursor_sources_and_outflow_survive_every_composition_rebuild() {
    let transport = PrecursorTransport {
        transfers: vec![],
        outflows: vec![PrecursorOutflow {
            cell: 0,
            volumetric_flow_m3_s: 0.7,
        }],
    };
    let initial = capture(0.0, 0.0, Some(2.0), transport.clone());
    let sources = DiffusionBurnupSources {
        neutron_rates_s: vec![0.0],
        precursor_rates_s: vec![3.0],
        nuclide_rates_s: vec![vec![4.0, 0.0]],
    };
    let duration = 0.4;
    let (next, report) = initial
        .trial_with_sources(duration, &[300.0], &sources, budget(), || Ok(()))
        .unwrap();
    assert_eq!(next.neutrons().model().precursor_transport(), &transport);
    assert!(next.neutrons().model().reactions()[0].macroscopic_cross_section_m_inv > 0.0);
    let removal = 1.1;
    let loss_fraction = -(-removal * duration).exp_m1();
    let precursor = 2.0 * (1.0 - loss_fraction) + 3.0 / removal * loss_fraction;
    let residence =
        2.0 * loss_fraction / removal + 3.0 / removal * (duration - loss_fraction / removal);
    close(next.neutrons().precursor_counts()[0], precursor, 2e-9);
    close(report.neutrons.precursor_decays, 0.4 * residence, 2e-9);
    close(report.neutrons.escaped_precursors, 0.7 * residence, 2e-9);
    close(report.neutrons.injected_precursors, 3.0 * duration, 1e-14);
    close(
        next.neutrons().neutron_counts()[0] + report.burnup[0].captures,
        report.neutrons.precursor_decays,
        2e-9,
    );
    close(
        next.inventories()[0].populations()[0] + report.burnup[0].captures,
        4.0 * duration,
        2e-9,
    );
}

#[test]
fn g4_invalid_sources_nonconvergence_and_late_cancellation_preserve_all_histories() {
    let (initial, sources) = equilibrium();
    let mut bad_sources = Vec::new();
    let mut bad = sources.clone();
    bad.neutron_rates_s.clear();
    bad_sources.push(bad);
    let mut bad = sources.clone();
    bad.precursor_rates_s.push(0.0);
    bad_sources.push(bad);
    let mut bad = sources.clone();
    bad.nuclide_rates_s.clear();
    bad_sources.push(bad);
    let mut bad = sources.clone();
    bad.nuclide_rates_s[0].pop();
    bad_sources.push(bad);
    let mut bad = sources.clone();
    bad.neutron_rates_s[0] = -1.0;
    bad_sources.push(bad);
    let mut bad = sources.clone();
    bad.nuclide_rates_s[0][0] = f64::NAN;
    bad_sources.push(bad);
    let mut bad = sources.clone();
    bad.nuclide_rates_s[0][0] = f64::INFINITY;
    bad_sources.push(bad);
    for bad in bad_sources {
        assert!(bad.validate(&initial).is_err());
        let mut attempted = initial.clone();
        assert!(
            attempted
                .advance_with_sources(0.0, &[300.0], &bad, budget(), || Ok(()))
                .is_err()
        );
        unchanged(&attempted, &initial);
    }
    let mut strict = budget();
    strict.max_iterations = 1;
    let mut attempted = initial.clone();
    assert!(matches!(
        attempted.advance_with_sources(1.25, &[300.0], &sources, strict, || Ok(())),
        Err(DiffusionBurnupError::Nonconvergence { .. })
    ));
    unchanged(&attempted, &initial);
    let mut calls = 0;
    initial
        .trial_with_sources(1.25, &[300.0], &sources, budget(), || {
            calls += 1;
            Ok(())
        })
        .unwrap();
    for stop in [1, calls / 2, calls] {
        let mut attempted = initial.clone();
        let mut call = 0;
        assert!(
            attempted
                .advance_with_sources(1.25, &[300.0], &sources, budget(), || {
                    call += 1;
                    if call == stop {
                        Err(FissionError::Cancelled)
                    } else {
                        Ok(())
                    }
                })
                .is_err()
        );
        unchanged(&attempted, &initial);
    }
    let (same, zero) = initial
        .trial_with_sources(0.0, &[300.0], &sources, budget(), || Ok(()))
        .unwrap();
    assert_eq!(
        same.neutrons().neutron_counts(),
        initial.neutrons().neutron_counts()
    );
    assert_eq!(same.inventories(), initial.inventories());
    assert_eq!(zero.produced_nuclei_by_cell, [0.0]);
    assert_eq!(zero.neutrons.injected_neutrons, 0.0);
}

#[test]
fn g4_scalar_supply_exports_charge_lost_small_feeds_before_thermal_publication() {
    let large = 2.0_f64.powi(60);
    let initial = state(
        vec![nuclide(0, 0.0), nuclide(1, 0.0)],
        vec![],
        vec![0.0, 0.0],
        vec![],
        0.0,
        None,
        PrecursorTransport::default(),
    );
    let sources = DiffusionBurnupSources {
        neutron_rates_s: vec![0.0],
        precursor_rates_s: vec![],
        nuclide_rates_s: vec![vec![large, 1.0]],
    };
    let mut allowance = budget();
    allowance.burnup.nuclei_tolerance = 1.0;
    let (next, report) = initial
        .trial_with_sources(1.0, &[300.0], &sources, allowance, || Ok(()))
        .unwrap();
    assert_eq!(next.inventories()[0].populations(), [large, 1.0]);
    assert_eq!(report.produced_nuclei, [vec![large, 1.0]]);
    assert_eq!(report.produced_nuclei_by_cell, [large]);
    assert_eq!(report.burnup[0].population_roundoff_nuclei, 1.0);
    allowance.burnup.nuclei_tolerance = 0.5;
    let mut refused = initial.clone();
    assert!(matches!(
        refused.advance_with_sources(1.0, &[300.0], &sources, allowance, || Ok(())),
        Err(DiffusionBurnupError::Burnup(BurnupError::Budget(_)))
    ));
    unchanged(&refused, &initial);
}
