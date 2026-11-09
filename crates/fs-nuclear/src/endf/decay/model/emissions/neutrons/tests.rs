//! Manufactured ENDF and analytic source/absorption references, not isotope data.

use super::super::super::{
    DecayEnergyFractions, DecayEnergyRouting, DecayModelEntry, EvaluatedDecayAssembly,
};
use super::super::DecayEmissionLimits;
use super::*;
use crate::EnergyChannels;
use crate::burnup::{
    BurnupBudget, BurnupData, BurnupInventory, BurnupLimits, BurnupRate, BurnupReaction,
};
use crate::decay::{DecayBudget, DecayInventory};
use crate::diffusion::decay::{DecayNeutronBudget, DecayNeutronSource, DecayNeutronState};
use crate::diffusion::{
    DiffusionBudget, DiffusionData, DiffusionLimits, DiffusionModel, DiffusionReaction,
    DiffusionReactionKind, DiffusionState,
};
use crate::endf::decay::{DecayEvaluation, EndfDecayReadLimits, EndfDecaySelection};
use crate::irradiation::{GroupFlux, NeutronGroups};
use fs_la::decay::DecayLimits;

const EV: f64 = 1.602_176_634e-19;

fn close(actual: f64, expected: f64, relative: f64) {
    assert!(
        (actual - expected).abs() <= relative * expected.abs().max(1e-300),
        "{actual:.17e} != {expected:.17e}"
    );
}

// All fixture numbers fit directly in an eleven-column ENDF numeric field.
fn row(values: [f64; 6], mt: u16) -> String {
    let fields: String = values
        .into_iter()
        .map(|value| format!("{value:>11}"))
        .collect();
    assert_eq!(fields.len(), 66);
    format!("{fields}{:4}{:2}{mt:3}{:5}", 777, 8, 1)
}

fn evaluation(
    stable: bool,
    spectrum: bool,
    neutron_mode: bool,
    continuum: bool,
    heavy_ev: f64,
) -> Arc<DecayEvaluation> {
    let id = if stable {
        NuclideId::new(11, if neutron_mode { 29 } else { 30 }, 0).unwrap()
    } else {
        NuclideId::new(10, 30, 0).unwrap()
    };
    let mode = if neutron_mode { 1.5 } else { 1.0 };
    let mut rows = vec![
        row(
            [
                1000.0 * f64::from(id.atomic_number()) + f64::from(id.mass_number()),
                f64::from(id.mass_number()),
                0.0,
                0.0,
                if stable { 1.0 } else { 0.0 },
                if spectrum && !stable { 1.0 } else { 0.0 },
            ],
            457,
        ),
        row(
            [if stable { 0.0 } else { 2.0 }, 0.0, 0.0, 0.0, 6.0, 0.0],
            457,
        ),
        row(
            [
                0.0,
                0.0,
                if stable { 0.0 } else { 10.0 },
                0.0,
                if stable { 0.0 } else { heavy_ev },
                0.0,
            ],
            457,
        ),
        row(
            [
                0.0,
                0.0,
                0.0,
                0.0,
                if stable { 6.0 } else { 12.0 },
                if stable { 0.0 } else { 2.0 },
            ],
            457,
        ),
    ];
    if stable {
        rows.push(row([0.0; 6], 457));
    } else {
        // Repeated residual modes exercise parent-wide, not per-branch, spectra.
        rows.push(row([mode, 0.0, 20.0, 0.0, 0.25, 0.0], 457));
        rows.push(row([mode, 0.0, 20.0, 0.0, 0.75, 0.0], 457));
        if spectrum {
            rows.push(row(
                [
                    0.0,
                    5.0,
                    if continuum { 1.0 } else { 0.0 },
                    0.0,
                    6.0,
                    if continuum { 0.0 } else { 1.0 },
                ],
                457,
            ));
            rows.push(row(
                [
                    if continuum { 0.0 } else { 1.0 },
                    0.0,
                    3.0,
                    0.0,
                    if continuum { 1.0 } else { 0.0 },
                    0.0,
                ],
                457,
            ));
            if continuum {
                // RP(E)=(E-1)/2 per eV on [1,3] eV: count=1, energy=7/3 eV.
                rows.push(row([mode, 0.0, 0.0, 0.0, 1.0, 2.0], 457));
                rows.push(row([2.0, 2.0, 0.0, 0.0, 0.0, 0.0], 457));
                rows.push(row([1.0, 0.0, 3.0, 1.0, 0.0, 0.0], 457));
            } else {
                rows.push(row([3.0, 0.0, 0.0, 0.0, 6.0, 0.0], 457));
                rows.push(row([mode, 0.0, 1.0, 0.0, 0.0, 0.0], 457));
            }
        }
    }
    rows.push(row([0.0; 6], 0));
    Arc::new(
        DecayEvaluation::parse(
            "manufactured-neutron-MF8".into(),
            EndfDecaySelection {
                material_number: 777,
                nuclide: id,
            },
            &rows.join("\n"),
            EndfDecayReadLimits::default(),
            || Ok(()),
        )
        .unwrap(),
    )
}

fn source(
    spectrum: bool,
    neutron_mode: bool,
    continuum: bool,
    heavy_ev: f64,
    heavy_outward: bool,
) -> Arc<PreparedDecayEmission> {
    let outward = DecayEnergyFractions {
        deposited: 0.0,
        radiation: 1.0,
        escaped: 0.0,
    };
    let routing = DecayEnergyRouting {
        light_particles: outward,
        electromagnetic: outward,
        heavy_particles: if heavy_outward {
            outward
        } else {
            DecayEnergyFractions {
                deposited: 1.0,
                radiation: 0.0,
                escaped: 0.0,
            }
        },
        include_reported_neutrinos: false,
    };
    // Intentionally retain stable first and parent second, rather than sorting.
    let assembly = Arc::new(
        EvaluatedDecayAssembly::new(
            "manufactured-neutron-assembly".into(),
            [true, false]
                .into_iter()
                .map(|stable| DecayModelEntry {
                    evaluation: evaluation(stable, spectrum, neutron_mode, continuum, heavy_ev),
                    routing,
                })
                .collect(),
            DecayLimits::default(),
            || Ok(()),
        )
        .unwrap(),
    );
    Arc::new(
        PreparedDecayEmission::new(
            assembly,
            DecayRadiation::Neutron,
            vec![0.0, 2.0 * EV, 4.0 * EV],
            DecayEmissionLimits::default(),
            || Ok(()),
        )
        .unwrap(),
    )
}

fn kinetic(energies: &[f64], capture: f64) -> Arc<NeutronKineticModel> {
    let mut edges = vec![0.0];
    for pair in energies.windows(2) {
        edges.push(0.5 * (pair[0] + pair[1]));
    }
    edges.push(energies[energies.len() - 1] * 2.0);
    let diffusion = Arc::new(
        DiffusionModel::new(
            DiffusionData {
                data_id: "manufactured-neutron-absorber".into(),
                groups: Arc::new(
                    NeutronGroups::new("manufactured-neutron-grid".into(), edges).unwrap(),
                ),
                speeds_m_s: vec![1.0; energies.len()],
                cell_volumes_m3: vec![1.0],
                delayed_groups: vec![],
                reactions: (0..energies.len())
                    .map(|incident_group| DiffusionReaction {
                        data_id: format!("capture-{incident_group}"),
                        cell: 0,
                        incident_group,
                        macroscopic_cross_section_m_inv: capture,
                        kind: DiffusionReactionKind::Capture,
                    })
                    .collect(),
                links: vec![],
                leakage: vec![],
            },
            DiffusionLimits::default(),
        )
        .unwrap(),
    );
    Arc::new(
        NeutronKineticModel::new(
            "manufactured-neutron-energy".into(),
            diffusion,
            energies.to_vec(),
            vec![1.0; energies.len()],
        )
        .unwrap(),
    )
}

fn bind(
    source: &Arc<PreparedDecayEmission>,
    model: Arc<NeutronKineticModel>,
    zero: &[usize],
) -> Result<EvaluatedDecayNeutrons, DecayNeutronBindingError> {
    source.bind_neutron_group_means(
        "manufactured-projection-and-zero-policy".into(),
        model,
        zero,
        || Ok(()),
    )
}

#[test]
fn parsed_delayed_neutrons_preserve_parent_counts_energy_and_true_birth_history() {
    let original = source(true, true, false, 4.0, true);
    let transport = kinetic(&[2.0 * EV, 4.0 * EV], 2.0);
    let binding = bind(&original, Arc::clone(&transport), &[]).unwrap();
    assert!(Arc::ptr_eq(binding.source(), &original));
    assert!(Arc::ptr_eq(binding.kinetic_model(), &transport));
    assert!(Arc::ptr_eq(
        binding.emission().model(),
        original.assembly().model()
    ));
    assert!(Arc::ptr_eq(
        binding.emission().groups(),
        transport.diffusion().groups()
    ));
    assert_eq!(binding.emission().neutron_yields()[0], vec![0.0; 2]);
    for yield_ in &binding.emission().neutron_yields()[1] {
        close(*yield_, 0.5, 2e-15);
    }
    close(binding.projection()[1].particles_per_decay, 1.0, 2e-15);
    close(
        binding.projection()[1].kinetic_energy_j_per_decay,
        3.0 * EV,
        2e-15,
    );
    let initial = DecayNeutronState::new(
        DiffusionState::new(Arc::clone(transport.diffusion()), vec![0.0; 2], vec![]).unwrap(),
        vec![DecayNeutronSource {
            cell: 0,
            inventory: DecayInventory::new(
                Arc::clone(original.assembly().model()),
                vec![0.0, 100.0],
            )
            .unwrap(),
            emission: Arc::clone(binding.emission()),
        }],
    )
    .unwrap();
    let budget = DecayNeutronBudget {
        diffusion: DiffusionBudget {
            limits: DiffusionLimits::default(),
            population_tolerance: 1e-8,
        },
        decay: DecayBudget {
            propagation: DecayLimits::default(),
            nuclei_tolerance: 1e-8,
            energy_tolerance_j: 1e-25,
        },
    };
    let t = 2.0_f64;
    let lambda = std::f64::consts::LN_2 / 2.0;
    let (next, report) = initial.trial(t, budget, || Ok(())).unwrap();
    close(report.decay[0].decays[1], 50.0, 3e-11);
    let surviving = 100.0 * lambda / (2.0 - lambda) * ((-lambda * t).exp() - (-2.0 * t).exp());
    for &count in next.neutrons().neutron_counts() {
        close(count, 0.5 * surviving, 3e-11);
    }
    close(report.neutrons.captures, 50.0 - surviving, 3e-11);
    close(report.emitted_neutrons[0].iter().sum(), 50.0, 3e-11);
}

#[test]
fn continuum_bin_moments_condense_without_midpoint_or_second_branch_weight() {
    let original = source(true, true, true, 4.0, true);
    let bound = bind(&original, kinetic(&[1.0 * EV, 3.0 * EV], 0.0), &[]).unwrap();
    let weights = &bound.emission().neutron_yields()[1];
    close(weights[0], 1.0 / 3.0, 3e-15);
    close(weights[1], 2.0 / 3.0, 3e-15);
    close(
        bound.projection()[1].kinetic_energy_j_per_decay,
        (7.0 / 3.0) * EV,
        3e-15,
    );
    assert!(bound.projection()[1].energy_projection_residual_j.abs() < 1e-32);
}

#[test]
fn unreported_is_not_zero_and_confirmation_cannot_contradict_decay_modes() {
    let missing = source(false, false, false, 0.0, true);
    let model = kinetic(&[2.0 * EV, 4.0 * EV], 0.0);
    assert!(matches!(
        bind(&missing, Arc::clone(&model), &[]),
        Err(DecayNeutronBindingError::Unreported { parent: 1, .. })
    ));
    let confirmed = bind(&missing, Arc::clone(&model), &[1]).unwrap();
    assert_eq!(confirmed.confirmed_nonemitting_parents(), &[1]);
    assert!(
        confirmed
            .emission()
            .neutron_yields()
            .iter()
            .flatten()
            .all(|&x| x == 0.0)
    );
    for invalid in [&[0][..], &[2][..], &[1, 1][..]] {
        assert!(matches!(
            bind(&missing, Arc::clone(&model), invalid),
            Err(DecayNeutronBindingError::Invalid(_))
        ));
    }
    let neutron_missing = source(false, true, false, 4.0, true);
    assert!(matches!(
        bind(&neutron_missing, model, &[1]),
        Err(DecayNeutronBindingError::Invalid(_))
    ));
}

#[test]
fn projection_refuses_wrong_species_outside_means_and_reused_heat_or_gamma_energy() {
    let original = source(true, true, false, 4.0, true);
    assert!(matches!(
        bind(&original, kinetic(&[4.0 * EV, 5.0 * EV], 0.0), &[]),
        Err(DecayNeutronBindingError::OutsideRepresentativeEnergies { parent: 1, .. })
    ));
    for original in [
        source(true, true, false, 2.0, true),
        source(true, true, false, 4.0, false),
    ] {
        // Ten eV of gamma radiation cannot fund a three-eV neutron.
        assert!(matches!(
            bind(&original, kinetic(&[2.0 * EV, 4.0 * EV], 0.0), &[]),
            Err(DecayNeutronBindingError::RadiationBudget { parent: 1, .. })
        ));
    }
    let photons = Arc::new(
        PreparedDecayEmission::new(
            Arc::clone(original.assembly()),
            DecayRadiation::Gamma,
            vec![0.0, 20.0 * EV],
            DecayEmissionLimits::default(),
            || Ok(()),
        )
        .unwrap(),
    );
    assert!(matches!(
        bind(&photons, kinetic(&[2.0 * EV, 4.0 * EV], 0.0), &[1]),
        Err(DecayNeutronBindingError::Invalid(_))
    ));
}

#[test]
fn preparation_is_bounded_and_late_cancellation_preserves_immutable_inputs() {
    let original = source(true, true, false, 4.0, true);
    let model = kinetic(&[2.0 * EV, 4.0 * EV], 0.0);
    let mut calls = 0;
    let cancelled =
        original.bind_neutron_group_means("cancelled".into(), Arc::clone(&model), &[], || {
            calls += 1;
            if calls >= 3 {
                Err(FissionError::Cancelled)
            } else {
                Ok(())
            }
        });
    assert!(matches!(
        cancelled,
        Err(DecayNeutronBindingError::Source(EndfError::Physics(
            FissionError::Cancelled
        )))
    ));
    assert!(bind(&original, Arc::clone(&model), &[]).is_ok());
    let mut bounded = original.as_ref().clone();
    bounded.limits.max_group_values = 1;
    assert!(matches!(
        bind(&Arc::new(bounded), model, &[]),
        Err(DecayNeutronBindingError::Source(EndfError::Physics(
            FissionError::Budget(_)
        )))
    ));
}

fn burnup_data(source: &EvaluatedDecayNeutrons) -> BurnupData {
    source
        .source()
        .assembly()
        .burnup_data(
            Arc::clone(source.kinetic_model().diffusion().groups()),
            250.0,
            1500.0,
            || Ok(()),
        )
        .unwrap()
}

#[test]
fn burnup_binding_retains_evaluation_with_reordered_decay_and_finite_capture() {
    let original = source(true, true, false, 4.0, true);
    let projected = Arc::new(bind(&original, kinetic(&[2.0 * EV, 4.0 * EV], 0.0), &[]).unwrap());
    let mut data = burnup_data(&projected);
    data.reactions.push(BurnupReaction {
        data_id: "manufactured-emitter-capture".into(),
        parent: 1,
        rate: BurnupRate::Capture {
            cross_sections_m2: vec![0.2, 0.4],
        },
        // Capture products explicitly leave this finite tracked network.
        products: vec![],
        energy_per_event: EnergyChannels::default(),
    });
    data.reactions.reverse();
    let model = Arc::new(BurnupModel::new(data, BurnupLimits::default()).unwrap());
    let bound = projected
        .bind_burnup(
            "evaluated-finite-emitter".into(),
            Arc::clone(&model),
            || Ok(()),
        )
        .unwrap();
    assert!(Arc::ptr_eq(bound.source(), &projected));
    assert!(Arc::ptr_eq(bound.burnup_binding().model(), &model));
    assert!(Arc::ptr_eq(bound.emission().model(), &model));
    assert_eq!(bound.emission().data_id(), "evaluated-finite-emitter");
    assert_eq!(
        bound.emission().neutron_yields(),
        projected.emission().neutron_yields()
    );
    let initial = BurnupInventory::new(Arc::clone(&model), vec![0.0, 100.0]).unwrap();
    let flux = GroupFlux::new(Arc::clone(model.groups()), vec![1.0, 1.0]).unwrap();
    let (_, report) = initial
        .trial(
            1.0,
            &flux,
            300.0,
            BurnupBudget {
                limits: BurnupLimits::default(),
                nuclei_tolerance: 1e-8,
                energy_tolerance_j: 1e-25,
            },
            || Ok(()),
        )
        .unwrap();
    let emitted = bound
        .burnup_binding()
        .from_burnup_step(&model, &report, || Ok(()))
        .unwrap();
    let lambda = std::f64::consts::LN_2 / 2.0;
    let expected_decays = 100.0 * lambda / (lambda + 0.6) * -(-(lambda + 0.6)).exp_m1();
    close(report.decays, expected_decays, 3e-11);
    assert!(report.captures > 0.0);
    close(
        emitted.emitted_particles.iter().sum(),
        expected_decays,
        3e-11,
    );
    close(
        emitted.emitted_energy_j.iter().sum(),
        3.0 * EV * expected_decays,
        3e-11,
    );
}

#[test]
fn burnup_binding_refuses_changed_identity_energy_grid_and_late_cancellation() {
    let original = source(true, true, false, 4.0, true);
    let projected = Arc::new(bind(&original, kinetic(&[2.0 * EV, 4.0 * EV], 0.0), &[]).unwrap());
    let data = burnup_data(&projected);
    let mut wrong_identity = data.clone();
    wrong_identity.nuclides[0].id = NuclideId::new(11, 28, 0).unwrap();
    let mut wrong_energy = data.clone();
    wrong_energy.reactions[0].energy_per_event.radiation_j += EV;
    for changed in [wrong_identity, wrong_energy] {
        let model = Arc::new(BurnupModel::new(changed, BurnupLimits::default()).unwrap());
        assert!(matches!(
            projected.bind_burnup("incompatible".into(), model, || Ok(())),
            Err(DecayNeutronBindingError::Source(_))
        ));
    }
    let mut wrong_grid = data.clone();
    wrong_grid.groups = Arc::new(
        NeutronGroups::new(
            "different-neutron-grid".into(),
            data.groups.edges_j().to_vec(),
        )
        .unwrap(),
    );
    let model = Arc::new(BurnupModel::new(wrong_grid, BurnupLimits::default()).unwrap());
    assert!(matches!(
        projected.bind_burnup("incompatible-grid".into(), model, || Ok(())),
        Err(DecayNeutronBindingError::Invalid(_))
    ));
    let model = Arc::new(BurnupModel::new(data, BurnupLimits::default()).unwrap());
    let mut calls = 0;
    let cancelled = projected.bind_burnup("cancelled".into(), Arc::clone(&model), || {
        calls += 1;
        if calls >= 4 {
            Err(FissionError::Cancelled)
        } else {
            Ok(())
        }
    });
    assert_eq!(calls, 4);
    assert!(matches!(
        cancelled,
        Err(DecayNeutronBindingError::Source(EndfError::Physics(
            FissionError::Cancelled
        )))
    ));
    assert!(
        projected
            .bind_burnup("accepted".into(), model, || Ok(()))
            .is_ok()
    );
}
