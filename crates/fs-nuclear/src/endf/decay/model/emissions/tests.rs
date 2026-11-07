//! G1 analytic spectra/populations and G3 interval/group/owner agreement.
//! Except for the pinned Co60/Ni60 records, all coefficients are manufactured.

use super::super::{
    DecayEnergyFractions, DecayEnergyRouting, DecayModelEntry, DecayProductRetention,
};
use super::*;
use crate::burnup::transport::{
    BurnupTransportBudget, BurnupTransportLimits, BurnupTransportState, NuclideTransfer,
    NuclideTransport,
};
use crate::burnup::{BurnupBudget, BurnupInventory, BurnupLimits, BurnupReaction};
use crate::decay::{DecayBudget, DecayInventory};
use crate::endf::decay::spectra::{
    DecayContinuum, DecayInterpolation, DiscreteDecayEmission, SpectrumValue,
};
use crate::endf::decay::{
    DecayEnergy, DecayEvaluation, DecayMetadata, DecayMode, EndfDecayReadLimits, EndfDecaySelection,
};
use crate::irradiation::{GroupFlux, NeutronGroups};
use fs_la::decay::DecayLimits;
use fs_qty::Frequency;
use fs_thermochem::decay::DecayLawV1;

const U: f64 = 1e-13;

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:.17e} != {expected:.17e}, tolerance {tolerance:.3e}"
    );
}

fn value(value: f64) -> SpectrumValue {
    SpectrumValue {
        value,
        uncertainty: 0.0,
    }
}

fn line(energy_j: f64, intensity: f64) -> DiscreteDecayEmission {
    DiscreteDecayEmission {
        energy_j: value(energy_j),
        decay_mode: 1.0,
        transition_type: 0.0,
        intensity: value(intensity),
        pair_or_positron: None,
        conversion_total: None,
        conversion_k: None,
        conversion_l: None,
    }
}

fn gamma(energy_j: f64) -> DecaySpectrum {
    DecaySpectrum {
        radiation: DecayRadiation::Gamma,
        discrete_normalization: value(1.0),
        average_energy_j: DecayEnergy {
            value_j: Some(energy_j),
            uncertainty_j: 0.0,
        },
        continuum_normalization: value(0.0),
        discrete: vec![line(energy_j, 1.0)].into(),
        continuum: None,
        covariance_flag: 0,
    }
}

fn routing() -> DecayEnergyRouting {
    DecayEnergyRouting {
        light_particles: DecayEnergyFractions {
            deposited: 1.0,
            radiation: 0.0,
            escaped: 0.0,
        },
        electromagnetic: DecayEnergyFractions {
            deposited: 0.0,
            radiation: 1.0,
            escaped: 0.0,
        },
        heavy_particles: DecayEnergyFractions {
            deposited: 0.0,
            radiation: 0.0,
            escaped: 1.0,
        },
        include_reported_neutrinos: false,
    }
}

fn assemble(entries: Vec<DecayModelEntry>) -> Arc<EvaluatedDecayAssembly> {
    Arc::new(
        EvaluatedDecayAssembly::new(
            "manufactured-emission-assembly".into(),
            entries,
            DecayLimits::default(),
            || Ok(()),
        )
        .unwrap(),
    )
}

/// A beta-minus chain with equal one-second half-lives and a stable endpoint.
fn chain(spectra: Vec<Vec<DecaySpectrum>>, split_first: bool) -> Arc<EvaluatedDecayAssembly> {
    let n = spectra.len();
    assemble(
        spectra
            .into_iter()
            .enumerate()
            .map(|(index, spectra)| {
                let stable = index + 1 == n;
                let fractions = if stable {
                    vec![]
                } else if index == 0 && split_first {
                    vec![0.25, 0.75]
                } else {
                    vec![1.0]
                };
                let evaluation = DecayEvaluation {
                    data_id: format!("manufactured-emission-parent-{index}").into(),
                    metadata: DecayMetadata {
                        material_number: 700 + index as u16,
                        nuclide: NuclideId::new(10 + index as u16, 30, 0).unwrap(),
                        atomic_weight_ratio: 30.0,
                        excited_state: 0,
                        spectrum_count: spectra.len(),
                        spin: None,
                        parity: None,
                    },
                    law: if stable {
                        DecayLawV1::stable()
                    } else {
                        DecayLawV1::from_decay_constant(Frequency::new(std::f64::consts::LN_2))
                            .unwrap()
                    },
                    half_life_uncertainty_s: 0.0,
                    average_energies: [0.0, 2.0 * U, 0.0]
                        .map(|v| DecayEnergy {
                            value_j: Some(v),
                            uncertainty_j: 0.0,
                        })
                        .into(),
                    modes: fractions
                        .into_iter()
                        .map(|probability| DecayMode {
                            rtyp: 1.0,
                            daughter_isomer: 0,
                            q_j: 5.0 * U,
                            q_uncertainty_j: 0.0,
                            probability,
                            probability_uncertainty: 0.0,
                        })
                        .collect::<Vec<_>>()
                        .into(),
                    branch_fraction_sum: if stable { 0.0 } else { 1.0 },
                    spectra: spectra.into(),
                };
                DecayModelEntry {
                    evaluation: Arc::new(evaluation),
                    routing: routing(),
                }
            })
            .collect(),
    )
}

fn prepare(
    assembly: Arc<EvaluatedDecayAssembly>,
    kind: DecayRadiation,
    edges: Vec<f64>,
) -> PreparedDecayEmission {
    PreparedDecayEmission::new(
        assembly,
        kind,
        edges,
        DecayEmissionLimits::default(),
        || Ok(()),
    )
    .unwrap()
}

fn decay_budget() -> DecayBudget {
    DecayBudget {
        propagation: DecayLimits::default(),
        nuclei_tolerance: 1e-6,
        energy_tolerance_j: 1e-15,
    }
}

fn inventory(assembly: &EvaluatedDecayAssembly) -> DecayInventory {
    let mut populations = vec![0.0; assembly.nuclides().len()];
    populations[0] = 1000.0;
    DecayInventory::new(Arc::clone(assembly.model()), populations).unwrap()
}

fn group_grid() -> Arc<NeutronGroups> {
    Arc::new(NeutronGroups::new("manufactured-incident-grid".into(), vec![0.0, 1.0]).unwrap())
}

fn burnup(assembly: &EvaluatedDecayAssembly, capture: bool) -> Arc<BurnupModel> {
    let mut data = assembly
        .burnup_data(group_grid(), 200.0, 1000.0, || Ok(()))
        .unwrap();
    data.reactions.reverse();
    if capture {
        // An explicit untracked capture product: removal competes with decay,
        // but neither the removal count nor its independently supplied energy
        // belongs to this evaluated radioactive photon source.
        data.reactions.insert(
            0,
            BurnupReaction {
                data_id: "manufactured-competing-capture".into(),
                parent: 0,
                rate: BurnupRate::Capture {
                    cross_sections_m2: vec![std::f64::consts::LN_2],
                },
                products: vec![],
                energy_per_event: crate::EnergyChannels::default(),
            },
        );
    }
    Arc::new(BurnupModel::new(data, BurnupLimits::default()).unwrap())
}

fn burnup_budget() -> BurnupBudget {
    BurnupBudget {
        limits: BurnupLimits::default(),
        nuclei_tolerance: 1e-6,
        energy_tolerance_j: 1e-15,
    }
}

#[test]
fn analytic_lines_linear_and_histogram_continua_preserve_counts_and_energy_moments() {
    let mut spectrum = gamma(U);
    spectrum.discrete_normalization = value(0.5);
    spectrum.discrete = vec![line(U, 2.0), line(2.0 * U, 4.0), line(4.0 * U, 2.0)].into();
    spectrum.continuum_normalization = value(2.0);
    spectrum.continuum = Some(DecayContinuum {
        decay_mode: 1.0,
        points: [
            (0.0, 0.0),
            (2.0 * U, 1.0 / U),
            (2.0 * U, 2.0 / U),
            (4.0 * U, 0.0),
        ]
        .map(|(energy_j, density_per_j)| DecaySpectrumPoint {
            energy_j,
            density_per_j,
        })
        .into(),
        interpolation: vec![
            DecayInterpolation {
                end_point: 2,
                law: 2,
            },
            DecayInterpolation {
                end_point: 4,
                law: 1,
            },
        ]
        .into(),
    });
    let assembly = chain(vec![vec![spectrum], vec![]], true);
    let coarse = prepare(
        Arc::clone(&assembly),
        DecayRadiation::Gamma,
        vec![0.0, U, 2.0 * U, 4.0 * U],
    );
    for (&actual, expected) in coarse
        .particles_per_decay(0)
        .unwrap()
        .iter()
        .zip([0.5, 2.5, 11.0])
    {
        near(actual, expected, 1e-13);
    }
    for (&actual, expected) in
        coarse
            .energy_per_decay_j(0)
            .unwrap()
            .iter()
            .zip([U / 3.0, 10.0 * U / 3.0, 32.0 * U])
    {
        near(actual, expected, 1e-26);
    }
    // Refined groups sum back to the coarse count AND first energy moment.
    let fine = prepare(
        Arc::clone(&assembly),
        DecayRadiation::Gamma,
        vec![0.0, U, 2.0 * U, 3.0 * U, 4.0 * U],
    );
    near(
        fine.particles_per_decay(0).unwrap()[2] + fine.particles_per_decay(0).unwrap()[3],
        11.0,
        1e-13,
    );
    near(
        fine.energy_per_decay_j(0).unwrap()[2] + fine.energy_per_decay_j(0).unwrap()[3],
        32.0 * U,
        1e-26,
    );
    let mut state = inventory(&assembly);
    let step = state.advance(1.0, decay_budget(), || Ok(())).unwrap();
    let emission = coarse
        .from_decay_step(state.model(), &step, || Ok(()))
        .unwrap();
    // Two physical daughter branches share parent-wide intensities. They are
    // counted once; neither BR nor its square is applied to the photon spectrum.
    near(step.branch_decays[0], 125.0, 1e-8);
    near(step.branch_decays[1], 375.0, 1e-8);
    for (&actual, expected) in emission
        .emitted_particles
        .iter()
        .zip([250.0, 1250.0, 5500.0])
    {
        near(actual, expected, 1e-8);
    }
    assert!(emission.unreported_decays.is_empty());
    assert_eq!(step.released, state.cumulative_released());
}

#[test]
fn daughter_emissions_and_partitioned_histories_use_integrated_counts() {
    let assembly = chain(vec![vec![gamma(U)], vec![gamma(3.0 * U)], vec![]], false);
    let source = prepare(
        Arc::clone(&assembly),
        DecayRadiation::Gamma,
        vec![0.0, 2.0 * U, 4.0 * U],
    );
    let mut whole = inventory(&assembly);
    let step = whole.advance(1.0, decay_budget(), || Ok(())).unwrap();
    let result = source
        .from_decay_step(whole.model(), &step, || Ok(()))
        .unwrap();
    let daughter = 1000.0 * (1.0 - 0.5 * (1.0 + std::f64::consts::LN_2));
    near(result.emitted_particles[0], 500.0, 1e-8);
    near(result.emitted_particles[1], daughter, 1e-8);
    near(result.emitted_energy_j[1], 3.0 * U * daughter, 1e-21);
    let mut split = inventory(&assembly);
    let first = split.advance(0.4, decay_budget(), || Ok(())).unwrap();
    let a = source
        .from_decay_step(split.model(), &first, || Ok(()))
        .unwrap();
    let second = split.advance(0.6, decay_budget(), || Ok(())).unwrap();
    let b = source
        .from_decay_step(split.model(), &second, || Ok(()))
        .unwrap();
    for g in 0..2 {
        near(
            a.emitted_particles[g] + b.emitted_particles[g],
            result.emitted_particles[g],
            1e-8,
        );
        near(
            a.emitted_energy_j[g] + b.emitted_energy_j[g],
            result.emitted_energy_j[g],
            1e-21,
        );
    }
}

#[test]
fn evaluated_cobalt_half_life_emits_published_photon_intensities() {
    let parse = |text, material_number, z| {
        Arc::new(
            DecayEvaluation::parse(
                "pinned ENDF/B-VIII.0 spectral regression".into(),
                EndfDecaySelection {
                    material_number,
                    nuclide: NuclideId::new(z, 60, 0).unwrap(),
                },
                text,
                EndfDecayReadLimits::default(),
                || Ok(()),
            )
            .unwrap(),
        )
    };
    let cobalt = parse(
        include_str!("../../../../../tests/data/endf_b_viii_0/dec-027_Co_060.endf"),
        558,
        27,
    );
    let nickel = parse(
        include_str!("../../../../../tests/data/endf_b_viii_0/dec-028_Ni_060.endf"),
        590,
        28,
    );
    let half_life = cobalt.law().half_life().unwrap().value();
    let assembly = assemble(vec![
        DecayModelEntry {
            evaluation: cobalt,
            routing: routing(),
        },
        DecayModelEntry {
            evaluation: nickel,
            routing: routing(),
        },
    ]);
    let ev = crate::endf::EV_J;
    let source = prepare(
        Arc::clone(&assembly),
        DecayRadiation::Gamma,
        vec![0.0, 1.2e6 * ev, 1.4e6 * ev, 3.0e6 * ev],
    );
    let mut state = inventory(&assembly);
    let step = state.advance(half_life, decay_budget(), || Ok(())).unwrap();
    let photons = source
        .from_decay_step(state.model(), &step, || Ok(()))
        .unwrap();
    near(step.decays[0], 500.0, 1e-8);
    // The small additional gamma lines remain included rather than discarded.
    assert!(photons.emitted_particles[0] >= 500.0 * 0.9985);
    assert!(photons.emitted_particles[0] < 500.0 * 0.999);
    near(photons.emitted_particles[1], 500.0 * 0.999826, 1e-8);
    near(
        photons.emitted_energy_j[1],
        500.0 * 0.999826 * 1.332492e6 * ev,
        1e-21,
    );
    assert_eq!(
        source.coverage(),
        &[
            DecayEmissionCoverage::Reported,
            DecayEmissionCoverage::Stable
        ]
    );
}

#[test]
fn decay_reordered_burnup_competition_and_compartment_sources_agree() {
    let assembly = chain(vec![vec![gamma(2.0 * U)], vec![]], true);
    let source = prepare(
        Arc::clone(&assembly),
        DecayRadiation::Gamma,
        vec![0.0, 4.0 * U],
    );
    let model = burnup(&assembly, false);
    let bound = source.bind_burnup(Arc::clone(&model), || Ok(())).unwrap();
    let flux = GroupFlux::new(Arc::clone(model.groups()), vec![0.0]).unwrap();
    let state = BurnupInventory::new(Arc::clone(&model), vec![1000.0, 0.0]).unwrap();
    let (_, step) = state
        .trial(1.0, &flux, 300.0, burnup_budget(), || Ok(()))
        .unwrap();
    let result = bound.from_burnup_step(&model, &step, || Ok(())).unwrap();
    let mut decay_state = inventory(&assembly);
    let decay_step = decay_state.advance(1.0, decay_budget(), || Ok(())).unwrap();
    let decay = source
        .from_decay_step(decay_state.model(), &decay_step, || Ok(()))
        .unwrap();
    near(
        result.emitted_particles[0],
        decay.emitted_particles[0],
        1e-8,
    );
    near(result.emitted_energy_j[0], decay.emitted_energy_j[0], 1e-21);
    let capture_model = burnup(&assembly, true);
    let capture_source = source
        .bind_burnup(Arc::clone(&capture_model), || Ok(()))
        .unwrap();
    let capture_state =
        BurnupInventory::new(Arc::clone(&capture_model), vec![1000.0, 0.0]).unwrap();
    let flux = GroupFlux::new(Arc::clone(capture_model.groups()), vec![1.0]).unwrap();
    let (_, capture_step) = capture_state
        .trial(1.0, &flux, 300.0, burnup_budget(), || Ok(()))
        .unwrap();
    let photons = capture_source
        .from_burnup_step(&capture_model, &capture_step, || Ok(()))
        .unwrap();
    near(capture_step.captures, 375.0, 1e-8);
    near(photons.emitted_particles[0], 375.0, 1e-8);

    let limits = BurnupTransportLimits::default();
    let transported = BurnupTransportState::new(
        vec![
            state,
            BurnupInventory::new(Arc::clone(&model), vec![0.0, 0.0]).unwrap(),
        ],
        vec![1.0, 2.0],
        NuclideTransport {
            transfers: vec![NuclideTransfer {
                from_compartment: 0,
                to_compartment: 1,
                rates_per_s: vec![std::f64::consts::LN_2, 0.0],
            }],
            outflows: vec![],
        },
        limits,
        || Ok(()),
    )
    .unwrap();
    let zero = GroupFlux::new(Arc::clone(model.groups()), vec![0.0]).unwrap();
    let (_, moved) = transported
        .trial(
            1.0,
            &[zero.clone(), zero],
            &[300.0, 300.0],
            BurnupTransportBudget {
                limits,
                nuclei_tolerance: 1e-6,
                energy_tolerance_j: 1e-15,
            },
            || Ok(()),
        )
        .unwrap();
    let local = bound
        .from_transport_step(transported.model(), &moved, || Ok(()))
        .unwrap();
    near(local[0].emitted_particles[0], 375.0, 1e-8);
    near(local[1].emitted_particles[0], 125.0, 1e-8);
    near(
        local.iter().map(|r| r.emitted_energy_j[0]).sum(),
        decay.emitted_energy_j[0],
        1e-21,
    );
}

#[test]
fn retained_alpha_binding_checks_exact_products_and_counts_emission_energy_once() {
    let mut spectrum = gamma(2.0 * U);
    spectrum.discrete = vec![DiscreteDecayEmission {
        decay_mode: 4.0,
        ..line(2.0 * U, 1.0)
    }]
    .into();
    let template = chain(vec![vec![spectrum], vec![]], false);
    let mut parent = (*template.evaluations()[0]).clone();
    parent.metadata.nuclide = NuclideId::new(4, 8, 0).unwrap();
    parent.metadata.atomic_weight_ratio = 8.0;
    parent.modes = vec![DecayMode {
        rtyp: 4.0,
        q_j: 12.0 * U,
        ..parent.modes()[0].clone()
    }]
    .into();
    parent.average_energies = [0.0, 2.0 * U, 5.0 * U]
        .map(|energy| DecayEnergy {
            value_j: Some(energy),
            uncertainty_j: 0.0,
        })
        .into();
    let mut helium = (*template.evaluations()[1]).clone();
    helium.metadata.nuclide = NuclideId::new(2, 4, 0).unwrap();
    helium.metadata.atomic_weight_ratio = 4.0;
    let original = assemble(vec![
        DecayModelEntry {
            evaluation: Arc::new(parent),
            routing: routing(),
        },
        DecayModelEntry {
            evaluation: Arc::new(helium),
            routing: routing(),
        },
    ]);
    let assembly = Arc::new(
        (*original)
            .clone()
            .with_burnup_product_retention(
                DecayProductRetention {
                    alpha_particles: true,
                    protons: false,
                },
                || Ok(()),
            )
            .unwrap(),
    );
    let source = prepare(
        Arc::clone(&assembly),
        DecayRadiation::Gamma,
        vec![0.0, 4.0 * U],
    );
    let data = assembly
        .burnup_data(group_grid(), 200.0, 1000.0, || Ok(()))
        .unwrap();
    assert_eq!(data.reactions[0].products.len(), 1);
    assert_eq!(data.reactions[0].products[0].yield_per_event, 2.0);
    for change in 0..4 {
        let mut altered = data.clone();
        match change {
            0 => altered.reactions[0].products[0].yield_per_event = 1.0,
            1 => altered.reactions[0].products.clear(),
            2 => altered.reactions[0].products[0].nuclide = 0,
            _ => altered.reactions[0].energy_per_event.escaped_j *= 2.0,
        }
        let wrong = Arc::new(BurnupModel::new(altered, BurnupLimits::default()).unwrap());
        assert!(source.bind_burnup(wrong, || Ok(())).is_err());
    }
    let model = Arc::new(BurnupModel::new(data, BurnupLimits::default()).unwrap());
    let residual_only = prepare(original, DecayRadiation::Gamma, vec![0.0, 4.0 * U]);
    assert!(
        residual_only
            .bind_burnup(Arc::clone(&model), || Ok(()))
            .is_err()
    );
    let bound = source.bind_burnup(Arc::clone(&model), || Ok(())).unwrap();
    let flux = GroupFlux::new(Arc::clone(model.groups()), vec![0.0]).unwrap();
    let state = BurnupInventory::new(Arc::clone(&model), vec![1000.0, 0.0]).unwrap();
    let (after, step) = state
        .trial(1.0, &flux, 300.0, burnup_budget(), || Ok(()))
        .unwrap();
    let photons = bound.from_burnup_step(&model, &step, || Ok(())).unwrap();
    near(after.populations()[0], 500.0, 1e-8);
    near(after.populations()[1], 1000.0, 1e-8);
    near(step.decays, 500.0, 1e-8);
    near(photons.emitted_particles[0], 500.0, 1e-8);
    near(photons.emitted_energy_j[0], 500.0 * 2.0 * U, 1e-21);
    near(step.released.radiation_j, 500.0 * 2.0 * U, 1e-21);
    near(step.released.escaped_j, 500.0 * 5.0 * U, 1e-21);
    near(step.generated_energy_j, 500.0 * 7.0 * U, 1e-21);
    assert_eq!(step.released.deposited_j, 0.0);
}

#[test]
fn missing_spectra_are_reported_and_endpoint_records_are_not_monoenergetic() {
    let assembly = chain(vec![vec![gamma(U)], vec![], vec![]], false);
    let source = prepare(
        Arc::clone(&assembly),
        DecayRadiation::Gamma,
        vec![0.0, 2.0 * U],
    );
    let mut state = inventory(&assembly);
    let step = state.advance(1.0, decay_budget(), || Ok(())).unwrap();
    let emission = source
        .from_decay_step(state.model(), &step, || Ok(()))
        .unwrap();
    assert_eq!(
        source.coverage(),
        &[
            DecayEmissionCoverage::Reported,
            DecayEmissionCoverage::Unreported,
            DecayEmissionCoverage::Stable
        ]
    );
    assert_eq!(
        emission.unreported_decays,
        vec![UnreportedDecayEmission {
            nuclide: assembly.nuclides()[1],
            decays: step.decays[1]
        }]
    );
    let mut beta = gamma(U);
    beta.radiation = DecayRadiation::BetaMinus;
    let beta_assembly = chain(vec![vec![beta.clone()], vec![]], false);
    assert!(matches!(
        PreparedDecayEmission::new(
            beta_assembly,
            DecayRadiation::BetaMinus,
            vec![0.0, 2.0 * U],
            DecayEmissionLimits::default(),
            || Ok(())
        ),
        Err(EndfError::Unsupported(_))
    ));
    // A supplied beta continuum is usable without inventing a shape from Q.
    beta.discrete = vec![].into();
    beta.continuum_normalization = value(1.0);
    beta.continuum = Some(DecayContinuum {
        decay_mode: 1.0,
        points: [(0.0, 1.0 / U), (U, 0.0)]
            .map(|(energy_j, density_per_j)| DecaySpectrumPoint {
                energy_j,
                density_per_j,
            })
            .into(),
        interpolation: vec![DecayInterpolation {
            end_point: 2,
            law: 1,
        }]
        .into(),
    });
    let continuum = prepare(
        chain(vec![vec![beta], vec![]], false),
        DecayRadiation::BetaMinus,
        vec![0.0, U],
    );
    near(continuum.particles_per_decay(0).unwrap()[0], 1.0, 1e-15);
    near(continuum.energy_per_decay_j(0).unwrap()[0], 0.5 * U, 1e-28);
}

#[test]
fn mismatched_models_channels_domains_limits_and_cancellation_refuse() {
    let assembly = chain(vec![vec![gamma(2.0 * U)], vec![]], false);
    let source = prepare(
        Arc::clone(&assembly),
        DecayRadiation::Gamma,
        vec![0.0, 4.0 * U],
    );
    let mut state = inventory(&assembly);
    let step = state.advance(1.0, decay_budget(), || Ok(())).unwrap();
    let clone_model = Arc::new((**state.model()).clone());
    assert!(
        source
            .from_decay_step(&clone_model, &step, || Ok(()))
            .is_err()
    );
    let mut bad = step.clone();
    bad.decays[1] = 1.0;
    assert!(
        source
            .from_decay_step(state.model(), &bad, || Ok(()))
            .is_err()
    );
    for limits in [
        DecayEmissionLimits {
            max_groups: 0,
            ..DecayEmissionLimits::default()
        },
        DecayEmissionLimits {
            max_group_values: 3,
            ..DecayEmissionLimits::default()
        },
        DecayEmissionLimits {
            max_work: 1,
            ..DecayEmissionLimits::default()
        },
    ] {
        assert!(matches!(
            PreparedDecayEmission::new(
                Arc::clone(&assembly),
                DecayRadiation::Gamma,
                vec![0.0, 4.0 * U],
                limits,
                || Ok(())
            ),
            Err(EndfError::Physics(FissionError::Budget(_)))
        ));
    }
    assert!(
        PreparedDecayEmission::new(
            Arc::clone(&assembly),
            DecayRadiation::Gamma,
            vec![0.0, U],
            DecayEmissionLimits::default(),
            || Ok(())
        )
        .is_err()
    );
    let custom = FissionError::InvalidInput("caller checkpoint sentinel");
    assert!(
        matches!(PreparedDecayEmission::new(Arc::clone(&assembly), DecayRadiation::Gamma, vec![0.0, 4.0 * U], DecayEmissionLimits::default(), || Err(custom.clone())), Err(EndfError::Physics(e)) if e == custom)
    );
    let mut calls = 0;
    let cancelled = source.from_decay_step(state.model(), &step, || {
        calls += 1;
        if calls == 2 {
            Err(custom.clone())
        } else {
            Ok(())
        }
    });
    assert!(matches!(cancelled, Err(EndfError::Physics(e)) if e == custom));
    assert_eq!(state.cumulative_decays(), step.decays);
    let mut data = assembly
        .burnup_data(group_grid(), 200.0, 1000.0, || Ok(()))
        .unwrap();
    data.reactions[0].energy_per_event.radiation_j *= 2.0;
    let wrong_route = Arc::new(BurnupModel::new(data, BurnupLimits::default()).unwrap());
    assert!(source.bind_burnup(wrong_route, || Ok(())).is_err());
    let model = burnup(&assembly, false);
    let bound = source.bind_burnup(Arc::clone(&model), || Ok(())).unwrap();
    assert!(bound.check_model(&Arc::new((*model).clone())).is_err());
}

#[test]
fn positron_continua_do_not_guess_capture_or_positron_normalization() {
    let mut spectrum = gamma(U);
    spectrum.radiation = DecayRadiation::ElectronCapturePositron;
    spectrum.discrete = vec![].into();
    spectrum.continuum_normalization = value(1.0);
    spectrum.continuum = Some(DecayContinuum {
        decay_mode: 2.0,
        points: [(0.0, 1.0 / U), (U, 0.0)]
            .map(|(energy_j, density_per_j)| DecaySpectrumPoint {
                energy_j,
                density_per_j,
            })
            .into(),
        interpolation: vec![DecayInterpolation {
            end_point: 2,
            law: 1,
        }]
        .into(),
    });
    // No discrete endpoint is needed to encounter the unresolved STYP=2
    // positron-intensity normalization of ENDF §8.4.2 rule 5(c).
    let continuum_only = chain(vec![vec![spectrum.clone()], vec![]], false);
    assert!(matches!(
        PreparedDecayEmission::new(
            continuum_only,
            DecayRadiation::ElectronCapturePositron,
            vec![0.0, 2.0 * U],
            DecayEmissionLimits::default(),
            || Ok(())
        ),
        Err(EndfError::Unsupported(_))
    ));
    // A zero FD must not bypass this refusal. A retained RIS on an endpoint
    // does not supply an admitted normalization for the continuum as a whole.
    spectrum.discrete_normalization = value(0.0);
    let mut endpoint = line(U, 1.0);
    endpoint.pair_or_positron = Some(value(0.4));
    spectrum.discrete = vec![endpoint].into();
    let retained_endpoint = chain(vec![vec![spectrum], vec![]], false);
    assert!(matches!(
        PreparedDecayEmission::new(
            retained_endpoint,
            DecayRadiation::ElectronCapturePositron,
            vec![0.0, 2.0 * U],
            DecayEmissionLimits::default(),
            || Ok(())
        ),
        Err(EndfError::Unsupported(_))
    ));
}

#[test]
fn subnormal_per_decay_rounding_is_not_amplified_by_large_inventories() {
    let mut spectrum = gamma(1.0);
    spectrum.discrete_normalization = value(1.6e-162);
    spectrum.discrete = vec![line(1.0, 2.0e-162)].into();
    let source = prepare(
        chain(vec![vec![spectrum], vec![]], false),
        DecayRadiation::Gamma,
        vec![0.0, 2.0],
    );
    let rounded = source.particles_per_decay(0).unwrap()[0];
    assert!(rounded.is_subnormal());
    // The exact coefficient is approximately 3.2e-324. Its nearest f64 is
    // 4.94e-324: multiplying that rounded value by 1e300 creates a normal-range
    // answer with more than 50% relative error, even without output overflow.
    assert!(matches!(
        source.apply(1.0, &[1e300, 0.0], &mut 0, &mut || Ok(())),
        Err(EndfError::Physics(FissionError::Arithmetic))
    ));
    let single = source
        .apply(1.0, &[1.0, 0.0], &mut 0, &mut || Ok(()))
        .unwrap();
    assert_eq!(single.emitted_particles, vec![rounded]);
    assert_eq!(single.emitted_energy_j, vec![rounded]);
    let idle = source
        .apply(0.0, &[0.0, 0.0], &mut 0, &mut || Ok(()))
        .unwrap();
    assert_eq!(idle.emitted_particles, vec![0.0]);
}

#[test]
fn unsupported_positive_continuum_is_not_silently_linearized() {
    let mut spectrum = gamma(U);
    spectrum.discrete = vec![].into();
    spectrum.continuum_normalization = value(1.0);
    spectrum.continuum = Some(DecayContinuum {
        decay_mode: 1.0,
        points: [(U, 1.0 / U), (2.0 * U, 2.0 / U)]
            .map(|(energy_j, density_per_j)| DecaySpectrumPoint {
                energy_j,
                density_per_j,
            })
            .into(),
        interpolation: vec![DecayInterpolation {
            end_point: 2,
            law: 6,
        }]
        .into(),
    });
    let assembly = chain(vec![vec![spectrum], vec![]], false);
    assert!(matches!(
        PreparedDecayEmission::new(
            assembly,
            DecayRadiation::Gamma,
            vec![0.0, 3.0 * U],
            DecayEmissionLimits::default(),
            || Ok(())
        ),
        Err(EndfError::Unsupported(_))
    ));
}

fn continuum_spectrum(law: u8, points: &[(f64, f64)], normalization: f64) -> DecaySpectrum {
    let mut spectrum = gamma(U);
    spectrum.discrete = vec![].into();
    spectrum.continuum_normalization = value(normalization);
    spectrum.continuum = Some(DecayContinuum {
        decay_mode: 1.0,
        points: points
            .iter()
            .map(|&(energy_j, density_per_j)| DecaySpectrumPoint {
                energy_j,
                density_per_j,
            })
            .collect::<Vec<_>>()
            .into(),
        interpolation: vec![DecayInterpolation {
            end_point: points.len(),
            law,
        }]
        .into(),
    });
    spectrum
}

fn prepare_continuum(spectrum: DecaySpectrum, edges: Vec<f64>) -> PreparedDecayEmission {
    prepare(
        chain(vec![vec![spectrum], vec![]], false),
        DecayRadiation::Gamma,
        edges,
    )
}

fn relative(actual: f64, expected: f64) {
    near(actual, expected, 3e-12 * expected.abs());
}

// G1: manufactured functions use the exact INT definitions in ENDF-102 (2023)
// §0.5.2.1. Native libm/closed antiderivatives form an independent oracle for
// the strict-arithmetic production operation tree.
#[test]
fn logarithmic_laws_retain_analytic_particle_and_energy_integrals() {
    let e2 = 2.0_f64.exp();
    let logarithmic = prepare_continuum(
        continuum_spectrum(3, &[(1.0, 1.0), (e2, 3.0)], 2.5),
        vec![1.0, 1.0_f64.exp(), e2],
    );
    let edges = logarithmic.group_edges_j();
    for (group, interval) in edges.windows(2).enumerate() {
        let [a, b] = [interval[0], interval[1]];
        let count = b * b.ln() - a * a.ln();
        let moment = |x: f64| x * x * (0.5 * x.ln() + 0.25);
        relative(
            logarithmic.particles_per_decay(0).unwrap()[group],
            2.5 * count,
        );
        relative(
            logarithmic.energy_per_decay_j(0).unwrap()[group],
            2.5 * (moment(b) - moment(a)),
        );
    }
    for sign in [-1.0, 1.0] {
        let exponential = prepare_continuum(
            continuum_spectrum(4, &[(0.0, 1.0), (2.0, (2.0_f64 * sign).exp())], 1.0),
            vec![0.0, 0.5, 1.25, 2.0],
        );
        for (group, interval) in exponential.group_edges_j().windows(2).enumerate() {
            let [a, b] = [interval[0], interval[1]];
            let count = ((sign * b).exp() - (sign * a).exp()) / sign;
            let moment = |x: f64| (sign * x - 1.0) * (sign * x).exp();
            relative(exponential.particles_per_decay(0).unwrap()[group], count);
            relative(
                exponential.energy_per_decay_j(0).unwrap()[group],
                moment(b) - moment(a),
            );
        }
    }
    for power in [
        -3.0,
        -2.0 - 1e-12,
        -2.0,
        -2.0 + 1e-12,
        -1.0 - 1e-12,
        -1.0,
        -1.0 + 1e-12,
        0.0,
        2.5,
    ] {
        let power_law = prepare_continuum(
            continuum_spectrum(5, &[(1.0, 1.0), (4.0, 4.0_f64.powf(power))], 1.0),
            vec![1.0, 1.125, 2.0, 3.5, 4.0],
        );
        for (group, interval) in power_law.group_edges_j().windows(2).enumerate() {
            let [a, b] = [interval[0], interval[1]];
            let moment = |order: f64| {
                let k = power + order + 1.0;
                if k == 0.0 {
                    (b / a).ln()
                } else {
                    a.powf(k) * (k * (b / a).ln()).exp_m1() / k
                }
            };
            relative(
                power_law.particles_per_decay(0).unwrap()[group],
                moment(0.0),
            );
            relative(power_law.energy_per_decay_j(0).unwrap()[group], moment(1.0));
        }
    }
}

#[test]
fn logarithmic_spectra_preserve_group_refinement_discontinuities_and_unit_scaling() {
    for law in 3..=5 {
        let coarse = prepare_continuum(
            continuum_spectrum(law, &[(1.0, 0.5), (4.0, 8.0)], 0.75),
            vec![1.0, 4.0],
        );
        for scale in [1e-180, U, 1.0, 1e180] {
            let refined = prepare_continuum(
                continuum_spectrum(
                    law,
                    &[(scale, 0.5 / scale), (4.0 * scale, 8.0 / scale)],
                    0.75,
                ),
                [1.0, 1.5, 2.25, 3.0, 4.0].map(|e| e * scale).to_vec(),
            );
            relative(
                refined.particles_per_decay(0).unwrap().iter().sum(),
                coarse.particles_per_decay(0).unwrap()[0],
            );
            relative(
                refined.energy_per_decay_j(0).unwrap().iter().sum::<f64>() / scale,
                coarse.energy_per_decay_j(0).unwrap()[0],
            );
        }
    }
    let mut mixed = continuum_spectrum(3, &[(1.0, 1.0), (2.0, 3.0), (2.0, 5.0), (4.0, 10.0)], 1.0);
    mixed.continuum.as_mut().unwrap().interpolation = vec![
        DecayInterpolation {
            end_point: 2,
            law: 3,
        },
        DecayInterpolation {
            end_point: 4,
            law: 5,
        },
    ]
    .into();
    let source = prepare_continuum(mixed, vec![1.0, 2.0, 3.0, 4.0]);
    relative(
        source.particles_per_decay(0).unwrap().iter().sum(),
        20.0 - 2.0 / 2.0_f64.ln(),
    );
    relative(
        source.energy_per_decay_j(0).unwrap().iter().sum(),
        5.5 - 1.5 / 2.0_f64.ln() + 140.0 / 3.0,
    );
}

#[test]
fn logarithmic_moments_keep_flat_and_nearly_coincident_limits() {
    for law in 3..=5 {
        for right in [1.0_f64.next_up(), 1.0 + 1e-10] {
            let source = prepare_continuum(
                continuum_spectrum(law, &[(1.0, 2.0), (right, 2.0)], 1.0),
                vec![1.0, right],
            );
            relative(
                source.particles_per_decay(0).unwrap()[0] / (right - 1.0),
                2.0,
            );
            relative(
                source.energy_per_decay_j(0).unwrap()[0] / (right - 1.0),
                1.0 + right,
            );
        }
        let delta = (1.0 + 1e-12) - 1.0;
        let source = prepare_continuum(
            continuum_spectrum(law, &[(1.0, 1.0), (4.0, 1.0 + delta)], 1.0),
            vec![1.0, 4.0],
        );
        let (count, energy) = if law == 4 {
            (3.0 + 1.5 * delta, 7.5 + 4.5 * delta)
        } else {
            (
                3.0 + delta * (4.0 - 3.0 / 4.0_f64.ln()),
                7.5 + delta * (8.0 - 3.75 / 4.0_f64.ln()),
            )
        };
        near(
            source.particles_per_decay(0).unwrap()[0],
            count,
            4e-15 * count,
        );
        near(
            source.energy_per_decay_j(0).unwrap()[0],
            energy,
            4e-15 * energy,
        );
    }
}

#[test]
fn logarithmic_moments_preserve_representable_extreme_scale_products() {
    // Neither E² nor the full endpoint density ratio is representable, but
    // the physical zeroth and first moments both are. No saturation is needed.
    let source = prepare_continuum(
        continuum_spectrum(5, &[(1e-200, 1e200), (1e200, 1e-200)], 1.0),
        vec![1e-200, 1.0, 1e200],
    );
    for &count in source.particles_per_decay(0).unwrap() {
        relative(count, 200.0 * 10.0_f64.ln());
    }
    relative(source.energy_per_decay_j(0).unwrap()[0], 1.0);
    relative(source.energy_per_decay_j(0).unwrap()[1], 1e200);
    let source = prepare_continuum(
        continuum_spectrum(4, &[(0.0, 1e300), (1.0, 1e-300)], 1e-300),
        vec![0.0, 0.5, 1.0],
    );
    let contrast = 600.0 * 10.0_f64.ln();
    relative(source.particles_per_decay(0).unwrap()[0], 1.0 / contrast);
    relative(
        source.energy_per_decay_j(0).unwrap()[0],
        1.0 / (contrast * contrast),
    );
    relative(source.particles_per_decay(0).unwrap()[1], 1e-300 / contrast);
    relative(
        source.energy_per_decay_j(0).unwrap()[1],
        1e-300 * (0.5 + 1.0 / contrast) / contrast,
    );
}

#[test]
fn logarithmic_domain_rules_reject_zero_logarithms_without_guessing() {
    for (law, points) in [
        (3, [(0.0, 1.0), (1.0, 2.0)]),
        (3, [(0.0, 0.0), (1.0, 0.0)]),
        (4, [(0.0, 0.0), (1.0, 1.0)]),
        (4, [(0.0, 0.0), (1.0, 0.0)]),
        (5, [(0.0, 1.0), (1.0, 2.0)]),
        (5, [(1.0, 1.0), (2.0, 0.0)]),
    ] {
        let assembly = chain(
            vec![vec![continuum_spectrum(law, &points, 1.0)], vec![]],
            false,
        );
        assert!(matches!(
            PreparedDecayEmission::new(
                assembly,
                DecayRadiation::Gamma,
                vec![0.0, 3.0],
                DecayEmissionLimits::default(),
                || Ok(())
            ),
            Err(EndfError::Physics(FissionError::InvalidInput(_)))
        ));
    }
    // INT=3 logs energy only: a true zero density endpoint remains valid.
    let permitted = prepare_continuum(
        continuum_spectrum(3, &[(1.0, 0.0), (2.0, 1.0)], 1.0),
        vec![1.0, 2.0],
    );
    relative(
        permitted.particles_per_decay(0).unwrap()[0],
        2.0 - 1.0 / 2.0_f64.ln(),
    );
    // A derived subnormal density would be rounded before a large supplied
    // normalization could restore its scale; refuse instead of amplifying it.
    let assembly = chain(
        vec![
            vec![continuum_spectrum(
                5,
                &[(1.0, 1e-320), (4.0, 2e-320)],
                1e300,
            )],
            vec![],
        ],
        false,
    );
    assert!(matches!(
        PreparedDecayEmission::new(
            assembly,
            DecayRadiation::Gamma,
            vec![1.0, 2.0, 4.0],
            DecayEmissionLimits::default(),
            || Ok(())
        ),
        Err(EndfError::Physics(FissionError::Arithmetic))
    ));
}
