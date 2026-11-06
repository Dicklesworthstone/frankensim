//! G1 parent/daughter energy accounts and G3 original-owner agreement.

use super::*;
use crate::burnup::{BurnupBudget, BurnupInventory, BurnupLimits, BurnupReaction};
use crate::decay::stochastic::{StochasticDecayBudget, StochasticDecayInventory, StreamKey};
use crate::decay::{DecayBudget, DecayInventory};
use crate::endf::decay::model::{
    DecayEnergyFractions, DecayEnergyRouting, DecayModelEntry, EvaluatedDecayAssembly,
};
use crate::endf::decay::{DecayEvaluation, EndfDecayReadLimits, EndfDecaySelection};
use crate::irradiation::{GroupFlux, NeutronGroups, NuclideId};
use fs_la::decay::DecayLimits;

const EV_J: f64 = 1.602_176_634e-19;

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:.17e} != {expected:.17e}"
    );
}

fn real(value: f64) -> String {
    for precision in (0..=6).rev() {
        let text = format!("{value:.precision$e}");
        let (mantissa, exponent) = text.split_once('e').unwrap();
        let exponent: i32 = exponent.parse().unwrap();
        let field = format!("{mantissa}{exponent:+}");
        if field.len() <= 11 {
            return format!("{field:>11}");
        }
    }
    panic!("manufactured ENDF field does not fit");
}

fn record(fields: [String; 6], mt: i32) -> String {
    format!("{}{:4}{:2}{mt:3}{:5}", fields.concat(), 777, 8, 1)
}

fn values(values: [f64; 6]) -> String {
    record(values.map(real), 457)
}

fn control(c1: f64, c2: f64, fields: [i32; 4], mt: i32) -> String {
    record(
        [
            real(c1),
            real(c2),
            format!("{:11}", fields[0]),
            format!("{:11}", fields[1]),
            format!("{:11}", fields[2]),
            format!("{:11}", fields[3]),
        ],
        mt,
    )
}

pub(crate) fn evaluation(
    z: u16,
    stable: bool,
    averages_ev: [f64; 3],
    line: Option<(f64, f64)>,
) -> Arc<DecayEvaluation> {
    let id = NuclideId::new(z, 30, 0).unwrap();
    let mut text = vec![
        control(
            f64::from(z) * 1000.0 + 30.0,
            30.0,
            [0, 0, i32::from(stable), i32::from(line.is_some())],
            457,
        ),
        control(if stable { 0.0 } else { 1.0 }, 0.0, [0, 0, 6, 0], 457),
        values([
            averages_ev[0],
            0.0,
            averages_ev[1],
            0.0,
            averages_ev[2],
            0.0,
        ]),
        control(0.0, 1.0, [0, 0, 6, i32::from(!stable)], 457),
        values(if stable {
            [0.0; 6]
        } else {
            [
                1.0,
                0.0,
                averages_ev.iter().sum::<f64>() + 1000.0,
                0.0,
                1.0,
                0.0,
            ]
        }),
    ];
    if let Some((energy_ev, intensity)) = line {
        text.extend([
            control(0.0, 0.0, [0, 0, 6, 1], 457),
            values([1.0, 0.0, energy_ev * intensity, 0.0, 0.0, 0.0]),
            control(energy_ev, 0.0, [0, 0, 6, 0], 457),
            values([1.0, 0.0, intensity, 0.0, 0.0, 0.0]),
        ]);
    }
    text.push(control(0.0, 0.0, [0; 4], 0));
    Arc::new(
        DecayEvaluation::parse(
            format!("manufactured-radiation-Z{z}"),
            EndfDecaySelection {
                material_number: 777,
                nuclide: id,
            },
            &text.join("\n"),
            EndfDecayReadLimits::default(),
            || Ok(()),
        )
        .unwrap(),
    )
}

pub(crate) fn routing(outward: f64) -> DecayEnergyRouting {
    let local = DecayEnergyFractions {
        deposited: 1.0,
        radiation: 0.0,
        escaped: 0.0,
    };
    DecayEnergyRouting {
        light_particles: local,
        heavy_particles: local,
        electromagnetic: DecayEnergyFractions {
            deposited: 1.0 - outward,
            radiation: outward,
            escaped: 0.0,
        },
        include_reported_neutrinos: false,
    }
}

fn assembly(
    parent: Arc<DecayEvaluation>,
    parent_routing: DecayEnergyRouting,
) -> Arc<EvaluatedDecayAssembly> {
    Arc::new(
        EvaluatedDecayAssembly::new(
            "manufactured radioactive radiation chain".into(),
            vec![
                DecayModelEntry {
                    evaluation: parent,
                    routing: parent_routing,
                },
                DecayModelEntry {
                    evaluation: evaluation(11, false, [0.0, 200.0, 0.0], Some((200.0, 1.0))),
                    routing: routing(0.25),
                },
                DecayModelEntry {
                    evaluation: evaluation(12, true, [0.0; 3], None),
                    routing: routing(0.0),
                },
            ],
            DecayLimits::default(),
            || Ok(()),
        )
        .unwrap(),
    )
}

/// Shared manufactured public-parser data only; consumers supply independent
/// expected outcomes. Parent/daughter half-lives are one second, with a stable
/// endpoint. The order is Z10, Z11, Z12 at fixed mass number 30.
pub(crate) fn gamma_chain(
    parent_ev: f64,
    daughter_ev: f64,
    parent_outward: f64,
    daughter_outward: f64,
) -> Arc<EvaluatedDecayAssembly> {
    Arc::new(
        EvaluatedDecayAssembly::new(
            "manufactured shared gamma chain".into(),
            vec![
                DecayModelEntry {
                    evaluation: evaluation(
                        10,
                        false,
                        [0.0, parent_ev, 0.0],
                        Some((parent_ev, 1.0)),
                    ),
                    routing: routing(parent_outward),
                },
                DecayModelEntry {
                    evaluation: evaluation(
                        11,
                        false,
                        [0.0, daughter_ev, 0.0],
                        Some((daughter_ev, 1.0)),
                    ),
                    routing: routing(daughter_outward),
                },
                DecayModelEntry {
                    evaluation: evaluation(12, true, [0.0; 3], None),
                    routing: routing(0.0),
                },
            ],
            DecayLimits::default(),
            || Ok(()),
        )
        .unwrap(),
    )
}

fn ordinary() -> Arc<EvaluatedDecayAssembly> {
    gamma_chain(100.0, 200.0, 0.5, 0.25)
}

fn source(
    assembly: Arc<EvaluatedDecayAssembly>,
    radiation: DecayRadiation,
) -> Arc<PreparedDecayEmission> {
    Arc::new(
        PreparedDecayEmission::new(
            assembly,
            radiation,
            vec![0.0, 150.0 * EV_J, 300.0 * EV_J],
            DecayEmissionLimits::default(),
            || Ok(()),
        )
        .unwrap(),
    )
}

fn prepare(assembly: Arc<EvaluatedDecayAssembly>) -> PreparedDecayRadiation {
    PreparedDecayRadiation::new(
        source(assembly, DecayRadiation::Gamma),
        DecayRadiationClosure {
            absolute_j_per_decay: 0.0,
            relative: 1e-12,
        },
        DecayEmissionLimits::default(),
        || Ok(()),
    )
    .unwrap()
}

fn budget() -> DecayBudget {
    DecayBudget {
        propagation: DecayLimits::default(),
        nuclei_tolerance: 1e-7,
        energy_tolerance_j: 1e-22,
    }
}

fn inventory(assembly: &EvaluatedDecayAssembly) -> DecayInventory {
    DecayInventory::new(Arc::clone(assembly.model()), vec![1000.0, 0.0, 0.0]).unwrap()
}

#[test]
fn g1_per_parent_routing_and_daughter_ingrowth_preserve_original_radiation() {
    let assembly = ordinary();
    let source = prepare(Arc::clone(&assembly));
    assert!(Arc::ptr_eq(source.model(), assembly.model()));
    assert_eq!(source.radiation(), DecayRadiation::Gamma);
    near(source.energy_per_decay_j(0).unwrap()[0], 50.0 * EV_J, 1e-31);
    near(source.energy_per_decay_j(1).unwrap()[1], 50.0 * EV_J, 1e-31);
    assert_eq!(source.energy_per_decay_j(2).unwrap(), &[0.0, 0.0]);
    assert!(source.energy_per_decay_j(3).is_none());
    let state = inventory(&assembly);
    let (_, step) = state.trial(1.0, budget(), || Ok(())).unwrap();
    let result = source
        .from_decay_step(state.model(), &step, || Ok(()))
        .unwrap();
    let daughter_decays = 1000.0 * (1.0 - 0.5 * (1.0 + std::f64::consts::LN_2));
    near(result.emitted_energy_j[0], 500.0 * 50.0 * EV_J, 1e-25);
    near(
        result.emitted_energy_j[1],
        daughter_decays * 50.0 * EV_J,
        1e-25,
    );
    near(
        result.emitted_energy_j.iter().sum(),
        step.released.radiation_j,
        1e-25,
    );
    near(
        result.original_radiation_j,
        step.released.radiation_j,
        1e-25,
    );
    assert_eq!(result.reconciled_spectral_energy_j, 0.0);
    let (middle, first) = state.trial(0.5, budget(), || Ok(())).unwrap();
    let (_, second) = middle.trial(0.5, budget(), || Ok(())).unwrap();
    let first = source
        .from_decay_step(state.model(), &first, || Ok(()))
        .unwrap();
    let second = source
        .from_decay_step(state.model(), &second, || Ok(()))
        .unwrap();
    for group in 0..2 {
        near(
            first.emitted_energy_j[group] + second.emitted_energy_j[group],
            result.emitted_energy_j[group],
            1e-25,
        );
    }
}

#[test]
fn g1_explicit_mean_discrepancy_is_reported_and_does_not_create_heat() {
    let assembly = assembly(
        evaluation(10, false, [0.0, 100.0001, 0.0], Some((100.0, 1.0))),
        routing(0.5),
    );
    let raw = source(Arc::clone(&assembly), DecayRadiation::Gamma);
    assert!(
        PreparedDecayRadiation::new(
            Arc::clone(&raw),
            DecayRadiationClosure {
                absolute_j_per_decay: 0.0,
                relative: 0.0
            },
            DecayEmissionLimits::default(),
            || Ok(())
        )
        .is_err()
    );
    let source = PreparedDecayRadiation::new(
        raw,
        DecayRadiationClosure {
            absolute_j_per_decay: 0.0,
            relative: 2e-6,
        },
        DecayEmissionLimits::default(),
        || Ok(()),
    )
    .unwrap();
    near(source.discrepancy_per_decay_j()[0], 0.00005 * EV_J, 1e-30);
    near(
        source.energy_per_decay_j(0).unwrap().iter().sum(),
        50.00005 * EV_J,
        1e-31,
    );
    let state = inventory(&assembly);
    let (_, step) = state.trial(1.0, budget(), || Ok(())).unwrap();
    let result = source
        .from_decay_step(state.model(), &step, || Ok(()))
        .unwrap();
    near(
        result.emitted_energy_j.iter().sum(),
        step.released.radiation_j,
        1e-25,
    );
    near(
        result.reconciled_spectral_energy_j,
        500.0 * source.discrepancy_per_decay_j()[0],
        1e-30,
    );
}

#[test]
fn g0_missing_spectra_and_other_outward_categories_are_not_normalized_away() {
    let loose = DecayRadiationClosure {
        absolute_j_per_decay: 1.0,
        relative: 1.0,
    };
    for line in [None, Some((100.0, 0.0))] {
        let assembly = assembly(evaluation(10, false, [0.0, 100.0, 0.0], line), routing(0.5));
        assert!(
            PreparedDecayRadiation::new(
                source(assembly, DecayRadiation::Gamma),
                loose,
                DecayEmissionLimits::default(),
                || Ok(())
            )
            .is_err()
        );
    }
    let mut routes = routing(0.5);
    routes.light_particles = DecayEnergyFractions {
        deposited: 0.5,
        radiation: 0.5,
        escaped: 0.0,
    };
    let mixed = assembly(
        evaluation(10, false, [10.0, 100.0, 0.0], Some((100.0, 1.0))),
        routes,
    );
    assert!(
        PreparedDecayRadiation::new(
            source(mixed, DecayRadiation::Gamma),
            loose,
            DecayEmissionLimits::default(),
            || Ok(())
        )
        .is_err()
    );
    let assembly = assembly(evaluation(10, false, [0.0, 100.0, 0.0], None), routing(0.0));
    let prepared = prepare(Arc::clone(&assembly));
    assert_eq!(prepared.original_radiation_per_decay_j()[0], 0.0);
    assert_eq!(prepared.energy_per_decay_j(0).unwrap(), &[0.0, 0.0]);
    let state = inventory(&assembly);
    let (_, step) = state.trial(1.0, budget(), || Ok(())).unwrap();
    let output = prepared
        .from_decay_step(state.model(), &step, || Ok(()))
        .unwrap();
    assert_eq!(output.emitted_energy_j[0], 0.0);
    assert!(output.emitted_energy_j[1] > 0.0);
    for radiation in [
        DecayRadiation::Other(99),
        DecayRadiation::FissionFragment,
        DecayRadiation::ElectronNeutrino,
    ] {
        assert!(matches!(
            PreparedDecayRadiation::new(
                source(Arc::clone(&assembly), radiation),
                loose,
                DecayEmissionLimits::default(),
                || Ok(())
            ),
            Err(EndfError::Unsupported(_))
        ));
    }
}

fn burnup(assembly: &EvaluatedDecayAssembly, capture: bool) -> Arc<BurnupModel> {
    let grid =
        Arc::new(NeutronGroups::new("manufactured neutron grid".into(), vec![0.0, 1.0]).unwrap());
    let mut data = assembly
        .burnup_data(grid, 200.0, 1000.0, || Ok(()))
        .unwrap();
    data.reactions.reverse();
    if capture {
        data.reactions.insert(
            0,
            BurnupReaction {
                data_id: "independent prompt capture radiation".into(),
                parent: 0,
                rate: BurnupRate::Capture {
                    cross_sections_m2: vec![std::f64::consts::LN_2],
                },
                products: vec![],
                energy_per_event: crate::EnergyChannels {
                    deposited_j: 0.0,
                    radiation_j: 7.0 * EV_J,
                    escaped_j: 0.0,
                },
            },
        );
    }
    Arc::new(BurnupModel::new(data, BurnupLimits::default()).unwrap())
}

fn burnup_step(model: Arc<BurnupModel>, flux: f64) -> BurnupStep {
    let flux = GroupFlux::new(Arc::clone(model.groups()), vec![flux]).unwrap();
    BurnupInventory::new(model, vec![1000.0, 0.0, 0.0])
        .unwrap()
        .trial(
            1.0,
            &flux,
            300.0,
            BurnupBudget {
                limits: BurnupLimits::default(),
                nuclei_tolerance: 1e-7,
                energy_tolerance_j: 1e-22,
            },
            || Ok(()),
        )
        .unwrap()
        .1
}

#[test]
fn g3_burnup_uses_original_decay_events_and_excludes_prompt_capture_radiation() {
    let assembly = ordinary();
    let partition = prepare(Arc::clone(&assembly));
    let model = burnup(&assembly, false);
    let bound = partition
        .bind_burnup(Arc::clone(&model), || Ok(()))
        .unwrap();
    assert!(Arc::ptr_eq(bound.model(), &model));
    assert!(Arc::ptr_eq(bound.source().model(), assembly.model()));
    let result = bound
        .from_burnup_step(&model, &burnup_step(Arc::clone(&model), 0.0), || Ok(()))
        .unwrap();
    let state = inventory(&assembly);
    let (_, step) = state.trial(1.0, budget(), || Ok(())).unwrap();
    let pure = partition
        .from_decay_step(state.model(), &step, || Ok(()))
        .unwrap();
    for group in 0..2 {
        near(
            result.emitted_energy_j[group],
            pure.emitted_energy_j[group],
            1e-25,
        );
    }
    let model = burnup(&assembly, true);
    let bound = partition
        .bind_burnup(Arc::clone(&model), || Ok(()))
        .unwrap();
    let step = burnup_step(Arc::clone(&model), 1.0);
    let result = bound.from_burnup_step(&model, &step, || Ok(())).unwrap();
    near(step.captures, 375.0, 1e-8);
    near(result.emitted_energy_j[0], 375.0 * 50.0 * EV_J, 1e-25);
    near(
        result.original_radiation_j + step.captures * 7.0 * EV_J,
        step.released.radiation_j,
        1e-25,
    );
    near(
        result.emitted_energy_j.iter().sum(),
        result.original_radiation_j,
        1e-25,
    );
}

#[test]
fn g0_identity_budget_duration_and_cancellation_refusals_preserve_source() {
    let assembly = ordinary();
    let partition = prepare(Arc::clone(&assembly));
    let state = inventory(&assembly);
    let (_, step) = state.trial(1.0, budget(), || Ok(())).unwrap();
    assert!(
        partition
            .from_decay_step(&Arc::new((**state.model()).clone()), &step, || Ok(()))
            .is_err()
    );
    let model = burnup(&assembly, false);
    let bound = partition
        .bind_burnup(Arc::clone(&model), || Ok(()))
        .unwrap();
    let report = burnup_step(Arc::clone(&model), 0.0);
    assert!(
        bound
            .from_burnup_step(&Arc::new((*model).clone()), &report, || Ok(()))
            .is_err()
    );
    let mut bad = step.clone();
    bad.duration_s = 0.0;
    assert!(
        partition
            .from_decay_step(state.model(), &bad, || Ok(()))
            .is_err()
    );
    let mut bad = step.clone();
    bad.decays[2] = 1.0;
    assert!(
        partition
            .from_decay_step(state.model(), &bad, || Ok(()))
            .is_err()
    );
    for constrained in [
        DecayEmissionLimits {
            max_work: 0,
            ..DecayEmissionLimits::default()
        },
        DecayEmissionLimits {
            max_groups: 1,
            ..DecayEmissionLimits::default()
        },
        DecayEmissionLimits {
            max_group_values: 1,
            ..DecayEmissionLimits::default()
        },
    ] {
        assert!(
            PreparedDecayRadiation::new(
                Arc::clone(partition.source()),
                DecayRadiationClosure {
                    absolute_j_per_decay: 0.0,
                    relative: 1e-12
                },
                constrained,
                || Ok(())
            )
            .is_err()
        );
    }
    for closure in [
        DecayRadiationClosure {
            absolute_j_per_decay: -1.0,
            relative: 0.0,
        },
        DecayRadiationClosure {
            absolute_j_per_decay: 0.0,
            relative: f64::NAN,
        },
    ] {
        assert!(
            PreparedDecayRadiation::new(
                Arc::clone(partition.source()),
                closure,
                DecayEmissionLimits::default(),
                || Ok(())
            )
            .is_err()
        );
    }
    let expected = partition
        .from_decay_step(state.model(), &step, || Ok(()))
        .unwrap();
    let mut visits = 0;
    partition
        .from_decay_step(state.model(), &step, || {
            visits += 1;
            Ok(())
        })
        .unwrap();
    for stop in [1, visits] {
        let mut seen = 0;
        assert_eq!(
            partition.from_decay_step(state.model(), &step, || {
                seen += 1;
                if seen == stop {
                    Err(FissionError::Cancelled)
                } else {
                    Ok(())
                }
            }),
            Err(EndfError::Physics(FissionError::Cancelled))
        );
    }
    assert_eq!(
        partition
            .from_decay_step(state.model(), &step, || Ok(()))
            .unwrap(),
        expected
    );
}

fn stochastic(assembly: &EvaluatedDecayAssembly) -> StochasticDecayInventory {
    StochasticDecayInventory::new(
        Arc::clone(assembly.model()),
        vec![128, 0, 0],
        StreamKey {
            seed: 77,
            kernel: 0xdec0_5ec7,
            tile: 31,
        },
    )
    .unwrap()
}

fn stochastic_budget() -> StochasticDecayBudget {
    StochasticDecayBudget {
        maximum_events: 256,
        energy_tolerance_j: 1e-22,
    }
}

#[test]
fn g3_integer_trajectory_radiation_preserves_replay_partition_and_random_stream() {
    let assembly = ordinary();
    let partition = prepare(Arc::clone(&assembly));
    let state = stochastic(&assembly);
    let before = state.clone();
    let (whole, step) = state.trial_to(2.0, stochastic_budget(), || Ok(())).unwrap();
    let projected = partition
        .from_stochastic_step(state.model(), &step, || Ok(()))
        .unwrap();
    assert!(step.decays[0] > 0 && step.decays[1] > 0);
    for parent in 0..2 {
        near(
            projected.emitted_energy_j[parent],
            step.decays[parent] as f64 * 50.0 * EV_J,
            1e-28,
        );
    }
    near(
        projected.original_radiation_j,
        step.released.radiation_j,
        1e-28,
    );
    near(
        projected.emitted_energy_j.iter().sum(),
        projected.original_radiation_j,
        1e-28,
    );
    assert_eq!(projected.reconciled_spectral_energy_j, 0.0);
    assert_eq!(
        partition
            .from_stochastic_step(state.model(), &step, || Ok(()))
            .unwrap(),
        projected
    );
    assert_eq!(state, before);
    let (middle, first) = state.trial_to(0.7, stochastic_budget(), || Ok(())).unwrap();
    let (split, second) = middle
        .trial_to(2.0, stochastic_budget(), || Ok(()))
        .unwrap();
    let first = partition
        .from_stochastic_step(state.model(), &first, || Ok(()))
        .unwrap();
    let second = partition
        .from_stochastic_step(state.model(), &second, || Ok(()))
        .unwrap();
    assert_eq!(whole, split);
    assert_eq!(whole.rng_checkpoint(), split.rng_checkpoint());
    for group in 0..2 {
        near(
            first.emitted_energy_j[group] + second.emitted_energy_j[group],
            projected.emitted_energy_j[group],
            1e-28,
        );
    }
    near(
        first.original_radiation_j + second.original_radiation_j,
        projected.original_radiation_j,
        1e-28,
    );
}

#[test]
fn g0_integer_radiation_requires_consistent_counts_and_exact_conversion_ceiling() {
    let assembly = ordinary();
    let partition = prepare(Arc::clone(&assembly));
    let state = stochastic(&assembly);
    let (_, step) = state.trial_to(2.0, stochastic_budget(), || Ok(())).unwrap();
    assert!(
        partition
            .from_stochastic_step(&Arc::new((**state.model()).clone()), &step, || Ok(()))
            .is_err()
    );
    let mut inconsistent = Vec::new();
    let mut bad = step.clone();
    bad.events += 1;
    inconsistent.push(bad);
    let mut bad = step.clone();
    bad.branch_decays[0] += 1;
    inconsistent.push(bad);
    let mut bad = step.clone();
    bad.decays[0] += 1;
    inconsistent.push(bad);
    let mut bad = step.clone();
    bad.decays.pop();
    inconsistent.push(bad);
    let mut bad = step.clone();
    bad.branch_decays.pop();
    inconsistent.push(bad);
    let mut bad = step.clone();
    bad.duration_s = 0.0;
    bad.end_time_s = bad.start_time_s;
    inconsistent.push(bad);
    let mut bad = step.clone();
    bad.start_time_s = -1.0;
    inconsistent.push(bad);
    let mut bad = step.clone();
    bad.events = MAX_STOCHASTIC_NUCLEI + 1;
    inconsistent.push(bad);
    let mut bad = step.clone();
    bad.branch_decays[0] = MAX_STOCHASTIC_NUCLEI + 1;
    inconsistent.push(bad);
    for bad in inconsistent {
        assert!(
            partition
                .from_stochastic_step(state.model(), &bad, || Ok(()))
                .is_err()
        );
    }
    // This is a public report-admission boundary, not a simulated 2^53-event run.
    let mut ceiling = step;
    ceiling.events = MAX_STOCHASTIC_NUCLEI;
    ceiling.decays = vec![MAX_STOCHASTIC_NUCLEI, 0, 0];
    ceiling.branch_decays = vec![MAX_STOCHASTIC_NUCLEI, 0];
    let accepted = partition
        .from_stochastic_step(state.model(), &ceiling, || Ok(()))
        .unwrap();
    near(
        accepted.original_radiation_j,
        MAX_STOCHASTIC_NUCLEI as f64 * 50.0 * EV_J,
        1e-16,
    );
    assert_eq!(accepted.emitted_energy_j[1], 0.0);
    assert_eq!(
        partition.from_stochastic_step(state.model(), &ceiling, || Err(FissionError::Cancelled)),
        Err(EndfError::Physics(FissionError::Cancelled))
    );
}
