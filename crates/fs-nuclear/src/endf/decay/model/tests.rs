//! G1 analytic inventories/calorimetry and G3 assembly-order/partition checks.
//! H3/He3 records are source-pinned evaluated data; other fixtures are explicitly
//! manufactured and make no empirical isotope or cross-section accuracy claim.

use super::*;
use crate::burnup::{BurnupBudget, BurnupInventory, BurnupLimits, BurnupModel};
use crate::decay::{DecayBudget, DecayInventory};
use crate::endf::decay::{EndfDecayReadLimits, EndfDecaySelection};
use crate::irradiation::GroupFlux;

const EV_J: f64 = 1.602_176_634e-19;

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:.17e} != {expected:.17e}, allowance {tolerance:.3e}"
    );
}

fn identity(z: u16, a: u16, isomer: u16) -> NuclideId {
    NuclideId::new(z, a, isomer).unwrap()
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

fn entry(evaluation: Arc<DecayEvaluation>) -> DecayModelEntry {
    DecayModelEntry {
        evaluation,
        routing: routing(),
    }
}

fn budget() -> DecayBudget {
    DecayBudget {
        propagation: DecayLimits::default(),
        nuclei_tolerance: 1e-7,
        energy_tolerance_j: 1e-20,
    }
}

fn assembly(entries: Vec<DecayModelEntry>) -> EvaluatedDecayAssembly {
    EvaluatedDecayAssembly::new(
        "explicit-source-assembly-and-routing".into(),
        entries,
        DecayLimits::default(),
        || Ok(()),
    )
    .unwrap()
}

fn real_records() -> (Arc<DecayEvaluation>, Arc<DecayEvaluation>) {
    let hydrogen = DecayEvaluation::parse(
        "ENDF/B-VIII.0 H3; Sandia xsec aaafe23d5f80a5eefbf49421a8e7ff591676c13a".into(),
        EndfDecaySelection {
            material_number: 4,
            nuclide: identity(1, 3, 0),
        },
        include_str!("../../../../tests/data/endf_b_viii_0/dec-001_H_003.endf"),
        EndfDecayReadLimits::default(),
        || Ok(()),
    )
    .unwrap();
    let helium = DecayEvaluation::parse(
        "ENDF/B-VIII.0 He3; Sandia xsec aaafe23d5f80a5eefbf49421a8e7ff591676c13a".into(),
        EndfDecaySelection {
            material_number: 9,
            nuclide: identity(2, 3, 0),
        },
        include_str!("../../../../tests/data/endf_b_viii_0/dec-002_He_003.endf"),
        EndfDecayReadLimits::default(),
        || Ok(()),
    )
    .unwrap();
    (Arc::new(hydrogen), Arc::new(helium))
}

// A short independent ENDF fixture encoder exercises the public reader/adapter
// boundary. All represented coefficients below are manufactured.
fn real(value: f64) -> String {
    for precision in (0..=6).rev() {
        let scientific = format!("{value:.precision$e}");
        let (mantissa, exponent) = scientific.split_once('e').unwrap();
        let exponent: i32 = exponent.parse().unwrap();
        let text = format!("{mantissa}{exponent:+}");
        if text.len() <= 11 {
            return format!("{text:>11}");
        }
    }
    panic!("fixture number too wide");
}

fn row(fields: [String; 6], mt: i32) -> String {
    format!("{}{:4}{:2}{mt:3}{:5}", fields.concat(), 777, 8, 1)
}

fn control(first: f64, second: f64, integers: [i32; 4], mt: i32) -> String {
    row(
        [
            real(first),
            real(second),
            format!("{:11}", integers[0]),
            format!("{:11}", integers[1]),
            format!("{:11}", integers[2]),
            format!("{:11}", integers[3]),
        ],
        mt,
    )
}

fn evaluation(
    id: NuclideId,
    half_life_s: Option<f64>,
    energies_ev: &[Option<f64>],
    modes: &[(f64, u16, f64, f64)],
) -> Arc<DecayEvaluation> {
    let mut lines = vec![control(
        1000.0 * f64::from(id.atomic_number()) + f64::from(id.mass_number()),
        f64::from(id.mass_number()),
        [
            i32::from(id.isomer()),
            i32::from(id.isomer()),
            i32::from(half_life_s.is_none()),
            0,
        ],
        457,
    )];
    let energies: Vec<_> = if half_life_s.is_none() {
        vec![Some(0.0); 3]
    } else {
        energies_ev.to_vec()
    };
    lines.push(control(
        half_life_s.unwrap_or(0.0),
        0.0,
        [0, 0, 2 * energies.len() as i32, 0],
        457,
    ));
    let values: Vec<_> = energies
        .iter()
        .flat_map(|energy| [energy.unwrap_or(-1.0), 0.0])
        .collect();
    for values in values.chunks(6) {
        lines.push(row(
            std::array::from_fn(|i| real(values.get(i).copied().unwrap_or(0.0))),
            457,
        ));
    }
    lines.push(control(
        0.0,
        1.0,
        [
            0,
            0,
            if modes.is_empty() {
                6
            } else {
                6 * modes.len() as i32
            },
            modes.len() as i32,
        ],
        457,
    ));
    if modes.is_empty() {
        lines.push(row(std::array::from_fn(|_| real(0.0)), 457));
    } else {
        for &(rtyp, isomer, probability, q) in modes {
            lines.push(row(
                [rtyp, f64::from(isomer), q, 0.0, probability, 0.0].map(real),
                457,
            ));
        }
    }
    lines.push(control(0.0, 0.0, [0; 4], 0));
    Arc::new(
        DecayEvaluation::parse(
            "manufactured-MF8-coefficients".into(),
            EndfDecaySelection {
                material_number: 777,
                nuclide: id,
            },
            &lines.join("\n"),
            EndfDecayReadLimits::default(),
            || Ok(()),
        )
        .unwrap(),
    )
}

fn stable(id: NuclideId) -> Arc<DecayEvaluation> {
    evaluation(id, None, &[], &[])
}

#[test]
fn evaluated_h3_half_life_and_mean_electron_energy_drive_original_inventory() {
    let (hydrogen, helium) = real_records();
    let half_life = hydrogen.law().half_life().unwrap().value();
    assert_eq!(half_life, 3.887896e8);
    let assembled = assembly(vec![entry(helium), entry(Arc::clone(&hydrogen))]);
    assert_eq!(assembled.index(identity(1, 3, 0)), Some(1));
    assert_eq!(assembled.index(identity(2, 4, 0)), None);
    assert!(Arc::ptr_eq(&assembled.evaluations()[1], &hydrogen));
    let initial = DecayInventory::new(Arc::clone(assembled.model()), vec![0.0, 1000.0]).unwrap();
    let (after, report) = initial.trial(half_life, budget(), || Ok(())).unwrap();
    close(after.populations()[0], 500.0, 1e-9);
    close(after.populations()[1], 500.0, 1e-9);
    close(report.decays[1], 500.0, 1e-9);
    close(report.released.deposited_j, 500.0 * 5690.0 * EV_J, 1e-24);
    assert_eq!(report.released.radiation_j, 0.0);
    assert_eq!(report.released.escaped_j, 0.0);
    assert!(
        hydrogen.modes()[0].q_j > 3.0 * assembled.model().branches()[0].energy_per_decay.total_j()
    );
    assert_eq!(initial.populations(), &[0.0, 1000.0]);
}

#[test]
fn equal_half_life_chain_resolves_original_order_and_partitioned_evolution() {
    let first = evaluation(
        identity(3, 10, 0),
        Some(2.0),
        &[Some(10.0), Some(0.0), Some(0.0)],
        &[(1.0, 0, 1.0, 50.0)],
    );
    let second = evaluation(
        identity(4, 10, 0),
        Some(2.0),
        &[Some(0.0), Some(20.0), Some(0.0)],
        &[(1.0, 0, 1.0, 100.0)],
    );
    let final_node = stable(identity(5, 10, 0));
    let assembled = assembly(vec![entry(final_node), entry(first), entry(second)]);
    let initial =
        DecayInventory::new(Arc::clone(assembled.model()), vec![0.0, 1000.0, 0.0]).unwrap();
    let (whole, report) = initial.trial(2.0, budget(), || Ok(())).unwrap();
    let intermediate = 500.0 * core::f64::consts::LN_2;
    close(whole.populations()[1], 500.0, 1e-9);
    close(whole.populations()[2], intermediate, 1e-9);
    close(whole.populations()[0], 500.0 - intermediate, 1e-9);
    close(report.released.deposited_j, 500.0 * 10.0 * EV_J, 1e-26);
    close(
        report.released.radiation_j,
        (500.0 - intermediate) * 20.0 * EV_J,
        1e-26,
    );
    let (half, _) = initial.trial(1.0, budget(), || Ok(())).unwrap();
    let (split, _) = half.trial(1.0, budget(), || Ok(())).unwrap();
    for (actual, expected) in split.populations().iter().zip(whole.populations()) {
        close(*actual, *expected, 1e-9);
    }
    close(
        split.cumulative_released().total_j(),
        whole.cumulative_released().total_j(),
        1e-26,
    );
}

#[test]
fn decay_mode_sequences_resolve_alpha_isomers_and_delayed_particles() {
    for (parent, rtyp, state, expected) in [
        (identity(92, 238, 0), 4.0, 0, identity(90, 234, 0)),
        (identity(95, 242, 2), 3.0, 1, identity(95, 242, 1)),
        (identity(3, 9, 0), 1.5, 0, identity(4, 8, 0)),
        (identity(7, 16, 0), 1.4, 0, identity(6, 12, 0)),
        (identity(5, 12, 0), 1.44, 0, identity(2, 4, 0)),
        (identity(93, 238, 0), 2.0, 0, identity(92, 238, 0)),
        (identity(16, 31, 0), 7.0, 0, identity(15, 30, 0)),
        (identity(10, 22, 0), 5.0, 0, identity(10, 21, 0)),
    ] {
        assert_eq!(daughter(parent, rtyp, state).unwrap(), expected);
    }
    let parent = identity(95, 242, 2);
    let record = evaluation(
        parent,
        Some(1.0),
        &[Some(0.0), Some(4.0), Some(0.0)],
        &[(3.0, 1, 1.0, 4.0)],
    );
    let assembled = assembly(vec![entry(record), entry(stable(identity(95, 242, 1)))]);
    assert_eq!(assembled.model().branches()[0].daughter, 1);
    for mode in [6.0, 1.6, 10.0, 1.05] {
        assert!(matches!(
            daughter(parent, mode, 0),
            Err(DecayAssemblyError::UnsupportedDecayMode { .. })
        ));
    }
}

#[test]
fn parent_mean_emissions_are_routed_once_and_optional_neutrinos_only_escape() {
    let mut means = vec![Some(0.0); 17];
    means[0] = Some(10.0);
    means[1] = Some(20.0);
    means[2] = Some(30.0);
    means[3] = Some(10.0);
    means[7] = Some(20.0);
    means[11] = Some(30.0);
    means[16] = Some(40.0);
    let parent = evaluation(
        identity(10, 20, 0),
        Some(1.0),
        &means,
        &[(1.0, 0, 0.25, 1000.0), (2.0, 0, 0.75, 100.0)],
    );
    let mut with_neutrinos = entry(parent);
    with_neutrinos.routing.include_reported_neutrinos = true;
    let assembled = assembly(vec![
        with_neutrinos,
        entry(stable(identity(11, 20, 0))),
        entry(stable(identity(9, 20, 0))),
    ]);
    assert_eq!(
        assembled.model().branches()[0].energy_per_decay,
        assembled.model().branches()[1].energy_per_decay
    );
    let initial =
        DecayInventory::new(Arc::clone(assembled.model()), vec![1000.0, 0.0, 0.0]).unwrap();
    let (after, report) = initial.trial(1.0, budget(), || Ok(())).unwrap();
    close(after.populations()[1], 125.0, 1e-9);
    close(after.populations()[2], 375.0, 1e-9);
    close(report.released.deposited_j, 500.0 * 10.0 * EV_J, 1e-26);
    close(report.released.radiation_j, 500.0 * 20.0 * EV_J, 1e-26);
    close(report.released.escaped_j, 500.0 * 70.0 * EV_J, 1e-26);
}

#[test]
fn incomplete_data_unknown_energy_invalid_routing_and_budget_refuse_atomically() {
    let (hydrogen, helium) = real_records();
    let construct = |entries, limits| {
        EvaluatedDecayAssembly::new("refusal-fixture".into(), entries, limits, || Ok(()))
    };
    assert!(matches!(
        construct(vec![entry(Arc::clone(&hydrogen))], DecayLimits::default()),
        Err(DecayAssemblyError::MissingDaughter { .. })
    ));
    assert!(matches!(
        construct(
            vec![entry(Arc::clone(&hydrogen)), entry(Arc::clone(&hydrogen))],
            DecayLimits::default()
        ),
        Err(DecayAssemblyError::DuplicateNuclide(_))
    ));
    let unknown = evaluation(
        identity(1, 3, 0),
        Some(1.0),
        &[Some(1.0), None, Some(0.0)],
        &[(1.0, 0, 1.0, 2.0)],
    );
    assert!(matches!(
        construct(
            vec![entry(unknown), entry(Arc::clone(&helium))],
            DecayLimits::default()
        ),
        Err(DecayAssemblyError::UnknownEnergy {
            component: "EEM",
            ..
        })
    ));
    let mut absent_neutrino = entry(Arc::clone(&hydrogen));
    absent_neutrino.routing.include_reported_neutrinos = true;
    assert!(matches!(
        construct(
            vec![absent_neutrino, entry(Arc::clone(&helium))],
            DecayLimits::default()
        ),
        Err(DecayAssemblyError::UnknownEnergy { .. })
    ));
    let mut invalid = entry(Arc::clone(&hydrogen));
    invalid.routing.light_particles.deposited = 0.5;
    assert!(matches!(
        construct(
            vec![invalid, entry(Arc::clone(&helium))],
            DecayLimits::default()
        ),
        Err(DecayAssemblyError::Invalid(_))
    ));
    let limits = DecayLimits {
        maximum_edges: 0,
        ..DecayLimits::default()
    };
    assert!(matches!(
        construct(
            vec![entry(Arc::clone(&hydrogen)), entry(Arc::clone(&helium))],
            limits
        ),
        Err(DecayAssemblyError::Budget(_))
    ));
    let bad_sum = evaluation(
        identity(1, 3, 0),
        Some(1.0),
        &[Some(1.0), Some(0.0), Some(0.0)],
        &[(1.0, 0, 0.5, 2.0)],
    );
    assert!(matches!(
        construct(
            vec![entry(bad_sum), entry(Arc::clone(&helium))],
            DecayLimits::default()
        ),
        Err(DecayAssemblyError::Decay(_))
    ));
    let mut calls = 0;
    let cancelled = EvaluatedDecayAssembly::new(
        "cancel".into(),
        vec![entry(Arc::clone(&hydrogen)), entry(helium)],
        DecayLimits::default(),
        || {
            calls += 1;
            if calls == 5 {
                Err(FissionError::Cancelled)
            } else {
                Ok(())
            }
        },
    );
    assert_eq!(
        cancelled.unwrap_err(),
        DecayAssemblyError::Physics(FissionError::Cancelled)
    );
    assert_eq!(hydrogen.law().half_life().unwrap().value(), 3.887896e8);
}

#[test]
fn imported_decay_and_capture_share_original_cyclic_burnup_and_calorimetry() {
    let (hydrogen, helium) = real_records();
    let half_life = hydrogen.law().half_life().unwrap().value();
    let lambda = hydrogen.law().decay_constant().value();
    let assembled = assembly(vec![entry(helium), entry(hydrogen)]);
    let groups = Arc::new(
        NeutronGroups::new("manufactured-one-group-exposure".into(), vec![0.0, 1e-12]).unwrap(),
    );
    let mut data = assembled
        .burnup_data(Arc::clone(&groups), 200.0, 1000.0, || Ok(()))
        .unwrap();
    assert!(data.reactions[0].data_id.contains("ENDF/B-VIII.0 H3"));
    let h3 = assembled.index(identity(1, 3, 0)).unwrap();
    let he3 = assembled.index(identity(2, 3, 0)).unwrap();
    data.reactions.push(BurnupReaction {
        data_id: "manufactured He3(n,p)H3 partial cross section and prompt energy".into(),
        parent: he3,
        rate: BurnupRate::Capture {
            cross_sections_m2: vec![1e-28],
        },
        products: vec![BurnupProduct {
            nuclide: h3,
            yield_per_event: 1.0,
        }],
        energy_per_event: EnergyChannels {
            deposited_j: 100.0 * EV_J,
            ..EnergyChannels::default()
        },
    });
    let limits = BurnupLimits::default();
    let model = Arc::new(BurnupModel::new(data, limits).unwrap());
    let initial = BurnupInventory::new(model, vec![0.0, 1000.0]).unwrap();
    let flux = GroupFlux::new(groups, vec![lambda / 1e-28]).unwrap();
    let (after, report) = initial
        .trial(
            half_life,
            &flux,
            300.0,
            BurnupBudget {
                limits,
                nuclei_tolerance: 1e-7,
                energy_tolerance_j: 1e-20,
            },
            || Ok(()),
        )
        .unwrap();
    // Equal forward/backward rates: H3 = N0(1+exp(-2 lambda t))/2.
    let decays = 500.0 * core::f64::consts::LN_2 + 187.5;
    let captures = 500.0 * core::f64::consts::LN_2 - 187.5;
    close(after.populations()[h3], 625.0, 1e-8);
    close(after.populations()[he3], 375.0, 1e-8);
    close(report.decays, decays, 1e-8);
    close(report.captures, captures, 1e-8);
    close(
        report.released.deposited_j,
        (decays * 5690.0 + captures * 100.0) * EV_J,
        1e-23,
    );
    assert_eq!(initial.populations(), &[0.0, 1000.0]);
}
