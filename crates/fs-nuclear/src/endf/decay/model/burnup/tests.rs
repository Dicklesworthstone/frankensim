//! Manufactured ENDF records and independent analytic event/energy solutions.
//! The coefficients stress multi-product semantics, not an empirical isotope
//! or cross-section validation. The original native inventory owns evolution.

use super::*;
use crate::EnergyChannels;
use crate::burnup::{BurnupBudget, BurnupInventory};
use crate::burnup::transport::{
    BurnupTransportBudget, BurnupTransportLimits, BurnupTransportState,
    NuclideTransfer, NuclideTransport,
};
use crate::endf::decay::{EndfDecayReadLimits, EndfDecaySelection};
use crate::endf::decay::model::{DecayEnergyFractions, EvaluatedDecayAssembly};
use crate::endf::fission_yields::{EndfFissionYieldReadLimits, EndfFissionYieldSelection};
use crate::irradiation::GroupFlux;

const EV_J: f64 = 1.602_176_634e-19;

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:.17e} != {expected:.17e}; allowance {tolerance:.3e}"
    );
}

fn identity(z: u16, a: u16, isomer: u16) -> NuclideId {
    NuclideId::new(z, a, isomer).unwrap()
}

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

fn row(fields: [String; 6], mf: i32, mt: i32) -> String {
    format!("{}{:4}{mf:2}{mt:3}{:5}", fields.concat(), 777, 1)
}

fn control(first: f64, second: f64, integers: [i32; 4], mf: i32, mt: i32) -> String {
    row(
        [
            real(first),
            real(second),
            format!("{:11}", integers[0]),
            format!("{:11}", integers[1]),
            format!("{:11}", integers[2]),
            format!("{:11}", integers[3]),
        ],
        mf,
        mt,
    )
}

fn za(id: NuclideId) -> f64 {
    1000.0 * f64::from(id.atomic_number()) + f64::from(id.mass_number())
}

fn evaluation(
    id: NuclideId,
    half_life_s: Option<f64>,
    energies_ev: [f64; 3],
    modes: &[(f64, f64)],
) -> Arc<DecayEvaluation> {
    evaluation_with_rfs(id, half_life_s, energies_ev, modes, 0)
}

fn evaluation_with_rfs(
    id: NuclideId,
    half_life_s: Option<f64>,
    energies_ev: [f64; 3],
    modes: &[(f64, f64)],
    rfs: u16,
) -> Arc<DecayEvaluation> {
    let mut lines = vec![
        control(
            za(id),
            f64::from(id.mass_number()),
            [
                i32::from(id.isomer()),
                i32::from(id.isomer()),
                i32::from(half_life_s.is_none()),
                0,
            ],
            8,
            457,
        ),
        control(half_life_s.unwrap_or(0.0), 0.0, [0, 0, 6, 0], 8, 457),
        row(
            [energies_ev[0], 0.0, energies_ev[1], 0.0, energies_ev[2], 0.0].map(real),
            8,
            457,
        ),
        control(
            0.0,
            1.0,
            [0, 0, if modes.is_empty() { 6 } else { 6 * modes.len() as i32 }, modes.len() as i32],
            8,
            457,
        ),
    ];
    if modes.is_empty() {
        lines.push(row([0.0; 6].map(real), 8, 457));
    } else {
        for &(rtyp, probability) in modes {
            // Q intentionally differs from the averaged emitted energies.
            lines.push(row([rtyp, f64::from(rfs), 1000.0, 0.0, probability, 0.0].map(real), 8, 457));
        }
    }
    lines.push(control(0.0, 0.0, [0; 4], 8, 0));
    Arc::new(
        DecayEvaluation::parse(
            "manufactured-decay-MF8".into(),
            EndfDecaySelection { material_number: 777, nuclide: id },
            &lines.join("\n"),
            EndfDecayReadLimits::default(),
            || Ok(()),
        )
        .unwrap(),
    )
}

fn stable(id: NuclideId) -> Arc<DecayEvaluation> {
    evaluation(id, None, [0.0; 3], &[])
}

fn yields(parent: NuclideId) -> Arc<SpontaneousFissionYieldEvaluation> {
    let mut lines = vec![
        control(za(parent), f64::from(parent.mass_number()), [-1, 0, 0, 3], 1, 451),
        control(0.0, 0.0, [i32::from(parent.isomer()), i32::from(parent.isomer()), 0, 6], 1, 451),
        control(0.0, 0.0, [0, 0, 5, 8], 1, 451),
        control(0.0, 0.0, [0; 4], 1, 451),
        control(0.0, 0.0, [0; 4], 1, 0),
        control(za(parent), f64::from(parent.mass_number()), [1, 0, 0, 0], 8, 454),
        control(0.0, 0.0, [0, 0, 8, 2], 8, 454),
    ];
    // A deliberately nonunit/non-two sum exposes illicit normalization.
    let values = [
        za(identity(39, 99, 0)), 0.0, 0.5, 0.01,
        za(identity(41, 101, 0)), 0.0, 1.75, 0.02,
    ];
    for chunk in values.chunks(6) {
        lines.push(row(
            std::array::from_fn(|index| real(chunk.get(index).copied().unwrap_or(0.0))),
            8,
            454,
        ));
    }
    lines.push(control(0.0, 0.0, [0; 4], 8, 0));
    Arc::new(
        SpontaneousFissionYieldEvaluation::parse(
            "manufactured-independent-SF-MF8".into(),
            EndfFissionYieldSelection { material_number: 777, nuclide: parent },
            &lines.join("\n"),
            EndfFissionYieldReadLimits::default(),
            || Ok(()),
        )
        .unwrap(),
    )
}

fn routing() -> DecayEnergyRouting {
    DecayEnergyRouting {
        light_particles: DecayEnergyFractions { deposited: 1.0, radiation: 0.0, escaped: 0.0 },
        electromagnetic: DecayEnergyFractions { deposited: 0.0, radiation: 1.0, escaped: 0.0 },
        heavy_particles: DecayEnergyFractions { deposited: 0.0, radiation: 0.0, escaped: 1.0 },
        include_reported_neutrinos: false,
    }
}

fn entry(evaluation: Arc<DecayEvaluation>) -> DecayModelEntry {
    DecayModelEntry { evaluation, routing: routing() }
}

fn groups() -> Arc<NeutronGroups> {
    Arc::new(NeutronGroups::new("manufactured-one-group".into(), vec![0.0, 4.0 * EV_J]).unwrap())
}

fn entries(retain_second_fragment: bool, retain_helium: bool) -> Vec<DecayModelEntry> {
    let mut entries = vec![
        entry(evaluation(
            identity(80, 200, 0), Some(2.0), [1.0, 2.0, 3.0], &[(6.0, 0.25), (4.0, 0.75)],
        )),
        entry(stable(identity(78, 196, 0))),
        entry(evaluation(
            identity(39, 99, 0), Some(1.0), [2.0, 1.0, 0.0], &[(1.0, 1.0)],
        )),
        entry(stable(identity(40, 99, 0))),
    ];
    if retain_second_fragment {
        entries.push(entry(stable(identity(41, 101, 0))));
    }
    if retain_helium {
        entries.push(entry(stable(identity(2, 4, 0))));
    }
    entries
}

fn assembly(export_second_fragment: bool, retain_helium: bool) -> EvaluatedBurnupDecayAssembly {
    let parent = identity(80, 200, 0);
    EvaluatedBurnupDecayAssembly::new(
        "manufactured-evaluated-SF".into(),
        entries(!export_second_fragment, retain_helium),
        vec![SpontaneousFissionProductBinding {
            parent,
            yields: yields(parent),
            exported_products: if export_second_fragment { vec![identity(41, 101, 0)] } else { vec![] },
        }],
        groups(),
        [300.0, 600.0],
        BurnupLimits::default(),
        || Ok(()),
    )
    .unwrap()
}

fn budget() -> BurnupBudget {
    BurnupBudget {
        limits: BurnupLimits::default(),
        nuclei_tolerance: 1e-8,
        energy_tolerance_j: 1e-24,
    }
}

fn source_budget() -> SpontaneousFissionReportBudget {
    SpontaneousFissionReportBudget { max_values: 64, max_work: 4096, nuclei_tolerance: 1e-8 }
}

fn initial(model: &Arc<BurnupModel>) -> BurnupInventory {
    let mut populations = vec![0.0; model.nuclides().len()];
    populations[0] = 100.0;
    BurnupInventory::new(Arc::clone(model), populations).unwrap()
}

fn flux(model: &Arc<BurnupModel>, value: f64) -> GroupFlux {
    GroupFlux::new(Arc::clone(model.groups()), vec![value]).unwrap()
}

#[test]
fn independent_fragments_decay_in_the_same_interval_with_parent_average_energy_once() {
    let source = Arc::new(assembly(false, false));
    let binding = source.bind_burnup(Arc::clone(source.model()), || Ok(())).unwrap();
    let original = initial(source.model());
    let (after, step) = original.trial(2.0, &flux(source.model(), 0.0), 450.0, budget(), || Ok(())).unwrap();
    let report = binding.from_burnup_step(source.model(), &step, source_budget(), || Ok(())).unwrap();

    // lambda_P=ln(2)/2, lambda_A=2*lambda_P. A(t) is independently
    // 100*(0.25*0.5)*(exp(-lambda_P*t)-exp(-2*lambda_P*t)).
    for (&actual, expected) in after.populations().iter().zip([50.0, 37.5, 3.125, 3.125, 21.875]) {
        close(actual, expected, 1e-9);
    }
    close(report.fissions[0], 12.5, 1e-9);
    assert!(report.exported_product_births[0].is_empty());
    close(after.populations().iter().sum(), 115.625, 1e-9);
    close(step.decays, 53.125, 1e-9);
    assert_eq!(step.fissions, 0.0); // Native fissions labels neutron-induced events.
    close(step.released.deposited_j, 56.25 * EV_J, 1e-27);
    close(step.released.radiation_j, 103.125 * EV_J, 1e-27);
    close(step.released.escaped_j, 150.0 * EV_J, 1e-27);
    assert_eq!(source.model().reactions()[0].energy_per_event, source.model().reactions()[1].energy_per_event);
    assert_eq!(source.model().reactions()[0].products[0].yield_per_event, 0.5);
    assert_eq!(source.model().reactions()[0].products[1].yield_per_event, 1.75);
    assert_eq!(source.spontaneous_fission_products()[0].yields.values()[0].uncertainty, 0.01);
    assert_eq!(source.mode_sources()[0].spontaneous_fission_binding, Some(0));
    assert_eq!(source.mode_sources()[1].spontaneous_fission_binding, None);
    assert_eq!(source.index(identity(41, 101, 0)), Some(4));
    assert_eq!(original.populations(), &[100.0, 0.0, 0.0, 0.0, 0.0]);

    let (half, _) = original.trial(1.0, &flux(source.model(), 0.0), 450.0, budget(), || Ok(())).unwrap();
    let (split, _) = half.trial(1.0, &flux(source.model(), 0.0), 450.0, budget(), || Ok(())).unwrap();
    for (&actual, &expected) in split.populations().iter().zip(after.populations()) {
        close(actual, expected, 1e-9);
    }
    close(split.cumulative_released().total_j(), after.cumulative_released().total_j(), 1e-27);
}

#[test]
fn explicit_exports_account_for_every_unrepresented_independent_fragment() {
    let source = Arc::new(assembly(true, false));
    let bound = source.bind_burnup(Arc::clone(source.model()), || Ok(())).unwrap();
    let (after, step) = initial(source.model())
        .trial(2.0, &flux(source.model(), 0.0), 450.0, budget(), || Ok(()))
        .unwrap();
    let report = bound.from_burnup_step(source.model(), &step, source_budget(), || Ok(())).unwrap();
    close(report.fissions[0], 12.5, 1e-9);
    close(report.exported_product_births[0][0], 21.875, 1e-9);
    close(after.populations().iter().sum::<f64>() + report.exported_product_births[0][0], 115.625, 1e-9);
    assert_eq!(source.index(identity(41, 101, 0)), None);
    assert_eq!(source.model().reactions()[0].products.len(), 1);

    let construct = |exported_products| {
        EvaluatedBurnupDecayAssembly::new(
            "export-boundary".into(),
            entries(false, false),
            vec![SpontaneousFissionProductBinding {
                parent: identity(80, 200, 0),
                yields: yields(identity(80, 200, 0)),
                exported_products,
            }],
            groups(), [300.0, 600.0], BurnupLimits::default(), || Ok(()),
        )
    };
    assert!(matches!(construct(vec![]), Err(DecayAssemblyError::MissingDaughter { .. })));
    assert!(construct(vec![identity(41, 101, 0), identity(41, 101, 0)]).is_err());
    assert!(construct(vec![identity(39, 99, 0)]).is_err());
    assert!(construct(vec![identity(42, 101, 0)]).is_err());
}

#[test]
fn reordered_sources_channels_and_products_with_capture_use_actual_spontaneous_events() {
    let first_parent = identity(80, 200, 0);
    let second_parent = identity(82, 204, 0);
    let mut evaluations = entries(true, false);
    evaluations.push(entry(evaluation(second_parent, Some(1.0), [4.0, 0.0, 0.0], &[(6.0, 1.0)])));
    // Supplied yield bindings intentionally reverse their evaluated parent order.
    let source = Arc::new(EvaluatedBurnupDecayAssembly::new(
        "reordered-SF-sources".into(), evaluations,
        vec![
            SpontaneousFissionProductBinding { parent: second_parent, yields: yields(second_parent), exported_products: vec![] },
            SpontaneousFissionProductBinding { parent: first_parent, yields: yields(first_parent), exported_products: vec![] },
        ],
        groups(), [300.0, 600.0], BurnupLimits::default(), || Ok(()),
    ).unwrap());
    let mut data = source.burnup_data(|| Ok(())).unwrap();
    data.reactions[0].products.reverse();
    data.reactions[3].products.reverse();
    data.reactions.push(BurnupReaction {
        data_id: "manufactured-competing-capture".into(),
        parent: 0,
        rate: BurnupRate::Capture { cross_sections_m2: vec![0.02] },
        products: vec![BurnupProduct { nuclide: 1, yield_per_event: 1.0 }],
        energy_per_event: EnergyChannels::default(),
    });
    data.reactions.reverse();
    let model = Arc::new(BurnupModel::new(data, BurnupLimits::default()).unwrap());
    let bound = source.bind_burnup(Arc::clone(&model), || Ok(())).unwrap();
    assert_eq!(bound.reaction_indices(), &[4, 3, 2, 1]);
    assert_eq!(source.mode_sources()[0].spontaneous_fission_binding, Some(1));
    assert_eq!(source.mode_sources()[3].spontaneous_fission_binding, Some(0));
    let mut populations = vec![0.0; source.nuclides().len()];
    populations[0] = 100.0;
    populations[source.index(second_parent).unwrap()] = 40.0;
    let original = BurnupInventory::new(Arc::clone(&model), populations).unwrap();
    let (_, step) = original.trial(2.0, &flux(&model, 2.0), 450.0, budget(), || Ok(())).unwrap();
    let report = bound.from_burnup_step(&model, &step, source_budget(), || Ok(())).unwrap();
    let lambda = core::f64::consts::LN_2 / 2.0;
    let expected = 100.0 * (0.25 * lambda) / (lambda + 0.04) * (-((-(lambda + 0.04) * 2.0).exp_m1()));
    close(report.fissions[0], 30.0, 1e-9);
    close(report.fissions[1], expected, 1e-9);
    assert!(step.captures > 0.0);
    assert!(report.fissions[1] < 12.5);

    let mut altered = source.burnup_data(|| Ok(())).unwrap();
    altered.reactions[0].products[0].yield_per_event += 0.125;
    let altered = Arc::new(BurnupModel::new(altered, BurnupLimits::default()).unwrap());
    assert!(source.bind_burnup(altered, || Ok(())).is_err());
    let mut altered = source.burnup_data(|| Ok(())).unwrap();
    altered.reactions[0].energy_per_event.deposited_j += EV_J;
    let altered = Arc::new(BurnupModel::new(altered, BurnupLimits::default()).unwrap());
    assert!(source.bind_burnup(altered, || Ok(())).is_err());
    assert!(bound.from_burnup_step(&Arc::new((*model).clone()), &step, source_budget(), || Ok(())).is_err());
}

#[test]
fn fully_exported_fission_products_need_no_invented_stable_inventory_or_future_heat() {
    let parent = identity(80, 200, 0);
    let source = Arc::new(EvaluatedBurnupDecayAssembly::new(
        "all-fragments-explicitly-exterior".into(),
        vec![
            entry(evaluation(parent, Some(2.0), [1.0, 2.0, 3.0], &[(6.0, 0.25), (4.0, 0.75)])),
            entry(stable(identity(78, 196, 0))),
        ],
        vec![SpontaneousFissionProductBinding {
            parent,
            yields: yields(parent),
            // Preserve declared exterior order independently of source ZA order.
            exported_products: vec![identity(41, 101, 0), identity(39, 99, 0)],
        }],
        groups(), [300.0, 600.0], BurnupLimits::default(), || Ok(()),
    ).unwrap());
    assert!(source.model().reactions()[0].products.is_empty());
    assert_eq!(source.nuclides().len(), 2);
    let bound = source.bind_burnup(Arc::clone(source.model()), || Ok(())).unwrap();
    let (after, step) = initial(source.model())
        .trial(2.0, &flux(source.model(), 0.0), 450.0, budget(), || Ok(()))
        .unwrap();
    let report = bound.from_burnup_step(source.model(), &step, source_budget(), || Ok(())).unwrap();
    close(report.fissions[0], 12.5, 1e-9);
    close(report.exported_product_births[0][0], 21.875, 1e-9);
    close(report.exported_product_births[0][1], 6.25, 1e-9);
    close(after.populations()[0], 50.0, 1e-9);
    close(after.populations()[1], 37.5, 1e-9);
    close(step.released.deposited_j, 50.0 * EV_J, 1e-27);
    close(step.released.radiation_j, 100.0 * EV_J, 1e-27);
    close(step.released.escaped_j, 150.0 * EV_J, 1e-27);
}

#[test]
fn ordinary_charged_product_lowering_is_shared_without_changing_fission_yields() {
    let policy = DecayProductRetention { alpha_particles: true, protons: false };
    let source = assembly(false, true).with_burnup_product_retention(policy, || Ok(())).unwrap();
    assert_eq!(source.burnup_product_retention(), policy);
    assert_eq!(source.model().reactions()[0].products.len(), 2);
    assert_eq!(source.model().reactions()[0].products[0].yield_per_event, 0.5);
    assert_eq!(source.model().reactions()[0].products[1].yield_per_event, 1.75);
    let (after, step) = initial(source.model())
        .trial(2.0, &flux(source.model(), 0.0), 450.0, budget(), || Ok(()))
        .unwrap();
    close(after.populations()[source.index(identity(2, 4, 0)).unwrap()], 37.5, 1e-9);
    close(step.released.deposited_j, 56.25 * EV_J, 1e-27);
    assert!(matches!(
        assembly(false, false).with_burnup_product_retention(policy, || Ok(())),
        Err(DecayAssemblyError::MissingRetainedProduct(_))
    ));

    // Ordinary-only assemblies emit exactly the existing lowering's native data.
    let original_entries = vec![
        entry(evaluation(identity(10, 20, 0), Some(2.0), [1.0, 2.0, 3.0], &[(4.0, 1.0)])),
        entry(stable(identity(8, 16, 0))),
        entry(stable(identity(2, 4, 0))),
    ];
    let groups = groups();
    let ordinary = EvaluatedDecayAssembly::new(
        "shared-ordinary".into(), original_entries.clone(), fs_la::decay::DecayLimits::default(), || Ok(()),
    )
    .unwrap()
    .with_burnup_product_retention(policy, || Ok(()))
    .unwrap()
    .burnup_data(Arc::clone(&groups), 300.0, 600.0, || Ok(()))
    .unwrap();
    let direct = EvaluatedBurnupDecayAssembly::new(
        "shared-ordinary".into(), original_entries, vec![], groups, [300.0, 600.0],
        BurnupLimits::default(), || Ok(()),
    )
    .unwrap()
    .with_burnup_product_retention(policy, || Ok(()))
    .unwrap();
    assert_eq!(ordinary, direct.burnup_data(|| Ok(())).unwrap());
}

#[test]
fn moving_parent_fissions_in_both_cells_with_local_exported_fragment_births() {
    let source = Arc::new(assembly(true, false));
    let bound = source.bind_burnup(Arc::clone(source.model()), || Ok(())).unwrap();
    let n = source.nuclides().len();
    let mut rates = vec![0.0; n];
    rates[0] = core::f64::consts::LN_2 / 2.0;
    let state = BurnupTransportState::new(
        vec![
            initial(source.model()),
            BurnupInventory::new(Arc::clone(source.model()), vec![0.0; n]).unwrap(),
        ],
        vec![1.0, 2.0],
        NuclideTransport {
            transfers: vec![NuclideTransfer { from_compartment: 0, to_compartment: 1, rates_per_s: rates }],
            outflows: vec![],
        },
        BurnupTransportLimits::default(),
        || Ok(()),
    )
    .unwrap();
    let (after, step) = state.trial(
        2.0, &[flux(source.model(), 0.0), flux(source.model(), 0.0)], &[450.0; 2],
        BurnupTransportBudget {
            limits: BurnupTransportLimits::default(), nuclei_tolerance: 1e-8, energy_tolerance_j: 1e-24,
        },
        || Ok(()),
    )
    .unwrap();
    let reports = bound.from_transport_step(source.model(), &step, source_budget(), || Ok(())).unwrap();
    close(after.inventories()[0].populations()[0], 25.0, 1e-9);
    close(after.inventories()[1].populations()[0], 25.0, 1e-9);
    close(step.transfer_counts[0][0], 37.5, 1e-9);
    // Donor total-decay events are 100/2*(1-exp(-2*lambda*2))=37.5;
    // receiver events are 50-37.5=12.5; each has the SF branching 0.25.
    close(reports[0].fissions[0], 9.375, 1e-9);
    close(reports[1].fissions[0], 3.125, 1e-9);
    close(reports[0].exported_product_births[0][0], 16.40625, 1e-9);
    close(reports[1].exported_product_births[0][0], 5.46875, 1e-9);
    assert!(bound.from_transport_step(
        source.model(), &step, SpontaneousFissionReportBudget { max_values: 3, ..source_budget() }, || Ok(()),
    ).is_err());
    assert_eq!(state.inventories()[0].populations()[0], 100.0);
}

#[test]
fn missing_sources_compound_fission_budgets_and_late_cancellation_refuse() {
    let construct = |entries, bindings, limits| {
        EvaluatedBurnupDecayAssembly::new(
            "admission-refusal".into(), entries, bindings, groups(), [300.0, 600.0], limits, || Ok(()),
        )
    };
    assert!(construct(entries(true, false), vec![], BurnupLimits::default()).is_err());
    let compound = vec![entry(evaluation(
        identity(80, 200, 0), Some(2.0), [1.0, 2.0, 3.0], &[(1.6, 1.0)],
    ))];
    assert!(matches!(
        construct(compound, vec![], BurnupLimits::default()),
        Err(DecayAssemblyError::UnsupportedDecayMode { .. })
    ));
    let bindings = || vec![SpontaneousFissionProductBinding {
        parent: identity(80, 200, 0), yields: yields(identity(80, 200, 0)), exported_products: vec![],
    }];
    let mut rfs_entries = entries(true, false);
    rfs_entries[0] = entry(evaluation_with_rfs(
        identity(80, 200, 0), Some(2.0), [1.0, 2.0, 3.0], &[(6.0, 1.0)], 1,
    ));
    assert!(matches!(
        construct(rfs_entries, bindings(), BurnupLimits::default()),
        Err(DecayAssemblyError::Invalid(_))
    ));
    let mut wrong = bindings();
    wrong[0].parent = identity(80, 200, 1);
    assert!(construct(entries(true, false), wrong, BurnupLimits::default()).is_err());
    assert!(matches!(
        construct(entries(true, false), bindings(), BurnupLimits { max_products: 1, ..BurnupLimits::default() }),
        Err(DecayAssemblyError::Budget(_))
    ));
    let source = Arc::new(assembly(true, false));
    let bound = source.bind_burnup(Arc::clone(source.model()), || Ok(())).unwrap();
    let original = initial(source.model());
    let (_, step) = original.trial(2.0, &flux(source.model(), 0.0), 450.0, budget(), || Ok(())).unwrap();
    assert!(bound.from_burnup_step(
        source.model(), &step, SpontaneousFissionReportBudget { max_work: 0, ..source_budget() }, || Ok(()),
    ).is_err());
    assert!(bound.from_burnup_step(
        source.model(), &step, SpontaneousFissionReportBudget { max_values: 1, ..source_budget() }, || Ok(()),
    ).is_err());
    let mut polls = 0;
    let accepted = bound.from_burnup_step(source.model(), &step, source_budget(), || { polls += 1; Ok(()) }).unwrap();
    let mut attempt = 0;
    assert_eq!(
        bound.from_burnup_step(source.model(), &step, source_budget(), || {
            attempt += 1;
            if attempt == polls { Err(FissionError::Cancelled) } else { Ok(()) }
        }),
        Err(DecayAssemblyError::Physics(FissionError::Cancelled))
    );
    assert_eq!(bound.from_burnup_step(source.model(), &step, source_budget(), || Ok(())).unwrap(), accepted);
    let mut invalid_step = step.clone();
    invalid_step.reaction_counts[0] = f64::NAN;
    assert!(bound.from_burnup_step(source.model(), &invalid_step, source_budget(), || Ok(())).is_err());
    let (_, zero_step) = original.trial(0.0, &flux(source.model(), 0.0), 450.0, budget(), || Ok(())).unwrap();
    let zero = bound.from_burnup_step(source.model(), &zero_step, source_budget(), || Ok(())).unwrap();
    assert_eq!(zero.fissions, vec![0.0]);
    assert_eq!(zero.exported_product_births, vec![vec![0.0]]);
    assert_eq!(original.populations(), &[100.0, 0.0, 0.0, 0.0]);
}
