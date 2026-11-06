//! Manufactured ENDF records and independent reaction/decay solutions.
//! These coefficients establish equations and wiring, not isotope accuracy.

use super::model::FissionYieldCollapseLimits;
use super::*;
use crate::EnergyChannels;
use crate::burnup::{
    BurnupBudget, BurnupData, BurnupInventory, BurnupLimits, BurnupModel, BurnupNuclide,
    BurnupProduct, BurnupRate, BurnupReaction,
};
use crate::cross_sections::NeutronReaction;
use crate::endf::spectrum::{FluxSpectrum, FluxSpectrumData, FluxSpectrumLimits};
use crate::endf::{EndfInterpolationRegion, EndfMetadata, PointwiseCrossSection};
use crate::irradiation::{GroupFlux, NeutronGroups};
use fs_qty::Frequency;
use fs_thermochem::decay::DecayLawV1;

fn parent() -> NuclideId {
    NuclideId::new(10, 20, 1).unwrap()
}
fn product(state: u16) -> NuclideId {
    NuclideId::new(4, 8, state).unwrap()
}
fn daughter() -> NuclideId {
    NuclideId::new(5, 8, 0).unwrap()
}
fn selection() -> EndfFissionYieldSelection {
    EndfFissionYieldSelection {
        material_number: 777,
        nuclide: parent(),
    }
}

fn real(value: f64) -> String {
    for precision in (0..=6).rev() {
        let scientific = format!("{value:.precision$e}");
        let (mantissa, exponent) = scientific.split_once('e').unwrap();
        let exponent: i32 = exponent.parse().unwrap();
        let field = format!("{mantissa}{exponent:+}");
        if field.len() <= 11 {
            return format!("{field:>11}");
        }
    }
    panic!("fixture value does not fit ENDF");
}

fn record(fields: [String; 6], mf: i32, mt: i32) -> String {
    assert!(fields.iter().all(|field| field.len() == 11));
    format!("{}{:4}{mf:2}{mt:3}{:5}", fields.concat(), 777, 1)
}

fn control(c1: f64, c2: f64, integers: [i32; 4], mf: i32, mt: i32) -> String {
    record(
        [
            real(c1),
            real(c2),
            format!("{:11}", integers[0]),
            format!("{:11}", integers[1]),
            format!("{:11}", integers[2]),
            format!("{:11}", integers[3]),
        ],
        mf,
        mt,
    )
}

fn fixture(rows: &[(f64, f64, f64)]) -> Vec<String> {
    let mut lines = vec![
        control(10020.0, 20.0, [-1, 0, 0, 3], 1, 451),
        control(0.0, 0.0, [1, 1, 0, 6], 1, 451),
        control(1.0, 3.0, [0, 0, 11, 8], 1, 451),
        control(0.0, 0.0, [0, 0, 0, 0], 1, 451),
        control(0.0, 0.0, [0, 0, 0, 0], 1, 0),
        control(10020.0, 20.0, [rows.len() as i32, 0, 0, 0], 8, 454),
    ];
    for (index, &(energy, first, second)) in rows.iter().enumerate() {
        let interpolation = if index == 0 { rows.len() as i32 - 1 } else { 2 };
        lines.push(control(energy, 0.0, [interpolation, 0, 8, 2], 8, 454));
        lines.push(record(
            [4008.0, 0.0, first, 0.01, 4008.0, 1.0].map(real),
            8,
            454,
        ));
        lines.push(record([second, 0.02, 0.0, 0.0, 0.0, 0.0].map(real), 8, 454));
    }
    lines.push(control(0.0, 0.0, [0, 0, 0, 0], 8, 0));
    lines
}

fn standard() -> Vec<String> {
    fixture(&[(1.0, 0.5, 1.5), (3.0, 1.5, 0.5)])
}
fn parse(lines: &[String]) -> Result<FissionYieldEvaluation, EndfError> {
    FissionYieldEvaluation::parse(
        "manufactured independent yields".into(),
        selection(),
        &lines.join("\n"),
        EndfFissionYieldReadLimits::default(),
        || Ok(()),
    )
}
fn field(line: &mut String, index: usize, value: &str) {
    assert!(value.len() <= 11);
    line.replace_range(index * 11..(index + 1) * 11, &format!("{value:>11}"));
}
fn near(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 3e-11 * expected.abs().max(1.0),
        "{actual:.17e} != {expected:.17e}"
    );
}

#[test]
fn independent_yields_retain_isomers_multiplicity_uncertainty_and_source_energy() {
    let evaluation = parse(&standard()).unwrap();
    assert_eq!(evaluation.metadata().nuclide, parent());
    assert_eq!(evaluation.metadata().library_version, 8);
    assert_eq!(evaluation.metadata().material_revision, 3);
    assert_eq!(evaluation.products(), &[product(0), product(1)]);
    assert_eq!(evaluation.points()[0].values[1].yield_per_fission, 1.5);
    assert_eq!(evaluation.points()[1].values[0].yield_per_fission, 1.5);
    assert_eq!(evaluation.points()[1].values[0].uncertainty, 0.01);
    assert_eq!(evaluation.points()[1].values[1].uncertainty, 0.02);
    assert_eq!(evaluation.points()[0].energy_j, EV_J);
    assert_eq!(evaluation.points()[1].energy_j, 3.0 * EV_J);
    assert!(!evaluation.is_energy_independent());
}

#[test]
fn format_identity_cumulative_and_incomplete_product_inputs_refuse() {
    let original = standard();
    for (line, index, value) in [
        (0, 0, "10021"),
        (1, 3, "0"),
        (2, 4, "5"),
        (5, 1, "21"),
        (6, 2, "0"),
        (6, 4, "7"),
        (7, 2, "-1"),
        (7, 3, "-0.1"),
        (7, 5, "0"),
        (8, 2, "1"),
        (9, 0, "1"),
        (9, 2, "5"),
        (10, 0, "4009"),
        (12, 0, "1"),
    ] {
        let mut changed = original.clone();
        field(&mut changed[line], index, value);
        assert!(
            parse(&changed).is_err(),
            "admitted changed line {line} field {index}"
        );
    }
    assert!(parse(&original[..original.len() - 1]).is_err());
    assert!(parse(&original[5..]).is_err()); // no parent isomer/sublibrary metadata
    let mut cumulative = original.clone();
    for line in &mut cumulative {
        if &line[70..75] == " 8454" {
            line.replace_range(72..75, "459");
        }
    }
    assert_eq!(parse(&cumulative), Err(EndfError::MissingSelection));
    let mut duplicated = original.clone();
    duplicated.extend(original.clone());
    assert_eq!(parse(&duplicated), Err(EndfError::AmbiguousSelection));
    let mut unequal_set = original.clone();
    field(&mut unequal_set[9], 4, "4");
    field(&mut unequal_set[9], 5, "1");
    assert!(parse(&unequal_set).is_err());
}

#[test]
fn source_work_limits_and_cancellation_are_enforced() {
    let text = standard().join("\n");
    let defaults = EndfFissionYieldReadLimits::default();
    for limits in [
        EndfFissionYieldReadLimits {
            max_bytes: text.len() - 1,
            ..defaults
        },
        EndfFissionYieldReadLimits {
            max_records: 12,
            ..defaults
        },
        EndfFissionYieldReadLimits {
            max_energy_points: 1,
            ..defaults
        },
        EndfFissionYieldReadLimits {
            max_products_per_energy: 1,
            ..defaults
        },
        EndfFissionYieldReadLimits {
            max_product_values: 3,
            ..defaults
        },
    ] {
        assert!(matches!(
            FissionYieldEvaluation::parse("fixture".into(), selection(), &text, limits, || Ok(())),
            Err(EndfError::Physics(FissionError::Budget(_)))
        ));
    }
    let mut polls = 0;
    assert_eq!(
        FissionYieldEvaluation::parse("fixture".into(), selection(), &text, defaults, || {
            polls += 1;
            if polls == 12 {
                Err(FissionError::Cancelled)
            } else {
                Ok(())
            }
        }),
        Err(EndfError::Physics(FissionError::Cancelled))
    );
}

fn grid(edges_ev: &[f64]) -> Arc<NeutronGroups> {
    Arc::new(
        NeutronGroups::new(
            "manufactured fission grid".into(),
            edges_ev.iter().map(|energy| energy * EV_J).collect(),
        )
        .unwrap(),
    )
}
fn cross(bounds_ev: [f64; 2], values: [f64; 2]) -> PointwiseCrossSection {
    PointwiseCrossSection {
        data_id: "manufactured processed fission".into(),
        metadata: EndfMetadata {
            material_number: 777,
            nuclide: parent(),
            temperature_k: 300.0,
            atomic_weight_ratio: 20.0,
            library_id: 0,
            library_version: 8,
            library_release: 0,
            material_revision: 3,
            resonance_flag: 2,
            derivation_flag: 1,
        },
        reaction: NeutronReaction::Fission,
        energies_j: bounds_ev.map(|energy| energy * EV_J).into(),
        cross_sections_m2: values.into(),
        interpolation_regions: vec![EndfInterpolationRegion {
            end_point: 2,
            interpolation: EndfInterpolation::Linear,
        }]
        .into(),
        mass_difference_q_j: 0.0,
        reaction_q_j: 0.0,
    }
}
fn spectrum(bounds_ev: [f64; 2]) -> FluxSpectrum {
    FluxSpectrum::new(
        FluxSpectrumData {
            data_id: "manufactured flat incident spectrum".into(),
            energies_j: bounds_ev.map(|energy| energy * EV_J).to_vec(),
            flux_density_m2_s_j: vec![1.0 / EV_J; 2],
            interpolation_regions: vec![EndfInterpolationRegion {
                end_point: 2,
                interpolation: EndfInterpolation::Linear,
            }],
        },
        FluxSpectrumLimits::default(),
    )
    .unwrap()
}

#[test]
fn product_condensation_preserves_reaction_integrals_and_declared_constant_laws() {
    let evaluation = parse(&standard()).unwrap();
    // With x=E/eV, sigma=0.1x and y=x/2, integral(sigma*y)/integral(sigma)=13/12.
    let weighted = evaluation
        .collapse(
            &cross([1.0, 3.0], [0.1, 0.3]),
            &spectrum([1.0, 3.0]),
            grid(&[1.0, 3.0]),
            FissionYieldCollapseLimits::default(),
            || Ok(()),
        )
        .unwrap();
    near(
        weighted.reaction_weighted().yields[0].values[0].unwrap(),
        13.0 / 12.0,
    );
    near(
        weighted.reaction_weighted().yields[1].values[0].unwrap(),
        11.0 / 12.0,
    );
    near(weighted.reaction_weighted().cross_sections_m2[0], 0.2);
    near(weighted.reaction_weighted().flux.values_m2_s()[0], 2.0);
    let constant = parse(&fixture(&[(1.0, 3.2, 0.0)])).unwrap();
    assert!(constant.is_energy_independent());
    let collapsed = constant
        .collapse(
            &cross([0.0, 5.0], [0.2, 0.2]),
            &spectrum([0.0, 5.0]),
            grid(&[0.0, 2.0, 5.0]),
            FissionYieldCollapseLimits::default(),
            || Ok(()),
        )
        .unwrap();
    for value in &collapsed.reaction_weighted().yields[0].values {
        near(value.unwrap(), 3.2);
    }
    assert!(
        evaluation
            .collapse(
                &cross([0.0, 5.0], [0.2, 0.2]),
                &spectrum([0.0, 5.0]),
                grid(&[0.0, 5.0]),
                FissionYieldCollapseLimits::default(),
                || Ok(())
            )
            .is_err()
    );
}

fn data(groups: Arc<NeutronGroups>) -> BurnupData {
    let node = |id, rate| BurnupNuclide {
        id,
        law: DecayLawV1::from_decay_constant(Frequency::new(rate)).unwrap(),
    };
    BurnupData {
        data_id: "manufactured product-resolved model".into(),
        groups,
        temperature_min_k: 300.0,
        temperature_max_k: 300.0,
        nuclides: vec![
            node(parent(), 0.0),
            node(product(0), 0.3),
            node(product(1), 0.0),
            node(daughter(), 0.0),
        ],
        reactions: vec![BurnupReaction {
            data_id: "manufactured product decay".into(),
            parent: 1,
            rate: BurnupRate::Decay {
                branch_fraction: 1.0,
            },
            products: vec![BurnupProduct {
                nuclide: 3,
                yield_per_event: 1.0,
            }],
            energy_per_event: EnergyChannels {
                deposited_j: 0.4,
                radiation_j: 0.1,
                escaped_j: 0.2,
            },
        }],
    }
}

fn prompt() -> [EnergyChannels; 2] {
    [EnergyChannels {
        deposited_j: 3.0,
        radiation_j: 0.0,
        escaped_j: 0.0,
    }; 2]
}

fn collapsed() -> super::model::FissionYieldCollapse {
    parse(&standard())
        .unwrap()
        .collapse(
            &cross([1.0, 3.0], [0.2, 0.2]),
            &spectrum([1.0, 3.0]),
            grid(&[1.0, 2.0, 3.0]),
            FissionYieldCollapseLimits::default(),
            || Ok(()),
        )
        .unwrap()
}

#[test]
fn evaluated_products_drive_actual_burnup_spectrum_response_and_shutdown_heat() {
    let collapsed = collapsed();
    let groups = Arc::clone(collapsed.reaction_weighted().flux.groups());
    let limits = BurnupLimits::default();
    let assembled = collapsed
        .append_to_burnup(
            &data(Arc::clone(&groups)),
            &prompt(),
            &[],
            limits,
            || Ok(()),
        )
        .unwrap();
    assert_eq!(assembled.reactions.len(), 3);
    for (group, reaction) in assembled.reactions[1..].iter().enumerate() {
        let BurnupRate::Fission { cross_sections_m2 } = &reaction.rate else {
            panic!()
        };
        near(cross_sections_m2[group], 0.2);
        assert_eq!(cross_sections_m2[1 - group], 0.0);
    }
    let model = Arc::new(BurnupModel::new(assembled, limits).unwrap());
    let budget = BurnupBudget {
        limits,
        nuclei_tolerance: 1e-8,
        energy_tolerance_j: 1e-8,
    };
    for (flux_values, yield_) in [(vec![2.0, 0.0], 0.75), (vec![0.0, 2.0], 1.25)] {
        let mut inventory =
            BurnupInventory::new(Arc::clone(&model), vec![100.0, 0.0, 0.0, 0.0]).unwrap();
        let flux = GroupFlux::new(Arc::clone(&groups), flux_values).unwrap();
        let report = inventory
            .advance(1.25, &flux, 300.0, budget, || Ok(()))
            .unwrap();
        let fissions = 100.0 * -(-0.5_f64).exp_m1();
        let surviving_product =
            100.0 * yield_ * 0.4 * ((-0.5_f64).exp() - (-0.375_f64).exp()) / (0.3 - 0.4);
        let decays = fissions * yield_ - surviving_product;
        near(report.fissions, fissions);
        near(inventory.populations()[0], 100.0 - fissions);
        near(inventory.populations()[1], surviving_product);
        near(inventory.populations()[2], fissions * (2.0 - yield_));
        near(inventory.populations()[3], decays);
        near(report.released.deposited_j, 3.0 * fissions + 0.4 * decays);
        near(report.released.radiation_j, 0.1 * decays);
        near(report.released.escaped_j, 0.2 * decays);
        let shutdown = GroupFlux::new(Arc::clone(&groups), vec![0.0; 2]).unwrap();
        let cooling = inventory
            .advance(2.0, &shutdown, 300.0, budget, || Ok(()))
            .unwrap();
        let later_decays = surviving_product * -(-0.6_f64).exp_m1();
        assert_eq!(cooling.fissions, 0.0);
        near(cooling.decays, later_decays);
        near(cooling.released.deposited_j, 0.4 * later_decays);
    }
}

#[test]
fn assembly_requires_explicit_products_avoids_duplicate_reactions_and_preserves_input() {
    let collapsed = collapsed();
    let groups = Arc::clone(collapsed.reaction_weighted().flux.groups());
    let source = data(groups);
    let original = source.clone();
    let limits = BurnupLimits::default();
    let assembled = collapsed
        .append_to_burnup(&source, &prompt(), &[], limits, || Ok(()))
        .unwrap();
    assert_eq!(source, original);
    assert!(
        collapsed
            .append_to_burnup(&assembled, &prompt(), &[], limits, || Ok(()))
            .is_err()
    );
    let mut missing = source.clone();
    missing.nuclides[2].id = NuclideId::new(6, 12, 0).unwrap();
    assert!(
        collapsed
            .append_to_burnup(&missing, &prompt(), &[], limits, || Ok(()))
            .is_err()
    );
    let exported = collapsed
        .append_to_burnup(&missing, &prompt(), &[product(1)], limits, || Ok(()))
        .unwrap();
    assert!(
        exported.reactions[1..]
            .iter()
            .all(|reaction| reaction.products.len() == 1)
    );
    assert!(
        collapsed
            .append_to_burnup(&source, &prompt(), &[product(1)], limits, || Ok(()))
            .is_err()
    );
    let mut temperature = source.clone();
    temperature.temperature_max_k = 500.0;
    assert!(
        collapsed
            .append_to_burnup(&temperature, &prompt(), &[], limits, || Ok(()))
            .is_err()
    );
    assert!(
        collapsed
            .append_to_burnup(
                &source,
                &prompt(),
                &[],
                BurnupLimits {
                    max_products: 2,
                    ..limits
                },
                || Ok(())
            )
            .is_err()
    );
    let mut polls = 0;
    assert_eq!(
        collapsed.append_to_burnup(&source, &prompt(), &[], limits, || {
            polls += 1;
            if polls == 15 {
                Err(FissionError::Cancelled)
            } else {
                Ok(())
            }
        }),
        Err(EndfError::Physics(FissionError::Cancelled))
    );
    assert_eq!(source, original);
}
