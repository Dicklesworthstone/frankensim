//! Independent temperature moments and actual finite-inventory evolution.
//! The manufactured coefficients test equations, not an evaluated isotope.

use super::*;
use crate::burnup::{BurnupBudget, BurnupInventory, BurnupModel, BurnupNuclide, BurnupReaction};
use crate::endf::fission_yields::model::FissionYieldCollapseLimits;
use crate::endf::fission_yields::{FissionYieldMetadata, FissionYieldPoint, FissionYieldValue};
use crate::endf::spectrum::{FluxSpectrumData, FluxSpectrumLimits};
use crate::endf::{EV_J, EndfInterpolation, EndfInterpolationRegion, PointwiseCrossSection};
use crate::irradiation::GroupFlux;
use fs_qty::Frequency;
use fs_thermochem::decay::DecayLawV1;

fn parent() -> NuclideId {
    NuclideId::new(10, 20, 1).unwrap()
}
fn product(index: usize) -> NuclideId {
    NuclideId::new(4 + index as u16, 8, 0).unwrap()
}
fn grid(edges_ev: &[f64]) -> Arc<NeutronGroups> {
    Arc::new(
        NeutronGroups::new(
            "temperature-series incident groups".into(),
            edges_ev.iter().map(|energy| energy * EV_J).collect(),
        )
        .unwrap(),
    )
}
fn evaluation() -> FissionYieldEvaluation {
    FissionYieldEvaluation {
        data_id: "temperature-series independent evaluation".into(),
        metadata: FissionYieldMetadata {
            material_number: 777,
            nuclide: parent(),
            atomic_weight_ratio: 20.0,
            library_id: 0,
            library_version: 8,
            library_release: 0,
            material_revision: 3,
        },
        products: vec![product(0), product(1)].into(),
        points: [(1.0, [0.5, 1.5]), (3.0, [1.5, 0.5])]
            .into_iter()
            .enumerate()
            .map(|(index, (energy, yields))| FissionYieldPoint {
                energy_j: energy * EV_J,
                interpolation_from_previous: (index > 0).then_some(EndfInterpolation::Linear),
                values: yields
                    .map(|yield_per_fission| FissionYieldValue {
                        yield_per_fission,
                        uncertainty: 0.01,
                    })
                    .into(),
            })
            .collect::<Vec<_>>()
            .into(),
    }
}
fn cross(temperature_k: f64, coefficients: [f64; 2]) -> PointwiseCrossSection {
    PointwiseCrossSection {
        data_id: format!("temperature-series processed fission at {temperature_k} K").into(),
        metadata: EndfMetadata {
            material_number: 777,
            nuclide: parent(),
            temperature_k,
            atomic_weight_ratio: 20.0,
            library_id: 0,
            library_version: 8,
            library_release: 0,
            material_revision: 3,
            resonance_flag: 2,
            derivation_flag: 1,
        },
        reaction: NeutronReaction::Fission,
        energies_j: [EV_J, 3.0 * EV_J].into(),
        cross_sections_m2: coefficients.into(),
        interpolation_regions: vec![EndfInterpolationRegion {
            end_point: 2,
            interpolation: EndfInterpolation::Linear,
        }]
        .into(),
        mass_difference_q_j: 0.0,
        reaction_q_j: 0.0,
    }
}
fn spectrum(ordinates: [f64; 2]) -> FluxSpectrum {
    FluxSpectrum::new(
        FluxSpectrumData {
            data_id: "temperature-series original spectrum".into(),
            energies_j: vec![EV_J, 3.0 * EV_J],
            flux_density_m2_s_j: ordinates.map(|value| value / EV_J).to_vec(),
            interpolation_regions: vec![EndfInterpolationRegion {
                end_point: 2,
                interpolation: EndfInterpolation::Linear,
            }],
        },
        FluxSpectrumLimits::default(),
    )
    .unwrap()
}
fn collapse(
    evaluation: &FissionYieldEvaluation,
    cross: &PointwiseCrossSection,
    spectrum: &FluxSpectrum,
    groups: &Arc<NeutronGroups>,
) -> FissionYieldCollapse {
    evaluation
        .collapse(
            cross,
            spectrum,
            Arc::clone(groups),
            FissionYieldCollapseLimits::default(),
            || Ok(()),
        )
        .unwrap()
}
fn source_data(groups: &Arc<NeutronGroups>) -> BurnupData {
    BurnupData {
        data_id: "temperature-series finite inventory".into(),
        groups: Arc::clone(groups),
        temperature_min_k: 300.0,
        temperature_max_k: 600.0,
        nuclides: [parent(), product(0), product(1)]
            .into_iter()
            .map(|id| BurnupNuclide {
                id,
                law: DecayLawV1::from_decay_constant(Frequency::new(0.0)).unwrap(),
            })
            .collect(),
        reactions: vec![BurnupReaction {
            data_id: "original competing capture".into(),
            parent: 0,
            rate: BurnupRate::Capture {
                cross_sections_m2: vec![0.05; groups.len()],
            },
            products: vec![],
            energy_per_event: energy(0.1, 0.02, 0.03),
        }],
    }
}
fn energy(deposited_j: f64, radiation_j: f64, escaped_j: f64) -> EnergyChannels {
    EnergyChannels {
        deposited_j,
        radiation_j,
        escaped_j,
    }
}
fn budget() -> BurnupBudget {
    BurnupBudget {
        limits: BurnupLimits::default(),
        nuclei_tolerance: 1e-8,
        energy_tolerance_j: 1e-8,
    }
}
fn near(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 3e-11 * expected.abs().max(1e-12),
        "{actual:.17e} != {expected:.17e}"
    );
}
fn standard(groups: &Arc<NeutronGroups>) -> [FissionYieldCollapse; 2] {
    let evaluation = evaluation();
    let spectrum = spectrum([1.0, 1.0]);
    [
        collapse(&evaluation, &cross(300.0, [0.1, 0.3]), &spectrum, groups),
        collapse(&evaluation, &cross(600.0, [0.8, 0.0]), &spectrum, groups),
    ]
}
fn knots<'a>(
    yields: &'a [FissionYieldCollapse; 2],
    prompt: &'a [Vec<EnergyChannels>; 2],
) -> [FissionYieldTemperatureKnot<'a>; 2] {
    [
        FissionYieldTemperatureKnot {
            yields: &yields[0],
            prompt_energy_per_group: &prompt[0],
        },
        FissionYieldTemperatureKnot {
            yields: &yields[1],
            prompt_energy_per_group: &prompt[1],
        },
    ]
}

#[test]
fn temperature_moments_drive_actual_removal_products_and_all_prompt_channels() {
    let groups = grid(&[1.0, 3.0]);
    let yields = standard(&groups);
    let prompt = [vec![energy(2.0, 0.4, 0.2)], vec![energy(5.0, 1.0, 0.8)]];
    let source = source_data(&groups);
    let original = source.clone();
    let assembled = append_temperature_series_to_burnup(
        &knots(&yields, &prompt),
        &source,
        &[],
        FissionYieldTemperatureLimits::default(),
        || Ok(()),
    )
    .unwrap();
    assert_eq!(assembled.reactions.len(), 3);
    assert_eq!(assembled.reactions[0], source.reactions[0]);
    assert_eq!(source, original);
    for reaction in &assembled.reactions[1..] {
        let BurnupRate::TabulatedFission { cross_sections } = &reaction.rate else {
            panic!("temperature assembly did not use original tabulated fission")
        };
        assert_eq!(cross_sections.temperatures_k(), &[300.0, 600.0]);
    }
    let model = Arc::new(BurnupModel::new(assembled, budget().limits).unwrap());
    let flux = GroupFlux::new(Arc::clone(&groups), vec![2.0]).unwrap();
    // Independent integrals on 1..3 eV give (sigma,Y0)=(.2,13/12),(.4,5/6).
    // At the middle temperature sigma*Y and sigma*Q interpolate, not Y and Q.
    for (temperature, sigma, first_yield, q) in [
        (300.0, 0.2, 13.0 / 12.0, energy(2.0, 0.4, 0.2)),
        (450.0, 0.3, 11.0 / 12.0, energy(4.0, 0.8, 0.6)),
        (600.0, 0.4, 5.0 / 6.0, energy(5.0, 1.0, 0.8)),
    ] {
        let mut inventory =
            BurnupInventory::new(Arc::clone(&model), vec![100.0, 0.0, 0.0]).unwrap();
        let duration = 0.7;
        let total_rate: f64 = 2.0 * sigma + 0.1;
        let removals = 100.0 * -(-total_rate * duration).exp_m1();
        let fissions = removals * (2.0 * sigma) / total_rate;
        let captures = removals * 0.1 / total_rate;
        let report = inventory
            .advance(duration, &flux, temperature, budget(), || Ok(()))
            .unwrap();
        near(report.fissions, fissions);
        near(report.captures, captures);
        near(inventory.populations()[0], 100.0 - removals);
        near(inventory.populations()[1], fissions * first_yield);
        near(inventory.populations()[2], fissions * (2.0 - first_yield));
        near(
            report.released.deposited_j,
            fissions * q.deposited_j + 0.1 * captures,
        );
        near(
            report.released.radiation_j,
            fissions * q.radiation_j + 0.02 * captures,
        );
        near(
            report.released.escaped_j,
            fissions * q.escaped_j + 0.03 * captures,
        );
    }
}

#[test]
fn zero_reaction_temperature_is_a_retained_anchor_not_an_interpolation_gap() {
    let groups = grid(&[1.0, 3.0]);
    let evaluation = evaluation();
    let spectrum = spectrum([1.0, 1.0]);
    let low = collapse(&evaluation, &cross(300.0, [0.2, 0.2]), &spectrum, &groups);
    let zero = collapse(&evaluation, &cross(450.0, [0.0, 0.0]), &spectrum, &groups);
    let high = collapse(&evaluation, &cross(600.0, [0.4, 0.4]), &spectrum, &groups);
    let prompt = [energy(1.0, 0.0, 0.0)];
    let knots = [&low, &zero, &high].map(|yields| FissionYieldTemperatureKnot {
        yields,
        prompt_energy_per_group: &prompt,
    });
    let mut data = source_data(&groups);
    data.reactions.clear();
    let assembled = append_temperature_series_to_burnup(
        &knots,
        &data,
        &[],
        FissionYieldTemperatureLimits::default(),
        || Ok(()),
    )
    .unwrap();
    assert_eq!(assembled.reactions.len(), 2);
    let mut coefficients = vec![0.0; 5];
    for reaction in &assembled.reactions {
        let BurnupRate::TabulatedFission { cross_sections } = &reaction.rate else {
            panic!()
        };
        assert_eq!(cross_sections.temperatures_k(), &[300.0, 450.0, 600.0]);
        for (sum, temperature) in coefficients
            .iter_mut()
            .zip([300.0, 375.0, 450.0, 525.0, 600.0])
        {
            *sum += cross_sections
                .cross_section_at_group(temperature, 0)
                .unwrap();
        }
    }
    for (actual, expected) in coefficients.into_iter().zip([0.2, 0.1, 0.0, 0.2, 0.4]) {
        near(actual, expected);
    }
    let model = Arc::new(BurnupModel::new(assembled, budget().limits).unwrap());
    let mut inventory = BurnupInventory::new(model, vec![100.0, 0.0, 0.0]).unwrap();
    let flux = GroupFlux::new(groups, vec![2.0]).unwrap();
    let report = inventory
        .advance(1.0, &flux, 450.0, budget(), || Ok(()))
        .unwrap();
    assert_eq!(report.fissions, 0.0);
    assert_eq!(report.released.total_j(), 0.0);
    assert_eq!(inventory.populations(), &[100.0, 0.0, 0.0]);
}

#[test]
fn temperature_channels_keep_incident_groups_disjoint_under_spectrum_changes() {
    let groups = grid(&[1.0, 2.0, 3.0]);
    let evaluation = evaluation();
    let spectrum = spectrum([1.0, 1.0]);
    let yields = [
        collapse(&evaluation, &cross(300.0, [0.2, 0.2]), &spectrum, &groups),
        collapse(&evaluation, &cross(600.0, [0.4, 0.4]), &spectrum, &groups),
    ];
    let prompt = [
        vec![energy(1.0, 0.0, 0.0), energy(2.0, 0.0, 0.0)],
        vec![energy(3.0, 0.0, 0.0), energy(4.0, 0.0, 0.0)],
    ];
    let mut data = source_data(&groups);
    data.reactions.clear();
    let assembled = append_temperature_series_to_burnup(
        &knots(&yields, &prompt),
        &data,
        &[],
        FissionYieldTemperatureLimits::default(),
        || Ok(()),
    )
    .unwrap();
    assert_eq!(assembled.reactions.len(), 4);
    let model = Arc::new(BurnupModel::new(assembled, budget().limits).unwrap());
    for (flux, yield_, prompt) in [
        (vec![2.0, 0.0], 0.75, 7.0 / 3.0),
        (vec![0.0, 2.0], 1.25, 10.0 / 3.0),
    ] {
        let mut inventory =
            BurnupInventory::new(Arc::clone(&model), vec![100.0, 0.0, 0.0]).unwrap();
        let flux = GroupFlux::new(Arc::clone(&groups), flux).unwrap();
        let report = inventory
            .advance(0.5, &flux, 450.0, budget(), || Ok(()))
            .unwrap();
        let fissions = 100.0 * -(-0.3_f64).exp_m1();
        near(report.fissions, fissions);
        near(inventory.populations()[1], fissions * yield_);
        near(report.released.deposited_j, fissions * prompt);
    }
}

#[test]
fn temperature_domain_original_evaluation_and_full_spectrum_must_match() {
    let groups = grid(&[1.0, 3.0]);
    let yields = standard(&groups);
    let prompt = std::array::from_fn(|_| vec![energy(1.0, 0.0, 0.0)]);
    let source = source_data(&groups);
    let limits = FissionYieldTemperatureLimits::default();
    let inputs = knots(&yields, &prompt);
    for invalid_knots in [
        vec![],
        vec![inputs[1], inputs[0]],
        vec![inputs[0], inputs[0]],
    ] {
        assert!(
            append_temperature_series_to_burnup(&invalid_knots, &source, &[], limits, || Ok(()))
                .is_err()
        );
    }
    for (minimum, maximum) in [
        (299.0, 600.0),
        (300.0, 601.0),
        (500.0, 400.0),
        (f64::NAN, 600.0),
    ] {
        let mut data = source.clone();
        data.temperature_min_k = minimum;
        data.temperature_max_k = maximum;
        assert!(
            append_temperature_series_to_burnup(&inputs, &data, &[], limits, || Ok(())).is_err()
        );
    }
    let mut subdomain = source.clone();
    subdomain.temperature_min_k = 350.0;
    subdomain.temperature_max_k = 550.0;
    assert!(
        append_temperature_series_to_burnup(&inputs, &subdomain, &[], limits, || Ok(())).is_ok()
    );
    subdomain.temperature_min_k = 300.0;
    subdomain.temperature_max_k = 300.0;
    assert!(
        append_temperature_series_to_burnup(&inputs[..1], &subdomain, &[], limits, || Ok(()))
            .is_ok()
    );

    let evaluated = evaluation();
    let mut changed_evaluation = evaluated.clone();
    let mut points = changed_evaluation.points().to_vec();
    let mut values = points[0].values.to_vec();
    values[0].yield_per_fission = 0.6;
    points[0].values = values.into();
    changed_evaluation.points = points.into();
    let mut wrong_revision = cross(600.0, [0.8, 0.0]);
    wrong_revision.metadata.material_revision += 1;
    let full_spectrum = spectrum([1.0, 1.0]);
    let alternatives = [
        collapse(
            &changed_evaluation,
            &cross(600.0, [0.8, 0.0]),
            &full_spectrum,
            &groups,
        ),
        collapse(&evaluated, &wrong_revision, &full_spectrum, &groups),
        // Same identity and integrated group flux; different within-group shape.
        collapse(
            &evaluated,
            &cross(600.0, [0.8, 0.0]),
            &spectrum([0.0, 2.0]),
            &groups,
        ),
    ];
    near(
        alternatives[2].reaction_weighted().flux.values_m2_s()[0],
        yields[0].reaction_weighted().flux.values_m2_s()[0],
    );
    for alternative in alternatives {
        let incompatible = [
            inputs[0],
            FissionYieldTemperatureKnot {
                yields: &alternative,
                prompt_energy_per_group: &prompt[1],
            },
        ];
        assert!(
            append_temperature_series_to_burnup(&incompatible, &source, &[], limits, || Ok(()))
                .is_err()
        );
    }
    // Identical independently allocated source records remain admissible.
    let independently_equal = collapse(
        &evaluation(),
        &cross(600.0, [0.8, 0.0]),
        &spectrum([1.0, 1.0]),
        &groups,
    );
    let equal = [
        inputs[0],
        FissionYieldTemperatureKnot {
            yields: &independently_equal,
            prompt_energy_per_group: &prompt[1],
        },
    ];
    assert!(append_temperature_series_to_burnup(&equal, &source, &[], limits, || Ok(())).is_ok());
}

#[test]
fn original_product_export_and_duplicate_fission_admission_are_retained() {
    let groups = grid(&[1.0, 3.0]);
    let yields = standard(&groups);
    let prompt = std::array::from_fn(|_| vec![energy(1.0, 0.0, 0.0)]);
    let knots = knots(&yields, &prompt);
    let limits = FissionYieldTemperatureLimits::default();
    let source = source_data(&groups);
    let mut missing = source.clone();
    missing.nuclides[2].id = NuclideId::new(6, 12, 0).unwrap();
    assert!(append_temperature_series_to_burnup(&knots, &missing, &[], limits, || Ok(())).is_err());
    let exported =
        append_temperature_series_to_burnup(&knots, &missing, &[product(1)], limits, || Ok(()))
            .unwrap();
    assert!(
        exported.reactions[1..]
            .iter()
            .all(|reaction| reaction.products.len() == 1)
    );
    assert!(
        append_temperature_series_to_burnup(&knots, &source, &[product(1)], limits, || Ok(()))
            .is_err()
    );
    let assembled =
        append_temperature_series_to_burnup(&knots, &source, &[], limits, || Ok(())).unwrap();
    assert!(
        append_temperature_series_to_burnup(&knots, &assembled, &[], limits, || Ok(())).is_err()
    );
    let malformed = [
        FissionYieldTemperatureKnot {
            yields: &yields[0],
            prompt_energy_per_group: &[],
        },
        knots[1],
    ];
    assert!(
        append_temperature_series_to_burnup(&malformed, &source, &[], limits, || Ok(())).is_err()
    );
}

#[test]
fn shared_table_and_global_channel_product_limits_are_enforced() {
    let groups = grid(&[1.0, 3.0]);
    let yields = standard(&groups);
    let prompt = std::array::from_fn(|_| vec![energy(1.0, 0.0, 0.0)]);
    let knots = knots(&yields, &prompt);
    let source = source_data(&groups);
    let limits = FissionYieldTemperatureLimits::default();
    let mut cases = Vec::new();
    let mut limited = limits;
    limited.burnup.max_reactions = 2; // each knot alone fits; their sum does not
    cases.push(limited);
    limited = limits;
    limited.burnup.max_products = 3; // two products per knot, four in the series
    cases.push(limited);
    limited = limits;
    limited.burnup.propagation.maximum_edges = 3;
    cases.push(limited);
    limited = limits;
    limited.cross_sections.max_values = 1;
    cases.push(limited);
    limited = limits;
    limited.cross_sections.max_temperatures = 1;
    cases.push(limited);
    limited = limits;
    limited.burnup.max_neutron_groups = 0;
    cases.push(limited);
    for limits in cases {
        assert!(matches!(
            append_temperature_series_to_burnup(&knots, &source, &[], limits, || Ok(())),
            Err(EndfError::Physics(FissionError::Budget(_)))
        ));
    }
    // Two one-group temperature rows are retained once, despite two channels.
    let exact = FissionYieldTemperatureLimits {
        cross_sections: CrossSectionLimits {
            max_values: 2,
            ..limits.cross_sections
        },
        ..limits
    };
    assert!(append_temperature_series_to_burnup(&knots, &source, &[], exact, || Ok(())).is_ok());
}

#[test]
fn cancellation_never_publishes_part_of_a_temperature_series() {
    let groups = grid(&[1.0, 3.0]);
    let yields = standard(&groups);
    let prompt = std::array::from_fn(|_| vec![energy(1.0, 0.0, 0.0)]);
    let knots = knots(&yields, &prompt);
    let source = source_data(&groups);
    let original = source.clone();
    let mut total = 0;
    append_temperature_series_to_burnup(
        &knots,
        &source,
        &[],
        FissionYieldTemperatureLimits::default(),
        || {
            total += 1;
            Ok(())
        },
    )
    .unwrap();
    for stop in [1, total / 2, total] {
        let mut polls = 0;
        assert_eq!(
            append_temperature_series_to_burnup(
                &knots,
                &source,
                &[],
                FissionYieldTemperatureLimits::default(),
                || {
                    polls += 1;
                    if polls == stop {
                        Err(FissionError::Cancelled)
                    } else {
                        Ok(())
                    }
                }
            ),
            Err(EndfError::Physics(FissionError::Cancelled))
        );
        assert_eq!(source, original);
    }
}
