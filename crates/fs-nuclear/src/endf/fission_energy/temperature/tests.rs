//! Manufactured evaluated laws and independent production/calorimetry checks.
//! These establish temperature interpolation and ownership, not isotope accuracy.

use super::*;
use crate::EnergyChannels;
use crate::burnup::{
    BurnupBudget, BurnupInventory, BurnupLimits, BurnupModel, BurnupNuclide, BurnupProduct,
    BurnupRate, BurnupReaction,
};
use crate::cross_sections::{CrossSectionLimits, NeutronReaction};
use crate::endf::fission_energy::model::FissionEnergyPreparationLimits;
use crate::endf::fission_energy::tests::{self as source, TableFixture};
use crate::endf::fission_yields::model::FissionYieldCollapseLimits;
use crate::endf::fission_yields::{
    EndfFissionYieldReadLimits, EndfFissionYieldSelection, FissionYieldEvaluation,
};
use crate::endf::neutron_yield::{
    EndfNeutronYieldReadLimits, EndfNeutronYieldSelection, FissionNeutronKind, FissionNeutronYield,
};
use crate::endf::reaction_weighted::ReactionWeightedLimits;
use crate::endf::spectrum::{FluxSpectrum, FluxSpectrumData, FluxSpectrumLimits};
use crate::endf::{EV_J, EndfInterpolation, EndfInterpolationRegion, PointwiseCrossSection};
use crate::irradiation::{GroupFlux, NeutronGroups};
use fs_qty::Frequency;
use fs_thermochem::decay::DecayLawV1;

fn parent() -> NuclideId {
    source::selection().nuclide
}

fn product(state: u16) -> NuclideId {
    NuclideId::new(4, 8, state).unwrap()
}

fn energy(
    fragment_coefficients: &[f64],
    uncertainty: f64,
    tables: &[TableFixture],
) -> Arc<FissionEnergyRelease> {
    let fragments = fragment_coefficients[0];
    let mut orders = vec![source::coefficients([
        fragments,
        2.0,
        3.0,
        1.0,
        5.0,
        7.0,
        11.0,
        fragments + 18.0,
        fragments + 29.0,
    ])];
    orders[0][0].uncertainty = uncertainty;
    for &coefficient in &fragment_coefficients[1..] {
        let mut values = [0.0; 9];
        values[0] = coefficient;
        values[7] = coefficient;
        values[8] = coefficient;
        orders.push(source::coefficients(values));
    }
    Arc::new(source::parse(&source::fixture(&orders, tables)).unwrap())
}

fn as_mf8(mut record: String) -> String {
    record.replace_range(70..72, " 8");
    record
}

fn products() -> FissionYieldEvaluation {
    let mut records = vec![
        source::control(10020.0, 19.8, [1, 1, 0, 3], 451),
        source::control(1.0, 0.0, [1, 1, 0, 6], 451),
        source::control(1.0, 2e7, [2, 0, 11, 8], 451),
        source::control(0.0, 0.0, [0; 4], 451),
        source::control(0.0, 0.0, [0; 4], 0),
        as_mf8(source::control(10020.0, 19.8, [2, 0, 0, 0], 454)),
    ];
    for (index, (incident, first, second)) in
        [(0.0, 0.5, 1.5), (2.0, 1.5, 0.5)].into_iter().enumerate()
    {
        records.push(as_mf8(source::control(
            incident,
            0.0,
            [if index == 0 { 1 } else { 2 }, 0, 8, 2],
            454,
        )));
        records.push(as_mf8(source::record(
            [4008.0, 0.0, first, 0.0, 4008.0, 1.0].map(source::real),
            454,
        )));
        records.push(as_mf8(source::record(
            [second, 0.0, 0.0, 0.0, 0.0, 0.0].map(source::real),
            454,
        )));
    }
    records.push(as_mf8(source::control(0.0, 0.0, [0; 4], 0)));
    FissionYieldEvaluation::parse(
        "manufactured-temperature-products".into(),
        EndfFissionYieldSelection {
            material_number: source::selection().material_number,
            nuclide: parent(),
        },
        &records.join("\n"),
        EndfFissionYieldReadLimits::default(),
        || Ok(()),
    )
    .unwrap()
}

fn region() -> Vec<EndfInterpolationRegion> {
    vec![EndfInterpolationRegion {
        end_point: 2,
        interpolation: EndfInterpolation::Linear,
    }]
}

fn cross(
    evaluation: &FissionEnergyRelease,
    temperature: f64,
    cross_sections: [f64; 2],
) -> PointwiseCrossSection {
    let mut metadata = evaluation.metadata().clone();
    metadata.temperature_k = temperature;
    metadata.resonance_flag = 2;
    metadata.derivation_flag = 1;
    PointwiseCrossSection {
        data_id: format!("manufactured fission at {temperature} K").into(),
        metadata,
        reaction: NeutronReaction::Fission,
        energies_j: vec![0.0, 2.0 * EV_J].into(),
        cross_sections_m2: cross_sections.to_vec().into(),
        interpolation_regions: region().into(),
        mass_difference_q_j: 0.0,
        reaction_q_j: 0.0,
    }
}

fn spectrum() -> FluxSpectrum {
    FluxSpectrum::new(
        FluxSpectrumData {
            data_id: "original flat temperature-series spectrum".into(),
            energies_j: vec![0.0, 2.0 * EV_J],
            flux_density_m2_s_j: vec![1.0 / EV_J; 2],
            interpolation_regions: region(),
        },
        FluxSpectrumLimits::default(),
    )
    .unwrap()
}

fn prepare(
    evaluation: Arc<FissionEnergyRelease>,
    reference: ReferenceEnergyTreatment,
) -> PreparedFissionEnergy {
    PreparedFissionEnergy::new(
        evaluation,
        [0.0, 2.0 * EV_J],
        reference,
        FissionEnergyPreparationLimits::default(),
        || Ok(()),
    )
    .unwrap()
}

struct Fixture {
    energy: [CollapsedFissionEnergy; 2],
    yields: [FissionYieldCollapse; 2],
    data: BurnupData,
}

impl Fixture {
    fn knots(&self) -> [FissionEnergyTemperatureKnot<'_>; 2] {
        std::array::from_fn(|index| FissionEnergyTemperatureKnot {
            energy: &self.energy[index],
            yields: &self.yields[index],
        })
    }

    fn recollapse(&self, prepared: &PreparedFissionEnergy, index: usize) -> CollapsedFissionEnergy {
        prepared
            .collapse(
                self.energy[index].source_cross_section(),
                self.energy[index].source_spectrum(),
                Arc::clone(&self.data.groups),
                ReactionWeightedLimits::default(),
                || Ok(()),
            )
            .unwrap()
    }
}

fn fixture_with(prepared: PreparedFissionEnergy) -> Fixture {
    let groups = Arc::new(
        NeutronGroups::new("temperature fission group".into(), vec![0.0, 2.0 * EV_J]).unwrap(),
    );
    let spectrum = spectrum();
    let products = products();
    let [cold, hot] = [(300.0, [0.1, 0.1]), (900.0, [0.0, 0.6])].map(|(temperature, sigma)| {
        let cross = cross(prepared.evaluation(), temperature, sigma);
        let energy = prepared
            .collapse(
                &cross,
                &spectrum,
                Arc::clone(&groups),
                ReactionWeightedLimits::default(),
                || Ok(()),
            )
            .unwrap();
        let yields = products
            .collapse(
                &cross,
                &spectrum,
                Arc::clone(&groups),
                FissionYieldCollapseLimits::default(),
                || Ok(()),
            )
            .unwrap();
        (energy, yields)
    });
    let node = |id, rate| BurnupNuclide {
        id,
        law: DecayLawV1::from_decay_constant(Frequency::new(rate)).unwrap(),
    };
    Fixture {
        energy: [cold.0, hot.0],
        yields: [cold.1, hot.1],
        data: BurnupData {
            data_id: "manufactured evaluated temperature feedback".into(),
            groups,
            temperature_min_k: 300.0,
            temperature_max_k: 900.0,
            nuclides: vec![
                node(parent(), 0.0),
                node(product(0), 0.3),
                node(product(1), 0.0),
                node(NuclideId::new(5, 8, 0).unwrap(), 0.0),
            ],
            reactions: vec![BurnupReaction {
                data_id: "original product decay".into(),
                parent: 1,
                rate: BurnupRate::Decay {
                    branch_fraction: 1.0,
                },
                products: vec![BurnupProduct {
                    nuclide: 3,
                    yield_per_event: 1.0,
                }],
                energy_per_event: EnergyChannels {
                    deposited_j: 4.0 * EV_J,
                    radiation_j: 3.0 * EV_J,
                    escaped_j: 6.0 * EV_J,
                },
            }],
        },
    }
}

fn fixture() -> Fixture {
    fixture_with(prepare(
        energy(&[4.0, 2.0], 0.0, &[]),
        ReferenceEnergyTreatment::Reject,
    ))
}

fn limits() -> FissionYieldTemperatureLimits {
    FissionYieldTemperatureLimits {
        burnup: BurnupLimits::default(),
        cross_sections: CrossSectionLimits::default(),
    }
}

fn routing() -> PromptFissionRouting {
    use super::super::burnup::FissionEnergyDestination;
    PromptFissionRouting {
        fragments: FissionEnergyDestination::Deposited,
        prompt_neutrons: FissionEnergyDestination::Radiation,
        prompt_photons: FissionEnergyDestination::Radiation,
    }
}

fn assemble(
    fixture: &Fixture,
    checkpoint: impl FnMut() -> Result<(), FissionError>,
) -> Result<BurnupData, EndfError> {
    append_prompt_temperature_series_to_burnup(
        &fixture.knots(),
        &fixture.data,
        routing(),
        &[],
        limits(),
        1e-30,
        checkpoint,
    )
}

/// Actual parse/collapse/series assembly shared by the original thermal tests.
/// Nuclides: parent, decaying product, stable product, decay daughter. At
/// x=(T-300)/600, the one-group moments are sigma=.1+.2x,
/// sigma*Y0=.1+.25x, sigma*EFR=(.6+1.4x)eV, sigma*Eprompt=(.9+2x)eV.
#[cfg(feature = "thermal")]
pub(crate) fn assembled_fixture() -> BurnupData {
    assemble(&fixture(), || Ok(())).unwrap()
}

fn near(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 8e-10 * expected.abs().max(1e-30),
        "{actual:.17e} != {expected:.17e}"
    );
}

#[test]
fn processed_temperature_rates_preserve_production_and_prompt_heat_at_knots_and_interiors() {
    let fixture = fixture();
    let source = fixture.data.clone();
    let assembled = assemble(&fixture, || Ok(())).unwrap();
    assert_eq!(fixture.data, source);
    assert_eq!(assembled.reactions[0], source.reactions[0]);
    assert!(assembled.reactions[1..].iter().all(|reaction| {
        reaction.data_id.contains("manufactured-energy-release")
            && matches!(reaction.rate, BurnupRate::TabulatedFission { .. })
    }));
    for (index, (sigma, product_yield, fragment)) in [(0.1, 1.0, 6.0), (0.3, 7.0 / 6.0, 20.0 / 3.0)]
        .into_iter()
        .enumerate()
    {
        near(
            fixture.energy[index].reaction_weighted().cross_sections_m2[0],
            sigma,
        );
        near(
            fixture.yields[index].reaction_weighted().yields[0].values[0].unwrap(),
            product_yield,
        );
        near(
            fixture.energy[index].component_j(FissionEnergyComponent::Fragments)[0].unwrap(),
            fragment * EV_J,
        );
    }
    let model = Arc::new(BurnupModel::new(assembled, limits().burnup).unwrap());
    let flux = GroupFlux::new(Arc::clone(&source.groups), vec![2.0]).unwrap();
    let budget = BurnupBudget {
        limits: limits().burnup,
        nuclei_tolerance: 1e-8,
        energy_tolerance_j: 1e-24,
    };
    for temperature in [300.0_f64, 450.0, 600.0, 750.0, 900.0] {
        let x = (temperature - 300.0) / 600.0;
        let sigma = 0.1 + 0.2 * x;
        let rate = 2.0 * sigma;
        let yield0 = (0.1 + 0.25 * x) / sigma;
        let fragment = (0.6 + 1.4 * x) / sigma;
        let duration = 1.25_f64;
        let fissions = 100.0 * -(-rate * duration).exp_m1();
        let survivor = if (rate - 0.3).abs() < 1e-12 {
            100.0 * yield0 * rate * duration * (-rate * duration).exp()
        } else {
            100.0 * yield0 * rate * ((-rate * duration).exp() - (-0.3 * duration).exp())
                / (0.3 - rate)
        };
        let decays = yield0 * fissions - survivor;
        let mut inventory =
            BurnupInventory::new(Arc::clone(&model), vec![100.0, 0.0, 0.0, 0.0]).unwrap();
        let step = inventory
            .advance(duration, &flux, temperature, budget, || Ok(()))
            .unwrap();
        near(step.fissions, fissions);
        near(step.decays, decays);
        near(inventory.populations()[1], survivor);
        near(inventory.populations()[2], (2.0 - yield0) * fissions);
        near(inventory.populations()[3], decays);
        near(
            step.released.deposited_j,
            (fragment * fissions + 4.0 * decays) * EV_J,
        );
        near(
            step.released.radiation_j,
            (3.0 * fissions + 3.0 * decays) * EV_J,
        );
        near(step.released.escaped_j, 6.0 * decays * EV_J);
        near(step.fission_energy_j, (fragment + 3.0) * fissions * EV_J);
        near(step.decay_energy_j, 13.0 * decays * EV_J);

        let shutdown = GroupFlux::new(Arc::clone(&source.groups), vec![0.0]).unwrap();
        let cooling = inventory
            .advance(2.0, &shutdown, 725.0, budget, || Ok(()))
            .unwrap();
        let later_decays = survivor * -(-0.6_f64).exp_m1();
        assert_eq!(cooling.fissions, 0.0);
        near(cooling.decays, later_decays);
        near(cooling.released.deposited_j, 4.0 * later_decays * EV_J);
        near(cooling.released.radiation_j, 3.0 * later_decays * EV_J);
        near(cooling.released.escaped_j, 6.0 * later_decays * EV_J);
    }
}

#[test]
fn original_evaluation_preparation_and_each_knot_source_must_match() {
    let fixture = fixture();
    let unchanged = fixture.data.clone();
    for change in 0..5 {
        let changed = match change {
            0 => prepare(
                energy(&[5.0, 2.0], 0.0, &[]),
                ReferenceEnergyTreatment::Reject,
            ),
            1 => prepare(
                energy(&[4.0, 2.0], 0.1, &[]),
                ReferenceEnergyTreatment::Reject,
            ),
            2 => prepare(
                energy(
                    &[4.0],
                    0.0,
                    &[TableFixture {
                        component: FissionEnergyComponent::Fragments,
                        derivation: 2,
                        points: vec![(0.0, 4.0), (2.0, 8.0)],
                        regions: vec![(2, 2)],
                    }],
                ),
                ReferenceEnergyTreatment::ConstantReference,
            ),
            3 => PreparedFissionEnergy::new(
                Arc::clone(fixture.energy[0].evaluation()),
                [0.0, 3.0 * EV_J],
                ReferenceEnergyTreatment::Reject,
                FissionEnergyPreparationLimits::default(),
                || Ok(()),
            )
            .unwrap(),
            _ => fixture.energy[0].source_preparation().clone(),
        };
        let hot = if change == 4 {
            // Same processed identity and group-average sigma, different shape.
            let mut altered = fixture.energy[1].source_cross_section().clone();
            altered.cross_sections_m2 = vec![0.6, 0.0].into();
            changed
                .collapse(
                    &altered,
                    fixture.energy[1].source_spectrum(),
                    Arc::clone(&fixture.data.groups),
                    ReactionWeightedLimits::default(),
                    || Ok(()),
                )
                .unwrap()
        } else {
            fixture.recollapse(&changed, 1)
        };
        let mut knots = fixture.knots();
        knots[1].energy = &hot;
        assert!(
            append_prompt_temperature_series_to_burnup(
                &knots,
                &fixture.data,
                routing(),
                &[],
                limits(),
                1e-30,
                || Ok(())
            )
            .is_err(),
            "admitted incompatible temperature energy source {change}"
        );
        assert_eq!(fixture.data, unchanged);
    }
    // Separate equal evaluations/preparations are valid; Arc identity is only
    // a fast path and never a substitute for the original numerical contents.
    let repeated = prepare(
        energy(&[4.0, 2.0], 0.0, &[]),
        ReferenceEnergyTreatment::Reject,
    );
    let hot = fixture.recollapse(&repeated, 1);
    let mut knots = fixture.knots();
    knots[1].energy = &hot;
    assert!(
        append_prompt_temperature_series_to_burnup(
            &knots,
            &fixture.data,
            routing(),
            &[],
            limits(),
            1e-30,
            || Ok(())
        )
        .is_ok()
    );
}

fn neutron_yield(outside_value: f64) -> Arc<FissionNeutronYield> {
    let records = vec![
        source::control(10020.0, 19.8, [1, 1, 0, 3], 451),
        source::control(1.0, 0.0, [1, 1, 0, 6], 451),
        source::control(1.0, 2e7, [2, 0, 10, 8], 451),
        source::control(0.0, 0.0, [0; 4], 451),
        source::control(0.0, 0.0, [0; 4], 0),
        source::control(10020.0, 19.8, [0, 2, 0, 0], 456),
        source::control(0.0, 0.0, [0, 0, 1, 3], 456),
        source::record([3, 2, 0, 0, 0, 0].map(|value| format!("{value:11}")), 456),
        source::record(
            [0.0, 2.4, 2.0, 2.4, 3.0, outside_value].map(source::real),
            456,
        ),
        source::control(0.0, 0.0, [0; 4], 0),
    ];
    Arc::new(
        FissionNeutronYield::parse(
            "same declared prompt-neutron source".into(),
            EndfNeutronYieldSelection {
                material_number: source::selection().material_number,
                nuclide: parent(),
                temperature_k: 0.0,
                kind: FissionNeutronKind::Prompt,
            },
            &records.join("\n"),
            EndfNeutronYieldReadLimits::default(),
            || Ok(()),
        )
        .unwrap(),
    )
}

#[test]
fn continuation_choices_and_full_sher_beck_sources_cannot_be_temperature_dependence() {
    let evaluated = energy(&[4.0], 0.0, &[]);
    let constant = prepare(
        Arc::clone(&evaluated),
        ReferenceEnergyTreatment::ConstantReference,
    );
    let sher_beck = |outside_value| {
        prepare(
            Arc::clone(&evaluated),
            ReferenceEnergyTreatment::SherBeck2023 {
                neutron_yield: neutron_yield(outside_value),
                reference_nu_at_zero: 2.4,
            },
        )
    };
    let reference = sher_beck(2.4);
    let fixture = fixture_with(reference);
    for prepared in [constant, sher_beck(2.5)] {
        let hot = fixture.recollapse(&prepared, 1);
        let mut knots = fixture.knots();
        knots[1].energy = &hot;
        assert!(
            append_prompt_temperature_series_to_burnup(
                &knots,
                &fixture.data,
                routing(),
                &[],
                limits(),
                1e-30,
                || Ok(())
            )
            .is_err()
        );
    }
    let same = sher_beck(2.4);
    let hot = fixture.recollapse(&same, 1);
    let mut knots = fixture.knots();
    knots[1].energy = &hot;
    assert!(
        append_prompt_temperature_series_to_burnup(
            &knots,
            &fixture.data,
            routing(),
            &[],
            limits(),
            1e-30,
            || Ok(())
        )
        .is_ok()
    );
}

#[test]
fn interpolation_approximation_cannot_change_between_temperature_knots() {
    let evaluated = energy(&[4.0, 2.0, 1.0], 0.0, &[]);
    let prepare_with = |absolute_tolerance_j| {
        PreparedFissionEnergy::new(
            Arc::clone(&evaluated),
            [0.0, 2.0 * EV_J],
            ReferenceEnergyTreatment::Reject,
            FissionEnergyPreparationLimits {
                absolute_tolerance_j,
                relative_tolerance: 0.0,
                ..FissionEnergyPreparationLimits::default()
            },
            || Ok(()),
        )
        .unwrap()
    };
    let coarse = prepare_with(EV_J);
    let fine = prepare_with(0.01 * EV_J);
    let fixture = fixture_with(coarse);
    let hot = fixture.recollapse(&fine, 1);
    let mut knots = fixture.knots();
    knots[1].energy = &hot;
    assert!(
        append_prompt_temperature_series_to_burnup(
            &knots,
            &fixture.data,
            routing(),
            &[],
            limits(),
            1e-30,
            || Ok(())
        )
        .is_err()
    );
}

#[test]
fn temperature_assembly_limits_allowances_and_cancellation_preserve_source_data() {
    let fixture = fixture();
    let unchanged = fixture.data.clone();
    assert!(
        append_prompt_temperature_series_to_burnup(
            &[],
            &fixture.data,
            routing(),
            &[],
            limits(),
            1e-30,
            || Ok(())
        )
        .is_err()
    );
    for change in 0..3 {
        let mut restricted = limits();
        match change {
            0 => restricted.cross_sections.max_temperatures = 1,
            1 => restricted.cross_sections.max_values = 1,
            _ => restricted.burnup.max_neutron_groups = 0,
        }
        assert!(matches!(
            append_prompt_temperature_series_to_burnup(
                &fixture.knots(),
                &fixture.data,
                routing(),
                &[],
                restricted,
                1e-30,
                || Ok(())
            ),
            Err(EndfError::Physics(FissionError::Budget(_)))
        ));
    }
    for allowance in [-1.0, f64::NAN, f64::INFINITY] {
        assert!(
            append_prompt_temperature_series_to_burnup(
                &fixture.knots(),
                &fixture.data,
                routing(),
                &[],
                limits(),
                allowance,
                || Ok(())
            )
            .is_err()
        );
    }
    let mut polls = 0;
    let assembled = assemble(&fixture, || {
        polls += 1;
        Ok(())
    })
    .unwrap();
    assert!(polls > 20);
    for stop in [1, polls / 2, polls] {
        let mut current = 0;
        assert!(matches!(
            assemble(&fixture, || {
                current += 1;
                if current == stop {
                    Err(FissionError::Cancelled)
                } else {
                    Ok(())
                }
            }),
            Err(EndfError::Physics(FissionError::Cancelled))
        ));
        assert_eq!(fixture.data, unchanged);
    }
    let model = Arc::new(BurnupModel::new(assembled, limits().burnup).unwrap());
    let mut inventory = BurnupInventory::new(model, vec![100.0, 0.0, 0.0, 0.0]).unwrap();
    let original = inventory.clone();
    let flux = GroupFlux::new(Arc::clone(&fixture.data.groups), vec![2.0]).unwrap();
    let budget = BurnupBudget {
        limits: limits().burnup,
        nuclei_tolerance: 1e-8,
        energy_tolerance_j: 1e-24,
    };
    for outside in [299.0, 901.0] {
        assert!(
            inventory
                .advance(1.0, &flux, outside, budget, || Ok(()))
                .is_err()
        );
        assert_eq!(inventory, original);
    }
}
