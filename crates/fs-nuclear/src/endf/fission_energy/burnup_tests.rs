//! Manufactured evaluations and independent inventory/calorimetry equations.
//! These tests establish component ownership, not evaluated-isotope accuracy.

use super::*;
use crate::burnup::{
    BurnupBudget, BurnupInventory, BurnupModel, BurnupNuclide, BurnupProduct, BurnupRate,
    BurnupReaction,
};
use crate::endf::fission_energy::model::{
    FissionEnergyPreparationLimits, PreparedFissionEnergy, ReferenceEnergyTreatment,
};
use crate::endf::fission_energy::{
    EndfFissionEnergyReadLimits, EndfFissionEnergySelection, FissionEnergyRelease,
};
use crate::endf::fission_yields::model::FissionYieldCollapseLimits;
use crate::endf::fission_yields::{
    EndfFissionYieldReadLimits, EndfFissionYieldSelection, FissionYieldEvaluation,
};
use crate::endf::reaction_weighted::ReactionWeightedLimits;
use crate::endf::spectrum::{FluxSpectrum, FluxSpectrumData, FluxSpectrumLimits};
use crate::endf::{EV_J, EndfInterpolation, EndfInterpolationRegion, PointwiseCrossSection};
use crate::irradiation::{GroupFlux, NeutronGroups};
use fs_qty::Frequency;
use fs_thermochem::decay::DecayLawV1;
use std::sync::Arc;

fn parent() -> NuclideId {
    NuclideId::new(10, 20, 1).unwrap()
}

fn product(isomer: u16) -> NuclideId {
    NuclideId::new(4, 8, isomer).unwrap()
}

fn real(value: f64) -> String {
    for precision in (0..=6).rev() {
        let value = format!("{value:.precision$e}");
        let (mantissa, exponent) = value.split_once('e').unwrap();
        let exponent: i32 = exponent.parse().unwrap();
        let field = format!("{mantissa}{exponent:+}");
        if field.len() <= 11 {
            return format!("{field:>11}");
        }
    }
    panic!("manufactured value does not fit ENDF");
}

fn record(fields: [String; 6], mf: i32, mt: i32) -> String {
    format!("{}{:4}{mf:2}{mt:3}{:5}", fields.concat(), 777, 1)
}

fn control(c1: f64, c2: f64, ints: [i32; 4], mf: i32, mt: i32) -> String {
    record(
        [
            real(c1),
            real(c2),
            format!("{:11}", ints[0]),
            format!("{:11}", ints[1]),
            format!("{:11}", ints[2]),
            format!("{:11}", ints[3]),
        ],
        mf,
        mt,
    )
}

fn metadata(sublibrary: i32) -> Vec<String> {
    vec![
        control(10020.0, 20.0, [1, 1, 0, 3], 1, 451),
        control(1.0, 0.0, [1, 1, 0, 6], 1, 451),
        control(1.0, 2e7, [0, 0, sublibrary, 8], 1, 451),
        control(0.0, 0.0, [0, 0, 0, 0], 1, 451),
        control(0.0, 0.0, [0; 4], 1, 0),
    ]
}

fn energy(fragments_ev: f64) -> Arc<FissionEnergyRelease> {
    // Delayed values deliberately differ from the explicit daughter decay law.
    // They remain observations and must never be added at the fission event.
    let components = [
        fragments_ev,
        2.0,
        3.0,
        1.0,
        5.0,
        7.0,
        11.0,
        fragments_ev + 18.0,
        fragments_ev + 29.0,
    ];
    let mut lines = metadata(10);
    lines.push(control(10020.0, 20.0, [0, 0, 0, 0], 1, 458));
    lines.push(control(0.0, 0.0, [0, 0, 18, 9], 1, 458));
    let values: Vec<_> = components.into_iter().flat_map(|v| [v, 0.0]).collect();
    for values in values.chunks_exact(6) {
        lines.push(record(std::array::from_fn(|i| real(values[i])), 1, 458));
    }
    lines.push(control(0.0, 0.0, [0; 4], 1, 0));
    Arc::new(
        FissionEnergyRelease::parse(
            "manufactured prompt energy".into(),
            EndfFissionEnergySelection {
                material_number: 777,
                nuclide: parent(),
                temperature_k: 0.0,
            },
            &lines.join("\n"),
            EndfFissionEnergyReadLimits::default(),
            || Ok(()),
        )
        .unwrap(),
    )
}

fn products() -> FissionYieldEvaluation {
    let mut lines = metadata(11);
    lines.push(control(10020.0, 20.0, [1, 0, 0, 0], 8, 454));
    lines.push(control(1.0, 0.0, [0, 0, 8, 2], 8, 454));
    lines.push(record(
        [4008.0, 0.0, 1.0, 0.0, 4008.0, 1.0].map(real),
        8,
        454,
    ));
    lines.push(record([1.0, 0.0, 0.0, 0.0, 0.0, 0.0].map(real), 8, 454));
    lines.push(control(0.0, 0.0, [0; 4], 8, 0));
    FissionYieldEvaluation::parse(
        "manufactured immediate products".into(),
        EndfFissionYieldSelection {
            material_number: 777,
            nuclide: parent(),
        },
        &lines.join("\n"),
        EndfFissionYieldReadLimits::default(),
        || Ok(()),
    )
    .unwrap()
}

fn grid(edges_ev: &[f64]) -> Arc<NeutronGroups> {
    Arc::new(
        NeutronGroups::new(
            "manufactured prompt grid".into(),
            edges_ev.iter().map(|v| v * EV_J).collect(),
        )
        .unwrap(),
    )
}

fn cross(
    evaluation: &FissionEnergyRelease,
    edges_ev: &[f64],
    sigma: &[f64],
) -> PointwiseCrossSection {
    let mut metadata = evaluation.metadata().clone();
    metadata.temperature_k = 300.0;
    metadata.resonance_flag = 2;
    metadata.derivation_flag = 1;
    PointwiseCrossSection {
        data_id: "manufactured fission reaction".into(),
        metadata,
        reaction: NeutronReaction::Fission,
        energies_j: edges_ev.iter().map(|v| v * EV_J).collect::<Vec<_>>().into(),
        cross_sections_m2: sigma.to_vec().into(),
        interpolation_regions: vec![EndfInterpolationRegion {
            end_point: edges_ev.len(),
            interpolation: EndfInterpolation::Linear,
        }]
        .into(),
        mass_difference_q_j: 0.0,
        reaction_q_j: 0.0,
    }
}

fn spectrum(bounds_ev: [f64; 2], amplitude: f64) -> FluxSpectrum {
    FluxSpectrum::new(
        FluxSpectrumData {
            data_id: "manufactured flat incident spectrum".into(),
            energies_j: bounds_ev.map(|v| v * EV_J).to_vec(),
            flux_density_m2_s_j: vec![amplitude / EV_J; 2],
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
    evaluation: Arc<FissionEnergyRelease>,
    cross: &PointwiseCrossSection,
    spectrum: &FluxSpectrum,
    groups: Arc<NeutronGroups>,
) -> CollapsedFissionEnergy {
    let edges = groups.edges_j();
    PreparedFissionEnergy::new(
        evaluation,
        [edges[0], edges[edges.len() - 1]],
        ReferenceEnergyTreatment::ConstantReference,
        FissionEnergyPreparationLimits::default(),
        || Ok(()),
    )
    .unwrap()
    .collapse(
        cross,
        spectrum,
        groups,
        ReactionWeightedLimits::default(),
        || Ok(()),
    )
    .unwrap()
}

fn routing() -> PromptFissionRouting {
    PromptFissionRouting {
        fragments: FissionEnergyDestination::Deposited,
        prompt_neutrons: FissionEnergyDestination::Radiation,
        prompt_photons: FissionEnergyDestination::Radiation,
    }
}

fn data(groups: Arc<NeutronGroups>) -> BurnupData {
    let node = |id, rate| BurnupNuclide {
        id,
        law: DecayLawV1::from_decay_constant(Frequency::new(rate)).unwrap(),
    };
    BurnupData {
        data_id: "manufactured evaluated-energy burnup".into(),
        groups,
        temperature_min_k: 300.0,
        temperature_max_k: 300.0,
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
    }
}

fn near(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 5e-11 * expected.abs().max(1e-30),
        "{actual:.17e} != {expected:.17e}"
    );
}

#[test]
fn prompt_routing_keeps_delayed_energy_and_incident_kinetic_input_distinct() {
    let evaluation = energy(8.0);
    let cross = cross(&evaluation, &[1.0, 3.0], &[0.2, 0.2]);
    let spectrum = spectrum([1.0, 3.0], 1.0);
    let collapsed = collapse(evaluation, &cross, &spectrum, grid(&[1.0, 2.0, 3.0]));
    let mut split = routing();
    split.prompt_photons = FissionEnergyDestination::Escaped;
    for (group, routed) in collapsed
        .prompt_products(split, || Ok(()))
        .unwrap()
        .iter()
        .enumerate()
    {
        let routed = routed.unwrap();
        near(routed.channels.deposited_j, 8.0 * EV_J);
        near(routed.channels.radiation_j, 2.0 * EV_J);
        near(routed.channels.escaped_j, EV_J);
        near(routed.outgoing_product_energy_j, 11.0 * EV_J);
        near(
            routed.incident_neutron_energy_j,
            (1.5 + group as f64) * EV_J,
        );
        near(routed.prompt_q_j, (9.5 - group as f64) * EV_J);
        assert_eq!(routed.interpolation_error_j, 0.0);
        assert!(routed.roundoff_j < 1e-30);
    }
    let local = PromptFissionRouting {
        fragments: FissionEnergyDestination::Deposited,
        prompt_neutrons: FissionEnergyDestination::Deposited,
        prompt_photons: FissionEnergyDestination::Deposited,
    };
    for group in collapsed
        .prompt_products(local, || Ok(()))
        .unwrap()
        .into_iter()
        .flatten()
    {
        near(group.channels.deposited_j, 11.0 * EV_J);
        assert_eq!(group.channels.radiation_j, 0.0);
        assert_eq!(group.channels.escaped_j, 0.0);
    }
    near(
        collapsed.component_j(FissionEnergyComponent::DelayedPhotons)[0].unwrap(),
        5.0 * EV_J,
    );
    near(
        collapsed.component_j(FissionEnergyComponent::Neutrinos)[0].unwrap(),
        11.0 * EV_J,
    );
}

#[test]
fn evaluated_prompt_products_and_original_decay_heat_are_credited_once_through_shutdown() {
    let evaluation = energy(8.0);
    let cross = cross(&evaluation, &[1.0, 3.0], &[0.2, 0.2]);
    let spectrum = spectrum([1.0, 3.0], 1.0);
    let groups = grid(&[1.0, 2.0, 3.0]);
    let collapsed = collapse(evaluation, &cross, &spectrum, Arc::clone(&groups));
    let yields = products()
        .collapse(
            &cross,
            &spectrum,
            Arc::clone(&groups),
            FissionYieldCollapseLimits::default(),
            || Ok(()),
        )
        .unwrap();
    let original = data(Arc::clone(&groups));
    let limits = BurnupLimits::default();
    let assembled = collapsed
        .append_prompt_products_to_burnup(&yields, &original, routing(), &[], limits, 1e-30, || {
            Ok(())
        })
        .unwrap();
    assert_eq!(original, data(Arc::clone(&groups)));
    assert_eq!(assembled.reactions.len(), 3);
    for (group, reaction) in assembled.reactions[1..].iter().enumerate() {
        let BurnupRate::Fission { cross_sections_m2 } = &reaction.rate else {
            panic!()
        };
        near(cross_sections_m2[group], 0.2);
        assert_eq!(cross_sections_m2[1 - group], 0.0);
        near(reaction.energy_per_event.deposited_j, 8.0 * EV_J);
        near(reaction.energy_per_event.radiation_j, 3.0 * EV_J);
        assert_eq!(reaction.energy_per_event.escaped_j, 0.0);
        assert!(reaction.data_id.contains("manufactured prompt energy"));
    }
    let model = Arc::new(BurnupModel::new(assembled, limits).unwrap());
    let mut inventory = BurnupInventory::new(model, vec![100.0, 0.0, 0.0, 0.0]).unwrap();
    let budget = BurnupBudget {
        limits,
        nuclei_tolerance: 1e-8,
        energy_tolerance_j: 1e-24,
    };
    let flux = GroupFlux::new(Arc::clone(&groups), vec![2.0, 0.0]).unwrap();
    let step = inventory
        .advance(1.25, &flux, 300.0, budget, || Ok(()))
        .unwrap();
    let fissions = 100.0 * -(-0.5_f64).exp_m1();
    let surviving_product = 100.0 * 0.4 * ((-0.5_f64).exp() - (-0.375_f64).exp()) / (0.3 - 0.4);
    let decays = fissions - surviving_product;
    near(step.fissions, fissions);
    near(step.decays, decays);
    near(inventory.populations()[1], surviving_product);
    near(
        step.released.deposited_j,
        (8.0 * fissions + 4.0 * decays) * EV_J,
    );
    near(
        step.released.radiation_j,
        (3.0 * fissions + 3.0 * decays) * EV_J,
    );
    near(step.released.escaped_j, 6.0 * decays * EV_J);
    near(step.fission_energy_j, 11.0 * fissions * EV_J);
    near(step.decay_energy_j, 13.0 * decays * EV_J);
    let shutdown = GroupFlux::new(groups, vec![0.0, 0.0]).unwrap();
    let cooling = inventory
        .advance(2.0, &shutdown, 300.0, budget, || Ok(()))
        .unwrap();
    let later_decays = surviving_product * -(-0.6_f64).exp_m1();
    assert_eq!(cooling.fissions, 0.0);
    near(cooling.decays, later_decays);
    near(cooling.released.deposited_j, 4.0 * later_decays * EV_J);
    near(cooling.released.radiation_j, 3.0 * later_decays * EV_J);
    near(cooling.released.escaped_j, 6.0 * later_decays * EV_J);
}

#[test]
fn bridge_checks_actual_condensation_and_ids_before_combining_prompt_energy_and_products() {
    let evaluation = energy(8.0);
    let original_cross = cross(&evaluation, &[1.0, 3.0], &[0.2, 0.2]);
    let original_spectrum = spectrum([1.0, 3.0], 1.0);
    let groups = grid(&[1.0, 2.0, 3.0]);
    let collapsed = collapse(
        evaluation,
        &original_cross,
        &original_spectrum,
        Arc::clone(&groups),
    );
    let original = data(Arc::clone(&groups));
    for change in 0..6 {
        let mut changed_cross = original_cross.clone();
        let mut changed_spectrum = original_spectrum.clone();
        let mut changed_groups = Arc::clone(&groups);
        match change {
            0 => changed_cross.cross_sections_m2 = vec![0.3, 0.3].into(),
            1 => changed_spectrum = spectrum([1.0, 3.0], 2.0),
            2 => changed_cross.data_id = "different processed source".into(),
            3 => changed_cross.metadata.material_revision += 1,
            4 => {
                changed_groups = Arc::new(
                    NeutronGroups::new(
                        "different named group grid".into(),
                        groups.edges_j().to_vec(),
                    )
                    .unwrap(),
                )
            }
            _ => changed_cross.metadata.temperature_k = 400.0,
        }
        let yields = products()
            .collapse(
                &changed_cross,
                &changed_spectrum,
                changed_groups,
                FissionYieldCollapseLimits::default(),
                || Ok(()),
            )
            .unwrap();
        assert!(
            collapsed
                .append_prompt_products_to_burnup(
                    &yields,
                    &original,
                    routing(),
                    &[],
                    BurnupLimits::default(),
                    1e-30,
                    || Ok(())
                )
                .is_err(),
            "admitted incompatible condensation case {change}"
        );
        assert_eq!(original, data(Arc::clone(&groups)));
    }
}

#[test]
fn equal_group_means_and_reused_ids_do_not_hide_different_incident_energy_shapes() {
    let evaluation = energy(8.0);
    let original_cross = cross(&evaluation, &[0.0, 2.0], &[0.25, 0.25]);
    let shaped_spectrum = |values| {
        FluxSpectrum::new(
            FluxSpectrumData {
                data_id: "same declared source identity".into(),
                energies_j: vec![0.0, 2.0 * EV_J],
                flux_density_m2_s_j: values,
                interpolation_regions: vec![EndfInterpolationRegion {
                    end_point: 2,
                    interpolation: EndfInterpolation::Linear,
                }],
            },
            FluxSpectrumLimits::default(),
        )
        .unwrap()
    };
    let original_spectrum = shaped_spectrum(vec![1.0, 1.0]);
    let groups = grid(&[0.0, 2.0]);
    let collapsed = collapse(
        evaluation,
        &original_cross,
        &original_spectrum,
        Arc::clone(&groups),
    );
    let original = data(Arc::clone(&groups));
    for change in 0..2 {
        let mut changed_cross = original_cross.clone();
        let mut changed_spectrum = original_spectrum.clone();
        if change == 0 {
            changed_cross.cross_sections_m2 = vec![0.0, 0.5].into();
        } else {
            changed_spectrum = shaped_spectrum(vec![0.0, 2.0]);
        }
        let yields = products()
            .collapse(
                &changed_cross,
                &changed_spectrum,
                Arc::clone(&groups),
                FissionYieldCollapseLimits::default(),
                || Ok(()),
            )
            .unwrap();
        let energy = collapsed.reaction_weighted();
        let products = yields.reaction_weighted();
        assert_eq!(energy.cross_sections_m2, products.cross_sections_m2);
        assert_eq!(energy.flux, products.flux);
        assert_eq!(energy.metadata, products.metadata);
        assert_eq!(energy.cross_section_data_id, products.cross_section_data_id);
        assert_eq!(energy.spectrum_data_id, products.spectrum_data_id);
        assert!(
            collapsed
                .append_prompt_products_to_burnup(
                    &yields,
                    &original,
                    routing(),
                    &[],
                    BurnupLimits::default(),
                    1e-30,
                    || Ok(())
                )
                .is_err(),
            "group averages concealed a changed source shape in case {change}"
        );
        assert_eq!(original, data(Arc::clone(&groups)));
    }
}

#[test]
fn zero_reaction_energy_remains_undefined_and_adds_no_fission_channel() {
    let evaluation = energy(8.0);
    let cross = cross(&evaluation, &[0.0, 1.0, 1.0, 2.0], &[0.0, 0.0, 0.2, 0.2]);
    let spectrum = spectrum([0.0, 2.0], 1.0);
    let groups = grid(&[0.0, 1.0, 2.0]);
    let collapsed = collapse(evaluation, &cross, &spectrum, Arc::clone(&groups));
    let prompt = collapsed.prompt_products(routing(), || Ok(())).unwrap();
    assert_eq!(prompt[0], None);
    assert!(prompt[1].is_some());
    let yields = products()
        .collapse(
            &cross,
            &spectrum,
            Arc::clone(&groups),
            FissionYieldCollapseLimits::default(),
            || Ok(()),
        )
        .unwrap();
    let assembled = collapsed
        .append_prompt_products_to_burnup(
            &yields,
            &data(groups),
            routing(),
            &[],
            BurnupLimits::default(),
            1e-30,
            || Ok(()),
        )
        .unwrap();
    assert_eq!(assembled.reactions.len(), 2);
    let BurnupRate::Fission { cross_sections_m2 } = &assembled.reactions[1].rate else {
        panic!()
    };
    assert_eq!(cross_sections_m2[0], 0.0);
    near(cross_sections_m2[1], 0.2);
}

#[test]
fn aggregation_allowance_limits_and_late_cancellation_preserve_original_data() {
    let evaluation = energy(1e30);
    let cross = cross(&evaluation, &[1.0, 3.0], &[0.2, 0.2]);
    let spectrum = spectrum([1.0, 3.0], 1.0);
    let groups = grid(&[1.0, 2.0, 3.0]);
    let collapsed = collapse(evaluation, &cross, &spectrum, Arc::clone(&groups));
    let yields = products()
        .collapse(
            &cross,
            &spectrum,
            Arc::clone(&groups),
            FissionYieldCollapseLimits::default(),
            || Ok(()),
        )
        .unwrap();
    let original = data(Arc::clone(&groups));
    let local = PromptFissionRouting {
        fragments: FissionEnergyDestination::Deposited,
        prompt_neutrons: FissionEnergyDestination::Deposited,
        prompt_photons: FissionEnergyDestination::Deposited,
    };
    let prompt = collapsed.prompt_products(local, || Ok(())).unwrap();
    let error = prompt
        .iter()
        .flatten()
        .map(|p| p.roundoff_j)
        .fold(0.0_f64, f64::max);
    assert!(error > 0.0);
    let limits = BurnupLimits::default();
    assert!(matches!(
        collapsed.append_prompt_products_to_burnup(
            &yields,
            &original,
            local,
            &[],
            limits,
            0.0,
            || Ok(())
        ),
        Err(EndfError::Physics(FissionError::Budget(_)))
    ));
    for allowance in [f64::NAN, -1.0, f64::INFINITY] {
        assert!(
            collapsed
                .append_prompt_products_to_burnup(
                    &yields,
                    &original,
                    local,
                    &[],
                    limits,
                    allowance,
                    || Ok(())
                )
                .is_err()
        );
    }
    assert!(
        collapsed
            .append_prompt_products_to_burnup(
                &yields,
                &original,
                local,
                &[],
                BurnupLimits {
                    max_neutron_groups: 1,
                    ..limits
                },
                error * 2.0,
                || Ok(())
            )
            .is_err()
    );
    let mut polls = 0;
    collapsed
        .append_prompt_products_to_burnup(
            &yields,
            &original,
            local,
            &[],
            limits,
            error * 2.0,
            || {
                polls += 1;
                Ok(())
            },
        )
        .unwrap();
    for cutoff in [1, polls / 2, polls] {
        let mut calls = 0;
        let result = collapsed.append_prompt_products_to_burnup(
            &yields,
            &original,
            local,
            &[],
            limits,
            error * 2.0,
            || {
                calls += 1;
                if calls == cutoff {
                    Err(FissionError::Cancelled)
                } else {
                    Ok(())
                }
            },
        );
        assert!(matches!(
            result,
            Err(EndfError::Physics(FissionError::Cancelled))
        ));
        assert_eq!(original, data(Arc::clone(&groups)));
    }
}
