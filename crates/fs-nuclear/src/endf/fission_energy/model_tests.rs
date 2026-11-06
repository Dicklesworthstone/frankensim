//! Independent component-law and reaction-integral checks with manufactured data.

use std::sync::Arc;

use super::model::{
    FissionEnergyPreparationLimits, PreparedFissionEnergy, ReferenceEnergyTreatment,
};
use super::tests::{
    TableFixture, coefficients, control, fixture, near, parse, real, record, selection,
};
use super::{FissionEnergyComponent as Component, FissionEnergyRelease};
use crate::FissionError;
use crate::cross_sections::NeutronReaction;
use crate::endf::neutron_yield::{
    EndfNeutronYieldReadLimits, EndfNeutronYieldSelection, FissionNeutronKind, FissionNeutronYield,
};
use crate::endf::reaction_weighted::ReactionWeightedLimits;
use crate::endf::spectrum::{FluxSpectrum, FluxSpectrumData, FluxSpectrumLimits};
use crate::endf::{
    EV_J, EndfError, EndfInterpolation, EndfInterpolationRegion, PointwiseCrossSection,
};
use crate::irradiation::NeutronGroups;

fn domain() -> [f64; 2] {
    [0.0, EV_J]
}

fn evaluated_reference() -> Arc<FissionEnergyRelease> {
    Arc::new(parse(&fixture(&[coefficients([100.0; 9])], &[])).unwrap())
}

fn nu(
    points: &[(f64, f64)],
    regions: &[(usize, i32)],
    kind: FissionNeutronKind,
    library_version: i32,
) -> Arc<FissionNeutronYield> {
    let mt = kind.mt();
    let mut lines = vec![
        control(10020.0, 19.8, [1, 1, 0, 3], 451),
        control(1.0, 0.0, [1, 1, 0, 6], 451),
        control(1.0, 2e7, [2, 0, 10, library_version], 451),
        control(0.0, 0.0, [0; 4], 451),
        control(0.0, 0.0, [0; 4], 0),
        control(10020.0, 19.8, [0, 2, 0, 0], mt),
    ];
    if kind == FissionNeutronKind::Delayed {
        lines.push(control(0.0, 0.0, [0, 0, 1, 0], mt));
        lines.push(record(
            [
                real(1.0),
                real(0.0),
                real(0.0),
                real(0.0),
                real(0.0),
                real(0.0),
            ],
            mt,
        ));
    }
    lines.push(control(
        0.0,
        0.0,
        [0, 0, regions.len() as i32, points.len() as i32],
        mt,
    ));
    for chunk in regions.chunks(3) {
        let mut fields = std::array::from_fn(|_| format!("{:11}", 0));
        for (index, &(end, law)) in chunk.iter().enumerate() {
            fields[index * 2] = format!("{end:11}");
            fields[index * 2 + 1] = format!("{law:11}");
        }
        lines.push(record(fields, mt));
    }
    for chunk in points.chunks(3) {
        let mut fields = std::array::from_fn(|_| real(0.0));
        for (index, &(energy_ev, multiplicity)) in chunk.iter().enumerate() {
            fields[index * 2] = real(energy_ev);
            fields[index * 2 + 1] = real(multiplicity);
        }
        lines.push(record(fields, mt));
    }
    lines.push(control(0.0, 0.0, [0; 4], 0));
    let selected = selection();
    Arc::new(
        FissionNeutronYield::parse(
            "manufactured-matching-nu".into(),
            EndfNeutronYieldSelection {
                material_number: selected.material_number,
                nuclide: selected.nuclide,
                temperature_k: selected.temperature_k,
                kind,
            },
            &(lines.join("\n") + "\n"),
            EndfNeutronYieldReadLimits::default(),
            || Ok(()),
        )
        .unwrap(),
    )
}

fn region(points: usize) -> Arc<[EndfInterpolationRegion]> {
    vec![EndfInterpolationRegion {
        end_point: points,
        interpolation: EndfInterpolation::Linear,
    }]
    .into()
}

fn cross(
    evaluation: &FissionEnergyRelease,
    energies_ev: &[f64],
    values_m2: &[f64],
) -> PointwiseCrossSection {
    let mut metadata = evaluation.metadata().clone();
    metadata.temperature_k = 300.0;
    metadata.resonance_flag = 2;
    metadata.derivation_flag = 1;
    PointwiseCrossSection {
        data_id: "manufactured-processed-fission".into(),
        metadata,
        reaction: NeutronReaction::Fission,
        energies_j: energies_ev
            .iter()
            .map(|energy| energy * EV_J)
            .collect::<Vec<_>>()
            .into(),
        cross_sections_m2: values_m2.to_vec().into(),
        interpolation_regions: region(energies_ev.len()),
        mass_difference_q_j: 0.0,
        reaction_q_j: 0.0,
    }
}

fn spectrum(values_per_ev: [f64; 2]) -> FluxSpectrum {
    FluxSpectrum::new(
        FluxSpectrumData {
            data_id: "manufactured-spectrum".into(),
            energies_j: domain().to_vec(),
            flux_density_m2_s_j: values_per_ev
                .into_iter()
                .map(|value| value / EV_J)
                .collect(),
            interpolation_regions: region(2).to_vec(),
        },
        FluxSpectrumLimits::default(),
    )
    .unwrap()
}

fn groups(edges_ev: &[f64]) -> Arc<NeutronGroups> {
    Arc::new(
        NeutronGroups::new(
            "manufactured-groups".into(),
            edges_ev.iter().map(|energy| energy * EV_J).collect(),
        )
        .unwrap(),
    )
}

#[test]
fn g0_reference_continuation_is_explicit_and_tables_need_no_fallback() {
    let evaluation = evaluated_reference();
    assert!(
        PreparedFissionEnergy::new(
            Arc::clone(&evaluation),
            domain(),
            ReferenceEnergyTreatment::Reject,
            FissionEnergyPreparationLimits::default(),
            || Ok(())
        )
        .is_err()
    );
    let prepared = PreparedFissionEnergy::new(
        Arc::clone(&evaluation),
        domain(),
        ReferenceEnergyTreatment::ConstantReference,
        FissionEnergyPreparationLimits::default(),
        || Ok(()),
    )
    .unwrap();
    assert!(Arc::ptr_eq(prepared.evaluation(), &evaluation));
    assert_eq!(prepared.domain_j(), domain());
    assert!(matches!(
        prepared.reference_treatment(),
        ReferenceEnergyTreatment::ConstantReference
    ));
    for component in Component::ALL {
        near(
            prepared.value_at_energy_j(component, 0.375 * EV_J).unwrap(),
            100.0 * EV_J,
        );
        assert_eq!(prepared.interpolation_error_j(component), 0.0);
    }
    let tables: Vec<_> = Component::ALL
        .into_iter()
        .map(|component| TableFixture {
            component,
            derivation: 2,
            points: vec![(0.0, 3.0), (1.0, 5.0)],
            regions: vec![(2, 2)],
        })
        .collect();
    let tabulated = Arc::new(parse(&fixture(&[coefficients([100.0; 9])], &tables)).unwrap());
    let prepared = PreparedFissionEnergy::new(
        tabulated,
        domain(),
        ReferenceEnergyTreatment::Reject,
        FissionEnergyPreparationLimits::default(),
        || Ok(()),
    )
    .unwrap();
    near(
        prepared
            .value_at_energy_j(Component::Fragments, 0.5 * EV_J)
            .unwrap(),
        4.0 * EV_J,
    );
}

#[test]
fn g1_sher_beck_uses_explicit_nu_zero_and_2023_neutrino_slope_with_table_precedence() {
    let reference = coefficients([160e6, 5e6, 3e4, 7e6, 6e6, 7e6, 9e6, 185e6, 194e6]);
    let photon_table = TableFixture {
        component: Component::PromptPhotons,
        derivation: 2,
        points: vec![(1e5, 4e6), (2e6, 8e6)],
        regions: vec![(2, 2)],
    };
    let evaluation = Arc::new(parse(&fixture(&[reference], &[photon_table])).unwrap());
    let yield_ = nu(
        &[(1e5, 2.45), (2e6, 2.60)],
        &[(2, 2)],
        FissionNeutronKind::Prompt,
        8,
    );
    // The original nu curve has no value at zero. Its explicit reference is
    // supplied independently, without inventing an extrapolated endpoint.
    assert!(yield_.curve().value_at_energy(0.0).is_err());
    let prepared = PreparedFissionEnergy::new(
        evaluation,
        [1e5 * EV_J, 2e6 * EV_J],
        ReferenceEnergyTreatment::SherBeck2023 {
            neutron_yield: yield_,
            reference_nu_at_zero: 2.4,
        },
        FissionEnergyPreparationLimits::default(),
        || Ok(()),
    )
    .unwrap();
    for (component, expected_ev) in [
        (Component::Fragments, 160e6),
        (Component::PromptNeutrons, 6e6),
        (Component::DelayedNeutrons, 3e4),
        (Component::PromptPhotons, 8e6),
        (Component::DelayedPhotons, 5.85e6),
        (Component::Betas, 6.85e6),
        (Component::Neutrinos, 8.8e6),
        (Component::Recoverable, 185.7e6),
        (Component::Total, 194.5e6),
    ] {
        near(
            prepared.value_at_energy_j(component, 2e6 * EV_J).unwrap(),
            expected_ev * EV_J,
        );
        assert_eq!(prepared.interpolation_error_j(component), 0.0);
    }
}

#[test]
fn g0_domain_and_evaluation_mismatch_do_not_infer_missing_physics() {
    let evaluation = Arc::new(parse(&fixture(&[coefficients([100e6; 9])], &[])).unwrap());
    let points = [(1e5, 2.4), (2e6, 2.6)];
    let make = |yield_, reference_nu_at_zero, domain| {
        PreparedFissionEnergy::new(
            Arc::clone(&evaluation),
            domain,
            ReferenceEnergyTreatment::SherBeck2023 {
                neutron_yield: yield_,
                reference_nu_at_zero,
            },
            FissionEnergyPreparationLimits::default(),
            || Ok(()),
        )
    };
    assert!(
        make(
            nu(&points, &[(2, 2)], FissionNeutronKind::Prompt, 8),
            2.4,
            [0.0, 2e6 * EV_J]
        )
        .is_err()
    );
    assert!(
        make(
            nu(&points, &[(2, 2)], FissionNeutronKind::Prompt, 9),
            2.4,
            [1e5 * EV_J, 2e6 * EV_J]
        )
        .is_err()
    );
    assert!(
        make(
            nu(&points, &[(2, 2)], FissionNeutronKind::Delayed, 8),
            2.4,
            [1e5 * EV_J, 2e6 * EV_J]
        )
        .is_err()
    );
    for invalid_reference in [-1.0, f64::NAN, f64::INFINITY] {
        assert!(
            make(
                nu(&points, &[(2, 2)], FissionNeutronKind::Total, 8),
                invalid_reference,
                [1e5 * EV_J, 2e6 * EV_J]
            )
            .is_err()
        );
    }
    let prepared = PreparedFissionEnergy::new(
        Arc::clone(&evaluation),
        domain(),
        ReferenceEnergyTreatment::ConstantReference,
        FissionEnergyPreparationLimits::default(),
        || Ok(()),
    )
    .unwrap();
    for outside in [-EV_J, 2.0 * EV_J, f64::NAN] {
        assert!(
            prepared
                .value_at_energy_j(Component::Fragments, outside)
                .is_err()
        );
    }
    let mut incompatible = cross(&evaluation, &[0.0, 1.0], &[1e-28, 1e-28]);
    incompatible.metadata.material_revision += 1;
    assert!(
        prepared
            .collapse(
                &incompatible,
                &spectrum([1.0, 1.0]),
                groups(&[0.0, 1.0]),
                ReactionWeightedLimits::default(),
                || Ok(())
            )
            .is_err()
    );
}

#[test]
fn g1_sher_beck_retains_histogram_boundaries_and_paired_nu_jumps() {
    let evaluation = Arc::new(parse(&fixture(&[coefficients([100e6; 9])], &[])).unwrap());
    for (points, regions, checks) in [
        (
            vec![(0.0, 2.4), (1e6, 2.6), (2e6, 2.8)],
            vec![(3, 1)],
            vec![(0.5e6, 2.4), (1e6, 2.6), (1.5e6, 2.6), (2e6, 2.8)],
        ),
        (
            vec![(0.0, 2.4), (1e6, 2.5), (1e6, 2.7), (2e6, 2.9)],
            vec![(4, 2)],
            vec![(0.5e6, 2.45), (1e6, 2.7), (1.5e6, 2.8), (2e6, 2.9)],
        ),
    ] {
        let prepared = PreparedFissionEnergy::new(
            Arc::clone(&evaluation),
            [0.0, 2e6 * EV_J],
            ReferenceEnergyTreatment::SherBeck2023 {
                neutron_yield: nu(&points, &regions, FissionNeutronKind::Total, 8),
                reference_nu_at_zero: 2.4,
            },
            FissionEnergyPreparationLimits::default(),
            || Ok(()),
        )
        .unwrap();
        for (energy_ev, multiplicity) in checks {
            let expected_ev = 100e6 + 1.307 * energy_ev - 8.07e6 * (multiplicity - 2.4);
            near(
                prepared
                    .value_at_energy_j(Component::PromptNeutrons, energy_ev * EV_J)
                    .unwrap(),
                expected_ev * EV_J,
            );
        }
    }
}

#[test]
fn g1_positive_polynomial_refinement_bounds_independent_values_and_group_means() {
    let evaluation = Arc::new(
        parse(&fixture(
            &[
                coefficients([2.0; 9]),
                coefficients([3.0; 9]),
                coefficients([4.0; 9]),
            ],
            &[],
        ))
        .unwrap(),
    );
    let xs = cross(&evaluation, &[0.0, 1.0], &[0.0, 1e-28]);
    let mut bounds = Vec::new();
    let mut errors = Vec::new();
    for tolerance_ev in [0.2, 0.002] {
        let prepared = PreparedFissionEnergy::new(
            Arc::clone(&evaluation),
            domain(),
            ReferenceEnergyTreatment::Reject,
            FissionEnergyPreparationLimits {
                absolute_tolerance_j: tolerance_ev * EV_J,
                relative_tolerance: 0.0,
                ..FissionEnergyPreparationLimits::default()
            },
            || Ok(()),
        )
        .unwrap();
        let bound = prepared.interpolation_error_j(Component::Fragments);
        assert!(bound > 0.0 && bound <= tolerance_ev * EV_J);
        for energy_ev in [0.13, 0.37, 0.61, 0.93] {
            let exact = (2.0 + 3.0 * energy_ev + 4.0 * energy_ev * energy_ev) * EV_J;
            let actual = prepared
                .value_at_energy_j(Component::Fragments, energy_ev * EV_J)
                .unwrap();
            assert!((actual - exact).abs() <= bound + 1e-12 * EV_J);
        }
        let collapsed = prepared
            .collapse(
                &xs,
                &spectrum([0.0, 1.0]),
                groups(&[0.0, 1.0]),
                ReactionWeightedLimits::default(),
                || Ok(()),
            )
            .unwrap();
        // Weight E^2: <E>=3/4 and <E^2>=3/5, so <2+3E+4E^2>=6.65 eV.
        let error = (collapsed.component_j(Component::Fragments)[0].unwrap() - 6.65 * EV_J).abs();
        assert!(error <= bound + 1e-12 * EV_J);
        assert_eq!(collapsed.interpolation_error_j(Component::Fragments), bound);
        bounds.push(bound);
        errors.push(error);
    }
    assert!(bounds[1] < bounds[0] / 10.0, "bounds: {bounds:?}");
    assert!(errors[1] < errors[0] / 10.0, "mean errors: {errors:?}");
}

#[test]
fn g1_reaction_weighted_component_and_incident_energies_use_the_full_integrals() {
    let evaluation = Arc::new(
        parse(&fixture(
            &[coefficients([2.0; 9]), coefficients([4.0; 9])],
            &[],
        ))
        .unwrap(),
    );
    let prepared = PreparedFissionEnergy::new(
        Arc::clone(&evaluation),
        domain(),
        ReferenceEnergyTreatment::Reject,
        FissionEnergyPreparationLimits::default(),
        || Ok(()),
    )
    .unwrap();
    let collapsed = prepared
        .collapse(
            &cross(&evaluation, &[0.0, 1.0], &[0.0, 1e-28]),
            &spectrum([0.0, 1.0]),
            groups(&[0.0, 1.0]),
            ReactionWeightedLimits::default(),
            || Ok(()),
        )
        .unwrap();
    // sigma=E barn, phi=E per eV, q=2+4E eV. The fission-weighted
    // q is 2+4*(3/4)=5 eV, not its midpoint or flux-only weighted value.
    for component in Component::ALL {
        near(collapsed.component_j(component)[0].unwrap(), 5.0 * EV_J);
        assert_eq!(collapsed.interpolation_error_j(component), 0.0);
    }
    near(collapsed.incident_energy_j()[0].unwrap(), 0.75 * EV_J);
    let weighted = collapsed.reaction_weighted();
    let events_per_parent_s = weighted.cross_sections_m2[0] * weighted.flux.values_m2_s()[0];
    assert!((events_per_parent_s - 1e-28 / 3.0).abs() < 1e-40);
    assert!(
        (events_per_parent_s * collapsed.component_j(Component::Fragments)[0].unwrap()
            - (5.0 / 3.0) * EV_J * 1e-28)
            .abs()
            < 1e-58
    );
    assert_eq!(collapsed.evaluation().metadata().temperature_k, 0.0);
    assert_eq!(weighted.metadata.temperature_k, 300.0);
}

#[test]
fn g0_zero_fission_groups_have_no_invented_component_or_incident_mean() {
    let evaluation = evaluated_reference();
    let prepared = PreparedFissionEnergy::new(
        Arc::clone(&evaluation),
        domain(),
        ReferenceEnergyTreatment::ConstantReference,
        FissionEnergyPreparationLimits::default(),
        || Ok(()),
    )
    .unwrap();
    let collapsed = prepared
        .collapse(
            &cross(
                &evaluation,
                &[0.0, 0.5, 0.5, 1.0],
                &[0.0, 0.0, 1e-28, 1e-28],
            ),
            &spectrum([1.0, 1.0]),
            groups(&[0.0, 0.5, 1.0]),
            ReactionWeightedLimits::default(),
            || Ok(()),
        )
        .unwrap();
    for component in Component::ALL {
        assert_eq!(collapsed.component_j(component)[0], None);
        near(collapsed.component_j(component)[1].unwrap(), 100.0 * EV_J);
    }
    assert_eq!(collapsed.incident_energy_j()[0], None);
    near(collapsed.incident_energy_j()[1].unwrap(), 0.75 * EV_J);
    assert_eq!(collapsed.reaction_weighted().cross_sections_m2[0], 0.0);
    assert!(collapsed.reaction_weighted().flux.values_m2_s()[0] > 0.0);
}

#[test]
fn g4_refinement_budgets_negative_interiors_and_cancellation_refuse_atomically() {
    let evaluation = Arc::new(
        parse(&fixture(
            &[
                coefficients([2.0; 9]),
                coefficients([3.0; 9]),
                coefficients([4.0; 9]),
            ],
            &[],
        ))
        .unwrap(),
    );
    let controls = FissionEnergyPreparationLimits {
        absolute_tolerance_j: 0.01 * EV_J,
        relative_tolerance: 0.0,
        ..FissionEnergyPreparationLimits::default()
    };
    for limits in [
        FissionEnergyPreparationLimits {
            max_depth: 0,
            ..controls
        },
        FissionEnergyPreparationLimits {
            max_evaluations: 1,
            ..controls
        },
        FissionEnergyPreparationLimits {
            max_curve_points: 3,
            ..controls
        },
    ] {
        assert!(
            PreparedFissionEnergy::new(
                Arc::clone(&evaluation),
                domain(),
                ReferenceEnergyTreatment::Reject,
                limits,
                || Ok(())
            )
            .is_err()
        );
    }
    // The endpoints are positive but q=(E-1/2)^2-0.1 is negative inside.
    // A loose interpolation allowance must not admit a negative physical law.
    let negative_inside = Arc::new(
        parse(&fixture(
            &[
                coefficients([0.15; 9]),
                coefficients([-1.0; 9]),
                coefficients([1.0; 9]),
            ],
            &[],
        ))
        .unwrap(),
    );
    assert!(
        PreparedFissionEnergy::new(
            negative_inside,
            domain(),
            ReferenceEnergyTreatment::Reject,
            FissionEnergyPreparationLimits {
                absolute_tolerance_j: EV_J,
                relative_tolerance: 0.0,
                ..controls
            },
            || Ok(())
        )
        .is_err()
    );
    let mut polls = 0;
    PreparedFissionEnergy::new(
        Arc::clone(&evaluation),
        domain(),
        ReferenceEnergyTreatment::Reject,
        controls,
        || {
            polls += 1;
            Ok(())
        },
    )
    .unwrap();
    for cutoff in [1, polls / 2, polls] {
        let mut calls = 0;
        let result = PreparedFissionEnergy::new(
            Arc::clone(&evaluation),
            domain(),
            ReferenceEnergyTreatment::Reject,
            controls,
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
    }
    assert_eq!(evaluation.polynomial_order(), 2);
    near(
        evaluation
            .value_at_energy_j(Component::Fragments, EV_J)
            .unwrap(),
        9.0 * EV_J,
    );
}
