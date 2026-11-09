//! Original parcel enthalpy, synchronous donor chronology and independent amounts.

use super::*;
use crate::miscibility::multicomponent::{
    MulticomponentMiscibilityError,
    advection::{MulticomponentMiscibilityAdvection, MulticomponentMiscibilityAdvectionFace},
    nonisothermal::{
        MulticomponentNonisothermalError, MulticomponentNonisothermalField,
        MulticomponentNonisothermalSolve,
        tests::{Case, enthalpy, entropy, near},
    },
    tests::matrix,
};

fn opposing(
    field: &MulticomponentNonisothermalField,
    rate: f64,
) -> MulticomponentMiscibilityAdvection {
    field
        .network()
        .prepare_advection(
            vec![MulticomponentMiscibilityAdvectionFace {
                forward_m3_s: rate,
                reverse_m3_s: rate,
            }],
            0.0,
            || false,
        )
        .unwrap()
}

fn advance(
    field: &MulticomponentNonisothermalField,
    flow: &MulticomponentMiscibilityAdvection,
    duration: f64,
) -> MulticomponentNonisothermalAdvectionCandidate {
    field
        .advect_candidate(
            flow,
            duration,
            MulticomponentNonisothermalAdvectionControls {
                maximum_courant: 0.25,
                ..Default::default()
            },
            MulticomponentNonisothermalSolve::default(),
            || false,
        )
        .unwrap()
}

#[test]
fn three_cell_cycle_transports_original_enthalpy_from_synchronous_unequal_donors() {
    let case = Case {
        a: matrix(3, &[(0, 1, 0.2), (0, 2, -0.1), (1, 2, 0.3)]),
        b: matrix(3, &[(0, 1, 120.0), (0, 2, -60.0), (1, 2, 80.0)]),
        volumes: vec![0.001, 0.002, 0.0005],
        fractions: vec![[0.8, 0.15, 0.05], [0.2, 0.5, 0.3], [0.1, 0.2, 0.7]],
        temperatures: vec![310.0, 430.0, 610.0],
        links: vec![[2, 1], [2, 0], [1, 0]],
        ..Case::default()
    };
    let initial = case.field();
    // Canonical face order is [0,1], [0,2], [1,2]; the cycle is 0->1->2->0.
    let rate = 0.00025;
    let flow = initial
        .network()
        .prepare_advection(
            vec![
                MulticomponentMiscibilityAdvectionFace {
                    forward_m3_s: rate,
                    reverse_m3_s: 0.0,
                },
                MulticomponentMiscibilityAdvectionFace {
                    forward_m3_s: 0.0,
                    reverse_m3_s: rate,
                },
                MulticomponentMiscibilityAdvectionFace {
                    forward_m3_s: rate,
                    reverse_m3_s: 0.0,
                },
            ],
            0.0,
            || false,
        )
        .unwrap();
    let candidate = advance(&initial, &flow, 0.4);
    assert_eq!(candidate.step().substeps, 1);
    near(candidate.step().maximum_courant, 0.2, 2e-15);
    let initial_h: Vec<_> = (0..3)
        .map(|cell| {
            enthalpy(
                initial.cell_component_moles(cell).unwrap(),
                case.temperatures[cell],
                &case,
            )
        })
        .collect();
    let mut total_h_change = 0.0;
    let mut total_s_change = 0.0;
    for (cell, donor) in [2, 0, 1].into_iter().enumerate() {
        let old = initial.cell_component_moles(cell).unwrap();
        let after = candidate.field().cell_component_moles(cell).unwrap();
        for species in 0..3 {
            near(
                after[species],
                old[species]
                    + 0.1 * (case.fractions[donor][species] - case.fractions[cell][species]),
                3e-15,
            );
        }
        let own_fraction = 0.0001 / case.volumes[cell];
        let donor_fraction = 0.0001 / case.volumes[donor];
        let expected_h = (1.0 - own_fraction) * initial_h[cell] + donor_fraction * initial_h[donor];
        let capacity: f64 = (0..3).map(|s| after[s] * case.cp[s]).sum();
        let expected_t = case.datum + (expected_h - enthalpy(after, case.datum, &case)) / capacity;
        let actual_t = candidate.field().temperatures()[cell].value();
        near(actual_t, expected_t, 1e-9);
        let final_h = enthalpy(after, actual_t, &case);
        near(final_h, expected_h, 1e-8);
        total_h_change += final_h - initial_h[cell];
        total_s_change +=
            entropy(after, actual_t, &case) - entropy(old, case.temperatures[cell], &case);
        near(after.iter().sum(), 1000.0 * case.volumes[cell], 3e-15);
    }
    for (face, (donor, sign)) in [(0, 1.0), (2, -1.0), (1, 1.0)].into_iter().enumerate() {
        for s in 0..3 {
            near(
                candidate.step().transferred_mol[face][s],
                sign * 0.1 * case.fractions[donor][s],
                3e-15,
            );
        }
        near(
            candidate.step().transferred_energy_j[face],
            sign * (0.0001 / case.volumes[donor]) * initial_h[donor],
            2e-10,
        );
    }
    near(total_h_change, 0.0, 1e-8);
    near(
        candidate.step().entropy_change_j_per_k,
        total_s_change,
        1e-9,
    );
    assert!(total_s_change > 0.0);
    assert_eq!(candidate.field().elapsed_s(), 0.4);
    assert_eq!(initial.elapsed_s(), 0.0);
}

#[test]
fn counterflow_relaxes_composition_and_enthalpy_toward_the_physical_thermal_equilibrium() {
    let case = Case {
        b: matrix(3, &[(0, 1, 120.0), (0, 2, -60.0), (1, 2, 80.0)]),
        temperatures: vec![350.0, 550.0],
        ..Case::default()
    };
    let initial = case.field();
    let flow = opposing(&initial, 0.0005);
    let initial_h: Vec<_> = (0..2)
        .map(|cell| {
            enthalpy(
                initial.cell_component_moles(cell).unwrap(),
                case.temperatures[cell],
                &case,
            )
        })
        .collect();
    for (duration, substeps) in [(2.0, 4), (20.0, 40)] {
        let mut controls = MulticomponentNonisothermalSolve::default();
        if substeps == 40 {
            // The long interval charges every sampled parcel's arithmetic
            // resolution; endpoint physics below still has tighter assertions.
            controls.energy_tolerance_j = 1e-7;
        }
        let candidate = initial
            .advect_candidate(
                &flow,
                duration,
                MulticomponentNonisothermalAdvectionControls {
                    maximum_courant: 0.25,
                    ..Default::default()
                },
                controls,
                || false,
            )
            .unwrap();
        assert_eq!(candidate.step().substeps, substeps);
        near(candidate.step().maximum_courant, 0.25, 2e-15);
        // Each synchronous half-second exchange contracts every conserved
        // concentration contrast by 1-Q*dt*(1/Vl+1/Vr)=0.625, including H/V.
        let decay = 0.625_f64.powi(substeps as i32);
        let mean_h = (initial_h[0] + initial_h[1]) / 3.0;
        let contrast_h = (initial_h[1] / 2.0 - initial_h[0]) * decay;
        let expected_h = [
            mean_h - 2.0 * contrast_h / 3.0,
            2.0 * (mean_h + contrast_h / 3.0),
        ];
        for cell in 0..2 {
            let amounts = candidate.field().cell_component_moles(cell).unwrap();
            for s in 0..3 {
                let mean_x = (case.fractions[0][s] + 2.0 * case.fractions[1][s]) / 3.0;
                let contrast_x = (case.fractions[1][s] - case.fractions[0][s]) * decay;
                let expected_x = if cell == 0 {
                    mean_x - 2.0 * contrast_x / 3.0
                } else {
                    mean_x + contrast_x / 3.0
                };
                near(
                    candidate.field().cell_mole_fractions(cell).unwrap()[s],
                    expected_x,
                    3e-14,
                );
                let signed_receipt = if cell == 0 { -1.0 } else { 1.0 };
                near(
                    amounts[s] - initial.cell_component_moles(cell).unwrap()[s],
                    signed_receipt * candidate.step().transferred_mol[0][s],
                    3e-14,
                );
            }
            let capacity: f64 = (0..3).map(|s| amounts[s] * case.cp[s]).sum();
            let expected_t =
                case.datum + (expected_h[cell] - enthalpy(amounts, case.datum, &case)) / capacity;
            near(
                candidate.field().temperatures()[cell].value(),
                expected_t,
                2e-9,
            );
        }
        near(
            candidate.step().transferred_energy_j[0],
            expected_h[1] - initial_h[1],
            1e-8,
        );
        near(candidate.step().energy_residual_j, 0.0, 1e-8);
        assert!(candidate.step().entropy_change_j_per_k > 0.0);
        if duration == 20.0 {
            // The equilibrium is dictated by conserved H and all inventories;
            // unequal species Cp rules out averaging the parcel temperatures.
            let total: Vec<_> = (0..3)
                .map(|s| initial.component_moles()[s] + initial.component_moles()[3 + s])
                .collect();
            let capacity: f64 = (0..3).map(|s| total[s] * case.cp[s]).sum();
            let equilibrium = case.datum
                + (initial_h.iter().sum::<f64>() - enthalpy(&total, case.datum, &case)) / capacity;
            assert!((equilibrium - (350.0 + 2.0 * 550.0) / 3.0).abs() > 1.0);
            for temperature in candidate.field().temperatures() {
                near(temperature.value(), equilibrium, 3e-6);
            }
        }
    }
}

#[test]
fn stirring_keeps_tiny_positive_traces_and_accumulates_sub_ulp_species_receipts() {
    let initial = Case {
        fractions: vec![[1e-280, 0.3, 0.7], [2e-280, 0.3, 0.7]],
        temperatures: vec![300.0; 2],
        ..Case::default()
    }
    .field();
    let flow = opposing(&initial, 0.0005);
    let candidate = advance(&initial, &flow, 0.2);
    near(candidate.step().transferred_mol[0][0] / 1e-281, -1.0, 2e-14);
    near(candidate.field().component_moles()[0] / 1e-280, 1.1, 2e-14);
    near(candidate.field().component_moles()[3] / 1e-280, 3.9, 2e-14);
    assert!(candidate.field().component_moles().iter().all(|&n| n > 0.0));

    let initial = Case {
        temperatures: vec![300.0; 2],
        ..Case::default()
    }
    .field();
    let flow = opposing(&initial, 1e-20);
    let first = advance(&initial, &flow, 1.0);
    let second = advance(first.field(), &flow, 1.0);
    assert_eq!(second.field().component_moles(), initial.component_moles());
    assert!(
        second
            .field()
            .state
            .amount_remainders
            .iter()
            .any(|&r| r != 0.0)
    );
    for s in 0..3 {
        let receipt = first.step().transferred_mol[0][s] + second.step().transferred_mol[0][s];
        near(second.field().state.amount_remainders[s], -receipt, 1e-32);
        near(
            second.field().state.amount_remainders[3 + s],
            receipt,
            1e-32,
        );
    }
    assert_eq!(second.field().elapsed_s(), 2.0);
}

#[test]
fn work_limits_source_domain_and_late_cancellation_leave_the_original_thermal_owner_untouched() {
    let initial = Case::default().field();
    let saved = initial.clone();
    let flow = opposing(&initial, 0.0005);
    let controls = MulticomponentNonisothermalAdvectionControls {
        maximum_courant: 0.1,
        maximum_substeps: 1,
        ..Default::default()
    };
    assert!(matches!(
        initial.advect_candidate(
            &flow,
            2.0,
            controls,
            MulticomponentNonisothermalSolve::default(),
            || false
        ),
        Err(MulticomponentNonisothermalError::Spatial(
            MulticomponentMiscibilityError::AdvectionWorkLimit { .. }
        ))
    ));
    let mut polls = 0usize;
    let accepted = initial
        .advect_candidate(&flow, 0.2, Default::default(), Default::default(), || {
            polls += 1;
            false
        })
        .unwrap();
    assert!(polls > 3);
    let mut seen = 0usize;
    assert!(matches!(
        initial.advect_candidate(&flow, 0.2, Default::default(), Default::default(), || {
            seen += 1;
            seen == polls
        }),
        Err(MulticomponentNonisothermalError::Cancelled)
    ));
    assert_eq!(initial, saved);
    assert_eq!(
        accepted,
        initial
            .advect_candidate(&flow, 0.2, Default::default(), Default::default(), || false)
            .unwrap()
    );

    let cold = Case {
        b: matrix(3, &[(0, 1, 120.0)]),
        cp: [80.0; 3],
        volumes: vec![0.001; 2],
        fractions: vec![[0.8, 0.1, 0.1], [0.1, 0.8, 0.1]],
        temperatures: vec![250.0; 2],
        ..Case::default()
    }
    .field();
    let saved = cold.clone();
    let flow = opposing(&cold, 0.0005);
    // Endothermic homogenization would cool these original parcels below
    // their declared lower temperature bracket. No clamped endpoint is legal.
    assert!(
        cold.advect_candidate(&flow, 0.5, Default::default(), Default::default(), || false)
            .is_err()
    );
    assert_eq!(cold, saved);
}

#[test]
fn unresolvable_absolute_parcel_enthalpy_cannot_create_a_reference_dependent_temperature() {
    let mut case = Case {
        volumes: vec![0.001; 2],
        fractions: vec![[0.25, 0.25, 0.5]; 2],
        ..Case::default()
    };
    let initial = case.field();
    let flow = opposing(&initial, 0.0005);
    let ordinary = advance(&initial, &flow, 0.5);
    assert_eq!(
        ordinary.field().component_moles(),
        initial.component_moles()
    );
    near(ordinary.step().transferred_energy_j[0], -4125.0, 1e-10);
    near(ordinary.field().temperatures()[0].value(), 375.0, 1e-9);
    near(ordinary.field().temperatures()[1].value(), 525.0, 1e-9);

    case.datum = 1e15;
    let initial = case.field();
    let saved = initial.clone();
    let flow = opposing(&initial, 0.0005);
    // All inventories, including low words, stay fixed. Rounding each huge
    // source parcel H before subtracting the counterstreams can nevertheless
    // shift the represented heat from 4125 J to 4126 J. Global H/S alone cannot
    // detect the resulting spurious 0.018 K redistribution.
    assert!(matches!(
        initial.advect_candidate(&flow, 0.5, Default::default(), Default::default(), || false,),
        Err(MulticomponentNonisothermalError::UnresolvedParcelEnthalpy { .. })
    ));
    assert_eq!(initial, saved);
}
