//! Independent parcel-energy balances, donor chronology and atomic refusal.

use super::*;
use crate::miscibility::advection::BinaryMiscibilityAdvectionFace;
use crate::miscibility::{BinaryMiscibilityCell, BinaryMiscibilityFace};
use fs_qty::AmountConcentration;
use fs_thermochem::solution::caloric::PureLiquidHeatCapacityV1;
use fs_thermochem::solution::{BinaryRegularSolutionV1, gradient::BinarySquareGradientV1};
use fs_thermochem::{MolarThermalQuantityV1, SpeciesId, UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R};

pub(in crate::miscibility::nonisothermal) fn fixture(
    cp: [f64; 2],
    interaction_kelvin: f64,
    diffusivity: f64,
    conductivity: f64,
    cells: &[[f64; 3]], // [volume m³, x_A, temperature K]
    faces: &[[usize; 2]],
) -> BinaryNonisothermalField {
    let bracket = [Temperature::new(250.0), Temperature::new(700.0)];
    let species = [SpeciesId::new("A").unwrap(), SpeciesId::new("B").unwrap()];
    let pair = BinaryRegularSolutionV1::new(
        species.clone(),
        0.0,
        interaction_kelvin,
        bracket[0],
        bracket[1],
    )
    .unwrap();
    let capacities = core::array::from_fn(|s| {
        PureLiquidHeatCapacityV1::new(
            species[s].clone(),
            Temperature::new(300.0),
            bracket
                .map(|t| (t, MolarThermalQuantityV1::new(cp[s])))
                .to_vec(),
        )
        .unwrap()
    });
    let caloric = BinaryLiquidCaloricModelV1::new(pair.clone(), capacities).unwrap();
    let model = BinarySquareGradientV1::new(
        pair,
        Temperature::new(300.0),
        AmountConcentration::new(1000.0),
        0.0,
    )
    .unwrap();
    let network = BinaryMiscibilityNetwork::new(
        model,
        cells
            .iter()
            .map(|cell| BinaryMiscibilityCell { volume_m3: cell[0] })
            .collect(),
        faces
            .iter()
            .map(|&cells| BinaryMiscibilityFace {
                cells,
                area_m2: 1.0,
                distance_m: 1.0,
                diffusivity_m2_s: diffusivity,
            })
            .collect(),
        1024 * 1024,
        || false,
    )
    .unwrap();
    let state = network
        .state_from_mole_fractions(cells.iter().map(|c| [c[1], 1.0 - c[1]]).collect(), || false)
        .unwrap();
    network
        .bind_nonisothermal_field(
            &state,
            caloric,
            cells.iter().map(|c| Temperature::new(c[2])).collect(),
            vec![
                BinaryNonisothermalFace {
                    thermal_conductivity_w_m_k: conductivity,
                    heat_of_transport_j_per_mol: 0.0,
                };
                faces.len()
            ],
            bracket,
            1024 * 1024,
            || false,
        )
        .unwrap()
}

pub(in crate::miscibility::nonisothermal) fn opposing_flow(
    field: &BinaryNonisothermalField,
    rate: f64,
) -> BinaryMiscibilityAdvection {
    field
        .network()
        .prepare_advection(
            vec![BinaryMiscibilityAdvectionFace {
                forward_m3_s: rate,
                reverse_m3_s: rate,
            }],
            0.0,
            || false,
        )
        .unwrap()
}

pub(in crate::miscibility::nonisothermal) fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual {actual:.17e}, expected {expected:.17e}, tolerance {tolerance:e}"
    );
}

// Independent constant-pure-Cp, chi=B/T parcel formula with the 300 K datum.
pub(in crate::miscibility::nonisothermal) fn enthalpy(
    amounts: [f64; 2],
    temperature: f64,
    cp: [f64; 2],
    interaction_kelvin: f64,
) -> f64 {
    (amounts[0] * cp[0] + amounts[1] * cp[1]) * (temperature - 300.0)
        + R * interaction_kelvin * amounts[0] * amounts[1] / (amounts[0] + amounts[1])
}

fn advect(
    initial: &BinaryNonisothermalField,
    flow: &BinaryMiscibilityAdvection,
    duration: f64,
) -> BinaryNonisothermalAdvectionCandidate {
    initial
        .advect_candidate(
            flow,
            duration,
            BinaryMiscibilityAdvectionControls::default(),
            BinaryNonisothermalSolve::default(),
            || false,
        )
        .unwrap()
}

#[test]
fn g1_unequal_heat_capacities_advect_enthalpy_instead_of_temperature() {
    let cp = [20.0, 60.0];
    let initial = fixture(
        cp,
        0.0,
        0.0,
        0.0,
        &[[0.001, 0.8, 350.0], [0.001, 0.2, 550.0]],
        &[[0, 1]],
    );
    let flow = opposing_flow(&initial, 0.0005);
    let candidate = advect(&initial, &flow, 0.5);
    let final_field = candidate.field();
    // A quarter of each donor parcel is exchanged. New capacities are 34 and
    // 46 J/K; new enthalpies are 4300 and 10100 J. Averaging T gives 400/500 K.
    for (cell, (x, h, capacity)) in [(0.65, 4300.0, 34.0), (0.35, 10100.0, 46.0)]
        .into_iter()
        .enumerate()
    {
        near(final_field.mole_fractions()[cell][0], x, 2e-14);
        near(
            final_field.temperatures()[cell].value(),
            300.0 + h / capacity,
            1e-8,
        );
        near(
            enthalpy(
                final_field.component_moles()[cell],
                final_field.temperatures()[cell].value(),
                cp,
                0.0,
            ),
            h,
            1e-8,
        );
    }
    near(
        candidate.step().transferred_component_moles[0][0],
        0.15,
        2e-14,
    );
    near(
        candidate.step().transferred_component_moles[0][1],
        -0.15,
        2e-14,
    );
    near(candidate.step().transferred_energy_j[0], -2900.0, 1e-10);
    assert!(candidate.step().entropy_change_j_per_k > 0.0);
    near(candidate.step().energy_residual_j, 0.0, 1e-8);
    assert_eq!(candidate.step().substeps, 1);
    assert_eq!(final_field.elapsed_s(), 0.5);
}

#[test]
fn g1_nonideal_mixing_heat_is_recovered_once_from_original_enthalpy() {
    for interaction_kelvin in [-600.0, 300.0] {
        let initial = fixture(
            [40.0; 2],
            interaction_kelvin,
            0.0,
            0.0,
            &[[0.001, 0.9, 350.0], [0.001, 0.1, 350.0]],
            &[[0, 1]],
        );
        let flow = opposing_flow(&initial, 0.0005);
        let candidate = advect(&initial, &flow, 0.5);
        // Both original parcels have H_ex=R*B*.09. Mixing gives x=.7/.3,
        // hence H_ex=R*B*.21; the exact sensible compensation is -R*B*.12.
        let expected_t = 350.0 + R * interaction_kelvin * (0.09 - 0.21) / 40.0;
        for cell in 0..2 {
            near(
                candidate.field().temperatures()[cell].value(),
                expected_t,
                1e-8,
            );
            near(
                candidate.field().mole_fractions()[cell][0],
                [0.7, 0.3][cell],
                2e-14,
            );
        }
        near(candidate.step().transferred_energy_j[0], 0.0, 1e-10);
        near(candidate.step().energy_residual_j, 0.0, 1e-8);
        assert!(candidate.step().entropy_change_j_per_k > 0.0);
    }
}

#[test]
fn g1_directed_circulation_samples_all_original_donors_synchronously() {
    let cells = [
        [0.001, 0.8, 310.0],
        [0.001, 0.35, 430.0],
        [0.001, 0.15, 610.0],
    ];
    let cp = [20.0, 60.0];
    // Deliberately reversed/reordered input faces still produce canonical
    // rates [0->1, 2->0, 1->2]. No just-received parcel may be sent onward.
    let initial = fixture(cp, 0.0, 0.0, 0.0, &cells, &[[2, 1], [2, 0], [1, 0]]);
    let q = 0.0005;
    let flow = initial
        .network()
        .prepare_advection(
            vec![
                BinaryMiscibilityAdvectionFace {
                    forward_m3_s: q,
                    reverse_m3_s: 0.0,
                },
                BinaryMiscibilityAdvectionFace {
                    forward_m3_s: 0.0,
                    reverse_m3_s: q,
                },
                BinaryMiscibilityAdvectionFace {
                    forward_m3_s: q,
                    reverse_m3_s: 0.0,
                },
            ],
            0.0,
            || false,
        )
        .unwrap();
    let candidate = advect(&initial, &flow, 0.5);
    for (cell, donor) in [2, 0, 1].into_iter().enumerate() {
        let x = 0.75 * cells[cell][1] + 0.25 * cells[donor][1];
        let h = 0.75 * enthalpy(initial.component_moles()[cell], cells[cell][2], cp, 0.0)
            + 0.25 * enthalpy(initial.component_moles()[donor], cells[donor][2], cp, 0.0);
        near(candidate.field().mole_fractions()[cell][0], x, 2e-14);
        near(
            candidate.field().temperatures()[cell].value(),
            300.0 + h / (cp[0] * x + cp[1] * (1.0 - x)),
            1e-8,
        );
    }
    for (face, (donor, sign)) in [(0, 1.0), (2, -1.0), (1, 1.0)].into_iter().enumerate() {
        for s in 0..2 {
            near(
                candidate.step().transferred_component_moles[face][s],
                sign * 0.25 * initial.component_moles()[donor][s],
                2e-14,
            );
        }
    }
    near(candidate.step().energy_residual_j, 0.0, 1e-8);
}

#[test]
fn g0_unequal_volumes_subcycle_without_losing_species_or_enthalpy() {
    let cp = [25.0, 45.0];
    let b = -100.0;
    let initial = fixture(
        cp,
        b,
        0.0,
        0.0,
        &[[0.001, 0.25, 320.0], [0.002, 0.75, 600.0]],
        &[[0, 1]],
    );
    let flow = opposing_flow(&initial, 0.0005);
    let candidate = initial
        .advect_candidate(
            &flow,
            2.0,
            BinaryMiscibilityAdvectionControls {
                maximum_courant: 0.25,
                ..Default::default()
            },
            BinaryNonisothermalSolve::default(),
            || false,
        )
        .unwrap();
    assert_eq!(candidate.step().substeps, 4);
    near(candidate.step().maximum_courant, 0.25, 2e-15);
    let mean = 7.0 / 12.0;
    let contrast = 0.5 * 0.625_f64.powi(4);
    for (cell, x) in [mean - 2.0 * contrast / 3.0, mean + contrast / 3.0]
        .into_iter()
        .enumerate()
    {
        near(candidate.field().mole_fractions()[cell][0], x, 2e-14);
        let total: f64 = candidate.field().component_moles()[cell].iter().sum();
        near(total, [1.0, 2.0][cell], 2e-14);
        let sign = if cell == 0 { -1.0 } else { 1.0 };
        for s in 0..2 {
            near(
                candidate.field().component_moles()[cell][s] - initial.component_moles()[cell][s],
                sign * candidate.step().transferred_component_moles[0][s],
                2e-13,
            );
        }
        let delta_h = enthalpy(
            candidate.field().component_moles()[cell],
            candidate.field().temperatures()[cell].value(),
            cp,
            b,
        ) - enthalpy(
            initial.component_moles()[cell],
            initial.temperatures()[cell].value(),
            cp,
            b,
        );
        near(
            delta_h,
            sign * candidate.step().transferred_energy_j[0],
            1e-8,
        );
    }
    assert!(candidate.step().entropy_change_j_per_k > 0.0);
    near(candidate.step().energy_residual_j, 0.0, 1e-8);
}

#[test]
fn g0_work_limits_and_foreign_flows_refuse_without_mutating_the_field() {
    let initial = fixture(
        [30.0; 2],
        0.0,
        0.0,
        0.0,
        &[[0.001, 0.2, 350.0], [0.001, 0.7, 500.0]],
        &[[0, 1]],
    );
    let original = initial.clone();
    let flow = opposing_flow(&initial, 0.001);
    let limited = initial.advect_candidate(
        &flow,
        2.0,
        BinaryMiscibilityAdvectionControls {
            maximum_substeps: 1,
            ..Default::default()
        },
        BinaryNonisothermalSolve::default(),
        || false,
    );
    assert!(matches!(
        limited,
        Err(BinaryNonisothermalError::Spatial(
            BinaryMiscibilityError::AdvectionWorkLimit { .. }
        ))
    ));
    let other = fixture(
        [30.0; 2],
        0.0,
        0.0,
        0.0,
        &[[0.002, 0.2, 350.0], [0.001, 0.7, 500.0]],
        &[[0, 1]],
    );
    let foreign = opposing_flow(&other, 0.001);
    assert!(matches!(
        initial.advect_candidate(
            &foreign,
            0.1,
            BinaryMiscibilityAdvectionControls::default(),
            BinaryNonisothermalSolve::default(),
            || false
        ),
        Err(BinaryNonisothermalError::Spatial(
            BinaryMiscibilityError::Invalid(_)
        ))
    ));
    assert_eq!(initial, original);
}

#[test]
fn g4_late_cancellation_is_atomic_and_replay_is_identical() {
    let initial = fixture(
        [20.0, 60.0],
        -100.0,
        0.0,
        0.0,
        &[[0.001, 0.3, 350.0], [0.001, 0.7, 500.0]],
        &[[0, 1]],
    );
    let original = initial.clone();
    let flow = opposing_flow(&initial, 0.0004);
    let mut polls = 0;
    let expected = initial
        .advect_candidate(
            &flow,
            0.7,
            BinaryMiscibilityAdvectionControls::default(),
            BinaryNonisothermalSolve::default(),
            || {
                polls += 1;
                false
            },
        )
        .unwrap();
    let mut cancelled_polls = 0;
    assert!(matches!(
        initial.advect_candidate(
            &flow,
            0.7,
            BinaryMiscibilityAdvectionControls::default(),
            BinaryNonisothermalSolve::default(),
            || {
                cancelled_polls += 1;
                cancelled_polls == polls
            }
        ),
        Err(BinaryNonisothermalError::Cancelled)
    ));
    assert_eq!(initial, original);
    assert_eq!(advect(&initial, &flow, 0.7), expected);
}
