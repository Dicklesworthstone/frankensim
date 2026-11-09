//! Independent friction, caloric, and whole-field entropy oracles.

use super::*;
use crate::miscibility::multicomponent::{
    MulticomponentMiscibilityCell, MulticomponentMiscibilityFace, MulticomponentMiscibilityField,
    MulticomponentMiscibilityNetwork, MulticomponentMiscibilitySolve,
    tests::{matrix, uniform_matrix},
};
use fs_qty::{AmountConcentration, MolarMass, Temperature};
use fs_thermochem::solution::caloric::PureLiquidHeatCapacityV1;
use fs_thermochem::solution::multicomponent::{
    MulticomponentRegularSolutionV1, caloric::MulticomponentLiquidCaloricModelV1,
};
use fs_thermochem::transport::maxwell_stefan::MaxwellStefanV1;
use fs_thermochem::{MolarThermalQuantityV1, SpeciesId, UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R};

const WORKSPACE: usize = 16 * 1024 * 1024;
const CONCENTRATION: f64 = 1000.0;

pub(super) fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual {actual:.17e}, expected {expected:.17e}, tolerance {tolerance:e}"
    );
}

fn species() -> Vec<SpeciesId> {
    (0..3)
        .map(|i| SpeciesId::new(format!("S{i:03}")).unwrap())
        .collect()
}

/// Small original-law fixture shared with the adaptive numerical tests.
pub(super) struct Case {
    pub(super) a: Vec<f64>,
    pub(super) b: Vec<f64>,
    pub(super) cp: [f64; 3],
    pub(super) datum: f64,
    pub(super) volumes: Vec<f64>,
    pub(super) fractions: Vec<[f64; 3]>,
    pub(super) temperatures: Vec<f64>,
    pub(super) links: Vec<[usize; 2]>,
    pub(super) diffusivity: Vec<f64>,
    pub(super) conductivity: f64,
    pub(super) heat: [f64; 3],
}

impl Default for Case {
    fn default() -> Self {
        Self {
            a: vec![0.0; 9],
            b: vec![0.0; 9],
            cp: [20.0, 40.0, 80.0],
            datum: 300.0,
            volumes: vec![0.001, 0.002],
            fractions: vec![[0.2, 0.3, 0.5], [0.4, 0.4, 0.2]],
            temperatures: vec![300.0, 600.0],
            links: vec![[0, 1]],
            diffusivity: uniform_matrix(3, 1e-4),
            conductivity: 0.0,
            heat: [0.0; 3],
        }
    }
}

impl Case {
    fn parts(
        &self,
    ) -> (
        MulticomponentMiscibilityNetwork,
        MulticomponentMiscibilityField,
        MulticomponentLiquidCaloricModelV1,
    ) {
        let law = MulticomponentRegularSolutionV1::new(
            species(),
            self.a.clone(),
            self.b.clone(),
            Temperature::new(100.0),
            Temperature::new(1000.0),
        )
        .unwrap();
        let caloric = MulticomponentLiquidCaloricModelV1::new(
            law.clone(),
            species()
                .into_iter()
                .zip(self.cp)
                .map(|(id, cp)| {
                    PureLiquidHeatCapacityV1::new(
                        id,
                        Temperature::new(self.datum),
                        vec![
                            (Temperature::new(100.0), MolarThermalQuantityV1::new(cp)),
                            (Temperature::new(1e15), MolarThermalQuantityV1::new(cp)),
                        ],
                    )
                    .unwrap()
                })
                .collect(),
        )
        .unwrap();
        let transport = MaxwellStefanV1::new(
            species()
                .into_iter()
                .zip([0.002, 0.028, 0.044])
                .map(|(id, mass)| (id, MolarMass::new(mass)))
                .collect(),
        )
        .unwrap();
        let network = MulticomponentMiscibilityNetwork::new(
            law,
            transport,
            Temperature::new(300.0),
            AmountConcentration::new(CONCENTRATION),
            self.volumes
                .iter()
                .map(|&volume_m3| MulticomponentMiscibilityCell { volume_m3 })
                .collect(),
            self.links
                .iter()
                .map(|&cells| MulticomponentMiscibilityFace {
                    cells,
                    area_m2: 1.0,
                    distance_m: 1.0,
                    diffusivity_matrix: 0,
                })
                .collect(),
            vec![self.diffusivity.clone()],
            WORKSPACE,
            || false,
        )
        .unwrap();
        let initial = network
            .bind(
                self.fractions
                    .iter()
                    .zip(network.capacities_mol())
                    .map(|(x, &n)| x.iter().map(|&x| n * x).collect())
                    .collect(),
                || false,
            )
            .unwrap();
        (network, initial, caloric)
    }

    fn faces(&self) -> Vec<MulticomponentNonisothermalFace> {
        self.links
            .iter()
            .map(|_| MulticomponentNonisothermalFace {
                thermal_conductivity_w_m_k: self.conductivity,
                heat_of_transport_j_per_mol: self.heat.to_vec(),
            })
            .collect()
    }

    pub(super) fn field(&self) -> MulticomponentNonisothermalField {
        let (network, initial, caloric) = self.parts();
        network
            .bind_nonisothermal_field(
                &initial,
                caloric,
                self.temperatures
                    .iter()
                    .copied()
                    .map(Temperature::new)
                    .collect(),
                self.faces(),
                [Temperature::new(250.0), Temperature::new(700.0)],
                WORKSPACE,
                || false,
            )
            .unwrap()
    }
}

fn enthalpy(amounts: &[f64], temperature: f64, case: &Case) -> f64 {
    let n: f64 = amounts.iter().sum();
    let sensible: f64 = (0..3)
        .map(|i| amounts[i] * case.cp[i] * (temperature - case.datum))
        .sum();
    sensible
        + R * (0..3)
            .flat_map(|i| (i + 1..3).map(move |j| (i, j)))
            .map(|(i, j)| case.b[3 * i + j] * amounts[i] * amounts[j] / n)
            .sum::<f64>()
}

fn entropy(amounts: &[f64], temperature: f64, case: &Case) -> f64 {
    let n: f64 = amounts.iter().sum();
    let pure: f64 = (0..3)
        .map(|i| amounts[i] * case.cp[i] * (temperature / case.datum).ln())
        .sum();
    let ideal: f64 = amounts.iter().map(|&a| a * (a / n).ln()).sum();
    let excess: f64 = (0..3)
        .flat_map(|i| (i + 1..3).map(move |j| (i, j)))
        .map(|(i, j)| case.a[3 * i + j] * amounts[i] * amounts[j] / n)
        .sum();
    pure - R * (ideal + excess)
}

// Independent 2x2 elimination of the original ternary friction equations.
// The complete force includes the authored Soret term, not a guessed Fick D.
fn face_oracle(case: &Case, left: usize, right: usize) -> ([f64; 3], f64, f64) {
    let [xl, xr] = [case.fractions[left], case.fractions[right]];
    let [tl, tr] = [case.temperatures[left], case.temperatures[right]];
    let x: [f64; 3] = core::array::from_fn(|s| 0.5 * (xl[s] + xr[s]));
    let beta_jump = 1.0 / tr - 1.0 / tl;
    let reduced: [f64; 3] = core::array::from_fn(|s| {
        (xr[s] / xl[s]).ln()
            + (0..3)
                .map(|j| {
                    (case.a[3 * s + j] + 0.5 * (1.0 / tl + 1.0 / tr) * case.b[3 * s + j])
                        * (xr[j] - xl[j])
                })
                .sum::<f64>()
            - case.heat[s] * beta_jump / R
    });
    let common: f64 = (0..3).map(|s| x[s] * reduced[s]).sum();
    let driving: [f64; 3] = core::array::from_fn(|s| x[s] * (reduced[s] - common));
    let [ab, ac, bc] = [
        case.diffusivity[1],
        case.diffusivity[2],
        case.diffusivity[5],
    ];
    let aa = x[1] / ab + (x[0] + x[2]) / ac;
    let bb = x[0] * (1.0 / ac - 1.0 / ab);
    let cc = x[1] * (1.0 / bc - 1.0 / ab);
    let dd = x[0] / ab + (x[1] + x[2]) / bc;
    let det = aa * dd - bb * cc;
    let ja = CONCENTRATION * (-driving[0] * dd + driving[1] * bb) / det;
    let jb = CONCENTRATION * (-aa * driving[1] + cc * driving[0]) / det;
    let flux = [ja, jb, -ja - jb];
    let thermal_mean = if tl == tr {
        tl
    } else {
        tl * tr * (tr / tl).ln() / (tr - tl)
    };
    let carried: f64 = (0..3)
        .map(|s| {
            let h = case.cp[s] * (thermal_mean - case.datum)
                + R * (0..3).map(|j| case.b[3 * s + j] * x[j]).sum::<f64>();
            (h + case.heat[s]) * flux[s]
        })
        .sum();
    let energy = carried - case.conductivity * (tr - tl);
    let production = -R * (0..3).map(|s| flux[s] * reduced[s]).sum::<f64>()
        + case.conductivity * (tr - tl).powi(2) / (tl * tr);
    (flux, energy, production)
}

#[test]
fn uniform_composition_has_only_fourier_heat_despite_unequal_cp_and_excess_enthalpy() {
    let case = Case {
        b: matrix(3, &[(0, 1, 120.0), (0, 2, -60.0), (1, 2, 80.0)]),
        fractions: vec![[0.2, 0.3, 0.5]; 2],
        conductivity: 2.0,
        ..Case::default()
    };
    let initial = case.field();
    let dt = 0.01;
    let candidate = initial
        .step_candidate(dt, MulticomponentNonisothermalSolve::default(), || false)
        .unwrap();
    let next = candidate.field();
    assert_eq!(next.component_moles(), initial.component_moles());
    assert_eq!(candidate.step().transferred_mol, vec![vec![0.0; 3]]);
    near(candidate.step().transferred_energy_j[0], -6.0, 1e-13);
    near(next.temperatures()[0].value(), 300.0 + 6.0 / 56.0, 1e-9);
    near(next.temperatures()[1].value(), 600.0 - 6.0 / 112.0, 1e-9);
    near(candidate.step().entropy_production_rate_w_per_k, 1.0, 2e-13);
    assert!(candidate.step().entropy_change_j_per_k > 0.0);
}

#[test]
fn unequal_cp_and_nonideal_ternary_transport_match_independent_energy_and_friction() {
    let case = Case {
        a: matrix(3, &[(0, 1, 0.2), (0, 2, -0.1), (1, 2, 0.3)]),
        b: matrix(3, &[(0, 1, 120.0), (0, 2, -60.0), (1, 2, 80.0)]),
        diffusivity: matrix(3, &[(0, 1, 1e-5), (0, 2, 3e-5), (1, 2, 7e-5)]),
        conductivity: 0.1,
        ..Case::default()
    };
    let initial = case.field();
    let dt = 0.02;
    let (flux, energy, production) = face_oracle(&case, 0, 1);
    let candidate = initial
        .step_candidate(dt, MulticomponentNonisothermalSolve::default(), || false)
        .unwrap();
    for s in 0..3 {
        near(candidate.step().transferred_mol[0][s], dt * flux[s], 3e-16);
    }
    near(candidate.step().transferred_energy_j[0], dt * energy, 1e-10);
    near(
        candidate.step().entropy_production_rate_w_per_k,
        production,
        2e-12,
    );
    let mut energy_change = 0.0;
    let mut entropy_change = 0.0;
    for cell in 0..2 {
        let before = initial.cell_component_moles(cell).unwrap();
        let after = candidate.field().cell_component_moles(cell).unwrap();
        let old_h = enthalpy(before, case.temperatures[cell], &case);
        let exchange = if cell == 0 { -dt * energy } else { dt * energy };
        let capacity: f64 = (0..3).map(|s| after[s] * case.cp[s]).sum();
        let mixing_h = enthalpy(after, case.datum, &case);
        let expected_t = case.datum + (old_h + exchange - mixing_h) / capacity;
        let actual_t = candidate.field().temperatures()[cell].value();
        near(actual_t, expected_t, 1e-9);
        energy_change += enthalpy(after, actual_t, &case) - old_h;
        entropy_change +=
            entropy(after, actual_t, &case) - entropy(before, case.temperatures[cell], &case);
    }
    near(energy_change, 0.0, 1e-8);
    near(candidate.step().energy_residual_j, energy_change, 1e-8);
    near(
        candidate.step().entropy_change_j_per_k,
        entropy_change,
        2e-10,
    );
    assert!(entropy_change > 0.0);
}

#[test]
fn equal_temperature_species_rates_reduce_to_the_original_isothermal_owner() {
    let case = Case {
        a: matrix(3, &[(0, 1, 0.2), (0, 2, -0.3), (1, 2, 0.1)]),
        b: matrix(3, &[(0, 1, 90.0), (0, 2, -30.0)]),
        temperatures: vec![300.0; 2],
        heat: [80.0, -40.0, 25.0],
        ..Case::default()
    };
    let (_, isothermal, _) = case.parts();
    let dt = 0.01;
    let original = isothermal
        .step_candidate(dt, MulticomponentMiscibilitySolve::default(), || false)
        .unwrap();
    let thermal = case
        .field()
        .step_candidate(dt, MulticomponentNonisothermalSolve::default(), || false)
        .unwrap();
    for s in 0..3 {
        near(
            thermal.step().transferred_mol[0][s],
            original.step().transferred_mol[0][s],
            2e-16,
        );
    }
    // The thermal owner still recovers local T from energy; an isothermal
    // constitutive rate does not prescribe an isothermal endpoint.
    assert!(
        thermal
            .field()
            .temperatures()
            .iter()
            .any(|t| t.value() != 300.0)
    );
}

#[test]
fn authored_soret_heat_has_reciprocal_energy_and_a_common_gauge_cancels() {
    let mut case = Case {
        fractions: vec![[0.2, 0.3, 0.5]; 2],
        heat: [1500.0, -800.0, 400.0],
        conductivity: 0.2,
        ..Case::default()
    };
    let dt = 0.02;
    let (flux, energy, production) = face_oracle(&case, 0, 1);
    assert!(flux.iter().any(|j| j.abs() > 1e-5));
    let result = case
        .field()
        .step_candidate(dt, MulticomponentNonisothermalSolve::default(), || false)
        .unwrap();
    for s in 0..3 {
        near(result.step().transferred_mol[0][s], dt * flux[s], 2e-16);
    }
    near(result.step().transferred_energy_j[0], dt * energy, 2e-10);
    near(
        result.step().entropy_production_rate_w_per_k,
        production,
        2e-12,
    );
    let reciprocal: f64 = (0..3).map(|s| case.heat[s] * flux[s]).sum();
    assert!(reciprocal.abs() > 0.1);
    for q in &mut case.heat {
        *q += 1234.0;
    }
    let shifted = case
        .field()
        .step_candidate(dt, MulticomponentNonisothermalSolve::default(), || false)
        .unwrap();
    assert_eq!(
        shifted.step().transferred_mol,
        result.step().transferred_mol
    );
    assert_eq!(
        shifted.step().transferred_energy_j,
        result.step().transferred_energy_j
    );
    for s in 0..3 {
        near(
            shifted.step().transferred_mol[0][s],
            result.step().transferred_mol[0][s],
            3e-16,
        );
    }
    near(
        shifted.step().transferred_energy_j[0],
        result.step().transferred_energy_j[0],
        2e-10,
    );
    for cell in 0..2 {
        near(
            shifted.field().temperatures()[cell].value(),
            result.field().temperatures()[cell].value(),
            1e-9,
        );
    }
}

#[test]
fn original_variable_heat_capacities_recover_heat_without_averaging_temperatures() {
    let case = Case {
        fractions: vec![[0.2, 0.3, 0.5]; 2],
        temperatures: vec![300.0, 550.0],
        conductivity: 3.0,
        ..Case::default()
    };
    let (network, initial, _) = case.parts();
    let intercept = [5.0, 10.0, 20.0];
    let slope = [0.05, 0.1, 0.04];
    let caloric = MulticomponentLiquidCaloricModelV1::new(
        network.model().clone(),
        species()
            .into_iter()
            .enumerate()
            .map(|(s, id)| {
                PureLiquidHeatCapacityV1::new(
                    id,
                    Temperature::new(300.0),
                    [100.0, 1000.0]
                        .into_iter()
                        .map(|t| {
                            (
                                Temperature::new(t),
                                MolarThermalQuantityV1::new(intercept[s] + slope[s] * t),
                            )
                        })
                        .collect(),
                )
                .unwrap()
            })
            .collect(),
    )
    .unwrap();
    let initial = network
        .bind_nonisothermal_field(
            &initial,
            caloric,
            case.temperatures
                .iter()
                .copied()
                .map(Temperature::new)
                .collect(),
            case.faces(),
            [Temperature::new(250.0), Temperature::new(700.0)],
            WORKSPACE,
            || false,
        )
        .unwrap();
    let dt = 0.1;
    let result = initial
        .step_candidate(dt, MulticomponentNonisothermalSolve::default(), || false)
        .unwrap();
    let mut expected_entropy = 0.0;
    for cell in 0..2 {
        let n = initial.cell_component_moles(cell).unwrap();
        let a: f64 = (0..3).map(|s| n[s] * intercept[s]).sum();
        let b: f64 = (0..3).map(|s| n[s] * slope[s]).sum();
        let old_t = case.temperatures[cell];
        let heat = if cell == 0 { 75.0 } else { -75.0 };
        let expected_t = ((a + b * old_t).powi(2) + 2.0 * b * heat).sqrt() / b - a / b;
        near(
            result.field().temperatures()[cell].value(),
            expected_t,
            1e-9,
        );
        expected_entropy += a * (expected_t / old_t).ln() + b * (expected_t - old_t);
    }
    assert_eq!(result.field().component_moles(), initial.component_moles());
    near(
        result.step().entropy_change_j_per_k,
        expected_entropy,
        1e-10,
    );
}

#[test]
fn all_faces_use_the_same_old_state_and_pair_every_species_and_energy_transfer() {
    let case = Case {
        volumes: vec![0.001, 0.002, 0.0005],
        fractions: vec![[0.2, 0.3, 0.5], [0.4, 0.4, 0.2], [0.3, 0.2, 0.5]],
        temperatures: vec![300.0, 600.0, 450.0],
        links: vec![[0, 1], [1, 2]],
        conductivity: 0.1,
        ..Case::default()
    };
    let initial = case.field();
    let dt = 0.05;
    let result = initial
        .step_candidate(dt, MulticomponentNonisothermalSolve::default(), || false)
        .unwrap();
    let mut expected = initial.component_moles().to_vec();
    let mut energy_changes = [0.0; 3];
    for (face, &[left, right]) in case.links.iter().enumerate() {
        let (flux, energy, _) = face_oracle(&case, left, right);
        for s in 0..3 {
            near(result.step().transferred_mol[face][s], dt * flux[s], 3e-16);
            expected[3 * left + s] -= dt * flux[s];
            expected[3 * right + s] += dt * flux[s];
        }
        near(result.step().transferred_energy_j[face], dt * energy, 2e-10);
        energy_changes[left] -= dt * energy;
        energy_changes[right] += dt * energy;
    }
    for (&actual, expected) in result.field().component_moles().iter().zip(expected) {
        near(actual, expected, 3e-15);
    }
    for cell in 0..3 {
        let before = initial.cell_component_moles(cell).unwrap();
        let after = result.field().cell_component_moles(cell).unwrap();
        near(
            after.iter().sum(),
            CONCENTRATION * case.volumes[cell],
            3e-15,
        );
        let difference = enthalpy(after, result.field().temperatures()[cell].value(), &case)
            - enthalpy(before, case.temperatures[cell], &case);
        near(difference, energy_changes[cell], 1e-8);
    }
    assert!(
        result
            .step()
            .component_residual_mol
            .iter()
            .all(|r| r.abs() < 1e-13)
    );
    near(result.step().energy_residual_j, 0.0, 1e-8);
    assert!(result.step().entropy_change_j_per_k > 0.0);
    assert_eq!(result.field().elapsed_s(), dt);
}

#[test]
fn enormous_sensible_reference_does_not_erase_small_fourier_heat() {
    let mut case = Case {
        fractions: vec![[0.2, 0.3, 0.5]; 2],
        conductivity: 0.2,
        ..Case::default()
    };
    let dt = 0.01;
    let ordinary = case
        .field()
        .step_candidate(dt, MulticomponentNonisothermalSolve::default(), || false)
        .unwrap();
    case.datum = 1e15;
    let large = case
        .field()
        .step_candidate(dt, MulticomponentNonisothermalSolve::default(), || false)
        .unwrap();
    for cell in 0..2 {
        near(
            large.field().temperatures()[cell].value(),
            ordinary.field().temperatures()[cell].value(),
            1e-9,
        );
    }
    near(large.step().energy_residual_j, 0.0, 1e-8);
    near(
        large.step().entropy_change_j_per_k,
        ordinary.step().entropy_change_j_per_k,
        1e-10,
    );
    near(large.step().transferred_energy_j[0], -0.6, 1e-14);
}

#[test]
fn tiny_positive_trace_and_adjacent_temperatures_keep_finite_independent_inventories() {
    let t = 400.0_f64;
    let case = Case {
        fractions: vec![[1e-280, 0.3, 0.7], [2e-280, 0.3, 0.7]],
        temperatures: vec![t, f64::from_bits(t.to_bits() + 1)],
        ..Case::default()
    };
    let initial = case.field();
    let result = initial
        .step_candidate(0.001, MulticomponentNonisothermalSolve::default(), || false)
        .unwrap();
    assert!(result.step().transferred_mol[0][0] < 0.0);
    assert!(result.field().component_moles()[0] > initial.component_moles()[0]);
    assert!(result.field().component_moles()[3] < initial.component_moles()[3]);
    assert!(
        result
            .field()
            .component_moles()
            .iter()
            .all(|&n| n > 0.0 && n.is_finite())
    );
    assert!(result.step().transferred_energy_j[0].is_finite());
    assert!(
        result
            .field()
            .temperatures()
            .iter()
            .all(|t| t.value().is_finite())
    );
    let trace_before = initial.component_moles()[0] + initial.component_moles()[3];
    let trace_after = result.field().component_moles()[0] + result.field().component_moles()[3];
    near(trace_after / trace_before, 1.0, 3e-15);
}

#[test]
fn qualification_covers_the_whole_temperature_bracket_and_explicit_heat_data() {
    let case = Case {
        b: uniform_matrix(3, 550.0),
        temperatures: vec![300.0, 350.0],
        ..Case::default()
    };
    let (network, initial, caloric) = case.parts();
    // chi=550/300 is globally admissible at the original reference point;
    // chi=550/250 fails the required whole-simplex lower bound at the cold end.
    assert!(
        network
            .bind_nonisothermal_field(
                &initial,
                caloric.clone(),
                vec![Temperature::new(300.0), Temperature::new(350.0)],
                case.faces(),
                [Temperature::new(250.0), Temperature::new(700.0)],
                WORKSPACE,
                || false,
            )
            .is_err()
    );
    let good = network.bind_nonisothermal_field(
        &initial,
        caloric.clone(),
        vec![Temperature::new(300.0), Temperature::new(350.0)],
        case.faces(),
        [Temperature::new(300.0), Temperature::new(700.0)],
        WORKSPACE,
        || false,
    );
    assert!(good.is_ok());
    let mut incomplete = case.faces();
    incomplete[0].heat_of_transport_j_per_mol.pop();
    assert!(
        network
            .bind_nonisothermal_field(
                &initial,
                caloric.clone(),
                vec![Temperature::new(300.0), Temperature::new(350.0)],
                incomplete,
                [Temperature::new(300.0), Temperature::new(700.0)],
                WORKSPACE,
                || false,
            )
            .is_err()
    );
    assert!(
        network
            .bind_nonisothermal_field(
                &initial,
                caloric,
                vec![Temperature::new(300.0), Temperature::new(350.0)],
                case.faces(),
                [Temperature::new(300.0), Temperature::new(700.0)],
                1,
                || false,
            )
            .is_err()
    );
}

#[test]
fn finite_step_entropy_positivity_and_late_cancellation_refuse_atomic_candidates() {
    let heat = Case {
        volumes: vec![0.001; 2],
        fractions: vec![[0.2, 0.3, 0.5]; 2],
        conductivity: 100.0,
        ..Case::default()
    }
    .field();
    // This explicit heat step remains inside [250,700] K but overshoots
    // thermal equilibration far enough to decrease the actual total entropy.
    assert!(matches!(
        heat.step_candidate(0.6, MulticomponentNonisothermalSolve::default(), || false),
        Err(MulticomponentNonisothermalError::Entropy { .. })
    ));
    let initial = Case::default().field();
    let snapshot = initial.clone();
    assert!(
        initial
            .step_candidate(1e6, MulticomponentNonisothermalSolve::default(), || false)
            .is_err()
    );
    let mut polls = 0usize;
    let accepted = initial
        .step_candidate(0.02, MulticomponentNonisothermalSolve::default(), || {
            polls += 1;
            false
        })
        .unwrap();
    assert!(polls > 3);
    let mut visits = 0usize;
    let cancelled =
        initial.step_candidate(0.02, MulticomponentNonisothermalSolve::default(), || {
            visits += 1;
            visits >= polls - 1
        });
    assert!(matches!(
        cancelled,
        Err(MulticomponentNonisothermalError::Cancelled)
    ));
    assert_eq!(initial, snapshot);
    assert_eq!(
        accepted,
        initial
            .step_candidate(0.02, MulticomponentNonisothermalSolve::default(), || false)
            .unwrap()
    );
}

#[test]
fn a_sensible_reference_cannot_turn_sub_ulp_species_transport_into_heat() {
    let initial = Case {
        datum: 1e15,
        temperatures: vec![300.0; 2],
        ..Case::default()
    }
    .field();
    let result =
        initial.step_candidate(1e-17, MulticomponentNonisothermalSolve::default(), || false);
    // Omitting these low-word species would turn reference enthalpy into a
    // resolvable, unphysical 0.0003 K thermal gradient. Refuse explicitly.
    assert!(matches!(
        result,
        Err(MulticomponentNonisothermalError::UnresolvedCaloricInventory { .. })
    ));
    let initial = Case {
        temperatures: vec![300.0; 2],
        ..Case::default()
    }
    .field();
    let candidate = initial
        .step_candidate(1e-17, MulticomponentNonisothermalSolve::default(), || false)
        .unwrap();
    assert_eq!(
        candidate.field().component_moles(),
        initial.component_moles()
    );
    assert!(
        candidate.step().transferred_mol[0]
            .iter()
            .any(|&n| n != 0.0)
    );
    // For ideal constant-Cp liquids at one initial T, transported partial
    // enthalpy exactly matches the inventory change. No temperature changes:
    // an arbitrary source reference must not create a thermal gradient.
    for temperature in candidate.field().temperatures() {
        near(temperature.value(), 300.0, 1e-9);
    }
}
