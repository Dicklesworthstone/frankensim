//! Independent caloric, full-friction and nonlocal interface-energy oracles.

use super::*;
use super::super::tests::{matrix, uniform_matrix};
use crate::miscibility::multicomponent::{
    MulticomponentMiscibilityCell, MulticomponentMiscibilityFace, MulticomponentMiscibilityNetwork,
    nonisothermal::{
        MulticomponentNonisothermalError, MulticomponentNonisothermalFace,
        MulticomponentNonisothermalSolve,
    },
    phase_field::{MulticomponentPhaseField, MulticomponentPhaseFieldNetwork},
};
use fs_qty::{AmountConcentration, MolarMass, Temperature};
use fs_thermochem::solution::caloric::PureLiquidHeatCapacityV1;
use fs_thermochem::solution::multicomponent::{
    MulticomponentRegularSolutionV1, caloric::MulticomponentLiquidCaloricModelV1,
};
use fs_thermochem::transport::maxwell_stefan::MaxwellStefanV1;
use fs_thermochem::{MolarThermalQuantityV1, SpeciesId, UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R};

pub(super) const WORKSPACE: usize = 32 * 1024 * 1024;
pub(super) const CONCENTRATION: f64 = 1000.0;
const DATUM: f64 = 300.0;

pub(super) fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual {actual:.17e}, expected {expected:.17e}, tolerance {tolerance:e}"
    );
}

pub(super) fn link(i: usize, j: usize) -> MulticomponentMiscibilityFace {
    MulticomponentMiscibilityFace {
        cells: [i, j],
        area_m2: 1.0,
        distance_m: 1.0,
        diffusivity_matrix: 0,
    }
}

#[derive(Clone)]
pub(super) struct Case {
    pub(super) a: Vec<f64>,
    pub(super) b: Vec<f64>,
    pub(super) cp: Vec<f64>,
    pub(super) cp_slope: Vec<f64>,
    pub(super) reference_temperature: f64,
    pub(super) volumes: Vec<f64>,
    pub(super) fractions: Vec<Vec<f64>>,
    pub(super) temperatures: Vec<f64>,
    pub(super) links: Vec<MulticomponentMiscibilityFace>,
    pub(super) diffusivities: Vec<Vec<f64>>,
    pub(super) conductivity: f64,
    pub(super) heat: Vec<f64>,
    pub(super) kappa: f64,
}

impl Default for Case {
    fn default() -> Self {
        Self {
            a: uniform_matrix(3, 4.0),
            b: vec![0.0; 9],
            cp: vec![25.0, 40.0, 70.0],
            cp_slope: vec![0.02, 0.05, 0.01],
            reference_temperature: 300.0,
            volumes: vec![0.001, 0.0015],
            fractions: vec![vec![0.2, 0.3, 0.5], vec![0.4, 0.4, 0.2]],
            temperatures: vec![350.0, 500.0],
            links: vec![link(0, 1)],
            diffusivities: vec![uniform_matrix(3, 1e-4)],
            conductivity: 0.2,
            heat: vec![0.0; 3],
            kappa: 2000.0,
        }
    }
}

impl Case {
    fn parts(
        &self,
    ) -> (
        MulticomponentPhaseFieldNetwork,
        MulticomponentPhaseField,
        MulticomponentLiquidCaloricModelV1,
    ) {
        let n = self.cp.len();
        let species = (0..n)
            .map(|s| SpeciesId::new(format!("S{s:03}")).unwrap())
            .collect::<Vec<_>>();
        let law = MulticomponentRegularSolutionV1::new(
            species.clone(),
            self.a.clone(),
            self.b.clone(),
            Temperature::new(100.0),
            Temperature::new(1000.0),
        )
        .unwrap();
        let caloric = MulticomponentLiquidCaloricModelV1::new(
            law.clone(),
            species
                .iter()
                .enumerate()
                .map(|(s, id)| {
                    PureLiquidHeatCapacityV1::new(
                        id.clone(),
                        Temperature::new(DATUM),
                        [100.0, 1000.0]
                            .into_iter()
                            .map(|temperature| {
                                (
                                    Temperature::new(temperature),
                                    MolarThermalQuantityV1::new(
                                        self.cp[s] + self.cp_slope[s] * (temperature - DATUM),
                                    ),
                                )
                            })
                            .collect(),
                    )
                    .unwrap()
                })
                .collect(),
        )
        .unwrap();
        let transport = MaxwellStefanV1::new(
            species
                .into_iter()
                .enumerate()
                .map(|(s, id)| (id, MolarMass::new(0.018 + 0.01 * s as f64)))
                .collect(),
        )
        .unwrap();
        let network = MulticomponentPhaseFieldNetwork::new(
            law,
            transport,
            Temperature::new(self.reference_temperature),
            AmountConcentration::new(CONCENTRATION),
            self.volumes
                .iter()
                .map(|&volume_m3| MulticomponentMiscibilityCell { volume_m3 })
                .collect(),
            self.links.clone(),
            self.diffusivities.clone(),
            self.kappa,
            WORKSPACE,
            || false,
        )
        .unwrap();
        let initial = network
            .bind(
                self.fractions
                    .iter()
                    .zip(network.capacities_mol())
                    .map(|(x, &capacity)| x.iter().map(|&x| capacity * x).collect())
                    .collect(),
                || false,
            )
            .unwrap();
        (network, initial, caloric)
    }

    fn thermal_faces(&self) -> Vec<MulticomponentNonisothermalFace> {
        self.links
            .iter()
            .map(|_| MulticomponentNonisothermalFace {
                thermal_conductivity_w_m_k: self.conductivity,
                heat_of_transport_j_per_mol: self.heat.clone(),
            })
            .collect()
    }

    fn field_with_budget(
        &self,
        maximum_workspace_bytes: usize,
    ) -> Result<MulticomponentNonisothermalPhaseField> {
        let (network, initial, caloric) = self.parts();
        network.bind_nonisothermal_phase_field(
            &initial,
            caloric,
            self.temperatures.iter().copied().map(Temperature::new).collect(),
            self.thermal_faces(),
            [Temperature::new(250.0), Temperature::new(700.0)],
            maximum_workspace_bytes,
            || false,
        )
    }

    pub(super) fn field(&self) -> MulticomponentNonisothermalPhaseField {
        self.field_with_budget(WORKSPACE).unwrap()
    }
}

// Independent integrations of the authored piecewise-linear Cp laws. The
// chosen fixtures remain inside a single linear segment, so these are exact.
pub(super) fn enthalpy(case: &Case, amounts: &[f64], temperature: f64) -> f64 {
    let n = amounts.len();
    let total: f64 = amounts.iter().sum();
    let delta = temperature - DATUM;
    let mut value = amounts
        .iter()
        .enumerate()
        .map(|(s, &amount)| {
            amount * (case.cp[s] * delta + 0.5 * case.cp_slope[s] * delta * delta)
        })
        .sum::<f64>();
    for s in 0..n {
        for t in s + 1..n {
            value += R * case.b[s * n + t] * amounts[s] * amounts[t] / total;
        }
    }
    value
}

pub(super) fn entropy(case: &Case, amounts: &[f64], temperature: f64) -> f64 {
    let n = amounts.len();
    let total: f64 = amounts.iter().sum();
    let logarithm = (temperature / DATUM).ln();
    let mut value = amounts
        .iter()
        .enumerate()
        .map(|(s, &amount)| {
            amount
                * ((case.cp[s] - DATUM * case.cp_slope[s]) * logarithm
                    + case.cp_slope[s] * (temperature - DATUM)
                    - R * (amount / total).ln())
        })
        .sum::<f64>();
    for s in 0..n {
        for t in s + 1..n {
            value -= R * case.a[s * n + t] * amounts[s] * amounts[t] / total;
        }
    }
    value
}

pub(super) fn rows(field: &MulticomponentNonisothermalPhaseField) -> Vec<Vec<f64>> {
    field
        .component_moles()
        .chunks_exact(field.network().species_count())
        .map(<[f64]>::to_vec)
        .collect()
}

pub(super) fn interface_shares(case: &Case, amounts: &[Vec<f64>]) -> Vec<f64> {
    let fractions = amounts
        .iter()
        .map(|row| {
            let total: f64 = row.iter().sum();
            row.iter().map(|amount| amount / total).collect::<Vec<_>>()
        })
        .collect::<Vec<_>>();
    let mut shares = vec![0.0; amounts.len()];
    for face in &case.links {
        let [i, j] = face.cells;
        let squared_jump = fractions[i]
            .iter()
            .zip(&fractions[j])
            .map(|(left, right)| (right - left) * (right - left))
            .sum::<f64>();
        let half_edge = case.kappa * face.area_m2 * squared_jump / (8.0 * face.distance_m);
        shares[i] += half_edge;
        shares[j] += half_edge;
    }
    shares
}

pub(super) fn local_energies(
    case: &Case,
    field: &MulticomponentNonisothermalPhaseField,
) -> Vec<f64> {
    let amounts = rows(field);
    interface_shares(case, &amounts)
        .into_iter()
        .enumerate()
        .map(|(cell, share)| {
            enthalpy(case, &amounts[cell], field.temperatures()[cell].value()) + share
        })
        .collect()
}

fn recover_temperature(case: &Case, amounts: &[f64], target_enthalpy: f64) -> f64 {
    let capacity = amounts.iter().enumerate().map(|(s, n)| n * case.cp[s]).sum::<f64>();
    let slope = amounts
        .iter()
        .enumerate()
        .map(|(s, n)| n * case.cp_slope[s])
        .sum::<f64>();
    let sensible = target_enthalpy - enthalpy(case, amounts, DATUM);
    // Rationalized positive-Cp root of C*delta + B*delta^2/2 = sensible.
    let delta = if slope == 0.0 {
        sensible / capacity
    } else {
        2.0 * sensible / (capacity + (capacity * capacity + 2.0 * slope * sensible).sqrt())
    };
    DATUM + delta
}

pub(super) fn audit(
    case: &Case,
    initial: &MulticomponentNonisothermalPhaseField,
    candidate: &MulticomponentNonisothermalPhaseFieldCandidate,
) {
    let next = candidate.field();
    let n = initial.network().species_count();
    let before = rows(initial);
    let after = rows(next);
    let old_shares = interface_shares(case, &before);
    let new_shares = interface_shares(case, &after);
    let mut energy_receipts = vec![0.0; before.len()];
    let mut species_receipts = vec![vec![0.0; n]; before.len()];
    for (face_index, face) in initial.network().faces().iter().enumerate() {
        let [i, j] = face.cells;
        energy_receipts[i] -= candidate.step().transferred_energy_j[face_index];
        energy_receipts[j] += candidate.step().transferred_energy_j[face_index];
        for s in 0..n {
            species_receipts[i][s] -= candidate.step().transferred_mol[face_index][s];
            species_receipts[j][s] += candidate.step().transferred_mol[face_index][s];
        }
    }
    let mut enthalpy_change = 0.0;
    let mut interface_change = 0.0;
    let mut entropy_change = 0.0;
    for cell in 0..before.len() {
        let old_t = initial.temperatures()[cell].value();
        let new_t = next.temperatures()[cell].value();
        let old_h = enthalpy(case, &before[cell], old_t);
        let new_h = enthalpy(case, &after[cell], new_t);
        let change_e = new_shares[cell] - old_shares[cell];
        near(new_h - old_h + change_e, energy_receipts[cell], 5e-8);
        near(
            new_t,
            recover_temperature(
                case,
                &after[cell],
                old_h + old_shares[cell] + energy_receipts[cell] - new_shares[cell],
            ),
            3e-9,
        );
        for s in 0..n {
            near(after[cell][s] - before[cell][s], species_receipts[cell][s], 3e-13);
            assert!(after[cell][s] > 0.0);
        }
        near(after[cell].iter().sum(), initial.network().capacities_mol()[cell], 1e-12);
        enthalpy_change += new_h - old_h;
        interface_change += change_e;
        entropy_change += entropy(case, &after[cell], new_t) - entropy(case, &before[cell], old_t);
    }
    for s in 0..n {
        near(
            after.iter().map(|row| row[s]).sum(),
            before.iter().map(|row| row[s]).sum(),
            1e-12,
        );
    }
    near(candidate.step().homogeneous_enthalpy_change_j, enthalpy_change, 5e-8);
    near(candidate.step().interfacial_energy_change_j, interface_change, 1e-9);
    near(candidate.step().energy_residual_j, enthalpy_change + interface_change, 5e-8);
    near(enthalpy_change + interface_change, 0.0, 5e-8);
    near(candidate.step().entropy_change_j_per_k, entropy_change, 5e-10);
    assert!(entropy_change >= -5e-10);
}

pub(super) fn spinodal_case() -> Case {
    let amplitude = 0.01;
    Case {
        a: vec![0.0; 9],
        b: uniform_matrix(3, 1600.0),
        volumes: vec![0.001; 2],
        fractions: vec![
            vec![1.0 / 3.0 - amplitude, 1.0 / 3.0 + amplitude, 1.0 / 3.0],
            vec![1.0 / 3.0 + amplitude, 1.0 / 3.0 - amplitude, 1.0 / 3.0],
        ],
        temperatures: vec![400.0; 2],
        conductivity: 0.0,
        kappa: 200.0,
        ..Case::default()
    }
}

#[test]
fn ternary_spinodal_growth_pays_interface_energy_and_heats_original_variable_cp_liquids() {
    let case = spinodal_case();
    let initial = case.field();
    let candidate = initial
        .step_candidate(0.2, MulticomponentNonisothermalSolve::default(), || false)
        .unwrap();
    let next = candidate.field();
    assert!(
        next.cell_mole_fractions(0).unwrap()[0]
            < initial.cell_mole_fractions(0).unwrap()[0]
    );
    assert!(
        next.cell_mole_fractions(1).unwrap()[0]
            > initial.cell_mole_fractions(1).unwrap()[0]
    );
    assert!(next.temperatures().iter().all(|temperature| temperature.value() > 400.0));
    assert!(candidate.step().interfacial_energy_change_j > 0.0);
    assert!(candidate.step().homogeneous_enthalpy_change_j < 0.0);
    assert!(candidate.step().entropy_change_j_per_k > 0.0);
    audit(&case, &initial, &candidate);
}

#[test]
fn an_unchanged_composition_neighbor_heats_or_cools_to_pay_its_interface_share() {
    // Strictly positive molecular diffusivities are retained on both faces.
    // Choose x0 analytically so the initial total force on face1-2 vanishes,
    // while face0-1 moves and changes the spectator's interface allocation.
    for (kappa, heating) in [(2000.0, false), (4000.0, true)] {
        let temperature: f64 = 450.0;
        let center: f64 = 1.0 / 3.0;
        let amplitude: f64 = 0.02;
        let p1 = center + amplitude;
        let p2 = center - amplitude;
        let reduced_gradient = kappa / (2.0 * R * temperature);
        let bulk_jump = (p2 / p1).ln() + 8.0 * amplitude;
        let p0 = center + 5.0 * amplitude - bulk_jump / reduced_gradient;
        let case = Case {
            cp: vec![30.0; 3],
            cp_slope: vec![0.0; 3],
            volumes: vec![0.001; 3],
            fractions: [p0, p1, p2]
                .into_iter()
                .map(|p| vec![p, 2.0 * center - p, center])
                .collect(),
            temperatures: vec![temperature; 3],
            links: vec![link(0, 1), link(1, 2)],
            conductivity: 0.0,
            kappa,
            ..Case::default()
        };
        let initial = case.field();
        let candidate = initial
            .step_candidate(0.05, MulticomponentNonisothermalSolve::default(), || false)
            .unwrap();
        let next = candidate.field();
        assert_eq!(
            next.cell_component_moles(2).unwrap(),
            initial.cell_component_moles(2).unwrap(),
        );
        for &amount in &candidate.step().transferred_mol[1] {
            near(amount, 0.0, 1e-16);
        }
        assert!(candidate.step().transferred_energy_j.iter().all(|&energy| energy == 0.0));
        let old_share = interface_shares(&case, &rows(&initial))[2];
        let new_share = interface_shares(&case, &rows(next))[2];
        let temperature_change = next.temperatures()[2].value() - temperature;
        assert_eq!(temperature_change > 0.0, heating);
        assert!(temperature_change.abs() > 1e-6);
        near(30.0 * temperature_change, old_share - new_share, 3e-8);
        audit(&case, &initial, &candidate);
    }
}

fn gradient_entropy_potentials(
    case: &Case,
    field: &MulticomponentNonisothermalPhaseField,
) -> Vec<Vec<f64>> {
    let n = field.network().species_count();
    let mut potentials = vec![vec![0.0; n]; field.network().cells().len()];
    for face in field.network().faces() {
        let [i, j] = face.cells;
        let left = field.cell_mole_fractions(i).unwrap();
        let right = field.cell_mole_fractions(j).unwrap();
        let left_moles: f64 = field.cell_component_moles(i).unwrap().iter().sum();
        let right_moles: f64 = field.cell_component_moles(j).unwrap().iter().sum();
        let beta_sum = 1.0 / field.temperatures()[i].value()
            + 1.0 / field.temperatures()[j].value();
        let coefficient = case.kappa * face.area_m2 * beta_sum / (4.0 * face.distance_m);
        for s in 0..n {
            potentials[i][s] += coefficient * (left[s] - right[s]) / left_moles;
            potentials[j][s] -= coefficient * (left[s] - right[s]) / right_moles;
        }
    }
    potentials
}

fn transported_enthalpy(case: &Case, mean: &[f64], left_t: f64, right_t: f64, s: usize) -> f64 {
    let thermal_mean = if left_t == right_t {
        left_t
    } else {
        left_t * right_t * (right_t / left_t).ln() / (right_t - left_t)
    };
    let pure = case.cp[s] * (thermal_mean - DATUM)
        + 0.5 * case.cp_slope[s]
            * (left_t * right_t - 2.0 * DATUM * thermal_mean + DATUM * DATUM);
    // The omitted excess partial-enthalpy term is common to all species and
    // cancels in the zero-net-molar frame.
    pure + R * mean
        .iter()
        .enumerate()
        .map(|(t, &fraction)| case.b[s * mean.len() + t] * fraction)
        .sum::<f64>()
}

fn ternary_face_oracle(
    case: &Case,
    field: &MulticomponentNonisothermalPhaseField,
) -> ([f64; 3], f64, f64) {
    let face = field.network().faces()[0];
    let [i, j] = face.cells;
    let left = field.cell_mole_fractions(i).unwrap();
    let right = field.cell_mole_fractions(j).unwrap();
    let left_t = field.temperatures()[i].value();
    let right_t = field.temperatures()[j].value();
    let beta_jump = 1.0 / right_t - 1.0 / left_t;
    let beta_mean = 0.5 * (1.0 / left_t + 1.0 / right_t);
    let mean: [f64; 3] = core::array::from_fn(|s| 0.5 * (left[s] + right[s]));
    let gradient = gradient_entropy_potentials(case, field);
    let effort: [f64; 3] = core::array::from_fn(|s| {
        (right[s] / left[s]).ln()
            + (0..3)
                .map(|t| (case.a[s * 3 + t] + case.b[s * 3 + t] * beta_mean)
                    * (right[t] - left[t]))
                .sum::<f64>()
            + (gradient[j][s] - gradient[i][s] - case.heat[s] * beta_jump) / R
    });
    let common = (0..3).map(|s| mean[s] * effort[s]).sum::<f64>();
    let driving: [f64; 3] =
        core::array::from_fn(|s| mean[s] * (effort[s] - common) / face.distance_m);
    let d = &case.diffusivities[face.diffusivity_matrix];
    let aa = mean[1] / d[1] + (mean[0] + mean[2]) / d[2];
    let ab = mean[0] * (1.0 / d[2] - 1.0 / d[1]);
    let ba = mean[1] * (1.0 / d[5] - 1.0 / d[1]);
    let bb = mean[0] / d[1] + (mean[1] + mean[2]) / d[5];
    let determinant = aa * bb - ab * ba;
    let ja = -CONCENTRATION * (bb * driving[0] - ab * driving[1]) / determinant;
    let jb = -CONCENTRATION * (aa * driving[1] - ba * driving[0]) / determinant;
    let flux = [ja, jb, -ja - jb];
    let energy = (0..3)
        .map(|s| (transported_enthalpy(case, &mean, left_t, right_t, s) + case.heat[s]) * flux[s])
        .sum::<f64>()
        - case.conductivity * (right_t - left_t) / face.distance_m;
    let production = -R * (0..3).map(|s| flux[s] * effort[s]).sum::<f64>()
        + case.conductivity * (right_t - left_t).powi(2)
            / (face.distance_m * left_t * right_t);
    (flux, energy, production)
}

#[test]
fn unequal_temperature_full_friction_and_reciprocal_heat_match_total_entropy_direction() {
    let case = Case {
        a: matrix(3, &[(0, 1, 3.2), (0, 2, 2.7), (1, 2, 3.4)]),
        b: matrix(3, &[(0, 1, 700.0), (0, 2, -100.0), (1, 2, 400.0)]),
        volumes: vec![0.0007, 0.0013],
        fractions: vec![vec![0.25, 0.30, 0.45], vec![0.42, 0.33, 0.25]],
        temperatures: vec![340.0, 520.0],
        links: vec![MulticomponentMiscibilityFace {
            area_m2: 0.7,
            distance_m: 0.4,
            ..link(0, 1)
        }],
        diffusivities: vec![matrix(3, &[(0, 1, 1e-5), (0, 2, 3e-5), (1, 2, 7e-5)])],
        heat: vec![1500.0, -800.0, 400.0],
        kappa: 400.0,
        ..Case::default()
    };
    let initial = case.field();
    let (flux, energy, production) = ternary_face_oracle(&case, &initial);
    let area = case.links[0].area_m2;
    let duration = 0.005;
    let candidate = initial
        .step_candidate(duration, MulticomponentNonisothermalSolve::default(), || false)
        .unwrap();
    for s in 0..3 {
        near(candidate.step().transferred_mol[0][s], duration * area * flux[s], 3e-14);
    }
    near(candidate.step().transferred_energy_j[0], duration * area * energy, 3e-10);
    near(candidate.step().entropy_production_rate_w_per_k, area * production, 3e-11);

    // Differentiate the independent total entropy while changing each local
    // U=H+half-incident-E and every original inventory along the predicted
    // physical direction. This validates the temperature-weighted nonlocal
    // gradient effort without invoking a production potential evaluator.
    let original_energy = local_energies(&case, &initial);
    let perturbed_entropy = |dt: f64| {
        let mut amounts = rows(&initial);
        for s in 0..3 {
            amounts[0][s] -= dt * area * flux[s];
            amounts[1][s] += dt * area * flux[s];
        }
        let shares = interface_shares(&case, &amounts);
        (0..2)
            .map(|cell| {
                let direction = if cell == 0 { -area * energy } else { area * energy };
                let target_h = original_energy[cell] + dt * direction - shares[cell];
                let t = recover_temperature(&case, &amounts[cell], target_h);
                entropy(&case, &amounts[cell], t)
            })
            .sum::<f64>()
    };
    let epsilon = 1e-3;
    let entropy_direction =
        (perturbed_entropy(epsilon) - perturbed_entropy(-epsilon)) / (2.0 * epsilon);
    near(entropy_direction, area * production, 5e-8);
    assert!(production > 0.0);
    audit(&case, &initial, &candidate);

    let mut shifted = case.clone();
    for heat in &mut shifted.heat {
        *heat += 1234.0;
    }
    let shifted_candidate = shifted
        .field()
        .step_candidate(duration, MulticomponentNonisothermalSolve::default(), || false)
        .unwrap();
    for s in 0..3 {
        near(
            shifted_candidate.step().transferred_mol[0][s],
            candidate.step().transferred_mol[0][s],
            3e-14,
        );
    }
    near(
        shifted_candidate.step().transferred_energy_j[0],
        candidate.step().transferred_energy_j[0],
        3e-10,
    );
    for cell in 0..2 {
        near(
            shifted_candidate.field().temperatures()[cell].value(),
            candidate.field().temperatures()[cell].value(),
            3e-9,
        );
    }
}

#[test]
fn binary_limit_preserves_original_kappa_and_endpoint_inverse_temperature_weights() {
    let case = Case {
        a: uniform_matrix(2, 2.5),
        b: uniform_matrix(2, 400.0),
        cp: vec![35.0, 60.0],
        cp_slope: vec![0.02, 0.04],
        volumes: vec![0.0009, 0.0018],
        fractions: vec![vec![0.31, 0.69], vec![0.64, 0.36]],
        temperatures: vec![350.0, 520.0],
        links: vec![MulticomponentMiscibilityFace {
            area_m2: 0.8,
            distance_m: 0.6,
            ..link(0, 1)
        }],
        diffusivities: vec![uniform_matrix(2, 7e-5)],
        conductivity: 0.15,
        heat: vec![1300.0, -500.0],
        kappa: 750.0,
        ..Case::default()
    };
    let initial = case.field();
    let face = case.links[0];
    let [left, right] = [0.31_f64, 0.64_f64];
    let [tl, tr] = [350.0, 520.0];
    let mean = 0.5 * (left + right);
    let beta_sum = 1.0 / tl + 1.0 / tr;
    let beta_jump = 1.0 / tr - 1.0 / tl;
    let chi = 2.5 + 200.0 * beta_sum;
    let gradient_exchange = case.kappa * face.area_m2 * beta_sum * (right - left)
        / (2.0 * face.distance_m)
        * (1.0 / initial.network().capacities_mol()[0]
            + 1.0 / initial.network().capacities_mol()[1]);
    let effort = (right / (1.0 - right)).ln() - (left / (1.0 - left)).ln()
        - 2.0 * chi * (right - left)
        + (gradient_exchange - (case.heat[0] - case.heat[1]) * beta_jump) / R;
    let flux = -CONCENTRATION * 7e-5 * mean * (1.0 - mean) * effort / face.distance_m;
    let carried = transported_enthalpy(&case, &[mean, 1.0 - mean], tl, tr, 0)
        - transported_enthalpy(&case, &[mean, 1.0 - mean], tl, tr, 1)
        + case.heat[0] - case.heat[1];
    let energy = carried * flux - case.conductivity * (tr - tl) / face.distance_m;
    let duration = 0.003;
    let candidate = initial
        .step_candidate(duration, MulticomponentNonisothermalSolve::default(), || false)
        .unwrap();
    near(candidate.step().transferred_mol[0][0], duration * face.area_m2 * flux, 3e-14);
    near(candidate.step().transferred_mol[0][1], -duration * face.area_m2 * flux, 3e-14);
    near(candidate.step().transferred_energy_j[0], duration * face.area_m2 * energy, 3e-10);
    near(
        interface_shares(&case, &rows(&initial)).iter().sum(),
        case.kappa * face.area_m2 * (right - left).powi(2) / (2.0 * face.distance_m),
        1e-10,
    );
    audit(&case, &initial, &candidate);
}

#[test]
fn every_species_and_total_local_energy_are_conserved_through_the_species_bound() {
    for n in [4, 17, 128] {
        let fractions = [1.0, -1.0]
            .into_iter()
            .map(|sign| {
                (0..n)
                    .map(|s| {
                        (1.0 + sign * 0.01
                            * (2.0 * std::f64::consts::PI * s as f64 / n as f64).cos())
                            / n as f64
                    })
                    .collect()
            })
            .collect();
        let case = Case {
            a: uniform_matrix(n, n as f64 + 1.0),
            b: uniform_matrix(n, 200.0),
            cp: vec![30.0; n],
            cp_slope: vec![0.02; n],
            volumes: vec![0.001, 0.0013],
            fractions,
            temperatures: vec![400.0, 420.0],
            diffusivities: vec![uniform_matrix(n, 1e-5)],
            conductivity: 0.05,
            heat: vec![0.0; n],
            kappa: 150.0,
            ..Case::default()
        };
        let initial = case.field();
        let candidate = initial
            .step_candidate(0.01, MulticomponentNonisothermalSolve::default(), || false)
            .unwrap();
        assert_eq!(candidate.step().transferred_mol[0].len(), n);
        assert_eq!(candidate.step().component_residual_mol.len(), n);
        assert!(candidate.step().component_residual_mol.iter().all(|r| r.abs() <= 1e-12));
        near(candidate.step().transferred_mol[0].iter().sum(), 0.0, 1e-12);
        audit(&case, &initial, &candidate);
    }
}

#[test]
fn thermal_bracket_gate_stays_isolated_and_refusals_preserve_all_nonlocal_state() {
    // The reference graph is convex at 500K (chi=1.8), but the thermal bracket
    // contains spinodal liquid at 250K (chi=3.6). Only the interface owner may
    // admit that bracket.
    let case = Case {
        a: vec![0.0; 9],
        b: uniform_matrix(3, 900.0),
        reference_temperature: 500.0,
        temperatures: vec![400.0, 500.0],
        ..Case::default()
    };
    let initial = case.field();
    let phase_graph = initial.network();
    let ordinary = MulticomponentMiscibilityNetwork::new(
        phase_graph.model().clone(),
        phase_graph.transport().clone(),
        phase_graph.temperature(),
        phase_graph.concentration(),
        phase_graph.cells().to_vec(),
        phase_graph.faces().to_vec(),
        case.diffusivities.clone(),
        WORKSPACE,
        || false,
    )
    .unwrap();
    let ordinary_initial = ordinary.bind(rows(&initial), || false).unwrap();
    assert!(
        ordinary
            .bind_nonisothermal_field(
                &ordinary_initial,
                initial.caloric_model().clone(),
                initial.temperatures().to_vec(),
                case.thermal_faces(),
                [Temperature::new(250.0), Temperature::new(700.0)],
                WORKSPACE,
                || false,
            )
            .is_err()
    );
    assert!(case.field_with_budget(initial.workspace_bytes() - 1).is_err());
    assert!(case.field_with_budget(initial.workspace_bytes()).is_ok());
    let snapshot = initial.clone();
    let controls = MulticomponentNonisothermalSolve::default();
    let mut polls = 0;
    let accepted = initial
        .step_candidate(0.01, controls, || {
            polls += 1;
            false
        })
        .unwrap();
    assert!(polls > 8);
    assert_ne!(accepted.field().temperatures(), initial.temperatures());
    let mut cancelled_polls = 0;
    assert!(matches!(
        initial.step_candidate(0.01, controls, || {
            cancelled_polls += 1;
            cancelled_polls == polls
        }),
        Err(MulticomponentNonisothermalError::Cancelled)
    ));
    assert_eq!(initial, snapshot);
    assert_eq!(
        accepted,
        initial.step_candidate(0.01, controls, || false).unwrap()
    );
    assert!(
        initial
            .step_candidate(
                0.01,
                MulticomponentNonisothermalSolve {
                    maximum_caloric_iterations: 0,
                    ..controls
                },
                || false,
            )
            .is_err()
    );
    assert_eq!(initial, snapshot);
}
