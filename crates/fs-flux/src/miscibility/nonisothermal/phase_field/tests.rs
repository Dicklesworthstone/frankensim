//! Independent caloric/entropy oracles for thermally evolving diffuse interfaces.

use super::*;
use crate::miscibility::{BinaryMiscibilityCell, BinaryMiscibilityFace};
use fs_qty::AmountConcentration;
use fs_thermochem::solution::caloric::PureLiquidHeatCapacityV1;
use fs_thermochem::solution::{BinaryRegularSolutionV1, gradient::BinarySquareGradientV1};
use fs_thermochem::{MolarThermalQuantityV1, SpeciesId, UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R};

const WORKSPACE: usize = 16 * 1024 * 1024;

#[derive(Clone, Copy)]
pub(super) struct Law {
    pub(super) a: f64,
    pub(super) b: f64,
    pub(super) excess_cp: f64,
    pub(super) cp: [f64; 2],
    pub(super) cp_slope: [f64; 2],
    pub(super) kappa: f64,
}

impl Default for Law {
    fn default() -> Self {
        Self {
            a: 3.0,
            b: 0.0,
            excess_cp: 0.0,
            cp: [40.0; 2],
            cp_slope: [0.0; 2],
            kappa: 200.0,
        }
    }
}

#[derive(Clone, Copy)]
pub(super) struct Link {
    pub(super) cells: [usize; 2],
    pub(super) area: f64,
    pub(super) distance: f64,
    pub(super) diffusivity: f64,
    pub(super) conductivity: f64,
}

pub(super) fn link(i: usize, j: usize, diffusivity: f64, conductivity: f64) -> Link {
    Link {
        cells: [i, j],
        area: 1.0,
        distance: 1.0,
        diffusivity,
        conductivity,
    }
}

fn prepare(
    law: Law,
    cells: &[[f64; 3]], // volume [m³], x_A, T [K]
    links: &[Link],
) -> (
    BinaryMiscibilityNetwork,
    BinaryLiquidCaloricModelV1,
    BinaryMiscibilityState,
) {
    let bracket = [Temperature::new(250.0), Temperature::new(700.0)];
    let species = [SpeciesId::new("A").unwrap(), SpeciesId::new("B").unwrap()];
    let pair = BinaryRegularSolutionV1::new(species.clone(), law.a, law.b, bracket[0], bracket[1])
        .unwrap()
        .with_excess_heat_capacity(law.excess_cp, Temperature::new(300.0))
        .unwrap();
    let pure = core::array::from_fn(|s| {
        PureLiquidHeatCapacityV1::new(
            species[s].clone(),
            Temperature::new(300.0),
            bracket
                .map(|t| {
                    (
                        t,
                        MolarThermalQuantityV1::new(
                            law.cp[s] + law.cp_slope[s] * (t.value() - 300.0),
                        ),
                    )
                })
                .to_vec(),
        )
        .unwrap()
    });
    let caloric = BinaryLiquidCaloricModelV1::new(pair.clone(), pure).unwrap();
    let model = BinarySquareGradientV1::new(
        pair,
        Temperature::new(300.0),
        AmountConcentration::new(1000.0),
        law.kappa,
    )
    .unwrap();
    let network = BinaryMiscibilityNetwork::new(
        model,
        cells
            .iter()
            .map(|c| BinaryMiscibilityCell { volume_m3: c[0] })
            .collect(),
        links
            .iter()
            .map(|l| BinaryMiscibilityFace {
                cells: l.cells,
                area_m2: l.area,
                distance_m: l.distance,
                diffusivity_m2_s: l.diffusivity,
            })
            .collect(),
        WORKSPACE,
        || false,
    )
    .unwrap();
    let state = network
        .state_from_mole_fractions(cells.iter().map(|c| [c[1], 1.0 - c[1]]).collect(), || false)
        .unwrap();
    (network, caloric, state)
}

pub(super) fn fixture(
    law: Law,
    cells: &[[f64; 3]],
    links: &[Link],
) -> BinaryNonisothermalPhaseField {
    let (network, caloric, state) = prepare(law, cells, links);
    network
        .bind_nonisothermal_phase_field(
            &state,
            caloric,
            cells.iter().map(|c| Temperature::new(c[2])).collect(),
            links
                .iter()
                .map(|l| BinaryNonisothermalFace {
                    thermal_conductivity_w_m_k: l.conductivity,
                    heat_of_transport_j_per_mol: 0.0,
                })
                .collect(),
            [Temperature::new(250.0), Temperature::new(700.0)],
            WORKSPACE,
            || false,
        )
        .unwrap()
}

pub(super) fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual {actual:.17e}, expected {expected:.17e}, tolerance {tolerance:e}"
    );
}

// Independent integration of Cp(T), plus the authored regular-solution excess.
fn bulk_enthalpy(law: Law, amounts: [f64; 2], t: f64) -> f64 {
    let delta = t - 300.0;
    let product = amounts[0] * amounts[1] / (amounts[0] + amounts[1]);
    (0..2)
        .map(|s| amounts[s] * (law.cp[s] * delta + 0.5 * law.cp_slope[s] * delta * delta))
        .sum::<f64>()
        + product * (R * law.b + law.excess_cp * delta)
}

fn bulk_entropy(law: Law, amounts: [f64; 2], t: f64) -> f64 {
    let total = amounts[0] + amounts[1];
    let logarithm = (t / 300.0).ln();
    (0..2)
        .map(|s| {
            amounts[s]
                * ((law.cp[s] - 300.0 * law.cp_slope[s]) * logarithm
                    + law.cp_slope[s] * (t - 300.0)
                    - R * (amounts[s] / total).ln())
        })
        .sum::<f64>()
        + amounts[0] * amounts[1] / total * (-R * law.a + law.excess_cp * logarithm)
}

fn interface_shares(law: Law, amounts: &[[f64; 2]], links: &[Link]) -> Vec<f64> {
    let mut shares = vec![0.0; amounts.len()];
    for link in links {
        let [i, j] = link.cells;
        let x = |cell: usize| amounts[cell][0] / (amounts[cell][0] + amounts[cell][1]);
        let difference = x(i) - x(j);
        let half_edge = 0.25 * law.kappa * link.area / link.distance * difference * difference;
        shares[i] += half_edge;
        shares[j] += half_edge;
    }
    shares
}

fn cell_energies(law: Law, field: &BinaryNonisothermalPhaseField, links: &[Link]) -> Vec<f64> {
    interface_shares(law, field.component_moles(), links)
        .into_iter()
        .enumerate()
        .map(|(i, interface)| {
            interface
                + bulk_enthalpy(
                    law,
                    field.component_moles()[i],
                    field.temperatures()[i].value(),
                )
        })
        .collect()
}

fn entropy(law: Law, field: &BinaryNonisothermalPhaseField) -> f64 {
    field
        .component_moles()
        .iter()
        .zip(field.temperatures())
        .map(|(&n, t)| bulk_entropy(law, n, t.value()))
        .sum()
}

fn audit(
    law: Law,
    links: &[Link],
    initial: &BinaryNonisothermalPhaseField,
    candidate: &BinaryNonisothermalPhaseFieldCandidate,
) {
    let next = candidate.field();
    let before = cell_energies(law, initial, links);
    let after = cell_energies(law, next, links);
    let mut heat = vec![0.0; before.len()];
    let mut species = vec![0.0; before.len()];
    for (face, link) in links.iter().enumerate() {
        let [i, j] = link.cells;
        let q = candidate.step().transferred_energy_j[face];
        heat[i] -= q;
        heat[j] += q;
        let n = candidate.step().transferred_a_mol[face];
        species[i] -= n;
        species[j] += n;
    }
    for i in 0..before.len() {
        near(after[i] - before[i], heat[i], 3e-8);
        near(
            next.component_moles()[i][0] - initial.component_moles()[i][0],
            species[i],
            3e-13,
        );
        near(
            next.component_moles()[i][1] - initial.component_moles()[i][1],
            -species[i],
            3e-13,
        );
    }
    near(after.iter().sum::<f64>(), before.iter().sum::<f64>(), 5e-8);
    let actual_entropy = entropy(law, next) - entropy(law, initial);
    assert!(actual_entropy >= -1e-10);
    near(
        candidate.step().entropy_change_j_per_k,
        actual_entropy,
        2e-10,
    );
}

#[test]
fn spinodal_growth_converts_excess_enthalpy_into_heat_and_interface_energy() {
    let law = Law {
        a: 0.0,
        b: 1200.0,
        ..Law::default()
    };
    let links = [link(0, 1, 1e-4, 0.0)];
    let initial = fixture(law, &[[0.001, 0.49, 400.0], [0.001, 0.51, 400.0]], &links);
    let candidate = initial
        .step_candidate(0.2, BinaryNonisothermalSolve::default(), || false)
        .unwrap();
    let next = candidate.field();
    assert!(next.mole_fractions()[0][0] < 0.49);
    assert!(next.mole_fractions()[1][0] > 0.51);
    assert!(next.temperatures().iter().all(|t| t.value() > 400.0));
    let interface_before: f64 = interface_shares(law, initial.component_moles(), &links)
        .iter()
        .sum();
    let interface_after: f64 = interface_shares(law, next.component_moles(), &links)
        .iter()
        .sum();
    assert!(interface_after > interface_before);
    near(
        candidate.step().interfacial_energy_change_j,
        interface_after - interface_before,
        2e-10,
    );
    near(
        candidate.step().homogeneous_enthalpy_change_j,
        -(interface_after - interface_before),
        3e-8,
    );
    audit(law, &links, &initial, &candidate);
}

#[test]
fn a_closed_neighbor_changes_temperature_to_pay_its_half_of_a_changing_interface() {
    let law = Law {
        kappa: 300.0,
        cp: [30.0; 2],
        ..Law::default()
    };
    let links = [
        link(0, 1, 1e-4, 0.0),
        Link {
            area: 0.7,
            distance: 0.4,
            ..link(1, 2, 0.0, 0.0)
        },
    ];
    let initial = fixture(
        law,
        &[
            [0.001, 0.25, 450.0],
            [0.001, 0.55, 450.0],
            [0.001, 0.8, 450.0],
        ],
        &links,
    );
    let candidate = initial
        .step_candidate(0.05, BinaryNonisothermalSolve::default(), || false)
        .unwrap();
    let next = candidate.field();
    assert_eq!(next.component_moles()[2], initial.component_moles()[2]);
    assert_eq!(candidate.step().transferred_a_mol[1], 0.0);
    assert_eq!(candidate.step().transferred_energy_j[1], 0.0);
    assert!(next.temperatures()[2].value() > initial.temperatures()[2].value());
    let initial_share = interface_shares(law, initial.component_moles(), &links)[2];
    let final_share = interface_shares(law, next.component_moles(), &links)[2];
    near(
        30.0 * (next.temperatures()[2].value() - 450.0),
        initial_share - final_share,
        3e-8,
    );
    audit(law, &links, &initial, &candidate);
}

#[test]
fn unequal_temperature_gradient_force_matches_total_entropy_at_fixed_local_energy() {
    // Identical pure calorics and B=C=Q=k=0 make transported enthalpy exactly
    // zero. The chemical force is therefore the derivative at EACH fixed E_i.
    let law = Law {
        kappa: 2000.0,
        ..Law::default()
    };
    let links = [
        link(0, 1, 1e-4, 0.0),
        Link {
            area: 0.7,
            distance: 0.4,
            ..link(1, 2, 0.0, 0.0)
        },
    ];
    let initial = fixture(
        law,
        &[
            [0.001, 0.25, 350.0],
            [0.001, 0.55, 450.0],
            [0.001, 0.8, 550.0],
        ],
        &links,
    );
    let fixed_energy = cell_energies(law, &initial, &links);
    let perturbed_entropy = |q: f64| {
        let mut amounts = initial.component_moles().to_vec();
        amounts[0][0] -= q;
        amounts[0][1] += q;
        amounts[1][0] += q;
        amounts[1][1] -= q;
        let shares = interface_shares(law, &amounts, &links);
        amounts
            .iter()
            .enumerate()
            .map(|(i, &n)| {
                let t = 300.0 + (fixed_energy[i] - shares[i]) / (40.0 * (n[0] + n[1]));
                bulk_entropy(law, n, t)
            })
            .sum::<f64>()
    };
    let delta = 1e-6;
    let entropy_gradient = (perturbed_entropy(delta) - perturbed_entropy(-delta)) / (2.0 * delta);
    let face_x = 0.5 * (0.25 + 0.55);
    let mobility = 1000.0 * 1e-4 * face_x * (1.0 - face_x) / R;
    let expected_rate = mobility * entropy_gradient;
    assert!(expected_rate.abs() > 1e-4);
    let dt = 1e-4;
    let candidate = initial
        .step_candidate(dt, BinaryNonisothermalSolve::default(), || false)
        .unwrap();
    let actual_rate = candidate.step().transferred_a_mol[0] / dt;
    near(actual_rate, expected_rate, expected_rate.abs() * 2e-3);
    audit(law, &links, &initial, &candidate);
}

#[test]
fn fourier_redistribution_uses_excess_heat_capacity_and_converges_to_the_thermal_mode() {
    let law = Law {
        cp: [20.0, 60.0],
        excess_cp: 40.0,
        ..Law::default()
    };
    let links = [link(0, 1, 0.0, 2.0)];
    let initial = fixture(law, &[[0.001, 0.5, 350.0], [0.001, 0.5, 550.0]], &links);
    // Each cell's ORIGINAL Cp is .5*20+.5*60+.25*40 = 50 J/K.
    let exact_left = 450.0 - 100.0 * (-2.0_f64 * 2.0 / 50.0).exp();
    let integrate = |steps: usize| {
        let mut field = initial.clone();
        for _ in 0..steps {
            field = field
                .step_candidate(
                    1.0 / steps as f64,
                    BinaryNonisothermalSolve::default(),
                    || false,
                )
                .unwrap()
                .into_parts()
                .0;
        }
        assert_eq!(field.component_moles(), initial.component_moles());
        near(field.elapsed_s(), 1.0, 1e-15);
        near(
            field.temperatures()[0].value() + field.temperatures()[1].value(),
            900.0,
            2e-9,
        );
        (field.temperatures()[0].value() - exact_left).abs()
    };
    let coarse = integrate(1);
    let medium = integrate(4);
    let fine = integrate(16);
    assert!(medium < 0.4 * coarse, "{coarse:e} -> {medium:e}");
    assert!(fine < 0.4 * medium, "{medium:e} -> {fine:e}");
    assert!(fine < 0.04);
}

#[test]
fn unequal_cells_retain_original_varying_pure_and_negative_excess_calorics() {
    let law = Law {
        a: 2.8,
        b: 300.0,
        excess_cp: -10.0,
        cp: [35.0, 60.0],
        cp_slope: [0.02, 0.04],
        kappa: 500.0,
    };
    let links = [
        Link {
            area: 0.7,
            distance: 0.4,
            ..link(0, 1, 1e-4, 0.2)
        },
        Link {
            area: 1.3,
            distance: 0.7,
            ..link(1, 2, 1e-4, 0.3)
        },
    ];
    let initial = fixture(
        law,
        &[
            [0.0007, 0.2, 340.0],
            [0.0013, 0.65, 520.0],
            [0.0009, 0.4, 460.0],
        ],
        &links,
    );
    let candidate = initial
        .step_candidate(0.001, BinaryNonisothermalSolve::default(), || false)
        .unwrap();
    assert_ne!(
        candidate.field().component_moles(),
        initial.component_moles()
    );
    assert_ne!(candidate.field().temperatures(), initial.temperatures());
    audit(law, &links, &initial, &candidate);
}

#[test]
fn zero_kappa_refuses_a_temperature_bracket_that_crosses_into_nonconvexity() {
    // chi(300 K)=2 admits the reference graph, but chi(250 K)=2.4.
    let law = Law {
        a: 0.0,
        b: 600.0,
        kappa: 0.0,
        ..Law::default()
    };
    let cells = [[0.001, 0.4, 400.0], [0.001, 0.6, 450.0]];
    let links = [link(0, 1, 1e-4, 0.0)];
    let (network, caloric, state) = prepare(law, &cells, &links);
    assert!(
        network
            .bind_nonisothermal_phase_field(
                &state,
                caloric,
                cells.iter().map(|c| Temperature::new(c[2])).collect(),
                vec![BinaryNonisothermalFace {
                    thermal_conductivity_w_m_k: 0.0,
                    heat_of_transport_j_per_mol: 0.0
                }],
                [Temperature::new(250.0), Temperature::new(700.0)],
                WORKSPACE,
                || false,
            )
            .is_err()
    );
    assert_eq!(state.component_moles(), &[[0.4, 0.6], [0.6, 0.4]]);
}

#[test]
fn final_poll_cancellation_preserves_every_state_and_replay_is_exact() {
    let links = [link(0, 1, 1e-4, 0.1), link(1, 2, 1e-4, 0.2)];
    let initial = fixture(
        Law::default(),
        &[
            [0.001, 0.3, 380.0],
            [0.001, 0.5, 430.0],
            [0.001, 0.7, 480.0],
        ],
        &links,
    );
    let before = initial.clone();
    let controls = BinaryNonisothermalSolve::default();
    let mut polls = 0usize;
    let first = initial
        .step_candidate(0.01, controls, || {
            polls += 1;
            false
        })
        .unwrap();
    assert!(polls > 3);
    let mut cancelled_polls = 0usize;
    let failed = initial.step_candidate(0.01, controls, || {
        cancelled_polls += 1;
        cancelled_polls == polls
    });
    assert_eq!(failed.unwrap_err(), BinaryNonisothermalError::Cancelled);
    assert_eq!(initial, before);
    let replay = initial.step_candidate(0.01, controls, || false).unwrap();
    assert_eq!(replay, first);
}
