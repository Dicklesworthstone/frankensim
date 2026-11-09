//! Independent finite-jump friction, cross-diffusion and conservative FV oracles.

use super::*;
use fs_qty::{AmountConcentration, MolarMass, Temperature};
use fs_thermochem::SpeciesId;
use fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R;
use fs_thermochem::solution::multicomponent::MulticomponentRegularSolutionV1;
use fs_thermochem::transport::maxwell_stefan::MaxwellStefanV1;

const WORKSPACE: usize = 16 * 1024 * 1024;
const CONCENTRATION: f64 = 1000.0;
const TEMPERATURE: f64 = 300.0;

pub(super) fn matrix(n: usize, pairs: &[(usize, usize, f64)]) -> Vec<f64> {
    let mut matrix = vec![0.0; n * n];
    for &(i, j, value) in pairs {
        matrix[i * n + j] = value;
        matrix[j * n + i] = value;
    }
    matrix
}

pub(super) fn uniform_matrix(n: usize, value: f64) -> Vec<f64> {
    let mut matrix = vec![value; n * n];
    for i in 0..n {
        matrix[i * n + i] = 0.0;
    }
    matrix
}

fn ids(n: usize) -> Vec<SpeciesId> {
    (0..n)
        .map(|i| SpeciesId::new(format!("S{i:03}")).unwrap())
        .collect()
}

fn model(n: usize, chi: Vec<f64>) -> MulticomponentRegularSolutionV1 {
    MulticomponentRegularSolutionV1::new(
        ids(n),
        chi,
        vec![0.0; n * n],
        Temperature::new(250.0),
        Temperature::new(700.0),
    )
    .unwrap()
}

fn transport(masses: &[f64]) -> MaxwellStefanV1 {
    MaxwellStefanV1::new(
        ids(masses.len())
            .into_iter()
            .zip(masses)
            .map(|(id, &m)| (id, MolarMass::new(m)))
            .collect(),
    )
    .unwrap()
}

pub(super) fn face(i: usize, j: usize) -> MulticomponentMiscibilityFace {
    MulticomponentMiscibilityFace {
        cells: [i, j],
        area_m2: 1.0,
        distance_m: 1.0,
        diffusivity_matrix: 0,
    }
}

pub(super) fn network(
    masses: &[f64],
    chi: Vec<f64>,
    volumes: &[f64],
    faces: Vec<MulticomponentMiscibilityFace>,
    diffusivity_matrices: Vec<Vec<f64>>,
) -> MulticomponentMiscibilityNetwork {
    MulticomponentMiscibilityNetwork::new(
        model(masses.len(), chi),
        transport(masses),
        Temperature::new(TEMPERATURE),
        AmountConcentration::new(CONCENTRATION),
        volumes
            .iter()
            .map(|&volume_m3| MulticomponentMiscibilityCell { volume_m3 })
            .collect(),
        faces,
        diffusivity_matrices,
        WORKSPACE,
        || false,
    )
    .unwrap()
}

pub(super) fn field(
    network: &MulticomponentMiscibilityNetwork,
    fractions: &[&[f64]],
) -> MulticomponentMiscibilityField {
    network
        .bind(
            fractions
                .iter()
                .zip(network.capacities_mol())
                .map(|(x, n)| x.iter().map(|x| x * n).collect())
                .collect(),
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

pub(super) fn gibbs_energy(field: &MulticomponentMiscibilityField, n: usize, chi: &[f64]) -> f64 {
    field
        .component_moles()
        .chunks_exact(n)
        .map(|amounts| {
            let total: f64 = amounts.iter().sum();
            let mut reduced = amounts.iter().map(|a| a * (a / total).ln()).sum::<f64>();
            for i in 0..n {
                for j in i + 1..n {
                    reduced += chi[i * n + j] * amounts[i] * amounts[j] / total;
                }
            }
            R * TEMPERATURE * reduced
        })
        .sum()
}

// The common excess-Gibbs term in each chemical potential is a gauge: its
// finite jump cancels in the balanced driving vector and in sum(J_i mu_i).
fn chemical_jumps(left: [f64; 3], right: [f64; 3], chi: &[f64]) -> [f64; 3] {
    core::array::from_fn(|i| {
        (right[i] / left[i]).ln()
            + (0..3)
                .map(|j| chi[3 * i + j] * (right[j] - left[j]))
                .sum::<f64>()
    })
}

// Eliminate J_C=-J_A-J_B in the ORIGINAL friction equations, then solve the
// resulting 2x2 system by its determinant. This oracle calls no MS evaluator.
fn ternary_flux(
    left: [f64; 3],
    right: [f64; 3],
    chi: &[f64],
    d: &[f64],
    distance: f64,
) -> [f64; 3] {
    let x: [f64; 3] = core::array::from_fn(|s| 0.5 * left[s] + 0.5 * right[s]);
    let mu = chemical_jumps(left, right, chi);
    let common: f64 = (0..3).map(|s| x[s] * mu[s]).sum();
    let driving: [f64; 3] = core::array::from_fn(|s| x[s] * (mu[s] - common) / distance);
    let [ab, ac, bc] = [d[1], d[2], d[5]];
    let a = x[1] / ab + (x[0] + x[2]) / ac;
    let b = x[0] * (1.0 / ac - 1.0 / ab);
    let c = x[1] * (1.0 / bc - 1.0 / ab);
    let e = x[0] / ab + (x[1] + x[2]) / bc;
    let determinant = a * e - b * c;
    let ja = CONCENTRATION * (-driving[0] * e + driving[1] * b) / determinant;
    let jb = CONCENTRATION * (-a * driving[1] + c * driving[0]) / determinant;
    [ja, jb, -ja - jb]
}

#[test]
fn binary_reduction_uses_exact_nonideal_chemical_jumps() {
    for interaction in [-1.5, 0.0, 1.5] {
        let chi = matrix(2, &[(0, 1, interaction)]);
        let network = network(
            &[0.018, 0.046],
            chi.clone(),
            &[0.001; 2],
            vec![face(0, 1)],
            vec![uniform_matrix(2, 1e-4)],
        );
        let initial = field(&network, &[&[0.2, 0.8], &[0.8, 0.2]]);
        let duration = 0.02;
        let candidate = initial
            .step_candidate(duration, MulticomponentMiscibilitySolve::default(), || {
                false
            })
            .unwrap();
        let exchange_jump = 2.0 * 4.0_f64.ln() - 1.2 * interaction;
        let expected = -duration * CONCENTRATION * 1e-4 * 0.25 * exchange_jump;
        near(candidate.step().transferred_mol[0][0], expected, 2e-16);
        near(candidate.step().transferred_mol[0][1], -expected, 2e-16);
        near(
            candidate.field().component_moles()[0],
            0.2 - expected,
            2e-15,
        );
        near(
            candidate.field().component_moles()[2],
            0.8 + expected,
            2e-15,
        );
        let actual_change =
            gibbs_energy(candidate.field(), 2, &chi) - gibbs_energy(&initial, 2, &chi);
        assert!(actual_change < 0.0);
        near(candidate.step().free_energy_change_j, actual_change, 3e-9);
    }
}

#[test]
fn ternary_friction_drives_one_species_uphill_while_total_free_energy_falls() {
    let chi = matrix(3, &[(0, 1, 0.2), (0, 2, -0.3), (1, 2, 0.1)]);
    let diffusivity = matrix(3, &[(0, 1, 1e-6), (0, 2, 1e-4), (1, 2, 1e-4)]);
    let network = network(
        &[0.002, 0.028, 0.044],
        chi.clone(),
        &[0.001; 2],
        vec![face(0, 1)],
        vec![diffusivity.clone()],
    );
    let [left, right] = [[0.3, 0.3, 0.4], [0.31, 0.1, 0.59]];
    let initial = field(&network, &[&left, &right]);
    let duration = 0.1;
    let candidate = initial
        .step_candidate(duration, MulticomponentMiscibilitySolve::default(), || {
            false
        })
        .unwrap();
    let expected = ternary_flux(left, right, &chi, &diffusivity, 1.0);
    assert!(right[0] > left[0] && expected[0] > 0.0);
    assert!(candidate.field().component_moles()[0] < initial.component_moles()[0]);
    for s in 0..3 {
        near(
            candidate.step().transferred_mol[0][s],
            duration * expected[s],
            3e-15,
        );
    }
    let jumps = chemical_jumps(left, right, &chi);
    let dissipated = -R * TEMPERATURE * (0..3).map(|s| expected[s] * jumps[s]).sum::<f64>();
    assert!(dissipated > 0.0);
    near(candidate.step().dissipation_rate_w, dissipated, 2e-10);
    let actual_change = gibbs_energy(candidate.field(), 3, &chi) - gibbs_energy(&initial, 3, &chi);
    assert!(actual_change < 0.0);
    near(candidate.step().free_energy_change_j, actual_change, 3e-9);
}

#[test]
fn unequal_volumes_and_masses_share_one_synchronous_conservative_face_update() {
    let masses = [0.002, 0.028, 0.044];
    let volumes = [0.0007, 0.0013, 0.0009];
    let chi = matrix(3, &[(0, 1, 0.2), (0, 2, -0.3), (1, 2, 0.1)]);
    let d = [
        matrix(3, &[(0, 1, 1e-6), (0, 2, 1e-4), (1, 2, 1e-4)]),
        matrix(3, &[(0, 1, 2e-4), (0, 2, 5e-5), (1, 2, 3e-4)]),
    ];
    let faces = vec![
        MulticomponentMiscibilityFace {
            area_m2: 0.7,
            distance_m: 0.4,
            ..face(0, 1)
        },
        MulticomponentMiscibilityFace {
            area_m2: 1.3,
            distance_m: 0.7,
            diffusivity_matrix: 1,
            ..face(1, 2)
        },
    ];
    let network = network(&masses, chi.clone(), &volumes, faces.clone(), d.to_vec());
    let fractions = [[0.3, 0.3, 0.4], [0.31, 0.1, 0.59], [0.2, 0.5, 0.3]];
    let initial = field(&network, &[&fractions[0], &fractions[1], &fractions[2]]);
    let dt = 0.01;
    let candidate = initial
        .step_candidate(dt, MulticomponentMiscibilitySolve::default(), || false)
        .unwrap();
    let mut expected = initial.component_moles().to_vec();
    for (f, face) in faces.iter().enumerate() {
        let [i, j] = face.cells;
        let rate = ternary_flux(fractions[i], fractions[j], &chi, &d[f], face.distance_m);
        for s in 0..3 {
            let transferred = dt * face.area_m2 * rate[s];
            near(candidate.step().transferred_mol[f][s], transferred, 3e-15);
            expected[3 * i + s] -= transferred;
            expected[3 * j + s] += transferred;
        }
    }
    for (actual, expected) in candidate.field().component_moles().iter().zip(expected) {
        near(*actual, expected, 2e-14);
    }
    for s in 0..3 {
        let total = |field: &MulticomponentMiscibilityField| {
            field
                .component_moles()
                .chunks_exact(3)
                .map(|n| n[s])
                .sum::<f64>()
        };
        near(total(candidate.field()), total(&initial), 3e-14);
    }
    for (i, &volume) in volumes.iter().enumerate() {
        near(
            candidate
                .field()
                .cell_component_moles(i)
                .unwrap()
                .iter()
                .sum(),
            CONCENTRATION * volume,
            3e-14,
        );
    }
    let mass = |field: &MulticomponentMiscibilityField, cell: usize| {
        field
            .cell_component_moles(cell)
            .unwrap()
            .iter()
            .zip(masses)
            .map(|(n, m)| n * m)
            .sum::<f64>()
    };
    near(
        (0..3).map(|i| mass(candidate.field(), i)).sum(),
        (0..3).map(|i| mass(&initial, i)).sum(),
        3e-15,
    );
    assert!((mass(candidate.field(), 0) - mass(&initial, 0)).abs() > 1e-7);
    assert!(candidate.step().free_energy_change_j < 0.0);
}

#[test]
fn authored_face_order_and_orientation_do_not_change_canonical_results() {
    let make = |faces| {
        network(
            &[0.002, 0.028, 0.044],
            vec![0.0; 9],
            &[0.001; 3],
            faces,
            vec![uniform_matrix(3, 1e-4)],
        )
    };
    let a = make(vec![face(0, 1), face(1, 2)]);
    let b = make(vec![face(2, 1), face(1, 0)]);
    let x: &[&[f64]] = &[&[0.2, 0.3, 0.5], &[0.4, 0.2, 0.4], &[0.3, 0.5, 0.2]];
    let first = field(&a, x)
        .step_candidate(0.02, MulticomponentMiscibilitySolve::default(), || false)
        .unwrap();
    let reordered = field(&b, x)
        .step_candidate(0.02, MulticomponentMiscibilitySolve::default(), || false)
        .unwrap();
    assert_eq!(first, reordered);
}

#[test]
fn positive_trace_inventory_moves_without_a_concentration_floor() {
    let network = network(
        &[0.002, 0.028],
        vec![0.0; 4],
        &[0.001; 2],
        vec![face(0, 1)],
        vec![uniform_matrix(2, 1e-4)],
    );
    let initial = field(&network, &[&[1e-120, 1.0], &[2e-120, 1.0]]);
    let dt = 0.1;
    let candidate = initial
        .step_candidate(dt, MulticomponentMiscibilitySolve::default(), || false)
        .unwrap();
    let transferred = -dt * CONCENTRATION * 1e-4 * 1.5e-120 * 2.0_f64.ln();
    near(
        candidate.step().transferred_mol[0][0],
        transferred,
        transferred.abs() * 2e-12,
    );
    near(
        candidate.step().transferred_mol[0][1],
        -transferred,
        transferred.abs() * 2e-12,
    );
    near(
        candidate.field().component_moles()[0],
        1e-120 - transferred,
        3e-134,
    );
    near(
        candidate.field().component_moles()[2],
        2e-120 + transferred,
        3e-134,
    );
    assert!(candidate.field().component_moles().iter().all(|v| *v > 0.0));
    assert_eq!(candidate.field().component_moles()[1], 1.0);
}

#[test]
fn admissions_reject_foreign_species_workspace_absence_and_non_global_convexity() {
    let make = |law, closure, budget| {
        MulticomponentMiscibilityNetwork::new(
            law,
            closure,
            Temperature::new(TEMPERATURE),
            AmountConcentration::new(CONCENTRATION),
            vec![MulticomponentMiscibilityCell { volume_m3: 0.001 }; 2],
            vec![face(0, 1)],
            vec![uniform_matrix(3, 1e-4)],
            budget,
            || false,
        )
    };
    let wrong = MaxwellStefanV1::new(vec![
        (SpeciesId::new("S000").unwrap(), MolarMass::new(0.002)),
        (SpeciesId::new("S001").unwrap(), MolarMass::new(0.028)),
        (SpeciesId::new("foreign").unwrap(), MolarMass::new(0.044)),
    ])
    .unwrap();
    assert!(make(model(3, vec![0.0; 9]), wrong, WORKSPACE).is_err());
    assert!(matches!(
        make(model(3, vec![0.0; 9]), transport(&[0.002, 0.028, 0.044]), 0),
        Err(MulticomponentMiscibilityError::Workspace { .. })
    ));
    // Every binary pair is convex and the UNIFORM ternary state is locally
    // stable (lambda_min=1-8.5/9>0), yet this model is not globally convex:
    // at x=[1/4,1/4,1/2], direction [1,1,-2] gives 16-17=-1.
    let attractive = model(3, matrix(3, &[(0, 1, -8.5)]));
    let unstable = attractive
        .evaluate(Temperature::new(TEMPERATURE), &[0.25, 0.25, 0.5], || false)
        .unwrap();
    near(
        unstable.reduced_curvature(&[1.0, 1.0, -2.0]).unwrap(),
        -1.0,
        2e-13,
    );
    assert!(matches!(
        make(attractive, transport(&[0.002, 0.028, 0.044]), WORKSPACE),
        Err(MulticomponentMiscibilityError::Convexity { .. })
    ));
    let admitted = make(
        model(3, vec![0.0; 9]),
        transport(&[0.002, 0.028, 0.044]),
        WORKSPACE,
    )
    .unwrap();
    assert!(
        admitted
            .bind(vec![vec![0.0, 0.4, 0.6], vec![0.2, 0.3, 0.5]], || false)
            .is_err()
    );
}

#[test]
fn large_steps_refuse_negative_inventory_and_positive_but_increasing_free_energy() {
    let network = network(
        &[0.018, 0.046],
        vec![0.0; 4],
        &[0.001; 2],
        vec![face(0, 1)],
        vec![uniform_matrix(2, 1e-4)],
    );
    let initial = field(&network, &[&[0.2, 0.8], &[0.8, 0.2]]);
    let before = initial.clone();
    let rate = CONCENTRATION * 1e-4 * 0.5 * 4.0_f64.ln();
    assert!(0.2 + 20.0 * rate > 1.0);
    assert!(matches!(
        initial.step_candidate(20.0, MulticomponentMiscibilitySolve::default(), || false),
        Err(MulticomponentMiscibilityError::Positivity { .. })
    ));
    // At 10 s every Euler inventory remains positive, but the profile has
    // overshot past its reflected initial state, increasing the convex G.
    let overshoot = 0.2 + 10.0 * rate;
    assert!(overshoot > 0.8 && overshoot < 1.0);
    assert!(matches!(
        initial.step_candidate(10.0, MulticomponentMiscibilitySolve::default(), || false),
        Err(MulticomponentMiscibilityError::FreeEnergy { .. })
    ));
    assert_eq!(initial, before);
}

#[test]
fn final_cancellation_is_atomic_and_replay_retains_all_original_receipts() {
    let network = network(
        &[0.002, 0.028, 0.044],
        vec![0.0; 9],
        &[0.001; 3],
        vec![face(0, 1), face(1, 2)],
        vec![uniform_matrix(3, 1e-4)],
    );
    let initial = field(
        &network,
        &[&[0.2, 0.3, 0.5], &[0.4, 0.2, 0.4], &[0.3, 0.5, 0.2]],
    );
    let saved = initial.clone();
    let controls = MulticomponentMiscibilitySolve::default();
    let mut polls = 0usize;
    let expected = initial
        .step_candidate(0.02, controls, || {
            polls += 1;
            false
        })
        .unwrap();
    assert!(polls > 3);
    let mut seen = 0usize;
    let result = initial.step_candidate(0.02, controls, || {
        seen += 1;
        seen == polls
    });
    assert_eq!(
        result.unwrap_err(),
        MulticomponentMiscibilityError::Cancelled
    );
    assert_eq!(initial, saved);
    assert_eq!(
        initial.step_candidate(0.02, controls, || false).unwrap(),
        expected
    );
}
