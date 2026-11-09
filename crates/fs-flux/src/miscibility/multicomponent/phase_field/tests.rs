//! Independent graph-dispersion, full-friction and inventory phase-field oracles.

use super::*;
use fs_qty::{AmountConcentration, MolarMass, Temperature};
use fs_thermochem::SpeciesId;
use fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R;
use fs_thermochem::solution::multicomponent::MulticomponentRegularSolutionV1;
use fs_thermochem::transport::maxwell_stefan::MaxwellStefanV1;

pub(super) const WORKSPACE: usize = 16 * 1024 * 1024;
pub(super) const CONCENTRATION: f64 = 1000.0;
pub(super) const TEMPERATURE: f64 = 300.0;

pub(super) fn uniform_matrix(n: usize, value: f64) -> Vec<f64> {
    let mut result = vec![value; n * n];
    for s in 0..n {
        result[s * n + s] = 0.0;
    }
    result
}

pub(super) fn matrix(n: usize, pairs: &[(usize, usize, f64)]) -> Vec<f64> {
    let mut result = vec![0.0; n * n];
    for &(i, j, value) in pairs {
        result[i * n + j] = value;
        result[j * n + i] = value;
    }
    result
}

fn species(n: usize) -> Vec<SpeciesId> {
    (0..n)
        .map(|s| SpeciesId::new(format!("S{s:03}")).unwrap())
        .collect()
}

#[allow(clippy::too_many_arguments)]
fn construct(
    masses: &[f64],
    chi: Vec<f64>,
    volumes: &[f64],
    faces: Vec<MulticomponentMiscibilityFace>,
    diffusivity_matrices: Vec<Vec<f64>>,
    kappa: f64,
    maximum_workspace_bytes: usize,
    cancelled: bool,
) -> Result<MulticomponentPhaseFieldNetwork> {
    let n = masses.len();
    let model = MulticomponentRegularSolutionV1::new(
        species(n),
        chi,
        vec![0.0; n * n],
        Temperature::new(250.0),
        Temperature::new(700.0),
    )
    .unwrap();
    let transport = MaxwellStefanV1::new(
        species(n)
            .into_iter()
            .zip(masses)
            .map(|(id, &mass)| (id, MolarMass::new(mass)))
            .collect(),
    )
    .unwrap();
    MulticomponentPhaseFieldNetwork::new(
        model,
        transport,
        Temperature::new(TEMPERATURE),
        AmountConcentration::new(CONCENTRATION),
        volumes
            .iter()
            .map(|&volume_m3| MulticomponentMiscibilityCell { volume_m3 })
            .collect(),
        faces,
        diffusivity_matrices,
        kappa,
        maximum_workspace_bytes,
        || cancelled,
    )
}

pub(super) fn network(
    masses: &[f64],
    chi: Vec<f64>,
    volumes: &[f64],
    faces: Vec<MulticomponentMiscibilityFace>,
    diffusivity_matrices: Vec<Vec<f64>>,
    kappa: f64,
) -> MulticomponentPhaseFieldNetwork {
    construct(
        masses,
        chi,
        volumes,
        faces,
        diffusivity_matrices,
        kappa,
        WORKSPACE,
        false,
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

pub(super) fn field(
    network: &MulticomponentPhaseFieldNetwork,
    fractions: Vec<Vec<f64>>,
) -> MulticomponentPhaseField {
    network
        .bind(
            fractions
                .into_iter()
                .zip(network.capacities_mol())
                .map(|(x, capacity)| x.into_iter().map(|x| x * capacity).collect())
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

// These energy oracles evaluate the declared original potential directly.
// They do not call the production energy or chemical-potential evaluators.
pub(super) fn bulk_energy(field: &MulticomponentPhaseField, chi: &[f64]) -> f64 {
    let n = field.network().species_count();
    field
        .component_moles()
        .chunks_exact(n)
        .map(|amounts| {
            let total: f64 = amounts.iter().sum();
            let mut reduced: f64 = amounts.iter().map(|a| a * (a / total).ln()).sum();
            for s in 0..n {
                for t in s + 1..n {
                    reduced += chi[s * n + t] * amounts[s] * amounts[t] / total;
                }
            }
            R * TEMPERATURE * reduced
        })
        .sum()
}

pub(super) fn interface_energy(field: &MulticomponentPhaseField, kappa: f64) -> f64 {
    field
        .network()
        .faces()
        .iter()
        .map(|face| {
            let left = field.cell_mole_fractions(face.cells[0]).unwrap();
            let right = field.cell_mole_fractions(face.cells[1]).unwrap();
            kappa * face.area_m2 / (4.0 * face.distance_m)
                * left
                    .iter()
                    .zip(right)
                    .map(|(a, b)| (b - a) * (b - a))
                    .sum::<f64>()
        })
        .sum()
}

pub(super) fn assert_inventories(
    initial: &MulticomponentPhaseField,
    final_field: &MulticomponentPhaseField,
) {
    let n = initial.network().species_count();
    for s in 0..n {
        let before: f64 = initial.component_moles().iter().skip(s).step_by(n).sum();
        let after: f64 = final_field
            .component_moles()
            .iter()
            .skip(s)
            .step_by(n)
            .sum();
        near(after, before, 4e-12);
    }
    for (cell, &capacity) in initial.network().capacities_mol().iter().enumerate() {
        near(
            final_field.cell_component_moles(cell).unwrap().iter().sum(),
            capacity,
            4e-12,
        );
    }
    assert!(
        final_field
            .component_moles()
            .iter()
            .all(|amount| amount.is_finite() && *amount > 0.0)
    );
}

pub(super) fn ring(mode: usize, amplitude: f64) -> (MulticomponentPhaseField, Vec<f64>, f64) {
    let cells = 16;
    let chi = uniform_matrix(3, 4.0);
    let kappa = CONCENTRATION * R * TEMPERATURE;
    let mut faces = (0..cells - 1).map(|i| face(i, i + 1)).collect::<Vec<_>>();
    faces.push(face(0, cells - 1));
    let network = network(
        &[0.018, 0.046, 0.092],
        chi.clone(),
        &vec![1.0; cells],
        faces,
        vec![uniform_matrix(3, 0.2)],
        kappa,
    );
    let fractions = (0..cells)
        .map(|cell| {
            let phase = 2.0 * std::f64::consts::PI * mode as f64 * cell as f64 / cells as f64;
            let perturbation = amplitude * phase.cos();
            vec![1.0 / 3.0 + perturbation, 1.0 / 3.0 - perturbation, 1.0 / 3.0]
        })
        .collect();
    (field(&network, fractions), chi, kappa)
}

#[test]
fn graph_spinodal_dispersion_grows_long_waves_and_damps_short_waves() {
    let amplitude = 1e-4;
    let duration = 0.01;
    for mode in [2, 3, 5] {
        let (initial, chi, kappa) = ring(mode, amplitude);
        let candidate = initial
            .step_candidate(duration, MulticomponentMiscibilitySolve::default(), || false)
            .unwrap();
        let cells = initial.network().cells().len();
        let cosines = (0..cells)
            .map(|cell| {
                (2.0 * std::f64::consts::PI * mode as f64 * cell as f64 / cells as f64).cos()
            })
            .collect::<Vec<_>>();
        // Exact graph eigenvalue, not the continuum k^2 approximation.
        // At equimolar x and equal pair diffusivities, tangent mobility is D/n.
        let q = 4.0 * (std::f64::consts::PI * mode as f64 / cells as f64).sin().powi(2);
        let expected_growth = -0.2 * q * (-1.0 + 0.5 * q) / 3.0;
        // Project paired transfers directly, avoiding cancellation in almost
        // unchanged mole fractions and large absolute Gibbs energies.
        let observed_growth: f64 = initial
            .network()
            .faces()
            .iter()
            .zip(&candidate.step().transferred_mol)
            .map(|(face, transfer)| {
                (cosines[face.cells[1]] - cosines[face.cells[0]])
                    * (transfer[0] - transfer[1])
                    / (cells as f64 * CONCENTRATION * duration * amplitude)
            })
            .sum();
        near(observed_growth, expected_growth, 5e-8);
        assert_eq!(observed_growth > 0.0, mode < 4);
        let bulk_change = bulk_energy(candidate.field(), &chi) - bulk_energy(&initial, &chi);
        let interface_change =
            interface_energy(candidate.field(), kappa) - interface_energy(&initial, kappa);
        near(candidate.step().bulk_free_energy_change_j, bulk_change, 2e-7);
        near(candidate.step().interfacial_energy_change_j, interface_change, 2e-10);
        near(
            candidate.step().free_energy_change_j,
            bulk_change + interface_change,
            2e-7,
        );
        assert!(candidate.step().free_energy_change_j < 0.0);
        if mode == 5 {
            // Smoothing a spinodal composition raises homogeneous Gibbs
            // energy; the larger interface-energy reduction makes it admissible.
            assert!(candidate.step().bulk_free_energy_change_j > 0.0);
            assert!(candidate.step().interfacial_energy_change_j < 0.0);
        } else {
            assert!(candidate.step().bulk_free_energy_change_j < 0.0);
            assert!(candidate.step().interfacial_energy_change_j > 0.0);
        }
        assert_inventories(&initial, candidate.field());
    }
}

fn ternary_fixture() -> (MulticomponentPhaseField, Vec<f64>, Vec<f64>, f64) {
    let chi = matrix(3, &[(0, 1, 4.2), (0, 2, 1.7), (1, 2, 3.4)]);
    let diffusivities = matrix(3, &[(0, 1, 0.04), (0, 2, 0.013), (1, 2, 0.07)]);
    let kappa = 0.12 * CONCENTRATION * R * TEMPERATURE;
    let mut link = face(0, 1);
    link.area_m2 = 1.3;
    link.distance_m = 0.8;
    let network = network(
        &[0.018, 0.046, 0.092],
        chi.clone(),
        &[0.7, 1.9],
        vec![link],
        vec![diffusivities.clone()],
        kappa,
    );
    (
        field(&network, vec![vec![0.44, 0.27, 0.29], vec![0.22, 0.51, 0.27]]),
        chi,
        diffusivities,
        kappa,
    )
}

// Independent elimination of J_C=-J_A-J_B in the original friction equations.
// This never calls a production Maxwell-Stefan or gradient-force evaluator.
fn ternary_oracle(
    initial: &MulticomponentPhaseField,
    chi: &[f64],
    diffusivities: &[f64],
    kappa: f64,
) -> ([f64; 3], [f64; 3], [f64; 3]) {
    let face = initial.network().faces()[0];
    let left = initial.cell_mole_fractions(face.cells[0]).unwrap();
    let right = initial.cell_mole_fractions(face.cells[1]).unwrap();
    let capacities = initial.network().capacities_mol();
    let gradient = kappa * face.area_m2 / (2.0 * face.distance_m * R * TEMPERATURE)
        * (1.0 / capacities[face.cells[0]] + 1.0 / capacities[face.cells[1]]);
    let mean: [f64; 3] = core::array::from_fn(|s| 0.5 * (left[s] + right[s]));
    let jumps: [f64; 3] = core::array::from_fn(|s| {
        (right[s] / left[s]).ln()
            + (0..3).map(|t| chi[s * 3 + t] * (right[t] - left[t])).sum::<f64>()
            + gradient * (right[s] - left[s])
    });
    let common: f64 = (0..3).map(|s| mean[s] * jumps[s]).sum();
    let driving: [f64; 3] =
        core::array::from_fn(|s| mean[s] * (jumps[s] - common) / face.distance_m);
    let [ab, ac, bc] = [diffusivities[1], diffusivities[2], diffusivities[5]];
    let a = mean[1] / ab + (mean[0] + mean[2]) / ac;
    let b = mean[0] * (1.0 / ac - 1.0 / ab);
    let c = mean[1] * (1.0 / bc - 1.0 / ab);
    let d = mean[0] / ab + (mean[1] + mean[2]) / bc;
    let determinant = a * d - b * c;
    let ja = -CONCENTRATION * (d * driving[0] - b * driving[1]) / determinant;
    let jb = -CONCENTRATION * (a * driving[1] - c * driving[0]) / determinant;
    ([ja, jb, -ja - jb], jumps, mean)
}

#[test]
fn unequal_pair_friction_and_unequal_capacities_match_independent_gradient_oracle() {
    let (initial, chi, diffusivities, kappa) = ternary_fixture();
    let duration = 0.005;
    let candidate = initial
        .step_candidate(duration, MulticomponentMiscibilitySolve::default(), || false)
        .unwrap();
    let (fluxes, jumps, mean) = ternary_oracle(&initial, &chi, &diffusivities, kappa);
    let face = initial.network().faces()[0];
    for (s, &flux) in fluxes.iter().enumerate() {
        near(
            candidate.step().transferred_mol[0][s] / (duration * face.area_m2),
            flux,
            1e-9,
        );
    }
    let force_power = -R * TEMPERATURE * face.area_m2
        * (0..3).map(|s| fluxes[s] * jumps[s]).sum::<f64>();
    let mut friction_power = 0.0;
    for s in 0..3 {
        for t in s + 1..3 {
            let relative_flux = mean[t] * fluxes[s] - mean[s] * fluxes[t];
            friction_power += relative_flux * relative_flux
                / (diffusivities[s * 3 + t] * mean[s] * mean[t]);
        }
    }
    friction_power *= R * TEMPERATURE * face.area_m2 * face.distance_m / CONCENTRATION;
    near(force_power, friction_power, 1e-8);
    near(candidate.step().dissipation_rate_w, force_power, 1e-7);
    let bulk_change = bulk_energy(candidate.field(), &chi) - bulk_energy(&initial, &chi);
    let interface_change =
        interface_energy(candidate.field(), kappa) - interface_energy(&initial, kappa);
    near(candidate.step().bulk_free_energy_change_j, bulk_change, 2e-7);
    near(candidate.step().interfacial_energy_change_j, interface_change, 1e-8);
    assert!(candidate.step().free_energy_change_j < 0.0);
    assert_inventories(&initial, candidate.field());
}

#[test]
fn species_permutation_preserves_full_gradient_coupling_and_energy() {
    let (initial, chi, diffusivities, kappa) = ternary_fixture();
    let permutation = [2, 0, 1];
    let permute = |matrix: &[f64]| {
        permutation
            .iter()
            .flat_map(|&s| permutation.iter().map(move |&t| matrix[s * 3 + t]))
            .collect::<Vec<_>>()
    };
    let permuted_network = network(
        &[0.092, 0.018, 0.046],
        permute(&chi),
        &[0.7, 1.9],
        initial.network().faces().to_vec(),
        vec![permute(&diffusivities)],
        kappa,
    );
    let permuted = field(
        &permuted_network,
        (0..2)
            .map(|cell| {
                let original = initial.cell_mole_fractions(cell).unwrap();
                permutation.iter().map(|&s| original[s]).collect()
            })
            .collect(),
    );
    let original = initial
        .step_candidate(0.005, MulticomponentMiscibilitySolve::default(), || false)
        .unwrap();
    let changed = permuted
        .step_candidate(0.005, MulticomponentMiscibilitySolve::default(), || false)
        .unwrap();
    for (new_s, &old_s) in permutation.iter().enumerate() {
        near(
            changed.step().transferred_mol[0][new_s],
            original.step().transferred_mol[0][old_s],
            1e-11,
        );
        for cell in 0..2 {
            near(
                changed.field().cell_mole_fractions(cell).unwrap()[new_s],
                original.field().cell_mole_fractions(cell).unwrap()[old_s],
                2e-14,
            );
        }
    }
    near(changed.step().free_energy_change_j, original.step().free_energy_change_j, 1e-7);
    near(changed.step().dissipation_rate_w, original.step().dissipation_rate_w, 1e-7);
}

#[test]
fn binary_limit_preserves_the_original_scalar_kappa_without_a_factor_of_two() {
    let interaction = 3.2;
    let kappa = 0.3 * CONCENTRATION * R * TEMPERATURE;
    let diffusivity = 0.023;
    let mut link = face(0, 1);
    link.area_m2 = 0.7;
    link.distance_m = 1.4;
    let network = network(
        &[0.018, 0.092],
        uniform_matrix(2, interaction),
        &[0.9, 1.8],
        vec![link],
        vec![uniform_matrix(2, diffusivity)],
        kappa,
    );
    let initial = field(&network, vec![vec![0.38, 0.62], vec![0.61, 0.39]]);
    let duration = 0.01;
    let candidate = initial
        .step_candidate(duration, MulticomponentMiscibilitySolve::default(), || false)
        .unwrap();
    let [left, right] = [0.38_f64, 0.61_f64];
    let mean = 0.5 * (left + right);
    let reduced_exchange = (right / (1.0 - right)).ln() - (left / (1.0 - left)).ln()
        - 2.0 * interaction * (right - left)
        + kappa * link.area_m2 / (link.distance_m * R * TEMPERATURE)
            * (1.0 / network.capacities_mol()[0] + 1.0 / network.capacities_mol()[1])
            * (right - left);
    let flux = -CONCENTRATION * diffusivity * mean * (1.0 - mean) * reduced_exchange
        / link.distance_m;
    near(candidate.step().transferred_mol[0][0], flux * link.area_m2 * duration, 1e-11);
    near(candidate.step().transferred_mol[0][1], -flux * link.area_m2 * duration, 1e-11);
    near(
        interface_energy(&initial, kappa),
        kappa * link.area_m2 * (right - left).powi(2) / (2.0 * link.distance_m),
        1e-10,
    );
    assert_inventories(&initial, candidate.field());
}

#[test]
fn every_original_component_is_conserved_through_the_supported_species_bound() {
    for n in [4, 17, 128] {
        let network = network(
            &vec![0.04; n],
            uniform_matrix(n, n as f64 + 1.0),
            &[1.0, 1.0],
            vec![face(0, 1)],
            vec![uniform_matrix(n, 0.02)],
            0.25 * CONCENTRATION * R * TEMPERATURE,
        );
        let left = (0..n)
            .map(|s| {
                (1.0 + 0.01 * (2.0 * std::f64::consts::PI * s as f64 / n as f64).cos())
                    / n as f64
            })
            .collect();
        let right = (0..n)
            .map(|s| {
                (1.0 - 0.01 * (2.0 * std::f64::consts::PI * s as f64 / n as f64).cos())
                    / n as f64
            })
            .collect();
        let initial = field(&network, vec![left, right]);
        let candidate = initial
            .step_candidate(0.02, MulticomponentMiscibilitySolve::default(), || false)
            .unwrap();
        assert_eq!(candidate.step().transferred_mol[0].len(), n);
        assert_eq!(candidate.step().component_residual_mol.len(), n);
        assert!(
            candidate
                .step()
                .component_residual_mol
                .iter()
                .all(|residual| residual.abs() <= 1e-12)
        );
        near(candidate.step().transferred_mol[0].iter().sum(), 0.0, 1e-12);
        assert!(candidate.step().free_energy_change_j <= 1e-9);
        assert_inventories(&initial, candidate.field());
    }
}

#[test]
fn nonconvex_admission_is_confined_to_positive_gradient_owner_and_bounded_workspace() {
    let build = |kappa, workspace, cancelled| {
        construct(
            &[0.018, 0.046, 0.092],
            uniform_matrix(3, 4.0),
            &[1.0, 1.0],
            vec![face(0, 1)],
            vec![uniform_matrix(3, 0.02)],
            kappa,
            workspace,
            cancelled,
        )
    };
    for kappa in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(build(kappa, WORKSPACE, false).is_err());
    }
    let kappa = CONCENTRATION * R * TEMPERATURE;
    let network = build(kappa, WORKSPACE, false).unwrap();
    assert!(matches!(
        build(kappa, network.workspace_bytes() - 1, false),
        Err(MulticomponentMiscibilityError::Workspace { .. })
    ));
    assert!(build(kappa, network.workspace_bytes(), false).is_ok());
    assert!(matches!(
        build(kappa, WORKSPACE, true),
        Err(MulticomponentMiscibilityError::Cancelled)
    ));
    // The ordinary zero-gradient transport API must keep its global-convexity
    // requirement; a phase-field network must never leak that weaker contract.
    assert!(matches!(
        MulticomponentMiscibilityNetwork::new(
            network.model().clone(),
            network.transport().clone(),
            network.temperature(),
            network.concentration(),
            network.cells().to_vec(),
            network.faces().to_vec(),
            vec![uniform_matrix(3, 0.02)],
            WORKSPACE,
            || false,
        ),
        Err(MulticomponentMiscibilityError::Convexity { .. })
    ));
}

#[test]
fn oversized_or_cancelled_phase_field_steps_leave_every_inventory_unchanged() {
    let (initial, _, _, _) = ternary_fixture();
    let original = initial.clone();
    assert!(matches!(
        initial.step_candidate(1e12, MulticomponentMiscibilitySolve::default(), || false),
        Err(MulticomponentMiscibilityError::Positivity { .. })
    ));
    assert_eq!(initial, original);
    let mut polls = 0;
    assert!(matches!(
        initial.step_candidate(0.005, MulticomponentMiscibilitySolve::default(), || {
            polls += 1;
            polls > 8
        }),
        Err(MulticomponentMiscibilityError::Cancelled)
    ));
    assert!(polls > 8);
    assert_eq!(initial, original);
}
