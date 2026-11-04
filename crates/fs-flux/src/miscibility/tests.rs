//! Independent energy derivatives, spinodal rates and conservative time evolution.

use super::*;
use fs_qty::{AmountConcentration, Temperature};
use fs_thermochem::solution::BinaryRegularSolutionV1;
use fs_thermochem::{SpeciesId, UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R};

fn model(chi: f64, scaled_kappa: f64) -> BinarySquareGradientV1 {
    BinarySquareGradientV1::new(
        BinaryRegularSolutionV1::new(
            [SpeciesId::new("A").unwrap(), SpeciesId::new("B").unwrap()],
            chi,
            0.0,
            Temperature::new(250.0),
            Temperature::new(600.0),
        )
        .unwrap(),
        Temperature::new(300.0),
        AmountConcentration::new(1.0),
        scaled_kappa * R * 300.0,
    )
    .unwrap()
}

fn face(i: usize, j: usize, distance: f64, diffusivity: f64) -> BinaryMiscibilityFace {
    BinaryMiscibilityFace {
        cells: [i, j],
        area_m2: 1.0,
        distance_m: distance,
        diffusivity_m2_s: diffusivity,
    }
}

fn network(
    chi: f64,
    kappa: f64,
    volumes: &[f64],
    faces: Vec<BinaryMiscibilityFace>,
) -> BinaryMiscibilityNetwork {
    BinaryMiscibilityNetwork::new(
        model(chi, kappa),
        volumes
            .iter()
            .map(|&volume_m3| BinaryMiscibilityCell { volume_m3 })
            .collect(),
        faces,
        16 * 1024 * 1024,
        || false,
    )
    .unwrap()
}

fn state(network: &BinaryMiscibilityNetwork, x: &[f64]) -> BinaryMiscibilityState {
    network
        .state_from_mole_fractions(x.iter().map(|&x| [x, 1.0 - x]).collect(), || false)
        .unwrap()
}

fn near(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a:.17e} != {b:.17e}; tolerance {tolerance:e}"
    );
}

fn strict() -> BinaryMiscibilitySolve {
    BinaryMiscibilitySolve {
        absolute_transfer_tolerance_mol: 1e-18,
        relative_transfer_tolerance: 1e-11,
        amount_balance_tolerance_mol: 1e-13,
        free_energy_tolerance_j: 1e-9,
        ..BinaryMiscibilitySolve::default()
    }
}

#[test]
fn chemical_potential_is_the_full_weighted_graph_energy_derivative() {
    let net = network(
        3.0,
        0.2,
        &[0.7, 1.3, 0.5],
        vec![face(0, 1, 0.4, 1.0), face(1, 2, 0.7, 1.0)],
    );
    let initial = state(&net, &[0.15, 0.4, 0.8]);
    let potential = net
        .exchange_potentials_j_per_mol(&initial, || false)
        .unwrap();
    for [i, j] in [[0, 1], [1, 2]] {
        let h = 1e-7;
        let shifted = |q: f64| {
            let mut amounts = initial.component_moles().to_vec();
            amounts[i][0] -= q;
            amounts[i][1] += q;
            amounts[j][0] += q;
            amounts[j][1] -= q;
            net.state_from_component_moles(amounts, || false).unwrap()
        };
        let derivative = (net.free_energy_j(&shifted(h), || false).unwrap()
            - net.free_energy_j(&shifted(-h), || false).unwrap())
            / (2.0 * h);
        near(
            derivative,
            potential[j].unwrap() - potential[i].unwrap(),
            1e-5,
        );
    }
}

#[test]
fn pure_contact_mixes_implicitly_without_clipping_or_fictitious_trace() {
    let net = network(0.0, 0.0, &[1.0, 1.0], vec![face(0, 1, 1.0, 1.0)]);
    let initial = state(&net, &[1.0, 0.0]);
    assert_eq!(
        net.exchange_potentials_j_per_mol(&initial, || false)
            .unwrap(),
        [None, None]
    );
    let controls = BinaryMiscibilitySolve {
        relative_transfer_tolerance: 1e-10,
        ..BinaryMiscibilitySolve::default()
    };
    let candidate = net
        .step_candidate(&initial, 1000.0, controls, || false)
        .unwrap();
    for fraction in candidate.state().mole_fractions() {
        assert!(fraction.iter().all(|&x| x > 0.0 && x < 1.0));
        near(fraction[0], 0.5, 1e-3);
    }
    let step = candidate.step();
    assert!(step.transferred_a_mol[0] > 0.0);
    assert!(step.free_energy_change_j < -3000.0);
    assert!(step.convex_split_dissipation_j > 0.0);
    assert!(
        step.free_energy_change_j + step.convex_split_dissipation_j
            <= controls.free_energy_tolerance_j
    );
    for species in 0..2 {
        near(
            candidate
                .state()
                .component_moles()
                .iter()
                .map(|n| n[species])
                .sum(),
            1.0,
            1e-14,
        );
        near(step.component_residual_mol[species], 0.0, 1e-14);
    }
    assert_eq!(initial.mole_fractions(), &[[1.0, 0.0], [0.0, 1.0]]);
}

#[test]
fn a_discrete_interface_balances_bulk_and_all_gradient_chemical_potentials() {
    let x = 0.2_f64;
    let chi = 3.0;
    // For two equal cells, stationary mu requires
    // kappa/(N RT) = [ln(x/(1-x))+chi(1-2x)]/(1-2x).
    let kappa = ((x / (1.0 - x)).ln() + chi * (1.0 - 2.0 * x)) / (1.0 - 2.0 * x);
    let net = network(chi, kappa, &[1.0, 1.0], vec![face(0, 1, 1.0, 1.0)]);
    let initial = state(&net, &[x, 1.0 - x]);
    let potentials = net
        .exchange_potentials_j_per_mol(&initial, || false)
        .unwrap();
    near(potentials[0].unwrap(), 0.0, 2e-12);
    near(potentials[1].unwrap(), 0.0, 2e-12);
    let candidate = net
        .step_candidate(&initial, 1.0, BinaryMiscibilitySolve::default(), || false)
        .unwrap();
    assert_eq!(
        candidate.state().component_moles(),
        initial.component_moles()
    );
    assert_eq!(candidate.step().transferred_a_mol, [0.0]);
    assert!(candidate.step().maximum_root_residual_mol <= 1e-14);
}

fn periodic(n: usize, chi: f64, kappa: f64, diffusivity: f64) -> BinaryMiscibilityNetwork {
    let h = 1.0 / n as f64;
    network(
        chi,
        kappa,
        &vec![h; n],
        (0..n)
            .map(|i| face(i, (i + 1) % n, h, diffusivity))
            .collect(),
    )
}

fn cosine_state(
    net: &BinaryMiscibilityNetwork,
    mode: usize,
    amplitude: f64,
) -> BinaryMiscibilityState {
    let n = net.cells().len();
    state(
        net,
        &(0..n)
            .map(|i| {
                0.5 + amplitude * (core::f64::consts::TAU * mode as f64 * i as f64 / n as f64).cos()
            })
            .collect::<Vec<_>>(),
    )
}

fn projection(state: &BinaryMiscibilityState, mode: usize) -> f64 {
    let n = state.mole_fractions().len();
    2.0 / n as f64
        * state
            .mole_fractions()
            .iter()
            .enumerate()
            .map(|(i, x)| {
                (x[0] - 0.5) * (core::f64::consts::TAU * mode as f64 * i as f64 / n as f64).cos()
            })
            .sum::<f64>()
}

#[test]
fn spinodal_long_waves_grow_while_square_gradient_damps_short_waves() {
    let n = 32;
    let chi = 3.0;
    let kappa = 0.002;
    let net = periodic(n, chi, kappa, 1.0);
    let h = 1.0 / n as f64;
    let dt = 2e-7;
    for mode in [1, 12] {
        let initial = cosine_state(&net, mode, 1e-4);
        let candidate = net
            .step_candidate(&initial, dt, strict(), || false)
            .unwrap();
        assert!(candidate.step().root_iterations < 8 * n);
        let a0 = projection(&initial, mode);
        let a1 = projection(candidate.state(), mode);
        // Independent linearized discrete CH eigenvalue at x=1/2, c=1.
        let lambda = 4.0
            * (core::f64::consts::PI * mode as f64 / n as f64)
                .sin()
                .powi(2)
            / (h * h);
        let expected_rate = -0.25 * lambda * (4.0 - 2.0 * chi + kappa * lambda);
        let measured_rate = (a1 / a0 - 1.0) / dt;
        near(measured_rate, expected_rate, 0.01 * expected_rate.abs());
        assert_eq!(a1 > a0, mode == 1);
        assert!(candidate.step().free_energy_change_j <= strict().free_energy_tolerance_j);
    }
}

#[test]
fn whole_graph_gradient_cache_agrees_with_independent_pair_restarts() {
    let faces = vec![face(0, 1, 0.8, 0.3), face(1, 2, 1.2, 0.8)];
    let net = network(3.0, 0.1, &[0.7, 1.0, 1.3], faces.clone());
    let initial = state(&net, &[0.1, 0.55, 0.9]);
    let expected = net
        .step_candidate(&initial, 0.2, strict(), || false)
        .unwrap();
    let mut first_faces = faces.clone();
    first_faces[1].diffusivity_m2_s = 0.0;
    let first = network(3.0, 0.1, &[0.7, 1.0, 1.3], first_faces);
    let first_state = first
        .state_from_component_moles(initial.component_moles().to_vec(), || false)
        .unwrap();
    let middle = first
        .step_candidate(&first_state, 0.2, strict(), || false)
        .unwrap();
    let mut last_faces = faces;
    last_faces[0].diffusivity_m2_s = 0.0;
    let last = network(3.0, 0.1, &[0.7, 1.0, 1.3], last_faces);
    let middle = last
        .state_from_component_moles(middle.state().component_moles().to_vec(), || false)
        .unwrap();
    let result = last
        .step_candidate(&middle, 0.2, strict(), || false)
        .unwrap();
    for (a, b) in expected
        .state()
        .component_moles()
        .iter()
        .zip(result.state().component_moles())
    {
        for species in 0..2 {
            near(a[species], b[species], 2e-14);
        }
    }
}

#[test]
fn canonical_faces_component_exchange_and_closed_pure_regions_are_preserved() {
    let faces = vec![
        face(0, 1, 0.6, 0.7),
        face(1, 2, 0.9, 0.3),
        face(2, 3, 1.0, 0.0),
    ];
    let a = network(3.0, 0.03, &[0.7, 1.3, 0.8, 1.0], faces.clone());
    let b = network(
        3.0,
        0.03,
        &[0.7, 1.3, 0.8, 1.0],
        faces
            .iter()
            .rev()
            .map(|f| BinaryMiscibilityFace {
                cells: [f.cells[1], f.cells[0]],
                ..*f
            })
            .collect(),
    );
    let initial = state(&a, &[0.1, 0.7, 0.3, 1.0]);
    let first = a.step_candidate(&initial, 0.1, strict(), || false).unwrap();
    let reordered = b.step_candidate(&initial, 0.1, strict(), || false).unwrap();
    assert_eq!(first, reordered);
    assert_eq!(
        first.state().component_moles()[3],
        initial.component_moles()[3]
    );
    assert_eq!(first.step().transferred_a_mol[2], 0.0);
    let reflected = a
        .state_from_component_moles(
            initial
                .component_moles()
                .iter()
                .map(|n| [n[1], n[0]])
                .collect(),
            || false,
        )
        .unwrap();
    let reflected = a
        .step_candidate(&reflected, 0.1, strict(), || false)
        .unwrap();
    for (x, y) in first
        .state()
        .component_moles()
        .iter()
        .zip(reflected.state().component_moles())
    {
        near(x[0], y[1], 1e-12);
        near(x[1], y[0], 1e-12);
    }
    near(
        first.step().free_energy_change_j,
        reflected.step().free_energy_change_j,
        1e-9,
    );
}

#[test]
fn ideal_discrete_diffusion_has_first_order_temporal_convergence() {
    let n = 8;
    let diffusivity = 0.01;
    let net = periodic(n, 0.0, 0.0, diffusivity);
    let initial = cosine_state(&net, 1, 1e-3);
    let duration = 0.2;
    let lambda = 4.0 * (core::f64::consts::PI / n as f64).sin().powi(2) * (n * n) as f64;
    let exact = 1e-3 * (-diffusivity * lambda * duration).exp();
    let mut errors = Vec::new();
    for steps in [10, 20, 40] {
        let mut state = initial.clone();
        for _ in 0..steps {
            state = net
                .step_candidate(&state, duration / steps as f64, strict(), || false)
                .unwrap()
                .into_parts()
                .0;
        }
        errors.push((projection(&state, 1) - exact).abs());
    }
    for pair in errors.windows(2) {
        assert!(
            pair[0] / pair[1] > 1.7 && pair[0] / pair[1] < 2.4,
            "first-order errors: {errors:?}"
        );
    }
}

#[test]
fn repeated_spinodal_exchanges_preserve_original_capacity_and_species() {
    let net = network(3.0, 0.1, &[1.0, 1.0], vec![face(0, 1, 1.0, 1.0)]);
    let initial = state(&net, &[0.6, 0.4]);
    let mut current = initial.clone();
    // Small conservative exchanges must not accumulate independent inventory
    // rounding into a false capacity change over a physical evolution window.
    // This uses the original constitutive composition admission on every trial.
    let controls = BinaryMiscibilitySolve {
        amount_balance_tolerance_mol: 2.0 * f64::EPSILON,
        ..strict()
    };
    for _ in 0..4096 {
        let candidate = net
            .step_candidate(&current, 2.0 / 4096.0, controls, || false)
            .unwrap();
        let report = candidate.step();
        assert!(report.transferred_a_mol[0] < 0.0);
        assert!(report.free_energy_change_j < 0.0);
        assert!(report.maximum_capacity_residual_mol <= 2.0 * f64::EPSILON);
        current = candidate.into_parts().0;
        for pair in current.component_moles() {
            near(pair[0] + pair[1], 1.0, 2.0 * f64::EPSILON);
        }
        for species in 0..2 {
            near(
                current.component_moles()[0][species] + current.component_moles()[1][species],
                1.0,
                2.0 * f64::EPSILON,
            );
        }
    }
    assert!(current.mole_fractions()[0][0] > 0.7);
    assert_eq!(current.elapsed_s(), 2.0);
    assert_eq!(initial.component_moles(), &[[0.6, 0.4], [0.4, 0.6]]);
}

#[test]
fn cancellation_bounds_source_identity_and_root_budget_are_atomic() {
    let net = network(
        3.0,
        0.1,
        &[1.0, 1.0, 1.0],
        vec![face(0, 1, 1.0, 1.0), face(1, 2, 1.0, 1.0)],
    );
    let initial = state(&net, &[1.0, 0.0, 0.6]);
    let before = initial.clone();
    let mut calls = 0;
    net.step_candidate(&initial, 0.1, strict(), || {
        calls += 1;
        false
    })
    .unwrap();
    for stop in [1, calls / 2, calls] {
        let mut current = 0;
        assert_eq!(
            net.step_candidate(&initial, 0.1, strict(), || {
                current += 1;
                current == stop
            }),
            Err(BinaryMiscibilityError::Cancelled)
        );
        assert_eq!(initial, before);
    }
    let short = BinaryMiscibilitySolve {
        maximum_root_iterations: 1,
        ..strict()
    };
    assert!(matches!(
        net.step_candidate(&initial, 0.1, short, || false),
        Err(BinaryMiscibilityError::Root { .. })
    ));
    assert!(
        net.step_candidate(&initial, f64::NAN, strict(), || false)
            .is_err()
    );
    let other = network(
        3.0,
        0.2,
        &[1.0, 1.0, 1.0],
        vec![face(0, 1, 1.0, 1.0), face(1, 2, 1.0, 1.0)],
    );
    assert!(matches!(
        other.step_candidate(&initial, 0.1, strict(), || false),
        Err(BinaryMiscibilityError::Invalid(_))
    ));
    assert!(
        net.state_from_component_moles(vec![[1.0, 1.0]; 3], || false)
            .is_err()
    );
    assert!(matches!(
        BinaryMiscibilityNetwork::new(
            model(0.0, 0.0),
            vec![BinaryMiscibilityCell { volume_m3: 1.0 }; 2],
            vec![face(0, 1, 1.0, 1.0)],
            1,
            || false
        ),
        Err(BinaryMiscibilityError::Workspace { .. })
    ));
}
