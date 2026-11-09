//! Shared-core saturation in original Kirchhoff flux coordinates. The oracles
//! below follow a two-winding magnetic circuit with independent leakage energy;
//! they do not infer a field geometry or claim hysteresis/adjoint certification.

use fs_phs::circuit::inductance::{
    InductorCurve, InductorFluxCurve, MagneticStorageBlock, SharedCoreInductance,
    assemble_circuit_with_storage_blocks,
};
use fs_phs::circuit::{
    Branch, CircuitDae, CircuitGraph, MutualInductance, assemble_circuit_with_mutual_inductance,
};
use fs_phs::discrete_gradient;

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual={actual:e}, expected={expected:e}, tolerance={tolerance:e}"
    );
}

fn core_knots() -> Vec<(f64, f64)> {
    vec![
        (-2.0, -6.0),
        (-1.0, -1.0),
        (0.0, 0.0),
        (1.0, 1.0),
        (2.0, 6.0),
    ]
}

fn law() -> SharedCoreInductance {
    SharedCoreInductance::try_new(vec![2.0, -1.0], vec![0.5, 0.25], core_knots()).unwrap()
}

// MMF(phi)=phi inside the knee and sign(phi)*(5|phi|-4) outside.
// Ampere-turn balance with the two leakage stores gives
// MMF(phi)+12phi=4(lambda_0-lambda_1), a scalar physical reference.
fn reference(flux: [f64; 2]) -> (f64, [f64; 2], f64) {
    let t = 4.0 * (flux[0] - flux[1]);
    let phi = if t.abs() <= 13.0 {
        t / 13.0
    } else {
        (t + 4.0 * t.signum()) / 17.0
    };
    let current = [2.0 * flux[0] - 4.0 * phi, 4.0 * flux[1] + 4.0 * phi];
    let x = phi.abs();
    let core_energy = if x <= 1.0 {
        0.5 * x * x
    } else {
        2.5 * x * x - 4.0 * x + 2.0
    };
    let energy = core_energy + 0.25 * current[0].powi(2) + 0.125 * current[1].powi(2);
    (phi, current, energy)
}

// Antiderivative of the independently solved phi(t); a straight flux-linkage
// path makes t affine, so this is exact even across both saturation knees.
fn integrated_phi(t: f64) -> f64 {
    let x = t.abs();
    if x <= 13.0 {
        x * x / 26.0
    } else {
        (0.5 * x * x + 4.0 * x - 26.0) / 17.0
    }
}

fn mean_current(before: [f64; 2], after: [f64; 2]) -> [f64; 2] {
    let t0 = 4.0 * (before[0] - before[1]);
    let t1 = 4.0 * (after[0] - after[1]);
    let phi = if t0 == t1 {
        reference(before).0
    } else {
        (integrated_phi(t1) - integrated_phi(t0)) / (t1 - t0)
    };
    [
        before[0] + after[0] - 4.0 * phi,
        2.0 * (before[1] + after[1]) + 4.0 * phi,
    ]
}

fn driven_graph() -> CircuitGraph {
    CircuitGraph {
        node_count: 3,
        branches: vec![
            (1, 0, Branch::Inductor { henries: 4.5 }),
            (1, 0, Branch::VoltageSource { port: 0 }),
            (2, 0, Branch::Inductor { henries: 1.25 }),
            (2, 0, Branch::VoltageSource { port: 1 }),
        ],
        transformers: vec![],
    }
}

fn assemble(graph: &CircuitGraph, branches: Vec<usize>, law: SharedCoreInductance) -> CircuitDae {
    assemble_circuit_with_storage_blocks(
        graph,
        &[],
        &[],
        &[],
        &[MagneticStorageBlock { branches, law }],
        &[],
    )
    .unwrap()
}

fn state(dae: &CircuitDae, flux: &[f64], input: &[f64]) -> Vec<f64> {
    let mut seed = vec![0.0; dae.system.state_dim()];
    for (&index, &value) in dae.flux_index.iter().zip(flux) {
        seed[index] = value;
    }
    dae.consistent_initial_state(&seed, input).unwrap()
}

#[test]
fn opposed_ampere_turns_preserve_leakage_while_common_excitation_saturates() {
    let law = law();
    assert_eq!(law.turns(), &[2.0, -1.0]);
    assert_eq!(law.leakage_inductance_h(), &[0.5, 0.25]);
    assert_eq!(law.core_flux_mmf_knots(), core_knots().as_slice());
    // 2*1-1*2=0: a nonzero differential current stores leakage energy but
    // neither magnetizes nor saturates the shared core.
    let flux = law.flux_linkages_wb(&[1.0, 2.0]).unwrap();
    assert_eq!(flux, vec![0.5, 0.5]);
    near(law.common_flux_wb(&flux).unwrap(), 0.0, 1e-15);
    near(law.stored_energy_j(&flux).unwrap(), 0.75, 1e-14);
    // Large differential linkage remains supported: core support constrains
    // ampere-turns, not each winding's total flux independently.
    let differential = [100.0, 100.0];
    near(law.common_flux_wb(&differential).unwrap(), 0.0, 1e-13);
    assert_eq!(law.currents_a(&differential).unwrap(), vec![200.0, 400.0]);
    for flux in [
        [0.5, 0.5],
        [1.1, -0.525],
        [4.6, -2.05],
        [-4.6, 2.05],
        [3.0, -2.0],
    ] {
        let (phi, currents, energy) = reference(flux);
        near(law.common_flux_wb(&flux).unwrap(), phi, 3e-14);
        near(law.stored_energy_j(&flux).unwrap(), energy, 3e-13);
        let actual = law.currents_a(&flux).unwrap();
        for (i, expected) in actual.iter().zip(currents) {
            near(*i, expected, 3e-14);
        }
        for (actual, expected) in law.flux_linkages_wb(&actual).unwrap().iter().zip(flux) {
            near(*actual, expected, 3e-13);
        }
    }
    let unsaturated = law.differential_inductance_h(&[1.1, -0.525]).unwrap();
    let saturated = law.differential_inductance_h(&[4.6, -2.05]).unwrap();
    for (actual, expected) in unsaturated.iter().zip([4.5, -2.0, -2.0, 1.25]) {
        near(*actual, expected, 2e-14);
    }
    for (actual, expected) in saturated.iter().zip([1.3, -0.4, -0.4, 0.45]) {
        near(*actual, expected, 2e-14);
    }
}

#[test]
fn millihenry_and_microhenry_models_keep_currents_and_finite_step_magnetic_work() {
    let before = [-3.0, 2.0];
    let after = [4.5, -2.5];
    let expected_current = mean_current(before, after);
    let expected_work = reference(after).2 - reference(before).2;
    for scale in [1e-3, 1e-6] {
        let source = SharedCoreInductance::try_new(
            vec![2.0, -1.0],
            vec![0.5 * scale, 0.25 * scale],
            core_knots()
                .into_iter()
                .map(|(phi, mmf)| (scale * phi, mmf))
                .collect(),
        )
        .unwrap();
        let a = before.map(|v| v * scale);
        let b = after.map(|v| v * scale);
        for (actual, expected) in source
            .discrete_currents_a(&a, &b)
            .unwrap()
            .iter()
            .zip(expected_current)
        {
            near(*actual, expected, 3e-13);
        }
        near(
            source.energy_change_j(&a, &b).unwrap(),
            expected_work * scale,
            3e-13 * scale,
        );
        near(
            source.stored_energy_j(&b).unwrap(),
            reference(after).2 * scale,
            3e-13 * scale,
        );
    }
}

#[test]
fn tight_coupling_retains_common_mode_current_energy_work_and_native_tangent() {
    for turns in [1.0, 100.0] {
        for leakage_ratio in [1e-6, 1e-8, 1e-10] {
            let leakage = turns * turns * leakage_ratio;
            let source = SharedCoreInductance::try_new(
                vec![turns, turns],
                vec![leakage, leakage],
                vec![(-2.0, -2.0), (0.0, 0.0), (2.0, 2.0)],
            )
            .unwrap();
            let flux = [turns, turns];
            // The common eigenvalue of diag(L_leak)+N*N^T is L_leak+2N².
            // This reference never subtracts the nearly equal total/core fluxes.
            let expected = turns / (leakage + 2.0 * turns * turns);
            for current in source.currents_a(&flux).unwrap() {
                near(current, expected, 1e-11 * expected);
            }
            let energy = turns * expected;
            near(
                source.stored_energy_j(&flux).unwrap(),
                energy,
                1e-12 * energy,
            );
            let before = flux.map(|v| 0.8 * v);
            for mean in source.discrete_currents_a(&before, &flux).unwrap() {
                near(mean, 0.9 * expected, 1e-11 * expected);
            }
            near(
                source.energy_change_j(&before, &flux).unwrap(),
                0.36 * energy,
                1e-12 * energy,
            );
            let dae = assemble(&driven_graph(), vec![0, 2], source);
            let mut x = vec![0.0; dae.system.state_dim()];
            for &coordinate in &dae.flux_index {
                x[coordinate] = turns;
            }
            let mut action = vec![f64::NAN; x.len()];
            assert!(
                dae.system
                    .storage()
                    .continuous_tangent()
                    .unwrap()
                    .apply(&x, &x, &mut action)
            );
            for &coordinate in &dae.flux_index {
                near(action[coordinate], expected, 1e-11 * expected);
            }
        }
    }
}

#[test]
fn energy_gradient_and_reciprocal_tangent_follow_the_same_saturated_potential() {
    let law = law();
    let dae = assemble(&driven_graph(), vec![0, 2], law.clone());
    let storage = dae.system.storage();
    let tangent = storage
        .continuous_tangent()
        .expect("continuous magnetic tangent");
    let epsilon = 2e-6;
    for flux in [[1.1, -0.525], [4.6, -2.05], [-4.6, 2.05]] {
        let currents = law.currents_a(&flux).unwrap();
        let inductance = law.differential_inductance_h(&flux).unwrap();
        assert_eq!(inductance[1], inductance[2]);
        let mut x = vec![0.0; dae.system.state_dim()];
        for (&slot, value) in dae.flux_index.iter().zip(flux) {
            x[slot] = value;
        }
        for column in 0..2 {
            let (mut plus, mut minus) = (flux, flux);
            plus[column] += epsilon;
            minus[column] -= epsilon;
            near(
                (law.stored_energy_j(&plus).unwrap() - law.stored_energy_j(&minus).unwrap())
                    / (2.0 * epsilon),
                currents[column],
                2e-8,
            );
            let mut direction = vec![0.0; x.len()];
            direction[dae.flux_index[column]] = 1.0;
            let mut action = vec![f64::NAN; x.len()];
            assert!(tangent.apply(&x, &direction, &mut action));
            let cp = law.currents_a(&plus).unwrap();
            let cm = law.currents_a(&minus).unwrap();
            for row in 0..2 {
                near(
                    action[dae.flux_index[row]],
                    (cp[row] - cm[row]) / (2.0 * epsilon),
                    2e-9,
                );
            }
            assert!(
                action[dae.system.differential_dim()..]
                    .iter()
                    .all(|&v| v == 0.0)
            );
            let (mut ip, mut im) = (currents.clone(), currents.clone());
            ip[column] += epsilon;
            im[column] -= epsilon;
            let lp = law.flux_linkages_wb(&ip).unwrap();
            let lm = law.flux_linkages_wb(&im).unwrap();
            for row in 0..2 {
                near(
                    (lp[row] - lm[row]) / (2.0 * epsilon),
                    inductance[row * 2 + column],
                    2e-9,
                );
            }
        }
    }
}

#[test]
fn low_word_saturation_knee_selects_the_physical_side_and_refuses_only_the_exact_kink() {
    let law = law();
    let dae = assemble(&driven_graph(), vec![0, 2], law.clone());
    let tangent = dae.system.storage().continuous_tangent().unwrap();
    let offset = 2.0_f64.powi(-55);
    for (second_flux, inductance, action) in [
        (offset, [4.5, -2.0, -2.0, 1.25], [10.0 / 13.0, 16.0 / 13.0]),
        (-offset, [1.3, -0.4, -0.4, 0.45], [18.0 / 17.0, 16.0 / 17.0]),
    ] {
        let flux = [3.25, second_flux];
        // Phi rounds to the knee although the represented linkages' exact
        // ampere-turn drive lies on opposite sides: 13 minus/plus 2^-53.
        assert_eq!(law.common_flux_wb(&flux).unwrap(), 1.0);
        for (actual, expected) in law
            .differential_inductance_h(&flux)
            .unwrap()
            .iter()
            .zip(inductance)
        {
            near(*actual, expected, 2e-14);
        }
        let mut x = vec![0.0; dae.system.state_dim()];
        x[dae.flux_index[0]] = flux[0];
        x[dae.flux_index[1]] = flux[1];
        let mut direction = vec![0.0; x.len()];
        direction[dae.flux_index[0]] = 1.0;
        let mut actual = vec![f64::NAN; x.len()];
        assert!(tangent.apply(&x, &direction, &mut actual));
        for (&coordinate, expected) in dae.flux_index.iter().zip(action) {
            near(actual[coordinate], expected, 2e-14);
        }
        x[dae.flux_index[1]] = 0.0;
        assert!(!tangent.apply(&x, &direction, &mut actual));
    }
}

#[test]
fn linear_shared_core_matches_existing_mutual_rl_dynamics_and_losses() {
    let graph = CircuitGraph {
        node_count: 5,
        branches: vec![
            (1, 0, Branch::VoltageSource { port: 0 }),
            (1, 2, Branch::Resistor { ohms: 1.0 }),
            (2, 0, Branch::Inductor { henries: 4.5 }),
            (3, 0, Branch::VoltageSource { port: 1 }),
            (3, 4, Branch::Resistor { ohms: 2.0 }),
            (4, 0, Branch::Inductor { henries: 1.25 }),
        ],
        transformers: vec![],
    };
    let native = assemble_circuit_with_mutual_inductance(
        &graph,
        &[MutualInductance {
            branches: (2, 5),
            henries: -2.0,
        }],
    )
    .unwrap();
    let linear = SharedCoreInductance::try_new(
        vec![2.0, -1.0],
        vec![0.5, 0.25],
        vec![(-100.0, -100.0), (0.0, 0.0), (100.0, 100.0)],
    )
    .unwrap();
    let shared = assemble(&graph, vec![2, 5], linear);
    let input = [1.0, -0.5];
    let mut expected = state(&native, &[0.3, -0.2], &input);
    let mut actual = state(&shared, &[0.3, -0.2], &input);
    for _ in 0..20 {
        let (a, defect) = shared.step_audited(&actual, &input, 0.025).unwrap();
        let (b, _) = native.step_audited(&expected, &input, 0.025).unwrap();
        for (x, y) in a.x.iter().zip(&b.x) {
            near(*x, *y, 2e-9);
        }
        near(a.delta_h, b.delta_h, 2e-10);
        near(a.supplied, b.supplied, 2e-10);
        near(a.dissipated, b.dissipated, 2e-10);
        assert!(defect < 2e-9);
        actual = a.x;
        expected = b.x;
    }
}

#[test]
fn native_source_work_crosses_both_saturation_knees_beside_held_storage() {
    let mut graph = driven_graph();
    graph.node_count = 4;
    graph.branches.extend([
        (3, 0, Branch::Inductor { henries: 2.0 }),
        (3, 0, Branch::VoltageSource { port: 2 }),
    ]);
    let dae = assemble(&graph, vec![0, 2], law());
    let before = [-3.0, 2.0];
    let after = [4.5, -2.5];
    let dt = 0.5;
    let input = [
        (after[0] - before[0]) / dt,
        (after[1] - before[1]) / dt,
        0.0,
    ];
    let initial = state(&dae, &[before[0], before[1], 1e8], &input);
    let (step, defect) = dae.step_audited(&initial, &input, dt).unwrap();
    let mean = mean_current(before, after);
    for j in 0..2 {
        near(step.x[dae.flux_index[j]], after[j], 3e-9);
        near(step.y[j], mean[j], 3e-9);
    }
    assert_eq!(step.x[dae.flux_index[2]], initial[dae.flux_index[2]]);
    let work = mean[0] * (after[0] - before[0]) + mean[1] * (after[1] - before[1]);
    let energy_change = reference(after).2 - reference(before).2;
    near(work, energy_change, 5e-14);
    near(step.delta_h, energy_change, 3e-9);
    near(step.supplied, energy_change, 3e-9);
    near(step.dissipated, 0.0, 1e-14);
    assert!(defect < 3e-9);
    // Analytic finite-step Jacobian must differentiate this same branch work,
    // including the path's crossings of both core knees.
    let storage = dae.system.storage();
    let rule = storage.discrete_gradient_rule().unwrap();
    let mut direction = vec![0.0; initial.len()];
    direction[dae.flux_index[0]] = 0.7;
    direction[dae.flux_index[1]] = -0.3;
    let mut action = vec![f64::NAN; initial.len()];
    assert!(rule.end_jacobian_action(&initial, &step.x, &direction, &mut action));
    let epsilon = 2e-6;
    let plus: Vec<_> = step
        .x
        .iter()
        .zip(&direction)
        .map(|(x, d)| x + epsilon * d)
        .collect();
    let minus: Vec<_> = step
        .x
        .iter()
        .zip(&direction)
        .map(|(x, d)| x - epsilon * d)
        .collect();
    let gp = discrete_gradient(storage, &initial, &plus);
    let gm = discrete_gradient(storage, &initial, &minus);
    for j in 0..2 {
        let slot = dae.flux_index[j];
        near(action[slot], (gp[slot] - gm[slot]) / (2.0 * epsilon), 3e-9);
    }
}

#[test]
fn winding_order_and_terminal_orientation_change_coordinates_not_physics() {
    let graph = driven_graph();
    let original = assemble(&graph, vec![0, 2], law());
    let reordered = assemble(
        &graph,
        vec![2, 0],
        SharedCoreInductance::try_new(vec![-1.0, 2.0], vec![0.25, 0.5], core_knots()).unwrap(),
    );
    let mut reversed_graph = graph;
    reversed_graph.branches[0].0 = 0;
    reversed_graph.branches[0].1 = 1;
    let reversed = assemble(
        &reversed_graph,
        vec![0, 2],
        SharedCoreInductance::try_new(vec![-2.0, -1.0], vec![0.5, 0.25], core_knots()).unwrap(),
    );
    let input = [1.0, -0.5];
    let a = state(&original, &[4.0, -2.0], &input);
    let b = state(&reordered, &[4.0, -2.0], &input);
    let c = state(&reversed, &[-4.0, -2.0], &input);
    let (a, _) = original.step_audited(&a, &input, 0.1).unwrap();
    let (b, _) = reordered.step_audited(&b, &input, 0.1).unwrap();
    let (c, _) = reversed.step_audited(&c, &input, 0.1).unwrap();
    for (x, y) in a.x.iter().zip(&b.x) {
        near(*x, *y, 2e-9);
    }
    near(
        a.x[original.flux_index[0]],
        -c.x[reversed.flux_index[0]],
        2e-9,
    );
    near(
        a.x[original.flux_index[1]],
        c.x[reversed.flux_index[1]],
        2e-9,
    );
    for (x, y) in a.y.iter().zip(&c.y) {
        near(*x, *y, 2e-9);
    }
    near(a.delta_h, c.delta_h, 2e-9);
    near(a.supplied, c.supplied, 2e-9);
}

#[test]
fn unsupported_leakage_overlap_and_core_support_refuse_without_state_mutation() {
    for leakage in [
        vec![0.0, 0.25],
        vec![-0.1, 0.25],
        vec![0.5],
        vec![f64::NAN, 0.25],
    ] {
        assert!(SharedCoreInductance::try_new(vec![2.0, -1.0], leakage, core_knots()).is_err());
    }
    let law = law();
    assert!(law.currents_a(&[8.0, 0.0]).is_err());
    assert!(law.flux_linkages_wb(&[4.0, 0.0]).is_err());
    let graph = driven_graph();
    for branches in [vec![0, 0], vec![0], vec![0, 1], vec![0, 9]] {
        assert!(
            assemble_circuit_with_storage_blocks(
                &graph,
                &[],
                &[],
                &[],
                &[MagneticStorageBlock {
                    branches,
                    law: law.clone()
                }],
                &[],
            )
            .is_err()
        );
    }
    let block = MagneticStorageBlock {
        branches: vec![0, 2],
        law: law.clone(),
    };
    let scalar = InductorCurve {
        branch: 0,
        curve: InductorFluxCurve::try_new(core_knots()).unwrap(),
    };
    assert!(
        assemble_circuit_with_storage_blocks(
            &graph,
            &[],
            &[],
            &[scalar],
            std::slice::from_ref(&block),
            &[],
        )
        .is_err()
    );
    assert!(
        assemble_circuit_with_storage_blocks(
            &graph,
            &[],
            &[],
            &[],
            std::slice::from_ref(&block),
            &[MutualInductance {
                branches: (0, 2),
                henries: 0.1
            }],
        )
        .is_err()
    );
    assert!(
        assemble_circuit_with_storage_blocks(&graph, &[], &[], &[], &[block.clone(), block], &[],)
            .is_err()
    );
    let dae = assemble(&graph, vec![0, 2], law);
    let input = [2.0, 0.0];
    let before = state(&dae, &[7.4, 0.0], &input);
    let saved = before.clone();
    // Volt-seconds would move lambda_0 to 7.6, beyond the common core's
    // supported t=30 limit. Neither extrapolation nor silent clamping is valid.
    assert!(dae.step_audited(&before, &input, 0.1).is_err());
    assert_eq!(before, saved);
}
