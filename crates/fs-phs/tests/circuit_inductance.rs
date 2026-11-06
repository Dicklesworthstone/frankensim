//! Native nonlinear magnetic storage: independent physical trajectories,
//! branch work, reciprocal block composition, and fail-closed source support.

use fs_phs::circuit::capacitance::{CapacitanceBlock, CapacitorChargeCurve, CapacitorCurve};
use fs_phs::circuit::inductance::{
    InductorCurve, InductorFluxCurve, assemble_circuit_with_inductor_curves,
    assemble_circuit_with_storage_curves,
};
use fs_phs::circuit::{
    Branch, CircuitDae, CircuitGraph, MutualInductance, TransformerLink, assemble_circuit,
};
use fs_phs::discrete_gradient;

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual={actual:e}, expected={expected:e}, tolerance={tolerance:e}"
    );
}

fn knots(scale: f64) -> Vec<(f64, f64)> {
    [
        (-2.0, -6.0),
        (-1.0, -1.0),
        (0.0, 0.0),
        (1.0, 1.0),
        (2.0, 6.0),
    ]
    .into_iter()
    .map(|(flux, current)| (scale * flux, current))
    .collect()
}

fn law(scale: f64) -> InductorFluxCurve {
    InductorFluxCurve::try_new(knots(scale)).unwrap()
}

// Independent antiderivative of I(λ)=λ for |λ|<=1 and
// I(λ)=sign(λ)*(5|λ|-4) beyond the saturation knot.
fn energy(flux: f64) -> f64 {
    let x = flux.abs();
    if x <= 1.0 {
        0.5 * x * x
    } else {
        2.5 * x * x - 4.0 * x + 2.0
    }
}

#[test]
fn saturated_energy_is_the_current_integral_and_inductance_is_differential() {
    let curve = law(1.0);
    assert_eq!(curve.knots(), knots(1.0).as_slice());
    assert_eq!(curve.flux_range_wb(), (-2.0, 2.0));
    close(curve.current_a(2.0).unwrap(), 6.0, 0.0);
    close(curve.differential_inductance_h(0.5).unwrap(), 1.0, 0.0);
    close(curve.differential_inductance_h(1.5).unwrap(), 0.2, 0.0);
    close(curve.stored_energy_j(2.0).unwrap(), 4.0, 0.0);
    assert_ne!(curve.stored_energy_j(2.0).unwrap(), 2.0 * 6.0 / 2.0);
    for flux in [-2.0, -1.7, -1.0, -1e-20, 0.0, 1e-20, 0.4, 1.2, 2.0] {
        close(curve.stored_energy_j(flux).unwrap(), energy(flux), 2e-15);
        close(
            curve.flux_wb(curve.current_a(flux).unwrap()).unwrap(),
            flux,
            2e-15 * flux.abs(),
        );
    }
    for scale in [1.0, 1e-9, 1e-180] {
        let curve = law(scale);
        let before = 2.0 * scale;
        let after = 0.5 * scale;
        let mean = curve.discrete_current_a(before, after).unwrap();
        close(mean, 31.0 / 12.0, 2e-14);
        close(
            curve.discrete_current_a(after, before).unwrap(),
            mean,
            2e-14,
        );
        close(
            curve.energy_change_j(before, after).unwrap(),
            -3.875 * scale,
            2e-14 * scale,
        );
        close(mean * (after - before), -3.875 * scale, 2e-14 * scale);
    }
}

fn rl_graph(inductance_h: f64) -> CircuitGraph {
    CircuitGraph {
        node_count: 3,
        branches: vec![
            (1, 0, Branch::VoltageSource { port: 0 }),
            (1, 2, Branch::Resistor { ohms: 1.0 }),
            (
                2,
                0,
                Branch::Inductor {
                    henries: inductance_h,
                },
            ),
        ],
        transformers: vec![],
    }
}

fn initial(dae: &CircuitDae, flux: f64, input: &[f64]) -> Vec<f64> {
    let mut state = vec![0.0; dae.system.state_dim()];
    state[dae.flux_index[0]] = flux;
    dae.consistent_initial_state(&state, input).unwrap()
}

fn driven_rl(steps: usize, scale: f64) -> f64 {
    let curve = law(scale);
    let dae = assemble_circuit_with_inductor_curves(
        &rl_graph(scale),
        &[InductorCurve {
            branch: 2,
            curve: curve.clone(),
        }],
        &[],
    )
    .unwrap();
    let input = [3.0];
    let mut state = initial(&dae, 0.0, &input);
    let dt = scale / steps as f64;
    let mut supplied = 0.0;
    let mut heat = 0.0;
    for _ in 0..steps {
        let (record, defect) = dae.step_audited(&state, &input, dt).unwrap();
        let old_flux = state[dae.flux_index[0]];
        let next_flux = record.x[dae.flux_index[0]];
        let current = curve.discrete_current_a(old_flux, next_flux).unwrap();
        assert!(next_flux > old_flux);
        close(next_flux - old_flux, dt * (3.0 - current), 5e-10 * scale);
        close(record.supplied, dt * 3.0 * current, 5e-10 * scale);
        close(record.dissipated, dt * current * current, 5e-10 * scale);
        assert!(defect < 5e-10 * scale, "defect={defect:e}");
        supplied += record.supplied;
        heat += record.dissipated;
        state = record.x;
    }
    let flux = state[dae.flux_index[0]];
    close(energy(flux / scale) * scale + heat, supplied, 5e-9 * scale);
    close(
        dae.inductor_currents(&state).unwrap()[0],
        5.0 * flux / scale - 4.0,
        5e-10,
    );
    flux / scale
}

#[test]
fn driven_rl_crosses_saturation_with_second_order_accuracy_and_realistic_units() {
    // λ'=3-I(λ). Before tc=ln(3/2), λ=3(1-exp(-t)). After
    // crossing λ=1, λ=1.4-.4*exp(-5(t-tc)). No reference integrator.
    let crossing = 1.5_f64.ln();
    let exact = 1.4 - 0.4 * (-5.0 * (1.0 - crossing)).exp();
    let errors = [40, 80, 160].map(|n| (driven_rl(n, 1.0) - exact).abs());
    assert!(errors[2] < 1e-5, "{errors:?}");
    for pair in errors.windows(2) {
        let ratio = pair[0] / pair[1];
        assert!((3.7..4.3).contains(&ratio), "{errors:?}");
    }
    let reference = driven_rl(80, 1.0);
    // Scale flux and time together: the same physical transient at mH and µH.
    for scale in [1e-3, 1e-6] {
        close(driven_rl(80, scale), reference, 2e-7);
    }
}

#[test]
fn nonlinear_lc_preserves_each_components_work_through_saturation() {
    let flux_law = law(1.0);
    let charge_law = CapacitorChargeCurve::try_new(knots(1.0)).unwrap();
    let dae = assemble_circuit_with_storage_curves(
        &CircuitGraph {
            node_count: 2,
            branches: vec![
                (1, 0, Branch::Inductor { henries: 1.0 }),
                (1, 0, Branch::Capacitor { farads: 1.0 }),
            ],
            transformers: vec![],
        },
        &[],
        &[CapacitorCurve {
            branch: 1,
            curve: charge_law.clone(),
        }],
        &[InductorCurve {
            branch: 0,
            curve: flux_law.clone(),
        }],
        &[],
    )
    .unwrap();
    let mut state = vec![0.0; dae.system.state_dim()];
    state[dae.charge_index[0]] = 1.5;
    state = dae.consistent_initial_state(&state, &[]).unwrap();
    let initial_energy = energy(1.5);
    let mut crossed = false;
    for _ in 0..200 {
        let (record, defect) = dae.step_audited(&state, &[], 0.01).unwrap();
        let (a, b) = (state[dae.flux_index[0]], record.x[dae.flux_index[0]]);
        let (q, r) = (state[dae.charge_index[0]], record.x[dae.charge_index[0]]);
        let current = flux_law.discrete_current_a(a, b).unwrap();
        let voltage = charge_law.discrete_voltage_v(q, r).unwrap();
        close(b - a, 0.01 * voltage, 5e-10);
        close(r - q, -0.01 * current, 5e-10);
        close(current * (b - a) + voltage * (r - q), 0.0, 5e-10);
        close(energy(b) + energy(r), initial_energy, 2e-8);
        assert!(defect < 5e-10, "defect={defect:e}");
        crossed |= b.abs() > 1.0;
        state = record.x;
    }
    assert!(crossed);
}

#[test]
fn nonlinear_flux_and_charge_keep_reciprocal_blocks_and_their_analytic_actions() {
    let graph = CircuitGraph {
        node_count: 7,
        branches: vec![
            (1, 0, Branch::Inductor { henries: 1.0 }),
            (1, 0, Branch::Capacitor { farads: 1.0 }),
            (2, 0, Branch::Inductor { henries: 2.0 }),
            (3, 0, Branch::Inductor { henries: 3.0 }),
            (4, 0, Branch::Capacitor { farads: 1.0 }),
            (5, 0, Branch::Capacitor { farads: 1.0 }),
            (6, 0, Branch::Capacitor { farads: 1.0 }),
        ],
        transformers: vec![],
    };
    let dae = assemble_circuit_with_storage_curves(
        &graph,
        &[CapacitanceBlock::try_new(vec![4, 5], vec![2.0, 1.0, 1.0, 3.0]).unwrap()],
        &[CapacitorCurve {
            branch: 1,
            curve: CapacitorChargeCurve::try_new(knots(1.0)).unwrap(),
        }],
        &[InductorCurve {
            branch: 0,
            curve: law(1.0),
        }],
        &[MutualInductance {
            branches: (2, 3),
            henries: 1.0,
        }],
    )
    .unwrap();
    let mut a = vec![0.0; dae.system.state_dim()];
    let mut b = a.clone();
    for (coordinate, before, after) in [
        (dae.flux_index[0], 2.0, 0.5),
        (dae.flux_index[1], 1.0, 2.5),
        (dae.flux_index[2], -0.5, 1.5),
        (dae.charge_index[0], 0.25, 1.5),
        (dae.charge_index[1], -1.0, 2.0),
        (dae.charge_index[2], 2.0, 1.0),
        (dae.charge_index[3], 1e8, 1e8),
    ] {
        a[coordinate] = before;
        b[coordinate] = after;
    }
    let mean = discrete_gradient(dae.system.storage(), &a, &b);
    close(mean[dae.flux_index[0]], 31.0 / 12.0, 2e-14);
    close(mean[dae.charge_index[0]], 51.0 / 40.0, 2e-14);
    // Both reciprocal matrices have inverse [[3,-1],[-1,2]]/5. They act
    // on their OWN coordinate midpoints despite both changing nonlinear laws.
    close(mean[dae.flux_index[1]], 0.95, 2e-14);
    close(mean[dae.flux_index[2]], -0.15, 2e-14);
    close(mean[dae.charge_index[1]], 0.0, 2e-14);
    close(mean[dae.charge_index[2]], 0.5, 2e-14);
    assert_eq!(mean[dae.charge_index[3]], 1e8);
    close(
        mean.iter()
            .zip(a.iter().zip(&b))
            .map(|(e, (x, y))| e * (y - x))
            .sum(),
        -1.65625,
        4e-14,
    );
    for (actual, expected) in dae
        .inductor_currents(&b)
        .unwrap()
        .iter()
        .zip([0.5, 1.2, 0.1])
    {
        close(*actual, expected, 2e-14);
    }
    let mut direction = vec![0.0; a.len()];
    for (coordinate, value) in [
        (dae.flux_index[0], 0.3),
        (dae.flux_index[1], 0.2),
        (dae.flux_index[2], -0.4),
        (dae.charge_index[0], -0.7),
        (dae.charge_index[1], 0.4),
        (dae.charge_index[2], 0.1),
    ] {
        direction[coordinate] = value;
    }
    let rule = dae.system.storage().discrete_gradient_rule().unwrap();
    let mut action = vec![0.0; a.len()];
    assert!(rule.end_jacobian_action(&a, &b, &direction, &mut action));
    let epsilon = 2e-6;
    let plus: Vec<_> = b
        .iter()
        .zip(&direction)
        .map(|(x, d)| x + epsilon * d)
        .collect();
    let minus: Vec<_> = b
        .iter()
        .zip(&direction)
        .map(|(x, d)| x - epsilon * d)
        .collect();
    let gp = discrete_gradient(dae.system.storage(), &a, &plus);
    let gm = discrete_gradient(dae.system.storage(), &a, &minus);
    for ((p, m), expected) in gp.iter().zip(&gm).zip(&action) {
        close((p - m) / (2.0 * epsilon), *expected, 2e-8);
    }
}

#[test]
fn saturating_magnetizing_branch_and_loaded_ideal_transformer_share_source_work() {
    let dae = assemble_circuit_with_inductor_curves(
        &CircuitGraph {
            node_count: 3,
            branches: vec![
                (1, 0, Branch::VoltageSource { port: 0 }),
                (1, 0, Branch::Inductor { henries: 1.0 }),
                (2, 0, Branch::Resistor { ohms: 4.0 }),
            ],
            transformers: vec![TransformerLink {
                primary: (1, 0),
                secondary: (2, 0),
                ratio: 2.0,
            }],
        },
        &[InductorCurve {
            branch: 1,
            curve: law(1.0),
        }],
        &[],
    )
    .unwrap();
    let state = initial(&dae, 0.0, &[1.5]);
    let (record, defect) = dae.step_audited(&state, &[1.5], 1.0).unwrap();
    // Primary volt-seconds set λ=1.5; secondary voltage is .75 V. The
    // same source supplies the nonlinear magnetic integral and real load heat.
    close(record.x[dae.flux_index[0]], 1.5, 3e-11);
    close(dae.inductor_currents(&record.x).unwrap()[0], 3.5, 3e-11);
    close(record.delta_h, 1.625, 3e-11);
    close(record.dissipated, 0.75 * 0.75 / 4.0, 3e-11);
    close(record.supplied, 1.625 + 0.75 * 0.75 / 4.0, 3e-11);
    assert!(defect < 3e-11, "defect={defect:e}");
}

#[test]
fn constant_flux_law_matches_linear_dynamics_and_empty_bindings_keep_fast_path() {
    let graph = rl_graph(2.0);
    let original = assemble_circuit(&graph).unwrap();
    let empty = assemble_circuit_with_inductor_curves(&graph, &[], &[]).unwrap();
    assert_eq!(
        original.system.storage().exact_quadratic_matrix(),
        empty.system.storage().exact_quadratic_matrix(),
    );
    assert!(empty.system.storage().exact_quadratic_matrix().is_some());
    let bounded = assemble_circuit_with_inductor_curves(
        &graph,
        &[InductorCurve {
            branch: 2,
            curve: InductorFluxCurve::try_new(vec![(-20.0, -10.0), (0.0, 0.0), (20.0, 10.0)])
                .unwrap(),
        }],
        &[],
    )
    .unwrap();
    let mut reference = initial(&original, 0.4, &[3.0]);
    let mut state = initial(&bounded, 0.4, &[3.0]);
    for _ in 0..20 {
        let (expected, _) = original.step_audited(&reference, &[3.0], 0.05).unwrap();
        let (actual, _) = bounded.step_audited(&state, &[3.0], 0.05).unwrap();
        for (a, b) in actual.x.iter().zip(&expected.x) {
            close(*a, *b, 2e-10);
        }
        close(actual.supplied, expected.supplied, 2e-10);
        close(actual.dissipated, expected.dissipated, 2e-10);
        state = actual.x;
        reference = expected.x;
    }
}

#[test]
fn unsupported_scalar_mutual_overlap_invalid_branches_and_source_exits_refuse() {
    for bad in [
        vec![(-1.0, -1.0), (0.0, 0.0), (1.0, 0.0)],
        vec![(-1.0, -1.0), (0.5, 0.5), (1.0, 1.0)],
        vec![(0.0, 0.0), (1.0, 1.0), (2.0, 3.0)],
        vec![(-1.0, -1.0), (0.0, 0.0), (1.0, f64::INFINITY)],
    ] {
        assert!(InductorFluxCurve::try_new(bad).is_err());
    }
    let graph = rl_graph(1.0);
    for branches in [vec![0], vec![3], vec![2, 2]] {
        let bindings: Vec<_> = branches
            .into_iter()
            .map(|branch| InductorCurve {
                branch,
                curve: law(1.0),
            })
            .collect();
        assert!(assemble_circuit_with_inductor_curves(&graph, &bindings, &[]).is_err());
    }
    let graph = CircuitGraph {
        node_count: 3,
        branches: vec![
            (1, 0, Branch::Inductor { henries: 1.0 }),
            (2, 0, Branch::Inductor { henries: 1.0 }),
        ],
        transformers: vec![],
    };
    for pair in [(0, 1), (1, 0)] {
        assert!(
            assemble_circuit_with_inductor_curves(
                &graph,
                &[InductorCurve {
                    branch: 0,
                    curve: law(1.0)
                }],
                &[MutualInductance {
                    branches: pair,
                    henries: 0.2
                }],
            )
            .is_err()
        );
    }
    let curve = law(1.0);
    assert!(curve.current_a(2.01).is_err());
    assert!(curve.flux_wb(6.01).is_err());
    assert!(curve.discrete_current_a(0.0, f64::NAN).is_err());
    let dae = assemble_circuit_with_inductor_curves(
        &CircuitGraph {
            node_count: 2,
            branches: vec![
                (1, 0, Branch::Inductor { henries: 1.0 }),
                (1, 0, Branch::VoltageSource { port: 0 }),
            ],
            transformers: vec![],
        },
        &[InductorCurve { branch: 0, curve }],
        &[],
    )
    .unwrap();
    let state = initial(&dae, 1.9, &[2.0]);
    let saved = state.clone();
    // An ideal 2 V source requests λ=2.1, beyond the supplied flux support.
    assert!(dae.step_audited(&state, &[2.0], 0.1).is_err());
    assert_eq!(state, saved);
}
