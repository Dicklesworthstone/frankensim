//! Full host tests: original exponential mechanics, source integrals and charts.
use super::*;
use fs_couple::render::plate::impact::exponential::ExponentialImpactSystem;
use fs_couple::render::plate::impact::material::heating::exponential::HeatedExponentialImpactSystem;
use fs_phs::linear_exact::{DenseExponentialBudget, MAX_LINEAR_EXPONENTIAL_BLOCK_DIM};

fn budget() -> DenseExponentialBudget {
    DenseExponentialBudget {
        max_dimension: MAX_LINEAR_EXPONENTIAL_BLOCK_DIM,
        max_dense_bytes: usize::MAX,
        max_multiply_terms: u64::MAX,
    }
}
fn exact(damping: f64, gate: &CancelGate) -> ExponentialImpactSystem {
    let body = ImpactBody {
        potential: BodyPotential::Linear(vec![3.0]),
        initial: vec![ModalAcousticState {
            displacement_m_sqrt_kg: 0.3,
            velocity_m_sqrt_kg_per_s: 0.4,
        }],
        damping_per_s: vec![damping],
    };
    let contact = Obstacle::new(
        vec![1.0],
        1,
        1,
        vec![0.0],
        vec![1.0],
        20.0,
        3.0,
        "synthetic quartic obstacle".into(),
    )
    .unwrap();
    ExponentialImpactSystem::new(
        vec![body],
        vec![contact],
        vec![],
        ImpactConfig {
            dt_s: DT,
            max_steps: 1000,
            maximum_energy_j: 20.0,
            energy_absolute_tolerance_j: 1e-11,
            energy_relative_tolerance: 1e-8,
            maximum_generalized_force: 10.0,
        },
        gate,
    )
    .unwrap()
    .with_relaxation_branches(
        vec![
            RelaxationBranch {
                projection: vec![1.0, 0.0],
                stiffness: 4.0,
                relaxation_time_s: TAU[0],
            },
            RelaxationBranch {
                projection: vec![-1.0, 0.0],
                stiffness: 1.0,
                relaxation_time_s: TAU[1],
            },
        ],
        InitialMemory::Unrelaxed,
        2,
        gate,
    )
    .unwrap()
}
fn thermal<'a>(
    c: &'a EquilibriumEnthalpyPhaseCurve,
    gate: &CancelGate,
) -> HeatedExponentialImpactSystem<'a> {
    HeatedExponentialImpactSystem::new(
        exact(0.4, gate),
        stores(c),
        vec![route(1, "support"), route(0, "body")],
        1e-10,
        budget(),
        gate,
    )
    .unwrap()
}

// Independent accepted held-effort trajectory. It keeps all four physical
// coordinates and integrates both source squares; no production heat helper.
fn trajectory_heat(a: &[f64], b: &[f64]) -> [f64; 2] {
    let mid: [f64; 4] = std::array::from_fn(|i| 0.5 * (a[i] + b[i]));
    let delta: [f64; 4] = std::array::from_fn(|i| b[i] - a[i]);
    assert!(a[0] > 0.0 && b[0] > 0.0);
    let correction = 5.0 * mid[0] * delta[0].powi(3) / delta.iter().map(|d| d * d).sum::<f64>();
    let held: [f64; 4] = std::array::from_fn(|i| {
        correction * delta[i] + if i == 0 { 20.0 * mid[0].powi(3) } else { 0.0 }
    });
    let rhs = |x: &[f64; 6]| {
        let e = [
            14.0 * x[0] - 2.0 * x[2] + x[3] + held[0],
            x[1] + held[1],
            x[2] - 2.0 * x[0] + held[2],
            x[3] + x[0] + held[3],
        ];
        [
            e[1],
            -e[0] - 0.4 * e[1],
            -e[2] / TAU[0],
            -e[3] / TAU[1],
            e[2] * e[2] / TAU[0],
            e[3] * e[3] / TAU[1],
        ]
    };
    let mut state = [a[0], a[1], a[2], a[3], 0.0, 0.0];
    let dt = DT / 64.0;
    for _ in 0..64 {
        let k1 = rhs(&state);
        let k2 = rhs(&std::array::from_fn(|i| state[i] + 0.5 * dt * k1[i]));
        let k3 = rhs(&std::array::from_fn(|i| state[i] + 0.5 * dt * k2[i]));
        let k4 = rhs(&std::array::from_fn(|i| state[i] + dt * k3[i]));
        for i in 0..6 {
            state[i] += dt * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]) / 6.0;
        }
    }
    for i in 0..4 {
        near(state[i], b[i], 2e-12);
    }
    [state[4], state[5]]
}

#[test]
fn original_trajectory_heats_material_stores_without_changing_impact_motion() {
    let curve = curve(1000.0);
    let gate = CancelGate::new_clock_free();
    let mut model = thermal(&curve, &gate);
    let mut plain = exact(0.4, &gate);
    let initial = plain.stored_energy_j();
    let mut external = 0.0;
    let mut midpoint_difference = 0.0;
    for _ in 0..80 {
        let before = plain.state().to_vec();
        let received: [f64; 2] = std::array::from_fn(|i| model.reservoirs()[i].received_heat_j());
        let expected = plain.step(&[0.0], &gate).unwrap();
        let frame = model.step(&[0.0], &gate).unwrap();
        assert_eq!(model.mechanics().state(), plain.state());
        assert_eq!(
            frame.mechanical.dissipated_energy_j,
            expected.dissipated_energy_j
        );
        let heat = trajectory_heat(&before, plain.state());
        for i in 0..2 {
            near(
                model.reservoirs()[i].received_heat_j() - received[i],
                heat[i],
                5e-13,
            );
        }
        let midq = 0.5 * (before[0] + plain.state()[0]);
        let midz = 0.5 * (before[2] + plain.state()[2]);
        midpoint_difference += (heat[0] - DT * (midz - 2.0 * midq).powi(2) / TAU[0]).abs();
        near(frame.deposited_heat_j, heat.iter().sum(), 1e-12);
        near(frame.combined_energy_residual_j, 0.0, 1e-10);
        external += frame.remaining_mechanical_loss_j;
    }
    assert!(
        midpoint_difference > 1e-7,
        "must distinguish actual integrated heat from midpoint power"
    );
    assert!(external > 1e-3);
    let stored: f64 = model
        .reservoirs()
        .iter()
        .map(|r| r.storage_change_j())
        .sum();
    near(
        model.mechanics().stored_energy_j() + stored + external,
        initial,
        1e-9,
    );
    for store in model.reservoirs() {
        near(
            store.state().temperature_k() - 300.0,
            store.received_heat_j() / (MASS * 10.0),
            1e-11,
        );
    }
}

#[test]
fn heat_attachment_and_continuation_preserve_current_history_clock_and_routes() {
    let curve = curve(1000.0);
    let gate = CancelGate::new_clock_free();
    let mut a = thermal(&curve, &gate);
    let mut b = thermal(&curve, &gate);
    for _ in 0..13 {
        a.step(&[0.0], &gate).unwrap();
        b.step(&[0.0], &gate).unwrap();
    }
    let (mechanics, stores) = a.into_parts();
    assert_eq!(mechanics.samples(), 13);
    let mut a = HeatedExponentialImpactSystem::new(
        mechanics,
        stores,
        vec![route(0, "body"), route(1, "support")],
        1e-10,
        budget(),
        &gate,
    )
    .unwrap();
    for _ in 0..17 {
        a.step(&[0.0], &gate).unwrap();
        b.step(&[0.0], &gate).unwrap();
    }
    assert_eq!(a.mechanics().state(), b.mechanics().state());
    assert_eq!(a.reservoirs(), b.reservoirs());
    assert_eq!(a.mechanics().samples(), 30);
    let mut live = exact(0.4, &gate);
    for _ in 0..4 {
        live.step(&[0.0], &gate).unwrap();
    }
    let state = live.state().to_vec();
    let attached = HeatedExponentialImpactSystem::new(
        live,
        super::stores(&curve),
        vec![route(0, "body")],
        1e-10,
        budget(),
        &gate,
    )
    .unwrap();
    assert_eq!(attached.mechanics().state(), state);
    assert_eq!(attached.mechanics().samples(), 4);
    assert_eq!(attached.reservoirs()[0].received_heat_j(), 0.0);
}

#[test]
fn source_chart_solver_and_cancellation_refusals_are_atomic_and_retryable() {
    let curve = curve(1000.0);
    let tiny = super::curve(1e-14);
    let gate = CancelGate::new_clock_free();
    let mut bounded = thermal(&tiny, &gate);
    let motion = bounded.mechanics().state().to_vec();
    let heat = bounded.reservoirs().to_vec();
    assert!(bounded.step(&[0.0], &gate).is_err());
    assert_eq!(bounded.mechanics().state(), motion);
    assert_eq!(bounded.reservoirs(), heat);
    assert_eq!(bounded.mechanics().samples(), 0);
    let mut mechanics = exact(0.4, &gate);
    mechanics.set_iteration_limit(0).unwrap();
    let mut model = HeatedExponentialImpactSystem::new(
        mechanics,
        stores(&curve),
        vec![route(0, "body"), route(1, "support")],
        1e-10,
        budget(),
        &gate,
    )
    .unwrap();
    let mut fresh = thermal(&curve, &gate);
    assert!(
        model.step(&[0.0], &gate).is_err(),
        "attachment must preserve Newton allowance"
    );
    let cancelled = CancelGate::new_clock_free();
    cancelled.request();
    assert!(model.step(&[0.0], &cancelled).is_err());
    assert_eq!(model.mechanics().state(), fresh.mechanics().state());
    assert_eq!(model.reservoirs(), fresh.reservoirs());
    model.set_iteration_limit(50).unwrap();
    model.step(&[0.0], &gate).unwrap();
    fresh.step(&[0.0], &gate).unwrap();
    assert_eq!(model.mechanics().state(), fresh.mechanics().state());
    assert_eq!(model.reservoirs(), fresh.reservoirs());
    for routes in [
        vec![route(0, "body"), route(0, "support")],
        vec![route(2, "body")],
        vec![route(0, "missing")],
    ] {
        assert!(
            HeatedExponentialImpactSystem::new(
                exact(0.4, &gate),
                stores(&curve),
                routes,
                1e-10,
                budget(),
                &gate
            )
            .is_err()
        );
    }
    let mut small = budget();
    small.max_dense_bytes = 1;
    assert!(
        HeatedExponentialImpactSystem::new(
            exact(0.4, &gate),
            stores(&curve),
            vec![route(0, "body")],
            1e-10,
            small,
            &gate
        )
        .is_err()
    );
}
