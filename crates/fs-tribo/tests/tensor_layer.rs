//! Manufactured tensor materials verify operators, not measured rosin or wood.
use fs_tribo::layer_contact::tensor::{ShearStiffness, TensorLayerRelaxationArm};
use fs_tribo::layer_contact::{ThirdBodyLayerContact, ThirdBodyLayerRelaxationArm};
use fs_tribo::partial_slip::TangentFrame;
use fs_tribo::tangential::transport::{
    ObjectiveTangentialLaw, TangentialHistoryRotation, TangentialTransportConfig,
};
use fs_tribo::{TangentialLaw, TangentialSolveConfig, TangentialStepInput};

fn input(slip: [f64; 2], dt: f64, compliance: f64) -> TangentialStepInput {
    TangentialStepInput {
        normal_n: 10.0,
        slip_increment_m: slip,
        dt_s: dt,
        contact_temperature_k: Some(300.0),
        nominal_pressure_pa: None,
        compliance_m_per_n: compliance,
        solve: TangentialSolveConfig {
            max_iterations: 96,
            maximum_force_n: 10000.0,
            force_absolute_tolerance_n: 1e-12,
            force_relative_tolerance: 1e-13,
        },
    }
}
fn close(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() <= tolerance, "{a:e} != {b:e}");
}
fn product(k: [[f64; 2]; 2], x: [f64; 2]) -> [f64; 2] {
    [
        k[0][0] * x[0] + k[0][1] * x[1],
        k[1][0] * x[0] + k[1][1] * x[1],
    ]
}
fn solve(a: [[f64; 2]; 2], b: [f64; 2]) -> [f64; 2] {
    let d = a[0][0] * a[1][1] - a[0][1] * a[1][0];
    [
        (a[1][1] * b[0] - a[0][1] * b[1]) / d,
        (a[0][0] * b[1] - a[1][0] * b[0]) / d,
    ]
}
fn model(yield_force: f64) -> ThirdBodyLayerContact {
    ThirdBodyLayerContact::new_tensor(
        ShearStiffness::positive_definite([[800.0, 70.0], [70.0, 500.0]]).unwrap(),
        vec![
            TensorLayerRelaxationArm {
                stiffness: ShearStiffness::from_gram_rows(vec![[20.0, 4.0], [0.0, 7.0]]).unwrap(),
                relaxation_time_s: 0.02,
            },
            TensorLayerRelaxationArm {
                stiffness: ShearStiffness::from_gram_rows(vec![[8.0, -11.0]]).unwrap(),
                relaxation_time_s: 0.07,
            },
        ],
        yield_force,
        0.7,
        300.0,
        1.0,
    )
    .unwrap()
}

#[test]
fn g1_full_tensor_cross_terms_match_independent_equilibrium_with_scalar_body_compliance() {
    let k = [[800.0, 70.0], [70.0, 500.0]];
    let law = ThirdBodyLayerContact::new_tensor(
        ShearStiffness::positive_definite(k).unwrap(),
        vec![],
        1000.0,
        1.0,
        300.0,
        1.0,
    )
    .unwrap();
    let x = input([0.003, -0.001], 0.001, 0.0002);
    let f = product(k, x.slip_increment_m).map(|v| 0.5 * v);
    let a = [
        [
            1.0 + 0.5 * x.compliance_m_per_n * k[0][0],
            0.5 * x.compliance_m_per_n * k[0][1],
        ],
        [
            0.5 * x.compliance_m_per_n * k[1][0],
            1.0 + 0.5 * x.compliance_m_per_n * k[1][1],
        ],
    ];
    let expected = solve(a, f);
    let step = law.advance(&law.initial_state(), &x).unwrap();
    for i in 0..2 {
        close(step.traction_n[i], expected[i], 3e-15);
    }
    close(step.closure_residual_j, 0.0, 1e-17);
    assert_eq!(step.dissipated_j, 0.0);
    assert!(law.stiffness_n_m().is_none());
    assert!(law.relaxation().is_none());
}

#[test]
fn g1_tensor_isotropic_limit_preserves_original_return_and_all_memory_on_nonproportional_paths() {
    for y in [100.0, 0.1] {
        let scalar = ThirdBodyLayerContact::new(800.0, y, 0.7, 300.0, 1.0)
            .unwrap()
            .with_relaxation(vec![ThirdBodyLayerRelaxationArm {
                stiffness_n_m: 400.0,
                relaxation_time_s: 0.02,
            }])
            .unwrap();
        let tensor = ThirdBodyLayerContact::new_tensor(
            ShearStiffness::positive_definite([[800.0, 0.0], [0.0, 800.0]]).unwrap(),
            vec![TensorLayerRelaxationArm {
                stiffness: ShearStiffness::from_gram_rows(vec![[20.0, 0.0], [0.0, 20.0]]).unwrap(),
                relaxation_time_s: 0.02,
            }],
            y,
            0.7,
            300.0,
            1.0,
        )
        .unwrap();
        let (mut a, mut b) = (scalar.initial_state(), tensor.initial_state());
        for i in 0..100 {
            let angle = 0.17 * i as f64;
            let x = input([0.0004 * angle.cos(), 0.0003 * angle.sin()], 0.0005, 0.0001);
            let sa = scalar.advance(&a, &x).unwrap();
            let sb = tensor.advance(&b, &x).unwrap();
            for j in 0..2 {
                close(sa.traction_n[j], sb.traction_n[j], 4e-12);
                close(
                    sa.state.elastic_slip_m()[j],
                    sb.state.elastic_slip_m()[j],
                    5e-15,
                );
                close(
                    sa.state.relaxation_slip_m()[0][j],
                    sb.state.relaxation_slip_m()[0][j],
                    5e-15,
                );
            }
            close(sa.dissipated_j, sb.dissipated_j, 4e-15);
            close(sb.closure_residual_j, 0.0, 5e-15);
            a = sa.state;
            b = sb.state;
        }
    }
}

#[test]
fn g1_noncommuting_and_rank_deficient_memory_keep_each_original_tensor_force_and_heat() {
    let law = model(1000.0);
    let x = input([0.002, -0.001], 0.001, 0.0);
    let step = law.advance(&law.initial_state(), &x).unwrap();
    let matrices = [
        [[400.0, 80.0], [80.0, 65.0]],
        [[64.0, -88.0], [-88.0, 121.0]],
    ];
    let times = [0.02, 0.07];
    let mut expected = product([[800.0, 70.0], [70.0, 500.0]], x.slip_increment_m).map(|v| 0.5 * v);
    let mut heat = 0.0;
    for j in 0..2 {
        let w = 1.0 / (1.0 + x.dt_s / (2.0 * times[j]));
        let mid = x.slip_increment_m.map(|v| 0.5 * w * v);
        let force = product(matrices[j], mid);
        heat += x.dt_s / times[j] * (mid[0] * force[0] + mid[1] * force[1]);
        for i in 0..2 {
            expected[i] += force[i];
            close(step.state.relaxation_slip_m()[j][i], 2.0 * mid[i], 1e-18);
        }
    }
    for i in 0..2 {
        close(step.traction_n[i], expected[i], 2e-15);
    }
    close(step.dissipated_j, heat, 1e-18);
    close(step.closure_residual_j, 0.0, 2e-18);
    let run = |n: usize| {
        let mut state = step.state.clone();
        let mut loss = 0.0;
        for _ in 0..n {
            let next = law
                .advance(&state, &input([0.0; 2], 0.02 / n as f64, 0.0))
                .unwrap();
            loss += next.dissipated_j;
            state = next.state;
        }
        close(
            law.stored_energy_j(&state) + loss,
            law.stored_energy_j(&step.state),
            2e-16,
        );
        let mut error = 0.0_f64;
        for j in 0..2 {
            for i in 0..2 {
                error = error.max(
                    (state.relaxation_slip_m()[j][i]
                        - step.state.relaxation_slip_m()[j][i] * (-0.02 / times[j]).exp())
                    .abs(),
                );
            }
        }
        error
    };
    let errors = [run(20), run(40), run(80)];
    for pair in errors.windows(2) {
        assert!((3.8..4.2).contains(&(pair[0] / pair[1])));
    }
}

#[test]
fn g1_tensor_plastic_work_is_independent_of_closure_residual_and_passive_under_reversals() {
    let law = model(0.1);
    let mut state = law.initial_state();
    let mut work = 0.0;
    let mut heat = 0.0;
    for i in 0..160 {
        let theta = 0.29 * i as f64;
        let x = input([0.0005 * theta.cos(), 0.0004 * theta.sin()], 0.0005, 0.0001);
        let step = law.advance(&state, &x).unwrap();
        let dp = [
            step.state.plastic_slip_m()[0] - state.plastic_slip_m()[0],
            step.state.plastic_slip_m()[1] - state.plastic_slip_m()[1],
        ];
        let channels = law.step_dissipation(&state, &step.state, x.dt_s).unwrap();
        close(
            channels.plastic_j,
            step.traction_n[0] * dp[0] + step.traction_n[1] * dp[1],
            2e-15,
        );
        assert!(channels.plastic_j >= 0.0 && channels.relaxation_j >= 0.0);
        work += step.traction_n[0] * step.slip_increment_m[0]
            + step.traction_n[1] * step.slip_increment_m[1];
        heat += step.dissipated_j;
        close(step.closure_residual_j, 0.0, 4e-15);
        state = step.state;
    }
    close(law.stored_energy_j(&state) + heat, work, 2e-13);
}

#[test]
fn g3_tensor_and_history_rotate_together_and_fixed_law_transport_refuses() {
    let law = model(0.1);
    let old = law
        .advance(&law.initial_state(), &input([0.003, 0.001], 0.001, 0.0))
        .unwrap();
    let angle = 0.63_f64;
    let a = TangentFrame::new([0.0, 0.0, 1.0], [1.0, 0.0, 0.0]).unwrap();
    let b = TangentFrame::new([0.0, 0.0, 1.0], [angle.cos(), angle.sin(), 0.0]).unwrap();
    let rotation = TangentialHistoryRotation::between(a, b, 0.0).unwrap();
    assert!(law.transport_state(&old.state, rotation).is_err());
    let transfer = law
        .transport_tensor_state(
            &old.state,
            rotation,
            TangentialTransportConfig {
                absolute_energy_tolerance_j: 1e-16,
                relative_energy_tolerance: 1e-13,
            },
        )
        .unwrap();
    close(transfer.stored_energy_residual_j, 0.0, 1e-17);
    let slip = [-0.0007, 0.0011];
    let x = law
        .advance(&old.state, &input(slip, 0.001, 0.0001))
        .unwrap();
    let y = transfer
        .law
        .advance(
            &transfer.state,
            &input(rotation.apply(slip).unwrap(), 0.001, 0.0001),
        )
        .unwrap();
    let force = rotation.apply(x.traction_n).unwrap();
    for i in 0..2 {
        close(force[i], y.traction_n[i], 4e-12);
    }
    close(x.dissipated_j, y.dissipated_j, 4e-15);
    close(x.stored_change_j, y.stored_change_j, 4e-15);
}

#[test]
fn g0_invalid_tensor_thermal_law_and_exhausted_return_refuse_without_mutation() {
    assert!(ShearStiffness::positive_definite([[1.0, 2.0], [2.0, 1.0]]).is_err());
    assert!(ShearStiffness::positive_definite([[1.0, 0.1], [0.2, 1.0]]).is_err());
    assert!(ShearStiffness::from_gram_rows(vec![[1.0, 0.0]; 7]).is_err());
    let law = model(0.1);
    let state = law.initial_state();
    let saved = state.clone();
    assert!(law.thermodynamics(&state, 0.0).is_err());
    assert!(
        law.advance(&state, &input([0.01, 0.02], 0.041, 0.0))
            .is_err()
    );
    let mut budget = input([0.01, 0.02], 0.001, 0.0);
    budget.solve.max_iterations = 1;
    assert!(law.advance(&state, &budget).is_err());
    assert_eq!(state, saved);
    let changed = ThirdBodyLayerContact::new_tensor(
        ShearStiffness::positive_definite([[801.0, 70.0], [70.0, 500.0]]).unwrap(),
        law.tensor_relaxation().unwrap().to_vec(),
        0.1,
        0.7,
        300.0,
        1.0,
    )
    .unwrap();
    assert!(changed.continue_state(&state).is_err());
    let mut arms = law.tensor_relaxation().unwrap().to_vec();
    arms[0].relaxation_time_s = 0.03;
    let kinetics = ThirdBodyLayerContact::new_tensor(
        ShearStiffness::positive_definite([[800.0, 70.0], [70.0, 500.0]]).unwrap(),
        arms,
        0.2,
        0.8,
        310.0,
        1.0,
    )
    .unwrap();
    assert!(kinetics.continue_state(&state).is_ok());
}
