//! Load-controlled creep on the original spatial tensor-Prony owner.
use fs_exec::Cx;
use fs_material::visco::tensor::{Mandel6, PronyBranch, TensorProny};
use fs_solid::linear3::visco::*;
use fs_solid::linear3::{TetAssemblyBudget, TetElasticError, TetStaticSolveConfig};

const NODES: [[f64; 3]; 4] = [
    [0.0, 0.0, 0.0],
    [1.0, 0.0, 0.0],
    [0.0, 1.0, 0.0],
    [0.0, 0.0, 1.0],
];
const CELL: [[usize; 4]; 1] = [[0, 1, 2, 3]];

fn with_cx<R>(f: impl FnOnce(&Cx<'_>, &fs_exec::CancelGate) -> R) -> R {
    let gate = fs_exec::CancelGate::new_clock_free();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            fs_exec::StreamKey {
                seed: 610,
                kernel_id: 610,
                tile: 0,
                iteration: 0,
            },
            fs_exec::Budget::INFINITE,
            fs_exec::ExecMode::Deterministic,
        );
        f(&cx, &gate)
    })
}

fn config() -> TetStaticSolveConfig {
    TetStaticSolveConfig {
        relative_residual_tolerance: 2e-12,
        maximum_iterations: 100,
    }
}

fn diagonal(value: [f64; 6]) -> Mandel6 {
    core::array::from_fn(|i| core::array::from_fn(|j| if i == j { value[i] } else { 0.0 }))
}

fn rank_one(row: [f64; 6], scale: f64) -> Mandel6 {
    core::array::from_fn(|i| core::array::from_fn(|j| scale * row[i] * row[j]))
}

fn material(law: TensorProny) -> TetViscoMaterial {
    TetViscoMaterial {
        law,
        density_kg_m3: 1000.0,
        material_state_identity: [61; 32],
    }
}

fn system(law: &TetViscoMaterial, free: &[usize], cx: &Cx<'_>) -> ViscoElements {
    let fixed: Vec<_> = (0..12).filter(|d| !free.contains(d)).collect();
    ViscoElements::new(
        SolidMeshRef::Tet4 {
            nodes_m: &NODES,
            elements: &CELL,
        },
        TetViscoField::Uniform(law),
        &fixed,
        TetAssemblyBudget::standard(),
        ViscoStateBudget {
            max_internal_variables: 1000,
        },
        cx,
    )
    .unwrap()
}

fn close(a: f64, b: f64, relative: f64) {
    assert!(
        (a - b).abs() <= relative * a.abs().max(b.abs()).max(1e-20),
        "{a:e} != {b:e}"
    );
}

fn nonproportional() -> TetViscoMaterial {
    let mut eq = diagonal([90e3, 70e3, 50e3, 60e3, 30e3, 40e3]);
    eq[0][3] = 10e3;
    eq[3][0] = 10e3;
    eq[0][5] = -8e3;
    eq[5][0] = -8e3;
    material(
        TensorProny::new_mandel(
            eq,
            vec![
                PronyBranch::new(rank_one([2.0, 0.0, 1.0, 3.0, 0.0, -1.0], 1000.0), 0.1).unwrap(),
                PronyBranch::new(rank_one([1.0, -1.0, 0.0, -1.0, 2.0, 3.0], 500.0), 0.7).unwrap(),
            ],
            0.02,
        )
        .unwrap()
        .with_kelvin_viscosity_mandel(rank_one([1.0, 0.0, 0.0, -2.0, 0.0, 1.0], 100.0))
        .unwrap(),
    )
}

#[test]
fn standard_linear_solid_creep_converges_to_analytic_retardation() {
    with_cx(|cx, _| {
        let e0 = 10e3;
        let e1 = 30e3;
        let tau = 0.15;
        let law = material(
            TensorProny::new(
                diagonal([e0; 6]),
                None,
                vec![PronyBranch::new(diagonal([e1, 0.0, 0.0, 0.0, 0.0, 0.0]), tau).unwrap()],
                (0.0, 1e8),
            )
            .unwrap(),
        );
        let solid = system(&law, &[3], cx);
        let mut initial = solid.state(cx).unwrap();
        let mut displacement = [0.0; 12];
        displacement[3] = 0.001;
        let preload = initial
            .step_internal_force(&displacement, 1e-5, cx)
            .unwrap();
        let mut force = [0.0; 12];
        force[3] = preload.full_force_n[3];
        // Unit affine tetrahedron: volume=1/6 and epsilon_xx=u_1x.
        let stress = 6.0 * force[3];
        let relaxed = stress / e0;
        let retardation = tau * (e0 + e1) / e0;
        let end_time = 0.9;
        let exact = relaxed + (displacement[3] - relaxed) * (-end_time / retardation).exp();
        let mut errors = Vec::new();
        for count in [8, 16, 32] {
            let prepared = solid
                .quasistatic_operator(end_time / f64::from(count), cx)
                .unwrap();
            let mut state = initial.clone();
            let mut previous = displacement[3];
            for _ in 0..count {
                let step = state
                    .step_quasistatic(&prepared, &force, &[0.0; 12], config(), cx)
                    .unwrap();
                assert!(step.displacement_m[3] > previous);
                assert!(step.constitutive.step_dissipated_j > 0.0);
                assert!(step.true_relative_residual < 2e-12);
                let ledger = step.constitutive.energy;
                close(ledger.work_j, ledger.stored_j + ledger.dissipated_j, 3e-12);
                previous = step.displacement_m[3];
            }
            errors.push((state.displacement_m()[3] - exact).abs());
        }
        assert!(errors[1] < 0.3 * errors[0], "{errors:?}");
        assert!(errors[2] < 0.3 * errors[1], "{errors:?}");
        assert!(errors[2] < 2e-4 * exact, "{errors:?}");
    });
}

#[test]
fn kelvin_creep_reaches_analytic_limit_without_invented_dissipation() {
    with_cx(|cx, _| {
        let e = 20e3;
        let eta = 4000.0;
        let law = material(
            TensorProny::new(
                diagonal([e; 6]),
                Some(diagonal([eta, 0.0, 0.0, 0.0, 0.0, 0.0])),
                vec![],
                (0.0, 1e8),
            )
            .unwrap(),
        );
        let solid = system(&law, &[3], cx);
        let mut load = [0.0; 12];
        load[3] = 4.0;
        let end_time = 0.4;
        let exact = 6.0 * load[3] / e * (1.0 - (-e * end_time / eta).exp());
        let mut errors = Vec::new();
        for count in [32, 64] {
            let dt = end_time / f64::from(count);
            let prepared = solid.quasistatic_operator(dt, cx).unwrap();
            let mut state = solid.state(cx).unwrap();
            let mut previous = 0.0;
            let mut independent_heat = 0.0;
            for _ in 0..count {
                let step = state
                    .step_quasistatic(&prepared, &load, &[0.0; 12], config(), cx)
                    .unwrap();
                let delta = step.displacement_m[3] - previous;
                independent_heat += eta * delta * delta / (6.0 * dt);
                close(
                    step.constitutive.energy.dissipated_j,
                    independent_heat,
                    2e-12,
                );
                close(
                    step.constitutive.energy.work_j,
                    step.constitutive.energy.stored_j + independent_heat,
                    2e-12,
                );
                previous = step.displacement_m[3];
            }
            errors.push((state.displacement_m()[3] - exact).abs());
        }
        assert!(errors[1] < 0.52 * errors[0], "{errors:?}");
        assert!(errors[1] < 6e-6, "{errors:?}");
    });
}

#[test]
fn nonproportional_memory_and_moving_supports_recover_original_history() {
    with_cx(|cx, _| {
        let solid = system(&nonproportional(), &[3, 4, 5], cx);
        let prepared = solid.quasistatic_operator(0.035, cx).unwrap();
        assert!(prepared.stiffness().get(0, 1).abs() > 100.0);
        assert!(prepared.stiffness().get(0, 2).abs() > 100.0);
        let mut reference = solid.state(cx).unwrap();
        let mut solved = solid.state(cx).unwrap();
        for sign in [1.0, -0.6, 0.3, 0.3, -0.9] {
            let mut expected = [0.0; 12];
            expected[0] = sign * 0.0003;
            expected[3] = sign * 0.001;
            expected[4] = -sign * 0.0004;
            expected[5] = sign * 0.0007;
            let truth = reference
                .step_internal_force(&expected, prepared.duration_s(), cx)
                .unwrap();
            let mut load = [0.0; 12];
            for d in [3, 4, 5] {
                load[d] = truth.full_force_n[d];
            }
            let mut support = [0.0; 12];
            support[0] = expected[0];
            let step = solved
                .step_quasistatic(&prepared, &load, &support, config(), cx)
                .unwrap();
            for d in 0..12 {
                close(step.displacement_m[d], expected[d], 2e-11);
                if d != 3 && d != 4 && d != 5 {
                    close(step.reaction_force_n[d], truth.full_force_n[d], 2e-11);
                }
            }
            close(step.constitutive.energy.work_j, truth.energy.work_j, 2e-11);
            close(
                step.constitutive.energy.stored_j,
                truth.energy.stored_j,
                2e-11,
            );
            close(
                step.constitutive.energy.dissipated_j,
                truth.energy.dissipated_j,
                2e-11,
            );
        }
        // Independent central differences of the original warm constitutive
        // history recover every assembled column, including shear cross terms.
        let h = 1e-6;
        for j in 0..3 {
            let mut plus = solved.clone();
            let mut minus = solved.clone();
            let mut up = solved.displacement_m().to_vec();
            let mut um = up.clone();
            up[j + 3] += h;
            um[j + 3] -= h;
            let fp = plus
                .step_internal_force(&up, prepared.duration_s(), cx)
                .unwrap();
            let fm = minus
                .step_internal_force(&um, prepared.duration_s(), cx)
                .unwrap();
            for i in 0..3 {
                close(
                    (fp.reduced_force_n[i] - fm.reduced_force_n[i]) / (2.0 * h),
                    prepared.stiffness().get(i, j),
                    2e-9,
                );
            }
        }
    });
}

#[test]
fn failed_load_solve_preserves_warm_memory_and_exact_retry() {
    with_cx(|cx, _| {
        let solid = system(&nonproportional(), &[3, 4, 5], cx);
        let prepared = solid.quasistatic_operator(0.02, cx).unwrap();
        let mut state = solid.state(cx).unwrap();
        let mut load = [0.0; 12];
        load[3] = 5.0;
        load[4] = -2.0;
        load[5] = 3.0;
        state
            .step_quasistatic(&prepared, &load, &[0.0; 12], config(), cx)
            .unwrap();
        let mut retry = state.clone();
        load[3] = -3.0;
        let before_energy = state.energy();
        let before_displacement = state.displacement_m().to_vec();
        assert!(matches!(
            state.step_quasistatic(
                &prepared,
                &load,
                &[0.0; 12],
                TetStaticSolveConfig {
                    maximum_iterations: 1,
                    ..config()
                },
                cx
            ),
            Err(ViscoElementError::Elastic(
                TetElasticError::StaticSolveFailed { .. }
            ))
        ));
        assert_eq!(state.energy(), before_energy);
        assert_eq!(state.displacement_m(), before_displacement);
        let mut excessive = load;
        excessive[3] = 1e12;
        assert!(
            state
                .step_quasistatic(&prepared, &excessive, &[0.0; 12], config(), cx)
                .is_err()
        );
        assert_eq!(state.energy(), before_energy);
        assert_eq!(state.displacement_m(), before_displacement);
        let a = state
            .step_quasistatic(&prepared, &load, &[0.0; 12], config(), cx)
            .unwrap();
        let b = retry
            .step_quasistatic(&prepared, &load, &[0.0; 12], config(), cx)
            .unwrap();
        assert_eq!(a.displacement_m, b.displacement_m);
        assert_eq!(a.constitutive.energy, b.constitutive.energy);
    });
}

#[test]
fn free_bodies_refuse_load_control_but_keep_prescribed_histories() {
    with_cx(|cx, _| {
        let solid = system(&nonproportional(), &(0..12).collect::<Vec<_>>(), cx);
        assert!(matches!(
            solid.quasistatic_operator(0.01, cx),
            Err(ViscoElementError::Elastic(
                TetElasticError::UnconstrainedRigidModes { .. }
            ))
        ));
        assert!(
            solid
                .state(cx)
                .unwrap()
                .step_internal_force(&[0.0; 12], 0.01, cx)
                .is_ok()
        );
    });
}

#[test]
fn cancellation_and_foreign_operator_leave_the_original_history_retryable() {
    with_cx(|cx, gate| {
        let law = nonproportional();
        let solid = system(&law, &[3, 4, 5], cx);
        let other = system(&law, &[3, 4, 5], cx);
        let prepared = solid.quasistatic_operator(0.04, cx).unwrap();
        let foreign = other.quasistatic_operator(0.04, cx).unwrap();
        let mut state = solid.state(cx).unwrap();
        let mut load = [0.0; 12];
        load[3] = 4.0;
        load[4] = -1.0;
        state
            .step_quasistatic(&prepared, &load, &[0.0; 12], config(), cx)
            .unwrap();
        let mut checkpointed = state.clone();
        assert!(matches!(
            state.step_quasistatic(&foreign, &load, &[0.0; 12], config(), cx),
            Err(ViscoElementError::InvalidHistory { .. })
        ));
        assert_eq!(state.energy(), checkpointed.energy());
        assert_eq!(state.displacement_m(), checkpointed.displacement_m());
        gate.request();
        assert!(matches!(
            state.step_quasistatic(&prepared, &load, &[0.0; 12], config(), cx),
            Err(ViscoElementError::Elastic(TetElasticError::Cancelled))
        ));
        assert_eq!(state.energy(), checkpointed.energy());
        assert_eq!(state.displacement_m(), checkpointed.displacement_m());
        with_cx(|retry_cx, _| {
            let a = state
                .step_quasistatic(&prepared, &load, &[0.0; 12], config(), retry_cx)
                .unwrap();
            let b = checkpointed
                .step_quasistatic(&prepared, &load, &[0.0; 12], config(), retry_cx)
                .unwrap();
            assert_eq!(a.displacement_m, b.displacement_m);
            assert_eq!(a.constitutive.energy, b.constitutive.energy);
        });
    });
}

#[test]
fn narrow_band_material_cannot_be_extrapolated_to_sustained_load_creep() {
    with_cx(|cx, _| {
        let band = (10.0, 10_000.0);
        let law = material(TensorProny::new(diagonal([20e3; 6]), None, vec![], band).unwrap());
        let solid = system(&law, &[3], cx);
        assert!(matches!(solid.quasistatic_operator(0.01, cx),
            Err(ViscoElementError::Material(fs_material::visco::ViscoError::OutOfBand {
                omega: 0.0,
                band: refused,
            })) if refused == band));
    });
}
