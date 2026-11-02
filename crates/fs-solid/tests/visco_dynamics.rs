//! Exact constitutive work coupled to midpoint inertia on the spatial solid.
use fs_exec::Cx;
use fs_material::visco::tensor::{Mandel6, PronyBranch, TensorProny};
use fs_solid::linear3::visco::*;
use fs_solid::linear3::{TetAssemblyBudget, TetStaticSolveConfig};

const NODES: [[f64; 3]; 4] = [[0.0; 3], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
const CELLS: [[usize; 4]; 1] = [[0, 1, 2, 3]];
fn with_cx<R>(f: impl FnOnce(&Cx<'_>, &fs_exec::CancelGate) -> R) -> R {
    let gate = fs_exec::CancelGate::new_clock_free();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            fs_exec::StreamKey {
                seed: 611,
                kernel_id: 611,
                tile: 0,
                iteration: 0,
            },
            fs_exec::Budget::INFINITE,
            fs_exec::ExecMode::Deterministic,
        );
        f(&cx, &gate)
    })
}
fn diagonal(v: [f64; 6]) -> Mandel6 {
    core::array::from_fn(|i| core::array::from_fn(|j| if i == j { v[i] } else { 0.0 }))
}
fn outer(v: [f64; 6], scale: f64) -> Mandel6 {
    core::array::from_fn(|i| core::array::from_fn(|j| scale * v[i] * v[j]))
}
fn material(law: TensorProny) -> TetViscoMaterial {
    TetViscoMaterial {
        law,
        density_kg_m3: 1000.0,
        material_state_identity: [62; 32],
    }
}
fn system(law: &TetViscoMaterial, free: &[usize], cx: &Cx<'_>) -> ViscoElements {
    let fixed: Vec<_> = (0..12).filter(|d| !free.contains(d)).collect();
    ViscoElements::new(
        SolidMeshRef::Tet4 {
            nodes_m: &NODES,
            elements: &CELLS,
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
fn config() -> ViscoDynamicConfig {
    ViscoDynamicConfig {
        maximum_steps: 2000,
        maximum_free_dofs: 100,
        solve: TetStaticSolveConfig {
            maximum_iterations: 200,
            relative_residual_tolerance: 2e-12,
        },
        absolute_energy_tolerance_j: 1e-14,
        relative_energy_tolerance: 2e-11,
    }
}
fn close(a: f64, b: f64, rel: f64) {
    assert!(
        (a - b).abs() <= rel * a.abs().max(b.abs()).max(1e-20),
        "{a:e} != {b:e}"
    );
}
fn coupled() -> TetViscoMaterial {
    let mut eq = diagonal([90e3, 70e3, 50e3, 60e3, 30e3, 40e3]);
    eq[0][3] = 10e3;
    eq[3][0] = 10e3;
    eq[0][5] = -8e3;
    eq[5][0] = -8e3;
    material(
        TensorProny::new_mandel(
            eq,
            vec![
                PronyBranch::new(outer([2.0, 0.0, 1.0, 3.0, 0.0, -1.0], 1000.0), 0.003).unwrap(),
                PronyBranch::new(outer([1.0, -1.0, 0.0, -1.0, 2.0, 3.0], 500.0), 0.3).unwrap(),
            ],
            0.02,
        )
        .unwrap()
        .with_kelvin_viscosity_mandel(outer([1.0, 0.0, 0.0, -2.0, 0.0, 1.0], 100.0))
        .unwrap(),
    )
}

#[test]
fn mean_stress_and_tangent_are_the_original_ramp_work_pair() {
    let law = coupled().law;
    let mut warm = law.state();
    law.step(&mut warm, [0.001, -0.0002, 0.0, 0.0007, 0.0, -0.0004], 0.01)
        .unwrap();
    for dt in [1e-8, 0.003, 0.1, 10.0] {
        let tangent = law.mean_tangent_mandel(dt).unwrap();
        let target = [0.0003, 0.0001, 0.0, -0.0008, 0.0002, 0.0001];
        let mut state = warm.clone();
        let response = law.step_with_mean(&mut state, target, dt).unwrap();
        let (end, mean) = (response.stress, response.mean_stress);
        let work: f64 = mean
            .iter()
            .zip(target.iter().zip(warm.strain))
            .map(|(s, (a, b))| s * (a - b))
            .sum();
        close(state.work - warm.work, work, 3e-12);
        assert!(end.iter().zip(mean).any(|(a, b)| (a - b).abs() > 0.01));
        for j in 0..6 {
            let mut plus = warm.clone();
            let mut minus = warm.clone();
            let mut ep = target;
            let mut em = target;
            ep[j] += 1e-6;
            em[j] -= 1e-6;
            let sp = law.step_with_mean(&mut plus, ep, dt).unwrap().mean_stress;
            let sm = law.step_with_mean(&mut minus, em, dt).unwrap().mean_stress;
            for i in 0..6 {
                let error = ((sp[i] - sm[i]) / 2e-6 - tangent[i][j]).abs();
                assert!(
                    error <= 2e-9 * tangent[i][i].sqrt() * tangent[j][j].sqrt(),
                    "dt={dt} i={i} j={j}"
                );
            }
        }
    }
    // A large cumulative counter cannot erase this interval's physical heat.
    let mut ordinary = warm.clone();
    let mut accumulated = warm.clone();
    accumulated.work += 1e30;
    accumulated.dissipated += 1e30;
    let a = law.step_with_mean(&mut ordinary, [0.0; 6], 0.01).unwrap();
    let b = law
        .step_with_mean(&mut accumulated, [0.0; 6], 0.01)
        .unwrap();
    assert_eq!(a, b);
    assert!(b.step_dissipated > 0.0);
    assert_eq!(accumulated.dissipated, 1e30);
}

#[test]
fn elastic_midpoint_keeps_energy_without_numerical_heat() {
    with_cx(|cx, _| {
        let law = material(TensorProny::new_mandel(diagonal([1200.0; 6]), vec![], 0.02).unwrap());
        let solid = system(&law, &[3], cx);
        let prepared = solid.dynamic_operator(0.03, config(), cx).unwrap();
        let mut initial = solid.state(cx).unwrap();
        let mut u = [0.0; 12];
        u[3] = 0.003;
        initial.step_internal_force(&u, 0.01, cx).unwrap();
        let mut v = [0.0; 12];
        v[3] = 0.002;
        let mut state = prepared.start(initial, &v, cx).unwrap();
        let energy = state.kinetic_energy_j() + state.material_history().energy().stored_j;
        for _ in 0..1000 {
            let step = state.advance(&[0.0; 12], cx).unwrap();
            assert_eq!(step.constitutive.step_dissipated_j, 0.0);
            assert_eq!(step.support_work_j, 0.0);
            close(
                step.kinetic_energy_j + step.constitutive.energy.stored_j,
                energy,
                2e-11,
            );
        }
    });
}

#[test]
fn kelvin_ringdown_converges_to_the_damped_oscillator_and_independent_heat() {
    with_cx(|cx, _| {
        let e = 1200.0;
        let eta = 60.0;
        let u0 = 0.003;
        let v0 = -0.002;
        let end: f64 = 1.2;
        let law = material(
            TensorProny::new_mandel(diagonal([e; 6]), vec![], 0.02)
                .unwrap()
                .with_kelvin_viscosity_mandel(diagonal([eta, 0.0, 0.0, 0.0, 0.0, 0.0]))
                .unwrap(),
        );
        let solid = system(&law, &[3], cx);
        // Unit affine tet with node 1x free: m=rho/60, k=E/6, c=eta/6.
        let alpha: f64 = 5.0 * eta / law.density_kg_m3;
        let omega = (10.0 * e / law.density_kg_m3 - alpha * alpha).sqrt();
        let b = (v0 + alpha * u0) / omega;
        let oscillation = u0 * (omega * end).cos() + b * (omega * end).sin();
        let exact_u = (-alpha * end).exp() * oscillation;
        let exact_v = (-alpha * end).exp()
            * (-alpha * oscillation - u0 * omega * (omega * end).sin()
                + b * omega * (omega * end).cos());
        let mut errors = Vec::new();
        for count in [16, 32, 64] {
            let dt = end / f64::from(count);
            let prepared = solid.dynamic_operator(dt, config(), cx).unwrap();
            let mut initial = solid.state(cx).unwrap();
            let mut u = [0.0; 12];
            u[3] = u0;
            initial.step_internal_force(&u, 0.02, cx).unwrap();
            let prior_heat = initial.energy().dissipated_j;
            let mut v = [0.0; 12];
            v[3] = v0;
            let mut state = prepared.start(initial, &v, cx).unwrap();
            let initial_total =
                state.kinetic_energy_j() + state.material_history().energy().stored_j;
            let mut heat = 0.0;
            let mut previous = u0;
            for _ in 0..count {
                let step = state.advance(&[0.0; 12], cx).unwrap();
                let next = state.material_history().displacement_m()[3];
                heat += eta * (next - previous).powi(2) / (6.0 * dt);
                close(
                    step.constitutive.energy.dissipated_j - prior_heat,
                    heat,
                    3e-11,
                );
                close(
                    step.kinetic_energy_j + step.constitutive.energy.stored_j + heat,
                    initial_total,
                    3e-11,
                );
                previous = next;
            }
            errors.push(
                (state.material_history().displacement_m()[3] - exact_u)
                    .hypot((state.velocity_mps()[3] - exact_v) / omega),
            );
        }
        assert!(errors[1] < 0.28 * errors[0], "{errors:?}");
        assert!(errors[2] < 0.28 * errors[1], "{errors:?}");
        assert!(
            errors[2] < 0.002 * exact_u.hypot(exact_v / omega),
            "{errors:?}"
        );
    });
}

#[test]
fn coupled_tensor_memory_closes_load_work_and_passive_ringdown() {
    with_cx(|cx, _| {
        let solid = system(&coupled(), &[3, 4, 5], cx);
        let prepared = solid.dynamic_operator(0.002, config(), cx).unwrap();
        assert!(prepared.effective_mass().get(0, 1).abs() > 1e-4);
        let mut initial = solid.state(cx).unwrap();
        let mut u = [0.0; 12];
        u[0] = 0.0002;
        u[3] = 0.001;
        u[4] = -0.0004;
        u[5] = 0.0006;
        initial.step_internal_force(&u, 0.001, cx).unwrap();
        let prior_heat = initial.energy().dissipated_j;
        let mut state = prepared.start(initial, &[0.0; 12], cx).unwrap();
        let initial_total = state.kinetic_energy_j() + state.material_history().energy().stored_j;
        let mut work = 0.0;
        let mut previous = initial_total;
        for i in 0..240 {
            let mut load = [0.0; 12];
            if i < 40 {
                load[3] = 3.0;
                load[4] = -1.0;
                load[5] = 0.5;
            }
            let before = state.material_history().displacement_m().to_vec();
            let step = state.advance(&load, cx).unwrap();
            let after = state.material_history().displacement_m();
            let independent: f64 = load
                .iter()
                .zip(after.iter().zip(before))
                .map(|(f, (a, b))| f * (a - b))
                .sum();
            close(step.external_work_j, independent, 1e-12);
            assert_eq!(after[0], u[0]);
            work += independent;
            let stored = step.kinetic_energy_j + step.constitutive.energy.stored_j;
            let heat = step.constitutive.energy.dissipated_j - prior_heat;
            close(initial_total + work, stored + heat, 3e-10);
            if i >= 40 {
                assert!(stored <= previous + 1e-13);
            }
            assert!(step.constitutive.step_dissipated_j >= 0.0);
            previous = stored;
        }
    });
}

#[test]
fn free_body_consistent_mass_translates_without_artificial_clamps() {
    with_cx(|cx, _| {
        let law = coupled();
        let solid = system(&law, &(0..12).collect::<Vec<_>>(), cx);
        let mut controls = config();
        controls.solve.relative_residual_tolerance = 1e-13;
        let prepared = solid.dynamic_operator(0.01, controls, cx).unwrap();
        let speed = [0.1, -0.03, 0.02];
        let acceleration = [0.04, -0.01, 0.02];
        let v: Vec<_> = (0..4).flat_map(|_| speed).collect();
        let load: Vec<_> = (0..4)
            .flat_map(|_| acceleration.map(|a| law.density_kg_m3 * a / 24.0))
            .collect();
        let mut state = prepared.start(solid.state(cx).unwrap(), &v, cx).unwrap();
        for _ in 0..20 {
            state.advance(&load, cx).unwrap();
        }
        let time = state.time_s();
        for node in 0..4 {
            for k in 0..3 {
                close(
                    state.velocity_mps()[3 * node + k],
                    speed[k] + acceleration[k] * time,
                    1e-10,
                );
                close(
                    state.material_history().displacement_m()[3 * node + k],
                    speed[k] * time + 0.5 * acceleration[k] * time * time,
                    1e-10,
                );
            }
        }
        assert!(state.material_history().energy().dissipated_j < 1e-18);
    });
}

#[test]
fn dynamic_refusals_keep_motion_material_memory_and_clock_atomic() {
    with_cx(|cx, gate| {
        let solid = system(&coupled(), &[3, 4, 5], cx);
        let prepared = solid.dynamic_operator(0.004, config(), cx).unwrap();
        let limited = solid
            .dynamic_operator(
                0.004,
                ViscoDynamicConfig {
                    maximum_steps: 1,
                    ..config()
                },
                cx,
            )
            .unwrap();
        let mut short = limited
            .start(solid.state(cx).unwrap(), &[0.0; 12], cx)
            .unwrap();
        short.advance(&[0.0; 12], cx).unwrap();
        assert!(short.advance(&[0.0; 12], cx).is_err());
        assert_eq!(short.steps(), 1);
        assert_eq!(short.time_s(), 0.004);
        let mut moving = [0.0; 12];
        moving[0] = 0.1;
        assert!(
            prepared
                .start(solid.state(cx).unwrap(), &moving, cx)
                .is_err()
        );
        let mut state = prepared
            .start(solid.state(cx).unwrap(), &[0.0; 12], cx)
            .unwrap();
        let mut load = [0.0; 12];
        load[3] = 2.0;
        load[4] = -1.0;
        state.advance(&load, cx).unwrap();
        let mut clean = state.clone();
        let mut excessive = load;
        excessive[3] = 1e12;
        assert!(state.advance(&excessive, cx).is_err());
        assert_eq!(
            state.material_history().energy(),
            clean.material_history().energy()
        );
        assert_eq!(
            state.material_history().displacement_m(),
            clean.material_history().displacement_m()
        );
        assert_eq!(state.velocity_mps(), clean.velocity_mps());
        assert_eq!(state.time_s(), clean.time_s());
        gate.request();
        assert!(state.advance(&load, cx).is_err());
        with_cx(|retry, _| {
            state.advance(&load, retry).unwrap();
            clean.advance(&load, retry).unwrap();
            assert_eq!(
                state.material_history().energy(),
                clean.material_history().energy()
            );
            assert_eq!(
                state.material_history().displacement_m(),
                clean.material_history().displacement_m()
            );
            assert_eq!(state.velocity_mps(), clean.velocity_mps());
            assert_eq!(state.time_s(), clean.time_s());
        });
    });
}
