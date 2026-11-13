//! Actual spatial/material constructors and original MIMO histories.
use super::*;
use fs_material::visco::tensor::{Mandel6, PronyBranch, TensorProny};
use fs_math::c64::C64;
use fs_solid::linear3::visco::*;
use fs_solid::linear3::{TetAssemblyBudget, TetStaticSolveConfig};
use fs_vfit::{
    mimo::{MimoModel, load::MimoLoadLimits},
    model::PoleTerm,
    passivity::CertificateClass,
};

fn with_cx(f: impl FnOnce(&Cx<'_>, &fs_exec::CancelGate)) {
    let gate = fs_exec::CancelGate::new_clock_free();
    fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default()).scope(|arena| {
        f(
            &Cx::new(
                &gate,
                arena,
                fs_exec::StreamKey {
                    seed: 491,
                    kernel_id: 96,
                    tile: 0,
                    iteration: 0,
                },
                fs_exec::Budget::INFINITE,
                fs_exec::ExecMode::Deterministic,
            ),
            &gate,
        )
    });
}

fn diagonal(value: f64) -> Mandel6 {
    core::array::from_fn(|i| core::array::from_fn(|j| if i == j { value } else { 0.0 }))
}

// Unit tetrahedron, only node 1x free: M=1, K=2.5; optional original
// Kelvin C=0.2 and all six Maxwell strains with K=6.5, tau=0.25.
fn solid(memory: bool, cx: &Cx<'_>) -> ViscoElements {
    let branches = if memory {
        vec![PronyBranch::new(diagonal(39.0), 0.25).unwrap()]
    } else {
        vec![]
    };
    let mut law = TensorProny::new_mandel(diagonal(15.0), branches, 0.05).unwrap();
    if memory {
        law = law.with_kelvin_viscosity_mandel(diagonal(1.2)).unwrap();
    }
    let material = TetViscoMaterial {
        law,
        density_kg_m3: 60.0,
        material_state_identity: [96; 32],
    };
    let fixed: Vec<_> = (0..12).filter(|&i| i != 3).collect();
    ViscoElements::new(
        SolidMeshRef::Tet4 {
            nodes_m: &[[0.0; 3], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
            elements: &[[0, 1, 2, 3]],
        },
        TetViscoField::Uniform(&material),
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
        maximum_steps: 1000,
        maximum_free_dofs: 100,
        solve: TetStaticSolveConfig {
            maximum_iterations: 100,
            relative_residual_tolerance: 2e-12,
        },
        absolute_energy_tolerance_j: 1e-15,
        relative_energy_tolerance: 2e-11,
    }
}

fn budget() -> SolidMimoTransientBudget {
    SolidMimoTransientBudget {
        ports: ViscoDynamicPortBudget {
            maximum_ports: 8,
            maximum_payload_bytes: 1 << 20,
        },
        maximum_load_payload_bytes: 3 << 20,
        maximum_relative_exchange_error: 2e-11,
    }
}

fn load(model: &MimoModel, dt: f64) -> PreparedMimoLoad {
    PreparedMimoLoad::new(
        model,
        dt,
        MimoLoadLimits {
            band_rad_s: [0.01, core::f64::consts::PI / dt],
            sampled_passivity_points: 33,
            max_states: 32,
            max_memory_bytes: 1 << 20,
        },
    )
    .unwrap()
}

fn mutual_model() -> MimoModel {
    MimoModel {
        ports: 2,
        poles: vec![
            PoleTerm::Real {
                pole: -2.0,
                residue: 0.0,
            },
            PoleTerm::Pair {
                pole: C64::new(-0.4, 1.8),
                residue: C64::ZERO,
            },
        ],
        residues: vec![
            [0.8, -0.25, -0.25, 0.6].map(C64::from_re).to_vec(),
            [0.4, 0.12, 0.12, 0.3].map(C64::from_re).to_vec(),
        ],
        d: vec![0.9, 0.15, 0.15, 0.7],
        e: vec![0.0; 4],
        reciprocal: true,
    }
}

fn close(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance * a.abs().max(b.abs()).max(1e-16),
        "{a:e} != {b:e}"
    );
}

#[test]
fn g1_complete_mutual_memory_matches_independent_scalar_equations_and_material_heat() {
    with_cx(|cx, _| {
        let body = solid(true, cx);
        assert_eq!(body.internal_variables(), 6);
        let dt = 0.02;
        let operator = body.dynamic_operator(dt, config(), cx).unwrap();
        let mut material = body.state(cx).unwrap();
        let (mut u, mut v, fixed_u, warm) = (0.003, 0.004, 0.001, 0.03);
        let mut displacement = [0.0; 12];
        displacement[3] = u;
        displacement[7] = fixed_u;
        material
            .step_internal_force(&displacement, warm, cx)
            .unwrap();
        let mut velocity = [0.0; 12];
        velocity[3] = v;
        let trajectory = operator.start(material, &velocity, cx).unwrap();
        let mut original = load(&mutual_model(), dt);
        for flow in [[0.003, -0.002], [-0.001, 0.004], [0.002, 0.001]] {
            original.step_mean_into(&flow, &mut [0.0; 2]).unwrap();
        }
        let mut reference = original.clone();
        let initial_poles = original.pole_states().to_vec();
        let map = [1.0, -0.4];
        let mut coupled =
            PreparedSolidMimoTransient::new(trajectory, original, &map, budget(), cx).unwrap();
        assert_eq!(coupled.source().pole_states(), initial_poles);
        assert_eq!(coupled.source_reports().1.class, CertificateClass::GridOnly);
        let tau = 0.25;
        let warm_gain = tau * (1.0 - (-warm / tau).exp()) / warm;
        let (mut spring, mut fixed_spring) = (warm_gain * u, warm_gain * fixed_u);
        let decay = (-dt / tau).exp();
        let mean_gain = tau * (1.0 - decay) / dt;
        let mut total_source_work = 0.0;
        let mut total_heat = 0.0;
        let mut total_applied = 0.0;
        let initial_energy = coupled.trajectory().kinetic_energy_j()
            + coupled.trajectory().material_history().energy().stored_j;
        for step in 0..80 {
            let applied = if step < 25 { 0.0008 } else { 0.0 };
            let mut history = [0.0; 2];
            reference.mean_history_output_into(&mut history).unwrap();
            let d = reference.mean_direct_term();
            let drag =
                map[0] * (d[0] * map[0] + d[1] * map[1]) + map[1] * (d[2] * map[0] + d[3] * map[1]);
            let predicted_mean = (2.0 * v / dt + applied
                - 2.5 * u
                - 6.5 * mean_gain * spring
                - map[0] * history[0]
                - map[1] * history[1])
                / (2.0 / dt + 0.5 * dt * 2.5 + 0.2 + 6.5 * tau * (1.0 - mean_gain) + drag);
            let mut force = [0.0; 12];
            force[3] = applied;
            let report = coupled.advance(&force, cx).unwrap();
            let actual_u = coupled.trajectory().material_history().displacement_m()[3];
            let actual_mean = (actual_u - u) / dt;
            close(actual_mean, predicted_mean, 8e-10);
            let expected_u = u + dt * predicted_mean;
            close(actual_u, expected_u, 2e-10);
            close(
                coupled.trajectory().velocity_mps()[3],
                2.0 * predicted_mean - v,
                8e-10,
            );
            close(report.structure.port.mean_flows[0], actual_mean, 2e-12);
            close(
                report.structure.port.mean_flows[1],
                -0.4 * actual_mean,
                2e-12,
            );
            let steady = tau * actual_mean;
            let transient = spring - steady;
            let integral = steady * steady * dt
                + 2.0 * steady * transient * tau * (1.0 - decay)
                + transient * transient * 0.5 * tau * (1.0 - decay * decay);
            let expected_heat = 0.2 * actual_mean * actual_mean * dt
                + 6.5 * integral / tau
                + 0.5 * 6.5 * fixed_spring * fixed_spring * (1.0 - decay * decay);
            close(
                report.structure.body.constitutive.step_dissipated_j,
                expected_heat,
                2e-10,
            );
            spring = decay * spring + tau * (1.0 - decay) * actual_mean;
            fixed_spring *= decay;
            let expected_storage = 0.5 * 2.5 * (actual_u * actual_u + fixed_u * fixed_u)
                + 0.5 * 6.5 * (spring * spring + fixed_spring * fixed_spring);
            close(
                report.structure.body.constitutive.energy.stored_j,
                expected_storage,
                2e-11,
            );
            let mut mean_effort = [0.0; 2];
            reference
                .step_mean_into(&report.structure.port.mean_flows, &mut mean_effort)
                .unwrap();
            assert_eq!(coupled.source().pole_states(), reference.pole_states());
            for (actual, expected) in report.structure.port.mean_efforts.iter().zip(mean_effort) {
                close(*actual, expected, 2e-11);
            }
            total_source_work += report.source_work_j;
            total_heat += report.structure.body.constitutive.step_dissipated_j;
            total_applied += report.structure.port.applied_work_j;
            close(
                coupled.trajectory().kinetic_energy_j()
                    + expected_storage
                    + total_heat
                    + total_source_work,
                initial_energy + total_applied,
                2e-10,
            );
            u = actual_u;
            v = coupled.trajectory().velocity_mps()[3];
        }
        assert!(total_heat > 0.0);
        assert_eq!(coupled.trajectory().steps(), 80);
        assert_eq!(coupled.source().source_model(), &mutual_model());
    });
}

#[test]
fn g3_three_clocks_converge_to_complete_continuous_solid_and_original_load_memory() {
    with_cx(|cx, _| {
        let body = solid(false, cx);
        let model = MimoModel {
            ports: 1,
            poles: vec![PoleTerm::Real {
                pole: -2.0,
                residue: 0.0,
            }],
            residues: vec![vec![C64::from_re(0.8)]],
            d: vec![0.3],
            e: vec![0.0],
            reciprocal: true,
        };
        let initial = [0.003, 0.001, 0.0004];
        let duration = 0.8;
        // Independent analytic exp(A*t)*initial by its rapidly convergent
        // vector Taylor series; no tested propagator or port Schur is used.
        let a = [[0.0, 1.0, 0.0], [-2.5, -0.3, -0.8], [0.0, 1.0, -2.0]];
        let (mut exact, mut term) = (initial, initial);
        for order in 1..=80 {
            term = core::array::from_fn(|i| {
                (0..3).map(|j| a[i][j] * term[j]).sum::<f64>() * duration / f64::from(order)
            });
            for i in 0..3 {
                exact[i] += term[i];
            }
        }
        let mut errors = Vec::new();
        for intervals in [20, 40, 80] {
            let dt = duration / f64::from(intervals);
            let operator = body.dynamic_operator(dt, config(), cx).unwrap();
            let mut material = body.state(cx).unwrap();
            let mut u = [0.0; 12];
            u[3] = initial[0];
            material.step_internal_force(&u, 0.01, cx).unwrap();
            let mut v = [0.0; 12];
            v[3] = initial[1];
            let trajectory = operator.start(material, &v, cx).unwrap();
            let original = load(&model, dt)
                .continuous_operator()
                .candidate_from_coordinates(&initial[2..], &initial[1..2])
                .unwrap();
            let mut coupled =
                PreparedSolidMimoTransient::new(trajectory, original, &[1.0], budget(), cx)
                    .unwrap();
            for _ in 0..intervals {
                coupled.advance(&[0.0; 12], cx).unwrap();
            }
            let actual = [
                coupled.trajectory().material_history().displacement_m()[3],
                coupled.trajectory().velocity_mps()[3],
                coupled.source().pole_states()[0].re,
            ];
            errors.push(
                actual
                    .iter()
                    .zip(exact)
                    .map(|(a, b)| (a - b).powi(2))
                    .sum::<f64>()
                    .sqrt(),
            );
        }
        assert!(
            errors.windows(2).all(|pair| pair[0] / pair[1] > 3.8),
            "{errors:?}"
        );
        assert!(errors[2] < 2e-7, "{errors:?}");
    });
}

#[test]
fn g1_singular_and_zero_direct_terms_keep_signed_cancelling_ports() {
    with_cx(|cx, _| {
        let body = solid(false, cx);
        let operator = body.dynamic_operator(0.03, config(), cx).unwrap();
        for d in [vec![0.0; 4], vec![0.5; 4]] {
            let mut material = body.state(cx).unwrap();
            let mut displacement = [0.0; 12];
            displacement[3] = 0.003;
            material
                .step_internal_force(&displacement, 0.01, cx)
                .unwrap();
            let state = operator.start(material, &[0.0; 12], cx).unwrap();
            let mut control = state.clone();
            let model = MimoModel {
                ports: 2,
                poles: vec![],
                residues: vec![],
                d,
                e: vec![0.0; 4],
                reciprocal: true,
            };
            let mut coupled = PreparedSolidMimoTransient::new(
                state,
                load(&model, 0.03),
                &[1.0, -1.0],
                budget(),
                cx,
            )
            .unwrap();
            for _ in 0..12 {
                control.advance(&[0.0; 12], cx).unwrap();
                let step = coupled.advance(&[0.0; 12], cx).unwrap();
                close(
                    coupled.trajectory().material_history().displacement_m()[3],
                    control.material_history().displacement_m()[3],
                    2e-12,
                );
                assert_eq!(step.source_work_j, 0.0);
            }
        }
    });
}

#[test]
fn g4_clock_budget_and_refused_intervals_preserve_both_accepted_owners() {
    with_cx(|cx, gate| {
        let body = solid(true, cx);
        let operator = body.dynamic_operator(0.02, config(), cx).unwrap();
        let state = operator
            .start(body.state(cx).unwrap(), &[0.0; 12], cx)
            .unwrap();
        let model = mutual_model();
        let wrong_clock = f64::from_bits(0.02_f64.to_bits() + 1);
        assert!(matches!(
            PreparedSolidMimoTransient::new(
                state.clone(),
                load(&model, wrong_clock),
                &[1.0, -0.4],
                budget(),
                cx
            ),
            Err(SolidMimoTransientError::Invalid(_))
        ));
        assert!(
            PreparedSolidMimoTransient::new(
                state.clone(),
                load(&model, 0.02),
                &[1.0, -0.4],
                SolidMimoTransientBudget {
                    maximum_load_payload_bytes: 1,
                    ..budget()
                },
                cx
            )
            .is_err()
        );
        let mut coupled =
            PreparedSolidMimoTransient::new(state, load(&model, 0.02), &[1.0, -0.4], budget(), cx)
                .unwrap();
        let mut force = [0.0; 12];
        force[3] = 0.001;
        coupled.advance(&force, cx).unwrap();
        let accepted_u = coupled
            .trajectory()
            .material_history()
            .displacement_m()
            .to_vec();
        let accepted_v = coupled.trajectory().velocity_mps().to_vec();
        let accepted_energy = coupled.trajectory().material_history().energy();
        let accepted_poles = coupled.source().pole_states().to_vec();
        for rejected in [f64::NAN, 1e12] {
            force[3] = rejected;
            assert!(coupled.advance(&force, cx).is_err());
            assert_eq!(
                coupled.trajectory().material_history().displacement_m(),
                accepted_u
            );
            assert_eq!(coupled.trajectory().velocity_mps(), accepted_v);
            assert_eq!(
                coupled.trajectory().material_history().energy(),
                accepted_energy
            );
            assert_eq!(coupled.source().pole_states(), accepted_poles);
            assert_eq!(coupled.trajectory().steps(), 1);
        }
        gate.request();
        assert!(matches!(
            coupled.advance(&[0.0; 12], cx),
            Err(SolidMimoTransientError::Cancelled)
        ));
        assert_eq!(coupled.source().pole_states(), accepted_poles);
        assert_eq!(
            coupled.trajectory().material_history().displacement_m(),
            accepted_u
        );
    });
}
