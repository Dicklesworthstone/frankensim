//! Original constrained attachments: independent analytic force/pole/power
//! specimens, including zero-storage reactions beside hereditary material.
use fs_exec::Cx;
use fs_la::nonlinear_eig::PoleOptions;
use fs_material::visco::tensor::{Mandel6, PronyBranch, TensorProny};
use fs_math::c64::C64;
use fs_phs::{DescriptorPortHamiltonian, PortHamiltonian, QuadraticStorage};
use fs_solid::linear3::TetAssemblyBudget;
use fs_solid::linear3::visco::*;

fn with_cx(f: impl FnOnce(&Cx<'_>, &fs_exec::CancelGate)) {
    let gate = fs_exec::CancelGate::new_clock_free();
    fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default()).scope(|arena| {
        f(
            &Cx::new(
                &gate,
                arena,
                fs_exec::StreamKey {
                    seed: 419,
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
// The one free displacement has original M=1, K=2.5. With memory it also
// has Kelvin C=0.2 and six original Maxwell strains with k=6.5, tau=0.25.
fn source(memory: bool, cx: &Cx<'_>) -> ViscoElements {
    let branches = if memory {
        vec![PronyBranch::new(diagonal(39.0), 0.25).unwrap()]
    } else {
        vec![]
    };
    let mut law = TensorProny::new_mandel(diagonal(15.0), branches, 0.01).unwrap();
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
            nodes_m: &[
                [0.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
                [0.0, 0.0, 1.0],
            ],
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
fn preparation() -> ViscoPoleBudget {
    ViscoPoleBudget {
        maximum_coordinates: 100,
        maximum_payload_bytes: 8 << 20,
        maximum_assembly_work: 1_000_000,
        maximum_evaluations: 1024,
        reference_rate_per_s: 2.0,
        physical_residual_tolerance: 2e-11,
    }
}
fn dynamic_budget() -> ViscoPoleDynamicBudget {
    ViscoPoleDynamicBudget {
        maximum_states: 20,
        maximum_ports: 8,
        maximum_payload_bytes: 1 << 20,
        maximum_preparation_work: 1_000_000,
    }
}
fn response_budget() -> ViscoJointHarmonicBudget {
    ViscoJointHarmonicBudget {
        maximum_coordinates: 100,
        maximum_payload_bytes: 8 << 20,
        maximum_work: 100_000_000,
    }
}
fn pole_options() -> PoleOptions {
    PoleOptions {
        spectral_scale: 2.0,
        step_tolerance: 2e-12,
        residual_tolerance: 2e-12,
        max_iterations: 40,
        max_dimension: 100,
        max_dense_bytes: 8 << 20,
        max_work: 100_000_000,
    }
}
fn spec<'a>(
    system: &'a DescriptorPortHamiltonian,
    force_map: &'a [f64],
) -> ViscoPoleDescriptorLoad<'a> {
    ViscoPoleDescriptorLoad {
        system,
        force_map,
        band_rad_s: (0.0, 100.0),
    }
}
// p'=-(c/m)*p+lambda; 0=-p/m+v; y=lambda.
// p is the original momentum, lambda a zero-storage constraint reaction.
fn inertance(mass: f64, drag: f64) -> DescriptorPortHamiltonian {
    DescriptorPortHamiltonian::new(
        2,
        1,
        1,
        vec![0.0, 1.0, -1.0, 0.0],
        vec![drag, 0.0, 0.0, 0.0],
        vec![0.0, 1.0],
        Box::new(QuadraticStorage::new(vec![1.0 / mass], 1).unwrap()),
    )
    .unwrap()
}
fn near(a: C64, b: C64, scale: f64) {
    assert!(
        (a - b).abs() < 3e-10 * scale.max(a.abs()).max(b.abs()),
        "{a:?} != {b:?}"
    );
}
fn close(a: f64, b: f64, scale: f64) {
    assert!(
        (a - b).abs() < 3e-10 * scale.max(a.abs()).max(b.abs()),
        "{a:e} != {b:e}"
    );
}
fn check_power(r: &ViscoJointHarmonicResponse) {
    let load = r.dynamic_load.as_ref().unwrap();
    let scale = r.average_input_power_w.abs()
        + r.average_material_heat_w.abs()
        + load.average_dissipation_w.abs()
        + r.omega_rad_s * r.average_total_storage_j.unwrap()
        + 1e-30;
    assert!(r.body_power_balance_residual_w.abs() < 3e-10 * scale);
    assert!(r.total_power_balance_residual_w.unwrap().abs() < 3e-10 * scale);
    assert!(load.power_balance_residual_w.abs() < 3e-10 * scale);
    assert!(r.mechanical_residual < 2e-11);
    assert!(r.memory_residual < 2e-11);
    assert!(r.dynamic_load_residual < 2e-11);
}

#[test]
fn constrained_inertance_retains_original_reaction_in_forced_response_and_finite_poles() {
    with_cx(|cx, _| {
        let body = source(false, cx);
        let (mass, drag) = (2.0, 0.3);
        let load = inertance(mass, drag);
        let force = C64::new(3e-4, -2e-4);
        for b in [0.75, -0.75] {
            let image = body
                .prepare_resonance_poles_with_descriptor_load(
                    spec(&load, &[b]),
                    preparation(),
                    dynamic_budget(),
                    cx,
                )
                .unwrap();
            assert_eq!(image.dimension(), 3);
            assert_eq!(image.dynamic_load_range(), Some(1..3));
            assert_eq!(image.dynamic_load_differential_range(), Some(1..2));
            assert_eq!(image.dynamic_load_algebraic_range(), Some(2..3));
            let total_mass = 1.0 + b * b * mass;
            let total_drag = b * b * drag;
            for omega in [0.1, 0.9, 2.5_f64.sqrt(), 3.0, 12.0] {
                let s = C64::new(0.0, -omega);
                let u =
                    force / (C64::from_re(2.5) + (s * s).scale(total_mass) + s.scale(total_drag));
                let flow = (s * u).scale(b);
                let momentum = flow.scale(mass);
                let reaction = (s.scale(mass) + C64::from_re(drag)) * flow;
                let factor = image.harmonic_factor(omega, response_budget(), cx).unwrap();
                let r = factor.solve(&[force], cx).unwrap();
                let account = r.dynamic_load.as_ref().unwrap();
                near(r.displacement_m[0], u, force.abs());
                near(account.state[0], momentum, force.abs());
                near(account.state[1], reaction, force.abs());
                near(account.flow[0], flow, force.abs());
                near(account.effort[0], reaction, force.abs());
                near(
                    account.state[0].scale(1.0 / mass),
                    account.flow[0],
                    force.abs(),
                );
                close(
                    account.average_storage_j,
                    0.25 * momentum.norm_sq() / mass,
                    1e-25,
                );
                close(
                    account.average_dissipation_w,
                    0.5 * drag * flow.norm_sq(),
                    1e-25,
                );
                close(
                    account.average_absorbed_power_w,
                    0.5 * drag * flow.norm_sq(),
                    1e-25,
                );
                close(
                    r.average_body_storage_j,
                    0.25 * (2.5 * u.norm_sq() + (s * u).norm_sq()),
                    1e-25,
                );
                assert_eq!(r.average_material_heat_w, 0.0);
                check_power(&r);
                // Reuse the same joint LU with a genuinely different phase.
                let phase = C64::new(-0.4, 0.7);
                let second = factor.solve(&[force * phase], cx).unwrap();
                for (&actual, &first) in second.physical_vector.iter().zip(&r.physical_vector) {
                    near(actual, first * phase, force.abs());
                }
            }
            let decay = -total_drag / (2.0 * total_mass);
            let expected = C64::new(decay, (2.5 / total_mass - decay * decay).sqrt());
            let seed = expected + C64::new(-0.02, 0.03);
            let seed_flow = seed.scale(b);
            let state = [
                C64::ONE,
                seed_flow.scale(mass),
                (seed.scale(mass) + C64::from_re(drag)) * seed_flow,
            ];
            let pole = image
                .resonance_pole_from_state(seed, &state, &pole_options(), cx)
                .unwrap();
            near(pole.rate_per_s, expected, 1.0);
            let account = pole.dynamic_load.as_ref().unwrap();
            let flow = (pole.rate_per_s * pole.physical_vector[0]).scale(b);
            near(account.state[0], flow.scale(mass), 1.0);
            near(
                account.state[1],
                (pole.rate_per_s.scale(mass) + C64::from_re(drag)) * flow,
                1.0,
            );
            close(
                account.storage_j,
                0.5 * account.state[0].norm_sq() / mass,
                1e-25,
            );
            close(account.dissipation_w, drag * flow.norm_sq(), 1e-25);
            let scale = 2.0 * pole.rate_per_s.abs() * (account.body_storage_j + account.storage_j)
                + account.dissipation_w;
            assert!(account.signed_power_defect_w.abs() < 3e-10 * scale);
            assert!(account.joint_decay_power_defect_w.abs() < 3e-10 * scale);
            assert!(pole.mechanical_residual < 2e-11 && pole.dynamic_load_residual < 2e-11);
        }
    });
}

#[test]
fn purely_algebraic_mutual_resistance_keeps_zero_storage_and_all_maxwell_strains() {
    with_cx(|cx, _| {
        let body = source(true, cx);
        let load = DescriptorPortHamiltonian::new(
            2,
            0,
            2,
            vec![0.0; 4],
            vec![2.0, 0.5, 0.5, 1.0],
            vec![1.0, 0.0, 0.0, 1.0],
            Box::new(QuadraticStorage::new(vec![], 0).unwrap()),
        )
        .unwrap();
        let image = body
            .prepare_resonance_poles_with_descriptor_load(
                spec(&load, &[1.0, -0.5]),
                preparation(),
                dynamic_budget(),
                cx,
            )
            .unwrap();
        assert_eq!(image.dimension(), 9);
        assert_eq!(image.dynamic_load_differential_range(), Some(7..7));
        assert_eq!(image.dynamic_load_algebraic_range(), Some(7..9));
        let force = C64::new(1.2e-4, 0.9e-4);
        for omega in [0.2, 0.9, 2.0, 3.0, 12.0] {
            let s = C64::new(0.0, -omega);
            let st = s.scale(0.25);
            // R^-1 B^T = [5/7,-6/7]^T, B R^-1 B^T = 8/7.
            let original = C64::from_re(2.5)
                + s * s
                + s.scale(0.2 + 8.0 / 7.0)
                + st.scale(6.5) / (C64::ONE + st);
            let u = force / original;
            let velocity = s * u;
            let q = st * u / (C64::ONE + st);
            let lambda = [velocity.scale(5.0 / 7.0), velocity.scale(-6.0 / 7.0)];
            let r = image
                .harmonic_factor(omega, response_budget(), cx)
                .unwrap()
                .solve(&[force], cx)
                .unwrap();
            near(r.displacement_m[0], u, 1e-12);
            near(r.physical_vector[1], q, 1e-12);
            for &strain in &r.physical_vector[2..7] {
                assert_eq!(strain, C64::ZERO);
            }
            let account = r.dynamic_load.as_ref().unwrap();
            for i in 0..2 {
                near(account.state[i], lambda[i], 1e-12);
                near(account.effort[i], lambda[i], 1e-12);
                near(account.flow[i], velocity.scale([1.0, -0.5][i]), 1e-12);
            }
            // Check both uncondensed constraint rows, not only host impedance.
            near(
                account.state[0].scale(2.0) + account.state[1].scale(0.5),
                velocity,
                1e-12,
            );
            near(
                account.state[0].scale(0.5) + account.state[1],
                velocity.scale(-0.5),
                1e-12,
            );
            assert_eq!(account.average_storage_j, 0.0);
            close(
                account.average_dissipation_w,
                0.5 * (8.0 / 7.0) * velocity.norm_sq(),
                1e-25,
            );
            close(
                r.average_material_heat_w,
                0.5 * (0.2 * velocity.norm_sq() + (6.5 / 0.25) * q.norm_sq()),
                1e-25,
            );
            close(
                r.average_body_storage_j,
                0.25 * (2.5 * u.norm_sq() + velocity.norm_sq() + 6.5 * q.norm_sq()),
                1e-25,
            );
            check_power(&r);
        }
    });
}

#[test]
fn descriptor_with_no_algebraic_coordinates_preserves_existing_ode_response_exactly() {
    with_cx(|cx, _| {
        let body = source(true, cx);
        let j = vec![0.0, -1.0, 1.0, 0.0];
        let r = vec![0.0, 0.0, 0.0, 0.35];
        let g = vec![1.0, -0.6, 0.0, 0.0];
        let q = vec![2.0, 0.0, 0.0, 1.0];
        let ode = PortHamiltonian::new(
            2,
            2,
            j.clone(),
            r.clone(),
            g.clone(),
            Box::new(QuadraticStorage::new(q.clone(), 2).unwrap()),
        )
        .unwrap();
        let descriptor = DescriptorPortHamiltonian::new(
            2,
            2,
            2,
            j,
            r,
            g,
            Box::new(QuadraticStorage::new(q, 2).unwrap()),
        )
        .unwrap();
        let map = [0.8, 0.3];
        let original = body
            .prepare_resonance_poles_with_dynamic_load(
                ViscoPoleDynamicLoad {
                    system: &ode,
                    force_map: &map,
                    band_rad_s: (0.0, 100.0),
                },
                preparation(),
                dynamic_budget(),
                cx,
            )
            .unwrap();
        let constrained = body
            .prepare_resonance_poles_with_descriptor_load(
                spec(&descriptor, &map),
                preparation(),
                dynamic_budget(),
                cx,
            )
            .unwrap();
        assert_eq!(constrained.dynamic_load_differential_range(), Some(7..9));
        assert_eq!(constrained.dynamic_load_algebraic_range(), Some(9..9));
        assert_eq!(original.dynamic_load_algebraic_range(), Some(9..9));
        for omega in [0.2, 2.0, 12.0] {
            let force = [C64::new(1.2e-4, -0.9e-4)];
            let a = original
                .harmonic_factor(omega, response_budget(), cx)
                .unwrap()
                .solve(&force, cx)
                .unwrap();
            let b = constrained
                .harmonic_factor(omega, response_budget(), cx)
                .unwrap()
                .solve(&force, cx)
                .unwrap();
            assert_eq!(a.physical_vector, b.physical_vector);
            assert_eq!(a.average_total_storage_j, b.average_total_storage_j);
            assert_eq!(a.average_material_heat_w, b.average_material_heat_w);
            assert_eq!(
                a.dynamic_load.as_ref().unwrap().average_dissipation_w,
                b.dynamic_load.as_ref().unwrap().average_dissipation_w,
            );
        }
    });
}

#[test]
fn original_constraints_capacity_band_singularity_and_cancellation_refuse_without_regularization() {
    with_cx(|cx, gate| {
        let body = source(true, cx);
        let load = inertance(2.0, 0.3);
        let attach = spec(&load, &[0.75]);
        for limited in [
            ViscoPoleDynamicBudget {
                maximum_states: 1,
                ..dynamic_budget()
            },
            ViscoPoleDynamicBudget {
                maximum_ports: 0,
                ..dynamic_budget()
            },
            ViscoPoleDynamicBudget {
                maximum_payload_bytes: 0,
                ..dynamic_budget()
            },
            ViscoPoleDynamicBudget {
                maximum_preparation_work: 0,
                ..dynamic_budget()
            },
        ] {
            assert!(
                body.prepare_resonance_poles_with_descriptor_load(
                    attach,
                    preparation(),
                    limited,
                    cx,
                )
                .is_err()
            );
        }
        assert!(
            body.prepare_resonance_poles_with_descriptor_load(
                attach,
                ViscoPoleBudget {
                    maximum_coordinates: 8,
                    ..preparation()
                },
                dynamic_budget(),
                cx,
            )
            .is_err()
        );
        assert!(
            body.prepare_resonance_poles_with_descriptor_load(
                spec(&load, &[]),
                preparation(),
                dynamic_budget(),
                cx,
            )
            .is_err()
        );
        let mut limited = attach;
        limited.band_rad_s = (0.0, 1.0);
        let image = body
            .prepare_resonance_poles_with_descriptor_load(
                limited,
                preparation(),
                dynamic_budget(),
                cx,
            )
            .unwrap();
        assert!(image.harmonic_factor(2.0, response_budget(), cx).is_err());
        let factor = image.harmonic_factor(0.5, response_budget(), cx).unwrap();
        let force = [C64::new(1e-4, 2e-4)];
        let before = factor.solve(&force, cx).unwrap();
        assert!(factor.solve(&[C64::from_re(100.0)], cx).is_err());
        assert!(factor.solve(&[], cx).is_err());
        assert_eq!(
            before.physical_vector,
            factor.solve(&force, cx).unwrap().physical_vector
        );

        // Both original algebraic equations demand the same zero velocity;
        // the two reactions are indeterminate. They must stay singular.
        let redundant = DescriptorPortHamiltonian::new(
            2,
            0,
            1,
            vec![0.0; 4],
            vec![0.0; 4],
            vec![1.0, 1.0],
            Box::new(QuadraticStorage::new(vec![], 0).unwrap()),
        )
        .unwrap();
        match body.prepare_resonance_poles_with_descriptor_load(
            spec(&redundant, &[1.0]),
            preparation(),
            dynamic_budget(),
            cx,
        ) {
            Ok(singular) => assert!(
                singular
                    .harmonic_factor(0.5, response_budget(), cx)
                    .is_err()
            ),
            Err(ViscoPoleError::Invalid(_)) => {}
            Err(other) => panic!("unexpected refusal: {other:?}"),
        }
        gate.request();
        assert!(factor.solve(&force, cx).is_err());
        assert!(image.harmonic_factor(0.5, response_budget(), cx).is_err());
        assert!(
            body.prepare_resonance_poles_with_descriptor_load(
                attach,
                preparation(),
                dynamic_budget(),
                cx,
            )
            .is_err()
        );
    });
}
