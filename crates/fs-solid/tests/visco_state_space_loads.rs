//! G1/G3/G4 original continuous boundary states, signed exchange and inertance.
use fs_exec::Cx;
use fs_la::nonlinear_eig::PoleOptions;
use fs_material::visco::tensor::{Mandel6, PronyBranch, TensorProny};
use fs_math::c64::C64;
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

fn attachment() -> ViscoPoleStateSpaceLoad<'static> {
    // x''+4x=v, effort=x'+0.4v'. A third original state z'=-3z
    // is mechanically dark and must not disappear from the joint spectrum.
    ViscoPoleStateSpaceLoad {
        states: 3,
        ports: 1,
        a: &[0.0, 1.0, 0.0, -4.0, 0.0, 0.0, 0.0, 0.0, -3.0],
        b: &[0.0, 1.0, 0.0],
        c: &[0.0, 1.0, 0.0],
        d: &[0.0],
        e: &[0.4],
        force_map: &[1.0],
        band_rad_s: (0.0, 100.0),
    }
}
fn near(a: C64, b: C64, scale: f64) {
    assert!(
        (a - b).abs() < 4e-10 * scale.max(a.abs()).max(b.abs()),
        "{a:?} != {b:?}"
    );
}
fn close(a: f64, b: f64, scale: f64) {
    near(C64::from_re(a), C64::from_re(b), scale);
}
fn check_boundary_account(r: &ViscoJointHarmonicResponse) {
    assert_eq!(r.average_total_storage_j, None);
    assert_eq!(r.total_power_balance_residual_w, None);
    assert!(r.port_load.is_none() && r.dynamic_load.is_none());
    let load = r.state_space_load.as_ref().unwrap();
    let scale = r.average_input_power_w.abs()
        + r.average_material_heat_w.abs()
        + load.average_absorbed_power_w.abs()
        + r.omega_rad_s * r.average_body_storage_j
        + 1e-30;
    assert!(r.body_power_balance_residual_w.abs() < 4e-10 * scale);
    assert!(r.mechanical_residual < 2e-11);
    assert!(r.memory_residual < 2e-11);
    assert!(r.dynamic_load_residual < 2e-11);
}

#[test]
fn g1_original_oscillator_and_derivative_inertance_cross_bare_resonances() {
    with_cx(|cx, _| {
        let body = source(false, cx);
        let image = body
            .prepare_resonance_poles_with_state_space_load(
                attachment(),
                preparation(),
                dynamic_budget(),
                cx,
            )
            .unwrap();
        assert_eq!(image.dimension(), 4);
        assert_eq!(image.state_space_load_range(), Some(1..4));
        assert_eq!(image.dynamic_load_range(), None);
        let force = C64::new(2e-4, -3e-4);
        for omega in [0.2, (2.5_f64 / 1.4).sqrt(), 2.0, 3.0, 12.0] {
            let s = C64::new(0.0, -omega);
            let bare = s * s + C64::from_re(4.0);
            let denominator = (C64::from_re(2.5) + (s * s).scale(1.4)) * bare + s * s;
            let u = force * bare / denominator;
            let x = s * force / denominator;
            let dx = s * x;
            let velocity = s * u;
            let effort = dx + (s * velocity).scale(0.4);
            let factor = image.harmonic_factor(omega, response_budget(), cx).unwrap();
            let r = factor.solve(&[force], cx).unwrap();
            let load = r.state_space_load.as_ref().unwrap();
            near(r.displacement_m[0], u, force.abs());
            near(load.state[0], x, force.abs());
            near(load.state[1], dx, force.abs());
            assert_eq!(load.state[2], C64::ZERO);
            near(load.flow[0], velocity, force.abs());
            near(load.effort[0], effort, force.abs());
            near((s * s + C64::from_re(2.5)) * u + effort, force, force.abs());
            close(
                r.average_body_storage_j,
                0.25 * (2.5 * u.norm_sq() + velocity.norm_sq()),
                1e-25,
            );
            assert_eq!(r.average_material_heat_w, 0.0);
            check_boundary_account(&r);
            if omega == 2.0 {
                assert!(u.abs() < 1e-20);
                near(load.state[1], force, force.abs());
            }
            let phase = C64::new(-0.2, 0.6);
            let next = factor.solve(&[force * phase], cx).unwrap();
            for (a, b) in next.physical_vector.iter().zip(&r.physical_vector) {
                near(*a, *b * phase, force.abs());
            }
        }
        // Exact two-oscillator polynomial including the ORIGINAL E=0.4:
        // 1.4*w^4-(2.5+1.4*4+1)*w^2+2.5*4=0.
        let a = 1.4_f64;
        let b = 2.5 + a * 4.0 + 1.0;
        let discriminant = (b * b - 4.0 * a * 10.0).sqrt();
        for squared in [
            (b - discriminant) / (2.0 * a),
            (b + discriminant) / (2.0 * a),
        ] {
            let expected = C64::new(0.0, squared.sqrt());
            let seed = expected + C64::new(-0.01, 0.02);
            let x = seed / (seed * seed + C64::from_re(4.0));
            let pole = image
                .resonance_pole_from_state(
                    seed,
                    &[C64::ONE, x, seed * x, C64::ZERO],
                    &pole_options(),
                    cx,
                )
                .unwrap();
            near(pole.rate_per_s, expected, 1.0);
            let load = pole.state_space_load.as_ref().unwrap();
            let scale = 2.0 * pole.rate_per_s.abs() * load.body_storage_j
                + load.absorbed_power_w.abs()
                + 1e-30;
            assert!(load.body_decay_power_residual_w.abs() < 4e-10 * scale);
            assert!(pole.dynamic_load.is_none());
        }
    });
}

#[test]
fn g3_full_mutual_direct_and_derivative_terms_keep_signed_work_and_original_material_heat() {
    with_cx(|cx, _| {
        let body = source(true, cx);
        let load = ViscoPoleStateSpaceLoad {
            states: 0,
            ports: 2,
            a: &[],
            b: &[],
            c: &[],
            d: &[1.0, 2.0, 2.0, 4.0],
            e: &[0.4, 0.1, 0.1, 0.2],
            force_map: &[1.0, -0.25],
            band_rad_s: (0.0, 100.0),
        };
        let image = body
            .prepare_resonance_poles_with_state_space_load(
                load,
                preparation(),
                dynamic_budget(),
                cx,
            )
            .unwrap();
        assert_eq!(image.dimension(), 7);
        assert_eq!(image.state_space_load_range(), Some(7..7));
        let force = C64::new(1e-4, -2e-4);
        for omega in [0.2, 0.9, 2.0, 3.0, 12.0] {
            let s = C64::new(0.0, -omega);
            let st = s.scale(0.25);
            let denominator = C64::from_re(2.5)
                + (s * s).scale(1.0 + 0.3625)
                + s.scale(0.2 + 0.25)
                + st.scale(6.5) / (C64::ONE + st);
            let u = force / denominator;
            let v = s * u;
            let q = st * u / (C64::ONE + st);
            let r = image
                .harmonic_factor(omega, response_budget(), cx)
                .unwrap()
                .solve(&[force], cx)
                .unwrap();
            let l = r.state_space_load.as_ref().unwrap();
            assert!(l.state.is_empty());
            near(r.displacement_m[0], u, force.abs());
            near(r.physical_vector[1], q, force.abs());
            for strain in &r.physical_vector[2..7] {
                assert_eq!(*strain, C64::ZERO);
            }
            let flow = [v, v.scale(-0.25)];
            for port in 0..2 {
                let mut expected = C64::ZERO;
                for other in 0..2 {
                    expected = expected
                        + flow[other].scale(load.d[2 * port + other])
                        + (s * flow[other]).scale(load.e[2 * port + other]);
                }
                near(l.flow[port], flow[port], force.abs());
                near(l.effort[port], expected, force.abs());
                close(
                    l.average_absorbed_port_power_w[port],
                    0.5 * (flow[port].conj() * expected).re,
                    1e-25,
                );
            }
            assert!(l.average_absorbed_port_power_w[0] > 0.0);
            assert!(l.average_absorbed_port_power_w[1] < 0.0);
            close(l.average_absorbed_power_w, 0.5 * 0.25 * v.norm_sq(), 1e-25);
            close(
                r.average_material_heat_w,
                0.5 * (0.2 * v.norm_sq() + (6.5 / 0.25) * q.norm_sq()),
                1e-25,
            );
            close(
                r.average_body_storage_j,
                0.25 * (2.5 * u.norm_sq() + v.norm_sq() + 6.5 * q.norm_sq()),
                1e-25,
            );
            check_boundary_account(&r);
        }
    });
}

#[test]
fn g3_original_hidden_load_state_remains_in_the_finite_joint_spectrum() {
    with_cx(|cx, _| {
        let body = source(false, cx);
        let image = body
            .prepare_resonance_poles_with_state_space_load(
                attachment(),
                preparation(),
                dynamic_budget(),
                cx,
            )
            .unwrap();
        let pole = image
            .resonance_pole_from_state(
                C64::from_re(-3.01),
                &[C64::ZERO, C64::ZERO, C64::ZERO, C64::ONE],
                &pole_options(),
                cx,
            )
            .unwrap();
        near(pole.rate_per_s, C64::from_re(-3.0), 1.0);
        assert_eq!(pole.physical_vector.len(), 4);
        let l = pole.state_space_load.unwrap();
        assert!(l.state[2].abs() > 0.0);
        for i in 0..3 {
            assert_eq!(pole.physical_vector[i], C64::ZERO);
        }
        assert_eq!(l.absorbed_power_w, 0.0);
        assert_eq!(l.body_storage_j, 0.0);
        assert_eq!(l.body_decay_power_residual_w, 0.0);
    });
}

#[test]
fn g4_original_matrix_budget_band_strain_and_cancellation_refusals_preserve_factor() {
    with_cx(|cx, gate| {
        let body = source(false, cx);
        for limits in [
            ViscoPoleDynamicBudget {
                maximum_states: 2,
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
                body.prepare_resonance_poles_with_state_space_load(
                    attachment(),
                    preparation(),
                    limits,
                    cx,
                )
                .is_err()
            );
        }
        for invalid_load in [
            ViscoPoleStateSpaceLoad {
                e: &[],
                ..attachment()
            },
            ViscoPoleStateSpaceLoad {
                c: &[0.0, 1.0],
                ..attachment()
            },
            ViscoPoleStateSpaceLoad {
                d: &[f64::NAN],
                ..attachment()
            },
            ViscoPoleStateSpaceLoad {
                force_map: &[],
                ..attachment()
            },
            ViscoPoleStateSpaceLoad {
                band_rad_s: (2.0, 1.0),
                ..attachment()
            },
        ] {
            assert!(
                body.prepare_resonance_poles_with_state_space_load(
                    invalid_load,
                    preparation(),
                    dynamic_budget(),
                    cx,
                )
                .is_err()
            );
        }
        let image = body
            .prepare_resonance_poles_with_state_space_load(
                ViscoPoleStateSpaceLoad {
                    band_rad_s: (0.1, 8.0),
                    ..attachment()
                },
                preparation(),
                dynamic_budget(),
                cx,
            )
            .unwrap();
        assert!(image.harmonic_factor(0.05, response_budget(), cx).is_err());
        assert!(image.harmonic_factor(9.0, response_budget(), cx).is_err());
        let factor = image.harmonic_factor(0.5, response_budget(), cx).unwrap();
        let force = [C64::new(1e-4, -2e-4)];
        let before = factor.solve(&force, cx).unwrap();
        assert!(factor.solve(&[], cx).is_err());
        assert!(factor.solve(&[C64::from_re(1e4)], cx).is_err());
        assert_eq!(
            before.physical_vector,
            factor.solve(&force, cx).unwrap().physical_vector
        );
        gate.request();
        assert!(factor.solve(&force, cx).is_err());
        assert!(
            body.prepare_resonance_poles_with_state_space_load(
                attachment(),
                preparation(),
                dynamic_budget(),
                cx,
            )
            .is_err()
        );
    });
}
