//! G1/G3/G4: original joint forced states, analytical attachments and power.
use fs_exec::Cx;
use fs_material::visco::tensor::{Mandel6, PronyBranch, TensorProny};
use fs_math::c64::C64;
use fs_phs::{PortHamiltonian, QuadraticStorage};
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
                    kernel_id: 95,
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
fn diagonal(v: f64) -> Mandel6 {
    core::array::from_fn(|i| core::array::from_fn(|j| if i == j { v } else { 0.0 }))
}
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
        material_state_identity: [95; 32],
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
fn response_budget() -> ViscoJointHarmonicBudget {
    ViscoJointHarmonicBudget {
        maximum_coordinates: 100,
        maximum_payload_bytes: 8 << 20,
        maximum_work: 100_000_000,
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
fn port_budget() -> ViscoPolePortBudget {
    ViscoPolePortBudget {
        maximum_ports: 8,
        maximum_payload_bytes: 1 << 20,
        maximum_preparation_work: 1_000_000,
    }
}
// State [host-minus-mass spring extension; independent mass momentum].
// Host velocity drives extension; original spring reaction opposes the host.
fn attachment(drag: f64, ports: usize, g: Vec<f64>) -> PortHamiltonian {
    PortHamiltonian::new(
        2,
        ports,
        vec![0.0, -1.0, 1.0, 0.0],
        vec![0.0, 0.0, 0.0, drag],
        g,
        Box::new(QuadraticStorage::new(vec![2.0, 0.0, 0.0, 1.0], 2).unwrap()),
    )
    .unwrap()
}
fn spec<'a>(system: &'a PortHamiltonian, b: &'a [f64]) -> ViscoPoleDynamicLoad<'a> {
    ViscoPoleDynamicLoad {
        system,
        force_map: b,
        band_rad_s: (0.0, 100.0),
    }
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
fn power(r: &ViscoJointHarmonicResponse) {
    let scale = r.average_input_power_w.abs()
        + r.average_material_heat_w.abs()
        + r.port_load
            .iter()
            .chain(&r.dynamic_load)
            .map(|l| l.average_dissipation_w.abs())
            .sum::<f64>()
        + r.omega_rad_s * r.average_total_storage_j.unwrap()
        + 1e-30;
    assert!(
        r.body_power_balance_residual_w.abs() < 3e-10 * scale,
        "{r:?}"
    );
    assert!(
        r.total_power_balance_residual_w.unwrap().abs() < 3e-10 * scale,
        "{r:?}"
    );
    assert!(
        r.mechanical_residual < 2e-11
            && r.memory_residual < 2e-11
            && r.dynamic_load_residual < 2e-11
    );
}

#[test]
fn g1_bare_attachment_antiresonance_and_bare_body_resonance_keep_regular_joint_response() {
    with_cx(|cx, _| {
        let body = source(false, cx);
        let load = attachment(0.0, 1, vec![1.0, 0.0]);
        let image = body
            .prepare_resonance_poles_with_dynamic_load(
                spec(&load, &[1.0]),
                preparation(),
                dynamic_budget(),
                cx,
            )
            .unwrap();
        let force = C64::new(3e-4, -2e-4);
        for omega in [0.1, 2_f64.sqrt(), 2.5_f64.sqrt(), 2.1, 7.0] {
            let factor = image.harmonic_factor(omega, response_budget(), cx).unwrap();
            let r = factor.solve(&[force], cx).unwrap();
            let s = C64::new(0.0, -omega);
            // Independent two-coordinate equation, with no division by the
            // standalone attachment denominator d=s²+2.
            let a = C64::from_re(2.5) + s * s;
            let d = C64::from_re(2.0) + s * s;
            let det = a * d + (s * s).scale(2.0);
            let u = force * d / det;
            let q = force * s * s / det;
            let p = s * (u - q);
            near(r.displacement_m[0], u, force.abs());
            assert_eq!(r.physical_vector.len(), 3);
            let account = r.dynamic_load.as_ref().unwrap();
            near(account.state[0], q, force.abs());
            near(account.state[1], p, force.abs());
            near(account.effort[0], q.scale(2.0), force.abs());
            close(
                r.average_body_storage_j,
                0.25 * (2.5 * u.norm_sq() + (s * u).norm_sq()),
                force.norm_sq(),
            );
            close(
                account.average_storage_j,
                0.25 * (2.0 * q.norm_sq() + p.norm_sq()),
                force.norm_sq(),
            );
            assert_eq!(r.average_material_heat_w, 0.0);
            assert_eq!(account.average_dissipation_w, 0.0);
            if omega == 2_f64.sqrt() {
                assert!(r.displacement_m[0].abs() < 1e-14 * force.abs());
                assert!(account.state[0].abs() > 1e-5);
            }
            power(&r);
            let phase = C64::new(-0.4, 0.7);
            let second = factor.solve(&[force * phase], cx).unwrap();
            for (x, y) in second.physical_vector.iter().zip(&r.physical_vector) {
                near(*x, *y * phase, force.abs());
            }
        }
    });
}

#[test]
fn g1_original_tensor_memory_kelvin_and_mutual_dynamic_load_have_independent_heat() {
    with_cx(|cx, _| {
        let body = source(true, cx);
        let drag = 0.35;
        let load = attachment(drag, 2, vec![1.0, -0.6, 0.0, 0.0]);
        let map = [0.8, 0.3];
        let h = map[0] - 0.6 * map[1];
        let image = body
            .prepare_resonance_poles_with_dynamic_load(
                spec(&load, &map),
                preparation(),
                dynamic_budget(),
                cx,
            )
            .unwrap();
        let force = C64::new(1.2e-4, 0.9e-4);
        for omega in [0.2, 0.9, 2_f64.sqrt(), 3.0, 12.0] {
            let s = C64::new(0.0, -omega);
            let st = s.scale(0.25);
            let original =
                C64::from_re(2.5) + s * s + s.scale(0.2) + st.scale(6.5) / (C64::ONE + st);
            let d = s * s + s.scale(drag) + C64::from_re(2.0);
            let excitation = s * (s + C64::from_re(drag));
            let det = original * d + excitation.scale(2.0 * h * h);
            let u = force * d / det;
            let q = force * excitation.scale(h) / det;
            let p = s * (u.scale(h) - q);
            let r = image
                .harmonic_factor(omega, response_budget(), cx)
                .unwrap()
                .solve(&[force], cx)
                .unwrap();
            assert_eq!(r.physical_vector.len(), 9);
            near(r.displacement_m[0], u, 1e-12);
            let memory = st * u / (C64::ONE + st);
            near(r.physical_vector[1], memory, 1e-12);
            for &x in &r.physical_vector[2..7] {
                assert_eq!(x, C64::ZERO);
            }
            let account = r.dynamic_load.as_ref().unwrap();
            near(account.state[0], q, 1e-12);
            near(account.state[1], p, 1e-12);
            for k in 0..2 {
                near(account.flow[k], (s * u).scale(map[k]), 1e-12);
                near(account.effort[k], q.scale(2.0 * [1.0, -0.6][k]), 1e-12);
            }
            let heat = 0.5 * (0.2 * (s * u).norm_sq() + 6.5 / 0.25 * memory.norm_sq());
            close(r.average_material_heat_w, heat, 1e-25);
            close(
                account.average_dissipation_w,
                0.5 * drag * p.norm_sq(),
                1e-25,
            );
            close(
                r.average_body_storage_j,
                0.25 * (2.5 * u.norm_sq() + (s * u).norm_sq() + 6.5 * memory.norm_sq()),
                1e-25,
            );
            assert!(r.average_material_heat_w > 0.0 && account.average_dissipation_w > 0.0);
            power(&r);
        }
    });
}

#[test]
fn g3_original_body_storage_is_resolved_beside_a_much_larger_attached_port_store() {
    with_cx(|cx, _| {
        let body = source(false, cx);
        let (k, c, m, b) = (1e12, 1e3, 1e8, 1.2);
        let image = body
            .prepare_resonance_poles_with_port_load(
                ViscoPolePortLoad {
                    ports: 1,
                    force_map: &[b],
                    stiffness: &[k],
                    damping: &[c],
                    inertance: &[m],
                    band_rad_s: (0.0, 100.0),
                },
                preparation(),
                port_budget(),
                cx,
            )
            .unwrap();
        let omega = 2.0;
        let force = C64::new(0.001, 0.0004);
        let r = image
            .harmonic_factor(omega, response_budget(), cx)
            .unwrap()
            .solve(&[force], cx)
            .unwrap();
        let u = force
            / C64::new(
                2.5 + b * b * k - omega * omega * (1.0 + b * b * m),
                -omega * b * b * c,
            );
        let v = C64::new(0.0, -omega) * u;
        near(r.displacement_m[0], u, 1e-30);
        let p = r.port_load.as_ref().unwrap();
        assert!(p.state.is_empty());
        close(
            r.average_body_storage_j,
            0.25 * (2.5 * u.norm_sq() + v.norm_sq()),
            1e-45,
        );
        close(
            p.average_storage_j,
            0.25 * b * b * (k * u.norm_sq() + m * v.norm_sq()),
            1e-40,
        );
        close(
            p.average_dissipation_w,
            0.5 * b * b * c * v.norm_sq(),
            1e-40,
        );
        assert!(p.average_storage_j > 1e10 * r.average_body_storage_j);
        assert!(r.average_body_storage_j > 0.0);
        assert_eq!(r.average_material_heat_w, 0.0);
        power(&r);
    });
}

#[test]
fn g4_full_cycle_strain_joint_singularity_budgets_and_cancellation_preserve_reusable_factor() {
    with_cx(|cx, gate| {
        let body = source(false, cx);
        let image = body.prepare_resonance_poles(preparation(), cx).unwrap();
        let factor = image.harmonic_factor(1.0, response_budget(), cx).unwrap();
        let force = C64::new(1e-4, 2e-4);
        let before = factor.solve(&[force], cx).unwrap();
        // Each component phase is below the source limit, but the full-cycle
        // strain amplitude sqrt(0.008²+0.008²) exceeds the declared 0.01.
        assert!(factor.solve(&[C64::new(0.012, 0.012)], cx).is_err());
        assert!(factor.solve(&[], cx).is_err());
        assert!(factor.solve(&[C64::new(f64::NAN, 0.0)], cx).is_err());
        for limited in [
            ViscoJointHarmonicBudget {
                maximum_coordinates: 0,
                ..response_budget()
            },
            ViscoJointHarmonicBudget {
                maximum_payload_bytes: 0,
                ..response_budget()
            },
            ViscoJointHarmonicBudget {
                maximum_work: 0,
                ..response_budget()
            },
        ] {
            assert!(image.harmonic_factor(1.0, limited, cx).is_err());
        }
        let after = factor.solve(&[force], cx).unwrap();
        assert_eq!(after.physical_vector, before.physical_vector);
        let zero = factor.solve(&[C64::ZERO], cx).unwrap();
        assert!(zero.physical_vector.iter().all(|&x| x == C64::ZERO));
        let singular = body
            .prepare_resonance_poles_with_port_load(
                ViscoPolePortLoad {
                    ports: 1,
                    force_map: &[1.0],
                    stiffness: &[1.5],
                    damping: &[0.0],
                    inertance: &[0.0],
                    band_rad_s: (0.0, 100.0),
                },
                preparation(),
                port_budget(),
                cx,
            )
            .unwrap();
        assert!(matches!(
            singular.harmonic_factor(2.0, response_budget(), cx),
            Err(ViscoJointHarmonicError::Linear(
                fs_la::eigen_complex::LuComplexError::Singular { .. }
            ))
        ));
        assert!(
            singular
                .harmonic_factor(101.0, response_budget(), cx)
                .is_err()
        );
        gate.request();
        assert!(factor.solve(&[force], cx).is_err());
        assert!(image.harmonic_factor(1.0, response_budget(), cx).is_err());
    });
}
