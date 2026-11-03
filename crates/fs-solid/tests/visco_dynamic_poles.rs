//! G1/G3/G4: independent attachment dynamics on original hereditary solids.
//! Analytical software specimens; no measured-material or spatial certification.
use fs_exec::Cx;
use fs_la::nonlinear_eig::PoleOptions;
use fs_la::nonlinear_eig::beyn::{BeynOptions, Ellipse};
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
                    kernel_id: 94,
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
    core::array::from_fn(|i| core::array::from_fn(|j| if i == j { value } else { 0. }))
}
fn solid(hereditary: bool, cx: &Cx<'_>) -> ViscoElements {
    let branches = if hereditary {
        vec![PronyBranch::new(diagonal(39.), 0.25).unwrap()]
    } else {
        vec![]
    };
    let material = TetViscoMaterial {
        law: TensorProny::new_mandel(diagonal(15.), branches, 0.01).unwrap(),
        density_kg_m3: 60.,
        material_state_identity: [94; 32],
    };
    let fixed: Vec<_> = (0..12).filter(|&i| i != 3).collect();
    ViscoElements::new(
        SolidMeshRef::Tet4 {
            nodes_m: &[[0., 0., 0.], [1., 0., 0.], [0., 1., 0.], [0., 0., 1.]],
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
fn budget() -> ViscoPoleBudget {
    ViscoPoleBudget {
        maximum_coordinates: 100,
        maximum_payload_bytes: 8 * 1024 * 1024,
        maximum_assembly_work: 1_000_000,
        maximum_evaluations: 1024,
        reference_rate_per_s: 2.,
        physical_residual_tolerance: 2e-11,
    }
}
fn load_budget() -> ViscoPoleDynamicBudget {
    ViscoPoleDynamicBudget {
        maximum_states: 20,
        maximum_ports: 8,
        maximum_payload_bytes: 1024 * 1024,
        maximum_preparation_work: 1_000_000,
    }
}
fn local() -> PoleOptions {
    PoleOptions {
        spectral_scale: 2.,
        step_tolerance: 2e-12,
        residual_tolerance: 2e-12,
        max_iterations: 40,
        max_dimension: 100,
        max_dense_bytes: 8 * 1024 * 1024,
        max_work: 100_000_000,
    }
}
fn contour(n: usize, nodes: usize) -> BeynOptions {
    BeynOptions {
        nodes,
        probe_columns: n,
        rank_relative_tolerance: 2e-9,
        count_tolerance: 1e-6,
        quadrature_tolerance: 2e-7,
        contour_min_reciprocal_condition: 1e-10,
        linear_residual_tolerance: 2e-11,
        svd_tolerance: 2e-14,
        svd_max_sweeps: 100,
        eigen_max_sweeps: 100,
        max_dimension: 100,
        max_dense_bytes: 16 * 1024 * 1024,
        max_work: 2_000_000_000,
        refinement: local(),
    }
}
// z=[relative spring extension; attached mass momentum], k=2, m=1.
// The attached mass has physical drag c to ground. Host velocity is the input,
// spring reaction the output; independent mass motion creates two resonances.
fn attachment(drag: f64, ports: usize, g: Vec<f64>) -> PortHamiltonian {
    PortHamiltonian::new(
        2,
        ports,
        vec![0., -1., 1., 0.],
        vec![0., 0., 0., drag],
        g,
        Box::new(QuadraticStorage::new(vec![2., 0., 0., 1.], 2).unwrap()),
    )
    .unwrap()
}
fn spec<'a>(system: &'a PortHamiltonian, force_map: &'a [f64]) -> ViscoPoleDynamicLoad<'a> {
    ViscoPoleDynamicLoad {
        system,
        force_map,
        band_rad_s: (0., 100.),
    }
}
fn near(a: C64, b: C64, tol: f64) {
    assert!(
        (a - b).abs() <= tol * b.abs().max(1.),
        "actual {a:?}; expected {b:?}"
    );
}
fn inner(a: C64, b: C64) -> f64 {
    (a.conj() * b).re
}

#[test]
fn g1_attached_mass_splits_modes_and_retains_states_at_the_bare_load_resonance() {
    with_cx(|cx, _| {
        let source = solid(false, cx);
        let load = attachment(0., 1, vec![1., 0.]);
        let prepared = source
            .prepare_resonance_poles_with_dynamic_load(
                spec(&load, &[1.]),
                budget(),
                load_budget(),
                cx,
            )
            .unwrap();
        assert_eq!(prepared.dimension(), 3);
        assert_eq!(prepared.dynamic_load_range(), Some(1..3));
        let frequencies = [
            ((6.5 - 22.25_f64.sqrt()) / 2.).sqrt(),
            ((6.5 + 22.25_f64.sqrt()) / 2.).sqrt(),
        ];
        for w in frequencies {
            let pole = prepared
                .resonance_pole(C64::new(-0.03, w * 1.04), &[C64::ONE], &local(), cx)
                .unwrap();
            near(pole.rate_per_s, C64::new(0., w), 3e-11);
            assert!(pole.mechanical_residual < 2e-11 && pole.dynamic_load_residual < 2e-11);
            let u = pole.physical_vector[0];
            let q = pole.physical_vector[1];
            let p = pole.physical_vector[2];
            near(
                q / u,
                (pole.rate_per_s * pole.rate_per_s)
                    / (pole.rate_per_s * pole.rate_per_s + C64::from_re(2.)),
                3e-11,
            );
            near(p, pole.rate_per_s * (u - q), 3e-11);
            let account = pole.dynamic_load.unwrap();
            assert_eq!(account.dissipation_w, 0.);
            assert_eq!(account.body_dissipation_w, 0.);
            assert!(account.storage_j > 0. && account.body_storage_j > 0.);
            assert!(
                account.joint_decay_power_defect_w.abs()
                    < 3e-11 * (account.storage_j + account.body_storage_j)
            );
        }
        // Evaluation at the standalone attachment's pole is regular for the
        // JOINT pencil. Its own resolvent would be singular at i*sqrt(2).
        let pole = prepared
            .resonance_pole_from_state(
                C64::new(0., 2_f64.sqrt()),
                &[C64::ONE, C64::ONE, C64::new(0., -2_f64.sqrt())],
                &local(),
                cx,
            )
            .unwrap();
        assert!(
            frequencies
                .iter()
                .any(|&w| (pole.rate_per_s - C64::new(0., w)).abs() < 1e-10)
        );
        let region = Ellipse {
            center: C64::new(0., 1.65),
            real_radius: 0.6,
            imaginary_radius: 0.9,
        };
        assert!(matches!(
            prepared.resonance_spectrum(region, &contour(3, 64), cx),
            Err(ViscoPoleError::Contour(
                fs_la::nonlinear_eig::beyn::BeynError::CountUnresolved { .. }
            ))
        ));
        for nodes in [128, 256] {
            let spectrum = prepared
                .resonance_spectrum(region, &contour(3, nodes), cx)
                .unwrap();
            assert_eq!(
                (spectrum.argument_count, spectrum.rank, spectrum.poles.len()),
                (2, 2, 2)
            );
            for w in frequencies {
                assert!(
                    spectrum
                        .poles
                        .iter()
                        .any(|p| (p.rate_per_s - C64::new(0., w)).abs() < 3e-10)
                );
            }
        }
    });
}

fn expected_hereditary(mut s: C64, drag: f64, h: [f64; 2]) -> C64 {
    // Independent scalar eliminant for checking only: the production joint
    // solver never performs either this material or attachment division.
    for _ in 0..50 {
        let t = C64::ONE + s.scale(0.25);
        let body = s * s + C64::from_re(2.5) + s.scale(1.625) / t;
        let body_d = s.scale(2.) + C64::from_re(1.625) / (t * t);
        let load = s * s + s.scale(drag) + C64::from_re(2.);
        let load_d = s.scale(2.) + C64::from_re(drag);
        let a = 2. * h[0] * h[0] + h[1] * h[1];
        let b = 2. * drag * h[0] * h[0];
        let f = body * load + (s * s).scale(a) + s.scale(b);
        let df = body_d * load + body * load_d + s.scale(2. * a) + C64::from_re(b);
        let step = f / df;
        s = s - step;
        if step.abs() < 1e-14 {
            break;
        }
    }
    s
}

#[test]
fn g1_hereditary_body_and_independent_attachment_have_separate_physical_losses() {
    with_cx(|cx, _| {
        let source = solid(true, cx);
        let load = attachment(0.3, 1, vec![1., 0.]);
        let prepared = source
            .prepare_resonance_poles_with_dynamic_load(
                spec(&load, &[1.]),
                budget(),
                load_budget(),
                cx,
            )
            .unwrap();
        assert_eq!(prepared.dimension(), 9);
        assert_eq!(prepared.memory_sites().len(), 1);
        assert_eq!(prepared.dynamic_load_range(), Some(7..9));
        let mut rates = Vec::new();
        for seed in [C64::new(-0.2, 1.), C64::new(-1., 2.7)] {
            let expected = expected_hereditary(seed, 0.3, [1., 0.]);
            let pole = prepared
                .resonance_pole(seed, &[C64::ONE], &local(), cx)
                .unwrap();
            near(pole.rate_per_s, expected, 5e-11);
            assert!(pole.rate_per_s.re < 0.);
            let s = pole.rate_per_s;
            let x = pole.physical_vector[0];
            let q = pole.physical_vector[1];
            near(q / x, s.scale(0.25) / (C64::ONE + s.scale(0.25)), 3e-11);
            let a = pole.dynamic_load.unwrap();
            let z = &a.state;
            let energy = 0.5 * (2. * z[0].norm_sq() + z[1].norm_sq());
            let heat = 0.3 * z[1].norm_sq();
            let body_energy = 0.5 * (2.5 * x.norm_sq() + (s * x).norm_sq() + 6.5 * q.norm_sq());
            let body_heat = 26. * q.norm_sq();
            assert!((a.storage_j - energy).abs() < 2e-13 * energy);
            assert!((a.dissipation_w - heat).abs() < 2e-13 * heat);
            assert!((a.body_storage_j - body_energy).abs() < 2e-13 * body_energy);
            assert!((a.body_dissipation_w - body_heat).abs() < 2e-13 * body_heat);
            assert!((a.absorbed_power_w - inner(s * x, z[0].scale(2.))).abs() < 1e-12);
            let scale = energy + body_energy + heat + body_heat;
            assert!(a.signed_power_defect_w.abs() < 5e-11 * scale);
            assert!(a.joint_decay_power_defect_w.abs() < 5e-11 * scale);
            rates.push(s);
        }
        assert!((rates[0] - rates[1]).abs() > 0.5);
    });
}

#[test]
fn g3_signed_mutual_ports_and_mechanically_dark_load_states_are_retained() {
    with_cx(|cx, _| {
        let source = solid(true, cx);
        let load = attachment(0.3, 2, vec![1., -1., 0.4, 0.6]);
        let b = [0.7, -0.2];
        let prepared = source
            .prepare_resonance_poles_with_dynamic_load(spec(&load, &b), budget(), load_budget(), cx)
            .unwrap();
        let h = [b[0] - b[1], 0.4 * b[0] + 0.6 * b[1]];
        let seed = C64::new(-0.2, 1.);
        let pole = prepared
            .resonance_pole(seed, &[C64::ONE], &local(), cx)
            .unwrap();
        near(pole.rate_per_s, expected_hereditary(seed, 0.3, h), 5e-11);
        let x = pole.physical_vector[0];
        let a = pole.dynamic_load.unwrap();
        for i in 0..2 {
            near(a.flow[i], pole.rate_per_s * x.scale(b[i]), 2e-12);
        }
        near(
            a.effort[0],
            a.state[0].scale(2.) + a.state[1].scale(0.4),
            2e-12,
        );
        near(
            a.effort[1],
            a.state[0].scale(-2.) + a.state[1].scale(0.6),
            2e-12,
        );
        assert!(a.signed_power_defect_w.abs() < 1e-10);
        // A disconnected original attachment still has its physical poles.
        // They are not an observable mechanical resonance of the body.
        let dark = source
            .prepare_resonance_poles_with_dynamic_load(
                spec(&load, &[0., 0.]),
                budget(),
                load_budget(),
                cx,
            )
            .unwrap();
        let expected = C64::new(-0.15, (2. - 0.15_f64.powi(2)).sqrt());
        let mut state = vec![C64::ZERO; 9];
        state[7] = C64::ONE;
        state[8] = -expected;
        let pole = dark
            .resonance_pole_from_state(expected + C64::new(0.03, 0.04), &state, &local(), cx)
            .unwrap();
        near(pole.rate_per_s, expected, 3e-11);
        assert!(pole.physical_vector[..7].iter().all(|z| z.abs() < 1e-12));
        assert!(
            pole.dynamic_load
                .unwrap()
                .flow
                .iter()
                .all(|z| z.abs() == 0.)
        );
    });
}

#[test]
fn g4_complete_state_capacity_source_band_and_cancellation_refuse_without_lossy_fallback() {
    with_cx(|cx, gate| {
        let source = solid(true, cx);
        let load = attachment(0.3, 1, vec![1., 0.]);
        let mut small = budget();
        small.maximum_coordinates = 8;
        assert!(matches!(
            source.prepare_resonance_poles_with_dynamic_load(
                spec(&load, &[1.]),
                small,
                load_budget(),
                cx
            ),
            Err(ViscoPoleError::Budget { .. })
        ));
        let mut small = load_budget();
        small.maximum_states = 1;
        assert!(matches!(
            source.prepare_resonance_poles_with_dynamic_load(
                spec(&load, &[1.]),
                budget(),
                small,
                cx
            ),
            Err(ViscoPoleError::Budget { .. })
        ));
        let mut small = load_budget();
        small.maximum_payload_bytes = 0;
        assert!(
            source
                .prepare_resonance_poles_with_dynamic_load(spec(&load, &[1.]), budget(), small, cx)
                .is_err()
        );
        assert!(
            source
                .prepare_resonance_poles_with_dynamic_load(
                    spec(&load, &[]),
                    budget(),
                    load_budget(),
                    cx
                )
                .is_err()
        );
        // The pHS storage interface can own nonlinear constitutive laws. They
        // must not silently become an invented quadratic attachment here.
        struct Nonlinear;
        impl fs_phs::Storage for Nonlinear {
            fn hamiltonian(&self, x: &[f64]) -> f64 {
                0.5 * x[0] * x[0] + 0.25 * x[0].powi(4)
            }
            fn gradient(&self, x: &[f64], out: &mut [f64]) {
                out[0] = x[0] + x[0].powi(3);
            }
        }
        let nonlinear =
            PortHamiltonian::new(1, 1, vec![0.], vec![0.2], vec![1.], Box::new(Nonlinear)).unwrap();
        assert!(matches!(
            source.prepare_resonance_poles_with_dynamic_load(
                spec(&nonlinear, &[1.]),
                budget(),
                load_budget(),
                cx
            ),
            Err(ViscoPoleError::Invalid(_))
        ));
        let mut limited = spec(&load, &[1.]);
        limited.band_rad_s = (0., 1.);
        let prepared = source
            .prepare_resonance_poles_with_dynamic_load(limited, budget(), load_budget(), cx)
            .unwrap();
        assert!(matches!(
            prepared.resonance_pole(C64::new(-1., 2.), &[C64::ONE], &local(), cx),
            Err(ViscoPoleError::SourceBand { .. })
        ));
        assert!(
            prepared
                .resonance_pole_from_state(C64::new(0., 0.5), &[C64::ONE], &local(), cx)
                .is_err()
        );
        gate.request();
        assert!(
            prepared
                .resonance_pole_from_state(C64::new(0., 0.5), &vec![C64::ONE; 9], &local(), cx)
                .is_err()
        );
        assert!(
            source
                .prepare_resonance_poles_with_dynamic_load(
                    spec(&load, &[1.]),
                    budget(),
                    load_budget(),
                    cx
                )
                .is_err()
        );
    });
}
