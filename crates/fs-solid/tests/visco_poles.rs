//! G1/G3: full spatial Prony poles versus independent analytic equations.
//! Constructed software specimens, not validation of measured material data.
use fs_exec::Cx;
use fs_la::nonlinear_eig::PoleOptions;
use fs_la::nonlinear_eig::beyn::{BeynError, BeynOptions, Ellipse};
use fs_material::visco::tensor::{Mandel6, PronyBranch, TensorProny};
use fs_math::c64::C64;
use fs_solid::linear3::TetAssemblyBudget;
use fs_solid::linear3::visco::*;

const NODES: [[f64; 3]; 4] = [[0., 0., 0.], [1., 0., 0.], [0., 1., 0.], [0., 0., 1.]];
const TET: [[usize; 4]; 1] = [[0, 1, 2, 3]];
fn with_cx(f: impl FnOnce(&Cx<'_>, &fs_exec::CancelGate)) {
    let gate = fs_exec::CancelGate::new_clock_free();
    fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default()).scope(|arena| {
        f(
            &Cx::new(
                &gate,
                arena,
                fs_exec::StreamKey {
                    seed: 413,
                    kernel_id: 91,
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
fn scalar_law() -> TensorProny {
    // M=1, K_eq=2.5, K_branch=6.5, tau=.25:
    // (1+.25s)(s²+2.5)+1.625s=.25(s+2)(s²+2s+5).
    TensorProny::new_mandel(
        diagonal(15.),
        vec![PronyBranch::new(diagonal(39.), 0.25).unwrap()],
        0.01,
    )
    .unwrap()
}
fn solid(law: TensorProny, free: &[usize], cx: &Cx<'_>) -> ViscoElements {
    let material = TetViscoMaterial {
        law,
        density_kg_m3: 60.,
        material_state_identity: [91; 32],
    };
    let fixed: Vec<_> = (0..12).filter(|i| !free.contains(i)).collect();
    ViscoElements::new(
        SolidMeshRef::Tet4 {
            nodes_m: &NODES,
            elements: &TET,
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
        maximum_evaluations: 512,
        reference_rate_per_s: 2.,
        physical_residual_tolerance: 2e-11,
    }
}
fn local() -> PoleOptions {
    PoleOptions {
        spectral_scale: 2.,
        step_tolerance: 2e-12,
        residual_tolerance: 2e-12,
        max_iterations: 32,
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
fn near(a: C64, b: C64, tol: f64) {
    assert!(
        (a - b).abs() <= tol * b.abs().max(1.),
        "actual {a:?}; expected {b:?}"
    );
}

#[test]
fn g1_original_maxwell_states_recover_analytic_complex_and_relaxation_poles() {
    with_cx(|cx, _| {
        let source = solid(scalar_law(), &[3], cx);
        let prepared = source.prepare_resonance_poles(budget(), cx).unwrap();
        assert_eq!(prepared.free_dofs(), &[3]);
        assert_eq!(prepared.dimension(), 7);
        assert_eq!(
            prepared.memory_sites(),
            &[ViscoMemorySite {
                element: 0,
                point: 0,
                branch: 0
            }]
        );
        for (seed, expected) in [
            (C64::new(-0.8, 1.8), C64::new(-1., 2.)),
            (C64::new(-0.8, -1.8), C64::new(-1., -2.)),
            (C64::from_re(-1.8), C64::from_re(-2.)),
        ] {
            let pole = prepared
                .resonance_pole(seed, &[C64::ONE], &local(), cx)
                .unwrap();
            near(pole.rate_per_s, expected, 2e-11);
            assert!(pole.mechanical_residual < 2e-11 && pole.memory_residual < 2e-11);
            assert!(pole.rate_per_s.re < 0.);
            let strain = pole.physical_vector[1] / pole.physical_vector[0];
            near(
                strain,
                expected.scale(0.25) / (C64::ONE + expected.scale(0.25)),
                2e-11,
            );
            assert!(pole.physical_vector[2..].iter().all(|z| z.abs() < 1e-12));
        }
    });
}

#[test]
fn g1_contour_discovers_actual_hereditary_resonance_with_live_quadrature_ladder() {
    with_cx(|cx, _| {
        let source = solid(scalar_law(), &[3], cx);
        let prepared = source.prepare_resonance_poles(budget(), cx).unwrap();
        let region = Ellipse {
            center: C64::new(-0.8, 2.1),
            real_radius: 0.7,
            imaginary_radius: 0.7,
        };
        let mut errors = Vec::new();
        for nodes in [32, 64, 128] {
            let found = prepared
                .resonance_spectrum(region, &contour(prepared.dimension(), nodes), cx)
                .unwrap();
            assert_eq!(found.argument_count, 1);
            assert_eq!(found.rank, 1);
            assert_eq!(found.poles.len(), 1);
            near(found.poles[0].rate_per_s, C64::new(-1., 2.), 2e-11);
            errors.push(found.moment_relative_change);
        }
        assert!(
            errors[1] < errors[0] && errors[2] < errors[0],
            "nested changes {errors:?}"
        );
    });
}

#[test]
fn g1_full_tensor_coupling_and_consistent_mass_match_independent_polynomial() {
    with_cx(|cx, _| {
        let mut eq = diagonal(12.);
        let mut branch = diagonal(18.);
        let mut viscosity = diagonal(0.6);
        eq[0][3] = 2.4;
        eq[3][0] = 2.4;
        branch[0][3] = -4.2;
        branch[3][0] = -4.2;
        viscosity[0][3] = 0.12;
        viscosity[3][0] = 0.12;
        let tau = 0.3;
        let law = TensorProny::new_mandel(eq, vec![PronyBranch::new(branch, tau).unwrap()], 0.02)
            .unwrap()
            .with_kelvin_viscosity_mandel(viscosity)
            .unwrap();
        // Free x displacement at node 1 and x displacement at node 2.
        // B columns are xx=1 and xy_Mandel=1/sqrt(2). M=[[1,.5],[.5,1]].
        let source = solid(law, &[3, 6], cx);
        let prepared = source.prepare_resonance_poles(budget(), cx).unwrap();
        let root2 = 2_f64.sqrt();
        let matrix = |a: &Mandel6| {
            [
                a[0][0] / 6.,
                a[0][3] / (6. * root2),
                a[3][0] / (6. * root2),
                a[3][3] / 12.,
            ]
        };
        let k = matrix(&eq);
        let b = matrix(&branch);
        let c = matrix(&viscosity);
        // Independent rational 2x2 determinant Newton, used only as an oracle.
        let eval = |s: C64| {
            let t = s.scale(tau);
            let d = C64::ONE + t;
            let mut a = [C64::ZERO; 4];
            let mut da = a;
            for i in 0..4 {
                let mass = if i == 0 || i == 3 { 1. } else { 0.5 };
                a[i] =
                    C64::from_re(k[i]) + s.scale(c[i]) + (s * s).scale(mass) + (t / d).scale(b[i]);
                da[i] = C64::from_re(c[i])
                    + s.scale(2. * mass)
                    + (C64::from_re(tau) / (d * d)).scale(b[i]);
            }
            (
                a[0] * a[3] - a[1] * a[2],
                da[0] * a[3] + a[0] * da[3] - da[1] * a[2] - a[1] * da[2],
                a,
            )
        };
        for seed in [C64::new(-0.2, 0.9), C64::new(-0.5, 2.2)] {
            let mut expected = seed;
            for _ in 0..30 {
                let (f, df, _) = eval(expected);
                let delta = f / df;
                expected = expected - delta;
                if delta.abs() < 1e-14 {
                    break;
                }
            }
            let (_, _, a) = eval(expected);
            let initial = [C64::ONE, -a[0] / a[1]];
            let pole = prepared
                .resonance_pole(seed, &initial, &local(), cx)
                .unwrap();
            near(pole.rate_per_s, expected, 4e-11);
            let (_, _, actual) = eval(pole.rate_per_s);
            let u = &pole.physical_vector[..2];
            assert!((actual[0] * u[0] + actual[1] * u[1]).abs() < 1e-10);
            assert!((actual[2] * u[0] + actual[3] * u[1]).abs() < 1e-10);
        }
    });
}

#[test]
fn g0_source_bands_complete_state_resources_and_cancellation_refuse() {
    with_cx(|cx, gate| {
        let source = solid(scalar_law().with_band((0.5, 4.)).unwrap(), &[3], cx);
        for b in [
            ViscoPoleBudget {
                maximum_coordinates: 6,
                ..budget()
            },
            ViscoPoleBudget {
                maximum_payload_bytes: 100,
                ..budget()
            },
            ViscoPoleBudget {
                maximum_assembly_work: 1,
                ..budget()
            },
        ] {
            assert!(matches!(
                source.prepare_resonance_poles(b, cx),
                Err(ViscoPoleError::Budget { .. })
            ));
        }
        let prepared = source.prepare_resonance_poles(budget(), cx).unwrap();
        assert!(matches!(
            prepared.resonance_pole(C64::new(-0.1, 0.2), &[C64::ONE], &local(), cx),
            Err(ViscoPoleError::SourceBand { .. })
        ));
        assert!(matches!(
            prepared.resonance_pole(C64::from_re(-4.), &[C64::ONE], &local(), cx),
            Err(ViscoPoleError::SeedAtRelaxationRate { .. })
        ));
        let outside = Ellipse {
            center: C64::new(-1., 2.),
            real_radius: 4.,
            imaginary_radius: 4.,
        };
        assert!(
            prepared
                .resonance_spectrum(outside, &contour(7, 32), cx)
                .is_err()
        );
        let short = source
            .prepare_resonance_poles(
                ViscoPoleBudget {
                    maximum_evaluations: 1,
                    ..budget()
                },
                cx,
            )
            .unwrap();
        assert!(matches!(
            short.resonance_pole(C64::new(-0.8, 1.8), &[C64::ONE], &local(), cx),
            Err(ViscoPoleError::Budget {
                quantity: "operator evaluations",
                ..
            })
        ));
        gate.request();
        assert!(
            prepared
                .resonance_pole(C64::new(-0.8, 1.8), &[C64::ONE], &local(), cx)
                .is_err()
        );
    });
}

#[test]
fn g3_lossless_poles_remain_lossless_and_oversized_contour_cannot_drop_roots() {
    with_cx(|cx, _| {
        let source = solid(
            TensorProny::new_mandel(diagonal(24.), vec![], 0.01).unwrap(),
            &[3],
            cx,
        );
        let prepared = source.prepare_resonance_poles(budget(), cx).unwrap();
        assert_eq!(prepared.dimension(), 1);
        assert!(prepared.memory_sites().is_empty());
        let pole = prepared
            .resonance_pole(C64::new(0., 1.9), &[C64::ONE], &local(), cx)
            .unwrap();
        near(pole.rate_per_s, C64::new(0., 2.), 2e-12);
        assert_eq!(pole.rate_per_s.re, 0.);
        assert_eq!(pole.memory_residual, 0.);
        let again = prepared
            .resonance_pole(C64::new(0., 1.9), &[C64::ONE], &local(), cx)
            .unwrap();
        assert_eq!(pole.rate_per_s.re.to_bits(), again.rate_per_s.re.to_bits());
        assert_eq!(pole.rate_per_s.im.to_bits(), again.rate_per_s.im.to_bits());
        let both = Ellipse {
            center: C64::ZERO,
            real_radius: 3.,
            imaginary_radius: 3.,
        };
        let refused = prepared.resonance_spectrum(both, &contour(1, 128), cx);
        match refused {
            Err(ViscoPoleError::Contour(BeynError::CountUnresolved { fine, coarse })) => {
                // The shared owner refuses an argument count larger than its
                // row/probe dimension even when both quadratures resolve it.
                near(fine, C64::from_re(2.), 1e-10);
                near(coarse, C64::from_re(2.), 1e-10);
            }
            other => panic!("two roots must not be returned as one: {other:?}"),
        }
    });
}
