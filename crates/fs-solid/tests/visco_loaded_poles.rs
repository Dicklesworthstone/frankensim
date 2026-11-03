//! G1/G3/G4: original hereditary equations with physical boundary K/C/M.
//! Software specimens and independent analytic equations, not material validation.
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
                    seed: 413,
                    kernel_id: 92,
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
fn solid(free: &[usize], cx: &Cx<'_>) -> ViscoElements {
    let material = TetViscoMaterial {
        law: TensorProny::new_mandel(
            diagonal(15.),
            vec![PronyBranch::new(diagonal(39.), 0.25).unwrap()],
            0.01,
        )
        .unwrap(),
        density_kg_m3: 60.,
        material_state_identity: [92; 32],
    };
    let fixed: Vec<_> = (0..12).filter(|i| !free.contains(i)).collect();
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
        maximum_evaluations: 512,
        reference_rate_per_s: 2.,
        physical_residual_tolerance: 2e-11,
    }
}
fn port_budget() -> ViscoPolePortBudget {
    ViscoPolePortBudget {
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
        max_iterations: 32,
        max_dimension: 100,
        max_dense_bytes: 8 * 1024 * 1024,
        max_work: 100_000_000,
    }
}
fn load() -> ViscoPolePortLoad<'static> {
    ViscoPolePortLoad {
        ports: 1,
        force_map: &[2.],
        stiffness: &[1.],
        damping: &[1.],
        inertance: &[0.25],
        band_rad_s: (0., 100.),
    }
}
fn near(a: C64, b: C64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance * b.abs().max(1.),
        "actual {a:?}; expected {b:?}"
    );
}
fn inner(a: C64, b: C64) -> f64 {
    a.re * b.re + a.im * b.im
}

#[test]
fn g1_passive_support_mass_stiffness_and_loss_recover_independent_loaded_cubic() {
    with_cx(|cx, _| {
        let source = solid(&[3], cx);
        // Bare M=1, Keq=2.5, Kbranch=6.5, tau=.25. B=2 adds
        // M=1, C=4, K=4: (1+.25s)(2s²+4s+6.5)+1.625s
        // = .5(s+2)(s²+4s+6.5). Original material coordinates stay present.
        let prepared = source
            .prepare_resonance_poles_with_port_load(load(), budget(), port_budget(), cx)
            .unwrap();
        assert_eq!(prepared.dimension(), 7);
        assert_eq!(prepared.memory_sites().len(), 1);
        for expected in [C64::new(-2., 2.5_f64.sqrt()), C64::from_re(-2.)] {
            let pole = prepared
                .resonance_pole(expected + C64::new(0.1, 0.05), &[C64::ONE], &local(), cx)
                .unwrap();
            near(pole.rate_per_s, expected, 3e-11);
            assert!(pole.mechanical_residual < 2e-11 && pole.memory_residual < 2e-11);
            let u = pole.physical_vector[0];
            near(
                pole.physical_vector[1] / u,
                expected.scale(0.25) / (C64::ONE + expected.scale(0.25)),
                3e-11,
            );
            let port = pole.port_load.unwrap();
            let d = u.scale(2.);
            let v = expected * d;
            let e = d + v + (expected * v).scale(0.25);
            near(port.displacement[0], d, 2e-11);
            near(port.flow[0], v, 3e-11);
            near(port.effort[0], e, 3e-11);
            let energy = 0.5 * (inner(d, d) + 0.25 * inner(v, v));
            let heat = inner(v, v);
            assert!((port.storage_j - energy).abs() < 3e-11 * energy);
            assert!((port.dissipation_w - heat).abs() < 3e-11 * heat);
            assert!((port.absorbed_power_w - inner(v, e)).abs() < 3e-11 * (heat + energy));
            assert!(port.signed_power_defect_w.abs() < 3e-13 * (heat + energy));
        }
        let bare = source
            .prepare_resonance_poles(budget(), cx)
            .unwrap()
            .resonance_pole(C64::new(-0.9, 1.9), &[C64::ONE], &local(), cx)
            .unwrap();
        near(bare.rate_per_s, C64::new(-1., 2.), 2e-11);
        assert!(bare.port_load.is_none());
    });
}

#[test]
fn g1_mutual_ports_preserve_full_coupled_mass_and_original_material_memory() {
    with_cx(|cx, _| {
        let source = solid(&[3, 6], cx);
        let b = [1., 0.3, -0.4, 1.2];
        let k = [3., -1., -1., 2.];
        let c = [0.5, 0.2, 0.2, 0.8];
        let m = [0.4, -0.1, -0.1, 0.7];
        let prepared = source
            .prepare_resonance_poles_with_port_load(
                ViscoPolePortLoad {
                    ports: 2,
                    force_map: &b,
                    stiffness: &k,
                    damping: &c,
                    inertance: &m,
                    band_rad_s: (0., 100.),
                },
                budget(),
                port_budget(),
                cx,
            )
            .unwrap();
        let project = |a: &[f64; 4]| {
            let mut out = [0.; 4];
            for i in 0..2 {
                for j in 0..2 {
                    for x in 0..2 {
                        for y in 0..2 {
                            out[i * 2 + j] += b[i * 2 + x] * a[x * 2 + y] * b[j * 2 + y];
                        }
                    }
                }
            }
            out
        };
        let (pk, pc, pm) = (project(&k), project(&c), project(&m));
        let eval = |s: C64| {
            let branch = [6.5, 0., 0., 3.25];
            let eq = [2.5, 0., 0., 1.25];
            let mass = [1., 0.5, 0.5, 1.];
            let t = s.scale(0.25);
            let den = C64::ONE + t;
            let a: [C64; 4] = core::array::from_fn(|i| {
                C64::from_re(eq[i] + pk[i])
                    + s.scale(pc[i])
                    + (s * s).scale(mass[i] + pm[i])
                    + (t / den).scale(branch[i])
            });
            let da: [C64; 4] = core::array::from_fn(|i| {
                C64::from_re(pc[i])
                    + s.scale(2. * (mass[i] + pm[i]))
                    + (C64::from_re(0.25) / (den * den)).scale(branch[i])
            });
            (
                a[0] * a[3] - a[1] * a[2],
                da[0] * a[3] + a[0] * da[3] - da[1] * a[2] - a[1] * da[2],
                a,
            )
        };
        let mut found = Vec::new();
        for seed in [C64::new(-0.5, 1.2), C64::new(-0.7, 2.8)] {
            let mut expected = seed;
            for _ in 0..40 {
                let (f, df, _) = eval(expected);
                let step = f / df;
                expected = expected - step;
                if step.abs() < 1e-14 {
                    break;
                }
            }
            let (_, _, a) = eval(expected);
            let pole = prepared
                .resonance_pole(seed, &[C64::ONE, -a[0] / a[1]], &local(), cx)
                .unwrap();
            near(pole.rate_per_s, expected, 4e-11);
            let (_, _, actual) = eval(pole.rate_per_s);
            let u = &pole.physical_vector[..2];
            for i in 0..2 {
                assert!((actual[2 * i] * u[0] + actual[2 * i + 1] * u[1]).abs() < 1e-10);
            }
            let port = pole.port_load.unwrap();
            for j in 0..2 {
                near(
                    port.displacement[j],
                    u[0].scale(b[j]) + u[1].scale(b[2 + j]),
                    3e-12,
                );
            }
            assert!(port.storage_j > 0. && port.dissipation_w > 0.);
            assert!(
                port.signed_power_defect_w.abs() < 5e-13 * (port.dissipation_w + port.storage_j)
            );
            found.push(pole.rate_per_s);
        }
        assert!(
            (found[0] - found[1]).abs() > 0.5,
            "distinct coupled roots {found:?}"
        );
    });
}

#[test]
fn g3_zero_boundary_preserves_original_pole_bits_and_band_is_not_extrapolated() {
    with_cx(|cx, _| {
        let source = solid(&[3], cx);
        let zero = ViscoPolePortLoad {
            stiffness: &[0.],
            damping: &[0.],
            inertance: &[0.],
            ..load()
        };
        let plain = source.prepare_resonance_poles(budget(), cx).unwrap();
        let loaded = source
            .prepare_resonance_poles_with_port_load(zero, budget(), port_budget(), cx)
            .unwrap();
        let seed = C64::new(-0.8, 1.8);
        let a = plain
            .resonance_pole(seed, &[C64::ONE], &local(), cx)
            .unwrap();
        let b = loaded
            .resonance_pole(seed, &[C64::ONE], &local(), cx)
            .unwrap();
        for (x, y) in a
            .physical_vector
            .iter()
            .zip(&b.physical_vector)
            .chain([(&a.rate_per_s, &b.rate_per_s)])
        {
            assert_eq!(
                (x.re.to_bits(), x.im.to_bits()),
                (y.re.to_bits(), y.im.to_bits())
            );
        }
        let boundary = b.port_load.unwrap();
        assert_eq!(boundary.storage_j, 0.);
        assert_eq!(boundary.dissipation_w, 0.);
        assert_eq!(boundary.absorbed_power_w, 0.);
        // A coupled semidefinite spring with a true common-motion nullspace
        // remains admitted; neither the ports nor their nullspace are removed.
        let singular = source
            .prepare_resonance_poles_with_port_load(
                ViscoPolePortLoad {
                    ports: 2,
                    force_map: &[1., 1.],
                    stiffness: &[1., -1., -1., 1.],
                    damping: &[0.; 4],
                    inertance: &[0.; 4],
                    ..load()
                },
                budget(),
                port_budget(),
                cx,
            )
            .unwrap()
            .resonance_pole(seed, &[C64::ONE], &local(), cx)
            .unwrap();
        near(singular.rate_per_s, a.rate_per_s, 1e-13);
        let port = singular.port_load.unwrap();
        assert_eq!(port.displacement.len(), 2);
        assert_eq!(port.storage_j, 0.);
        let narrow = source
            .prepare_resonance_poles_with_port_load(
                ViscoPolePortLoad {
                    band_rad_s: (0., 1.),
                    ..load()
                },
                budget(),
                port_budget(),
                cx,
            )
            .unwrap();
        assert!(matches!(
            narrow.resonance_pole(seed, &[C64::ONE], &local(), cx),
            Err(ViscoPoleError::SourceBand { .. })
        ));
    });
}

#[test]
fn g0_g4_active_malformed_overbudget_or_cancelled_loads_refuse() {
    with_cx(|cx, gate| {
        let source = solid(&[3], cx);
        for bad in [
            ViscoPolePortLoad {
                stiffness: &[-1e-100],
                ..load()
            },
            ViscoPolePortLoad {
                damping: &[-1e-100],
                ..load()
            },
            ViscoPolePortLoad {
                inertance: &[-1e-100],
                ..load()
            },
            ViscoPolePortLoad {
                force_map: &[f64::NAN],
                ..load()
            },
            ViscoPolePortLoad {
                inertance: &[],
                ..load()
            },
            ViscoPolePortLoad {
                band_rad_s: (2., 1.),
                ..load()
            },
            ViscoPolePortLoad {
                ports: 2,
                force_map: &[1., 1.],
                stiffness: &[1., 2., 2., 1.],
                damping: &[0.; 4],
                inertance: &[0.; 4],
                ..load()
            },
            // A zero diagonal cannot conceal an active cross coupling.
            ViscoPolePortLoad {
                ports: 2,
                force_map: &[1., 1.],
                stiffness: &[0., 1e-100, 1e-100, 1.],
                damping: &[0.; 4],
                inertance: &[0.; 4],
                ..load()
            },
            ViscoPolePortLoad {
                ports: 2,
                force_map: &[1., 1.],
                stiffness: &[1., 0., 0.1, 1.],
                damping: &[0.; 4],
                inertance: &[0.; 4],
                ..load()
            },
        ] {
            assert!(
                source
                    .prepare_resonance_poles_with_port_load(bad, budget(), port_budget(), cx)
                    .is_err()
            );
        }
        // The documented extra payload/work envelopes for n=p=1 are 360
        // bytes and 2940 scalar-work units. Exact caps admit, one below refuses.
        let exact = ViscoPolePortBudget {
            maximum_ports: 1,
            maximum_payload_bytes: 360,
            maximum_preparation_work: 2940,
        };
        assert!(
            source
                .prepare_resonance_poles_with_port_load(load(), budget(), exact, cx)
                .is_ok()
        );
        for cap in [
            ViscoPolePortBudget {
                maximum_ports: 0,
                ..exact
            },
            ViscoPolePortBudget {
                maximum_payload_bytes: 359,
                ..exact
            },
            ViscoPolePortBudget {
                maximum_preparation_work: 2939,
                ..exact
            },
        ] {
            assert!(matches!(
                source.prepare_resonance_poles_with_port_load(load(), budget(), cap, cx),
                Err(ViscoPoleError::Budget { .. })
            ));
        }
        let prepared = source
            .prepare_resonance_poles_with_port_load(load(), budget(), port_budget(), cx)
            .unwrap();
        gate.request();
        assert!(
            source
                .prepare_resonance_poles_with_port_load(load(), budget(), port_budget(), cx)
                .is_err()
        );
        assert!(
            prepared
                .resonance_pole(C64::new(-2., 1.6), &[C64::ONE], &local(), cx)
                .is_err()
        );
    });
}
