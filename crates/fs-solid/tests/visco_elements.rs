//! IRP M1.3-W2a: synthetic constitutive fixtures, not specimen validation.
use fs_exec::Cx;
use fs_material::visco::ViscoError;
use fs_material::visco::tensor::{Mandel6, PronyBranch, TensorProny};
use fs_math::det;
use fs_solid::linear3::tet10::{TetQuadrature, promote_p1_to_p2};
use fs_solid::linear3::visco::*;
use fs_solid::linear3::{TetAssemblyBudget, TetElasticMaterial};

const Q: TetQuadrature = TetQuadrature {
    stiffness_points_per_dir: 3,
    mass_points_per_dir: 4,
};

fn with_cx<R>(f: impl FnOnce(&Cx<'_>, &fs_exec::CancelGate) -> R) -> R {
    let gate = fs_exec::CancelGate::new_clock_free();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            fs_exec::StreamKey {
                seed: 60,
                kernel_id: 60,
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
fn fixture(n: usize) -> TetViscoMaterial {
    let eq = diagonal([11e9, 0.9e9, 0.5e9, 1.4e9, 0.12e9, 1.3e9]);
    let branch = PronyBranch::new(diagonal([2e8, 3e7, 2e7, 4e7, 1e7, 3e7]), 1e-3).unwrap();
    TetViscoMaterial {
        law: TensorProny::new(eq, None, vec![branch; n], (0.0, 1e12)).unwrap(),
        density_kg_m3: 450.0,
        material_state_identity: [60; 32],
    }
}
fn nodes() -> Vec<[f64; 3]> {
    vec![[0.0; 3], [0.1, 0.0, 0.0], [0.0, 0.1, 0.0], [0.0, 0.0, 0.1]]
}
fn displacement(nodes: &[[f64; 3]], strain: [f64; 6]) -> Vec<f64> {
    let s = core::f64::consts::FRAC_1_SQRT_2;
    nodes
        .iter()
        .flat_map(|p| {
            [
                strain[0] * p[0] + s * (strain[3] * p[1] + strain[5] * p[2]),
                strain[1] * p[1] + s * (strain[3] * p[0] + strain[4] * p[2]),
                strain[2] * p[2] + s * (strain[5] * p[0] + strain[4] * p[1]),
            ]
        })
        .collect()
}
fn quadratic(a: &fs_sparse::Csr, u: &[f64]) -> f64 {
    let mut y = vec![0.0; u.len()];
    a.spmv(u, &mut y);
    u.iter().zip(y).map(|(x, y)| x * y).sum()
}
fn q6(a: &Mandel6, u: &[f64; 6]) -> f64 {
    (0..6)
        .map(|i| u[i] * (0..6).map(|j| a[i][j] * u[j]).sum::<f64>())
        .sum()
}
fn close(a: f64, b: f64, rel: f64) {
    assert!(
        (a - b).abs() <= rel * a.abs().max(b.abs()).max(1e-25),
        "{a:e} != {b:e}, relative={:e}",
        (a - b).abs() / a.abs().max(b.abs())
    );
}
fn system(mesh: SolidMeshRef<'_>, m: &TetViscoMaterial, cx: &Cx<'_>) -> ViscoElements {
    ViscoElements::new(
        mesh,
        TetViscoField::Uniform(m),
        &[],
        TetAssemblyBudget::standard(),
        ViscoStateBudget {
            max_internal_variables: 100_000,
        },
        cx,
    )
    .unwrap()
}

#[test]
fn g0_uniform_strain_complex_energy_exact() {
    with_cx(|cx, _| {
        let v = nodes();
        let t = [[0, 1, 2, 3]];
        let (p, t10) = promote_p1_to_p2(&v, &t, None, cx).unwrap();
        let m = fixture(1);
        let strain = [2e-5, -1e-5, 3e-5, 4e-5, -2e-5, 1e-5];
        for (mesh, xyz) in [
            (
                SolidMeshRef::Tet4 {
                    nodes_m: &v,
                    elements: &t,
                },
                &v,
            ),
            (
                SolidMeshRef::Tet10 {
                    nodes_m: &p,
                    elements: &t10,
                    quadrature: Q,
                },
                &p,
            ),
        ] {
            let s = system(mesh, &m, cx);
            let u = displacement(xyz, strain);
            for w in [0.0, 1000.0, 1e12] {
                let k = s.complex_stiffness(w, cx).unwrap();
                let (r, i) = m.law.complex_modulus(w).unwrap();
                close(quadratic(&k.real, &u), q6(&r, &strain) * 0.001 / 6.0, 1e-13);
                if w > 0.0 {
                    close(quadratic(&k.imag, &u), q6(&i, &strain) * 0.001 / 6.0, 1e-13);
                }
                if w == 0.0 {
                    assert_eq!(k.real, s.equilibrium().stiffness);
                }
                if w == 1e12 {
                    let instant = core::array::from_fn(|i| {
                        core::array::from_fn(|j| {
                            m.law.c_eq()[i][j] + m.law.branches()[0].c_pa[i][j]
                        })
                    });
                    close(
                        quadratic(&k.real, &u),
                        q6(&instant, &strain) * 0.001 / 6.0,
                        1e-13,
                    );
                }
            }
        }
    });
}

#[test]
fn g0_empty_branches_bitwise_equal_elastic() {
    with_cx(|cx, _| {
        let v = nodes();
        let t = [[0, 1, 2, 3]];
        let (p, t10) = promote_p1_to_p2(&v, &t, None, cx).unwrap();
        let m = fixture(0);
        for mesh in [
            SolidMeshRef::Tet4 {
                nodes_m: &v,
                elements: &t,
            },
            SolidMeshRef::Tet10 {
                nodes_m: &p,
                elements: &t10,
                quadrature: Q,
            },
        ] {
            let s = system(mesh, &m, cx);
            let k = s.complex_stiffness(1000.0, cx).unwrap();
            assert_eq!(k.real, s.equilibrium().stiffness);
            assert!(k.imag.row(0).1.iter().all(|v| *v == 0.0));
        }
    });
}

#[test]
fn g0_passivity_under_random_multiaxial_histories() {
    with_cx(|cx, _| {
        let v = nodes();
        let t = [[0, 1, 2, 3]];
        let (p, t10) = promote_p1_to_p2(&v, &t, None, cx).unwrap();
        let m = fixture(3);
        for (mesh, n) in [
            (
                SolidMeshRef::Tet4 {
                    nodes_m: &v,
                    elements: &t,
                },
                v.len(),
            ),
            (
                SolidMeshRef::Tet10 {
                    nodes_m: &p,
                    elements: &t10,
                    quadrature: Q,
                },
                p.len(),
            ),
        ] {
            let s = system(mesh, &m, cx);
            let mut state = s.state(cx).unwrap();
            let mut rng = fs_rand::StreamKey {
                seed: 0x4d31_3357_3261,
                kernel: 60,
                tile: n as u32,
            }
            .stream();
            for _ in 0..4096 {
                let u: Vec<f64> = (0..3 * n).map(|_| 2e-7 * (rng.next_f64() - 0.5)).collect();
                let out = state
                    .step_internal_force(&u, 1e-5 + 1e-3 * rng.next_f64(), cx)
                    .unwrap();
                assert!(out.step_dissipated_j >= 0.0);
                close(
                    out.energy.work_j,
                    out.energy.stored_j + out.energy.dissipated_j,
                    1e-12,
                );
            }
        }
    });
}

#[test]
fn g1_time_domain_element_force_matches_complex_stiffness() {
    with_cx(|cx, _| {
        let v = nodes();
        let t = [[0, 1, 2, 3]];
        let m = fixture(1);
        let s = system(
            SolidMeshRef::Tet4 {
                nodes_m: &v,
                elements: &t,
            },
            &m,
            cx,
        );
        let w = 1000.0;
        let u = displacement(&v, [1e-5, -2e-5, 3e-5, 1e-5, 2e-5, -1e-5]);
        let k = s.complex_stiffness(w, cx).unwrap();
        let mut kr = vec![0.0; 12];
        let mut ki = kr.clone();
        k.real.spmv(&u, &mut kr);
        k.imag.spmv(&u, &mut ki);
        let mut errors = Vec::new();
        for n in [32, 64, 128] {
            let mut state = s.state(cx).unwrap();
            let dt = core::f64::consts::TAU / (w * f64::from(n));
            let mut err = 0.0;
            let mut norm = 0.0;
            for j in 1..=20 * n {
                let phase = w * dt * f64::from(j);
                let (a, b) = (det::sin(phase), det::cos(phase));
                let out = state
                    .step_internal_force(&u.iter().map(|v| v * a).collect::<Vec<_>>(), dt, cx)
                    .unwrap();
                if j > 19 * n {
                    for i in 0..12 {
                        let exact = kr[i] * a + ki[i] * b;
                        let difference = out.full_force_n[i] - exact;
                        err += difference * difference;
                        norm += exact * exact;
                    }
                }
            }
            errors.push(det::sqrt(err / norm));
        }
        for pair in errors.windows(2) {
            let order = det::ln(pair[0] / pair[1]) / det::ln(2.0);
            assert!(
                (order - 2.0).abs() < 0.2,
                "errors={errors:?}, order={order}"
            );
        }
        println!("G1 harmonic endpoint-force dt ladder [32,64,128]: {errors:?}");
    });
}

#[test]
fn g3_rotation_metamorphic() {
    with_cx(|cx, _| {
        let v = nodes();
        let t = [[0, 1, 2, 3]];
        let m = fixture(2);
        let r = fs_material::tensor::rotation([0.0, 0.0, 1.0], 0.37);
        let rot = |x: [f64; 3]| core::array::from_fn(|i| (0..3).map(|j| r[i][j] * x[j]).sum());
        let vr: Vec<_> = v.iter().copied().map(rot).collect();
        let mr = TetViscoMaterial {
            law: m.law.rotated(r).unwrap(),
            ..m.clone()
        };
        let a = system(
            SolidMeshRef::Tet4 {
                nodes_m: &v,
                elements: &t,
            },
            &m,
            cx,
        );
        let b = system(
            SolidMeshRef::Tet4 {
                nodes_m: &vr,
                elements: &t,
            },
            &mr,
            cx,
        );
        let u = displacement(&v, [1e-5, 2e-5, -1e-5, 3e-5, -2e-5, 4e-5]);
        let ur: Vec<f64> = u
            .chunks_exact(3)
            .flat_map(|p| rot([p[0], p[1], p[2]]))
            .collect();
        let ka = a.complex_stiffness(1000.0, cx).unwrap();
        let kb = b.complex_stiffness(1000.0, cx).unwrap();
        close(quadratic(&ka.real, &u), quadratic(&kb.real, &ur), 1e-12);
        close(
            modal_loss_energy_quotient(&ka, &u).unwrap(),
            modal_loss_energy_quotient(&kb, &ur).unwrap(),
            1e-12,
        );
    });
}

#[test]
fn g0_refuses_frequency_outside_law_band() {
    with_cx(|cx, _| {
        let v = nodes();
        let t = [[0, 1, 2, 3]];
        let mut m = fixture(1);
        m.law = m.law.with_band((100.0, 2000.0)).unwrap();
        let s = system(
            SolidMeshRef::Tet4 {
                nodes_m: &v,
                elements: &t,
            },
            &m,
            cx,
        );
        for w in [0.0, 99.0, 2001.0] {
            assert!(matches!(
                s.complex_stiffness(w, cx),
                Err(ViscoElementError::Material(ViscoError::OutOfBand { .. }))
            ));
        }
        assert!(s.complex_stiffness(f64::NAN, cx).is_err());
        assert!(
            modal_loss_energy_quotient(&s.complex_stiffness(1000.0, cx).unwrap(), &[0.0; 12])
                .is_err()
        );
    });
}

#[test]
fn g0_internal_state_budget_cap_and_cap_plus_one() {
    with_cx(|cx, _| {
        let v = nodes();
        let t = [[0, 1, 2, 3]];
        let (p, t10) = promote_p1_to_p2(&v, &t, None, cx).unwrap();
        let m = fixture(64);
        for (mesh, needed) in [
            (
                SolidMeshRef::Tet4 {
                    nodes_m: &v,
                    elements: &t,
                },
                384,
            ),
            (
                SolidMeshRef::Tet10 {
                    nodes_m: &p,
                    elements: &t10,
                    quadrature: Q,
                },
                27 * 384,
            ),
        ] {
            let make = |cap| {
                ViscoElements::new(
                    mesh,
                    TetViscoField::Uniform(&m),
                    &[],
                    TetAssemblyBudget::standard(),
                    ViscoStateBudget {
                        max_internal_variables: cap,
                    },
                    cx,
                )
            };
            assert_eq!(make(needed).unwrap().internal_variables(), needed);
            assert!(
                matches!(make(needed-1),Err(ViscoElementError::InternalStateBudget {needed:n,cap}) if n==needed && cap==needed-1)
            );
        }
        assert!(
            TensorProny::new(
                *m.law.c_eq(),
                None,
                vec![m.law.branches()[0]; 65],
                (0.0, 1e12)
            )
            .is_err()
        );
    });
}

#[test]
fn g0_invalid_history_and_cancellation_preserve_state() {
    with_cx(|cx, gate| {
        let v = nodes();
        let t = [[0, 1, 2, 3]];
        let m = fixture(1);
        let s = system(
            SolidMeshRef::Tet4 {
                nodes_m: &v,
                elements: &t,
            },
            &m,
            cx,
        );
        let mut state = s.state(cx).unwrap();
        state.step_internal_force(&[1e-7; 12], 1e-4, cx).unwrap();
        let energy = state.energy();
        for (u, dt) in [
            (vec![0.0; 11], 1e-3),
            (vec![f64::NAN; 12], 1e-3),
            (vec![0.0; 12], 0.0),
        ] {
            assert!(state.step_internal_force(&u, dt, cx).is_err());
            assert_eq!(state.energy(), energy);
        }
        // Overflows partway through a candidate; the original history survives.
        assert!(
            state
                .step_internal_force(&[f64::MAX; 12], 1e-6, cx)
                .is_err()
        );
        assert_eq!(state.energy(), energy);
        gate.request();
        assert!(state.step_internal_force(&[0.0; 12], 1e-3, cx).is_err());
        assert_eq!(state.energy(), energy);
    });
}

#[test]
fn g0_isotropic_zener_uniform_strain_loss() {
    with_cx(|cx, _| {
        // Uniaxial constrained strain: the independent P-wave Zener modulus is
        // lambda+2mu; this is not a free-side Poisson-relaxation beam oracle.
        let elastic = TetElasticMaterial::try_new(1200.0, 3e9, 0.25, [60; 32]).unwrap();
        let c = *elastic.stiffness_mandel_pa();
        let fraction = 0.2;
        let tau = 0.002;
        let branch = PronyBranch::new(c.map(|row| row.map(|v| fraction * v)), tau).unwrap();
        let m = TetViscoMaterial {
            law: TensorProny::new(c, None, vec![branch], (0.0, 1e12)).unwrap(),
            density_kg_m3: 1200.0,
            material_state_identity: [60; 32],
        };
        let v = nodes();
        let t = [[0, 1, 2, 3]];
        let s = system(
            SolidMeshRef::Tet4 {
                nodes_m: &v,
                elements: &t,
            },
            &m,
            cx,
        );
        let u = displacement(&v, [1e-5, 0.0, 0.0, 0.0, 0.0, 0.0]);
        for w in [100.0, 500.0, 2000.0] {
            let x = w * tau;
            let k = s.complex_stiffness(w, cx).unwrap();
            close(
                modal_loss_energy_quotient(&k, &u).unwrap(),
                fraction * x / (1.0 + (1.0 + fraction) * x * x),
                1e-13,
            );
        }
    });
}

#[test]
fn g0_heterogeneous_material_force_energy_and_support_reactions() {
    with_cx(|cx, _| {
        let mut xyz = nodes();
        xyz.extend(nodes().into_iter().map(|p| [p[0] + 0.2, p[1], p[2]]));
        let tets = [[0, 1, 2, 3], [4, 5, 6, 7]];
        let materials = [fixture(1), fixture(3)];
        let s = ViscoElements::new(
            SolidMeshRef::Tet4 {
                nodes_m: &xyz,
                elements: &tets,
            },
            TetViscoField::PerElement(&materials),
            &[0, 1, 2],
            TetAssemblyBudget::standard(),
            ViscoStateBudget {
                max_internal_variables: 24,
            },
            cx,
        )
        .unwrap();
        let strains = [
            [2e-5, -1e-5, 3e-5, 4e-5, -2e-5, 1e-5],
            [-1e-5, 3e-5, -2e-5, 2e-5, 1e-5, -3e-5],
        ];
        let mut u = displacement(&xyz[..4], strains[0]);
        u.extend(displacement(&xyz[4..], strains[1]));
        let mut state = s.state(cx).unwrap();
        let out = state.step_internal_force(&u, 0.0003, cx).unwrap();
        let mut expected_storage = 0.0;
        for e in 0..2 {
            let mut point = materials[e].law.state();
            let stress = materials[e]
                .law
                .step(&mut point, strains[e], 0.0003)
                .unwrap();
            expected_storage += point.stored * 0.001 / 6.0;
            // For this axis-aligned tetrahedron, the positive x/y/z corner
            // forces are its stress-tensor columns times V/h.
            let q = core::f64::consts::FRAC_1_SQRT_2;
            let columns = [
                [stress[0], q * stress[3], q * stress[5]],
                [q * stress[3], stress[1], q * stress[4]],
                [q * stress[5], q * stress[4], stress[2]],
            ];
            for a in 0..3 {
                for c in 0..3 {
                    close(
                        out.full_force_n[12 * e + 3 * (a + 1) + c],
                        columns[a][c] * 0.001 / 0.6,
                        1e-13,
                    );
                }
            }
            for c in 0..3 {
                close(
                    out.full_force_n[12 * e + c],
                    -columns.iter().map(|v| v[c]).sum::<f64>() * 0.001 / 0.6,
                    1e-13,
                );
            }
        }
        close(out.energy.stored_j, expected_storage, 1e-13);
        assert_eq!(out.reduced_force_n, out.full_force_n[3..]);
        // The first element can finish before a later element refuses. Verify
        // internal strains as well as the published ledger survive rollback.
        let mut control = state.clone();
        let mut invalid = u.clone();
        invalid[12..].fill(f64::MAX);
        assert!(state.step_internal_force(&invalid, 0.0002, cx).is_err());
        let next: Vec<_> = u.iter().map(|v| 0.7 * v).collect();
        let actual = state.step_internal_force(&next, 0.0002, cx).unwrap();
        let expected = control.step_internal_force(&next, 0.0002, cx).unwrap();
        assert_eq!(actual.full_force_n, expected.full_force_n);
        assert_eq!(actual.energy, expected.energy);
        let k = s.complex_stiffness(1000.0, cx).unwrap();
        assert_eq!(k.real.nrows(), 21);
        assert_eq!(s.equilibrium().free_dofs, (3..24).collect::<Vec<_>>());
        assert!(
            ViscoElements::new(
                SolidMeshRef::Tet4 {
                    nodes_m: &xyz,
                    elements: &tets
                },
                TetViscoField::PerElement(&materials[..1]),
                &[],
                TetAssemblyBudget::standard(),
                ViscoStateBudget {
                    max_internal_variables: 24
                },
                cx,
            )
            .is_err()
        );
    });
}

#[test]
fn g1_dispersive_zener_mode_matches_closed_form() {
    with_cx(|cx, _| {
        let v = nodes();
        let t = [[0, 1, 2, 3]];
        let (p, t10) = promote_p1_to_p2(&v, &t, None, cx).unwrap();
        for (mesh, n) in [
            (
                SolidMeshRef::Tet4 {
                    nodes_m: &v,
                    elements: &t,
                },
                12,
            ),
            (
                SolidMeshRef::Tet10 {
                    nodes_m: &p,
                    elements: &t10,
                    quadrature: Q,
                },
                30,
            ),
        ] {
            // One retained displacement coordinate yields an independent
            // analytic quadratic for s=omega². Both element orders use their
            // own consistent mass. This is a discrete Zener fixture, not a
            // continuum bar or an empirical wood validation.
            let fixed: Vec<_> = (0..n).filter(|&i| i != 3).collect();
            let elastic = fixture(0);
            let make = |m: &TetViscoMaterial| {
                ViscoElements::new(
                    mesh,
                    TetViscoField::Uniform(m),
                    &fixed,
                    TetAssemblyBudget::standard(),
                    ViscoStateBudget {
                        max_internal_variables: 1000,
                    },
                    cx,
                )
                .unwrap()
            };
            let base = make(&elastic);
            let lambda0 =
                base.equilibrium().stiffness.row(0).1[0] / base.equilibrium().mass.row(0).1[0];
            let omega0 = det::sqrt(lambda0);
            let alpha = 4.0;
            let tau = 1.0 / omega0;
            let c = *elastic.law.c_eq();
            let branch = PronyBranch::new(c.map(|r| r.map(|v| alpha * v)), tau).unwrap();
            let m = TetViscoMaterial {
                law: TensorProny::new(c, None, vec![branch], (0.0, 1e12)).unwrap(),
                ..elastic
            };
            let s = make(&m);
            let b = ViscoModeBudget {
                modes: 1,
                max_free_dofs: 1,
                max_iterations: 64,
                frequency_bracket_rad_s: (omega0 / 10.0, omega0 * 10.0),
                relative_frequency_tolerance: 1e-10,
            };
            let modes = s.modes(b, cx).unwrap();
            let mode = &modes[0];
            let expected = omega0 * det::sqrt(0.5 * (alpha + det::sqrt(alpha * alpha + 4.0)));
            close(mode.omega_rad_s, expected, 2e-10);
            let x = expected * tau;
            close(
                mode.loss_factor,
                alpha * x / (1.0 + (1.0 + alpha) * x * x),
                3e-10,
            );
            assert!(mode.relative_frequency_residual <= 4e-10);
            assert!(
                mode.omega_rad_s > 2.0 * omega0,
                "equilibrium-only modes must fail this fixture"
            );
            assert!(matches!(
                s.modes(
                    ViscoModeBudget {
                        max_iterations: 1,
                        ..b
                    },
                    cx
                ),
                Err(ViscoElementError::ModalConvergence { mode: 0 })
            ));
            assert!(matches!(
                s.modes(
                    ViscoModeBudget {
                        frequency_bracket_rad_s: (omega0 / 10.0, omega0 / 2.0),
                        ..b
                    },
                    cx
                ),
                Err(ViscoElementError::UnbracketedMode { mode: 0 })
            ));
            assert!(matches!(
                s.modes(
                    ViscoModeBudget {
                        max_free_dofs: 0,
                        ..b
                    },
                    cx
                ),
                Err(ViscoElementError::InvalidModalBudget)
            ));
        }
    });
}
