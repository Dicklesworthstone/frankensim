//! Physical CST/DKT forced motion with original tensor dispersion and ports.

use fs_exec::{CancelGate, Cx};
use fs_material::visco::tensor::{Mandel6, PronyBranch, TensorProny};
use fs_math::c64::C64;
use fs_modal::harmonic::{
    HarmonicParameterControls, HarmonicPortControls, HarmonicPortParameterDirection,
    HarmonicResponseControls,
};
use fs_plate::shell::viscosity::ShellViscosityBudget;
use fs_plate::shell::viscosity::harmonic::{
    TensorPronyShell, TensorPronyShellSection, TensorShellHarmonicError,
};
use fs_plate::shell::{CstDktMass, ShellMesh, ShellSupport, assemble_shell_with_mass};
use fs_plate::{PlateError, PlateSection, dkt_curvature_quadrature};

fn with_cx(run: impl FnOnce(&Cx<'_>, &CancelGate)) {
    let gate = CancelGate::new_clock_free();
    fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default()).scope(|arena| {
        run(
            &Cx::new(
                &gate,
                arena,
                fs_exec::StreamKey {
                    seed: 53,
                    kernel_id: 23,
                    tile: 0,
                    iteration: 0,
                },
                fs_exec::Budget::INFINITE,
                fs_exec::ExecMode::Deterministic,
            ),
            &gate,
        );
    });
}
fn mesh() -> ShellMesh {
    ShellMesh::new(
        vec![[0.0, 0.0, 0.0], [0.4, 0.0, 0.0], [0.0, 0.3, 0.0]],
        vec![[0, 1, 2]],
    )
    .unwrap()
}
fn budget() -> ShellViscosityBudget {
    ShellViscosityBudget {
        max_elements: 4,
        max_full_dofs: 36,
        max_dofs: 36,
        max_modes: 18,
        max_work: 20_000_000,
    }
}
fn controls() -> HarmonicResponseControls {
    HarmonicResponseControls {
        maximum_coordinates: 36,
        maximum_payload_bytes: 1 << 22,
        maximum_relative_residual: 2e-11,
    }
}
fn tangent() -> HarmonicParameterControls {
    HarmonicParameterControls {
        maximum_additional_payload_bytes: 1 << 20,
        maximum_work: 10_000_000,
    }
}
fn ports() -> HarmonicPortControls {
    HarmonicPortControls {
        maximum_ports: 2,
        maximum_projection_terms: 1_000_000,
        maximum_payload_bytes: 1 << 20,
    }
}
fn isotropic() -> Mandel6 {
    core::array::from_fn(|i| {
        core::array::from_fn(|j| {
            if i < 3 && j < 3 {
                4.0 + if i == j { 6.0 } else { 0.0 }
            } else if i == j {
                6.0
            } else {
                0.0
            }
        })
    })
}
fn scaled(c: Mandel6, scale: f64) -> Mandel6 {
    c.map(|row| row.map(|v| v * scale))
}
fn proportional(limit: f64) -> TensorProny {
    let c = isotropic();
    TensorProny::new_mandel(
        c,
        vec![
            PronyBranch::new(scaled(c, 0.8), 0.2).unwrap(),
            PronyBranch::new(scaled(c, 0.4), 2.0).unwrap(),
        ],
        limit,
    )
    .unwrap()
    .with_kelvin_viscosity_mandel(scaled(c, 0.03))
    .unwrap()
    .with_band((0.0, 100.0))
    .unwrap()
}
fn source(law: TensorProny, mass: CstDktMass, cx: &Cx<'_>) -> TensorPronyShell {
    TensorPronyShell::assemble(
        mesh(),
        vec![TensorPronyShellSection::new(law, 0.1, 4.0).unwrap()],
        &[0, 1],
        ShellSupport::Clamped,
        mass,
        0.0,
        budget(),
        cx,
    )
    .unwrap()
}
fn scalar(omega: f64) -> (C64, C64) {
    let mut value = C64::new(1.0, 0.03 * omega);
    let mut derivative = C64::new(0.0, 0.03);
    for (strength, tau) in [(0.8, 0.2), (0.4, 2.0)] {
        let x = omega * tau;
        let d = 1.0 + x * x;
        value = value + C64::new(strength * x * x / d, strength * x / d);
        derivative = derivative
            + C64::new(
                2.0 * strength * tau * x / (d * d),
                strength * tau * (1.0 - x * x) / (d * d),
            );
    }
    (value, derivative)
}
fn multiply(matrix: &fs_sparse::Csr, x: &[C64]) -> Vec<C64> {
    (0..x.len())
        .map(|i| (0..x.len()).fold(C64::ZERO, |sum, j| sum + x[j].scale(matrix.get(i, j))))
        .collect()
}
fn near(actual: C64, expected: C64, tolerance: f64) {
    assert!(
        (actual - expected).abs() < tolerance * expected.abs().max(1e-12),
        "{actual:?} != {expected:?}"
    );
}
fn real_near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() < tolerance * expected.abs().max(1e-15),
        "{actual:e} != {expected:e}"
    );
}
fn inner(x: &[C64], y: &[C64]) -> f64 {
    x.iter().zip(y).map(|(&a, &b)| (a.conj() * b).re).sum()
}

#[test]
fn manufactured_bending_has_analytic_motion_material_heat_and_frequency_slopes() {
    with_cx(|cx, _| {
        for mass in [CstDktMass::Lumped, CstDktMass::Consistent] {
            let source = source(proportional(0.02), mass, cx);
            let original = assemble_shell_with_mass(
                &mesh(),
                &PlateSection::isotropic(54.0 / 7.0, 2.0 / 7.0, 0.1, 4.0).unwrap(),
                &[0, 1],
                ShellSupport::Clamped,
                mass,
            )
            .unwrap();
            for i in 0..6 {
                for j in 0..6 {
                    real_near(source.model().m.get(i, j), original.m.get(i, j), 1e-12);
                }
            }
            let q0 = [
                C64::ZERO,
                C64::ZERO,
                C64::new(2e-5, -3e-6),
                C64::new(-1e-5, 2e-6),
                C64::new(3e-5, 1e-6),
                C64::ZERO,
            ];
            let dq = [
                C64::ZERO,
                C64::ZERO,
                C64::new(1e-6, -4e-7),
                C64::new(-3e-7, 2e-7),
                C64::new(2e-7, 8e-7),
                C64::ZERO,
            ];
            for omega in [0.3, 0.5, 2.0, 5.0, 12.0] {
                let q: Vec<_> = q0
                    .iter()
                    .zip(dq)
                    .map(|(&a, b)| a + b.scale(omega))
                    .collect();
                let (s, ds) = scalar(omega);
                let kq = multiply(&original.k, &q);
                let mq = multiply(&original.m, &q);
                let kdq = multiply(&original.k, &dq);
                let mdq = multiply(&original.m, &dq);
                // This pure bending subspace has exactly zero drilling strain.
                let f: Vec<_> = (0..6)
                    .map(|i| s.conj() * kq[i] - mq[i].scale(omega * omega))
                    .collect();
                let df: Vec<_> = (0..6)
                    .map(|i| {
                        ds.conj() * kq[i] + s.conj() * kdq[i]
                            - mq[i].scale(2.0 * omega)
                            - mdq[i].scale(omega * omega)
                    })
                    .collect();
                let prepared = source
                    .prepare_harmonic_response(omega, controls(), cx)
                    .unwrap()
                    .with_frequency_derivative(1 << 20, cx)
                    .unwrap();
                let result = prepared.solve(&f, Some(&df), None, tangent(), cx).unwrap();
                for i in 0..6 {
                    near(result.response.displacement[i], q[i], 2e-8);
                    near(result.displacement_derivative[i], dq[i], 2e-7);
                    near(
                        result.velocity_derivative[i],
                        C64::new(0.0, -1.0) * q[i] + C64::new(0.0, -omega) * dq[i],
                        2e-8,
                    );
                }
                let energy = inner(&q, &kq);
                let heat = 0.5 * omega * s.im * energy;
                let dheat = 0.5 * (s.im + omega * ds.im) * energy + omega * s.im * inner(&dq, &kq);
                real_near(result.response.average_viscous_power_w, heat, 1e-8);
                real_near(result.average_viscous_power_derivative, dheat, 1e-8);
                assert!(result.response.power_balance_residual_w.abs() < 1e-9 * heat);
                assert!(result.power_balance_residual_derivative.abs() < 1e-8 * dheat.abs());
                assert_eq!(
                    prepared.base().solve(&f, cx).unwrap().displacement,
                    result.response.displacement
                );
            }
        }
    });
}

fn anisotropic() -> TensorProny {
    let outer = |diagonal: f64, v: [f64; 6]| {
        core::array::from_fn(|i| {
            core::array::from_fn(|j| v[i] * v[j] + if i == j { diagonal } else { 0.0 })
        })
    };
    TensorProny::new_mandel(
        outer(5.0, [0.4, -0.5, 0.7, 0.2, -0.3, 0.6]),
        vec![
            PronyBranch::new(outer(0.6, [1.0, -2.0, 0.5, 0.25, -0.75, 1.5]), 0.125).unwrap(),
            PronyBranch::new(outer(0.9, [-0.5, 0.75, 1.0, 1.25, 0.5, -1.0]), 2.0).unwrap(),
        ],
        0.05,
    )
    .unwrap()
    .with_band((0.0, 100.0))
    .unwrap()
}

#[test]
fn full_anisotropic_dispersion_and_mutual_load_slopes_match_actual_frequency_changes() {
    with_cx(|cx, _| {
        let source = source(anisotropic(), CstDktMass::Consistent, cx);
        let db: Vec<_> = (0..12).map(|i| 0.0001 * ((i % 3) as f64 - 1.0)).collect();
        let dz = [
            C64::new(0.002, 0.001),
            C64::new(-0.0001, 0.0002),
            C64::new(-0.0001, 0.0002),
            C64::new(0.003, -0.001),
        ];
        let df: Vec<_> = (0..6)
            .map(|i| C64::new(1e-8 * (i + 1) as f64, -2e-8))
            .collect();
        let input = |w: f64| {
            let b: Vec<_> = (0..12)
                .map(|i| 0.02 * ((i % 4) as f64 - 1.5) + w * db[i])
                .collect();
            let z = [
                C64::new(0.2, 0.01) + dz[0].scale(w),
                C64::new(0.03, 0.005) + dz[1].scale(w),
                C64::new(0.03, 0.005) + dz[2].scale(w),
                C64::new(0.3, -0.02) + dz[3].scale(w),
            ];
            let f: Vec<_> = (0..6)
                .map(|i| C64::new(2e-7 * (i + 1) as f64, -1e-7) + df[i].scale(w))
                .collect();
            (b, z, f)
        };
        for omega in [0.5, 2.0, 8.0] {
            let (b, z, f) = input(omega);
            let prepared = source
                .prepare_harmonic_response_with_port_impedance(
                    omega,
                    &b,
                    2,
                    &z,
                    controls(),
                    ports(),
                    cx,
                )
                .unwrap()
                .with_frequency_derivative(1 << 20, cx)
                .unwrap();
            let result = prepared
                .solve(
                    &f,
                    Some(&df),
                    Some(HarmonicPortParameterDirection {
                        force_map: &db,
                        impedance: &dz,
                    }),
                    tangent(),
                    cx,
                )
                .unwrap();
            let h = 1e-5 * omega;
            let sample = |w| {
                let (b, z, f) = input(w);
                source
                    .prepare_harmonic_response_with_port_impedance(
                        w,
                        &b,
                        2,
                        &z,
                        controls(),
                        ports(),
                        cx,
                    )
                    .unwrap()
                    .solve(&f, cx)
                    .unwrap()
            };
            let plus = sample(omega + h);
            let minus = sample(omega - h);
            for i in 0..6 {
                near(
                    result.displacement_derivative[i],
                    (plus.displacement[i] - minus.displacement[i]).scale(0.5 / h),
                    3e-6,
                );
                near(
                    result.velocity_derivative[i],
                    (plus.velocity[i] - minus.velocity[i]).scale(0.5 / h),
                    3e-6,
                );
            }
            real_near(
                result.average_viscous_power_derivative,
                (plus.average_viscous_power_w - minus.average_viscous_power_w) / (2.0 * h),
                3e-6,
            );
            let dp = result
                .port_load
                .as_ref()
                .unwrap()
                .average_absorbed_power_derivative;
            real_near(
                dp,
                (plus.port_load.unwrap().average_absorbed_power_w
                    - minus.port_load.unwrap().average_absorbed_power_w)
                    / (2.0 * h),
                3e-6,
            );
            assert!(result.power_balance_residual_derivative.abs() < 1e-18);
        }
    });
}

#[test]
fn original_drilling_has_exactly_zero_material_loss_and_storage_dispersion() {
    with_cx(|cx, _| {
        let source = source(proportional(0.02), CstDktMass::Consistent, cx);
        let mut q = vec![C64::ZERO; 6];
        q[5] = C64::new(2e-5, -3e-5);
        let kq = multiply(&source.model().k, &q);
        let mq = multiply(&source.model().m, &q);
        for omega in [0.5, 3.0, 10.0] {
            let f: Vec<_> = kq
                .iter()
                .zip(&mq)
                .map(|(&k, &m)| k - m.scale(omega * omega))
                .collect();
            let r = source
                .prepare_harmonic_response(omega, controls(), cx)
                .unwrap()
                .solve(&f, cx)
                .unwrap();
            for i in 0..6 {
                near(r.displacement[i], q[i], 2e-7);
            }
            assert!(r.average_viscous_power_w.abs() < 1e-30);
        }
    });
}

#[test]
fn soft_frequency_storage_survives_a_stiff_reference_without_matrix_subtraction() {
    with_cx(|cx, _| {
        let c = isotropic();
        let law = TensorProny::new_mandel(
            c,
            vec![PronyBranch::new(scaled(c, 1e18), 1.0).unwrap()],
            0.02,
        )
        .unwrap()
        .with_band((0.0, 1000.0))
        .unwrap();
        let source = TensorPronyShell::assemble(
            mesh(),
            vec![TensorPronyShellSection::new(law, 0.1, 4.0).unwrap()],
            &[0, 1],
            ShellSupport::Clamped,
            CstDktMass::Consistent,
            1000.0,
            budget(),
            cx,
        )
        .unwrap();
        let elastic = assemble_shell_with_mass(
            &mesh(),
            &PlateSection::isotropic(54.0 / 7.0, 2.0 / 7.0, 0.1, 4.0).unwrap(),
            &[0, 1],
            ShellSupport::Clamped,
            CstDktMass::Consistent,
        )
        .unwrap();
        let q = [
            C64::ZERO,
            C64::ZERO,
            C64::new(1e-5, 2e-6),
            C64::new(3e-6, -1e-6),
            C64::new(-2e-6, 3e-6),
            C64::ZERO,
        ];
        let omega = 1e-20;
        let s = C64::new(1.0, 0.01);
        let kq = multiply(&elastic.k, &q);
        let mq = multiply(&elastic.m, &q);
        let f: Vec<_> = kq
            .iter()
            .zip(mq)
            .map(|(&k, m)| s.conj() * k - m.scale(omega * omega))
            .collect();
        let result = source
            .prepare_harmonic_response(omega, controls(), cx)
            .unwrap()
            .solve(&f, cx)
            .unwrap();
        for i in 0..6 {
            near(result.displacement[i], q[i], 2e-7);
        }
    });
}

#[test]
fn vertex_strain_admission_catches_extrema_beyond_original_quadrature_points() {
    with_cx(|cx, _| {
        let source = source(proportional(0.02), CstDktMass::Consistent, cx);
        let (b, _) = dkt_curvature_quadrature(&[0.0, 0.4, 0.0], &[0.0, 0.0, 0.3], 0).unwrap();
        // Only last node w moves. The DKT curvature varies linearly, so vertices
        // can exceed every area midpoint even for a centered homogeneous layer.
        let mids = b.map(|p| [p[0][6], p[1][6], p[2][6]]);
        let norm = |k: [f64; 3]| {
            let [x, y, s] = k.map(|v| v * 0.05);
            (x * x + y * y + 0.5 * s * s + 0.16 * (x + y) * (x + y)).sqrt()
        };
        let qm = mids.map(norm).into_iter().fold(0.0, f64::max);
        let qv = (0..3)
            .map(|v| {
                norm(core::array::from_fn(|i| {
                    mids[v][i] + mids[(v + 2) % 3][i] - mids[(v + 1) % 3][i]
                }))
            })
            .fold(0.0, f64::max);
        assert!(qv > 1.1 * qm);
        let amplitude = 0.02 / (0.5 * (qm + qv));
        assert!(amplitude * qm < 0.02 && amplitude * qv > 0.02);
        let mut q = vec![C64::ZERO; 6];
        q[2] = C64::from_re(amplitude);
        let omega = 2.0;
        let kq = multiply(&source.model().k, &q);
        let mq = multiply(&source.model().m, &q);
        let f: Vec<_> = kq
            .iter()
            .zip(mq)
            .map(|(&k, m)| scalar(omega).0.conj() * k - m.scale(omega * omega))
            .collect();
        assert!(matches!(
            source
                .prepare_harmonic_response(omega, controls(), cx)
                .unwrap()
                .solve(&f, cx),
            Err(TensorShellHarmonicError::Material { .. })
        ));
    });
}

#[test]
fn original_source_bands_budgets_and_cancellation_refuse_without_partial_response() {
    with_cx(|cx, gate| {
        let source = source(proportional(0.02), CstDktMass::Consistent, cx);
        assert!(matches!(
            source.prepare_harmonic_response(101.0, controls(), cx),
            Err(TensorShellHarmonicError::Material { .. })
        ));
        assert!(
            source
                .prepare_harmonic_response(0.0, controls(), cx)
                .is_err()
        );
        assert!(
            source
                .prepare_harmonic_response(
                    1.0,
                    HarmonicResponseControls {
                        maximum_payload_bytes: 1,
                        ..controls()
                    },
                    cx
                )
                .is_err()
        );
        assert!(
            source
                .prepare_harmonic_response(1.0, controls(), cx)
                .unwrap()
                .with_frequency_derivative(1, cx)
                .is_err()
        );
        let too_large = vec![C64::from_re(1e6); 6];
        assert!(matches!(
            source
                .prepare_harmonic_response(1.0, controls(), cx)
                .unwrap()
                .solve(&too_large, cx),
            Err(TensorShellHarmonicError::Material { .. })
        ));
        gate.request();
        assert!(matches!(
            source.prepare_harmonic_response(1.0, controls(), cx),
            Err(TensorShellHarmonicError::Plate(PlateError::Cancelled))
        ));
    });
}
