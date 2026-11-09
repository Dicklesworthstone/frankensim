//! Complete moving-support inertia, layer strain and physical boundary work.
use super::*;
use fs_fft::harmonic_window::HarmonicWindow;
use fs_modal::harmonic::HarmonicPortControls;
use fs_plate::shell::viscosity::harmonic::{
    PreparedTensorPronyShellOperator, TensorShellMaterialBudget,
};

fn ports() -> HarmonicPortControls {
    HarmonicPortControls {
        maximum_ports: 2,
        maximum_projection_terms: 100_000,
        maximum_payload_bytes: 1 << 20,
    }
}
fn impedance() -> [C64; 4] {
    [
        C64::new(0.8, 0.3),
        C64::new(0.1, -0.2),
        C64::new(-0.1, -0.2),
        C64::new(1.2, -0.4),
    ]
}
fn force_map(n: usize) -> Vec<f64> {
    (0..n)
        .flat_map(|i| [0.3 + 0.07 * (i % 4) as f64, -0.2 + 0.11 * (i % 5) as f64])
        .collect()
}
fn inner(a: C64, b: C64) -> f64 {
    a.re * b.re + a.im * b.im
}
fn scalar(a: f64, b: f64, scale: f64) {
    assert!(
        (a - b).abs() <= 3e-10 * scale.abs().max(1e-30),
        "{a:e} != {b:e}"
    );
}
fn action(
    operator: &PreparedTensorPronyShellOperator<'_>,
    q: &[C64],
    b: &[f64],
    z: &[C64; 4],
) -> Vec<C64> {
    let w = operator.omega_rad_s();
    (0..q.len())
        .map(|i| {
            (0..q.len()).fold(C64::ZERO, |sum, j| {
                let mut a = C64::new(
                    operator.stiffness().get(i, j) - w * w * operator.mass().get(i, j),
                    -w * operator.damping().get(i, j),
                );
                for p in 0..2 {
                    for r in 0..2 {
                        let load = z[2 * p + r].scale(b[2 * i + p] * b[2 * j + r]);
                        a = a + C64::new(load.im, -load.re).scale(w);
                    }
                }
                sum + a * q[j]
            })
        })
        .collect()
}
fn full_source(cx: &Cx<'_>) -> TensorPronyShell {
    TensorPronyShell::assemble(
        ShellMesh::new(
            vec![[0.0, 0.0, 0.0], [0.4, 0.0, 0.0], [0.0, 0.3, 0.0]],
            vec![[0, 1, 2]],
        )
        .unwrap(),
        vec![TensorPronyShellSection::from_laminate(laminate()).unwrap()],
        &[],
        ShellSupport::Free,
        CstDktMass::Consistent,
        0.0,
        budget(),
        cx,
    )
    .unwrap()
}
fn observed_heat(shell: &TensorPronyShell, q: &[C64], w: f64, cx: &Cx<'_>) -> f64 {
    let observer = shell
        .prepare_material_spectrum_observer(
            &[w],
            TensorShellMaterialBudget {
                maximum_points: 12,
                maximum_payload_bytes: 1 << 22,
                maximum_work: 1_000_000_000,
            },
            cx,
        )
        .unwrap();
    let period = core::f64::consts::TAU / w;
    observer
        .observe_window(&[q], HarmonicWindow::new(0.137, period).unwrap(), cx)
        .unwrap()
        .material_heat_j
        / period
}

#[test]
fn moving_translations_and_rotations_keep_full_loaded_reactions_and_original_layer_heat() {
    with_cx(|cx| {
        let shell = full_source(cx);
        let support: Vec<_> = (0..12).collect();
        let b = force_map(18);
        let z = impedance();
        let q: Vec<_> = (0..18)
            .map(|i| C64::new(1e-5 * (i as f64 - 6.0), 0.7e-5 * (3.0 - (i % 7) as f64)))
            .collect();
        for w in [0.4, 2.4, 9.0] {
            let original = shell.prepare_harmonic_operator(w, controls(), cx).unwrap();
            // Signed offset-layer translation/rotation inertia is physically present.
            assert!(original.mass().get(1, 3).abs() > 0.0);
            let aq = action(&original, &q, &b, &z);
            let mut force = aq.clone();
            for (i, f) in force[..12].iter_mut().enumerate() {
                *f = C64::new(2e-6 * (i as f64 - 4.0), -3e-6);
            }
            let prepared = shell
                .prepare_supported_harmonic_response_with_port_impedance(
                    &support,
                    w,
                    &b,
                    2,
                    &z,
                    controls(),
                    ports(),
                    cx,
                )
                .unwrap();
            assert_eq!(prepared.support_dofs(), &support);
            assert_eq!(prepared.free_dofs(), &[12, 13, 14, 15, 16, 17]);
            assert_eq!(prepared.represented_dofs(), &(0..18).collect::<Vec<_>>());
            assert!(core::ptr::eq(prepared.source(), &shell));
            let result = prepared.solve(&force, &q[..12], cx).unwrap();
            for i in 0..18 {
                near(result.displacement[i], q[i], 2e-9);
            }
            for i in 0..12 {
                near(result.support_reaction[i], aq[i] - force[i], 2e-9);
            }
            let velocity: Vec<_> = q.iter().map(|v| C64::new(v.im, -v.re).scale(w)).collect();
            let flow: [C64; 2] = core::array::from_fn(|p| {
                velocity
                    .iter()
                    .enumerate()
                    .fold(C64::ZERO, |v, (i, u)| v + u.scale(b[2 * i + p]))
            });
            let effort: [C64; 2] =
                core::array::from_fn(|p| z[2 * p] * flow[0] + z[2 * p + 1] * flow[1]);
            let boundary = 0.5 * (inner(flow[0], effort[0]) + inner(flow[1], effort[1]));
            let intrinsic = 0.5
                * (0..18)
                    .map(|i| {
                        (0..18)
                            .map(|j| original.damping().get(i, j) * inner(velocity[i], velocity[j]))
                            .sum::<f64>()
                    })
                    .sum::<f64>();
            let load = result.port_load.as_ref().unwrap();
            for p in 0..2 {
                near(load.flow[p], flow[p], 2e-9);
                near(load.effort[p], effort[p], 2e-9);
            }
            scalar(load.average_absorbed_power_w, boundary, boundary);
            scalar(result.average_viscous_power_w, intrinsic, intrinsic);
            scalar(
                observed_heat(&shell, &result.displacement, w, cx),
                intrinsic,
                intrinsic,
            );
            let support_work = 0.5
                * (0..12)
                    .map(|i| inner(aq[i] - force[i], velocity[i]))
                    .sum::<f64>();
            scalar(
                result.average_support_power_w,
                support_work,
                support_work.abs() + boundary + intrinsic,
            );
            scalar(result.power_balance_residual_w, 0.0, intrinsic + boundary);
            let repeated = prepared.solve(&force, &q[..12], cx).unwrap();
            assert_eq!(result.displacement, repeated.displacement);
            assert_eq!(result.support_reaction, repeated.support_reaction);
        }
    });
}

#[test]
fn original_labels_and_empty_all_prescribed_limits_preserve_the_same_shell() {
    with_cx(|cx| {
        let shell = source(laminate(), CstDktMass::Consistent, cx).unwrap();
        let b = force_map(6);
        let z = impedance();
        let w = 2.4;
        let force = [C64::new(2e-6, -3e-6); 6];
        let ordinary = shell
            .prepare_harmonic_response_with_port_impedance(w, &b, 2, &z, controls(), ports(), cx)
            .unwrap()
            .solve(&force, cx)
            .unwrap();
        let empty = shell
            .prepare_supported_harmonic_response_with_port_impedance(
                &[],
                w,
                &b,
                2,
                &z,
                controls(),
                ports(),
                cx,
            )
            .unwrap();
        assert_eq!(empty.represented_dofs(), &[12, 13, 14, 15, 16, 17]);
        let result = empty.solve(&force, &[], cx).unwrap();
        for i in 0..6 {
            near(result.displacement[i], ordinary.displacement[i], 3e-10);
        }
        scalar(
            result.average_viscous_power_w,
            ordinary.average_viscous_power_w,
            ordinary.average_viscous_power_w,
        );
        let labels = [12, 13, 14, 15, 16, 17];
        let q = [C64::new(0.0001, -0.00003); 6];
        let all = shell
            .prepare_supported_harmonic_response_with_port_impedance(
                &labels,
                w,
                &b,
                2,
                &z,
                controls(),
                ports(),
                cx,
            )
            .unwrap();
        assert!(all.free_dofs().is_empty());
        let result = all.solve(&force, &q, cx).unwrap();
        assert_eq!(result.displacement, q);
        assert_eq!(result.relative_residual, 0.0);
        let original = shell.prepare_harmonic_operator(w, controls(), cx).unwrap();
        let aq = action(&original, &q, &b, &z);
        for i in 0..6 {
            near(result.support_reaction[i], aq[i] - force[i], 3e-10);
        }
        let bare = shell
            .prepare_supported_harmonic_response(&labels, w, controls(), cx)
            .unwrap()
            .solve(&force, &q, cx)
            .unwrap();
        assert!(bare.port_load.is_none());
        scalar(
            bare.average_viscous_power_w,
            result.average_viscous_power_w,
            result.average_viscous_power_w,
        );
    });
}

#[test]
fn source_support_material_and_load_refusals_leave_cached_response_reusable() {
    with_cx(|cx| {
        let shell = source(laminate(), CstDktMass::Consistent, cx).unwrap();
        for labels in [&[0][..], &[18][..], &[13, 12][..], &[12, 12][..]] {
            assert!(
                shell
                    .prepare_supported_harmonic_response(labels, 2.4, controls(), cx)
                    .is_err()
            );
        }
        for w in [0.0, 21.0, f64::NAN] {
            assert!(
                shell
                    .prepare_supported_harmonic_response(&[12], w, controls(), cx)
                    .is_err()
            );
        }
        assert!(
            shell
                .prepare_supported_harmonic_response(
                    &[12],
                    2.4,
                    HarmonicResponseControls {
                        maximum_payload_bytes: 0,
                        ..controls()
                    },
                    cx
                )
                .is_err()
        );
        let b = force_map(6);
        let z = impedance();
        for (map, load, limits) in [
            (&b[..11], z, ports()),
            (
                &b[..],
                [C64::from_re(-1.0), C64::ZERO, C64::ZERO, C64::ONE],
                ports(),
            ),
            (
                &b[..],
                z,
                HarmonicPortControls {
                    maximum_projection_terms: 0,
                    ..ports()
                },
            ),
        ] {
            assert!(
                shell
                    .prepare_supported_harmonic_response_with_port_impedance(
                        &[12],
                        2.4,
                        map,
                        2,
                        &load,
                        controls(),
                        limits,
                        cx,
                    )
                    .is_err()
            );
        }
        let prepared = shell
            .prepare_supported_harmonic_response(&[12, 13, 14, 15, 16, 17], 2.4, controls(), cx)
            .unwrap();
        let force = [C64::ZERO; 6];
        let q = [C64::new(0.0001, 0.0002); 6];
        let expected = prepared.solve(&force, &q, cx).unwrap();
        assert!(prepared.solve(&force, &[C64::from_re(0.5); 6], cx).is_err());
        assert!(prepared.solve(&force, &q[..5], cx).is_err());
        assert_eq!(
            expected.displacement,
            prepared.solve(&force, &q, cx).unwrap().displacement
        );
    });
}

#[test]
fn support_solve_and_preparation_obey_the_original_cancellation_gate() {
    let gate = CancelGate::new_clock_free();
    fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default()).scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            fs_exec::StreamKey {
                seed: 59,
                kernel_id: 24,
                tile: 0,
                iteration: 0,
            },
            fs_exec::Budget::INFINITE,
            fs_exec::ExecMode::Deterministic,
        );
        let shell = source(laminate(), CstDktMass::Consistent, &cx).unwrap();
        let prepared = shell
            .prepare_supported_harmonic_response(&[12], 2.4, controls(), &cx)
            .unwrap();
        gate.request();
        assert!(prepared.solve(&[C64::ZERO; 6], &[C64::ZERO], &cx).is_err());
        assert!(
            shell
                .prepare_supported_harmonic_response(&[12], 2.4, controls(), &cx)
                .is_err()
        );
    });
}
