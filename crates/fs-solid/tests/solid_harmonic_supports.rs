//! G1/G3 support-driven solid motion with native consistent inertia.

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_math::c64::C64;
use fs_modal::harmonic::{HarmonicResponseControls, HarmonicResponseError};
use fs_solid::linear3::{
    TetAssemblyBudget, TetElasticAssembly, TetElasticMaterial, TetLinearElasticProblem,
    TetMaterialField,
};
use fs_sparse::{Coo, Csr};

fn with_cx(cancelled: bool, run: impl FnOnce(&Cx<'_>)) {
    let gate = CancelGate::new_clock_free();
    if cancelled {
        gate.request();
    }
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        run(&Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 62,
                kernel_id: 64,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        ));
    });
}

fn controls() -> HarmonicResponseControls {
    HarmonicResponseControls {
        maximum_coordinates: 32,
        maximum_payload_bytes: 1 << 20,
        maximum_relative_residual: 1e-11,
    }
}

fn csr(values: &[f64], n: usize) -> Csr {
    let mut matrix = Coo::new(n, n);
    for row in 0..n {
        for col in 0..n {
            matrix.push(row, col, values[row * n + col]);
        }
    }
    matrix.assemble()
}

fn close(actual: C64, expected: C64) {
    assert!(
        (actual - expected).abs() <= 2e-11 * expected.abs().max(1e-5),
        "actual={actual:?}, expected={expected:?}"
    );
}

fn native_tet(cx: &Cx<'_>) -> TetElasticAssembly {
    let material = TetElasticMaterial::try_new(120.0, 120.0, 0.0, [9; 32]).unwrap();
    TetLinearElasticProblem {
        nodes_m: &[
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ],
        tetrahedra: &[[0, 1, 2, 3]],
        materials: TetMaterialField::Uniform(&material),
        fixed_dofs: &[],
        budget: TetAssemblyBudget::standard(),
    }
    .assemble(cx)
    .unwrap()
}

#[test]
fn original_global_support_labels_retain_cross_mass_damping_and_support_loads() {
    with_cx(false, |cx| {
        let source = TetElasticAssembly {
            mass: csr(&[2.0, 0.3, 0.3, 1.5], 2),
            stiffness: csr(&[7.0, -7.0, -7.0, 7.0], 2),
            free_dofs: vec![1, 7],
            total_mass_kg: 4.1,
            minimum_scaled_jacobian: 1.0,
            element_volumes_m3: vec![1.0],
        };
        let damping = csr(&[0.4, -0.2, -0.2, 0.6], 2);
        let omega = 1.7;
        let prepared = source
            .prepare_supported_harmonic_response(&[1], Some(&damping), omega, controls(), cx)
            .unwrap();
        assert!(std::ptr::eq(prepared.source(), &source));
        assert_eq!(prepared.represented_dofs(), &[1, 7]);
        assert_eq!(prepared.support_dofs(), &[1]);
        assert_eq!(prepared.free_dofs(), &[7]);
        let support = C64::new(0.03, -0.04);
        let force = [C64::new(0.1, -0.4), C64::new(0.5, 0.2)];
        let actual = prepared.solve(&force, &[support], cx).unwrap();
        let a00 = C64::new(7.0 - omega * omega * 2.0, -omega * 0.4);
        let a01 = C64::new(-7.0 - omega * omega * 0.3, omega * 0.2);
        let a11 = C64::new(7.0 - omega * omega * 1.5, -omega * 0.6);
        let free = (force[1] - a01 * support) / a11;
        close(actual.displacement[0], support);
        close(actual.displacement[1], free);
        let reaction = a00 * support + a01 * free - force[0];
        close(actual.support_reaction[0], reaction);
        let velocity = [
            C64::new(support.im, -support.re).scale(omega),
            C64::new(free.im, -free.re).scale(omega),
        ];
        let support_power = 0.5 * (reaction.conj() * velocity[0]).re;
        let heat = 0.5
            * (velocity[0].conj() * (velocity[0].scale(0.4) - velocity[1].scale(0.2))
                + velocity[1].conj() * (velocity[1].scale(0.6) - velocity[0].scale(0.2)))
            .re;
        assert!((actual.average_support_power_w - support_power).abs() < 1e-13);
        assert!((actual.average_viscous_power_w - heat).abs() < 1e-13);
        assert!(actual.power_balance_residual_w.abs() < 1e-12);
    });
}

#[test]
fn native_tetrahedron_base_excitation_keeps_consistent_inertia() {
    with_cx(false, |cx| {
        let source = native_tet(cx);
        let supports: Vec<_> = (0..9).collect();
        for omega in [0.0, 0.7, 1.6] {
            let prepared = source
                .prepare_supported_harmonic_response(&supports, None, omega, controls(), cx)
                .unwrap();
            assert_eq!(prepared.free_dofs(), &[9, 10, 11]);
            let base = C64::new(0.02, if omega == 0.0 { 0.0 } else { -0.01 });
            let mut prescribed = vec![C64::ZERO; 9];
            for node in 0..3 {
                prescribed[3 * node] = base;
            }
            let actual = prepared.solve(&[C64::ZERO; 12], &prescribed, cx).unwrap();
            // Independent P1 tetrahedron integrals at nu=0: K_xx=E*V/2=10,
            // M_33=rho*V/10=2, and each of three support cross masses is 1.
            let tip = base.scale((10.0 + 3.0 * omega * omega) / (10.0 - 2.0 * omega * omega));
            close(actual.displacement[9], tip);
            close(actual.displacement[10], C64::ZERO);
            close(actual.displacement[11], C64::ZERO);
            assert!(actual.relative_residual < 1e-11);
            assert!(actual.power_balance_residual_w.abs() < 1e-12);
        }
    });
}

#[test]
fn slanted_native_tet_admits_only_assembly_roundoff_without_changing_the_source() {
    with_cx(false, |cx| {
        let material = TetElasticMaterial::try_new(120.0, 120.0, 0.27, [9; 32]).unwrap();
        let source = TetLinearElasticProblem {
            nodes_m: &[
                [0.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [0.2, 1.1, 0.0],
                [0.1, 0.3, 0.9],
            ],
            tetrahedra: &[[0, 1, 2, 3]],
            materials: TetMaterialField::Uniform(&material),
            fixed_dofs: &[],
            budget: TetAssemblyBudget::standard(),
        }
        .assemble(cx)
        .unwrap();
        // A real non-axis-aligned element has last-bit skew because the native
        // operator contracts B^T D B independently in the two triangles.
        assert!((0..12).any(|row| {
            (0..row).any(|col| source.stiffness.get(row, col) != source.stiffness.get(col, row))
        }));
        let original_mass = source.mass.clone();
        let original_stiffness = source.stiffness.clone();
        let supports: Vec<_> = (0..9).collect();
        for omega in [0.0, 1.3] {
            let prepared = source
                .prepare_supported_harmonic_response(&supports, None, omega, controls(), cx)
                .unwrap();
            let base = C64::new(0.02, if omega == 0.0 { 0.0 } else { -0.01 });
            let mut prescribed = vec![C64::ZERO; 9];
            for node in 0..3 {
                prescribed[3 * node] = base;
            }
            let actual = prepared.solve(&[C64::ZERO; 12], &prescribed, cx).unwrap();
            if omega == 0.0 {
                // With the base translated statically, the whole unconstrained
                // tip translates rigidly, independently of element geometry.
                close(actual.displacement[9], base);
                close(actual.displacement[10], C64::ZERO);
                close(actual.displacement[11], C64::ZERO);
            }
            for row in 0..12 {
                let original_balance = (0..12).fold(C64::ZERO, |sum, col| {
                    sum + actual.displacement[col].scale(
                        original_stiffness.get(row, col)
                            - omega * omega * original_mass.get(row, col),
                    )
                });
                let reaction = if row < 9 {
                    actual.support_reaction[row]
                } else {
                    C64::ZERO
                };
                assert!((original_balance - reaction).abs() < 1e-12);
            }
            assert!(actual.relative_residual < 1e-11);
            assert!(actual.power_balance_residual_w.abs() < 1e-12);
        }
        assert_eq!(source.mass, original_mass);
        assert_eq!(source.stiffness, original_stiffness);

        // A constitutive/nonconservative skew is not assembly roundoff and
        // must not be silently replaced by a reciprocal elastic operator.
        let mut invalid = source;
        let mut stiffness = vec![0.0; 144];
        for row in 0..12 {
            for col in 0..12 {
                stiffness[12 * row + col] = invalid.stiffness.get(row, col);
            }
        }
        stiffness[1] += 1e-5;
        invalid.stiffness = csr(&stiffness, 12);
        assert!(matches!(
            invalid.prepare_supported_harmonic_response(&supports, None, 1.3, controls(), cx),
            Err(HarmonicResponseError::Invalid(
                "solid stiffness asymmetry exceeds assembly roundoff"
            ))
        ));
    });
}

#[test]
fn fully_prescribed_rigid_translation_recovers_total_inertial_reaction() {
    with_cx(false, |cx| {
        let source = native_tet(cx);
        let supports: Vec<_> = (0..12).collect();
        let omega = 1.3;
        let prepared = source
            .prepare_supported_harmonic_response(&supports, None, omega, controls(), cx)
            .unwrap();
        assert!(prepared.free_dofs().is_empty());
        let amplitude = C64::new(0.1, -0.07);
        let mut motion = vec![C64::ZERO; 12];
        for node in 0..4 {
            motion[3 * node + 1] = amplitude;
        }
        let response = prepared.solve(&[C64::ZERO; 12], &motion, cx).unwrap();
        let total = (0..4).fold(C64::ZERO, |sum, node| {
            sum + response.support_reaction[3 * node + 1]
        });
        close(
            total,
            amplitude.scale(-omega * omega * source.total_mass_kg),
        );
        assert!(response.average_support_power_w.abs() < 1e-12);
    });
}

#[test]
fn absent_support_budget_and_cancel_refuse_before_publication() {
    with_cx(false, |cx| {
        let source = native_tet(cx);
        assert!(matches!(
            source.prepare_supported_harmonic_response(&[12], None, 1.0, controls(), cx),
            Err(HarmonicResponseError::Invalid(_))
        ));
        assert!(matches!(
            source.prepare_supported_harmonic_response(
                &[0],
                None,
                1.0,
                HarmonicResponseControls {
                    maximum_payload_bytes: 1,
                    ..controls()
                },
                cx
            ),
            Err(HarmonicResponseError::Budget { .. })
        ));
    });
    with_cx(true, |cx| {
        let source = TetElasticAssembly {
            mass: csr(&[1.0], 1),
            stiffness: csr(&[1.0], 1),
            free_dofs: vec![4],
            total_mass_kg: 1.0,
            minimum_scaled_jacobian: 1.0,
            element_volumes_m3: vec![1.0],
        };
        assert!(matches!(
            source.prepare_supported_harmonic_response(&[4], None, 1.0, controls(), cx),
            Err(HarmonicResponseError::Cancelled)
        ));
    });
}
