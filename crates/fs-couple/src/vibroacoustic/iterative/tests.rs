//! G1/G3: original structural reaction, outgoing acoustic fields and port work.
use super::*;
use fs_alloc::{ArenaConfig, ArenaPool};
use fs_bem::helmholtz::{PanelRule, Preconditioner, solve_radiation_with_rule};
use fs_exec::{Budget, ExecMode, StreamKey};
use fs_fmm::DirectionalOptions;
use fs_modal::harmonic::HarmonicPortControls;

fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new();
    let pool = ArenaPool::new(ArenaConfig::default());
    pool.scope(|arena| {
        f(&Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 17,
                kernel_id: 71,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        ))
    })
}
fn controls() -> StructuralRadiationControls {
    StructuralRadiationControls {
        mechanical: HarmonicResponseControls {
            maximum_coordinates: 64,
            maximum_payload_bytes: 8_000_000,
            maximum_relative_residual: 1e-10,
        },
        maximum_trace_entries: 1_000_000,
        maximum_mobility_entries: 4096,
        maximum_payload_bytes: 8_000_000,
        maximum_relative_velocity_defect: 1e-8,
        maximum_relative_power_defect: 1e-8,
    }
}
fn options() -> IterativeOptions {
    IterativeOptions {
        tol: 2e-11,
        restart: 40,
        max_iters: 240,
        fmm_order: 6,
        leaf_capacity: 256,
        panel_rule: PanelRule::TriangleNear {
            max_pairs: 2_000_000,
        },
        precond: Preconditioner::NearFieldBlockJacobi {
            max_block_panels: 20,
        },
        ..IterativeOptions::default()
    }
}
fn relative(a: &[C64], b: &[C64]) -> f64 {
    a.iter()
        .zip(b)
        .fold(0.0_f64, |s, (a, b)| s.hypot((*a - *b).abs()))
        / b.iter().fold(0.0_f64, |s, b| s.hypot(b.abs())).max(1e-300)
}
fn matrix_pencil<'a>(
    mass: &'a [f64],
    damping: &'a [f64],
    stiffness: &'a [f64],
    n: usize,
) -> StructuralRadiationPencil<'a> {
    StructuralRadiationPencil {
        mass,
        damping,
        stiffness,
        coordinates: n,
    }
}

#[test]
fn pulsating_structural_sphere_converges_to_independent_loaded_impedance() {
    let medium = Medium {
        density: 1.2,
        sound_speed: 3.0,
    };
    let (r, omega, m, damping, k) = (0.1, 24.0, 0.025, 0.2, 10.0);
    let ka = omega * r / medium.sound_speed;
    let exact_z = C64::new(ka * ka, -ka).scale(
        4.0 * core::f64::consts::PI * r * r * medium.density * medium.sound_speed / (1.0 + ka * ka),
    );
    let drive = C64::new(1.0, 0.3);
    let expected_velocity = drive / (C64::new(damping, k / omega - omega * m) + exact_z);
    let expected_power = 0.5 * exact_z.re * expected_velocity.norm_sq();
    with_cx(|cx| {
        let gate = CancelGate::new_clock_free();
        let mut previous = f64::INFINITY;
        for level in 1..=2 {
            let surface = SpherePanels::icosphere(r, level).unwrap();
            let traces = vec![vec![1.0; surface.areas().len()]];
            let source = StructuralRadiation::new(
                &surface,
                &traces,
                matrix_pencil(&[m], &[damping], &[k], 1),
                omega,
                medium,
                controls(),
                cx,
            )
            .unwrap();
            let prepared = source.prepare(&options(), &gate, cx).unwrap();
            let actual = prepared.solve(&[drive], &gate, cx).unwrap();
            let error =
                (actual.motion.velocity[0] - expected_velocity).abs() / expected_velocity.abs();
            let power_error =
                (actual.signed_radiation_power_w - expected_power).abs() / expected_power;
            eprintln!(
                "structural sphere level={level} panels={} velocity={error:e} power={power_error:e}",
                surface.areas().len()
            );
            assert!(error < previous);
            previous = error;
            assert!(actual.relative_interface_work_defect < 1e-10);
            assert!(actual.relative_total_work_defect < 1e-10);
            assert!(actual.structural_dissipation_w > 0.0 && actual.signed_radiation_power_w > 0.0);
            if level == 2 {
                assert!(error < 0.04 && power_error < 0.06);
            }
        }
    });
}

fn paired_spheres() -> SpherePanels {
    let sphere = SpherePanels::icosphere(0.01, 0).unwrap();
    let mut triangles = Vec::new();
    for z in [-1.0, 1.0] {
        triangles.extend(
            sphere
                .triangles()
                .unwrap()
                .iter()
                .map(|t| t.map(|p| [p[0], p[1], p[2] + z])),
        );
    }
    SpherePanels::from_triangles(triangles).unwrap()
}

#[test]
fn directional_full_structure_matches_original_loaded_dense_owner_with_gyroscopic_coupling() {
    let surface = paired_spheres();
    let medium = Medium {
        density: 1.2,
        sound_speed: 343.0,
    };
    let omega = 34.0 * medium.sound_speed;
    let mut traces = vec![vec![0.0; 40]; 2];
    for i in 0..40 {
        traces[usize::from(i >= 20)][i] = 1.0;
    }
    let m = [2e-6, 0.2e-6, 0.2e-6, 1.5e-6];
    let k = [400.0, -30.0, -30.0, 350.0];
    let c = [0.015, 0.004, -0.002, 0.01];
    let f = [C64::new(1.0, 0.3), C64::new(0.2, -0.1)];
    let mut settings = options();
    settings.panel_rule = PanelRule::Triangle {
        max_pairs: 2_000_000,
    };
    settings.directional = Some(DirectionalOptions {
        order: 8,
        leaf_capacity: 64,
        ..DirectionalOptions::for_tolerance(1e-2).unwrap()
    });
    with_cx(|cx| {
        let gate = CancelGate::new_clock_free();
        let source = StructuralRadiation::new(
            &surface,
            &traces,
            matrix_pencil(&m, &c, &k, 2),
            omega,
            medium,
            controls(),
            cx,
        )
        .unwrap();
        let prepared = source.prepare(&settings, &gate, cx).unwrap();
        let actual = prepared.solve(&f, &gate, cx).unwrap();
        assert!(actual.report.converged_euclidean());
        assert_eq!(actual.receipt.fmm_order, 8);
        assert!(
            source.prepare(&options(), &gate, cx).is_err(),
            "fixture must require high-kD route"
        );
        let mut impedance = vec![C64::ZERO; 4];
        let mut columns = Vec::new();
        for j in 0..2 {
            let velocity = traces[j]
                .iter()
                .map(|&v| C64::from_re(v))
                .collect::<Vec<_>>();
            let field = solve_radiation_with_rule(
                &surface,
                34.0,
                medium,
                &velocity,
                Formulation::BurtonMiller,
                settings.panel_rule,
            )
            .unwrap();
            for i in 0..2 {
                let mut z = ComplexSum::default();
                for ((&p, &shape), &area) in
                    field.pressure.iter().zip(&traces[i]).zip(surface.areas())
                {
                    z.scaled_product(p, shape, area);
                }
                impedance[i * 2 + j] = z.finish().unwrap();
            }
            columns.push(field);
        }
        let reference = PreparedHarmonicResponse::try_with_port_impedance(
            &m,
            &c,
            &k,
            2,
            &[1.0, 0.0, 0.0, 1.0],
            2,
            &impedance,
            omega,
            controls().mechanical,
            HarmonicPortControls {
                maximum_ports: 4,
                maximum_projection_terms: 1_000_000,
                maximum_payload_bytes: 8_000_000,
            },
            cx,
        )
        .unwrap()
        .solve(&f, cx)
        .unwrap();
        eprintln!(
            "directional original-coordinate velocity error={:e}",
            relative(&actual.motion.velocity, &reference.velocity)
        );
        assert!(relative(&actual.motion.velocity, &reference.velocity) < 2e-7);
        let pressure = (0..40)
            .map(|i| {
                columns[0].pressure[i] * reference.velocity[0]
                    + columns[1].pressure[i] * reference.velocity[1]
            })
            .collect::<Vec<_>>();
        eprintln!(
            "directional original-coordinate pressure error={:e} total-work={:e}",
            relative(&actual.field.pressure, &pressure),
            actual.relative_total_work_defect
        );
        assert!(relative(&actual.field.pressure, &pressure) < 2e-7);
        let load = reference.port_load.as_ref().unwrap();
        assert!(
            (load.average_absorbed_power_w - actual.signed_radiation_power_w).abs()
                < 2e-7 * load.average_absorbed_power_w.abs()
        );
        assert!(
            (actual.external_input_power_w - reference.average_input_power_w).abs()
                < 2e-7 * reference.average_input_power_w.abs()
        );
        assert!(actual.relative_total_work_defect < 1e-9);
    });
}

#[test]
fn original_coordinate_congruence_and_dark_coordinate_keep_fields_and_work() {
    let surface = SpherePanels::icosphere(0.1, 1).unwrap();
    let medium = Medium {
        density: 1.2,
        sound_speed: 3.0,
    };
    let (m, c, k) = (
        [0.02, 0.004, 0.004, 0.03],
        [0.3, 0.07, -0.03, 0.2],
        [12.0, -2.0, -2.0, 15.0],
    );
    let f = [C64::new(1.0, 0.4), C64::new(-0.3, 0.1)];
    let traces = vec![vec![1.0; 80], vec![0.0; 80]];
    let scale = [1000.0, 0.02];
    let congruent = |matrix: &[f64; 4]| {
        core::array::from_fn::<_, 4, _>(|i| matrix[i] * scale[i / 2] * scale[i % 2])
    };
    let (m2, c2, k2) = (congruent(&m), congruent(&c), congruent(&k));
    let traces2 = (0..2)
        .map(|r| traces[r].iter().map(|v| v * scale[r]).collect::<Vec<_>>())
        .collect::<Vec<_>>();
    with_cx(|cx| {
        let gate = CancelGate::new_clock_free();
        let a = StructuralRadiation::new(
            &surface,
            &traces,
            matrix_pencil(&m, &c, &k, 2),
            20.0,
            medium,
            controls(),
            cx,
        )
        .unwrap();
        let b = StructuralRadiation::new(
            &surface,
            &traces2,
            matrix_pencil(&m2, &c2, &k2, 2),
            20.0,
            medium,
            controls(),
            cx,
        )
        .unwrap();
        let x = a
            .prepare(&options(), &gate, cx)
            .unwrap()
            .solve(&f, &gate, cx)
            .unwrap();
        let y = b
            .prepare(&options(), &gate, cx)
            .unwrap()
            .solve(&[f[0].scale(scale[0]), f[1].scale(scale[1])], &gate, cx)
            .unwrap();
        assert!(
            x.motion.velocity[1].abs() > 0.01,
            "dark coordinates still move and couple structurally"
        );
        let restored = y
            .motion
            .velocity
            .iter()
            .zip(scale)
            .map(|(&v, s)| v.scale(s))
            .collect::<Vec<_>>();
        assert!(relative(&restored, &x.motion.velocity) < 2e-9);
        assert!(relative(&x.field.pressure, &y.field.pressure) < 2e-9);
        for (a, b) in [
            (x.external_input_power_w, y.external_input_power_w),
            (x.structural_dissipation_w, y.structural_dissipation_w),
            (x.signed_radiation_power_w, y.signed_radiation_power_w),
        ] {
            assert!((a - b).abs() < 2e-9 * a.abs().max(1e-20));
        }
        // Removing the full off-diagonal structural interaction would change
        // the observed bright coordinate even though the second trace is zero.
        let bare = a.mechanical.solve(&f, cx).unwrap();
        assert!(
            (bare.velocity[0] - x.motion.velocity[0]).abs() > 0.01 * x.motion.velocity[0].abs()
        );
    });
}

#[test]
fn nonlocal_surface_work_is_the_original_full_structural_dissipation() {
    let surface = SpherePanels::icosphere(0.1, 1).unwrap();
    let traces = vec![
        vec![1.0; 80],
        surface.normals().iter().map(|n| n[2]).collect(),
    ];
    let (m, c, k) = (
        [0.02, 0.004, 0.004, 0.03],
        [0.3, 0.07, -0.03, 0.2],
        [12.0, -2.0, -2.0, 15.0],
    );
    let medium = Medium {
        density: 1.2,
        sound_speed: 3.0,
    };
    with_cx(|cx| {
        let source = StructuralRadiation::new(
            &surface,
            &traces,
            matrix_pencil(&m, &c, &k, 2),
            20.0,
            medium,
            controls(),
            cx,
        )
        .unwrap();
        let pressure = (0..80)
            .map(|i| C64::new(0.7 + 0.13 * i as f64, (i as f64 * 0.4).sin()))
            .collect::<Vec<_>>();
        let force = source.traction(&pressure).unwrap();
        let response = source.mechanical.solve(&force, cx).unwrap();
        let mut normalized = vec![C64::ZERO; 80];
        source.apply(&pressure, &mut normalized).unwrap();
        let actual_velocity = source.trace(&response.velocity).unwrap();
        for i in 0..80 {
            assert!(
                (normalized[i] - actual_velocity[i].scale(medium.density * medium.sound_speed))
                    .abs()
                    < 1e-11
            );
        }
        let power = pressure
            .iter()
            .zip(&normalized)
            .zip(surface.areas())
            .map(|((&p, &v), &a)| {
                0.5 * (p * v.conj()).re * a / (medium.density * medium.sound_speed)
            })
            .sum::<f64>();
        assert!(power > 0.0);
        assert!((power - response.average_viscous_power_w).abs() < 2e-12 * power);
        let asymmetric = (source.mobility()[1] - source.mobility()[2]).abs();
        assert!(asymmetric > 0.01, "gyro coupling must not be symmetrized");
    });
}

#[test]
fn source_scope_singular_bare_resonance_budgets_and_retry_are_explicit() {
    let surface = SpherePanels::icosphere(0.1, 1).unwrap();
    let traces = vec![vec![1.0; 80]];
    let medium = Medium {
        density: 1.2,
        sound_speed: 3.0,
    };
    with_cx(|cx| {
        assert!(matches!(
            StructuralRadiation::new(
                &surface,
                &traces,
                matrix_pencil(&[1.0], &[0.0], &[400.0], 1),
                20.0,
                medium,
                controls(),
                cx
            ),
            Err(StructuralRadiationError::Mechanical(
                HarmonicResponseError::Singular { .. }
            ))
        ));
        assert!(
            StructuralRadiation::new(
                &surface,
                &traces,
                matrix_pencil(&[0.02], &[-0.1], &[12.0], 1),
                20.0,
                medium,
                controls(),
                cx
            )
            .is_err()
        );
        for control in [
            StructuralRadiationControls {
                maximum_trace_entries: 19,
                ..controls()
            },
            StructuralRadiationControls {
                maximum_mobility_entries: 0,
                ..controls()
            },
            StructuralRadiationControls {
                maximum_payload_bytes: 1,
                ..controls()
            },
        ] {
            assert!(
                StructuralRadiation::new(
                    &surface,
                    &traces,
                    matrix_pencil(&[0.02], &[0.1], &[12.0], 1),
                    20.0,
                    medium,
                    control,
                    cx
                )
                .is_err()
            );
        }
        let source = StructuralRadiation::new(
            &surface,
            &traces,
            matrix_pencil(&[0.02], &[0.1], &[12.0], 1),
            20.0,
            medium,
            controls(),
            cx,
        )
        .unwrap();
        let gate = CancelGate::new_clock_free();
        let prepared = source.prepare(&options(), &gate, cx).unwrap();
        let force = [C64::new(1.0, -0.2)];
        let first = prepared.solve(&force, &gate, cx).unwrap();
        assert!(
            prepared
                .solve(&[C64::new(f64::NAN, 0.0)], &gate, cx)
                .is_err()
        );
        let cancelled = CancelGate::new_clock_free();
        cancelled.request();
        assert!(matches!(
            prepared.solve(&force, &cancelled, cx),
            Err(StructuralRadiationError::Cancelled)
        ));
        let retry = prepared.solve(&force, &gate, cx).unwrap();
        assert_eq!(first.motion.velocity, retry.motion.velocity);
        assert_eq!(first.field.pressure, retry.field.pressure);
        let zero = prepared.solve(&[C64::ZERO], &gate, cx).unwrap();
        assert_eq!(zero.field.pressure, vec![C64::ZERO; 80]);
        assert_eq!(zero.external_input_power_w, 0.0);
    });
}
