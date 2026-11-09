//! G1/G3: actual radiation resolves a bare resonance in the full joint equations.
use super::*;
use fs_alloc::{ArenaConfig, ArenaPool};
use fs_bem::helmholtz::{PanelRule, Preconditioner, solve_radiation_with_rule};
use fs_exec::{Budget, ExecMode, StreamKey};
use fs_fmm::DirectionalOptions;
use fs_modal::harmonic::{HarmonicPortControls, PreparedHarmonicResponse};

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
fn controls() -> JointStructuralRadiationControls {
    JointStructuralRadiationControls {
        mechanical: HarmonicResponseControls {
            maximum_coordinates: 64,
            maximum_payload_bytes: 8_000_000,
            maximum_relative_residual: 1e-9,
        },
        maximum_trace_entries: 1_000_000,
        maximum_payload_bytes: 16_000_000,
        maximum_relative_power_defect: 1e-8,
    }
}
fn options() -> IterativeOptions {
    IterativeOptions {
        tol: 2e-10,
        restart: 60,
        max_iters: 300,
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
fn pencil<'a>(m: &'a [f64], c: &'a [f64], k: &'a [f64], n: usize) -> StructuralRadiationPencil<'a> {
    StructuralRadiationPencil {
        coordinates: n,
        mass: m,
        damping: c,
        stiffness: k,
    }
}
fn relative(a: &[C64], b: &[C64]) -> f64 {
    let difference = a
        .iter()
        .zip(b)
        .flat_map(|(a, b)| [a.re - b.re, a.im - b.im])
        .collect::<Vec<_>>();
    let reference = b.iter().flat_map(|b| [b.re, b.im]).collect::<Vec<_>>();
    fs_solver::norm2(&difference) / fs_solver::norm2(&reference).max(1e-300)
}

#[test]
fn actual_radiation_resolves_exact_bare_resonance_against_the_independent_sphere() {
    let medium = Medium {
        density: 1.2,
        sound_speed: 1.0,
    };
    let (radius, omega, mass, stiffness) = (0.1, 8.0, 0.25, 16.0);
    let ka = omega * radius / medium.sound_speed;
    let z = C64::new(ka * ka, -ka).scale(
        4.0 * core::f64::consts::PI * radius * radius * medium.density * medium.sound_speed
            / (1.0 + ka * ka),
    );
    let force = C64::new(1.0, 0.3);
    let exact = force / z;
    with_cx(|cx| {
        assert!(matches!(
            PreparedHarmonicResponse::try_new(
                &[mass],
                &[0.0],
                &[stiffness],
                1,
                omega,
                controls().mechanical,
                cx
            ),
            Err(HarmonicResponseError::Singular { .. })
        ));
        let gate = CancelGate::new_clock_free();
        let mut previous = f64::INFINITY;
        for level in 1..=2 {
            let surface = SpherePanels::icosphere(radius, level).unwrap();
            let traces = vec![vec![1.0; surface.areas().len()]];
            let owner = PreparedJointStructuralRadiation::new(
                &surface,
                &traces,
                pencil(&[mass], &[0.0], &[stiffness], 1),
                omega,
                medium,
                &options(),
                controls(),
                &gate,
                cx,
            )
            .unwrap();
            let solution = owner.solve(&[force], &gate, cx).unwrap();
            let error = (solution.velocity[0] - exact).abs() / exact.abs();
            eprintln!(
                "joint bare sphere level={level} velocity error={error:e} equation={:e} power={:e}",
                solution.relative_mechanical_equation_defect, solution.relative_total_work_defect
            );
            assert!(error < previous);
            previous = error;
            assert_eq!(solution.structural_dissipation_w, 0.0);
            assert!(solution.signed_radiation_power_w > 0.0);
            assert!(solution.relative_total_work_defect < 1e-9);
            assert!(relative(&solution.opposing_pressure_force, &[force]) < 1e-9);
            if level == 2 {
                assert!(error < 0.04);
            }
        }
    });
}

fn paired_spheres() -> SpherePanels {
    let one = SpherePanels::icosphere(0.01, 0).unwrap();
    let mut triangles = Vec::new();
    for z in [-1.0, 1.0] {
        triangles.extend(
            one.triangles()
                .unwrap()
                .iter()
                .map(|t| t.map(|p| [p[0], p[1], p[2] + z])),
        );
    }
    SpherePanels::from_triangles(triangles).unwrap()
}
#[test]
fn directional_joint_fields_match_dense_original_loaded_gyro_equations() {
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
    let force = [C64::new(1.0, 0.3), C64::new(0.2, -0.1)];
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
        let owner = PreparedJointStructuralRadiation::new(
            &surface,
            &traces,
            pencil(&m, &c, &k, 2),
            omega,
            medium,
            &settings,
            controls(),
            &gate,
            cx,
        )
        .unwrap();
        let solution = owner.solve(&force, &gate, cx).unwrap();
        let mut impedance = vec![C64::ZERO; 4];
        let mut columns = Vec::new();
        for j in 0..2 {
            let v = traces[j]
                .iter()
                .map(|&t| C64::from_re(t))
                .collect::<Vec<_>>();
            let field = solve_radiation_with_rule(
                &surface,
                34.0,
                medium,
                &v,
                Formulation::BurtonMiller,
                settings.panel_rule,
            )
            .unwrap();
            for i in 0..2 {
                let mut sum = ComplexSum::default();
                for ((&p, &t), &a) in field.pressure.iter().zip(&traces[i]).zip(surface.areas()) {
                    sum.scaled_product(p, t, a);
                }
                impedance[i * 2 + j] = sum.finish().unwrap();
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
        .solve(&force, cx)
        .unwrap();
        let pressure = (0..40)
            .map(|i| {
                columns[0].pressure[i] * reference.velocity[0]
                    + columns[1].pressure[i] * reference.velocity[1]
            })
            .collect::<Vec<_>>();
        let ev = relative(&solution.velocity, &reference.velocity);
        let ep = relative(&solution.field.pressure, &pressure);
        eprintln!(
            "joint directional velocity={ev:e} pressure={ep:e} total work={:e}",
            solution.relative_total_work_defect
        );
        assert!(ev < 2e-7 && ep < 2e-7);
        assert_eq!(solution.receipt.fmm_order, 8);
        assert!(solution.relative_total_work_defect < 1e-9);
        assert!(solution.structural_dissipation_w > 0.0);
    });
}

#[test]
fn original_coordinate_units_and_coupled_dark_motion_preserve_physical_fields() {
    let surface = SpherePanels::icosphere(0.1, 1).unwrap();
    let medium = Medium {
        density: 1.2,
        sound_speed: 3.0,
    };
    let m = [0.25, 0.02, 0.02, 0.2];
    let c = [0.1, 0.03, -0.01, 0.2];
    let k = [20.0, -2.0, -2.0, 25.0];
    let traces = vec![vec![1.0; 80], vec![0.0; 80]];
    let f = [C64::new(1.0, 0.2), C64::ZERO];
    let units = [1e-3, -7.0];
    let transform =
        |a: [f64; 4]| core::array::from_fn::<_, 4, _>(|i| a[i] * units[i / 2] * units[i % 2]);
    let (mm, cc, kk) = (transform(m), transform(c), transform(k));
    let tt = traces
        .iter()
        .zip(units)
        .map(|(r, s)| r.iter().map(|t| t * s).collect())
        .collect::<Vec<Vec<_>>>();
    let ff = [f[0].scale(units[0]), f[1].scale(units[1])];
    with_cx(|cx| {
        let gate = CancelGate::new_clock_free();
        let x = PreparedJointStructuralRadiation::new(
            &surface,
            &traces,
            pencil(&m, &c, &k, 2),
            8.0,
            medium,
            &options(),
            controls(),
            &gate,
            cx,
        )
        .unwrap()
        .solve(&f, &gate, cx)
        .unwrap();
        let y = PreparedJointStructuralRadiation::new(
            &surface,
            &tt,
            pencil(&mm, &cc, &kk, 2),
            8.0,
            medium,
            &options(),
            controls(),
            &gate,
            cx,
        )
        .unwrap()
        .solve(&ff, &gate, cx)
        .unwrap();
        let restored = y
            .velocity
            .iter()
            .zip(units)
            .map(|(v, s)| v.scale(s))
            .collect::<Vec<_>>();
        assert!(x.velocity[1].abs() > 1e-3);
        assert!(relative(&restored, &x.velocity) < 1e-8);
        assert!(relative(&y.field.pressure, &x.field.pressure) < 1e-8);
        assert!(
            (x.structural_dissipation_w - y.structural_dissipation_w).abs()
                < 1e-8 * x.structural_dissipation_w
        );
    });
}

#[test]
fn unresolved_dark_resonance_budgets_cancellation_and_reuse_remain_explicit() {
    let surface = SpherePanels::icosphere(0.1, 1).unwrap();
    let medium = Medium {
        density: 1.2,
        sound_speed: 3.0,
    };
    let traces = vec![vec![1.0; 80], vec![0.0; 80]];
    let m = [0.25, 0.0, 0.0, 0.25];
    let k = [17.0, 0.0, 0.0, 16.0];
    let c = [0.1, 0.0, 0.0, 0.0];
    with_cx(|cx| {
        let gate = CancelGate::new_clock_free();
        let owner = PreparedJointStructuralRadiation::new(
            &surface,
            &traces,
            pencil(&m, &c, &k, 2),
            8.0,
            medium,
            &options(),
            controls(),
            &gate,
            cx,
        )
        .unwrap();
        assert!(matches!(
            owner.solve(&[C64::ZERO, C64::ONE], &gate, cx),
            Err(JointStructuralRadiationError::NotConverged { .. })
        ));
        let force = [C64::new(1.0, 0.2), C64::ZERO];
        let first = owner.solve(&force, &gate, cx).unwrap();
        assert!(
            owner
                .solve_controlled(&force, &gate, cx, |count| count == 3)
                .is_err()
        );
        let retry = owner.solve(&force, &gate, cx).unwrap();
        assert_eq!(first.velocity, retry.velocity);
        assert_eq!(first.field.pressure, retry.field.pressure);
        let zero = owner.solve(&[C64::ZERO; 2], &gate, cx).unwrap();
        assert_eq!(zero.field.pressure, vec![C64::ZERO; 80]);
        assert_eq!(zero.velocity, vec![C64::ZERO; 2]);
        assert_eq!(zero.external_input_power_w, 0.0);
        assert!(
            owner
                .solve(&[C64::new(f64::NAN, 0.0), C64::ZERO], &gate, cx)
                .is_err()
        );
        let mut tiny = controls();
        tiny.maximum_payload_bytes = owner.added_payload_bytes() - 1;
        assert!(
            PreparedJointStructuralRadiation::new(
                &surface,
                &traces,
                pencil(&m, &c, &k, 2),
                8.0,
                medium,
                &options(),
                tiny,
                &gate,
                cx
            )
            .is_err()
        );
        let open =
            SpherePanels::from_triangles(surface.triangles().unwrap()[..79].to_vec()).unwrap();
        let incomplete_traces = vec![vec![1.0; 79], vec![0.0; 79]];
        assert!(
            PreparedJointStructuralRadiation::new(
                &open,
                &incomplete_traces,
                pencil(&m, &c, &k, 2),
                8.0,
                medium,
                &options(),
                controls(),
                &gate,
                cx,
            )
            .is_err()
        );
        let mut short = options();
        short.max_iters = 1;
        let insufficient = PreparedJointStructuralRadiation::new(
            &surface,
            &traces,
            pencil(&m, &c, &k, 2),
            8.0,
            medium,
            &short,
            controls(),
            &gate,
            cx,
        )
        .unwrap();
        assert!(matches!(
            insufficient.solve(&force, &gate, cx),
            Err(JointStructuralRadiationError::NotConverged { .. })
        ));
        let mut invalid_m = m;
        invalid_m[0] = -1.0;
        assert!(
            PreparedJointStructuralRadiation::new(
                &surface,
                &traces,
                pencil(&invalid_m, &c, &k, 2),
                8.0,
                medium,
                &options(),
                controls(),
                &gate,
                cx
            )
            .is_err()
        );
    });
}
