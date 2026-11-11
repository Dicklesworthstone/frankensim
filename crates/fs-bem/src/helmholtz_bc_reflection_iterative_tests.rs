use super::*;
use crate::helmholtz::bc::{IncidentSource, PlaneWave, PreparedReflectionBoundary};

const MEDIUM: Medium = Medium {
    density: 1.2,
    sound_speed: 340.0,
};

fn options() -> IterativeOptions {
    IterativeOptions {
        tol: 2e-11,
        restart: 60,
        max_iters: 120,
        fmm_order: 6,
        leaf_capacity: 64,
        near_field_cache_bytes: 0,
        precond: Preconditioner::NearFieldBlockJacobi {
            max_block_panels: 64,
        },
        ..IterativeOptions::default()
    }
}

fn relative(actual: &[C64], reference: &[C64]) -> f64 {
    let numerator: f64 = actual
        .iter()
        .zip(reference)
        .map(|(&a, &b)| (a - b).norm_sq())
        .sum();
    let denominator: f64 = reference.iter().map(|z| z.norm_sq()).sum();
    fs_math::det::sqrt(numerator / denominator.max(f64::MIN_POSITIVE))
}

fn incident(surface: &SpherePanels, k: f64) -> IncidentTraces {
    IncidentTraces::from_sources(
        surface,
        k,
        MEDIUM,
        &[IncidentSource::PlaneWave(PlaneWave {
            direction: [0.6, 0.0, 0.8],
            amplitude_pa: C64::new(0.7, -0.2),
        })],
        1e-4,
    )
    .unwrap()
}

fn original_equation_residual(surface: &SpherePanels, solution: &RadiationSolution) -> f64 {
    let n = surface.areas().len();
    let (a, b) = crate::helmholtz::assemble_dense_with_static(
        surface,
        solution.k,
        alpha_for(Formulation::BurtonMiller, solution.k),
        None,
    );
    let q_scale = C64::new(0.0, solution.k * MEDIUM.density * MEDIUM.sound_speed);
    let mut error = 0.0;
    let mut scale = 0.0;
    for i in 0..n {
        let mut left = C64::ZERO;
        let mut right = C64::ZERO;
        for j in 0..n {
            left = left + a[i * n + j] * solution.pressure[j];
            right = right + b[i * n + j] * q_scale * solution.velocity[j];
        }
        error += (left - right).norm_sq();
        scale += left.norm_sq() + right.norm_sq();
    }
    fs_math::det::sqrt(error / scale.max(f64::MIN_POSITIVE))
}

#[test]
fn mixed_reflection_matches_dense_and_the_original_boundary_equation() {
    let surface = SpherePanels::icosphere(0.1, 0).unwrap();
    let n = surface.areas().len();
    let reflection: Vec<_> = (0..n)
        .map(|i| match i % 4 {
            0 => C64::ONE,
            1 => -C64::ONE,
            2 => C64::new(0.3, -0.4),
            _ => C64::new(0.6, 0.8),
        })
        .collect();
    let gate = CancelGate::new_clock_free();
    let dense = PreparedReflectionBoundary::new(
        &surface,
        3.0,
        MEDIUM,
        &reflection,
        false,
        Formulation::BurtonMiller,
    )
    .unwrap();
    let prepared = PreparedIterativeReflectionBoundary::new(
        &surface,
        3.0,
        MEDIUM,
        &reflection,
        false,
        Formulation::BurtonMiller,
        &options(),
        &gate,
    )
    .unwrap();
    let drive: Vec<_> = (0..n)
        .map(|i| C64::new(0.01 * (i + 1) as f64, -0.02))
        .collect();
    let expected = dense.solve(&drive).unwrap();
    let (actual, report, receipt) = prepared.solve(&drive, &gate).unwrap();
    assert!(report.converged_euclidean() && report.rel_residual <= options().tol);
    assert!(receipt.matvecs > 0);
    assert!(relative(&actual.radiation.pressure, &expected.radiation.pressure) < 2e-9);
    assert!(relative(&actual.radiation.velocity, &expected.radiation.velocity) < 2e-9);
    assert!(original_equation_residual(&surface, &actual.radiation) < 2e-10);
    assert!((actual.absorbed_power_w / expected.absorbed_power_w - 1.0).abs() < 2e-9);
    let incident = incident(&surface, 3.0);
    let expected = dense.scatter(&incident).unwrap();
    let (actual, report, _) = prepared.scatter(&incident, &gate).unwrap();
    assert!(report.converged_euclidean());
    assert!(relative(&actual.scattered.pressure, &expected.scattered.pressure) < 2e-9);
    assert!(relative(&actual.scattered.velocity, &expected.scattered.velocity) < 2e-9);
    assert!(original_equation_residual(&surface, &actual.scattered) < 2e-10);
    assert!((actual.absorbed_power_w / expected.absorbed_power_w - 1.0).abs() < 2e-9);
    for (i, &r) in reflection.iter().enumerate() {
        if r == -C64::ONE {
            assert_eq!(actual.total_pressure_pa()[i], C64::ZERO);
        }
    }
}

#[test]
fn pressure_release_uses_its_own_flow_block_and_retains_zero_wall_work() {
    let surface = SpherePanels::icosphere(0.1, 0).unwrap();
    let gate = CancelGate::new_clock_free();
    let reflection = vec![-C64::ONE; surface.areas().len()];
    let prepared = PreparedIterativeReflectionBoundary::new(
        &surface,
        3.0,
        MEDIUM,
        &reflection,
        false,
        Formulation::BurtonMiller,
        &options(),
        &gate,
    )
    .unwrap();
    let (response, report, _) = prepared.scatter(&incident(&surface, 3.0), &gate).unwrap();
    assert!(report.converged_euclidean());
    // The preconditioner blocks are fs-fmm's leaf groups, and the tree has at least the eight octants of depth
    // 1, so this 20-panel sphere is several blocks and the transformed block-Jacobi is not the complete operator:
    // no iteration count is asserted. The solve converges to its own stated tolerance.
    assert!(
        report.rel_residual <= options().tol,
        "relative residual {:e} above the solve's tolerance {:e}",
        report.rel_residual,
        options().tol
    );
    assert!(response.total_pressure_pa().iter().all(|&p| p == C64::ZERO));
    assert_eq!(response.absorbed_power_w, 0.0);
    assert!(response.scattered.velocity.iter().any(|v| v.abs() > 1e-6));
    assert!(original_equation_residual(&surface, &response.scattered) < 2e-10);
}

#[test]
fn transformed_blocks_follow_each_existing_burton_miller_panel_rule() {
    let gate = CancelGate::new_clock_free();
    // The Triangle rule admits a cloud only if its leaf box spans CBIE_MIN_LEAF_RADII = 4 panel radii
    // (helmholtz_iterative.rs). The tree depth is clamped to at least 1, so a leaf is at most half the body, 0.1
    // here, and no leaf_capacity enlarges it. icosphere(0.1, 1) has panel radius 0.03568 and 4r = 0.1427, which is
    // refused; icosphere(0.1, 2) has radius 0.01876 and 4r = 0.0750, which is admitted. The Triangle arm therefore
    // runs on level 2; the other two arms stay on level 1. The Triangle rule's own ceiling, `max_pairs`, bounds
    // the near pairs it integrates exactly: the ordered pairs i != j with |c_i - c_j| <= 4 r_j
    // (helmholtz_iterative.rs, triangle_is_near). That count is a property of the mesh, 380, 3,600 and 12,840 on
    // icosphere(0.1, 0), (0.1, 1) and (0.1, 2), counted independently, so level 1 sits under 10000 and the
    // level-2 Triangle arm needs 12,840.
    for (rule, level) in [
        (PanelRule::Centroid, 1),
        (PanelRule::TriangleNear { max_pairs: 10000 }, 1),
        (PanelRule::Triangle { max_pairs: 12840 }, 2),
    ] {
        let surface = SpherePanels::icosphere(0.1, level).unwrap();
        let reflection: Vec<_> = (0..surface.areas().len())
            .map(|i| match i % 3 {
                0 => C64::ONE,
                1 => -C64::ONE,
                _ => C64::new(0.2, 0.5),
            })
            .collect();
        let mut config = options();
        config.leaf_capacity = 128;
        config.precond = Preconditioner::NearFieldBlockJacobi {
            max_block_panels: 128,
        };
        config.panel_rule = rule;
        let prepared = PreparedIterativeReflectionBoundary::new(
            &surface,
            3.0,
            MEDIUM,
            &reflection,
            false,
            Formulation::BurtonMiller,
            &config,
            &gate,
        )
        .unwrap();
        let ReflectionFactor::Local(operator) = &prepared.factor else {
            panic!("mixed boundary");
        };
        // The blocks are fs-fmm's leaf groups (at least the eight octants of depth 1), so no mesh is one block.
        // Each panel is in exactly one block, and the exactness of the transform is per block: for x supported
        // on one block, every pair inside a leaf is a direct near-field pair, so the image restricted to that
        // block is the block's own transformed matrix times x, and its stored inverse returns x there.
        assert!(!operator.blocks.is_empty());
        let mut covered = vec![0usize; reflection.len()];
        for block in &operator.blocks {
            for &i in &block.indices {
                covered[i] += 1;
            }
        }
        assert!(covered.iter().all(|&count| count == 1), "{rule:?}");
        for block in &operator.blocks {
            let mut x = vec![C64::ZERO; reflection.len()];
            for &i in &block.indices {
                x[i] = C64::new(0.03 * (i + 1) as f64, -0.04);
            }
            let pressure: Vec<_> = x
                .iter()
                .zip(&prepared.pressure_coefficient)
                .map(|(&x, &a)| a * x)
                .collect();
            let flow: Vec<_> = x
                .iter()
                .zip(&prepared.flow_coefficient)
                .map(|(&x, &b)| b * x)
                .collect();
            let actual = prepared.image(operator, &pressure, &flow, &gate).unwrap();
            let recovered = operator.inverse(&actual);
            let own: Vec<_> = block.indices.iter().map(|&i| recovered[i]).collect();
            let expected: Vec<_> = block.indices.iter().map(|&i| x[i]).collect();
            assert!(relative(&own, &expected) < 3e-10, "{rule:?}");
        }
    }
}

#[test]
fn rigid_arithmetic_and_actual_failure_diagnostics_are_preserved() {
    let surface = SpherePanels::icosphere(0.1, 0).unwrap();
    let n = surface.areas().len();
    let gate = CancelGate::new_clock_free();
    let drive = vec![C64::new(0.01, -0.02); n];
    let old = PreparedIterativeBoundary::new(
        &surface,
        3.0,
        MEDIUM,
        Admittance::Rigid,
        Formulation::BurtonMiller,
        &options(),
        &gate,
    )
    .unwrap()
    .solve(&drive, &gate)
    .unwrap()
    .0;
    let rigid = PreparedIterativeReflectionBoundary::new(
        &surface,
        3.0,
        MEDIUM,
        &vec![C64::ONE; n],
        false,
        Formulation::BurtonMiller,
        &options(),
        &gate,
    )
    .unwrap();
    let actual = rigid.solve(&drive, &gate).unwrap().0;
    assert_eq!(actual.radiation.pressure, old.radiation.pressure);
    assert_eq!(actual.radiation.velocity, old.radiation.velocity);
    let mut stopped = options();
    stopped.max_iters = 0;
    let release = PreparedIterativeReflectionBoundary::new(
        &surface,
        3.0,
        MEDIUM,
        &vec![-C64::ONE; n],
        false,
        Formulation::BurtonMiller,
        &stopped,
        &gate,
    )
    .unwrap();
    let incident = incident(&surface, 3.0);
    match release.scatter(&incident, &gate) {
        Err(IterativeError::NotConverged {
            pressure,
            report,
            receipt,
        }) => {
            assert!(!report.converged_euclidean());
            assert_eq!(receipt.iterations, 0);
            for (&p, &incoming) in pressure.iter().zip(incident.pressure_pa()) {
                assert_eq!(p, -incoming);
            }
        }
        _ => panic!("a pressure trace alone cannot certify the unsolved flow"),
    }
    let mut budget = options();
    budget.max_near_entries = 0;
    assert!(matches!(
        PreparedIterativeReflectionBoundary::new(
            &surface,
            3.0,
            MEDIUM,
            &vec![-C64::ONE; n],
            false,
            Formulation::BurtonMiller,
            &budget,
            &gate,
        ),
        Err(IterativeError::NearFieldBudget { .. })
    ));
    for r in [C64::from_re(1.01), C64::new(f64::NAN, 0.0)] {
        assert!(
            PreparedIterativeReflectionBoundary::new(
                &surface,
                3.0,
                MEDIUM,
                &vec![r; n],
                false,
                Formulation::BurtonMiller,
                &options(),
                &gate,
            )
            .is_err()
        );
    }
    gate.request();
    assert!(matches!(
        rigid.solve(&drive, &gate),
        Err(IterativeError::Cancelled { .. })
    ));
}
