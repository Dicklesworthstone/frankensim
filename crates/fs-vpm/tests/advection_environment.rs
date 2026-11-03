//! G2 incident-flow/image oracles, G3 Galilean covariance, G4 atomic refusals.

use fs_vpm::farfield3d::{FarField, WakeCoreEvolutionMode};
use fs_vpm::filament3d::{
    FilamentWake, Refusal, SelfInducedShed, WakeBoundary, WakeRateCertificate,
};
use fs_vpm::images3d::{ImagePlane, ImageSystem};

const CORE: f64 = 0.2;
const SPAN: f64 = 1.4;
const DT: f64 = 0.5;

fn ring(gamma: f64) -> FilamentWake {
    let mut wake = FilamentWake::new(
        WakeRateCertificate {
            shed_hz: 1.0 / DT,
            n_stations: 1,
            max_rows: 64,
        },
        vec![[0.0, -SPAN / 2.0, 1.0], [0.0, SPAN / 2.0, 1.0]],
    )
    .unwrap();
    wake.shed(&[gamma], [0.0; 3]).unwrap();
    for p in &mut wake.rows[0].nodes {
        p[0] = -1.0;
    }
    wake.rows[0].core2_m2 = CORE * CORE;
    wake
}

fn config(substeps: usize) -> SelfInducedShed {
    SelfInducedShed {
        ambient_velocity_mps: [0.0; 3],
        emitted_core_radius_m: CORE,
        core_mode: WakeCoreEvolutionMode::Frozen,
        substeps,
        max_displacement_fraction: 1.0,
        max_work: 10_000_000,
    }
}

fn distance(a: [f64; 3], b: [f64; 3]) -> f64 {
    (a[0] - b[0]).hypot(a[1] - b[1]).hypot(a[2] - b[2])
}

fn assert_close(a: [f64; 3], b: [f64; 3], tolerance: f64) {
    assert!(distance(a, b) <= tolerance, "{a:?} != {b:?}");
}

#[test]
fn explicit_default_environment_is_bit_identical_to_the_legacy_entry_point() {
    let mut legacy = ring(1.0);
    let mut explicit = legacy.clone();
    let mut settings = config(16);
    settings.core_mode = WakeCoreEvolutionMode::CoreSpreading { nu_eff_m2ps: 0.01 };
    settings.ambient_velocity_mps = [-0.1, 0.0, 0.05];
    for gamma in [1.5, 0.8] {
        legacy
            .shed_self_induced(&[gamma], settings, || false)
            .unwrap();
        explicit
            .shed_self_induced_with_environment(
                &[gamma],
                settings,
                WakeBoundary::default(),
                |_, _| Ok([-0.0, 0.0, -0.0]),
                || false,
            )
            .unwrap();
        assert_eq!(explicit, legacy);
    }
}

#[test]
fn spatial_incident_flow_converges_to_the_independent_affine_exponential_solution() {
    let original = ring(0.0);
    let (lambda, offset, decay) = (0.6_f64, 0.15_f64, 0.2_f64);
    let exact = original.rows[0].nodes[0];
    let exact = [
        (exact[0] + offset / lambda) * (lambda * DT).exp() - offset / lambda,
        exact[1] * (-decay * DT).exp(),
        exact[2],
    ];
    let errors: Vec<_> = [4, 8, 16]
        .into_iter()
        .map(|substeps| {
            let mut wake = original.clone();
            wake.shed_self_induced_with_environment(
                &[7.0],
                config(substeps),
                WakeBoundary::default(),
                |p, _| Ok([lambda * p[0] + offset, -decay * p[1], 0.0]),
                || false,
            )
            .unwrap();
            // The endpoint strength cannot affect this initially zero-Γ interval.
            assert_eq!(wake.rows[0].gamma, [0.0]);
            assert_eq!(wake.rows[1].gamma, [7.0]);
            distance(wake.rows[0].nodes[0], exact)
        })
        .collect();
    for pair in errors.windows(2) {
        assert!((3.7..4.1).contains(&(pair[0] / pair[1])), "{errors:?}");
    }
}

#[test]
fn time_dependent_incident_acceleration_uses_elapsed_seconds_at_every_stage() {
    let mut wake = ring(0.0);
    let before = wake.clone();
    let mut settings = config(3);
    settings.ambient_velocity_mps = [0.3, -0.2, 0.1];
    let acceleration = [0.4, 0.2, -0.3];
    let mut times = Vec::new();
    wake.shed_self_induced_with_environment(
        &[9.0],
        settings,
        WakeBoundary::default(),
        |_, elapsed| {
            times.push(elapsed);
            Ok(acceleration.map(|a| a * elapsed))
        },
        || false,
    )
    .unwrap();
    let expected_times: Vec<_> = (0..3)
        .flat_map(|step| {
            [step as f64 / 3.0 * DT, (step as f64 + 0.5) / 3.0 * DT]
                .into_iter()
                .flat_map(|time| [time; 2])
        })
        .collect();
    assert_eq!(times, expected_times);
    for (&a, &b) in before.rows[0].nodes.iter().zip(&wake.rows[0].nodes) {
        assert_close(
            b,
            std::array::from_fn(|i| {
                a[i] + settings.ambient_velocity_mps[i] * DT + 0.5 * acceleration[i] * DT * DT
            }),
            3e-15,
        );
    }
    assert_eq!(wake.rows[1].nodes, before.line_nodes);
    assert_eq!(wake.rows[1].age_ticks, 0.0);
}

#[test]
fn prescribed_translation_and_matching_ambient_flow_preserve_the_anchored_ring_motion() {
    let mut stationary = ring(4.0);
    let mut moving = stationary.clone();
    let drift = [0.2, -0.1, 0.15];
    let mut moving_settings = config(64);
    moving_settings.ambient_velocity_mps = drift;
    for interval in 1..=2 {
        let endpoint: Vec<_> = moving
            .line_nodes
            .iter()
            .map(|p| std::array::from_fn(|i| p[i] + DT * drift[i]))
            .collect();
        stationary
            .shed_self_induced(&[9.0], config(64), || false)
            .unwrap();
        moving
            .shed_self_induced_with_environment(
                &[9.0],
                moving_settings,
                WakeBoundary {
                    line_endpoint_m: Some(&endpoint),
                    image_plane: None,
                },
                |_, _| Ok([0.0; 3]),
                || false,
            )
            .unwrap();
        for (fixed, translated) in stationary.rows.iter().zip(&moving.rows) {
            for (&a, &b) in fixed.nodes.iter().zip(&translated.nodes) {
                assert_close(
                    b,
                    std::array::from_fn(|i| a[i] + interval as f64 * DT * drift[i]),
                    3e-14,
                );
            }
        }
        assert_eq!(moving.line_nodes, endpoint);
        assert_eq!(moving.rows.last().unwrap().nodes, endpoint);
        assert_eq!(moving.rows.last().unwrap().age_ticks, 0.0);
        if interval == 1 {
            let angle = 4.0 * SPAN * DT
                / (4.0 * std::f64::consts::PI * (1.0 + SPAN * SPAN + CORE * CORE).sqrt())
                * (1.0 / (1.0 + CORE * CORE) + 1.0 / (SPAN * SPAN + CORE * CORE));
            for p in &moving.rows[0].nodes {
                assert_close(
                    *p,
                    [
                        -angle.cos() + DT * drift[0],
                        p[1],
                        1.0 + angle.sin() + DT * drift[2],
                    ],
                    2e-7,
                );
            }
            assert_eq!(moving.rows[0].gamma, [4.0]);
            assert_eq!(moving.rows[1].gamma, [9.0]);
        }
    }
    // The row born on the line at the first endpoint is free during the next
    // interval, even though its initial coordinates match the prescribed line.
    assert!(distance(moving.rows[1].nodes[0], moving.line_nodes[0]) > 1e-3);
}

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

// Independent quadrature of Γ/(4π) ∫ dl×r / (|r|²+a²)^(3/2).
// No production segment kernel, source traversal, image builder or RK helper.
fn line_quadrature(p: [f64; 3], a: [f64; 3], b: [f64; 3], gamma: f64) -> [f64; 3] {
    const N: usize = 1024;
    let dl = std::array::from_fn(|i| b[i] - a[i]);
    let mut result = [0.0; 3];
    for j in 0..=N {
        let r: [f64; 3] = std::array::from_fn(|i| p[i] - a[i] - j as f64 / N as f64 * dl[i]);
        let rho2 = r.iter().map(|v| v * v).sum::<f64>() + CORE * CORE;
        let weight = if j == 0 || j == N {
            1.0
        } else if j % 2 == 0 {
            2.0
        } else {
            4.0
        };
        let numerator = cross(dl, r);
        for i in 0..3 {
            result[i] += weight * numerator[i] / (rho2 * rho2.sqrt());
        }
    }
    result.map(|v| gamma * v / (12.0 * std::f64::consts::PI * N as f64))
}

fn hand_image_ring(p: [f64; 3], line: [[f64; 3]; 2], free: [[f64; 3]; 2]) -> [f64; 3] {
    // Orientation follows one positive material quad, independent of the
    // production lattice's shared-edge traversal and summation order.
    let vertices = [line[0], line[1], free[1], free[0]];
    let mut velocity = [0.0; 3];
    for edge in 0..4 {
        let a = vertices[edge];
        let b = vertices[(edge + 1) % 4];
        let physical = line_quadrature(p, a, b, 1.0);
        let reflected = line_quadrature(p, [a[0], a[1], -a[2]], [b[0], b[1], -b[2]], -1.0);
        for i in 0..3 {
            velocity[i] += physical[i] + reflected[i];
        }
    }
    velocity
}

#[test]
fn image_evolution_matches_independent_two_stage_line_quadrature_and_slip_condition() {
    let mut wake = ring(1.0);
    wake.cert.shed_hz = 100.0;
    let dt = 0.01;
    let line = [wake.line_nodes[0], wake.line_nodes[1]];
    let initial = [wake.rows[0].nodes[0], wake.rows[0].nodes[1]];
    let ambient = [-0.2, 0.03, 0.0];
    let midpoint: [[f64; 3]; 2] = initial.map(|p| {
        let u = hand_image_ring(p, line, initial);
        std::array::from_fn(|i| p[i] + 0.5 * dt * (u[i] + ambient[i]))
    });
    let expected: [[f64; 3]; 2] = std::array::from_fn(|j| {
        let u = hand_image_ring(midpoint[j], line, midpoint);
        std::array::from_fn(|i| initial[j][i] + dt * (u[i] + ambient[i]))
    });
    let mut settings = config(1);
    settings.ambient_velocity_mps = ambient;
    wake.shed_self_induced_with_environment(
        &[3.0],
        settings,
        WakeBoundary {
            line_endpoint_m: None,
            image_plane: Some(ImagePlane { z0_m: 0.0 }),
        },
        |_, _| Ok([0.0; 3]),
        || false,
    )
    .unwrap();
    for (actual, expected) in wake.rows[0].nodes.iter().zip(expected) {
        assert_close(*actual, expected, 3e-12);
    }
    let far = FarField {
        mode: WakeCoreEvolutionMode::Frozen,
        shed_dt_s: dt,
        cells: Vec::new(),
    };
    let image = ImageSystem::build(&wake, &far, ImagePlane { z0_m: 0.0 }).unwrap();
    for p in [[-0.2, 0.1, 0.0], [0.3, -0.4, 0.0], [-1.3, 1.1, 0.0]] {
        let physical = wake.induced_velocity(p);
        let reflected = image.eval(p);
        assert!((physical[2] + reflected[2]).abs() < 2e-14);
    }
    assert_eq!(wake.rows.len(), 2);
    assert_eq!(wake.rows[0].gamma, [1.0]);
    assert_eq!(wake.rows[1].gamma, [3.0]);
}

#[test]
fn both_midpoint_and_endpoint_plane_crossings_refuse_without_publication() {
    for (vertical, expected_calls) in [(-0.5, 1), (-0.3, 3), (-0.2, 3)] {
        let mut wake = ring(0.0);
        for p in &mut wake.rows[0].nodes {
            p[2] = 0.1;
        }
        let before = wake.clone();
        let mut calls = 0;
        let error = wake
            .shed_self_induced_with_environment(
                &[1.0],
                config(1),
                WakeBoundary {
                    line_endpoint_m: None,
                    image_plane: Some(ImagePlane { z0_m: 0.0 }),
                },
                |_, _| {
                    calls += 1;
                    Ok([0.0, 0.0, vertical])
                },
                || false,
            )
            .unwrap_err();
        assert_eq!(error.code, "image-plane-invalid");
        assert_eq!(calls, expected_calls);
        assert_eq!(wake, before);
    }
}

fn unavailable() -> Refusal {
    Refusal {
        code: "test-incident-unavailable",
        message: "incident flow is unavailable at this stage".into(),
        ranked_repairs: vec!["supply the required physical field".into()],
    }
}

#[test]
fn nonfinite_and_early_or_late_incident_errors_roll_back_all_geometry_and_ages() {
    let original = ring(0.2);
    for case in 0..4 {
        let mut wake = original.clone();
        let mut calls = 0;
        let error = wake
            .shed_self_induced_with_environment(
                &[2.0],
                config(4),
                WakeBoundary::default(),
                |_, elapsed| {
                    calls += 1;
                    match case {
                        0 => Err(unavailable()),
                        1 => Ok([f64::NAN, 0.0, 0.0]),
                        2 => Ok([0.0, f64::INFINITY, 0.0]),
                        _ if elapsed >= 0.75 * DT => Err(unavailable()),
                        _ => Ok([0.02, 0.01, 0.0]),
                    }
                },
                || false,
            )
            .unwrap_err();
        if case == 0 || case == 3 {
            assert_eq!(error, unavailable());
        } else {
            assert_eq!(error.code, "wake-advection-nonfinite");
        }
        if case == 3 {
            assert_eq!(calls, 13);
        }
        assert_eq!(wake, original);
    }
}

#[test]
fn boundary_admission_line_displacement_and_image_work_limits_are_atomic() {
    let original = ring(0.2);
    for case in 0..5 {
        let mut wake = original.clone();
        let mut endpoint = original.line_nodes.clone();
        let mut settings = config(1);
        let mut plane = None;
        let expected = match case {
            0 => {
                endpoint.pop();
                "wake-advection-line-invalid"
            }
            1 => {
                endpoint[0][0] = f64::NAN;
                "wake-advection-line-invalid"
            }
            2 => {
                endpoint[0][2] = 0.0;
                plane = Some(ImagePlane { z0_m: 0.0 });
                "wake-advection-line-invalid"
            }
            3 => {
                endpoint[0][0] += 0.3;
                "wake-advection-displacement"
            }
            _ => {
                settings.max_work = 62;
                plane = Some(ImagePlane { z0_m: 0.0 });
                "wake-advection-work-exhausted"
            }
        };
        let error = wake
            .shed_self_induced_with_environment(
                &[1.0],
                settings,
                WakeBoundary {
                    line_endpoint_m: Some(&endpoint),
                    image_plane: plane,
                },
                |_, _| Ok([0.0; 3]),
                || false,
            )
            .unwrap_err();
        assert_eq!(error.code, expected, "case {case}");
        assert_eq!(wake, original);
    }
    // 2 stages × 3 (targets+1) × 9 (8 image+physical sources+1),
    // plus 9 for final source-geometry admission.
    let mut at_cap = original.clone();
    let mut settings = config(1);
    settings.max_work = 63;
    at_cap
        .shed_self_induced_with_environment(
            &[1.0],
            settings,
            WakeBoundary {
                line_endpoint_m: None,
                image_plane: Some(ImagePlane { z0_m: 0.0 }),
            },
            |_, _| Ok([0.0; 3]),
            || false,
        )
        .unwrap();
    settings.max_work = 35;
    let mut without_images = original;
    without_images
        .shed_self_induced(&[1.0], settings, || false)
        .unwrap();
}

#[test]
fn final_resolved_edge_collapse_refuses_and_finer_steps_resolve_incompressible_strain() {
    // The divergence-free field (2*t*x, -2*t*y, 0) has the exact trajectory
    // y(t) = y(0)*exp(-t²). Tiny nonzero Γ keeps the physical span active while
    // its self-induced displacement is below coordinate roundoff here.
    let mut original = ring(1e-30);
    original.cert.shed_hz = 1.0;
    original.line_nodes = vec![[1.0, -1.0, 2.0], [1.0, 1.0, 2.0]];
    original.rows[0].nodes = vec![[0.0, -1.0, 2.0], [0.0, 1.0, 2.0]];
    original.rows[0].core2_m2 = 4.0;
    let mut settings = config(1);
    settings.emitted_core_radius_m = 2.0;
    let mut wake = original.clone();
    let error = wake
        .shed_self_induced_with_environment(
            &[1e-30],
            settings,
            WakeBoundary::default(),
            |p, t| Ok([2.0 * t * p[0], -2.0 * t * p[1], 0.0]),
            || false,
        )
        .unwrap_err();
    // With one step both sampled source geometries remain regular, but the
    // final RK update sends both span endpoints to (0, 0, 2). It must refuse.
    assert_eq!(error.code, "wake-advection-core-invalid");
    assert_eq!(wake, original);

    let exact = (-1.0_f64).exp();
    let mut errors = Vec::new();
    for substeps in [16, 32, 64] {
        let mut wake = original.clone();
        settings.substeps = substeps;
        wake.shed_self_induced_with_environment(
            &[1e-30],
            settings,
            WakeBoundary::default(),
            |p, t| Ok([2.0 * t * p[0], -2.0 * t * p[1], 0.0]),
            || false,
        )
        .unwrap();
        errors.push((wake.rows[0].nodes[1][1] - exact).abs());
        if substeps == 64 {
            for (p, sign) in wake.rows[0].nodes.iter().zip([-1.0, 1.0]) {
                assert_close(*p, [0.0, sign * exact, 2.0], 2e-5);
            }
        }
    }
    for pair in errors.windows(2) {
        assert!((3.8..4.5).contains(&(pair[0] / pair[1])), "{errors:?}");
    }
}
