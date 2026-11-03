//! Normal-section lifting-line equations, frame covariance and physical loads.
use fs_ga::{So3, So3Tangent, Vec3};
use fs_wing::nonlinear::continuation::{
    ContinuationControls, solve_with_incident, solve_with_section_frames,
};
use fs_wing::nonlinear::{
    InfluenceOperator, SectionFrame, StripRegime, StripSpec, solve_nonlinear_with_section_frames,
};
use fs_wing::{Panel, SurfaceId, flat_surface, induced_velocity_free};

const RHO: f64 = 1.225;

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a.into_iter().zip(b).map(|(x, y)| x * y).sum()
}
fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    core::array::from_fn(|i| a[i] - b[i])
}
fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
fn rotate(r: So3, a: [f64; 3]) -> [f64; 3] {
    let value = r.rotate(Vec3::new(a[0], a[1], a[2])).unwrap();
    [value.x, value.y, value.z]
}
fn rotation(a: [f64; 3]) -> So3 {
    So3::exp(So3Tangent::new(Vec3::new(a[0], a[1], a[2]))).unwrap()
}
fn close(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() < tolerance, "{a:e} != {b:e}");
}
fn vector_close(a: [f64; 3], b: [f64; 3], tolerance: f64) {
    for i in 0..3 {
        close(a[i], b[i], tolerance);
    }
}
fn layout(rows: usize) -> (Vec<Panel>, Vec<StripSpec>) {
    let panels = flat_surface(SurfaceId::WingLower, 6.0, 1.3, 0.0, 0.0, 6, rows).unwrap();
    let strips = (0..6)
        .map(|s| StripSpec {
            panel_indices: (0..rows).map(|row| row * 6 + s).collect(),
            chord_m: 1.3,
            twist_rad: 0.0,
        })
        .collect();
    (panels, strips)
}
fn frames(panels: &[Panel], strips: &[StripSpec]) -> Vec<SectionFrame> {
    strips
        .iter()
        .map(|strip| SectionFrame::from_panel(&panels[strip.panel_indices[0]]).unwrap())
        .collect()
}
fn section(_: usize, alpha: f64) -> (f64, StripRegime) {
    (core::f64::consts::TAU * alpha, StripRegime::Attached)
}
fn controls() -> ContinuationControls {
    ContinuationControls {
        maximum_incidence_rad: 0.7,
        maximum_incidence_change_rad: 0.1,
        maximum_row_condition: 8.0,
        minimum_trailing_clearance_chords: 0.02,
        maximum_newton_iterations: 32,
        maximum_residual_evaluations: 600,
        maximum_jacobian_products: 4096,
        section_difference_rad: 1e-6,
    }
}

#[test]
fn nonlinear_sections_follow_vertical_canted_and_reversed_world_axes() {
    let (panels, strips) = layout(2);
    let original_frames = frames(&panels, &strips);
    let rotations = [
        rotation([core::f64::consts::FRAC_PI_2, 0.0, 0.0]),
        rotation([0.35, -0.7, 0.4]),
        rotation([0.0, core::f64::consts::PI, 0.0]),
    ];
    for sign in [-1.0, 1.0] {
        let free = [12.0, 0.3, sign * 1.1];
        let incident: Vec<[f64; 3]> = (0..strips.len())
            .map(|s| [12.0 + 0.15 * s as f64, 0.3, sign * (0.8 + 0.1 * s as f64)])
            .collect();
        let op = InfluenceOperator::build(&panels, free, RHO).unwrap();
        let original = solve_nonlinear_with_section_frames(
            &op,
            &panels,
            &strips,
            free,
            RHO,
            &section,
            &original_frames,
            &incident,
            None,
        )
        .unwrap();
        for (index, r) in rotations.into_iter().enumerate() {
            let offset = [1.7, -2.3, 0.8];
            let moved: Vec<Panel> = panels
                .iter()
                .map(|p| {
                    let point = |x| {
                        let value = rotate(r, x);
                        core::array::from_fn(|i| value[i] + offset[i])
                    };
                    Panel {
                        a: point(p.a),
                        b: point(p.b),
                        ctrl: point(p.ctrl),
                        normal: rotate(r, p.normal),
                        ..*p
                    }
                })
                .collect();
            let free = rotate(r, free);
            let local: Vec<[f64; 3]> = incident.iter().map(|v| rotate(r, *v)).collect();
            let local_frames = frames(&moved, &strips);
            let op = InfluenceOperator::build(&moved, free, RHO).unwrap();
            let result = solve_nonlinear_with_section_frames(
                &op,
                &moved,
                &strips,
                free,
                RHO,
                &section,
                &local_frames,
                &local,
                None,
            )
            .unwrap();
            for (&g, &reference) in result.gamma.iter().zip(&original.gamma) {
                close(g, reference, 3e-10);
            }
            for (&alpha, &reference) in result.alpha_eff_rad.iter().zip(&original.alpha_eff_rad) {
                close(alpha, reference, 2e-12);
            }
            vector_close(
                result.total_force_n,
                rotate(r, original.total_force_n),
                2e-8,
            );
            assert_eq!(result.branch_id, original.branch_id);
            if index == 0 {
                // The vertical surface produces a restoring sideforce for
                // either lateral incident-flow sign, not a spurious z lift.
                assert!(result.total_force_n[1] * free[1] < 0.0);
                assert!(result.total_force_n[1].abs() > 10.0);
                assert!(result.total_force_n[2].abs() < 1e-10);
            }
            if index == 2 {
                assert!(free[0] < 0.0);
            }
        }
    }
}

#[test]
fn swept_normal_section_area_gives_q_s_cl_and_excludes_spanwise_speed() {
    let sweep: f64 = 0.4;
    let streamwise_chord = 1.3;
    let projected_span = 6.0;
    let mut panels = Vec::new();
    let mut strips = Vec::new();
    for i in 0..6 {
        let y0 = -3.0 + i as f64;
        let y1 = y0 + 1.0;
        let xm = sweep.tan() * 0.5 * (y0 + y1);
        panels.push(Panel {
            surface: SurfaceId::WingLower,
            a: [sweep.tan() * y0, y0, 0.0],
            b: [sweep.tan() * y1, y1, 0.0],
            ctrl: [xm - 0.5 * streamwise_chord, 0.5 * (y0 + y1), 0.0],
            normal: [0.0, 0.0, -1.0],
            width_m: 1.0 / sweep.cos(),
        });
        strips.push(StripSpec {
            panel_indices: vec![i],
            chord_m: streamwise_chord * sweep.cos(),
            twist_rad: 0.0,
        });
    }
    let frames = frames(&panels, &strips);
    let free = [12.0, 0.0, 1.0];
    let op = InfluenceOperator::build(&panels, free, RHO).unwrap();
    let span = frames[0].span_direction();
    let span_gust = core::array::from_fn(|i| free[i] + 30.0 * span[i]);
    for cl in [-0.7, 0.7] {
        let closure = |_: usize, _: f64| (cl, StripRegime::Attached);
        let mut reference = None;
        for velocity in [free, span_gust] {
            let local = vec![velocity; strips.len()];
            let result = solve_nonlinear_with_section_frames(
                &op, &panels, &strips, free, RHO, &closure, &frames, &local, None,
            )
            .unwrap();
            let normal_speed_squared = (free[0] * sweep.cos()).powi(2) + free[2] * free[2];
            let expected_magnitude =
                0.5 * RHO * normal_speed_squared * streamwise_chord * projected_span * cl.abs();
            close(
                dot(result.total_force_n, result.total_force_n).sqrt(),
                expected_magnitude,
                2e-8,
            );
            assert!(result.total_force_n[2] * cl < 0.0);
            for (&gamma, strip) in result.gamma.iter().zip(&strips) {
                close(
                    gamma,
                    -0.5 * normal_speed_squared.sqrt() * strip.chord_m * cl,
                    2e-10,
                );
            }
            close(dot(result.total_force_n, velocity), 0.0, 1e-9);
            if let Some(old) = reference {
                vector_close(result.total_force_n, old, 2e-8);
            } else {
                reference = Some(result.total_force_n);
            }
        }
    }
}

/// Independent complete-horseshoe field minus finite bound-vortex quadrature.
fn trailing(point: [f64; 3], panels: &[Panel], gamma: &[f64], free: [f64; 3]) -> [f64; 3] {
    let mut result = induced_velocity_free(point, panels, gamma, free);
    for (p, &g) in panels.iter().zip(gamma) {
        let r1 = sub(point, p.a);
        let r2 = sub(point, p.b);
        let segment = sub(p.b, p.a);
        let product = cross(r1, r2);
        let squared = dot(product, product);
        let n1 = dot(r1, r1).sqrt();
        let n2 = dot(r2, r2).sqrt();
        if squared >= 1e-8 && n1 >= 1e-8 && n2 >= 1e-8 {
            let scale = g * (dot(segment, r1) / n1 - dot(segment, r2) / n2)
                / (4.0 * core::f64::consts::PI * squared);
            for i in 0..3 {
                result[i] -= product[i] * scale;
            }
        }
    }
    result
}

#[test]
fn induced_incidence_is_the_directional_angle_derivative_at_finite_incidence() {
    let (panels, strips) = layout(2);
    let frames = frames(&panels, &strips);
    let free = [10.0, 0.0, 2.5];
    let local = vec![free; strips.len()];
    let op = InfluenceOperator::build(&panels, free, RHO).unwrap();
    let result = solve_nonlinear_with_section_frames(
        &op, &panels, &strips, free, RHO, &section, &frames, &local, None,
    )
    .unwrap();
    let mut old_projection_error = 0.0_f64;
    for (s, strip) in strips.iter().enumerate() {
        let p = panels[strip.panel_indices[0]];
        let point = core::array::from_fn(|i| 0.5 * (p.a[i] + p.b[i]));
        let induced = trailing(point, &panels, &result.gamma, free);
        let u = free[0];
        let v = free[2];
        let derivative = (u * induced[2] - v * induced[0]) / (u * u + v * v);
        let base = v.atan2(u);
        close(result.alpha_eff_rad[s], base + derivative, 2e-12);
        let epsilon = 1e-5;
        let upper = (v + epsilon * induced[2]).atan2(u + epsilon * induced[0]);
        let lower = (v - epsilon * induced[2]).atan2(u - epsilon * induced[0]);
        close((upper - lower) / (2.0 * epsilon), derivative, 1e-9);
        old_projection_error =
            old_projection_error.max((derivative - induced[2] / (u * u + v * v).sqrt()).abs());
        let total: f64 = strip.panel_indices.iter().map(|&j| result.gamma[j]).sum();
        let target =
            -0.5 * (u * u + v * v).sqrt() * strip.chord_m * section(s, base + derivative).0;
        assert!((target - total).abs() / (total.abs() + 1e-9) < 2e-10);
    }
    assert!(
        old_projection_error > 1e-4,
        "finite incidence must make chord induction observable"
    );
}

#[test]
fn framed_continuation_reuses_the_physical_equation_and_refuses_legacy_anchors() {
    let (panels, strips) = layout(2);
    let frames = frames(&panels, &strips);
    let free = [12.0, 0.1, 0.7];
    let op = InfluenceOperator::build(&panels, free, RHO).unwrap();
    let local = vec![free; strips.len()];
    let initial = solve_with_section_frames(
        &op,
        &panels,
        &strips,
        free,
        RHO,
        &section,
        &frames,
        &local,
        None,
        controls(),
    )
    .unwrap();
    let next: Vec<[f64; 3]> = local
        .iter()
        .enumerate()
        .map(|(s, v)| [v[0] + 0.5, v[1], v[2] + 0.02 * s as f64])
        .collect();
    let continued = solve_with_section_frames(
        &op,
        &panels,
        &strips,
        free,
        RHO,
        &section,
        &frames,
        &next,
        Some(&initial.continuation),
        controls(),
    )
    .unwrap();
    let direct = solve_nonlinear_with_section_frames(
        &op,
        &panels,
        &strips,
        free,
        RHO,
        &section,
        &frames,
        &next,
        Some(initial.continuation.gamma()),
    )
    .unwrap();
    for (&a, &b) in continued.solution.gamma.iter().zip(&direct.gamma) {
        close(a, b, 3e-9);
    }
    vector_close(continued.solution.total_force_n, direct.total_force_n, 1e-7);
    assert!(continued.continuation.maximum_incidence_change_rad() > 0.0);
    let repeat = solve_with_section_frames(
        &op,
        &panels,
        &strips,
        free,
        RHO,
        &section,
        &frames,
        &next,
        Some(&initial.continuation),
        controls(),
    )
    .unwrap();
    assert_eq!(continued, repeat);
    let legacy = solve_with_incident(
        &op,
        &panels,
        &strips,
        free,
        RHO,
        &section,
        &local,
        None,
        None,
        controls(),
    )
    .unwrap();
    let saved = legacy.continuation.clone();
    assert_eq!(
        solve_with_section_frames(
            &op,
            &panels,
            &strips,
            free,
            RHO,
            &section,
            &frames,
            &local,
            Some(&legacy.continuation),
            controls(),
        )
        .unwrap_err()
        .code,
        "nonlinear-continuation-stale"
    );
    assert_eq!(legacy.continuation, saved);
    assert_eq!(
        solve_with_incident(
            &op,
            &panels,
            &strips,
            free,
            RHO,
            &section,
            &local,
            None,
            Some(&initial.continuation),
            controls(),
        )
        .unwrap_err()
        .code,
        "nonlinear-continuation-stale"
    );
}

#[test]
fn section_frames_admit_actual_geometry_and_reject_wrong_rows_or_reverse_flow() {
    let (mut panels, strips) = layout(2);
    let frames = frames(&panels, &strips);
    assert_eq!(frames[0], SectionFrame::BODY_XZ);
    let flow = frames[0].flow([10.0, 35.0, 2.0]).unwrap();
    close(flow.speed_mps, 104.0_f64.sqrt(), 1e-14);
    assert_eq!(flow.velocity_mps, [10.0, 0.0, 2.0]);
    assert_eq!(flow.spanwise_mps, 35.0);
    assert!(SectionFrame::new([1.0, 0.0, 0.0], [1.0, 0.0, 0.0]).is_err());
    assert!(SectionFrame::new([f64::NAN, 0.0, 0.0], [0.0, 1.0, 0.0]).is_err());
    assert!(frames[0].flow([-1.0, 0.0, 0.0]).is_err());
    assert!(frames[0].flow([0.0, 0.0, 1.0]).is_err());
    assert!(frames[0].flow([1.0, f64::NAN, 0.0]).is_err());
    let free = [12.0, 0.0, 0.7];
    let local = vec![free; strips.len()];
    // A panel normal is a BC orientation, so either sign defines the same plane.
    for p in &mut panels {
        p.normal = p.normal.map(|x| -x);
    }
    assert_eq!(SectionFrame::from_panel(&panels[0]).unwrap(), frames[0]);
    let op = InfluenceOperator::build(&panels, free, RHO).unwrap();
    assert!(
        solve_nonlinear_with_section_frames(
            &op, &panels, &strips, free, RHO, &section, &frames, &local, None
        )
        .is_ok()
    );
    assert_eq!(
        solve_nonlinear_with_section_frames(
            &op,
            &panels,
            &strips,
            free,
            RHO,
            &section,
            &frames[1..],
            &local,
            None
        )
        .unwrap_err()
        .code,
        "section-frame-invalid"
    );
    assert_eq!(
        solve_nonlinear_with_section_frames(
            &op,
            &panels,
            &strips,
            free,
            f64::NAN,
            &section,
            &frames,
            &local,
            None
        )
        .unwrap_err()
        .code,
        "strips-invalid"
    );
    let row = strips[0].panel_indices[1];
    let p = &mut panels[row];
    core::mem::swap(&mut p.a, &mut p.b);
    let op = InfluenceOperator::build(&panels, free, RHO).unwrap();
    assert_eq!(
        solve_nonlinear_with_section_frames(
            &op, &panels, &strips, free, RHO, &section, &frames, &local, None
        )
        .unwrap_err()
        .code,
        "section-frame-invalid"
    );
}
