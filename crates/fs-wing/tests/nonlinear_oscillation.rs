//! G0/G3: frozen native Huffman launch inputs that exposed an undamped
//! alternating strip mode. Geometry, closure, wash and tolerance are unchanged.
use fs_math::det;
use fs_wing::nonlinear::{
    GAMMA_TOL, InfluenceOperator, MAX_ITERATIONS, NonlinearReport, StripRegime, StripSpec,
    solve_nonlinear,
};
use fs_wing::{Panel, SurfaceId, flat_surface, induced_velocity_free};

fn layout() -> (Vec<Panel>, Vec<StripSpec>) {
    // The exact wright_openloop_v1 layout at the captured tick-two canard
    // setting, frozen here without a dependency from this L3 crate to L4.
    let mut panels = flat_surface(SurfaceId::WingLower, 12.29, 1.981, 0.0, 0.0, 8, 2).unwrap();
    panels.extend(flat_surface(SurfaceId::WingUpper, 12.29, 1.981, 0.0, -1.89, 8, 2).unwrap());
    panels.extend(flat_surface(SurfaceId::CanardLower, 3.66, 0.61, 2.9, 1.05, 4, 1).unwrap());
    panels.extend(flat_surface(SurfaceId::CanardUpper, 3.66, 0.61, 2.9, 0.35, 4, 1).unwrap());
    let mut strips = Vec::new();
    for plane in 0..2 {
        for span in 0..8 {
            strips.push(StripSpec {
                panel_indices: vec![plane * 16 + span, plane * 16 + 8 + span],
                chord_m: 1.981,
                twist_rad: 0.0,
            });
        }
    }
    for s in 0..8 {
        strips.push(StripSpec {
            panel_indices: vec![32 + s],
            chord_m: 0.61,
            twist_rad: 0.13025215778133248,
        });
    }
    (panels, strips)
}

fn closure(_: usize, alpha: f64) -> (f64, StripRegime) {
    let attached = core::f64::consts::TAU * (alpha + 0.1);
    let a = alpha.abs();
    if a <= 0.30 {
        (attached, StripRegime::Attached)
    } else if a < 0.45 {
        let t = (a - 0.30) / 0.15;
        let blend = t * t * (3.0 - 2.0 * t);
        let separated = 1.98 * det::sin(alpha) * det::cos(alpha);
        (
            attached * (1.0 - blend) + separated * blend,
            StripRegime::Blended,
        )
    } else {
        (
            1.98 * det::sin(alpha) * det::cos(alpha),
            StripRegime::Separated,
        )
    }
}

fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    core::array::from_fn(|i| a[i] - b[i])
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a.iter().zip(b).map(|(a, b)| a * b).sum()
}

// Independently remove the finite bound segment from the public complete
// horseshoe field. This does not reuse the nonlinear cached downwash matrix,
// iteration state, or reported residual.
fn bound_vertical(point: [f64; 3], panel: Panel) -> f64 {
    let r1 = sub(point, panel.a);
    let r2 = sub(point, panel.b);
    let cross = [
        r1[1] * r2[2] - r1[2] * r2[1],
        r1[2] * r2[0] - r1[0] * r2[2],
        r1[0] * r2[1] - r1[1] * r2[0],
    ];
    let cross_sq = dot(cross, cross);
    let norm1 = det::sqrt(dot(r1, r1));
    let norm2 = det::sqrt(dot(r2, r2));
    if cross_sq < 1e-8 || norm1 < 1e-8 || norm2 < 1e-8 {
        return 0.0;
    }
    let segment = sub(panel.b, panel.a);
    cross[2] * (dot(segment, r1) / norm1 - dot(segment, r2) / norm2)
        / (4.0 * core::f64::consts::PI * cross_sq)
}

fn final_mismatch(
    panels: &[Panel],
    strips: &[StripSpec],
    freestream: [f64; 3],
    wash: &[f64],
    result: &NonlinearReport,
) -> f64 {
    let speed = dot(freestream, freestream).sqrt();
    let incidence = det::atan2(freestream[2], freestream[0]);
    let mut worst = 0.0f64;
    for (s, strip) in strips.iter().enumerate() {
        let front = panels[strip.panel_indices[0]];
        let point = core::array::from_fn(|i| 0.5 * (front.a[i] + front.b[i]));
        let complete = induced_velocity_free(point, panels, &result.gamma, freestream)[2];
        let bound: f64 = panels
            .iter()
            .zip(&result.gamma)
            .map(|(&panel, &gamma)| gamma * bound_vertical(point, panel))
            .sum();
        let alpha = incidence + strip.twist_rad + (complete - bound) / speed;
        let (cl, regime) = closure(s, alpha);
        assert_eq!(regime, result.regimes[s]);
        let target = -0.5 * (speed + wash[s]) * strip.chord_m * cl;
        let actual: f64 = strip.panel_indices.iter().map(|&j| result.gamma[j]).sum();
        worst = worst.max((target - actual).abs() / (actual.abs() + 1e-9));
    }
    worst
}

fn verify_frozen_launch(freestream: [f64; 3], increment: f64) {
    let (panels, strips) = layout();
    let op = InfluenceOperator::build(&panels, freestream, 1.17).unwrap();
    let mut wash = vec![0.0; strips.len()];
    for strip in [2, 3, 4, 5, 10, 11, 12, 13] {
        wash[strip] = increment;
    }
    let run = |warm| {
        solve_nonlinear(
            &op,
            &panels,
            &strips,
            freestream,
            1.17,
            &closure,
            Some(&wash),
            warm,
        )
        .unwrap()
    };
    let solved = run(None);
    let recomputed = final_mismatch(&panels, &strips, freestream, &wash, &solved);
    println!(
        "frozen launch Vx={}: evaluations={} returned={} independently_recomputed={}",
        freestream[0], solved.iterations, solved.residual, recomputed
    );
    assert!(solved.iterations <= MAX_ITERATIONS);
    assert!(solved.residual < GAMMA_TOL);
    assert!(
        recomputed < GAMMA_TOL,
        "accepted circulation must satisfy the unchanged actual closure"
    );
    assert!(solved.regimes.iter().all(|&r| r == StripRegime::Attached));
    assert_eq!(
        solved,
        run(None),
        "same admitted branch and input replay bitwise"
    );
    let continued = run(Some(&solved.gamma));
    assert_eq!(continued.branch_id, solved.branch_id);
    assert!(final_mismatch(&panels, &strips, freestream, &wash, &continued) < GAMMA_TOL);
    for (a, b) in solved.gamma.iter().zip(continued.gamma) {
        assert!(
            (a - b).abs() < 1e-9,
            "continuation must retain the same attached solution"
        );
    }
}

#[test]
fn slowly_growing_alternating_mode_converges_within_original_budget() {
    // Original owner: all 600 evaluations kept omega=.7, late norm ratio
    // 1.060932378987 while direction cosine approached -1; residual 3.220064.
    verify_frozen_launch(
        [4.075381272311152, 0.0, -0.035457888543788126],
        0.6891081914105683,
    );
}

#[test]
fn alternating_residual_plateau_converges_on_the_attached_branch() {
    // Original owner: a two-cycle of relative residuals 6.113/26.739 with
    // fixed norm 26.83568 evaded the consecutive-growth safeguard to 600.
    verify_frozen_launch(
        [3.0759286923386306, 0.0, -0.02640662086636544],
        0.7534899483164144,
    );
}
