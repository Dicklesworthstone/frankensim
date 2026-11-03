//! G0/G3: frozen native flight inputs that exposed an undamped alternating
//! strip mode or excessive damping at zero circulation. Physics is unchanged.
use fs_math::det;
use fs_wing::nonlinear::{
    GAMMA_TOL, InfluenceOperator, MAX_ITERATIONS, NonlinearReport, StripRegime, StripSpec,
    solve_nonlinear, solve_nonlinear_with_incident,
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
    let mut worst = 0.0f64;
    for (s, strip) in strips.iter().enumerate() {
        let chordwise = freestream[0] + wash[s];
        let speed = chordwise.hypot(freestream[2]);
        let incidence = det::atan2(freestream[2], chordwise);
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
        let target = -0.5 * speed * strip.chord_m * cl;
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
        "accepted circulation must satisfy the closure in the actual washed local flow"
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

#[test]
fn sign_changing_wing_circulation_does_not_mistake_normalization_for_divergence() {
    // Actual assisted Historical(3) flight at tick 974. The upper inner wing
    // passes through zero circulation while the dimensional mismatch contracts.
    // The old relative-growth guard halved omega five times to .021875 and
    // exhausted 600 evaluations at 9.460544473899003e-7. No aircraft control,
    // section law, residual tolerance, or iteration budget is changed here.
    let (panels, mut strips) = layout();
    for strip in &mut strips[16..] {
        strip.twist_rad = -0.22515074957833842;
    }
    let freestream = [12.522663350147573, 0.0, -1.9512321004310316];
    let incident: Vec<_> = (0..24)
        .map(|s| match s {
            2..=5 => [11.770859874399594, 0.0, -1.4076438498700594],
            0..=7 => [11.581045305148415, 0.0, -1.4076438498700594],
            10..=13 => [14.736956716146942, 0.0, -1.4076438498700594],
            8..=15 => [14.547142146895762, 0.0, -1.4076438498700594],
            16..=19 => [9.933213726399888, 0.0, 3.2927957285101117],
            _ => [11.031768112232239, 0.0, 3.2927957285101117],
        })
        .collect();
    let op = InfluenceOperator::build(&panels, freestream, 1.17).unwrap();
    let run = |warm| {
        solve_nonlinear_with_incident(
            &op, &panels, &strips, freestream, 1.17, &closure, &incident, None, warm,
        )
        .unwrap()
    };
    let check_physics = |result: &NonlinearReport| {
        let mut worst = 0.0f64;
        for (s, strip) in strips.iter().enumerate() {
            let front = panels[strip.panel_indices[0]];
            let point = core::array::from_fn(|i| 0.5 * (front.a[i] + front.b[i]));
            let complete = induced_velocity_free(point, &panels, &result.gamma, freestream)[2];
            let bound: f64 = panels
                .iter()
                .zip(&result.gamma)
                .map(|(&panel, &gamma)| gamma * bound_vertical(point, panel))
                .sum();
            let v = incident[s];
            let speed = det::sqrt(v[0] * v[0] + v[2] * v[2]);
            let alpha = det::atan2(v[2], v[0]) + strip.twist_rad + (complete - bound) / speed;
            let (cl, regime) = closure(s, alpha);
            assert_eq!(regime, StripRegime::Attached);
            assert_eq!(result.regimes[s], regime);
            let target = -0.5 * speed * strip.chord_m * cl;
            let actual: f64 = strip.panel_indices.iter().map(|&j| result.gamma[j]).sum();
            worst = worst.max((target - actual).abs() / (actual.abs() + 1e-9));
        }
        assert!(
            worst < GAMMA_TOL,
            "independently reconstructed final physical mismatch {worst:e}"
        );
        worst
    };
    let solved = run(None);
    let recomputed = check_physics(&solved);
    println!(
        "frozen sign crossing: evaluations={} returned={} independently_recomputed={}",
        solved.iterations, solved.residual, recomputed
    );
    assert!(solved.iterations < MAX_ITERATIONS);
    assert!(solved.residual < GAMMA_TOL);
    assert_eq!(solved, run(None), "same input and branch replay bitwise");
    let continued = run(Some(&solved.gamma));
    check_physics(&continued);
    assert_eq!(continued.branch_id, solved.branch_id);
    for (cold, warm) in solved.gamma.iter().zip(&continued.gamma) {
        assert!((cold - warm).abs() < 1e-9);
    }
    // These opposing surface loads are physical: do not avoid the zero
    // crossing by clipping circulation or replacing it with an unsigned value.
    let total = |s: usize| {
        strips[s]
            .panel_indices
            .iter()
            .map(|&j| solved.gamma[j])
            .sum::<f64>()
    };
    assert!(total(0) > 0.0 && total(8) < 0.0);
    assert!(total(11) < 0.0 && total(11).abs() < 0.01 * total(8).abs());
}
