//! Frozen native branch ambiguity that formerly ended fixed-control flight.
use fs_flyer::aircraft::wright_openloop_v1;
use fs_flyer::propcoupling::{
    CANDIDATE_B, IncidentInflow, coupled_prop_airframe_step_with_incident, wash_factor,
};
use fs_math::det;
use fs_wing::nonlinear::{
    GAMMA_TOL, InfluenceOperator, StripRegime, solve_nonlinear_with_incident,
};
use fs_wing::{Panel, induced_velocity_free};

// The same Estimated thin-airfoil/flat-plate section fixture as the native
// aircraft. Both branches solve this identical nonlinear physical equation.
fn section(_: usize, alpha: f64) -> (f64, StripRegime) {
    let attached = core::f64::consts::TAU * (alpha + 0.1);
    let magnitude = alpha.abs();
    if magnitude <= 0.30 {
        (attached, StripRegime::Attached)
    } else {
        let separated = fs_airfoil::FLAT_PLATE_CD90 * det::sin(alpha) * det::cos(alpha);
        if magnitude >= 0.45 {
            (separated, StripRegime::Separated)
        } else {
            let t = (magnitude - 0.30) / 0.15;
            let weight = t * t * (3.0 - 2.0 * t);
            (
                attached * (1.0 - weight) + separated * weight,
                StripRegime::Blended,
            )
        }
    }
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a.iter().zip(b).map(|(a, b)| a * b).sum()
}

// Independent finite-bound-segment integral, removed from the public complete
// horseshoe field to reconstruct the actual trailing-sheet incidence.
fn bound_vertical(point: [f64; 3], panel: Panel) -> f64 {
    let r1: [f64; 3] = core::array::from_fn(|i| point[i] - panel.a[i]);
    let r2: [f64; 3] = core::array::from_fn(|i| point[i] - panel.b[i]);
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
    let segment = core::array::from_fn(|i| panel.b[i] - panel.a[i]);
    cross[2] * (dot(segment, r1) / norm1 - dot(segment, r2) / norm2)
        / (4.0 * core::f64::consts::PI * cross_sq)
}

#[test]
fn continued_native_branch_closes_both_wing_and_disk_momentum() {
    // Exact former fixed-control midpoint at tick 553. Cold wing restarts
    // jumped between G~4.361 and G~3.209 m/s across x~3.26, then refused the
    // second real residual growth. The original cap, tolerance and plant stay.
    let mut design = wright_openloop_v1();
    for disk in &mut design.disks {
        disk.omega_rad_s = 46.05132511391593;
    }
    let (panels, strips) = design.layout(0.13949685982663312).unwrap();
    let freestream = [4.645176920823272, 0.0, 2.2564945853080958];
    let rho = 1.294;
    let incident: Vec<_> = (0..24)
        .map(|s| match s {
            0..=7 => [4.9428401392783865, 0.0, 2.0846560898207787],
            8..=15 => [4.005201001144776, 0.0, 2.0846560898207787],
            16..=19 => [5.463750771574836, 0.0, 0.5987585111951548],
            _ => [5.116477016710537, 0.0, 0.5987585111951548],
        })
        .collect();
    let disk_incident = [[4.446734775186529, 0.0, 3.202071409267176]; 2];
    let initial_slip = [3.1820175601410963, 3.182017560141096];
    let op = InfluenceOperator::build(&panels, freestream, rho).unwrap();
    let wash: Vec<_> = strips
        .iter()
        .map(|strip| {
            let front = panels[strip.panel_indices[0]];
            let x = 0.5 * (front.a[0] + front.b[0]);
            let y = 0.5 * (front.a[1] + front.b[1]);
            design.disks.map(|disk| {
                wash_factor(
                    x - disk.center_m[0],
                    (y - disk.center_m[1]).abs(),
                    design.rotor.radius_m,
                )
            })
        })
        .collect();
    let initial_flow: Vec<_> = incident
        .iter()
        .zip(&wash)
        .map(|(&v, factors)| {
            let left = factors[0] * initial_slip[0];
            let right = factors[1] * initial_slip[1];
            [v[0] + (left + right), v[1], v[2]]
        })
        .collect();
    let first_wing = solve_nonlinear_with_incident(
        &op,
        &panels,
        &strips,
        freestream,
        rho,
        &section,
        &initial_flow,
        None,
        None,
    )
    .unwrap();
    let run = || {
        coupled_prop_airframe_step_with_incident(
            &op,
            &panels,
            &strips,
            &section,
            &design.rotor,
            &design.disks,
            freestream,
            rho,
            &CANDIDATE_B,
            Some(initial_slip),
            IncidentInflow {
                strip_velocity_mps: &incident,
                disk_velocity_mps: disk_incident,
            },
        )
        .unwrap()
    };
    let result = run();
    assert_eq!(result, run(), "same admitted branch replays bitwise");
    assert!(result.corrections <= CANDIDATE_B.cap);
    assert!(result.residual < CANDIDATE_B.tol);
    let flow = result.strip_velocity_mps.as_ref().unwrap();
    // At the same accepted outer iterate, a fresh cold solve chooses the
    // distinct low-induction root. Its lower strips 2 and 5 become attached;
    // continuation retains those strips on the original separated branch.
    let restarted = solve_nonlinear_with_incident(
        &op, &panels, &strips, freestream, rho, &section, flow, None, None,
    )
    .unwrap();
    for s in [2, 5] {
        assert_eq!(first_wing.regimes[s], StripRegime::Separated);
        assert_eq!(
            section(s, result.strip_alpha_eff_rad[s]).1,
            StripRegime::Separated
        );
        assert_eq!(restarted.regimes[s], StripRegime::Attached);
    }

    let mut worst_wing = 0.0f64;
    for (s, strip) in strips.iter().enumerate() {
        let front = panels[strip.panel_indices[0]];
        let probe = core::array::from_fn(|i| 0.5 * (front.a[i] + front.b[i]));
        let complete = induced_velocity_free(probe, &panels, &result.gamma, freestream)[2];
        let bound: f64 = panels
            .iter()
            .zip(&result.gamma)
            .map(|(&panel, &gamma)| gamma * bound_vertical(probe, panel))
            .sum();
        let v = flow[s];
        let speed = det::sqrt(v[0] * v[0] + v[2] * v[2]);
        let alpha = det::atan2(v[2], v[0]) + strip.twist_rad + (complete - bound) / speed;
        assert!((alpha - result.strip_alpha_eff_rad[s]).abs() < 1e-13);
        let target = -0.5 * speed * strip.chord_m * section(s, alpha).0;
        let actual: f64 = strip.panel_indices.iter().map(|&j| result.gamma[j]).sum();
        worst_wing = worst_wing.max((actual - target).abs() / (actual.abs() + 1e-9));
    }
    assert!(
        worst_wing < GAMMA_TOL,
        "actual section mismatch {worst_wing:e}"
    );

    let area = core::f64::consts::PI * design.rotor.radius_m * design.rotor.radius_m;
    let mut worst_disk = 0.0f64;
    for k in 0..2 {
        let induced =
            induced_velocity_free(design.disks[k].center_m, &panels, &result.gamma, freestream);
        let axial = disk_incident[k][0] + induced[0];
        let native =
            fs_airscrew::bemt_solve(&design.rotor, rho, axial, design.disks[k].omega_rad_s)
                .unwrap();
        assert_eq!(result.thrust_n[k], native.thrust_n);
        assert_eq!(result.torque_nm[k], native.torque_nm);
        let vi = result.w_slip[k];
        let thrust = 2.0 * rho * area * vi * (axial + vi);
        assert!((thrust - native.thrust_n).abs() < 1e-12 * native.thrust_n.abs());
        // Infer the evaluated candidate from its actual returned strip wash,
        // rather than confusing it with the next predictor `w_slip`.
        let s = wash
            .iter()
            .position(|w| w[k] > 0.0 && w[1 - k] == 0.0)
            .unwrap();
        let evaluated = (flow[s][0] - incident[s][0]) / wash[s][k];
        worst_disk = worst_disk.max((vi - evaluated).abs() / (evaluated.abs() + 1e-6));
    }
    assert!(worst_disk < CANDIDATE_B.tol);
    assert!((worst_disk - result.residual).abs() < 1e-14);
    println!(
        "continued frozen flight: corrections={} wing={worst_wing:e} disk={worst_disk:e} slip={:?}",
        result.corrections, result.w_slip
    );
}
