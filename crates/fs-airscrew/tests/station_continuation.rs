//! G1 annular momentum and G3 local continuation at a captured flight failure.

use fs_airscrew::{BladeStation, MAX_STATION_ITERS, Rotor, W_TOL, bemt_solve};
use fs_math::det;

const RHO: f64 = 1.294;
const AXIAL: f64 = 3.079326558429997;
const OMEGA: f64 = 45.4237104851136;

fn rotor() -> Rotor {
    Rotor {
        radius_m: 1.2954,
        n_blades: 2,
        camber_ratio: 0.04,
        stations: [
            (0.30, 0.13, 40.0),
            (0.45, 0.17, 30.0),
            (0.60, 0.20, 23.0),
            (0.75, 0.21, 18.5),
            (0.88, 0.20, 15.5),
            (0.96, 0.16, 14.0),
        ]
        .into_iter()
        .map(|(r, chord, beta)| BladeStation {
            r_over_r: r,
            chord_m: chord,
            beta_rad: beta * core::f64::consts::PI / 180.0,
        })
        .collect(),
    }
}

struct Loads {
    thrust: f64,
    torque: f64,
    momentum: f64,
    coefficient: f64,
    target_w: f64,
}

// Independently reconstruct dimensional blade loads and annular momentum
// from public geometry/receipts. No production station solve helper is used.
fn loads(rotor: &Rotor, index: usize, axial: f64, w: f64) -> Loads {
    let station = rotor.stations[index];
    let radius = rotor.radius_m * station.r_over_r;
    let tangent = OMEGA * radius;
    let phi = det::atan2(axial + w, tangent);
    let alpha = (station.beta_rad - phi).clamp(-1.3, 1.3);
    let attached_cl = 2.0 * core::f64::consts::PI * (alpha + 2.0 * rotor.camber_ratio);
    let attached_cd = 0.012 + 0.8 * alpha * alpha;
    let (cl, cd) = if alpha.abs() <= 0.25 {
        (attached_cl, attached_cd)
    } else {
        let plate = fs_airfoil::flat_plate_separated(alpha, 6.0).unwrap();
        let (separated_cl, separated_cd) = fs_airfoil::body_to_wind(plate.cn, plate.ca, alpha);
        let t = ((alpha.abs() - 0.25) / 0.2).min(1.0);
        let blend = t * t * (3.0 - 2.0 * t);
        (
            attached_cl * (1.0 - blend) + separated_cl * blend,
            attached_cd * (1.0 - blend) + separated_cd.abs() * blend,
        )
    };
    let sine = det::sin(phi);
    let cosine = det::cos(phi);
    let blade = 0.5
        * RHO
        * ((axial + w).powi(2) + tangent.powi(2))
        * station.chord_m
        * f64::from(rotor.n_blades);
    let thrust = blade * (cl * cosine - cd * sine);
    let torque = blade * radius * (cl * sine + cd * cosine);
    let prandtl = |arm: f64| {
        2.0 / core::f64::consts::PI
            * det::acos(det::exp(
                -0.5 * f64::from(rotor.n_blades) * arm / (radius * sine.abs()),
            ))
    };
    let factor = (prandtl(rotor.radius_m - radius)
        * prandtl(radius - 0.8 * rotor.stations[0].r_over_r * rotor.radius_m))
    .clamp(1e-3, 1.0);
    let coefficient = 4.0 * core::f64::consts::PI * radius * RHO * factor;
    let target_product = thrust / coefficient;
    Loads {
        thrust,
        torque,
        momentum: coefficient * (axial + w) * w,
        coefficient,
        target_w: 2.0 * target_product / (axial + det::sqrt(axial * axial + 4.0 * target_product)),
    }
}

fn reference_root(rotor: &Rotor, axial: f64) -> f64 {
    let mut lower = 4.15;
    let mut upper = 4.18;
    let balance = |w: f64| {
        let value = loads(rotor, 0, axial, w);
        value.thrust - value.momentum
    };
    assert!(balance(lower) > 0.0 && balance(upper) < 0.0);
    for _ in 0..64 {
        let middle = lower + 0.5 * (upper - lower);
        if balance(middle) > 0.0 {
            lower = middle;
        } else {
            upper = middle;
        }
    }
    lower + 0.5 * (upper - lower)
}

#[test]
fn captured_tick_577_point_closes_actual_momentum_within_the_original_budget() {
    let rotor = rotor();
    // This unchanged operating point previously exhausted 120 evaluations:
    // its local branch first qualifies at evaluation 113, after the old
    // one-time evaluation-96 gate. The bracket probe consumes evaluation 114.
    let solution = bemt_solve(&rotor, RHO, AXIAL, OMEGA).unwrap();
    assert_eq!(MAX_STATION_ITERS, 120);
    assert_eq!(solution.stations[0].iterations, 118);
    assert!((solution.stations[0].w_mps - reference_root(&rotor, AXIAL)).abs() < 2e-9);
    let mut thrust = 0.0;
    let mut torque = 0.0;
    for (index, receipt) in solution.stations.iter().enumerate() {
        assert!(receipt.iterations <= MAX_STATION_ITERS);
        let values = loads(&rotor, index, AXIAL, receipt.w_mps);
        let tolerance = W_TOL * (receipt.w_mps.abs() + 1e-9);
        assert!((values.target_w - receipt.w_mps).abs() <= tolerance + 2e-14);
        assert!(
            (values.thrust - values.momentum).abs()
                <= values.coefficient * tolerance * (AXIAL + 2.0 * receipt.w_mps + tolerance)
                    + 2e-12
        );
        assert!(AXIAL + 2.0 * receipt.w_mps > 0.0);
        let low = if index == 0 {
            0.8 * receipt.r_over_r
        } else {
            0.5 * (rotor.stations[index - 1].r_over_r + receipt.r_over_r)
        };
        let high = if index + 1 == rotor.stations.len() {
            1.0
        } else {
            0.5 * (receipt.r_over_r + rotor.stations[index + 1].r_over_r)
        };
        let width = rotor.radius_m * (high - low);
        thrust += values.thrust * width;
        torque += values.torque * width;
    }
    assert!((thrust - solution.thrust_n).abs() < 2e-11);
    assert!((torque - solution.torque_nm).abs() < 2e-11);
    assert!(solution.power.residual_w.abs() < 1e-10);
}

#[test]
fn nearby_operating_points_stay_on_the_same_locally_resolved_branch() {
    let rotor = rotor();
    let mut previous: Option<(f64, f64, f64, f64)> = None;
    for delta in [-1e-4, -1e-6, 0.0, 1e-6, 1e-4] {
        let axial = AXIAL + delta;
        let solution = bemt_solve(&rotor, RHO, axial, OMEGA).unwrap();
        let receipt = solution.stations[0];
        assert!(receipt.iterations <= MAX_STATION_ITERS);
        assert!((receipt.w_mps - reference_root(&rotor, axial)).abs() < 2e-9);
        if let Some((old_axial, old_w, old_thrust, old_torque)) = previous {
            let step = axial - old_axial;
            assert!(receipt.w_mps > old_w && receipt.w_mps - old_w < 2.0 * step);
            assert!((solution.thrust_n - old_thrust).abs() < 50.0 * step);
            assert!((solution.torque_nm - old_torque).abs() < 10.0 * step);
        }
        previous = Some((axial, receipt.w_mps, solution.thrust_n, solution.torque_nm));
    }
}

#[test]
fn a_late_bracket_still_refuses_when_the_same_budget_cannot_resolve_it() {
    let rotor = rotor();
    let error = bemt_solve(&rotor, RHO, AXIAL - 0.001, OMEGA).unwrap_err();
    assert_eq!(error.code, "station-did-not-converge");
    assert!(
        error
            .message
            .contains("120 station evaluations exhausted in the local bracket")
    );
    // This is a numerical refusal, not an invented domain exclusion: the
    // independent force equation still has a positive normal-flow root.
    assert!(reference_root(&rotor, AXIAL - 0.001) > 0.0);
}
