//! Actual annular balance, local stalled-section branches and bounded work.

use fs_airscrew::{
    BladeStation, MAX_BRAKING_INDUCTION, MAX_STATION_ITERS, Rotor, W_TOL, bemt_solve,
};
use fs_math::det;

const RHO: f64 = 1.294;
const FAILURE_OMEGA: f64 = 45.154808569972736;
const FAILURE_AXIAL: f64 = 3.1309275058577004;

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
        .map(|(r, chord, degrees)| BladeStation {
            r_over_r: r,
            chord_m: chord,
            beta_rad: degrees * std::f64::consts::PI / 180.0,
        })
        .collect(),
    }
}

/// Independent dimensional blade-element loading and annular momentum at
/// the PUBLIC receipt's induction. No production iteration helpers are used.
fn loads(rotor: &Rotor, i: usize, axial: f64, omega: f64, w: f64) -> [f64; 4] {
    let station = rotor.stations[i];
    let radius = station.r_over_r * rotor.radius_m;
    let tangent = omega * radius;
    let phi = det::atan2(axial + w, tangent);
    let alpha = (station.beta_rad - phi).clamp(-1.3, 1.3);
    let attached_cl = 2.0 * std::f64::consts::PI * (alpha + 2.0 * rotor.camber_ratio);
    let attached_cd = 0.012 + 0.8 * alpha * alpha;
    let (cl, cd) = if alpha.abs() <= 0.25 {
        (attached_cl, attached_cd)
    } else {
        let plate = fs_airfoil::flat_plate_separated(alpha, 6.0).unwrap();
        let (cl, cd) = fs_airfoil::body_to_wind(plate.cn, plate.ca, alpha);
        let t = ((alpha.abs() - 0.25) / 0.2).min(1.0);
        let blend = t * t * (3.0 - 2.0 * t);
        (
            attached_cl * (1.0 - blend) + cl * blend,
            attached_cd * (1.0 - blend) + cd.abs() * blend,
        )
    };
    let sine = det::sin(phi);
    let cosine = det::cos(phi);
    let blade_factor = 0.5
        * RHO
        * ((axial + w).powi(2) + tangent.powi(2))
        * station.chord_m
        * f64::from(rotor.n_blades);
    let thrust = blade_factor * (cl * cosine - cd * sine);
    let torque = blade_factor * radius * (cl * sine + cd * cosine);
    let factor = |arm: f64| {
        if sine.abs() <= 1e-9 {
            1.0
        } else {
            2.0 / std::f64::consts::PI
                * det::acos(det::exp(
                    -f64::from(rotor.n_blades) * arm / (2.0 * radius * sine.abs()),
                ))
        }
    };
    let prandtl = (factor(rotor.radius_m - radius)
        * factor(radius - 0.8 * rotor.stations[0].r_over_r * rotor.radius_m))
    .clamp(1e-3, 1.0);
    let coefficient = 4.0 * std::f64::consts::PI * radius * RHO * prandtl;
    [thrust, torque, coefficient * (axial + w) * w, coefficient]
}

fn independent_root(rotor: &Rotor, axial: f64, mut low: f64, mut high: f64) -> f64 {
    let residual = |w| {
        let value = loads(rotor, 0, axial, FAILURE_OMEGA, w);
        value[0] - value[2]
    };
    assert!(residual(low) * residual(high) < 0.0);
    let low_sign = residual(low).signum();
    for _ in 0..60 {
        let midpoint = low + 0.5 * (high - low);
        if residual(midpoint).signum() == low_sign {
            low = midpoint;
        } else {
            high = midpoint;
        }
    }
    low + 0.5 * (high - low)
}

#[test]
fn formerly_refused_native_flight_point_satisfies_actual_annular_momentum() {
    let rotor = rotor();
    let solution = bemt_solve(&rotor, RHO, FAILURE_AXIAL, FAILURE_OMEGA).unwrap();
    let reference = independent_root(&rotor, FAILURE_AXIAL, 4.1, 4.3);
    assert!((solution.stations[0].w_mps - reference).abs() < 2e-9);
    assert!(solution.stations[0].iterations > 96);
    for receipt in &solution.stations {
        assert!(receipt.iterations <= MAX_STATION_ITERS);
    }
}

#[test]
fn local_finish_retains_the_picard_branch_when_three_positive_roots_exist() {
    let rotor = rotor();
    let first = independent_root(&rotor, 3.0, 1.8, 2.1);
    let middle = independent_root(&rotor, 3.0, 2.3, 2.7);
    let last = independent_root(&rotor, 3.0, 4.0, 4.2);
    assert!(first < middle && middle < last);
    let solution = bemt_solve(&rotor, RHO, 3.0, FAILURE_OMEGA).unwrap();
    assert!((solution.stations[0].w_mps - first).abs() < 2e-9);
    assert!(solution.stations[0].iterations > 96);
    assert!(solution.stations[0].iterations <= MAX_STATION_ITERS);
}

#[test]
fn accepted_induction_closes_momentum_and_reconstructs_returned_loads() {
    let rotor = rotor();
    let omega_350 = 350.0 / 60.0 * std::f64::consts::TAU;
    for (axial, omega) in [
        (0.0, FAILURE_OMEGA),
        (3.0, FAILURE_OMEGA),
        (FAILURE_AXIAL, FAILURE_OMEGA),
        (5.0, FAILURE_OMEGA),
        (10.7, FAILURE_OMEGA),
        (13.6, FAILURE_OMEGA),
        (14.3, omega_350),
        (16.0, omega_350),
        (50.0, omega_350),
    ] {
        let solution = bemt_solve(&rotor, RHO, axial, omega).unwrap();
        if axial == 14.3 {
            assert!(solution.stations.iter().any(|station| station.w_mps < 0.0));
            assert!(solution.stations.iter().any(|station| station.w_mps > 0.0));
            assert!(
                solution.thrust_n > 0.0,
                "a powered rotor can have braking inner annuli"
            );
        }
        let mut thrust = 0.0;
        let mut torque = 0.0;
        for (i, receipt) in solution.stations.iter().enumerate() {
            let values = loads(&rotor, i, axial, omega, receipt.w_mps);
            let tolerance = values[3]
                * W_TOL
                * (receipt.w_mps.abs() + 1e-9)
                * (axial + 2.0 * receipt.w_mps + W_TOL);
            assert!(
                (values[0] - values[2]).abs() <= tolerance + 2e-12,
                "V={axial}, station={i}, dimensional momentum defect={}",
                values[0] - values[2]
            );
            let low = if i == 0 {
                0.8 * rotor.stations[0].r_over_r
            } else {
                0.5 * (rotor.stations[i - 1].r_over_r + receipt.r_over_r)
            };
            let high = if i + 1 == solution.stations.len() {
                1.0
            } else {
                0.5 * (receipt.r_over_r + rotor.stations[i + 1].r_over_r)
            };
            let width = rotor.radius_m * (high - low);
            thrust += width * values[0];
            torque += width * values[1];
            assert!(receipt.iterations <= MAX_STATION_ITERS);
            assert!(receipt.w_mps >= -MAX_BRAKING_INDUCTION * axial);
            assert!(axial + 2.0 * receipt.w_mps >= 0.2 * axial);
        }
        assert!((solution.thrust_n - thrust).abs() < 2e-12);
        assert!((solution.torque_nm - torque).abs() < 2e-12);
    }
}

#[test]
fn nonfinite_inputs_and_unresolved_nonreal_static_loading_refuse() {
    let rotor = rotor();
    for (rho, axial, omega) in [
        (f64::INFINITY, 1.0, FAILURE_OMEGA),
        (RHO, f64::NAN, FAILURE_OMEGA),
        (RHO, 1.0, f64::INFINITY),
    ] {
        assert_eq!(
            bemt_solve(&rotor, rho, axial, omega).unwrap_err().code,
            "operating-point-invalid"
        );
    }
    let mut invalid = rotor.clone();
    invalid.camber_ratio = f64::NAN;
    assert_eq!(
        bemt_solve(&invalid, RHO, 1.0, FAILURE_OMEGA)
            .unwrap_err()
            .code,
        "rotor-invalid"
    );
    let mut negative_static = rotor.clone();
    negative_static.camber_ratio = -0.1;
    negative_static.stations[0].beta_rad = 0.001;
    let error = bemt_solve(&negative_static, RHO, 0.0, FAILURE_OMEGA).unwrap_err();
    assert_eq!(error.code, "station-did-not-converge");
    assert!(error.message.contains("non-real trial target"));
}

/// Choose chord from an independently imposed dimensional balance at w=-aV.
/// This changes a test rotor, never a production fixture or force coefficient.
fn specified_braking_rotor(a: f64) -> (Rotor, f64, f64, f64) {
    let mut rotor = rotor();
    rotor.camber_ratio = 0.0;
    rotor.stations[0].beta_rad = 0.001;
    rotor.stations[0].chord_m = 1.0;
    let axial = 2.5;
    let omega = 20.0 / (rotor.radius_m * rotor.stations[0].r_over_r);
    let w = -a * axial;
    let unit_chord = loads(&rotor, 0, axial, omega, w);
    rotor.stations[0].chord_m = unit_chord[2] / unit_chord[0];
    assert!(rotor.stations[0].chord_m > 0.0);
    let exact = loads(&rotor, 0, axial, omega, w);
    assert!((exact[0] - exact[2]).abs() < 2e-12);
    (rotor, axial, omega, w)
}

#[test]
fn a_nonreal_trial_target_can_continue_to_a_valid_braking_root() {
    let (rotor, axial, omega, reference) = specified_braking_rotor(0.35);
    let initial = loads(&rotor, 0, axial, omega, 0.1 * axial);
    assert!(
        axial * axial + 4.0 * initial[0] / initial[3] < 0.0,
        "the original cold trial must have a non-real quadratic target"
    );
    let solution = bemt_solve(&rotor, RHO, axial, omega).unwrap();
    let first = solution.stations[0];
    assert!((first.w_mps - reference).abs() < 2e-9);
    assert!(first.iterations <= MAX_STATION_ITERS);
    let actual = loads(&rotor, 0, axial, omega, first.w_mps);
    assert!((actual[0] - actual[2]).abs() < 2e-8);
}

#[test]
fn a_real_normal_root_does_not_admit_empirically_unsupported_high_induction() {
    let (rotor, axial, omega, reference) = specified_braking_rotor(0.41);
    assert!(
        reference > -0.5 * axial,
        "normal far wake is still downstream"
    );
    let exact = loads(&rotor, 0, axial, omega, reference);
    assert!(axial * axial + 4.0 * exact[0] / exact[3] > 0.0);
    let error = bemt_solve(&rotor, RHO, axial, omega).unwrap_err();
    assert_eq!(error.code, "station-wake-domain");
    assert!(error.message.contains("a=0.4"));
}
