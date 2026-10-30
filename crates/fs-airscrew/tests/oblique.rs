//! G0/G1/G3: actual blade-force quadrature, physical power and frame parity.
use fs_airscrew::dynamic_inflow::UniformAxialInflow;
use fs_airscrew::oblique::{AzimuthQuadrature, DiskFrame, ObliqueLoads, SpinDirection};
use fs_airscrew::{BladeStation, Rotor};

const AXIS: [f64; 3] = [1.0, 0.0, 0.0];
const RADIAL: [f64; 3] = [0.0, 1.0, 0.0];
const ZERO: [f64; 3] = [0.0; 3];
const RHO: f64 = 1.225;
const OMEGA: f64 = 60.0;

fn rotor() -> Rotor {
    Rotor {
        radius_m: 1.0,
        n_blades: 3,
        camber_ratio: 0.02,
        stations: [0.25, 0.45, 0.65, 0.85]
            .map(|r| BladeStation {
                r_over_r: r,
                chord_m: 0.12 - 0.04 * r,
                beta_rad: 9.0_f64.atan2(OMEGA * r) + 0.06,
            })
            .to_vec(),
    }
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a.into_iter().zip(b).map(|(x, y)| x * y).sum()
}

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn near(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a:.16e} != {b:.16e}; defect {:.3e}",
        (a - b).abs()
    );
}

fn vector_near(a: [f64; 3], b: [f64; 3], tolerance: f64) {
    for j in 0..3 {
        near(a[j], b[j], tolerance);
    }
}

fn evaluate(
    points: usize,
    spin: SpinDirection,
    velocity: [f64; 3],
    rate: [f64; 3],
) -> ObliqueLoads {
    AzimuthQuadrature::new(points)
        .unwrap()
        .evaluate(
            &rotor(),
            RHO,
            DiskFrame::new(AXIS, RADIAL, spin).unwrap(),
            velocity,
            rate,
            2.0,
            OMEGA,
        )
        .unwrap()
}

#[test]
fn axial_limit_matches_native_sections_and_ignores_already_counted_axial_carrier_rate() {
    let axial = UniformAxialInflow::new(1.0, RHO)
        .unwrap()
        .evaluate(&rotor(), 7.0, 2.0, OMEGA)
        .unwrap();
    let loads = evaluate(32, SpinDirection::AlongThrust, [7.0, 0.0, 0.0], ZERO);
    near(loads.thrust_n, axial.thrust_n, 2e-12);
    near(loads.torque_nm, axial.torque_nm, 2e-12);
    vector_near(loads.force_n, [axial.thrust_n, 0.0, 0.0], 2e-12);
    vector_near(loads.hub_moment_nm, ZERO, 2e-13);
    near(loads.power_residual_w, 0.0, 2e-10);
    assert_eq!(
        loads,
        evaluate(
            32,
            SpinDirection::AlongThrust,
            [7.0, 0.0, 0.0],
            [11.0, 0.0, 0.0]
        )
    );
    assert_eq!(
        loads,
        evaluate(32, SpinDirection::AgainstThrust, [7.0, 0.0, 0.0], ZERO)
    );
}

// Independent attached-section integration. It resolves lift/drag from unit
// velocity vectors (not sin/cos(phi)), uses midpoint azimuths and explicit
// radial panel widths, and does not call the production blade-force owner.
fn attached_reference(velocity: [f64; 3], tilt: [f64; 3]) -> ([f64; 3], [f64; 3], f64, f64) {
    let rotor = rotor();
    let mut force = ZERO;
    let mut moment = ZERO;
    let mut torque = 0.0;
    let mut dissipation = 0.0;
    let n = 8192;
    for (station, dr) in rotor.stations.iter().zip([0.15, 0.2, 0.2, 0.25]) {
        let r = station.r_over_r;
        for k in 0..n {
            let angle = std::f64::consts::TAU * (k as f64 + 0.5) / n as f64;
            let radial = [0.0, angle.cos(), angle.sin()];
            let arm = radial.map(|x| r * x);
            let tangent = cross(AXIS, radial);
            let local = cross(tilt, arm);
            let incident = core::array::from_fn(|j| velocity[j] + local[j]);
            let ua = incident[0] + 2.0;
            let ut = OMEGA * r + dot(incident, tangent);
            let speed = ua.hypot(ut);
            let alpha = station.beta_rad - ua.atan2(ut);
            assert!(
                alpha.abs() < 0.25,
                "reference only applies to attached sections"
            );
            let cl = std::f64::consts::TAU * (alpha + 2.0 * rotor.camber_ratio);
            let cd = 0.012 + 0.8 * alpha * alpha;
            let q = 0.5 * RHO * speed * speed * station.chord_m * 3.0 * dr / n as f64;
            let thrust = q * (cl * ut - cd * ua) / speed;
            let drag = q * (cl * ua + cd * ut) / speed;
            let f = core::array::from_fn(|j| thrust * AXIS[j] - drag * tangent[j]);
            let m = cross(arm, f);
            for j in 0..3 {
                force[j] += f[j];
                moment[j] += m[j];
            }
            torque += r * drag;
            dissipation += q * cd * speed;
        }
    }
    (force, moment, torque, dissipation)
}

#[test]
fn oblique_loads_match_independent_section_integration_and_converge_in_azimuth() {
    let velocity = [7.0, 6.0, -1.5];
    let rate = [0.0, 0.4, -0.2];
    let (force, moment, torque, dissipation) = attached_reference(velocity, rate);
    let error = |loads: &ObliqueLoads| {
        (0..3)
            .map(|j| {
                (loads.force_n[j] - force[j]).powi(2)
                    + (loads.hub_moment_nm[j] - if j == 0 { 0.0 } else { moment[j] }).powi(2)
            })
            .sum::<f64>()
            .sqrt()
    };
    let a = evaluate(8, SpinDirection::AlongThrust, velocity, rate);
    let b = evaluate(12, SpinDirection::AlongThrust, velocity, rate);
    let c = evaluate(16, SpinDirection::AlongThrust, velocity, rate);
    let errors = [error(&a), error(&b), error(&c)];
    println!("oblique attached azimuth errors (8, 12, 16 points): {errors:?}");
    assert!(errors[1] < errors[0] * 0.05);
    assert!(errors[2] < errors[1] * 0.05);
    vector_near(c.force_n, force, 3e-10);
    vector_near(c.hub_moment_nm, [0.0, moment[1], moment[2]], 3e-10);
    near(c.torque_nm, torque, 3e-10);
    near(c.profile_dissipation_w, dissipation, 3e-9);
    assert!(c.force_n[1].abs() > 0.1);
    assert!(c.hub_moment_nm[1].abs() > 0.1);
    near(moment[0], -torque, 2e-10);
}

#[test]
fn handedness_reverses_transverse_hub_moment_while_crossflow_force_still_opposes_flow() {
    let a = evaluate(64, SpinDirection::AlongThrust, [7.0, 4.0, 0.0], ZERO);
    let b = evaluate(64, SpinDirection::AgainstThrust, [7.0, 4.0, 0.0], ZERO);
    vector_near(a.force_n, b.force_n, 2e-11);
    vector_near(a.hub_moment_nm, b.hub_moment_nm.map(|x| -x), 2e-11);
    near(a.thrust_n, b.thrust_n, 2e-11);
    near(a.torque_nm, b.torque_nm, 2e-11);
    assert!(a.force_n[1] < -0.1);
    assert!(a.hub_moment_nm[1].abs() > 0.1);
}

fn rotate(a: [f64; 3]) -> [f64; 3] {
    let (s, c) = 0.43_f64.sin_cos();
    let b = [c * a[0] - s * a[1], s * a[0] + c * a[1], a[2]];
    let (s, c) = 0.71_f64.sin_cos();
    [c * b[0] + s * b[2], b[1], -s * b[0] + c * b[2]]
}

#[test]
fn rotated_frame_and_physical_shaft_carrier_power_are_consistent() {
    let velocity = [7.0, 3.0, -1.0];
    let rate = [9.0, 0.4, -0.2];
    let a = evaluate(64, SpinDirection::AgainstThrust, velocity, rate);
    let b = AzimuthQuadrature::new(64)
        .unwrap()
        .evaluate(
            &rotor(),
            RHO,
            DiskFrame::new(rotate(AXIS), rotate(RADIAL), SpinDirection::AgainstThrust).unwrap(),
            rotate(velocity),
            rotate(rate),
            2.0,
            OMEGA,
        )
        .unwrap();
    vector_near(b.force_n, rotate(a.force_n), 4e-11);
    vector_near(b.hub_moment_nm, rotate(a.hub_moment_nm), 4e-11);
    near(b.torque_nm, a.torque_nm, 4e-11);
    // Negative spin is the physical vector -Ω n. The same axial air torque
    // acts on carrier absolute rotation and the relative shaft coordinate.
    let total_hub = core::array::from_fn(|j| a.hub_moment_nm[j] + AXIS[j] * a.torque_nm);
    let relative_spin_along_axis = -OMEGA - rate[0];
    let total_power =
        dot(a.force_n, velocity) + dot(total_hub, rate) + a.torque_nm * relative_spin_along_axis;
    near(total_power, a.kinematic_power_w, 3e-10);
    near(
        total_power + a.induced_work_w + a.profile_dissipation_w,
        0.0,
        3e-10,
    );
    near(a.power_residual_w, 0.0, 3e-10);
    assert!(a.profile_dissipation_w > 0.0);
    assert!(a.maximum_spanwise_speed_mps > 2.0);
}

#[test]
fn bounded_quadrature_and_unmodeled_section_regimes_refuse_without_clipping() {
    for count in [0, 7, 9, 511, 513, usize::MAX] {
        assert_eq!(
            AzimuthQuadrature::new(count).unwrap_err().code,
            "oblique-rotor-invalid"
        );
    }
    assert!(AzimuthQuadrature::new(512).is_ok());
    assert!(DiskFrame::new(AXIS, AXIS, SpinDirection::AlongThrust).is_err());
    let rule = AzimuthQuadrature::new(32).unwrap();
    let frame = DiskFrame::new(AXIS, RADIAL, SpinDirection::AlongThrust).unwrap();
    assert_eq!(
        rule.evaluate(&rotor(), RHO, frame, [7.0, 100.0, 0.0], ZERO, 2.0, OMEGA)
            .unwrap_err()
            .code,
        "oblique-rotor-section-domain"
    );
    let mut high_pitch = rotor();
    for station in &mut high_pitch.stations {
        station.beta_rad = 1.4;
    }
    assert_eq!(
        rule.evaluate(&high_pitch, RHO, frame, ZERO, ZERO, 0.0, OMEGA)
            .unwrap_err()
            .code,
        "oblique-rotor-section-domain"
    );
    assert!(
        rule.evaluate(&rotor(), f64::MAX, frame, [7.0, 1.0, 0.0], ZERO, 2.0, OMEGA)
            .is_err()
    );
    let expected = evaluate(32, SpinDirection::AlongThrust, [7.0, 2.0, 1.0], ZERO);
    assert_eq!(
        expected,
        evaluate(32, SpinDirection::AlongThrust, [7.0, 2.0, 1.0], ZERO)
    );
}
