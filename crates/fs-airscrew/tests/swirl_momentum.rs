//! Independent blade loads, original SI momentum balances and limiting cases.
use fs_airscrew::swirl::{self, Controls, Induction};
use fs_airscrew::{BladeStation, Rotor};

const RHO: f64 = 1.225;

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:e} != {expected:e}"
    );
}

fn prandtl(r: f64, phi: f64) -> f64 {
    let one = |distance: f64| {
        2.0 / std::f64::consts::PI * (-1.5 * distance / (r * phi.sin().abs())).exp().acos()
    };
    (one(1.0 - r) * one(r - 0.24)).clamp(1e-3, 1.0)
}

/// Construct a synthetic attached-flow rotor whose ORIGINAL native section
/// law has prescribed w,u at every annulus. No production coefficient changes.
/// For t=u/w, force directions require Cl/Cd=(cos(phi)+t sin(phi))/(t cos(phi)-sin(phi)).
/// Solve Cl=2 pi alpha and Cd=.012+.8 alpha², then choose the physical chord
/// from the axial momentum equation. The angular balance follows independently.
fn manufactured(axial: f64, omega: f64, states: [Induction; 4]) -> Rotor {
    let radii = [0.3, 0.5, 0.7, 0.9];
    let stations = radii
        .into_iter()
        .zip(states)
        .map(|(r, state)| {
            let va = axial + state.axial_mps;
            let vt = omega * r - state.tangential_mps;
            let phi = va.atan2(vt);
            let t = state.tangential_mps / state.axial_mps;
            let ratio = (phi.cos() + t * phi.sin()) / (t * phi.cos() - phi.sin());
            let b = std::f64::consts::TAU;
            let alpha = 0.024 * ratio / (b + (b * b - 4.0 * 0.8 * 0.012 * ratio * ratio).sqrt());
            assert!(alpha.abs() < 0.25);
            let cl = b * alpha;
            let cd = 0.012 + 0.8 * alpha * alpha;
            let thrust_per_chord =
                1.5 * RHO * (va * va + vt * vt) * (cl * phi.cos() - cd * phi.sin());
            let thrust =
                4.0 * std::f64::consts::PI * RHO * r * prandtl(r, phi) * va * state.axial_mps;
            BladeStation {
                r_over_r: r,
                chord_m: thrust / thrust_per_chord,
                beta_rad: phi + alpha,
            }
        })
        .collect();
    Rotor {
        radius_m: 1.0,
        n_blades: 3,
        camber_ratio: 0.0,
        stations,
    }
}

fn hover_fixture() -> (Rotor, [Induction; 4]) {
    let states = [0.3, 0.5, 0.7, 0.9].map(|r| Induction {
        axial_mps: 0.12 * 50.0 * r,
        tangential_mps: 0.03 * 50.0 * r,
    });
    (manufactured(0.0, 50.0, states), states)
}

fn width(rotor: &Rotor, index: usize) -> f64 {
    let low = if index == 0 {
        0.8 * rotor.stations[0].r_over_r
    } else {
        0.5 * (rotor.stations[index - 1].r_over_r + rotor.stations[index].r_over_r)
    };
    let high = if index + 1 == rotor.stations.len() {
        1.0
    } else {
        0.5 * (rotor.stations[index].r_over_r + rotor.stations[index + 1].r_over_r)
    };
    rotor.radius_m * (high - low)
}

/// Public-state reconstruction uses only ordinary trigonometry and the
/// declared attached section equations, independently of production helpers.
fn check(rotor: &Rotor, density: f64, axial: f64, omega: f64, solved: &swirl::Solution) {
    let mut thrust = 0.0;
    let mut torque = 0.0;
    let mut axial_power = 0.0;
    let mut swirl_power = 0.0;
    let mut drag_power = 0.0;
    for (i, result) in solved.stations.iter().enumerate() {
        let blade = rotor.stations[i];
        let r = blade.r_over_r * rotor.radius_m;
        let dr = width(rotor, i);
        let w = result.induction.axial_mps;
        let u = result.induction.tangential_mps;
        let va = axial + w;
        let vt = omega * r - u;
        let phi = va.atan2(vt);
        let alpha = blade.beta_rad - phi;
        assert!(alpha.abs() < 0.25, "oracle fixture must remain attached");
        let cl = std::f64::consts::TAU * (alpha + 2.0 * rotor.camber_ratio);
        let cd = 0.012 + 0.8 * alpha * alpha;
        let pressure =
            0.5 * density * (va * va + vt * vt) * blade.chord_m * f64::from(rotor.n_blades) * dr;
        let dt = pressure * (cl * phi.cos() - cd * phi.sin());
        let dq = pressure * r * (cl * phi.sin() + cd * phi.cos());
        let geometric = 4.0 * std::f64::consts::PI * density * r * prandtl(r, phi) * dr;
        let momentum_t = geometric * va * w;
        let momentum_q = geometric * r * va * u;
        near(result.prandtl_f, prandtl(r, phi), 2e-14);
        near(result.thrust_n, dt, 3e-12 * dt.abs().max(1.0));
        near(result.torque_nm, dq, 3e-12 * dq.abs().max(1.0));
        near(
            result.axial_momentum_flux_n,
            momentum_t,
            3e-12 * dt.abs().max(1.0),
        );
        near(
            result.angular_momentum_flux_nm,
            momentum_q,
            3e-12 * dq.abs().max(1.0),
        );
        let roundoff = 32.0 * f64::EPSILON * (dt.abs() + dq.abs() / r);
        near(dt, momentum_t, 3e-10 * dt.abs() + roundoff);
        near(dq, momentum_q, 3e-10 * dq.abs() + roundoff * r);
        near(
            result.axial_residual_n,
            result.thrust_n - result.axial_momentum_flux_n,
            0.0,
        );
        near(
            result.angular_residual_nm,
            result.torque_nm - result.angular_momentum_flux_nm,
            0.0,
        );
        assert!(result.evaluations <= swirl::MAX_SWIRL_EVALUATIONS);
        thrust += dt;
        torque += dq;
        axial_power += dt * w;
        swirl_power += dq * u / r;
        drag_power += pressure * cd * va.hypot(vt);
    }
    near(solved.thrust_n, thrust, 3e-12 * thrust.abs().max(1.0));
    near(solved.torque_nm, torque, 3e-12 * torque.abs().max(1.0));
    let tolerance = 3e-12 * solved.power.shaft_w.abs().max(1.0);
    near(solved.power.axial_induced_w, axial_power, tolerance);
    near(solved.power.swirl_w, swirl_power, tolerance);
    near(solved.power.profile_drag_w, drag_power, tolerance);
    near(
        solved.power.shaft_w,
        axial * thrust + axial_power + swirl_power + drag_power,
        tolerance,
    );
    near(solved.power.residual_w, 0.0, tolerance);
}

#[test]
fn hover_solves_actual_linear_and_angular_momentum_and_changes_blade_forces() {
    let (rotor, expected) = hover_fixture();
    let result = swirl::bemt_solve(&rotor, RHO, 0.0, 50.0, Controls::default(), None).unwrap();
    for (station, state) in result.stations.iter().zip(expected) {
        near(station.induction.axial_mps, state.axial_mps, 2e-8);
        near(station.induction.tangential_mps, state.tangential_mps, 2e-8);
    }
    check(&rotor, RHO, 0.0, 50.0, &result);
    assert!(result.power.swirl_w > 0.0);
    assert_eq!(result.power.thrust_w, 0.0);
    // Holding axial induction fixed while removing only the tangential
    // induction changes the ORIGINAL blade forces: this is a coupled solve,
    // not an angular-momentum receipt added to unchanged no-swirl forces.
    let mut without_swirl = [0.0; 2];
    for (i, state) in result.stations.iter().enumerate() {
        let blade = rotor.stations[i];
        let va = state.induction.axial_mps;
        let vt = 50.0 * blade.r_over_r;
        let phi = va.atan2(vt);
        let alpha = blade.beta_rad - phi;
        let cl = std::f64::consts::TAU * alpha;
        let cd = 0.012 + 0.8 * alpha * alpha;
        let pressure = 1.5 * RHO * (va * va + vt * vt) * blade.chord_m * width(&rotor, i);
        without_swirl[0] += pressure * (cl * phi.cos() - cd * phi.sin());
        without_swirl[1] += pressure * blade.r_over_r * (cl * phi.sin() + cd * phi.cos());
    }
    assert!((without_swirl[0] - result.thrust_n).abs() > 0.005 * result.thrust_n);
    assert!((without_swirl[1] - result.torque_nm).abs() > 0.005 * result.torque_nm);
}

#[test]
fn windmilling_keeps_signed_angular_momentum_and_rejects_high_induction() {
    for fraction in [0.2, 0.39, 0.4, 0.41] {
        let states = [Induction {
            axial_mps: -fraction * 10.0,
            tangential_mps: -0.15,
        }; 4];
        let rotor = manufactured(10.0, 60.0, states);
        let result = swirl::bemt_solve(&rotor, RHO, 10.0, 60.0, Controls::default(), Some(&states));
        if fraction > 0.4 {
            assert_eq!(result.unwrap_err().code, "swirl-wake-domain");
        } else {
            let result = result.unwrap();
            assert!(result.thrust_n < 0.0 && result.torque_nm < 0.0);
            assert!(result.power.shaft_w < 0.0 && result.power.swirl_w > 0.0);
            check(&rotor, RHO, 10.0, 60.0, &result);
            // Cold starts reach the same interior branch. At the exact 0.4
            // boundary the seed tests admission without asking a numerical
            // root slightly outside that boundary to bypass the physical gate.
            if fraction < 0.4 {
                let cold =
                    swirl::bemt_solve(&rotor, RHO, 10.0, 60.0, Controls::default(), None).unwrap();
                near(cold.thrust_n, result.thrust_n, 1e-8 * result.thrust_n.abs());
                near(
                    cold.torque_nm,
                    result.torque_nm,
                    1e-8 * result.torque_nm.abs(),
                );
            }
        }
    }
}

#[test]
fn weak_loading_recovers_the_no_swirl_limit_without_freezing_angular_transfer() {
    let (base, _) = hover_fixture();
    let mut errors = Vec::new();
    let mut tangential_fractions = Vec::new();
    for factor in [0.03, 0.003, 0.0003, 0.00003] {
        let mut rotor = base.clone();
        for station in &mut rotor.stations {
            station.chord_m *= factor;
            station.beta_rad += 0.35;
        }
        let swirl = swirl::bemt_solve(&rotor, RHO, 8.0, 50.0, Controls::default(), None).unwrap();
        let plain = fs_airscrew::bemt_solve(&rotor, RHO, 8.0, 50.0).unwrap();
        errors.push(
            ((swirl.thrust_n - plain.thrust_n) / plain.thrust_n)
                .abs()
                .max(((swirl.torque_nm - plain.torque_nm) / plain.torque_nm).abs()),
        );
        tangential_fractions.push(
            swirl
                .stations
                .iter()
                .map(|s| s.induction.tangential_mps.abs() / (50.0 * s.r_over_r))
                .fold(0.0_f64, f64::max),
        );
        assert!(swirl.power.swirl_w > 0.0);
    }
    assert!(
        tangential_fractions.windows(2).all(|p| p[0] > 5.0 * p[1]),
        "{tangential_fractions:?}"
    );
    assert!(tangential_fractions[3] < 1e-4, "{tangential_fractions:?}");
    assert!(errors.windows(2).all(|p| p[0] > 5.0 * p[1]), "{errors:?}");
    assert!(errors[3] < 1e-4, "{errors:?}");
}

#[test]
fn torque_free_windmilling_has_zero_swirl_without_false_nonconvergence() {
    let states = [Induction {
        axial_mps: -2.0,
        tangential_mps: 0.0,
    }; 4];
    let rotor = manufactured(10.0, 60.0, states);
    let exact =
        swirl::bemt_solve(&rotor, RHO, 10.0, 60.0, Controls::default(), Some(&states)).unwrap();
    let cold = swirl::bemt_solve(&rotor, RHO, 10.0, 60.0, Controls::default(), None).unwrap();
    for result in [&exact, &cold] {
        check(&rotor, RHO, 10.0, 60.0, result);
        assert!(result.thrust_n < 0.0);
        assert!(result.torque_nm.abs() < 1e-11);
        assert!(
            result
                .stations
                .iter()
                .all(|s| s.induction.tangential_mps.abs() < 1e-11)
        );
        near(result.power.swirl_w, 0.0, 1e-20);
        near(result.power.shaft_w, 0.0, 1e-9);
    }
}

#[test]
fn density_and_velocity_scaling_retain_the_actual_momentum_solution() {
    let (rotor, _) = hover_fixture();
    let base = swirl::bemt_solve(&rotor, RHO, 0.0, 50.0, Controls::default(), None).unwrap();
    let dense = swirl::bemt_solve(&rotor, 2.0 * RHO, 0.0, 50.0, Controls::default(), None).unwrap();
    let fast = swirl::bemt_solve(&rotor, RHO, 0.0, 100.0, Controls::default(), None).unwrap();
    for ((a, b), c) in base
        .stations
        .iter()
        .zip(&dense.stations)
        .zip(&fast.stations)
    {
        near(b.induction.axial_mps, a.induction.axial_mps, 2e-9);
        near(b.induction.tangential_mps, a.induction.tangential_mps, 2e-9);
        near(c.induction.axial_mps, 2.0 * a.induction.axial_mps, 2e-8);
        near(
            c.induction.tangential_mps,
            2.0 * a.induction.tangential_mps,
            2e-8,
        );
    }
    near(dense.thrust_n, 2.0 * base.thrust_n, 2e-8 * base.thrust_n);
    near(dense.torque_nm, 2.0 * base.torque_nm, 2e-8 * base.torque_nm);
    near(fast.thrust_n, 4.0 * base.thrust_n, 2e-8 * base.thrust_n);
    near(fast.torque_nm, 4.0 * base.torque_nm, 2e-8 * base.torque_nm);
    near(
        fast.power.shaft_w,
        8.0 * base.power.shaft_w,
        2e-8 * base.power.shaft_w,
    );
    check(&rotor, RHO, 0.0, 100.0, &fast);
}

#[test]
fn bounded_work_bad_seeds_and_unmodeled_flow_refuse_without_partial_results() {
    let (rotor, _) = hover_fixture();
    let limited = Controls {
        maximum_station_evaluations: 3,
        ..Controls::default()
    };
    assert_eq!(
        swirl::bemt_solve(&rotor, RHO, 0.0, 50.0, limited, None)
            .unwrap_err()
            .code,
        "swirl-nonconvergence"
    );
    for controls in [
        Controls {
            maximum_station_evaluations: 257,
            ..Controls::default()
        },
        Controls {
            relative_tolerance: f64::NAN,
            ..Controls::default()
        },
    ] {
        assert_eq!(
            swirl::bemt_solve(&rotor, RHO, 0.0, 50.0, controls, None)
                .unwrap_err()
                .code,
            "swirl-input-invalid"
        );
    }
    for (density, axial, omega) in [
        (f64::NAN, 0.0, 50.0),
        (RHO, -1.0, 50.0),
        (RHO, 0.0, 0.0),
        (RHO, 0.0, f64::INFINITY),
    ] {
        assert!(
            swirl::bemt_solve(&rotor, density, axial, omega, Controls::default(), None).is_err()
        );
    }
    let bad = [Induction {
        axial_mps: 1.0,
        tangential_mps: 1e3,
    }; 4];
    assert_eq!(
        swirl::bemt_solve(&rotor, RHO, 0.0, 50.0, Controls::default(), Some(&bad))
            .unwrap_err()
            .code,
        "swirl-wake-domain"
    );
    assert_eq!(
        swirl::bemt_solve(&rotor, RHO, 0.0, 50.0, Controls::default(), Some(&bad[..3]))
            .unwrap_err()
            .code,
        "swirl-input-invalid"
    );
}
