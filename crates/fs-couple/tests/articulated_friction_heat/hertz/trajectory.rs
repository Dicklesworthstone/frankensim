//! G1/G3/G4: native arrest, original finite heat and endpoint source feedback.
#[path = "trajectory/march.rs"]
mod march;

use super::*;
use fs_couple::friction_heat::hertz::trajectory::propose_rigid_hertz_trajectory_heat;
use fs_mbd::tangential_friction::hertz::trajectory::{
    RigidHertzTrajectoryError, RigidHertzTrajectoryReport,
};

fn run_trajectory<'m>(
    pair: &Pair,
    owner: &SourceThermalCoulomb,
    thermal: &DepositionSimulation<'m>,
    route: &FrictionHeatRoute,
    duration: f64,
    controls: RigidHertzHeatControls,
    resolve: impl Fn(&DepositionSimulation<'m>, f64) -> Result<RigidHertzHeatSource<f64>, &'static str>,
    gate: &CancelGate,
) -> Result<
    HeatedRigidHertzCandidate<'m, f64, RigidHertzTrajectoryReport>,
    RigidHertzHeatError<&'static str, RigidHertzTrajectoryError>,
> {
    propose_rigid_hertz_trajectory_heat(
        &pair.bodies,
        &pair.properties,
        pair.contact(owner),
        thermal,
        route,
        thermal.time_s(),
        duration,
        controls,
        &[0.0; 3],
        300.0,
        resolve,
        gate,
    )
}

// Eliminate the original three caloric coordinates and the unchanged ordinary
// link 10. Once arrest is reached, relative kinetic energy fixes the heat, while
// source temperature still determines the actual stopping time and travel.
fn arrested_reference(spin: bool, rate: f64, duration: f64) -> (f64, [f64; 3], f64, f64) {
    let q = 0.5 * rate * rate / mobility(spin);
    let d0 = duration * K[0] * AREA / A;
    let d1 = duration * K[1] * AREA / A;
    let d02 = duration * AREA / (A / K[0] + A / K[2]);
    let c0 = CAP[0] + d02 * CAP[2] / (CAP[2] + d02);
    let response = d0 * c0 / (c0 + d0) + d1 * CAP[1] / (CAP[1] + d1);
    let theta = q / response;
    let t0 = d0 * theta / (c0 + d0);
    let effort = effort_scale(spin) * (MU[0] - SLOPE[0] * theta);
    (
        300.0 + theta,
        [
            300.0 + t0,
            300.0 + d1 * theta / (CAP[1] + d1),
            300.0 + d02 * t0 / (CAP[2] + d02),
        ],
        q,
        rate / (mobility(spin) * effort),
    )
}

#[test]
fn finite_sliding_and_spin_arrest_heat_original_cells_without_dwell_loss() {
    for (spin, rate, durations) in [(false, 0.5, [0.3, 1.5]), (true, 1.0, [3.0, 12.0])] {
        let pair = Pair::new(if spin { 0.0 } else { rate }, if spin { rate } else { 0.0 });
        let owner = law(0);
        let original = owner.clone();
        let thermal = thermal([300.0; 3], So3::identity());
        let routes = routes(&thermal);
        let mut stopped: Option<[RigidBodyState; 2]> = None;
        for duration in durations {
            let candidate = run_trajectory(
                &pair,
                &owner,
                &thermal,
                &routes[0],
                duration,
                controls(),
                source,
                &CancelGate::new(),
            )
            .unwrap();
            let (interface, temperatures, heat, stop) = arrested_reference(spin, rate, duration);
            let motion = &candidate.mechanical;
            close(candidate.source, interface);
            close(
                candidate.heat.interfaces[0].interface_temperature_k,
                interface,
            );
            close(motion.relaxation.kinetic_phase_end_s.unwrap(), stop);
            assert!(stop < duration);
            close(motion.relaxation.moving_duration_s, stop);
            close(motion.relaxation.stationary_duration_s, duration - stop);
            close(
                if spin {
                    motion.relaxation.spin_increment_rad
                } else {
                    motion.relaxation.slip_increment_m[0]
                },
                0.5 * rate * stop,
            );
            close(motion.spin_after_rad_s, 0.0);
            close(motion.slip_velocity_after_m_s[0], 0.0);
            close(motion.relaxation.dissipated_energy_j, heat);
            close(candidate.heat.dissipated_j, heat);
            close(energy(&candidate.thermal) - energy(&thermal), heat);
            close(
                kinetic(&candidate.bodies_after) - kinetic(&pair.bodies),
                -heat,
            );
            close(candidate.combined_energy_residual_j, 0.0);
            for i in 0..3 {
                close(
                    candidate.thermal.temperature_k(CellId(i as u64)).unwrap(),
                    temperatures[i],
                );
            }
            assert!(
                temperatures[2] > 300.0,
                "the original other link must still conduct"
            );
            close(candidate.heat.interfaces[0].body_heat_j.iter().sum(), heat);
            assert_eq!(candidate.thermal.time_s(), duration);
            assert_eq!(candidate.friction_source.owner_identity(), owner.identity());
            assert_eq!(
                candidate
                    .friction_source
                    .friction_source()
                    .state_coordinate("T"),
                Some(candidate.source)
            );
            close(
                candidate.friction_source.law().parameters().kinetic_mu,
                MU[0] - SLOPE[0] * (interface - 300.0),
            );
            assert!(candidate.heat.entropy.unwrap().friction_generation_j_k > 0.0);
            if let Some(previous) = stopped {
                for i in 0..2 {
                    for (actual, expected) in [
                        (
                            candidate.bodies_after[i].linear_momentum_world(),
                            previous[i].linear_momentum_world(),
                        ),
                        (
                            candidate.bodies_after[i].angular_momentum_body(),
                            previous[i].angular_momentum_body(),
                        ),
                    ] {
                        close(actual.x, expected.x);
                        close(actual.y, expected.y);
                        close(actual.z, expected.z);
                    }
                }
            }
            stopped = Some(candidate.bodies_after);
        }
        assert_eq!(owner, original);
        assert_eq!(thermal.time_s(), 0.0);
        assert_eq!(energy(&thermal), 0.0);
    }
}

#[test]
fn stationary_remainder_retains_hot_to_cold_conduction_without_friction_heat() {
    let pair = Pair::new(0.0, 0.0);
    let owner = law(0);
    let initial = [310.0, 290.0, 300.0];
    let thermal = thermal(initial, So3::identity());
    let routes = routes(&thermal);
    let candidate = run_trajectory(
        &pair,
        &owner,
        &thermal,
        &routes[0],
        DT,
        controls(),
        source,
        &CancelGate::new(),
    )
    .unwrap();
    let g = [1, 2].map(|i| DT * AREA / (A / K[0] + A / K[i]));
    let reference = solve_dense(
        vec![
            vec![CAP[0] + g[0] + g[1], -g[0], -g[1]],
            vec![-g[0], CAP[1] + g[0], 0.0],
            vec![-g[1], 0.0, CAP[2] + g[1]],
        ],
        (0..3).map(|i| CAP[i] * initial[i]).collect(),
    );
    for i in 0..3 {
        close(
            candidate.thermal.temperature_k(CellId(i as u64)).unwrap(),
            reference[i],
        );
    }
    assert_eq!(
        candidate.mechanical.relaxation.kinetic_phase_end_s,
        Some(0.0)
    );
    assert_eq!(candidate.mechanical.relaxation.stationary_duration_s, DT);
    assert_eq!(candidate.heat.dissipated_j, 0.0);
    assert_eq!(candidate.bodies_after, pair.bodies);
    assert!(candidate.heat.interfaces[0].body_heat_j[0] < 0.0);
    assert!(candidate.heat.interfaces[0].body_heat_j[1] > 0.0);
    close(candidate.combined_energy_residual_j, 0.0);
    let entropy = candidate.heat.entropy.unwrap();
    assert_eq!(entropy.friction_generation_j_k, 0.0);
    assert!(entropy.interface_conduction_generation_j_k > 0.0);
}

// Independent continuous scalar force/thermal equations, with the massless
// face eliminated algebraically. State is speed, three temperature rises,
// actual travel and physical friction work. The moving branch is continued
// smoothly only to locate its zero; the rest branch has zero friction.
fn derivative(y: [f64; 6], moving: bool) -> [f64; 6] {
    let rate = if moving { y[0] } else { 0.0 };
    let g = [
        K[0] * AREA / A,
        K[1] * AREA / A,
        AREA / (A / K[0] + A / K[2]),
    ];
    let theta =
        (g[0] * y[1] + g[1] * y[2] + 10.0 * MU[0] * rate) / (g[0] + g[1] + 10.0 * SLOPE[0] * rate);
    let effort = if moving {
        10.0 * (MU[0] - SLOPE[0] * theta)
    } else {
        0.0
    };
    [
        -mobility(false) * effort,
        (g[0] * (theta - y[1]) + g[2] * (y[3] - y[1])) / CAP[0],
        g[1] * (theta - y[2]) / CAP[1],
        g[2] * (y[1] - y[3]) / CAP[2],
        rate,
        effort * rate,
    ]
}
fn rk4(y: [f64; 6], dt: f64, moving: bool) -> [f64; 6] {
    let a = derivative(y, moving);
    let b = derivative(std::array::from_fn(|i| y[i] + 0.5 * dt * a[i]), moving);
    let c = derivative(std::array::from_fn(|i| y[i] + 0.5 * dt * b[i]), moving);
    let d = derivative(std::array::from_fn(|i| y[i] + dt * c[i]), moving);
    std::array::from_fn(|i| y[i] + dt * (a[i] + 2.0 * b[i] + 2.0 * c[i] + d[i]) / 6.0)
}
fn continuous_reference(total: f64) -> [f64; 6] {
    let steps = 20_000;
    let dt = total / steps as f64;
    let mut y = [0.5, 0.0, 0.0, 0.0, 0.0, 0.0];
    let mut moving = true;
    for _ in 0..steps {
        let endpoint = rk4(y, dt, moving);
        if moving && endpoint[0] <= 0.0 {
            let (mut lower, mut upper) = (0.0, dt);
            for _ in 0..48 {
                let mid = 0.5 * (lower + upper);
                if rk4(y, mid, true)[0] > 0.0 {
                    lower = mid;
                } else {
                    upper = mid;
                }
            }
            let stop = 0.5 * (lower + upper);
            y = rk4(y, stop, true);
            y[0] = 0.0;
            y = rk4(y, dt - stop, false);
            moving = false;
        } else {
            y = endpoint;
        }
    }
    assert!(!moving);
    close(y[5], 0.5 * 0.5 * 0.5 / mobility(false));
    y
}

#[test]
fn thermal_windows_converge_to_continuous_fade_and_heat_release() {
    let total = 0.3;
    let reference = continuous_reference(total);
    let owner = law(0);
    let mut errors = Vec::new();
    for count in [8, 16, 32] {
        let mut pair = Pair::new(0.5, 0.0);
        let mut thermal = thermal([300.0; 3], So3::identity());
        let route = routes(&thermal)[0].clone();
        let (mut travel, mut heat) = (0.0, 0.0);
        for index in 0..count {
            let target = total * (index + 1) as f64 / count as f64;
            let candidate = run_trajectory(
                &pair,
                &owner,
                &thermal,
                &route,
                target - thermal.time_s(),
                controls(),
                source,
                &CancelGate::new(),
            )
            .unwrap();
            travel += candidate.mechanical.relaxation.slip_increment_m[0];
            heat += candidate.heat.dissipated_j;
            close(candidate.combined_energy_residual_j, 0.0);
            pair.bodies = candidate.bodies_after;
            thermal = candidate.thermal;
        }
        close(heat, reference[5]);
        close(energy(&thermal), reference[5]);
        let mut error = (travel - reference[4]).abs();
        for i in 0..3 {
            error = error.max(
                (thermal.temperature_k(CellId(i as u64)).unwrap() - 300.0 - reference[i + 1]).abs(),
            );
        }
        errors.push(error);
    }
    assert!(
        errors[1] < 0.7 * errors[0] && errors[2] < 0.7 * errors[1],
        "{errors:?}"
    );
    assert!(errors[2] < 2e-4, "{errors:?}");
}

#[test]
fn late_source_cancel_and_native_refusal_publish_no_partial_stop_or_heat() {
    let pair = Pair::new(0.5, 0.0);
    let saved = pair.bodies;
    let owner = law(0);
    let original = owner.clone();
    let thermal = thermal([300.0; 3], So3::identity());
    let saved_cells = thermal.cells().clone();
    let routes = routes(&thermal);
    let reference = run_trajectory(
        &pair,
        &owner,
        &thermal,
        &routes[0],
        0.3,
        controls(),
        source,
        &CancelGate::new(),
    )
    .unwrap();
    let calls = Cell::new(0);
    let failed = run_trajectory(
        &pair,
        &owner,
        &thermal,
        &routes[0],
        0.3,
        controls(),
        |t, temperature| {
            calls.set(calls.get() + 1);
            if calls.get() == reference.trials {
                Err("late trajectory thermal source")
            } else {
                source(t, temperature)
            }
        },
        &CancelGate::new(),
    );
    assert!(matches!(
        failed,
        Err(RigidHertzHeatError::Source(
            "late trajectory thermal source"
        ))
    ));
    assert_eq!(calls.get(), reference.trials);
    let gate = CancelGate::new();
    let calls = Cell::new(0);
    assert!(
        run_trajectory(
            &pair,
            &owner,
            &thermal,
            &routes[0],
            0.3,
            controls(),
            |t, temperature| {
                calls.set(calls.get() + 1);
                if calls.get() == reference.trials {
                    gate.request();
                }
                source(t, temperature)
            },
            &gate,
        )
        .is_err()
    );
    assert!(matches!(
        run_trajectory(
            &pair,
            &owner,
            &thermal,
            &routes[0],
            0.3,
            RigidHertzHeatControls {
                maximum_trials: reference.trials - 1,
                ..controls()
            },
            source,
            &CancelGate::new(),
        ),
        Err(RigidHertzHeatError::TrialBudget)
    ));
    let curved = Pair::new(0.1, 3.0);
    assert!(matches!(
        run_trajectory(
            &curved,
            &owner,
            &thermal,
            &routes[0],
            0.3,
            controls(),
            source,
            &CancelGate::new(),
        ),
        Err(RigidHertzHeatError::Mechanical(_))
    ));
    // Exact tangential relaxation alone exposes normal rates to its caller.
    // This fixed thermal face must refuse initial opening and normal recoil
    // generated by an oblique contact arm throughout an otherwise exact ray.
    let mut opening = pair;
    opening.bodies[0] = RigidBodyState::new(
        saved[0].pose(),
        saved[0]
            .linear_momentum_world()
            .add(Vec3::new(0.0, 0.0, 0.1)),
        saved[0].angular_momentum_body(),
    )
    .unwrap();
    let mut oblique = pair;
    oblique.binding.arm_a_body = Vec3::new(0.03, 0.0, ARMS[0]);
    oblique.bodies[0] = RigidBodyState::new(
        Pose::new(
            saved[0]
                .pose()
                .position_world()
                .add(Vec3::new(-0.03, 0.0, 0.0)),
            saved[0].pose().orientation(),
        )
        .unwrap(),
        saved[0].linear_momentum_world(),
        saved[0].angular_momentum_body(),
    )
    .unwrap();
    for unsupported in [opening, oblique] {
        assert!(matches!(
            run_trajectory(
                &unsupported,
                &owner,
                &thermal,
                &routes[0],
                0.3,
                controls(),
                source,
                &CancelGate::new(),
            ),
            Err(RigidHertzHeatError::Mechanical(_))
        ));
    }
    assert_eq!(pair.bodies, saved);
    assert_eq!(owner, original);
    assert_eq!(thermal.cells(), &saved_cells);
    assert_eq!(thermal.time_s(), 0.0);
    let retry = run_trajectory(
        &pair,
        &owner,
        &thermal,
        &routes[0],
        0.3,
        controls(),
        source,
        &CancelGate::new(),
    )
    .unwrap();
    assert_eq!(retry.bodies_after, reference.bodies_after);
    assert_eq!(retry.thermal.cells(), reference.thermal.cells());
    assert_eq!(retry.mechanical, reference.mechanical);
}
