//! G1/G3/G4: real motor work, retained static mode and source-resolved heat.
use super::*;
use fs_couple::friction_heat::hertz::driven::*;
use fs_couple::friction_heat::hertz::march::*;
use fs_mbd::Wrench;
use fs_tribo::sliding_patch::contact::driven::{HertzDrivenControls, HertzDrivenMode};
use fs_tribo::sliding_patch::{HertzSlidingPatch, HertzSlidingPatchParameters};

fn driven_controls(step: f64) -> RigidHertzDrivenHeatControls {
    let mut heat = controls();
    // The 12 s conduction block has entries near 2400. A 200 microkelvin
    // interface stencil avoids subtracting roundoff-dominated temperatures
    // in its directional Jacobian; all true residual/work gates are unchanged.
    heat.difference_step_j = 1e-5;
    RigidHertzDrivenHeatControls {
        heat,
        driven: HertzDrivenControls {
            maximum_step_s: step,
            minimum_step_s: 1e-14,
            maximum_attempts: 100_000,
            maximum_evaluations: 400_001,
            scaled_rate_tolerance_m_s: 1e-10,
            rest_scaled_rate_tolerance_m_s: 1e-10,
        },
    }
}
fn march_controls(step: f64) -> HertzHeatMarchControls {
    HertzHeatMarchControls {
        initial_step_s: step,
        minimum_step_s: 1e-9,
        maximum_step_s: step,
        maximum_attempts: 4_000,
        maximum_accepted_phases: 8_000,
        maximum_thermal_trials: 400_000,
        tolerance: HertzHeatMarchTolerance {
            linear_momentum_ns: 2e-6,
            angular_momentum_nms: 2e-7,
            cell_energy_j: 2e-7,
            scaled_travel_m: 2e-7,
            dissipated_j: 2e-7,
            relative: 0.0,
        },
        energy_tolerance_j: 5e-8,
    }
}
fn spin_motor(moment: f64) -> [Wrench; 2] {
    [1.0, -1.0].map(|sign| Wrench {
        force_world: Vec3::ZERO,
        torque_body: Vec3::new(0.0, 0.0, sign * moment),
    })
}
fn common_drive() -> [Wrench; 2] {
    [
        Wrench {
            force_world: Vec3::new(5.0, 2.5, 0.0),
            torque_body: Vec3::new(0.0, 0.0, 0.075),
        },
        Wrench::ZERO,
    ]
}
fn y_mobility() -> f64 {
    1.0 / MASS[0] + 1.0 / MASS[1] + ARMS[0].powi(2) / I[0][0] + ARMS[1].powi(2) / I[1][0]
}
fn components(v: Vec3) -> [f64; 3] {
    [v.x, v.y, v.z]
}

#[allow(clippy::too_many_arguments)]
fn phase<'m>(
    pair: &Pair,
    owner: &SourceThermalCoulomb,
    thermal: &DepositionSimulation<'m>,
    external: [Wrench; 2],
    mode: HertzDrivenMode,
    dt: f64,
    guess: f64,
) -> HeatedRigidHertzDrivenCandidate<'m, f64> {
    propose_rigid_hertz_driven_heat(
        &pair.bodies,
        &pair.properties,
        pair.contact(owner),
        external,
        mode,
        thermal,
        &routes(thermal)[0],
        thermal.time_s(),
        dt,
        driven_controls(0.001),
        &[0.0; 3],
        guess,
        source,
        &CancelGate::new(),
    )
    .unwrap_or_else(|error| panic!("driven Hertz thermal phase: start={}, dt={dt}, mode={mode:?}, temperature_guess={guess}: {error:?}", thermal.time_s()))
}

#[test]
fn retained_static_mode_transmits_joint_force_spin_and_common_motion_work_without_heat() {
    let pair = Pair::new(0.0, 0.0);
    let owner = law(0);
    let thermal = thermal([300.0; 3], So3::identity());
    let external = common_drive();
    let dt = 0.2;
    let out = propose_rigid_hertz_driven_heat_march(
        &pair.bodies,
        &pair.properties,
        pair.contact(&owner),
        external,
        HertzDrivenMode::Resting,
        &thermal,
        &routes(&thermal)[0],
        dt,
        driven_controls(0.001),
        march_controls(0.02),
        300.0,
        source,
        &CancelGate::new(),
    )
    .unwrap();
    assert_eq!(out.mode_after, HertzDrivenMode::Resting);
    let m = out.march;
    // Solve the three actual no-slip acceleration constraints independently.
    // The spin demand exceeds the kinetic capacity: reclassifying a tiny
    // numerical rate as sliding would lose this valid retained static state.
    let force = Vec3::new(-2.5 / mobility(false), -1.25 / y_mobility(), 0.0);
    let moment = -(0.075 / I[0][2]) / mobility(true);
    assert!(moment.abs() > effort_scale(true) * MU[0]);
    let mut kinetic = 0.0;
    let mut work = 0.0;
    for (i, sign) in [(0, 1.0), (1, -1.0)] {
        let resultant = external[i].force_world.add(force.scale(sign));
        let torque = external[i]
            .torque_body
            .add(Vec3::new(0.0, 0.0, ARMS[i]).cross(force.scale(sign)))
            .add(Vec3::new(0.0, 0.0, sign * moment));
        let expected_p = resultant.scale(dt);
        let expected_l = torque.scale(dt);
        for (actual, expected) in components(m.bodies_after[i].linear_momentum_world())
            .into_iter()
            .zip(components(expected_p))
        {
            close(actual, expected);
        }
        for (actual, expected) in components(m.bodies_after[i].angular_momentum_body())
            .into_iter()
            .zip(components(expected_l))
        {
            close(actual, expected);
        }
        kinetic += expected_p.norm_squared() / (2.0 * MASS[i]);
        work += 0.5 * dt * dt * external[i].force_world.dot(resultant) / MASS[i];
        for j in 0..3 {
            kinetic += components(expected_l)[j].powi(2) / (2.0 * I[i][j]);
            work += 0.5 * dt * dt * components(external[i].torque_body)[j] * components(torque)[j]
                / I[i][j];
        }
    }
    assert!(work > 0.1);
    close(work, kinetic);
    close(m.external_mechanical_work_j, work);
    close(m.kinetic_energy_change_j, kinetic);
    assert_eq!(m.dissipated_j, 0.0);
    assert_eq!(m.slip_increment_m, [0.0; 2]);
    assert_eq!(m.spin_increment_rad, 0.0);
    for receipt in &m.phases {
        let HertzHeatMarchMotion::Driven(native) = &receipt.mechanical else {
            panic!("wrong owner")
        };
        assert_eq!(native.driven.mode, HertzDrivenMode::Resting);
        assert_eq!(native.driven.drive_work_j, 0.0);
        assert!(native.external_work_j > 0.0);
        assert_eq!(receipt.heat.dissipated_j, 0.0);
    }
    assert_eq!(m.thermal.cells(), thermal.cells());
    for residual in m.energy_residual_j {
        close(residual, 0.0);
    }
}

// Elimination of the original three caloric cells, including unchanged link10.
fn thermal_response(dt: f64, theta: f64) -> ([f64; 3], f64) {
    let d0 = dt * K[0] * AREA / A;
    let d1 = dt * K[1] * AREA / A;
    let d02 = dt * AREA / (A / K[0] + A / K[2]);
    let c0 = CAP[0] + d02 * CAP[2] / (CAP[2] + d02);
    let t0 = d0 * theta / (c0 + d0);
    (
        [t0, d1 * theta / (CAP[1] + d1), d02 * t0 / (CAP[2] + d02)],
        (d0 * c0 / (c0 + d0) + d1 * CAP[1] / (CAP[1] + d1)) * theta,
    )
}
fn motor_motion(theta: f64, motor: f64, dt: f64) -> (f64, f64, f64, f64) {
    let effort = effort_scale(true) * (MU[0] - SLOPE[0] * theta);
    let deceleration = mobility(true) * (effort - motor);
    let stop = 3.0 / deceleration;
    let moving = dt.min(stop);
    let angle = 3.0 * moving - 0.5 * deceleration * moving * moving;
    (
        (3.0 - deceleration * moving).max(0.0),
        angle,
        effort * angle,
        stop,
    )
}

#[test]
fn motor_work_and_heat_match_independent_spin_thermal_root_through_exact_arrest() {
    let pair = Pair::new(0.0, 3.0);
    let owner = law(0);
    let thermal = thermal([300.0; 3], So3::identity());
    let motor = 0.005;
    for dt in [0.2, 12.0] {
        let (mut lower, mut upper) = (0.0, 20.0);
        for _ in 0..60 {
            let mid = 0.5 * (lower + upper);
            if thermal_response(dt, mid).1 < motor_motion(mid, motor, dt).2 {
                lower = mid;
            } else {
                upper = mid;
            }
        }
        let theta = 0.5 * (lower + upper);
        let (rate, angle, heat, stop) = motor_motion(theta, motor, dt);
        let out = phase(
            &pair,
            &owner,
            &thermal,
            spin_motor(motor),
            HertzDrivenMode::Sliding,
            dt,
            300.0,
        );
        close(out.source, 300.0 + theta);
        close(out.mechanical.driven.spin_rate_after_rad_s, rate);
        close(out.mechanical.driven.spin_increment_rad, angle);
        close(out.heat.dissipated_j, heat);
        close(out.mechanical.external_work_j, motor * angle);
        close(out.mechanical.kinetic_change_j + heat, motor * angle);
        close(out.combined_energy_residual_j, 0.0);
        for (i, t) in thermal_response(dt, theta).0.into_iter().enumerate() {
            close(
                out.thermal.temperature_k(CellId(i as u64)).unwrap(),
                300.0 + t,
            );
        }
        assert_eq!(out.friction_source.owner_identity(), owner.identity());
        assert_eq!(
            out.friction_source.friction_source().state_coordinate("T"),
            Some(out.source)
        );
        if dt > stop {
            assert_eq!(out.mechanical.driven.mode, HertzDrivenMode::Resting);
            close(out.mechanical.driven.first_arrest_time_s.unwrap(), stop);
            assert!(out.mechanical.driven.arrest_time_is_exact);
            let accepted = Pair {
                bodies: out.bodies_after,
                ..pair
            };
            let rest = phase(
                &accepted,
                &owner,
                &out.thermal,
                spin_motor(motor),
                HertzDrivenMode::Resting,
                1.0,
                out.source,
            );
            assert_eq!(rest.heat.dissipated_j, 0.0);
            close(rest.mechanical.external_work_j, 0.0);
            assert!(
                rest.thermal.temperature_k(CellId(0)).unwrap()
                    != out.thermal.temperature_k(CellId(0)).unwrap()
            );
        } else {
            assert_eq!(out.mechanical.driven.mode, HertzDrivenMode::Sliding);
        }
    }
    // At 7.5 s, the single-temperature scalar root is still sliding, while
    // the two evolving source snapshots have arrested. Deliberately permissive
    // continuous gates isolate the discrete hysteresis check: it must demand
    // refinement even when every momentum/work difference would otherwise fit.
    let mut unresolved = march_controls(7.5);
    unresolved.minimum_step_s = 7.5;
    unresolved.tolerance = HertzHeatMarchTolerance {
        linear_momentum_ns: 1.0,
        angular_momentum_nms: 1.0,
        cell_energy_j: 1.0,
        scaled_travel_m: 1.0,
        dissipated_j: 1.0,
        relative: 0.0,
    };
    let mismatch = propose_rigid_hertz_driven_heat_march(
        &pair.bodies,
        &pair.properties,
        pair.contact(&owner),
        spin_motor(motor),
        HertzDrivenMode::Sliding,
        &thermal,
        &routes(&thermal)[0],
        7.5,
        driven_controls(0.001),
        unresolved,
        300.0,
        source,
        &CancelGate::new(),
    );
    assert!(
        matches!(mismatch, Err(HertzHeatMarchError::Unresolved { defect, .. })
        if defect.mode_disagreement && defect.normalized == 2.0)
    );
    assert_eq!(thermal.time_s(), 0.0);
    assert_eq!(pair.bodies, Pair::new(0.0, 3.0).bodies);
}

#[test]
fn original_temperature_dependent_static_capacity_releases_the_same_motor() {
    let pair = Pair::new(0.0, 0.0);
    let owner = law(0);
    let motor = effort_scale(true) * 0.79;
    for initial in [300.0, 304.0] {
        let thermal = thermal([initial; 3], So3::identity());
        let out = phase(
            &pair,
            &owner,
            &thermal,
            spin_motor(motor),
            HertzDrivenMode::Resting,
            0.01,
            initial,
        );
        if initial == 300.0 {
            assert_eq!(out.mechanical.driven.mode, HertzDrivenMode::Resting);
            assert_eq!(out.heat.dissipated_j, 0.0);
        } else {
            assert_eq!(out.mechanical.driven.mode, HertzDrivenMode::Sliding);
            assert!(out.mechanical.driven.spin_rate_after_rad_s > 0.0);
            assert!(out.heat.dissipated_j > 0.0);
            assert!(out.mechanical.external_work_j > out.heat.dissipated_j);
            assert!(out.friction_source.law().parameters().static_mu < 0.79);
            close(
                out.mechanical.external_work_j,
                motor * out.mechanical.driven.spin_increment_rad,
            );
        }
        close(out.combined_energy_residual_j, 0.0);
        assert_eq!(out.friction_source.owner_identity(), owner.identity());
    }
}

// Independent continuous equations: the original pressure integral supplies
// instantaneous traction, while this RK4 integrates all three driven rates,
// the actual caloric capacities, friction heat and motor work. It uses no
// native driven trajectory, thermal Newton solve or adaptive phase controller.
fn continuous_reference(duration: f64, steps: usize) -> [f64; 11] {
    let unit = HertzSlidingPatch::new(HertzSlidingPatchParameters {
        normal_force_n: 10.0,
        contact_radius_m: radius(),
        friction_coefficient: 1.0,
    })
    .unwrap();
    let derivative = |y: [f64; 11]| {
        let wrench = unit.evaluate([y[0], y[1]], y[2]).unwrap();
        let g0 = K[0] * AREA / A;
        let g1 = K[1] * AREA / A;
        let g02 = AREA / (A / K[0] + A / K[2]);
        let p = wrench.dissipated_power_w;
        let theta = (g0 * y[3] + g1 * y[4] + MU[0] * p) / (g0 + g1 + SLOPE[0] * p);
        let mu = MU[0] - SLOPE[0] * theta;
        [
            mobility(false) * (0.5 + mu * wrench.force_n[0]),
            y_mobility() * (0.25 + mu * wrench.force_n[1]),
            mobility(true) * (0.02 + mu * wrench.spin_moment_nm),
            (g0 * (theta - y[3]) + g02 * (y[5] - y[3])) / CAP[0],
            g1 * (theta - y[4]) / CAP[1],
            g02 * (y[3] - y[5]) / CAP[2],
            y[0],
            y[1],
            y[2],
            mu * p,
            0.5 * y[0] + 0.25 * y[1] + 0.02 * y[2],
        ]
    };
    let mut y = [0.1, 0.0, 3.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0];
    let h = duration / steps as f64;
    for _ in 0..steps {
        let a = derivative(y);
        let b = derivative(std::array::from_fn(|i| y[i] + 0.5 * h * a[i]));
        let c = derivative(std::array::from_fn(|i| y[i] + 0.5 * h * b[i]));
        let d = derivative(std::array::from_fn(|i| y[i] + h * c[i]));
        y = std::array::from_fn(|i| y[i] + h * (a[i] + 2.0 * b[i] + 2.0 * c[i] + d[i]) / 6.0);
    }
    y
}

#[test]
fn mixed_driven_motion_heat_and_motor_work_converge_to_continuous_three_rate_reference() {
    let pair = Pair::new(0.1, 3.0);
    let owner = law(0);
    let thermal = thermal([300.0; 3], So3::identity());
    // Equal and opposite wrenches applied at the contact point make supplied
    // full-body work exactly F dot relative slip + M times relative angle.
    let external = std::array::from_fn(|i| {
        let sign = if i == 0 { 1.0 } else { -1.0 };
        let force = Vec3::new(0.5 * sign, 0.25 * sign, 0.0);
        Wrench {
            force_world: force,
            torque_body: Vec3::new(0.0, 0.0, ARMS[i]).cross(force).add(Vec3::new(
                0.0,
                0.0,
                0.02 * sign,
            )),
        }
    });
    let duration = 0.005;
    let reference = continuous_reference(duration, 8_000);
    let check = continuous_reference(duration, 4_000);
    for i in 0..11 {
        assert!((reference[i] - check[i]).abs() < 1e-10);
    }
    let mut errors = Vec::new();
    for step in [0.0005, 0.000125] {
        let mut control = march_controls(2.0 * step);
        // Explicit native/thermal resolution study. These finite local gates
        // also remain active, independently of the continuous global oracle.
        control.tolerance = HertzHeatMarchTolerance {
            linear_momentum_ns: 1e-4,
            angular_momentum_nms: 1e-5,
            cell_energy_j: 1e-5,
            scaled_travel_m: 1e-6,
            dissipated_j: 1e-6,
            relative: 0.0,
        };
        let out = propose_rigid_hertz_driven_heat_march(
            &pair.bodies,
            &pair.properties,
            pair.contact(&owner),
            external,
            HertzDrivenMode::Sliding,
            &thermal,
            &routes(&thermal)[0],
            duration,
            driven_controls(step),
            control,
            300.0,
            source,
            &CancelGate::new(),
        )
        .unwrap();
        assert_eq!(out.mode_after, HertzDrivenMode::Sliding);
        let m = out.march;
        let HertzHeatMarchMotion::Driven(last) = &m.phases.last().unwrap().mechanical else {
            panic!("wrong native owner")
        };
        let rate_error = (last.driven.velocity_after_m_s[0] - reference[0])
            .abs()
            .max((last.driven.velocity_after_m_s[1] - reference[1]).abs())
            .max(radius() * (last.driven.spin_rate_after_rad_s - reference[2]).abs());
        let temperature_error = (0..3)
            .map(|i| {
                (m.thermal.temperature_k(CellId(i as u64)).unwrap() - 300.0 - reference[3 + i])
                    .abs()
            })
            .fold(0.0_f64, f64::max);
        let heat_error = (m.dissipated_j - reference[9]).abs();
        let work_error = (m.external_mechanical_work_j - reference[10]).abs();
        errors.push([rate_error, temperature_error, heat_error, work_error]);
        assert!(last.driven.velocity_after_m_s[1] > 0.0);
        close(
            m.external_mechanical_work_j,
            0.5 * m.slip_increment_m[0]
                + 0.25 * m.slip_increment_m[1]
                + 0.02 * m.spin_increment_rad,
        );
        close(m.thermal_energy_change_j, m.dissipated_j);
        for residual in m.energy_residual_j {
            assert!(residual.abs() < control.energy_tolerance_j);
        }
        eprintln!(
            "driven Hertz thermal: step={step}, phases={}, trials={}, errors={:?}",
            m.phases.len(),
            m.thermal_trials,
            errors.last().unwrap()
        );
    }
    for i in 0..4 {
        assert!(errors[1][i] < 0.5 * errors[0][i], "{errors:?}");
    }
    assert!(
        errors[1][0] < 1e-5 && errors[1][1] < 2e-4 && errors[1][2] < 1e-7 && errors[1][3] < 1e-7,
        "{errors:?}"
    );
}

#[test]
fn late_driven_source_cancel_and_work_budget_refusals_leave_all_input_owners_unchanged() {
    let pair = Pair::new(0.0, 0.0);
    let owner = law(0);
    let thermal = thermal([300.0; 3], So3::identity());
    let cells = thermal.cells().clone();
    let route = routes(&thermal)[0].clone();
    let input_mode = HertzDrivenMode::Resting;
    let run = |duration,
               limit,
               resolve: &dyn Fn(
        &DepositionSimulation<'_>,
        f64,
    ) -> Result<RigidHertzHeatSource<f64>, &'static str>,
               gate: &CancelGate| {
        propose_rigid_hertz_driven_heat_march(
            &pair.bodies,
            &pair.properties,
            pair.contact(&owner),
            common_drive(),
            input_mode,
            &thermal,
            &route,
            duration,
            driven_controls(0.001),
            limit,
            300.0,
            resolve,
            gate,
        )
    };
    let limit = march_controls(0.02);
    let first = run(0.02, limit, &source, &CancelGate::new()).unwrap();
    let reference = run(0.06, limit, &source, &CancelGate::new()).unwrap();
    let calls = Cell::new(0);
    let failure = run(
        0.06,
        limit,
        &|state, t| {
            calls.set(calls.get() + 1);
            if calls.get() > first.march.thermal_trials {
                Err("late driven source")
            } else {
                source(state, t)
            }
        },
        &CancelGate::new(),
    );
    assert!(matches!(
        failure,
        Err(HertzHeatMarchError::Driven(RigidHertzHeatError::Source(
            "late driven source"
        )))
    ));
    assert_eq!(calls.get(), first.march.thermal_trials + 1);
    let cancel = CancelGate::new();
    calls.set(0);
    assert!(
        run(
            0.06,
            limit,
            &|state, t| {
                calls.set(calls.get() + 1);
                if calls.get() > first.march.thermal_trials {
                    cancel.request();
                }
                source(state, t)
            },
            &cancel
        )
        .is_err()
    );
    for budget in 0..3 {
        let mut bounded = limit;
        match budget {
            0 => bounded.maximum_attempts = 1,
            1 => bounded.maximum_accepted_phases = 2,
            _ => bounded.maximum_thermal_trials = first.march.thermal_trials,
        }
        assert!(matches!(
            run(0.06, bounded, &source, &CancelGate::new()),
            Err(HertzHeatMarchError::Budget(_))
        ));
    }
    assert_eq!(thermal.cells(), &cells);
    assert_eq!(thermal.time_s(), 0.0);
    assert_eq!(pair.bodies, Pair::new(0.0, 0.0).bodies);
    assert_eq!(input_mode, HertzDrivenMode::Resting);
    let retry = run(0.06, limit, &source, &CancelGate::new()).unwrap();
    assert_eq!(retry.mode_after, reference.mode_after);
    assert_eq!(retry.march.bodies_after, reference.march.bodies_after);
    assert_eq!(retry.march.thermal.cells(), reference.march.thermal.cells());
    assert_eq!(
        retry.march.external_mechanical_work_j,
        reference.march.external_mechanical_work_j
    );
}
