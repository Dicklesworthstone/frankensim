//! Physical rate-law arrest, finite-cell cooling, reversal and thermal release.
use super::*;

#[derive(Clone, Copy, Debug)]
struct Motion {
    rate: f64,
    angle: f64,
    heat: f64,
    stop: Option<f64>,
}

// Independently integrate each signed exponential branch. Heat uses A*|omega|
// and B*omega² separately; it is never recovered from a native work residual.
fn motion(theta: f64, initial: f64, moment: f64, duration: f64) -> Motion {
    let a = effort_scale(true) * (MU[0] - SLOPE[0] * theta);
    let b = 10.0 * 0.4 * radius().powi(2) * (VISC + VISC_T * theta);
    let decay = mobility(true) * b;
    let segment = |speed: f64, sign: f64, t: f64| {
        let terminal = (sign * moment - a) / b;
        let delta = speed - terminal;
        let first = -(-decay * t).exp_m1() / decay;
        let second = -(-2.0 * decay * t).exp_m1() / (2.0 * decay);
        let travel = terminal * t + delta * first;
        let square = terminal.powi(2) * t + 2.0 * terminal * delta * first + delta.powi(2) * second;
        Motion {
            rate: sign * (terminal + delta * (-decay * t).exp()),
            angle: sign * travel,
            heat: a * travel + b * square,
            stop: None,
        }
    };
    let sign = initial.signum();
    let terminal = (sign * moment - a) / b;
    let stop = if terminal < 0.0 {
        (initial.abs() / -terminal).ln_1p() / decay
    } else {
        f64::INFINITY
    };
    let mut out = segment(initial.abs(), sign, duration.min(stop));
    if stop <= duration {
        out.rate = 0.0;
        out.stop = Some(stop);
        let static_capacity = effort_scale(true) * (0.8 - 0.005 * theta);
        if moment.abs() > static_capacity && duration > stop {
            let reversed = segment(0.0, moment.signum(), duration - stop);
            out.rate = reversed.rate;
            out.angle += reversed.angle;
            out.heat += reversed.heat;
        }
    }
    out
}

#[test]
fn sourced_affine_stop_hold_and_reversal_close_independent_motion_heat_and_motor_work() {
    let pair = Pair::new(0.0, 3.0);
    let owner = rate_owner(false, 320.0);
    let thermal = thermal([300.0; 3], So3::identity());
    for (moment, duration) in [(0.005, 12.0), (-0.1, 3.0)] {
        let theta = temperature_root(duration, |t| motion(t, 3.0, moment, duration).heat);
        let expected = motion(theta, 3.0, moment, duration);
        let out = phase(
            &pair,
            &owner,
            &thermal,
            moment,
            HertzDrivenMode::Sliding,
            duration,
            source,
            &CancelGate::new(),
        )
        .unwrap();
        let path = &out.mechanical.mechanical.driven;
        close(out.source, 300.0 + theta);
        close(path.spin_rate_after_rad_s, expected.rate);
        close(path.spin_increment_rad, expected.angle);
        close(path.first_arrest_time_s.unwrap(), expected.stop.unwrap());
        let bracket = out.mechanical.stopping_time_bracket_s.unwrap();
        close(bracket[0], expected.stop.unwrap());
        assert_eq!(bracket[0], bracket[1]);
        assert!(out.mechanical.exact_invariant && path.arrest_time_is_exact);
        assert!(!out.mechanical.adaptive_error_estimates_available);
        close(out.heat.dissipated_j, expected.heat);
        close(
            out.mechanical.mechanical.external_work_j,
            moment * expected.angle,
        );
        close(
            out.mechanical.mechanical.kinetic_change_j,
            0.5 * (expected.rate.powi(2) - 9.0) / mobility(true),
        );
        close(out.combined_energy_residual_j, 0.0);
        for (i, t) in thermal_response(duration, theta).0.into_iter().enumerate() {
            close(
                out.thermal.temperature_k(CellId(i as u64)).unwrap(),
                300.0 + t,
            );
        }
        assert_eq!(out.friction_source.owner_identity(), owner.identity());
        if moment > 0.0 {
            assert_eq!(path.mode, HertzDrivenMode::Resting);
            let accepted = Pair {
                bodies: out.bodies_after,
                ..pair
            };
            let cooled = phase(
                &accepted,
                &owner,
                &out.thermal,
                moment,
                path.mode,
                1.0,
                source,
                &CancelGate::new(),
            )
            .unwrap();
            assert_eq!(cooled.heat.dissipated_j, 0.0);
            assert_eq!(
                cooled.mechanical.mechanical.driven.mode,
                HertzDrivenMode::Resting
            );
            close(cooled.mechanical.mechanical.external_work_j, 0.0);
            close(energy(&cooled.thermal), energy(&out.thermal));
            assert_ne!(cooled.thermal.cells(), out.thermal.cells());
        } else {
            assert_eq!(path.mode, HertzDrivenMode::Sliding);
            assert!(path.spin_rate_after_rad_s < 0.0);
            assert!(out.heat.dissipated_j > 0.0);
            assert!(out.mechanical.mechanical.external_work_j < 0.0);
        }
    }
    assert_eq!(thermal.time_s(), 0.0);
    assert_eq!(pair.bodies, Pair::new(0.0, 3.0).bodies);
}

#[test]
fn thermal_time_refinement_cannot_erase_a_reversal_when_endpoint_modes_agree() {
    let pair = Pair::new(0.0, 3.0);
    let owner = rate_owner(false, 320.0);
    let thermal = thermal([300.0; 3], So3::identity());
    let (moment, duration) = (-0.1, 1.66);
    // Independent branch/three-cell elimination places the full-phase event
    // after1.66s and the two-temperature event near1.65985765s. Both endpoints
    // are Sliding, but their physical event histories are different.
    let coarse = phase(
        &pair,
        &owner,
        &thermal,
        moment,
        HertzDrivenMode::Sliding,
        duration,
        source,
        &CancelGate::new(),
    )
    .unwrap();
    let first = phase(
        &pair,
        &owner,
        &thermal,
        moment,
        HertzDrivenMode::Sliding,
        0.5 * duration,
        source,
        &CancelGate::new(),
    )
    .unwrap();
    let advanced = Pair {
        bodies: first.bodies_after,
        ..pair
    };
    let second = phase(
        &advanced,
        &owner,
        &first.thermal,
        moment,
        first.mechanical.mechanical.driven.mode,
        0.5 * duration,
        source,
        &CancelGate::new(),
    )
    .unwrap();
    for candidate in [&coarse, &second] {
        assert_eq!(
            candidate.mechanical.mechanical.driven.mode,
            HertzDrivenMode::Sliding
        );
    }
    assert!(
        coarse
            .mechanical
            .mechanical
            .driven
            .first_arrest_time_s
            .is_none()
    );
    assert!(coarse.mechanical.mechanical.driven.spin_rate_after_rad_s > 0.0);
    assert!(
        second
            .mechanical
            .mechanical
            .driven
            .first_arrest_time_s
            .is_some()
    );
    assert!(second.mechanical.mechanical.driven.spin_rate_after_rad_s < 0.0);
    let mut control = march_controls();
    control.initial_step_s = duration;
    control.minimum_step_s = duration;
    control.maximum_step_s = duration;
    // Isolate the discrete physical-event gate from continuous norm gates.
    control.tolerance = HertzHeatMarchTolerance {
        linear_momentum_ns: 1.0,
        angular_momentum_nms: 1.0,
        cell_energy_j: 1.0,
        scaled_travel_m: 1.0,
        dissipated_j: 1.0,
        relative: 0.0,
    };
    let result = propose_rigid_hertz_rate_driven_heat_march(
        &pair.bodies,
        &pair.properties,
        contact(&pair, &owner),
        motor(moment),
        HertzDrivenMode::Sliding,
        &thermal,
        &routes(&thermal)[0],
        duration,
        driven_controls(),
        control,
        300.0,
        source,
        &CancelGate::new(),
    );
    assert!(
        matches!(result, Err(HertzHeatMarchError::Unresolved { defect, .. })
        if defect.arrest_disagreement && !defect.mode_disagreement && defect.normalized == 2.0)
    );
    assert_eq!(thermal.time_s(), 0.0);
}

#[test]
fn original_stribeck_travel_estimates_enter_the_external_work_comparison() {
    let pair = Pair::new(0.0, 3.0);
    let owner = rate_owner(true, 320.0);
    let thermal = thermal([300.0; 3], So3::identity());
    let (moment, duration) = (0.09, 0.2);
    let coarse = phase(
        &pair,
        &owner,
        &thermal,
        moment,
        HertzDrivenMode::Sliding,
        duration,
        source,
        &CancelGate::new(),
    )
    .unwrap();
    let first = phase(
        &pair,
        &owner,
        &thermal,
        moment,
        HertzDrivenMode::Sliding,
        0.5 * duration,
        source,
        &CancelGate::new(),
    )
    .unwrap();
    let advanced = Pair {
        bodies: first.bodies_after,
        ..pair
    };
    // Use the same accepted source temperature as the existing march's next
    // phase guess, so the comparison concerns these actual numerical receipts.
    let second = propose_rigid_hertz_rate_driven_heat(
        &advanced.bodies,
        &advanced.properties,
        contact(&advanced, &owner),
        motor(moment),
        first.mechanical.mechanical.driven.mode,
        &first.thermal,
        &routes(&thermal)[0],
        first.thermal.time_s(),
        0.5 * duration,
        driven_controls(),
        &[0.0; 3],
        first.source,
        source,
        &CancelGate::new(),
    )
    .unwrap();
    let reports = [&coarse.mechanical, &first.mechanical, &second.mechanical];
    let uncertainty: f64 = reports
        .into_iter()
        .map(|native| {
            assert!(native.adaptive_error_estimates_available && !native.exact_invariant);
            let g = native.mechanical.driven.driving_scaled_wrench_n;
            fs_math::det::hypot(fs_math::det::hypot(g[0], g[1]), g[2])
                * native.estimated_scaled_travel_error_m
        })
        .sum();
    let difference = (coarse.mechanical.mechanical.external_work_j
        - (first.mechanical.mechanical.external_work_j
            + second.mechanical.mechanical.external_work_j))
        .abs();
    assert!(uncertainty > 0.0);
    let expected = difference + uncertainty;
    assert!(expected > difference);
    let mut control = march_controls();
    control.initial_step_s = duration;
    control.minimum_step_s = duration;
    control.maximum_step_s = duration;
    control.tolerance = HertzHeatMarchTolerance {
        linear_momentum_ns: 1.0,
        angular_momentum_nms: 1.0,
        cell_energy_j: 1.0,
        scaled_travel_m: 1.0,
        dissipated_j: 0.5 * expected,
        relative: 0.0,
    };
    let result = propose_rigid_hertz_rate_driven_heat_march(
        &pair.bodies,
        &pair.properties,
        contact(&pair, &owner),
        motor(moment),
        HertzDrivenMode::Sliding,
        &thermal,
        &routes(&thermal)[0],
        duration,
        driven_controls(),
        control,
        300.0,
        source,
        &CancelGate::new(),
    );
    let Err(HertzHeatMarchError::Unresolved { defect, .. }) = result else {
        panic!("unresolved original work estimate was accepted")
    };
    assert!((defect.external_mechanical_work_j - expected).abs() <= 64.0 * f64::EPSILON * expected);
    assert!(defect.external_mechanical_work_j > difference);
    assert!(!defect.mode_disagreement && !defect.arrest_disagreement);
}

fn warming(
    state: &DepositionSimulation<'_>,
    temperature: f64,
) -> Result<RigidHertzHeatSource<f64>, &'static str> {
    let mut out = source(state, temperature)?;
    // Independently supplied1K/s in each original caloric cell. Internal
    // friction/conduction is still owned exclusively by the original network.
    out.transport.sources_w = CAP
        .into_iter()
        .enumerate()
        .map(|(i, capacity)| (CellId(i as u64), capacity))
        .collect();
    Ok(out)
}

#[test]
fn stopped_contact_holds_then_original_heating_releases_the_same_motor() {
    let pair = Pair::new(0.0, 3.0);
    let owner = rate_owner(false, 320.0);
    let thermal = thermal([300.0; 3], So3::identity());
    let moment = -0.08;
    let mut control = march_controls();
    control.initial_step_s = 0.2;
    control.maximum_step_s = 0.2;
    let out = propose_rigid_hertz_rate_driven_heat_march(
        &pair.bodies,
        &pair.properties,
        contact(&pair, &owner),
        motor(moment),
        HertzDrivenMode::Sliding,
        &thermal,
        &routes(&thermal)[0],
        9.0,
        driven_controls(),
        control,
        300.0,
        warming,
        &CancelGate::new(),
    )
    .unwrap();
    assert_eq!(out.mode_after, HertzDrivenMode::Sliding);
    let m = out.march;
    let mut first_stop = None;
    let mut held = 0;
    let mut released = None;
    for phase in &m.phases {
        let HertzHeatMarchMotion::RateDriven(native) = &phase.mechanical else {
            panic!("wrong owner")
        };
        let path = &native.mechanical.driven;
        if first_stop.is_none() {
            first_stop = path.first_arrest_time_s.map(|t| phase.start_s + t);
        }
        if path.mode == HertzDrivenMode::Resting {
            held += 1;
            if path.first_arrest_time_s.is_none() {
                assert_eq!(phase.heat.dissipated_j, 0.0);
                close(native.mechanical.external_work_j, 0.0);
            }
        } else if held > 0 && released.is_none() {
            released = Some(phase.end_s);
            assert!(path.spin_rate_after_rad_s < 0.0);
            assert!(
                effort_scale(true)
                    * phase
                        .friction_source
                        .law()
                        .friction_law()
                        .static_coefficient()
                        .unwrap()
                    < moment.abs()
            );
        }
        assert_eq!(phase.friction_source.owner_identity(), owner.identity());
    }
    assert!(first_stop.unwrap() < 3.0);
    assert!(held > 2 && released.unwrap() > first_stop.unwrap() + 1.0);
    close(m.external_heat_j, 9.0 * CAP.iter().sum::<f64>());
    close(m.external_mechanical_work_j, moment * m.spin_increment_rad);
    close(
        m.thermal_energy_change_j,
        m.dissipated_j + m.external_heat_j,
    );
    for residual in m.energy_residual_j {
        assert!(residual.abs() < control.energy_tolerance_j);
    }
    eprintln!(
        "rate thermal arrest/release: stop={first_stop:?}, release={released:?}, held={held}, phases={}, trials={}",
        m.phases.len(),
        m.thermal_trials
    );
    assert_eq!(thermal.time_s(), 0.0);
}

#[test]
fn late_refusal_after_private_arrest_publishes_no_body_mode_or_caloric_state() {
    let pair = Pair::new(0.0, 0.2);
    let owner = rate_owner(false, 320.0);
    let thermal = thermal([300.0; 3], So3::identity());
    let before = thermal.cells().clone();
    let run = |resolve: &dyn Fn(
        &DepositionSimulation<'_>,
        f64,
    ) -> Result<RigidHertzHeatSource<f64>, &'static str>,
               gate: &CancelGate| {
        let mut control = march_controls();
        control.initial_step_s = 0.02;
        control.maximum_step_s = 0.02;
        propose_rigid_hertz_rate_driven_heat_march(
            &pair.bodies,
            &pair.properties,
            contact(&pair, &owner),
            motor(-0.08),
            HertzDrivenMode::Sliding,
            &thermal,
            &routes(&thermal)[0],
            0.4,
            driven_controls(),
            control,
            300.0,
            resolve,
            gate,
        )
    };
    let reference = run(&source, &CancelGate::new()).unwrap();
    assert_eq!(reference.mode_after, HertzDrivenMode::Resting);
    assert!(reference.march.phases.iter().any(|p| {
        match &p.mechanical {
            HertzHeatMarchMotion::RateDriven(n) => n
                .mechanical
                .driven
                .first_arrest_time_s
                .is_some_and(|t| p.start_s + t < 0.2),
            _ => false,
        }
    }));
    let fail = run(
        &|state, t| {
            if state.time_s() > 0.2 {
                Err("after private arrest")
            } else {
                source(state, t)
            }
        },
        &CancelGate::new(),
    );
    assert!(matches!(
        fail,
        Err(HertzHeatMarchError::RateDriven(
            RigidHertzHeatError::Source("after private arrest")
        ))
    ));
    let gate = CancelGate::new();
    assert!(
        run(
            &|state, t| {
                if state.time_s() > 0.2 {
                    gate.request();
                }
                source(state, t)
            },
            &gate
        )
        .is_err()
    );
    assert_eq!(thermal.cells(), &before);
    assert_eq!(thermal.time_s(), 0.0);
    assert_eq!(pair.bodies, Pair::new(0.0, 0.2).bodies);
    let retry = run(&source, &CancelGate::new()).unwrap();
    assert_eq!(retry.mode_after, reference.mode_after);
    assert_eq!(retry.march.bodies_after, reference.march.bodies_after);
    assert_eq!(retry.march.thermal.cells(), reference.march.thermal.cells());
    assert_eq!(
        retry.march.external_mechanical_work_j,
        reference.march.external_mechanical_work_j
    );
}
