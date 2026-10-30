//! G1/G3/G4 event/thermal checks reuse the original two-finite-body fixture.
//! Independent references integrate the scalar force/torque ODE and then solve
//! the original two-capacity heat balance; no native trajectory report enters
//! their heat, stopping-time or endpoint reconstruction.

use super::*;
use fs_contact::interface_binding::bind_dry_interface_system_card;
use fs_contact::interface_binding::thermal_friction::ThermalFrictionAxes;
use fs_couple::friction_heat::hertz_rolling::collinear::*;
use fs_matdb::{
    ClaimSet, InterfaceSystemCard, MaterialStateId, QueryPoint, SurfaceSpec, SystemContext,
};
use fs_material::state_point::{
    MaterialPropertySelection, ScalarAdmissibility, ScalarPropertyRequirement,
};
use fs_mbd::tangential_friction::hertz_rolling::collinear::{
    RigidHertzRollingCollinearControls, RigidHertzRollingCollinearHistory,
};
use fs_qty::Dims;

fn original_source(static_mu: f64, kinetic_mu: f64, slope: f64) -> SourceThermalCoulomb {
    let surface = |name: &str| SurfaceSpec {
        material: MaterialStateId {
            chemistry: name.into(),
            phase: "solid".into(),
            process: "synthetic fixture".into(),
            revision: 0,
        },
        texture_frame: format!("texture-{name}"),
    };
    let mut claims = ClaimSet::new();
    claims
        .insert_claim(super::super::curve("static", 0, [static_mu; 2]))
        .unwrap();
    claims
        .insert_claim(super::super::curve(
            "kinetic",
            0,
            [kinetic_mu + 20.0 * slope, kinetic_mu - 20.0 * slope],
        ))
        .unwrap();
    let card = InterfaceSystemCard::assemble(
        surface("middle"),
        surface("outer-0"),
        SystemContext {
            medium: "dry".into(),
            third_body: None,
            environment: "synthetic".into(),
            history: "original".into(),
        },
        claims,
        vec![],
    )
    .unwrap();
    let requirement = |name| {
        ScalarPropertyRequirement::try_new(name, Dims::NONE, ScalarAdmissibility::NonNegative)
            .unwrap()
    };
    SourceThermalCoulomb::resolve(
        bind_dry_interface_system_card(&card, InputAuthority::SyntheticFixture).unwrap(),
        QueryPoint::new().with("T", 300.0).unwrap(),
        requirement("static"),
        requirement("kinetic"),
        MaterialPropertySelection::SingleClaimOnly,
        ThermalFrictionAxes {
            temperature: "T".into(),
            nominal_pressure: None,
        },
    )
    .unwrap()
}

fn event_controls() -> RigidHertzRollingCollinearHeatControls {
    let old = controls();
    RigidHertzRollingCollinearHeatControls {
        mechanical: RigidHertzRollingCollinearControls {
            mechanical: old.mechanical,
            maximum_phases: 8,
            maximum_root_iterations: 128,
            time_tolerance_s: 1e-12,
        },
        maximum_cells: old.maximum_cells,
        maximum_trials: 512,
        maximum_newton_iterations: 30,
        difference_step_j: old.difference_step_j,
        temperature_scale_j_per_k: old.temperature_scale_j_per_k,
        heat: old.heat,
        combined_energy_tolerance_j: old.combined_energy_tolerance_j,
        newton: old.newton,
    }
}

#[derive(Clone, Copy, Debug)]
struct Motion {
    stop: f64,
    surface: f64,
    bulk: f64,
    bulk_before_stop: f64,
    angle: f64,
    force_impulse: f64,
    moment_impulse: f64,
}

// Independent original planar acceleration equations:
// u_dot=a F+b M, omega_dot=b F+g M, M=-c omega.
// After an actual stop, u=0 needs F=b*c*omega/a; this fixture's independent
// static capacity admits that reaction throughout the remaining decay.
fn integrated_motion(f: &Fixture, kinetic_mu: f64, duration: f64) -> Motion {
    let ma = f.properties[0].mass();
    let ia = 0.4 * ma * R * R;
    let ib = 0.13;
    let a = ma.recip() + 1.0 / 7.0 + R * R / ia + D * D / ib;
    let b = -R / ia + D / ib;
    let g = ia.recip() + ib.recip();
    let c = f.tau * f.patch.normal_load_n() * R;
    let force = -kinetic_mu * f.patch.normal_load_n();
    let lambda = c * g;
    let slide = |t: f64| {
        if lambda == 0.0 {
            let u = f.slip + a * force * t;
            let omega = f.roll + b * force * t;
            let distance = f.slip * t + 0.5 * a * force * t * t;
            let angle = f.roll * t + 0.5 * b * force * t * t;
            (u, omega, distance, angle, 0.0)
        } else {
            let equilibrium = b * force / lambda;
            let difference = f.roll - equilibrium;
            let e1 = -(-lambda * t).exp_m1();
            let e2 = -(-2.0 * lambda * t).exp_m1();
            let u = f.slip + (a - b * b / g) * force * t - (b / g) * difference * e1;
            let omega = f.roll - difference * e1;
            let distance = f.slip * t + 0.5 * (a - b * b / g) * force * t * t
                - (b / g) * difference * (t - e1 / lambda);
            let angle = equilibrium * t + difference * e1 / lambda;
            let bulk = c
                * (equilibrium * equilibrium * t
                    + 2.0 * equilibrium * difference * e1 / lambda
                    + difference * difference * e2 / (2.0 * lambda));
            (u, omega, distance, angle, bulk)
        }
    };
    let stop = if f.slip == 0.0 {
        0.0
    } else {
        assert!(
            slide(duration).0 < 0.0,
            "reference window must include physical arrest"
        );
        let (mut lo, mut hi) = (0.0, duration);
        for _ in 0..100 {
            let mid = 0.5 * (lo + hi);
            if slide(mid).0 > 0.0 {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        0.5 * (lo + hi)
    };
    let (_, omega, distance, slide_angle, slide_bulk) = slide(stop);
    let remaining = duration - stop;
    let static_lambda = c * (g - b * b / a);
    let (static_angle, static_bulk) = if static_lambda == 0.0 {
        (omega * remaining, 0.0)
    } else {
        (
            omega * -(-static_lambda * remaining).exp_m1() / static_lambda,
            c * omega * omega * -(-2.0 * static_lambda * remaining).exp_m1()
                / (2.0 * static_lambda),
        )
    };
    Motion {
        stop,
        surface: -force * distance,
        bulk: slide_bulk + static_bulk,
        bulk_before_stop: slide_bulk,
        angle: slide_angle + static_angle,
        force_impulse: force * stop + (b * c / a) * static_angle,
        moment_impulse: -c * (slide_angle + static_angle),
    }
}

#[derive(Clone, Copy, Debug)]
struct ThermalReference {
    motion: Motion,
    face: f64,
    temperatures: [f64; 2],
}
fn event_reference(
    f: &Fixture,
    before: [f64; 2],
    outside: [f64; 2],
    duration: f64,
    mu: f64,
    slope: f64,
) -> ThermalReference {
    let capacity = [
        f.properties[0].mass() * CP[0],
        f.properties[1].mass() * CP[1],
    ];
    let at = |face: f64| {
        let motion = integrated_motion(f, mu - slope * (face - 300.0), duration);
        let temperatures = std::array::from_fn(|i| {
            (capacity[i] * before[i]
                + duration * (SIDE_G[i] * face + outside[i])
                + if i == 0 { motion.bulk } else { 0.0 })
                / (capacity[i] + duration * SIDE_G[i])
        });
        (
            duration
                * (0..2)
                    .map(|i| SIDE_G[i] * (face - temperatures[i]))
                    .sum::<f64>()
                - motion.surface,
            ThermalReference {
                motion,
                face,
                temperatures,
            },
        )
    };
    let (mut lo, mut hi) = (280.0, 320.0);
    assert!(at(lo).0 < 0.0 && at(hi).0 > 0.0);
    for _ in 0..100 {
        let mid = 0.5 * (lo + hi);
        if at(mid).0 > 0.0 {
            hi = mid;
        } else {
            lo = mid;
        }
    }
    at(0.5 * (lo + hi)).1
}

fn verify_event(
    f: &Fixture,
    before: &DepositionSimulation<'_>,
    result: &HeatedRigidHertzRollingCollinearCandidate<'_, [f64; 3]>,
    expected: ThermalReference,
    outside: [f64; 2],
    duration: f64,
) {
    let native = &result.mechanical;
    assert!(native.history.is_stopped());
    let stopped = native
        .report
        .phases
        .iter()
        .find(|phase| phase.reached_stop)
        .unwrap();
    near(
        stopped.start_time_s + stopped.duration_s,
        expected.motion.stop,
        5e-11,
    );
    near(
        native
            .report
            .phases
            .iter()
            .map(|phase| phase.duration_s)
            .sum(),
        duration,
        1e-14,
    );
    assert!(native.report.phases.len() >= 2);
    assert_eq!(native.report.phases.last().unwrap().surface_heat_j, 0.0);
    near(native.report.surface_heat_j, expected.motion.surface, 8e-11);
    near(native.report.bulk_heat_j, expected.motion.bulk, 8e-11);
    near(result.source[0], expected.face, 2e-7);
    for i in 0..2 {
        near(
            result.thermal.temperature_k(IDS[i]).unwrap(),
            expected.temperatures[i],
            2e-7,
        );
    }
    let initial = f.bodies[0].linear_momentum_world().x;
    near(
        native.bodies_after[0].linear_momentum_world().x,
        initial + expected.motion.force_impulse,
        3e-10,
    );
    near(
        native.bodies_after[1].linear_momentum_world().x,
        -expected.motion.force_impulse,
        3e-10,
    );
    near(
        native.bodies_after[0].angular_momentum_body().y,
        f.bodies[0].angular_momentum_body().y - R * expected.motion.force_impulse
            + expected.motion.moment_impulse,
        3e-11,
    );
    near(
        native.bodies_after[1].angular_momentum_body().y,
        -D * expected.motion.force_impulse - expected.motion.moment_impulse,
        3e-11,
    );
    near(
        pair_energy(f, &native.bodies_after) - pair_energy(f, &f.bodies)
            + energy_change(before, &result.thermal)
            - duration * outside.iter().sum::<f64>(),
        0.0,
        3e-10,
    );
    near(
        result.external_heat_j,
        duration * outside.iter().sum::<f64>(),
        1e-14,
    );
    near(
        native.history.rolling_history().state().cumulative_heat_j(),
        expected.motion.bulk,
        8e-11,
    );
    let radius = (f.patch.contact_area_m2() / core::f64::consts::PI).sqrt();
    assert_eq!(
        result.friction_source,
        f.source
            .hertz_at(result.source[0], None, f.patch.normal_load_n(), radius,)
            .unwrap()
    );
}

#[test]
fn exact_coulomb_stop_preserves_work_and_continues_original_conduction() {
    let mut f = Fixture::new(0.0, 0.2, 0.05);
    f.source = original_source(0.8, 0.4, 0.0);
    let initial = [296.0, 304.0];
    let duration = 0.5;
    let t = thermal(&f, initial, false);
    let route = route(&f, &t);
    let history = RigidHertzRollingCollinearHistory::default();
    let result = propose_rigid_hertz_rolling_collinear_heat(
        &f.bodies,
        &f.properties,
        &history,
        f.contact(),
        &t,
        &route,
        bulk_binding(),
        0.0,
        duration,
        event_controls(),
        &[0.0; 2],
        300.0,
        |state, face| source(&f, state, face, [0.0; 2]),
        &CancelGate::new(),
    )
    .unwrap();
    let expected = event_reference(&f, initial, [0.0; 2], duration, 0.4, 0.0);
    verify_event(&f, &t, &result, expected, [0.0; 2], duration);
    assert!(expected.motion.stop < duration / 10.0);
    assert_eq!(result.mechanical.report.bulk_heat_j, 0.0);
    assert!(result.thermal.temperature_k(IDS[0]).unwrap() > initial[0]);
    assert!(result.thermal.temperature_k(IDS[1]).unwrap() < initial[1]);
    assert_eq!(history.rolling_history(), &RigidRollingHistory::new());
}

#[test]
fn arbitrarily_large_static_capacity_does_not_erase_initial_sliding_or_stopping_heat() {
    let mut f = Fixture::new(0.0, 0.2, 0.0);
    let duration = 0.3;
    let mut results = Vec::new();
    for static_mu in [0.8, 20.0] {
        f.source = original_source(static_mu, 0.4, 0.0);
        let t = thermal(&f, [300.0; 2], false);
        let route = route(&f, &t);
        let result = propose_rigid_hertz_rolling_collinear_heat(
            &f.bodies,
            &f.properties,
            &RigidHertzRollingCollinearHistory::default(),
            f.contact(),
            &t,
            &route,
            bulk_binding(),
            0.0,
            duration,
            event_controls(),
            &[0.0; 2],
            300.0,
            |state, face| source(&f, state, face, [0.0; 2]),
            &CancelGate::new(),
        )
        .unwrap();
        let expected = event_reference(&f, [300.0; 2], [0.0; 2], duration, 0.4, 0.0);
        verify_event(&f, &t, &result, expected, [0.0; 2], duration);
        assert!(result.mechanical.report.surface_heat_j > 0.01);
        results.push(result);
    }
    assert_eq!(
        results[0].mechanical.bodies_after,
        results[1].mechanical.bodies_after
    );
    near(
        results[0].mechanical.report.surface_heat_j,
        results[1].mechanical.report.surface_heat_j,
        1e-14,
    );
}

#[test]
fn actual_face_temperature_changes_stopping_and_bulk_temperature_keeps_its_own_source() {
    // The original quasistatic gate compares relaxation to indentation/speed,
    // not contact-radius/speed. Keep its unmodified 0.1 limit throughout arrest.
    let mut f = Fixture::new(0.0005, 0.2, 0.05);
    f.source = original_source(0.8, 0.4, 0.01);
    let initial = [305.0, 301.0];
    let outside = [0.02, -0.01];
    let duration = 5.0;
    let t = thermal(&f, initial, false);
    let route = route(&f, &t);
    let result = propose_rigid_hertz_rolling_collinear_heat(
        &f.bodies,
        &f.properties,
        &RigidHertzRollingCollinearHistory::default(),
        f.contact(),
        &t,
        &route,
        bulk_binding(),
        0.0,
        duration,
        event_controls(),
        &[0.0; 2],
        300.0,
        |state, face| source(&f, state, face, outside),
        &CancelGate::new(),
    )
    .unwrap();
    let expected = event_reference(&f, initial, outside, duration, 0.4, 0.01);
    verify_event(&f, &t, &result, expected, outside, duration);
    assert!(expected.motion.stop > integrated_motion(&f, 0.4, duration).stop);
    near(
        result.mechanical.report.rolling_temperature_k,
        expected.temperatures[0],
        2e-7,
    );
    assert!((expected.face - expected.temperatures[0]).abs() > 0.05);
}

#[test]
fn exact_rolling_variance_and_post_stop_decay_heat_the_original_sphere_once() {
    let mut f = Fixture::new(0.0005, 0.2, 0.05);
    f.source = original_source(0.8, 0.4, 0.0);
    // A longer admitted window retains substantial physical rolling decay
    // without increasing material viscosity outside its original domain.
    let duration = 10.0;
    let initial = [300.0; 2];
    let t = thermal(&f, initial, false);
    let route = route(&f, &t);
    let result = propose_rigid_hertz_rolling_collinear_heat(
        &f.bodies,
        &f.properties,
        &RigidHertzRollingCollinearHistory::default(),
        f.contact(),
        &t,
        &route,
        bulk_binding(),
        0.0,
        duration,
        event_controls(),
        &[0.0; 2],
        300.0,
        |state, face| source(&f, state, face, [0.0; 2]),
        &CancelGate::new(),
    )
    .unwrap();
    let expected = event_reference(&f, initial, [0.0; 2], duration, 0.4, 0.0);
    verify_event(&f, &t, &result, expected, [0.0; 2], duration);
    let c = f.tau * f.patch.normal_load_n() * R;
    let incorrect_mean_heat = c * expected.motion.angle.powi(2) / duration;
    assert!(expected.motion.bulk > incorrect_mean_heat + 1e-5);
    assert!(expected.motion.bulk > 10.0 * expected.motion.bulk_before_stop);
    let caloric = energy_change(&t, &result.thermal);
    near(
        caloric,
        expected.motion.surface + expected.motion.bulk,
        3e-10,
    );
}

#[test]
fn opaque_stopped_history_continues_decay_and_cooling_without_restarting_slip() {
    let mut f = Fixture::new(0.0005, 0.2, 0.05);
    f.source = original_source(0.8, 0.4, 0.0);
    let t = thermal(&f, [300.0; 2], false);
    let route = route(&f, &t);
    let first = propose_rigid_hertz_rolling_collinear_heat(
        &f.bodies,
        &f.properties,
        &RigidHertzRollingCollinearHistory::default(),
        f.contact(),
        &t,
        &route,
        bulk_binding(),
        0.0,
        5.0,
        event_controls(),
        &[0.0; 2],
        300.0,
        |state, face| source(&f, state, face, [0.0; 2]),
        &CancelGate::new(),
    )
    .unwrap();
    let second = propose_rigid_hertz_rolling_collinear_heat(
        &first.mechanical.bodies_after,
        &f.properties,
        &first.mechanical.history,
        f.contact(),
        &first.thermal,
        &route,
        bulk_binding(),
        5.0,
        5.0,
        event_controls(),
        &[0.0; 2],
        first.source[0],
        |state, face| source(&f, state, face, [0.0; 2]),
        &CancelGate::new(),
    )
    .unwrap();
    assert!(second.mechanical.history.is_stopped());
    assert_eq!(second.mechanical.report.surface_heat_j, 0.0);
    assert!(second.mechanical.report.bulk_heat_j > 0.0);
    let full = integrated_motion(&f, 0.4, 10.0);
    near(
        first.mechanical.report.surface_heat_j + second.mechanical.report.surface_heat_j,
        full.surface,
        1e-10,
    );
    near(
        first.mechanical.report.bulk_heat_j + second.mechanical.report.bulk_heat_j,
        full.bulk,
        1e-10,
    );
    near(
        second.mechanical.bodies_after[0].linear_momentum_world().x,
        f.bodies[0].linear_momentum_world().x + full.force_impulse,
        3e-10,
    );
    near(
        pair_energy(&f, &second.mechanical.bodies_after) - pair_energy(&f, &f.bodies)
            + energy_change(&t, &second.thermal),
        0.0,
        5e-10,
    );
}

#[test]
fn discarded_stopping_trajectories_cancellation_and_late_sources_leave_exact_retry() {
    let mut f = Fixture::new(0.0005, 0.2, 0.05);
    f.source = original_source(0.8, 0.4, 0.0);
    let t = thermal(&f, [300.0; 2], false);
    let route = route(&f, &t);
    let history = RigidHertzRollingCollinearHistory::default();
    let run = |controls, gate: &CancelGate| {
        propose_rigid_hertz_rolling_collinear_heat(
            &f.bodies,
            &f.properties,
            &history,
            f.contact(),
            &t,
            &route,
            bulk_binding(),
            0.0,
            5.0,
            controls,
            &[0.0; 2],
            300.0,
            |state, face| source(&f, state, face, [0.0; 2]),
            gate,
        )
    };
    let baseline = run(event_controls(), &CancelGate::new()).unwrap();
    let calls = Cell::new(0usize);
    let failure = propose_rigid_hertz_rolling_collinear_heat(
        &f.bodies,
        &f.properties,
        &history,
        f.contact(),
        &t,
        &route,
        bulk_binding(),
        0.0,
        5.0,
        event_controls(),
        &[0.0; 2],
        300.0,
        |state, face| {
            calls.set(calls.get() + 1);
            if calls.get() == 3 {
                Err("late source after complete stop proposal")
            } else {
                source(&f, state, face, [0.0; 2])
            }
        },
        &CancelGate::new(),
    );
    assert!(matches!(
        failure,
        Err(RigidHertzRollingCollinearHeatError::Source(
            "late source after complete stop proposal"
        ))
    ));
    let mut bounded = event_controls();
    bounded.mechanical.maximum_phases = 1;
    assert!(run(bounded, &CancelGate::new()).is_err());
    let mut bounded = event_controls();
    bounded.maximum_trials = 1;
    assert!(matches!(
        run(bounded, &CancelGate::new()),
        Err(RigidHertzRollingCollinearHeatError::TrialBudget)
    ));
    let gate = CancelGate::new();
    gate.request();
    assert!(run(event_controls(), &gate).is_err());
    let late_gate = CancelGate::new();
    let late_calls = Cell::new(0usize);
    let cancelled = propose_rigid_hertz_rolling_collinear_heat(
        &f.bodies,
        &f.properties,
        &history,
        f.contact(),
        &t,
        &route,
        bulk_binding(),
        0.0,
        5.0,
        event_controls(),
        &[0.0; 2],
        300.0,
        |state, face| {
            late_calls.set(late_calls.get() + 1);
            if late_calls.get() == 3 {
                late_gate.request();
            }
            source(&f, state, face, [0.0; 2])
        },
        &late_gate,
    );
    assert!(cancelled.is_err());
    assert_eq!(late_calls.get(), 3);
    let retry = run(event_controls(), &CancelGate::new()).unwrap();
    assert_eq!(
        baseline.mechanical.bodies_after,
        retry.mechanical.bodies_after
    );
    assert_eq!(baseline.mechanical.history, retry.mechanical.history);
    assert_eq!(baseline.source, retry.source);
    assert_eq!(history.rolling_history(), &RigidRollingHistory::new());
    assert!(!history.is_stopped());
    assert_eq!(t.time_s(), 0.0);
}
