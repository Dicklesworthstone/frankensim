//! Integration tests for multirate stepping with conserved interface
//! accumulators (bead `frankensim-rc-1109-root-3gc12.5.9`, Package 5a).

#![allow(clippy::float_cmp)]

use core::cell::Cell;
use std::rc::Rc;

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_evidence::Color;
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_system::elements::{
    Damper, Force, GasContents, GasVolume, HeatCapacity, Mass, Motion, Restriction, SourceAccount,
    Spring, Waveform,
};
use fs_system::{
    EventDirection, EventSetConfig, EventSetStatus, MultirateConfig, MultirateError, NodeId,
    RouteTemperature, SlowInterpolation, System, SystemBuilder, SystemError,
};

fn with_cx<R>(gate: &CancelGate, f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let pool = ArenaPool::new(ArenaConfig::default());
    pool.scope(|arena| {
        f(&Cx::new(
            gate,
            arena,
            StreamKey {
                seed: 59,
                kernel_id: 9,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        ))
    })
}

fn close(actual: f64, expected: f64, tol: f64) {
    assert!(
        (actual - expected).abs() <= tol,
        "actual={actual:.17e} expected={expected:.17e} diff={:.3e} tol={tol:.3e}",
        (actual - expected).abs()
    );
}

/// Build a fast damped oscillator (`mass`, `spring`, `damper`) heating a slow
/// `HeatCapacity` (`body`) through `HeatRoute`.
fn build_oscillator_heating_body(
    route_temp: RouteTemperature,
) -> (System, NodeId, NodeId, NodeId, NodeId) {
    let mut b = SystemBuilder::new();
    let mass = b
        .add_with_tangent(
            "mass",
            Mass::new(0.02).unwrap(),
            Motion {
                position_m: 0.05,
                velocity_m_s: 0.0,
            },
        )
        .unwrap();
    let spring = b
        .add_with_tangent("spring", Spring::new(2.0e3, 0.0).unwrap(), 0.05)
        .unwrap();
    let damper = b
        .add_with_tangent("damper", Damper::new(0.4).unwrap(), ())
        .unwrap();
    let body = b
        .add_with_tangent("body", HeatCapacity::new(50.0).unwrap(), 300.0)
        .unwrap();

    let mp = b.port(mass, "p").unwrap();
    let sa = b.port(spring, "a").unwrap();
    let sb = b.port(spring, "b").unwrap();
    let da = b.port(damper, "a").unwrap();
    let db = b.port(damper, "b").unwrap();
    let dheat = b.port(damper, "heat").unwrap();
    let bheat = b.port(body, "heat").unwrap();

    b.common_flow(mp, &[sa, da]).unwrap();
    b.boundary(sb, Box::new(|_, _, out| out[0] = 0.0)).unwrap();
    b.boundary(db, Box::new(|_, _, out| out[0] = 0.0)).unwrap();
    b.heat_route(dheat, &[(bheat, 1.0)], None, route_temp)
        .unwrap();

    let mut sys = b.seal().unwrap();
    if matches!(route_temp, RouteTemperature::SinkEquivalent) {
        let _ = sys.guess(dheat, &[300.0]);
    }
    (sys, mass, spring, damper, body)
}

#[test]
fn rate_ratio_one_reproduces_single_rate_bitwise() {
    let (mut single, m1, _s1, _d1, b1) =
        build_oscillator_heating_body(RouteTemperature::SinkEquivalent);
    let (mut multi, m2, _s2, _d2, b2) =
        build_oscillator_heating_body(RouteTemperature::SinkEquivalent);

    let plan = multi
        .plan_multirate(1.0e-3, MultirateConfig::default())
        .unwrap();
    assert!(plan.is_single_rate());
    assert_eq!(plan.groups().len(), 1);
    assert_eq!(plan.groups()[0].substeps, 1);

    let gate = CancelGate::new_clock_free();
    with_cx(&gate, |cx| {
        for _ in 0..40 {
            let rs = single.step(1.0e-3, cx).unwrap();
            let rm = multi.step_multirate(1.0e-3, &plan, cx).unwrap();
            assert_eq!(rm.dt_s.to_bits(), rs.dt_s.to_bits());
            assert_eq!(rm.time_s.to_bits(), rs.time_s.to_bits());
            assert_eq!(
                rm.interpolation_numerical_loss_j.to_bits(),
                0.0_f64.to_bits()
            );
            assert_eq!(rm.cross_rate_iterations, 0);
        }
    });

    assert_eq!(single.time_s().to_bits(), multi.time_s().to_bits());
    assert_eq!(single.steps(), multi.steps());

    let ms1 = single.state::<Motion>(m1).unwrap();
    let ms2 = multi.state::<Motion>(m2).unwrap();
    assert_eq!(ms1.position_m.to_bits(), ms2.position_m.to_bits());
    assert_eq!(ms1.velocity_m_s.to_bits(), ms2.velocity_m_s.to_bits());

    let t1 = *single.state::<f64>(b1).unwrap();
    let t2 = *multi.state::<f64>(b2).unwrap();
    assert_eq!(t1.to_bits(), t2.to_bits());

    let ls = single.ledger();
    let lm = multi.ledger();
    assert_eq!(
        ls.stored_change_j().to_bits(),
        lm.stored_change_j().to_bits()
    );
    assert_eq!(
        ls.total_dissipated_j().to_bits(),
        lm.total_dissipated_j().to_bits()
    );
    assert_eq!(
        ls.energy_residual_j().to_bits(),
        lm.energy_residual_j().to_bits()
    );
    assert_eq!(
        ls.connection_defect_j().to_bits(),
        lm.connection_defect_j().to_bits()
    );
    assert_eq!(
        ls.entropy_generated_j_per_k().to_bits(),
        lm.entropy_generated_j_per_k().to_bits()
    );
}

#[test]
#[allow(clippy::too_many_lines)]
fn fast_oscillator_heating_slow_thermal_mass_conserves_energy_and_converges() {
    let gate = CancelGate::new_clock_free();
    let horizon_s = 0.08_f64;
    let dt_fast = 2.5e-4_f64;
    let total_fast_steps = (horizon_s / dt_fast).round() as usize;

    // Fine single-rate reference at `dt_fast` with two-way `SinkEquivalent` route.
    let (mut reference, mr, _sr, _dr, br) =
        build_oscillator_heating_body(RouteTemperature::SinkEquivalent);
    with_cx(&gate, |cx| {
        for _ in 0..total_fast_steps {
            reference.step(dt_fast, cx).unwrap();
        }
    });
    let ref_motion = *reference.state::<Motion>(mr).unwrap();
    let ref_temp = *reference.state::<f64>(br).unwrap();
    let ref_diss = reference.ledger().total_dissipated_j();
    assert!(ref_temp > 300.0 + 1.0e-3, "body must noticeably warm");

    // Test BOTH one-way (`Fixed(350 K)`) and two-way (`SinkEquivalent`, cyclic cross-rate Newton).
    for (label, route_temp) in [
        ("one-way", RouteTemperature::Fixed(350.0)),
        ("two-way", RouteTemperature::SinkEquivalent),
    ] {
        let mut prev_temp_err = f64::INFINITY;
        for macro_substeps in [16_usize, 8_usize] {
            let macro_dt = dt_fast * (macro_substeps as f64);
            let macro_steps = total_fast_steps / macro_substeps;
            let (mut sys, m, s, d, body) = build_oscillator_heating_body(route_temp);
            let config = MultirateConfig {
                node_substeps: vec![
                    (m, macro_substeps),
                    (s, macro_substeps),
                    (d, macro_substeps),
                    (body, 1),
                ],
                interpolation: SlowInterpolation::ZeroOrderHold,
                ..MultirateConfig::default()
            };
            let plan = sys.plan_multirate(macro_dt, config).unwrap();
            assert!(!plan.is_single_rate());
            assert_eq!(plan.groups().len(), 2);
            assert_eq!(plan.groups()[0].substeps, macro_substeps);
            assert_eq!(plan.groups()[1].substeps, 1);
            if matches!(route_temp, RouteTemperature::SinkEquivalent) {
                assert_eq!(plan.cyclic_stage_count(), 1);
            } else {
                assert_eq!(plan.cyclic_stage_count(), 0);
            }

            with_cx(&gate, |cx| {
                for _ in 0..macro_steps {
                    let rep = sys.step_multirate(macro_dt, &plan, cx).unwrap();
                    assert_eq!(rep.interpolation_numerical_loss_j, 0.0);
                    assert_eq!(rep.group_substeps, vec![macro_substeps, 1]);
                }
            });

            close(sys.time_s(), horizon_s, 1.0e-12);
            let mot = *sys.state::<Motion>(m).unwrap();
            let temp = *sys.state::<f64>(body).unwrap();
            let diss = sys.ledger().total_dissipated_j();

            // Damper force is temperature-independent, so fast mechanical
            // trajectory matches the fine single-rate reference closely.
            close(mot.position_m, ref_motion.position_m, 1.0e-11);
            close(mot.velocity_m_s, ref_motion.velocity_m_s, 1.0e-9);
            close(diss, ref_diss, 1.0e-10);
            let temp_err = (temp - ref_temp).abs();
            assert!(
                temp_err <= 1.0e-8,
                "{label} M={macro_substeps}: temp={temp} vs ref={ref_temp} err={temp_err:.3e}"
            );
            assert!(
                temp_err <= prev_temp_err + 1.0e-14,
                "{label}: refining macro step should not worsen temperature error"
            );
            prev_temp_err = temp_err;

            // Exact energy & entropy ledger closure across the rate boundary.
            let ledger = sys.ledger();
            close(ledger.energy_residual_j(), 0.0, 1.0e-13);
            close(ledger.connection_defect_j(), 0.0, 1.0e-13);
            close(ledger.numerical_loss_j(), 0.0, 1.0e-15);
            assert!(
                ledger.entropy_generated_j_per_k() >= -1.0e-15,
                "{label}: second law violated"
            );
            let outcome = ledger.outcome();
            assert!(
                matches!(outcome.energy, Color::Verified { .. }),
                "{label}: ZOH energy outcome must be Verified, got {:?}",
                outcome.energy
            );
            assert!(
                matches!(outcome.overall, Color::Verified { .. }),
                "{label}: ZOH overall outcome must be Verified, got {:?}",
                outcome.overall
            );
        }
    }
}

/// Build a two-way mechanically coupled fast+slow system:
/// - Fast oscillator (`fast_mass`, `coupling_spring`) attached to a moving
///   `slow_mass` driven by an external `Force`, coupled via lossless
///   `CommonFlow` across the rate boundary.
fn build_two_way_mechanical_pair() -> (System, NodeId, NodeId, NodeId, NodeId) {
    let mut b = SystemBuilder::new();
    let fast_mass = b
        .add_with_tangent(
            "fast_mass",
            Mass::new(0.05).unwrap(),
            Motion {
                position_m: 0.0,
                velocity_m_s: 0.4,
            },
        )
        .unwrap();
    let coupling_spring = b
        .add_with_tangent("coupling_spring", Spring::new(400.0, 0.0).unwrap(), 0.02)
        .unwrap();
    let slow_mass = b
        .add_with_tangent(
            "slow_mass",
            Mass::new(2.0).unwrap(),
            Motion {
                position_m: 0.0,
                velocity_m_s: 0.1,
            },
        )
        .unwrap();
    let drive = b
        .add_with_tangent(
            "drive",
            Force::new(Waveform::constant(6.0).unwrap()),
            SourceAccount::default(),
        )
        .unwrap();

    let fm_p = b.port(fast_mass, "p").unwrap();
    let cs_a = b.port(coupling_spring, "a").unwrap();
    let cs_b = b.port(coupling_spring, "b").unwrap();
    let sm_p = b.port(slow_mass, "p").unwrap();
    let dr_p = b.port(drive, "p").unwrap();

    b.common_flow(fm_p, &[cs_a]).unwrap();
    b.common_flow(sm_p, &[cs_b, dr_p]).unwrap();

    (
        b.seal().unwrap(),
        fast_mass,
        coupling_spring,
        slow_mass,
        drive,
    )
}

#[test]
#[allow(clippy::too_many_lines)]
fn slow_to_fast_linear_interpolation_books_exact_numerical_loss() {
    let gate = CancelGate::new_clock_free();
    let macro_dt = 4.0e-3;
    let steps = 25_usize;

    // 1. ZeroOrderHold: zero interpolation loss, exact conservation, Verified.
    let (mut sys_zoh, fm1, cs1, sm1, dr1) = build_two_way_mechanical_pair();
    let plan_zoh = sys_zoh
        .plan_multirate(
            macro_dt,
            MultirateConfig {
                node_substeps: vec![(fm1, 8), (cs1, 8), (sm1, 1), (dr1, 1)],
                interpolation: SlowInterpolation::ZeroOrderHold,
                ..MultirateConfig::default()
            },
        )
        .unwrap();
    assert_eq!(plan_zoh.cyclic_stage_count(), 1);
    let mut total_cross_iters = 0_u32;
    with_cx(&gate, |cx| {
        for _ in 0..steps {
            let rep = sys_zoh.step_multirate(macro_dt, &plan_zoh, cx).unwrap();
            assert_eq!(rep.interpolation_numerical_loss_j, 0.0);
            total_cross_iters += rep.cross_rate_iterations;
        }
    });
    assert!(
        total_cross_iters >= steps as u32,
        "cyclic two-way mechanical pair must exercise cross-rate Newton iterations"
    );
    let lz = sys_zoh.ledger();
    close(lz.numerical_loss_j(), 0.0, 1.0e-15);
    close(lz.energy_residual_j(), 0.0, 1.0e-12);
    close(lz.connection_defect_j(), 0.0, 1.0e-12);
    assert!(matches!(lz.outcome().energy, Color::Verified { .. }));

    // 2. Linear interpolation: nonzero quadrature discrepancy explicitly booked
    //    into `interpolation_numerical_loss_j` and `ledger.numerical_loss_j()`,
    //    keeping `energy_residual_j == 0` while refusing `Color::Verified`.
    let (mut sys_lin, fm2, cs2, sm2, dr2) = build_two_way_mechanical_pair();
    let plan_lin = sys_lin
        .plan_multirate(
            macro_dt,
            MultirateConfig {
                node_substeps: vec![(fm2, 8), (cs2, 8), (sm2, 1), (dr2, 1)],
                interpolation: SlowInterpolation::Linear,
                ..MultirateConfig::default()
            },
        )
        .unwrap();
    let mut booked_total = 0.0_f64;
    let mut nonzero_steps = 0_usize;
    with_cx(&gate, |cx| {
        for step_idx in 0..steps {
            let rep = sys_lin.step_multirate(macro_dt, &plan_lin, cx).unwrap();
            if step_idx == 0 {
                assert_eq!(rep.interpolation_numerical_loss_j, 0.0);
            } else if rep.interpolation_numerical_loss_j.abs() > 1.0e-12 {
                nonzero_steps += 1;
            }
            booked_total += rep.interpolation_numerical_loss_j;
        }
    });
    assert!(
        nonzero_steps >= 10,
        "linear slow-to-fast interpolation must produce measurable booked numerical loss, got {nonzero_steps} nonzero steps (total={booked_total:.3e})"
    );
    let ll = sys_lin.ledger();
    close(ll.numerical_loss_j(), booked_total, 1.0e-12);
    close(ll.energy_residual_j(), 0.0, 1.0e-12);
    close(ll.connection_defect_j(), 0.0, 1.0e-12);
    assert!(
        matches!(
            &ll.outcome().energy,
            Color::Estimated { estimator, .. } if estimator == "window-ledger:numerical-loss"
        ),
        "linear interpolation numerical loss must color window Estimated, got {:?}",
        ll.outcome().energy
    );
}

#[test]
fn cross_rate_stream_mass_enthalpy_and_carried_conservation() {
    let mut b = SystemBuilder::new();
    let pipe = b
        .add_with_tangent("pipe", Restriction::new(2.0e-8, 1).unwrap(), ())
        .unwrap();
    let gas_a = GasVolume::new(287.0, 1.4, 1).unwrap();
    let contents_a = gas_a.contents(2.2e5, 310.0, 0.005, &[0.35]).unwrap();
    let gas_b = GasVolume::new(287.0, 1.4, 1).unwrap();
    let contents_b = gas_b.contents(1.8e5, 300.0, 0.05, &[0.10]).unwrap();
    let total_mass_0 = contents_a.mass_kg + contents_b.mass_kg;
    let total_tracer_0 = contents_a.carried_kg[0] + contents_b.carried_kg[0];

    let va = b
        .add_with_tangent("small_vessel", gas_a, contents_a)
        .unwrap();
    let vb = b
        .add_with_tangent("large_vessel", gas_b, contents_b)
        .unwrap();

    b.bond(b.port(va, "stream").unwrap(), b.port(pipe, "a").unwrap())
        .unwrap();
    b.bond(b.port(vb, "stream").unwrap(), b.port(pipe, "b").unwrap())
        .unwrap();
    b.boundary(
        b.port(va, "heat").unwrap(),
        Box::new(|t0, dt, out| {
            let t = t0 + 0.5 * dt;
            out[0] = if t < 0.08 { 400.0 } else { -300.0 };
        }),
    )
    .unwrap();
    b.boundary(
        b.port(vb, "heat").unwrap(),
        Box::new(|_, _, out| out[0] = 0.0),
    )
    .unwrap();

    let mut sys = b.seal().unwrap();
    let macro_dt = 0.01;
    let plan = sys
        .plan_multirate(
            macro_dt,
            MultirateConfig {
                node_substeps: vec![(va, 8), (pipe, 8), (vb, 1)],
                interpolation: SlowInterpolation::ZeroOrderHold,
                ..MultirateConfig::default()
            },
        )
        .unwrap();
    assert_eq!(plan.groups().len(), 2);

    let gate = CancelGate::new_clock_free();
    with_cx(&gate, |cx| {
        for _ in 0..20 {
            sys.step_multirate(macro_dt, &plan, cx).unwrap();
        }
    });

    let ca = sys.state::<GasContents>(va).unwrap();
    let cb = sys.state::<GasContents>(vb).unwrap();
    close(ca.mass_kg + cb.mass_kg, total_mass_0, 1.0e-13);
    close(ca.carried_kg[0] + cb.carried_kg[0], total_tracer_0, 1.0e-13);

    let ledger = sys.ledger();
    close(ledger.mass_residual_kg().unwrap(), 0.0, 1.0e-13);
    close(ledger.carried_residual(0).unwrap(), 0.0, 1.0e-13);
    close(ledger.energy_residual_j(), 0.0, 1.0e-11);
    close(ledger.connection_defect_j(), 0.0, 1.0e-11);
    assert!(matches!(ledger.outcome().energy, Color::Verified { .. }));
    assert!(matches!(
        ledger.outcome().mass,
        Some(Color::Verified { .. })
    ));
}

#[test]
fn cross_rate_fast_and_slow_events_synchronize_all_clocks() {
    let (mut sys, mass, spring, damper, body) =
        build_oscillator_heating_body(RouteTemperature::SinkEquivalent);
    let macro_dt = 8.0e-3;
    let plan = sys
        .plan_multirate(
            macro_dt,
            MultirateConfig {
                node_substeps: vec![(mass, 8), (spring, 8), (damper, 8), (body, 1)],
                interpolation: SlowInterpolation::ZeroOrderHold,
                ..MultirateConfig::default()
            },
        )
        .unwrap();

    let gate = CancelGate::new_clock_free();
    let cfg = EventSetConfig {
        time_tolerance_s: 1.0e-7,
        max_bisections: 48,
    };

    // Guard 0: fast spring stretch crosses 0.0 (falling).
    // Guard 1: slow body crosses T = 300.015 K (rising).
    let mut located_fast = false;
    let mut located_slow = false;
    with_cx(&gate, |cx| {
        for _ in 0..30 {
            let rep = sys
                .step_multirate_to_events(
                    macro_dt,
                    &plan,
                    cfg,
                    &[EventDirection::Falling, EventDirection::Rising],
                    |s, vals| {
                        vals[0] = *s.state::<f64>(spring).unwrap();
                        vals[1] = *s.state::<f64>(body).unwrap() - 300.015;
                    },
                    cx,
                )
                .unwrap();
            assert_eq!(rep.time_s.to_bits(), sys.time_s().to_bits());
            if let Some(step) = &rep.step {
                assert_eq!(step.time_s.to_bits(), sys.time_s().to_bits());
            }
            if let EventSetStatus::Located(events) = rep.status {
                for ev in events {
                    assert!(
                        ev.bracket.upper_dt_s - ev.bracket.lower_dt_s <= cfg.time_tolerance_s,
                        "bracket wider than tolerance"
                    );
                    if ev.index == 0 {
                        located_fast = true;
                        assert!(ev.bracket.upper_value <= 0.0);
                    }
                    if ev.index == 1 {
                        located_slow = true;
                        assert!(ev.bracket.upper_value >= 0.0);
                    }
                }
                if located_fast && located_slow {
                    break;
                }
            }
        }
    });
    assert!(located_fast, "must locate fast-group position crossing");
    assert!(located_slow, "must locate slow-group temperature crossing");
    close(sys.ledger().energy_residual_j(), 0.0, 1.0e-12);
}

#[test]
fn cross_rate_slow_guard_alone_synchronizes_all_clocks() {
    let (mut sys, mass, spring, damper, body) =
        build_oscillator_heating_body(RouteTemperature::SinkEquivalent);
    let macro_dt = 8.0e-3;
    let plan = sys
        .plan_multirate(
            macro_dt,
            MultirateConfig {
                node_substeps: vec![(mass, 8), (spring, 8), (damper, 8), (body, 1)],
                interpolation: SlowInterpolation::ZeroOrderHold,
                ..MultirateConfig::default()
            },
        )
        .unwrap();

    let gate = CancelGate::new_clock_free();
    let cfg = EventSetConfig {
        time_tolerance_s: 1.0e-7,
        max_bisections: 48,
    };
    let mut located_slow = false;
    with_cx(&gate, |cx| {
        for _ in 0..30 {
            let rep = sys
                .step_multirate_to_events(
                    macro_dt,
                    &plan,
                    cfg,
                    &[EventDirection::Rising],
                    |s, vals| {
                        vals[0] = *s.state::<f64>(body).unwrap() - 300.015;
                    },
                    cx,
                )
                .unwrap();
            assert_eq!(rep.time_s.to_bits(), sys.time_s().to_bits());
            if let EventSetStatus::Located(events) = rep.status {
                assert_eq!(events.len(), 1);
                let ev = &events[0];
                assert!(
                    ev.bracket.upper_dt_s - ev.bracket.lower_dt_s <= cfg.time_tolerance_s,
                    "slow-guard bracket wider than tolerance"
                );
                assert!(ev.bracket.upper_value >= 0.0);
                located_slow = true;
                break;
            }
        }
    });
    assert!(located_slow, "must locate slow-group temperature crossing");
    close(sys.ledger().energy_residual_j(), 0.0, 1.0e-12);
}

#[test]
fn multirate_refusal_rolls_back_atomically() {
    // Mid-step refusal rollback: cyclic mechanical pair with `cross_rate_max_iterations = 0`
    // advances ALL fast substeps (mutating working `self.accepted`) before refusing
    // with `CrossRateNotConverged`; rollback must restore entry state bit-for-bit.
    let gate = CancelGate::new_clock_free();
    let (mut ref_sys, fm, cs, sm, dr) = build_two_way_mechanical_pair();
    let refusing_plan = ref_sys
        .plan_multirate(
            4.0e-3,
            MultirateConfig {
                node_substeps: vec![(fm, 8), (cs, 8), (sm, 1), (dr, 1)],
                cross_rate_max_iterations: 0,
                ..MultirateConfig::default()
            },
        )
        .unwrap();
    let fm_before = *ref_sys.state::<Motion>(fm).unwrap();
    let sm_before = *ref_sys.state::<Motion>(sm).unwrap();
    with_cx(&gate, |cx| {
        let err = ref_sys
            .step_multirate(4.0e-3, &refusing_plan, cx)
            .unwrap_err();
        assert!(matches!(
            err,
            MultirateError::CrossRateNotConverged { iterations: 0, .. }
        ));
    });
    let fm_after = *ref_sys.state::<Motion>(fm).unwrap();
    let sm_after = *ref_sys.state::<Motion>(sm).unwrap();
    assert_eq!(
        fm_after.position_m.to_bits(),
        fm_before.position_m.to_bits()
    );
    assert_eq!(
        fm_after.velocity_m_s.to_bits(),
        fm_before.velocity_m_s.to_bits()
    );
    assert_eq!(
        sm_after.position_m.to_bits(),
        sm_before.position_m.to_bits()
    );
    assert_eq!(
        sm_after.velocity_m_s.to_bits(),
        sm_before.velocity_m_s.to_bits()
    );
    assert_eq!(ref_sys.steps(), 0);
    assert_eq!(ref_sys.time_s().to_bits(), 0.0_f64.to_bits());

    // Foreign plan, RateRatioExceeded, and OverrideViolatesLimit refusals.
    let (mut sys, mass, ..) = build_oscillator_heating_body(RouteTemperature::SinkEquivalent);
    let macro_dt = 8.0e-3;
    let (other_sys, ..) = build_oscillator_heating_body(RouteTemperature::Fixed(300.0));
    let other_plan = other_sys
        .plan_multirate(macro_dt, MultirateConfig::default())
        .unwrap();
    with_cx(&gate, |cx| {
        let err = sys.step_multirate(macro_dt, &other_plan, cx).unwrap_err();
        assert!(matches!(err, MultirateError::ForeignPlan));
    });

    let err_ratio = sys
        .plan_multirate(
            1.0,
            MultirateConfig {
                max_substeps: 4,
                node_max_dt_s: vec![(mass, 1.0e-3)],
                ..MultirateConfig::default()
            },
        )
        .unwrap_err();
    assert!(matches!(
        err_ratio,
        MultirateError::RateRatioExceeded { .. }
    ));

    let err_override = sys
        .plan_multirate(
            0.01,
            MultirateConfig {
                max_substeps: 32,
                node_max_dt_s: vec![(mass, 1.0e-3)],
                node_substeps: vec![(mass, 2)],
                ..MultirateConfig::default()
            },
        )
        .unwrap_err();
    assert!(matches!(
        err_override,
        MultirateError::OverrideViolatesLimit { .. }
    ));
}

#[test]
fn multirate_cancellation_rolls_back_atomically() {
    // Mid-step cancellation rollback: boundary closure cancels the gate on fast
    // substep 4 AFTER substeps 0..3 have already mutated working fast states;
    // rollback must restore entry state, clock, and ledger bit-for-bit.
    let mid_gate = Rc::new(CancelGate::new_clock_free());
    let calls = Rc::new(Cell::new(0_usize));
    let mut b = SystemBuilder::new();
    let cm = b
        .add_with_tangent(
            "mass",
            Mass::new(0.02).unwrap(),
            Motion {
                position_m: 0.05,
                velocity_m_s: 0.0,
            },
        )
        .unwrap();
    let cspring = b
        .add_with_tangent("spring", Spring::new(2.0e3, 0.0).unwrap(), 0.05)
        .unwrap();
    let cd = b
        .add_with_tangent("damper", Damper::new(0.4).unwrap(), ())
        .unwrap();
    let cbody = b
        .add_with_tangent("body", HeatCapacity::new(50.0).unwrap(), 300.0)
        .unwrap();
    b.common_flow(
        b.port(cm, "p").unwrap(),
        &[b.port(cspring, "a").unwrap(), b.port(cd, "a").unwrap()],
    )
    .unwrap();
    let gate_hook = Rc::clone(&mid_gate);
    let calls_hook = Rc::clone(&calls);
    b.boundary(
        b.port(cspring, "b").unwrap(),
        Box::new(move |_, _, out| {
            let n = calls_hook.get() + 1;
            calls_hook.set(n);
            if n == 4 {
                gate_hook.request();
            }
            out[0] = 0.0;
        }),
    )
    .unwrap();
    b.boundary(b.port(cd, "b").unwrap(), Box::new(|_, _, out| out[0] = 0.0))
        .unwrap();
    b.heat_route(
        b.port(cd, "heat").unwrap(),
        &[(b.port(cbody, "heat").unwrap(), 1.0)],
        None,
        RouteTemperature::Fixed(300.0),
    )
    .unwrap();
    let mut cancel_sys = b.seal().unwrap();
    let cancel_plan = cancel_sys
        .plan_multirate(
            8.0e-3,
            MultirateConfig {
                node_substeps: vec![(cm, 8), (cspring, 8), (cd, 8), (cbody, 1)],
                ..MultirateConfig::default()
            },
        )
        .unwrap();
    let cm_before = *cancel_sys.state::<Motion>(cm).unwrap();
    with_cx(&mid_gate, |cx| {
        let err = cancel_sys
            .step_multirate(8.0e-3, &cancel_plan, cx)
            .unwrap_err();
        assert!(matches!(
            err,
            MultirateError::System(ref e) if matches!(**e, SystemError::Cancelled)
        ));
    });
    assert!(
        calls.get() >= 4,
        "cancellation must trigger mid-step after fast substeps already ran"
    );
    let cm_after = *cancel_sys.state::<Motion>(cm).unwrap();
    assert_eq!(
        cm_after.position_m.to_bits(),
        cm_before.position_m.to_bits()
    );
    assert_eq!(
        cm_after.velocity_m_s.to_bits(),
        cm_before.velocity_m_s.to_bits()
    );
    assert_eq!(cancel_sys.steps(), 0);
    assert_eq!(cancel_sys.time_s().to_bits(), 0.0_f64.to_bits());
}
