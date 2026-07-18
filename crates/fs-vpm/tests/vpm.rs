//! Battery for the vortex particle method (fs-vpm). Each test is an analytic
//! vortex fixture: a single vortex induces a tangential Gamma/(2 pi r) field and
//! stays put, a counter-rotating pair self-propels at Gamma/(2 pi d), a
//! co-rotating pair conserves total circulation and its centroid.

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, BudgetRefusal, CancelGate, Cx, ExecMode, StreamKey, VirtualClock};
use fs_vpm::{
    VortexParticle, VpmBudget, VpmError, induced_velocity, simulate, simulate_with_cx,
    total_circulation, vorticity_centroid,
};
use std::f64::consts::{PI, TAU};
use std::mem::size_of;

const TEST_STREAM: StreamKey = StreamKey {
    seed: 0x5650_4D5F_4741_554E,
    kernel_id: 0x5650_4D,
    tile: 0,
    iteration: 0,
};

fn with_cx<R>(gate: &CancelGate, budget: Budget, f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let pool = ArenaPool::new(ArenaConfig::default());
    let clock = VirtualClock::new();
    let result = pool.scope(|arena| {
        let cx = Cx::new(gate, arena, TEST_STREAM, budget, ExecMode::Deterministic)
            .with_time_source(&clock);
        f(&cx)
    });
    assert!(
        pool.stats().quiescent(),
        "VPM test Cx arena must be quiescent"
    );
    result
}

fn roomy_budget(particles: usize, steps: usize) -> VpmBudget {
    VpmBudget::new(particles, steps, u64::MAX, usize::MAX)
}

fn assert_particles_bitwise_equal(left: &[VortexParticle], right: &[VortexParticle]) {
    assert_eq!(left.len(), right.len());
    for (left, right) in left.iter().zip(right) {
        assert_eq!(left.pos.map(f64::to_bits), right.pos.map(f64::to_bits));
        assert_eq!(left.circulation.to_bits(), right.circulation.to_bits());
    }
}

fn mag(v: [f64; 2]) -> f64 {
    (v[0] * v[0] + v[1] * v[1]).sqrt()
}

#[test]
fn a_single_vortex_induces_the_analytic_tangential_field() {
    // Gamma = 2 pi at the origin -> |u| = Gamma/(2 pi r) = 1/r, purely tangential.
    let v = [VortexParticle::new([0.0, 0.0], TAU)];
    let u1 = induced_velocity(&v, [1.0, 0.0], 0.0);
    assert!((mag(u1) - 1.0).abs() < 1e-12); // 1/r at r=1
    assert!(u1[0].abs() < 1e-12 && (u1[1] - 1.0).abs() < 1e-12); // tangential (+y)
    let u2 = induced_velocity(&v, [2.0, 0.0], 0.0);
    assert!((mag(u2) - 0.5).abs() < 1e-12); // 1/r at r=2
    // tangential everywhere: u . r = 0.
    let p = [1.3, -0.7];
    let u = induced_velocity(&v, p, 0.0);
    assert!((u[0] * p[0] + u[1] * p[1]).abs() < 1e-12);
}

#[test]
fn a_single_vortex_does_not_move_itself() {
    let start = [VortexParticle::new([0.3, -0.2], TAU)];
    let end = simulate(&start, 0.01, 100, 0.0);
    assert!((end[0].pos[0] - 0.3).abs() < 1e-12);
    assert!((end[0].pos[1] - (-0.2)).abs() < 1e-12);
}

#[test]
fn a_counter_rotating_pair_self_propels_at_the_analytic_speed() {
    // +Gamma at (0, 0.5), -Gamma at (0, -0.5): d = 1, speed Gamma/(2 pi d) = 1 in +x.
    let start = [
        VortexParticle::new([0.0, 0.5], TAU),
        VortexParticle::new([0.0, -0.5], -TAU),
    ];
    let t = 0.4;
    let end = simulate(&start, 0.01, 40, 0.0);
    // translated +x by speed*t = 0.4, vertical positions unchanged.
    assert!((end[0].pos[0] - t).abs() < 1e-6, "x = {}", end[0].pos[0]);
    assert!((end[1].pos[0] - t).abs() < 1e-6);
    assert!((end[0].pos[1] - 0.5).abs() < 1e-9 && (end[1].pos[1] + 0.5).abs() < 1e-9);
    // separation preserved (rigid translation).
    let sep = mag([end[0].pos[0] - end[1].pos[0], end[0].pos[1] - end[1].pos[1]]);
    assert!((sep - 1.0).abs() < 1e-9);
}

#[test]
fn a_co_rotating_pair_conserves_circulation_and_centroid() {
    // two +Gamma vortices orbit their shared centroid at the origin.
    let start = [
        VortexParticle::new([0.5, 0.0], TAU),
        VortexParticle::new([-0.5, 0.0], TAU),
    ];
    let c0 = total_circulation(&start);
    let end = simulate(&start, 0.005, 200, 1e-6);
    // total circulation is invariant.
    assert!((total_circulation(&end) - c0).abs() < 1e-12);
    assert!((c0 - 2.0 * TAU).abs() < 1e-12);
    // the centroid of vorticity stays at the origin.
    let centroid = vorticity_centroid(&end).unwrap();
    assert!(mag(centroid) < 1e-3, "centroid drifted to {centroid:?}");
    // they actually moved (orbiting), so this is not a trivial fixed point.
    assert!(mag([end[0].pos[0] - 0.5, end[0].pos[1]]) > 0.05);
}

#[test]
fn a_symmetric_pair_has_no_defined_centroid() {
    // total circulation zero -> centroid undefined (guarded).
    let pair = [
        VortexParticle::new([0.0, 0.5], TAU),
        VortexParticle::new([0.0, -0.5], -TAU),
    ];
    assert!(total_circulation(&pair).abs() < 1e-12);
    assert!(vorticity_centroid(&pair).is_none());
}

#[test]
fn the_desingularized_core_bounds_the_velocity() {
    // with a finite core, the self-field is finite even at the particle.
    let v = [VortexParticle::new([0.0, 0.0], TAU)];
    let u = induced_velocity(&v, [0.0, 0.0], 0.1);
    assert!(u[0].abs() < 1e-12 && u[1].abs() < 1e-12); // exactly at the particle
    // just off-core the speed is bounded by the desingularized kernel.
    let near = induced_velocity(&v, [0.01, 0.0], 0.1);
    assert!(mag(near) < 1.0 / (2.0 * PI) * 10.0); // far below the singular 1/r
}

#[test]
fn the_method_is_deterministic() {
    let start = [
        VortexParticle::new([0.5, 0.0], TAU),
        VortexParticle::new([-0.5, 0.0], TAU),
    ];
    let a = simulate(&start, 0.005, 50, 1e-6);
    let b = simulate(&start, 0.005, 50, 1e-6);
    assert_eq!(a[0].pos[0].to_bits(), b[0].pos[0].to_bits());
}

#[test]
fn g0_checked_path_matches_legacy_bits_and_reports_exact_work() {
    let start = [
        VortexParticle::new([0.0, 0.5], TAU),
        VortexParticle::new([0.0, -0.5], -TAU),
    ];
    let steps = 40;
    let legacy = simulate(&start, 0.01, steps, 0.0);
    let gate = CancelGate::new_clock_free();
    let checked = with_cx(&gate, Budget::INFINITE, |cx| {
        simulate_with_cx(cx, &start, 0.01, steps, 0.0, roomy_budget(2, steps))
    })
    .expect("valid direct RK4 plan must admit");

    assert_particles_bitwise_equal(&checked.particles, &legacy);
    assert_eq!(checked.steps_completed, steps);
    assert_eq!(checked.pair_evaluations, 4 * 2 * 2 * steps as u64);
    assert_eq!(
        checked.work_units,
        2 + steps as u64 * (4 * 2 * 2 + 4 * 2 + 1)
    );
    assert_eq!(checked.budget.planned_cost, checked.work_units);
    assert_eq!(checked.budget.cost_charged, checked.work_units);
    assert_eq!(
        checked.peak_live_bytes,
        2 * (5 * size_of::<VortexParticle>() + 4 * size_of::<[f64; 2]>())
    );
}

#[test]
fn g0_local_caps_and_checked_overflow_refuse_before_work() {
    let start = [
        VortexParticle::new([0.0, 0.5], TAU),
        VortexParticle::new([0.0, -0.5], -TAU),
    ];
    let peak_live = 2 * (5 * size_of::<VortexParticle>() + 4 * size_of::<[f64; 2]>());
    let gate = CancelGate::new_clock_free();
    with_cx(&gate, Budget::INFINITE, |cx| {
        assert_eq!(
            simulate_with_cx(
                cx,
                &start,
                0.01,
                1,
                0.0,
                VpmBudget::new(1, 1, 16, peak_live),
            ),
            Err(VpmError::CapExceeded {
                resource: "particles",
                required: 2,
                maximum: 1,
            })
        );
        assert_eq!(
            simulate_with_cx(
                cx,
                &start,
                0.01,
                1,
                0.0,
                VpmBudget::new(2, 0, 16, peak_live),
            ),
            Err(VpmError::CapExceeded {
                resource: "steps",
                required: 1,
                maximum: 0,
            })
        );
        assert_eq!(
            simulate_with_cx(
                cx,
                &start,
                0.01,
                1,
                0.0,
                VpmBudget::new(2, 1, 15, peak_live),
            ),
            Err(VpmError::CapExceeded {
                resource: "pair evaluations",
                required: 16,
                maximum: 15,
            })
        );
        assert_eq!(
            simulate_with_cx(
                cx,
                &start,
                0.01,
                1,
                0.0,
                VpmBudget::new(2, 1, 16, peak_live - 1),
            ),
            Err(VpmError::CapExceeded {
                resource: "logical live bytes",
                required: peak_live as u128,
                maximum: (peak_live - 1) as u128,
            })
        );

        let singleton = [VortexParticle::new([0.0, 0.0], 1.0)];
        assert_eq!(
            simulate_with_cx(
                cx,
                &singleton,
                0.01,
                usize::MAX,
                0.0,
                VpmBudget::new(1, usize::MAX, u64::MAX, usize::MAX),
            ),
            Err(VpmError::PlanOverflow {
                resource: "pair evaluations",
            })
        );
    });
}

#[test]
fn g0_zero_step_identity_and_invalid_numbers_are_explicit() {
    let valid = [VortexParticle::new([0.25, -0.5], 2.0)];
    let gate = CancelGate::new_clock_free();
    with_cx(&gate, Budget::INFINITE, |cx| {
        let identity = simulate_with_cx(
            cx,
            &valid,
            -0.25,
            0,
            0.0,
            VpmBudget::new(1, 0, 0, size_of::<VortexParticle>()),
        )
        .expect("zero-step signed-time identity must admit");
        assert_particles_bitwise_equal(&identity.particles, &valid);
        assert_eq!(identity.work_units, 1);
        assert_eq!(identity.pair_evaluations, 0);

        assert_eq!(
            simulate_with_cx(cx, &valid, f64::NAN, 0, 0.0, roomy_budget(1, 0)),
            Err(VpmError::InvalidScalar {
                field: "dt",
                bits: f64::NAN.to_bits(),
                particle: None,
            })
        );
        assert_eq!(
            simulate_with_cx(cx, &valid, 0.1, 0, -1.0, roomy_budget(1, 0)),
            Err(VpmError::InvalidScalar {
                field: "core",
                bits: (-1.0_f64).to_bits(),
                particle: None,
            })
        );

        let invalid = [VortexParticle::new([0.0, f64::INFINITY], 1.0)];
        assert_eq!(
            simulate_with_cx(cx, &invalid, 0.1, 0, 0.0, roomy_budget(1, 0)),
            Err(VpmError::InvalidScalar {
                field: "position.y",
                bits: f64::INFINITY.to_bits(),
                particle: Some(0),
            })
        );
    });
}

#[test]
fn g0_nonfinite_intermediate_refuses_complete_state_publication() {
    let start = [
        VortexParticle::new([f64::MAX, 0.0], 1.0),
        VortexParticle::new([-f64::MAX, 0.0], 1.0),
    ];
    let gate = CancelGate::new_clock_free();
    let error = with_cx(&gate, Budget::INFINITE, |cx| {
        simulate_with_cx(cx, &start, 0.1, 1, 0.0, roomy_budget(2, 1))
    });
    assert_eq!(
        error,
        Err(VpmError::NonFiniteResult {
            phase: "fs-vpm.rk4-k1",
            particle: 0,
        })
    );
}

#[test]
fn g4_ambient_cost_poll_and_cancellation_refuse_without_partial_state() {
    let start = [
        VortexParticle::new([0.0, 0.5], TAU),
        VortexParticle::new([0.0, -0.5], -TAU),
    ];
    let exact_work = 2 + (4 * 2 * 2 + 4 * 2 + 1);
    let active = CancelGate::new_clock_free();
    let cost_error = with_cx(
        &active,
        Budget::INFINITE.with_cost_quota(exact_work - 1),
        |cx| simulate_with_cx(cx, &start, 0.01, 1, 0.0, roomy_budget(2, 1)),
    );
    assert_eq!(
        cost_error,
        Err(VpmError::ExecutionBudget(
            BudgetRefusal::CostPlanExceedsQuota {
                planned: exact_work,
                quota: exact_work - 1,
            }
        ))
    );

    let poll_error = with_cx(&active, Budget::INFINITE.with_poll_quota(0), |cx| {
        simulate_with_cx(cx, &start, 0.01, 1, 0.0, roomy_budget(2, 1))
    });
    assert_eq!(
        poll_error,
        Err(VpmError::ExecutionBudget(BudgetRefusal::PollsExhausted {
            phase: "fs-vpm.admission",
            quota: 0,
        }))
    );

    let cancelled = CancelGate::new_clock_free();
    cancelled.request();
    let cancelled_error = with_cx(&cancelled, Budget::INFINITE, |cx| {
        simulate_with_cx(cx, &start, 0.01, 1, 0.0, roomy_budget(2, 1))
    });
    assert_eq!(
        cancelled_error,
        Err(VpmError::ExecutionBudget(BudgetRefusal::Cancelled {
            phase: "fs-vpm.admission",
        }))
    );

    let stride_fixture = [VortexParticle::new([0.0, 0.0], 0.0); 17];
    let stride_error = with_cx(&active, Budget::INFINITE.with_poll_quota(4), |cx| {
        simulate_with_cx(cx, &stride_fixture, 0.01, 1, 0.0, roomy_budget(17, 1))
    });
    assert_eq!(
        stride_error,
        Err(VpmError::ExecutionBudget(BudgetRefusal::PollsExhausted {
            phase: "fs-vpm.rk4-k1",
            quota: 4,
        }))
    );
}
