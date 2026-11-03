//! G0/G1/G3 checks of the admitted momentum-preserving VPM transaction.

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_vpm::{
    VortexParticle, VpmBudget, VpmError, conservative::simulate_with_cx, pairwise::advance_pair,
};

fn with_cx<T>(gate: &CancelGate, f: impl FnOnce(&Cx<'_>) -> T) -> T {
    let pool = ArenaPool::new(ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            gate,
            arena,
            StreamKey {
                seed: 7,
                kernel_id: 8,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        f(&cx)
    })
}

fn caps(n: usize, steps: usize) -> VpmBudget {
    VpmBudget::new(n, steps, u64::MAX, usize::MAX)
}

#[test]
fn checked_pair_matches_analytic_flow_and_owns_one_particle_buffer() {
    let pair = [
        VortexParticle::new([0.5, 0.0], 2.0),
        VortexParticle::new([-0.5, 0.0], 1.0),
    ];
    let analytic = advance_pair(pair, 0.4, 0.1).unwrap();
    with_cx(&CancelGate::new_clock_free(), |cx| {
        let run = simulate_with_cx(cx, &pair, 0.4, 1, 0.1, caps(2, 1)).unwrap();
        for (a, b) in run.particles.iter().zip(analytic) {
            for axis in 0..2 {
                assert!((a.pos[axis] - b.pos[axis]).abs() < 4e-16);
            }
        }
        assert_eq!(run.pair_evaluations, 2);
        assert_eq!(run.work_units, 5);
        assert_eq!(
            run.peak_live_bytes,
            2 * std::mem::size_of::<VortexParticle>()
        );
        assert_eq!(run.budget.cost_charged, run.work_units);
        assert_eq!(run.budget.planned_cost, run.work_units);
        let identity = simulate_with_cx(cx, &pair, 0.0, 3, 0.1, caps(2, 3)).unwrap();
        assert_eq!(identity.particles, pair);
        assert_eq!(identity.pair_evaluations, 0);
        assert_eq!(identity.work_units, 5);
    });
}

#[test]
fn admitted_signed_cloud_retains_impulse_and_is_reversible() {
    let initial = [
        VortexParticle::new([-0.7, -0.1], 1.0),
        VortexParticle::new([0.5, -0.3], -0.5),
        VortexParticle::new([0.2, 0.8], -0.5),
    ];
    let moments = |state: &[VortexParticle]| -> [f64; 3] {
        state.iter().fold([0.0; 3], |mut sum, p| {
            sum[0] += p.circulation * p.pos[0];
            sum[1] += p.circulation * p.pos[1];
            sum[2] += p.circulation * (p.pos[0] * p.pos[0] + p.pos[1] * p.pos[1]);
            sum
        })
    };
    with_cx(&CancelGate::new_clock_free(), |cx| {
        let forward = simulate_with_cx(cx, &initial, 0.03, 2000, 0.1, caps(3, 2000)).unwrap();
        for (a, b) in moments(&forward.particles).iter().zip(moments(&initial)) {
            assert!((a - b).abs() < 3e-11);
        }
        let reverse =
            simulate_with_cx(cx, &forward.particles, -0.03, 2000, 0.1, caps(3, 2000)).unwrap();
        for (a, b) in reverse.particles.iter().zip(initial) {
            for axis in 0..2 {
                assert!((a.pos[axis] - b.pos[axis]).abs() < 4e-10);
            }
        }
        assert_eq!(forward.pair_evaluations, 12000);
        assert_eq!(forward.work_units, 14003);
    });
}

#[test]
fn failed_admission_or_nonrepresentable_pair_never_publishes_a_partial_cloud() {
    let pair = [
        VortexParticle::new([0.0, 0.5], 1.0),
        VortexParticle::new([0.0, -0.5], -1.0),
    ];
    with_cx(&CancelGate::new_clock_free(), |cx| {
        assert!(matches!(
            simulate_with_cx(cx, &pair, 0.1, 1, -0.1, caps(2, 1)),
            Err(VpmError::InvalidScalar { .. })
        ));
        let limited = VpmBudget::new(2, 1, 1, usize::MAX);
        assert!(matches!(
            simulate_with_cx(cx, &pair, 0.1, 1, 0.1, limited),
            Err(VpmError::CapExceeded { .. })
        ));
        let too_close = [
            VortexParticle::new([1e-200, 0.0], 1.0),
            VortexParticle::new([0.0, 0.0], 1.0),
        ];
        assert!(matches!(
            simulate_with_cx(cx, &too_close, 1.0, 1, 0.0, caps(2, 1)),
            Err(VpmError::NonFiniteResult { .. })
        ));
    });
    let gate = CancelGate::new_clock_free();
    gate.request();
    with_cx(&gate, |cx| {
        assert!(matches!(
            simulate_with_cx(cx, &pair, 0.1, 1, 0.1, caps(2, 1)),
            Err(VpmError::ExecutionBudget(_))
        ));
    });
    assert_eq!(pair[0].pos, [0.0, 0.5]);
    assert_eq!(pair[1].pos, [0.0, -0.5]);
}
