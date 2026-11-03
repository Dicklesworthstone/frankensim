//! G0/G1/G3 checks of the production exact-pair vortex map.

use fs_vpm::{
    VortexParticle,
    pairwise::{PairFlowError, advance_pair},
};
use std::f64::consts::{PI, TAU};

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual} != {expected}, tol={tolerance}"
    );
}

fn moments(particles: &[VortexParticle]) -> [f64; 3] {
    particles.iter().fold([0.0; 3], |mut sum, p| {
        sum[0] += p.circulation * p.pos[0];
        sum[1] += p.circulation * p.pos[1];
        sum[2] += p.circulation * (p.pos[0] * p.pos[0] + p.pos[1] * p.pos[1]);
        sum
    })
}

fn compose_pairs(state: &mut [VortexParticle], dt: f64, steps: usize, core: f64) {
    for _ in 0..steps {
        for i in 0..state.len() {
            for j in (i + 1)..state.len() {
                let pair = advance_pair([state[i], state[j]], 0.5 * dt, core).unwrap();
                state[i] = pair[0];
                state[j] = pair[1];
            }
        }
        for i in (0..state.len()).rev() {
            for j in ((i + 1)..state.len()).rev() {
                let pair = advance_pair([state[i], state[j]], 0.5 * dt, core).unwrap();
                state[i] = pair[0];
                state[j] = pair[1];
            }
        }
    }
}

#[test]
fn symmetric_composition_conserves_mixed_cloud_impulse_and_converges_quadratically() {
    let initial = [
        VortexParticle::new([-0.7, -0.1], 1.0),
        VortexParticle::new([0.5, -0.3], -0.5),
        VortexParticle::new([0.2, 0.8], -0.5),
    ];
    let reference = fs_vpm::simulate(&initial, 0.0001, 10000, 0.1);
    let mut coarse = initial;
    let mut fine = initial;
    compose_pairs(&mut coarse, 0.05, 20, 0.1);
    compose_pairs(&mut fine, 0.025, 40, 0.1);
    let error = |state: &[VortexParticle]| -> f64 {
        state
            .iter()
            .zip(&reference)
            .map(|(a, b)| (a.pos[0] - b.pos[0]).powi(2) + (a.pos[1] - b.pos[1]).powi(2))
            .sum::<f64>()
            .sqrt()
    };
    assert!((3.8..4.2).contains(&(error(&coarse) / error(&fine))));
    let before = moments(&initial);
    let mut long = initial;
    compose_pairs(&mut long, 0.03, 2000, 0.1);
    for (actual, expected) in moments(&long).iter().zip(before) {
        close(*actual, expected, 3e-11);
    }
    compose_pairs(&mut long, -0.03, 2000, 0.1);
    for (a, b) in long.iter().zip(initial) {
        for axis in 0..2 {
            close(a.pos[axis], b.pos[axis], 4e-10);
        }
    }
}

#[test]
fn equal_vortices_rotate_exactly_with_the_declared_softened_rate() {
    let pair = [
        VortexParticle::new([0.5, 0.0], TAU),
        VortexParticle::new([-0.5, 0.0], TAU),
    ];
    for core in [0.0, 0.2, 3.0] {
        let time = 0.5 * PI * (1.0 + core * core);
        let end = advance_pair(pair, time, core).unwrap();
        close(end[0].pos[0], -0.5, 6e-16);
        close(end[1].pos[0], 0.5, 6e-16);
        close(end[0].pos[1], 0.0, 6e-16);
        close(end[1].pos[1], 0.0, 6e-16);
        assert_eq!(pair[0].circulation.to_bits(), end[0].circulation.to_bits());
    }
}

#[test]
fn exact_dipole_and_nearly_cancelling_strengths_translate_continuously() {
    let pair = [
        VortexParticle::new([0.0, 0.5], TAU),
        VortexParticle::new([0.0, -0.5], -TAU),
    ];
    for core in [0.0, 0.25] {
        let end = advance_pair(pair, 0.4, core).unwrap();
        for (i, p) in end.iter().enumerate() {
            close(p.pos[0], 0.4 / (1.0 + core * core), 2e-16);
            close(p.pos[1], pair[i].pos[1], 0.0);
        }
        let mut nearby = pair;
        nearby[1].circulation *= 1.0 - 1e-13;
        let almost = advance_pair(nearby, 0.4, core).unwrap();
        for (a, b) in end.iter().zip(almost) {
            close(a.pos[0], b.pos[0], 1e-13);
            close(a.pos[1], b.pos[1], 1e-13);
        }
    }
}

#[test]
fn signed_pairs_preserve_all_three_impulse_moments_and_reverse() {
    for strengths in [
        [3.0, 7.0],
        [3.0, -7.0],
        [-3.0, -7.0],
        [3.0, -3.0],
        [0.0, 7.0],
    ] {
        let start = [
            VortexParticle::new([1.75, -0.25], strengths[0]),
            VortexParticle::new([-0.75, 1.25], strengths[1]),
        ];
        let before = moments(&start);
        let mut state = start;
        for _ in 0..400 {
            state = advance_pair(state, 0.12, 0.2).unwrap();
        }
        for (actual, expected) in moments(&state).iter().zip(before) {
            close(*actual, expected, 4e-11);
        }
        for _ in 0..400 {
            state = advance_pair(state, -0.12, 0.2).unwrap();
        }
        for (actual, expected) in state.iter().zip(start) {
            for axis in 0..2 {
                close(actual.pos[axis], expected.pos[axis], 3e-11);
            }
        }
    }
}

#[test]
fn translating_and_rotating_the_frame_commutes_with_pair_flow() {
    let start = [
        VortexParticle::new([0.3, -0.2], 2.0),
        VortexParticle::new([-0.7, 0.5], -0.8),
    ];
    let transform =
        |p: VortexParticle| VortexParticle::new([8.0 - p.pos[1], -4.0 + p.pos[0]], p.circulation);
    let base = advance_pair(start, 0.3, 0.1).unwrap().map(transform);
    let moved = advance_pair(start.map(transform), 0.3, 0.1).unwrap();
    for (a, b) in base.iter().zip(moved) {
        for axis in 0..2 {
            close(a.pos[axis], b.pos[axis], 3e-15);
        }
    }
}

#[test]
fn dimensional_scaling_survives_unrepresentable_squared_lengths() {
    let base = [
        VortexParticle::new([0.5, 0.0], 1.0),
        VortexParticle::new([-0.5, 0.0], -0.75),
    ];
    // Time scales by L² with circulation fixed. Squared lengths or a naive
    // quotient under/overflow, although angles and final positions are finite.
    for scale in [1e-160, 1e-150, 1e150] {
        let scaled = base.map(|p| VortexParticle::new(p.pos.map(|x| scale * x), p.circulation));
        let scaled_time = 0.4 * scale * scale;
        let result = advance_pair(scaled, scaled_time, 0.3 * scale).unwrap();
        // A subnormal time has coarser input quantization; compare its actual
        // represented physical duration, not the unrepresentable nominal one.
        let represented = advance_pair(base, scaled_time / scale / scale, 0.3).unwrap();
        for (a, b) in represented.iter().zip(result) {
            for axis in 0..2 {
                close(b.pos[axis] / scale, a.pos[axis], 2e-15);
            }
        }
    }
    // Γ rescales by L² instead. The summed strength exceeds f64::MAX while
    // the angular rate and each strength remain representable.
    let large = [
        VortexParticle::new([5e153, 0.0], 1e308),
        VortexParticle::new([-5e153, 0.0], 1e308),
    ];
    let evolved = advance_pair(large, 0.1, 0.0).unwrap();
    let angle: f64 = 0.1 / PI;
    close(evolved[0].pos[0] / 5e153, angle.cos(), 4e-16);
    close(evolved[0].pos[1] / 5e153, angle.sin(), 4e-16);
}

#[test]
fn pair_time_refinement_agrees_with_a_fine_independent_rk4_reference() {
    let initial = [
        VortexParticle::new([0.4, -0.2], 0.7),
        VortexParticle::new([-0.5, 0.6], -0.2),
    ];
    let exact = advance_pair(initial, 0.6, 0.15).unwrap();
    let reference = fs_vpm::simulate(&initial, 0.0005, 1200, 0.15);
    for (a, b) in exact.iter().zip(reference) {
        for axis in 0..2 {
            close(a.pos[axis], b.pos[axis], 4e-14);
        }
    }
}

#[test]
fn fixed_leaves_validate_identity_and_fail_transactionally() {
    let pair = [
        VortexParticle::new([0.0, 0.0], 1.0),
        VortexParticle::new([0.0, 1.0], -1.0),
    ];
    assert_eq!(advance_pair(pair, 0.0, 0.2).unwrap(), pair);
    let coincident = [pair[0], pair[0]];
    assert_eq!(advance_pair(coincident, 0.5, 0.0).unwrap(), coincident);
    for bad in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        assert_eq!(
            advance_pair(pair, bad, 0.1),
            Err(PairFlowError::InvalidInput)
        );
        assert_eq!(
            advance_pair(pair, 0.0, bad),
            Err(PairFlowError::InvalidInput)
        );
        let mut invalid = pair;
        invalid[1].circulation = bad;
        assert_eq!(
            advance_pair(invalid, 0.0, 0.1),
            Err(PairFlowError::InvalidInput)
        );
    }
    assert_eq!(
        advance_pair(pair, 0.2, -1.0),
        Err(PairFlowError::InvalidInput)
    );
    let overflow = [
        VortexParticle::new([f64::MAX, 0.0], 1.0),
        VortexParticle::new([-f64::MAX, 0.0], 1.0),
    ];
    assert_eq!(
        advance_pair(overflow, 0.2, 0.0),
        Err(PairFlowError::Unrepresentable)
    );
    let angle_overflow = [
        VortexParticle::new([1e-200, 0.0], 1.0),
        VortexParticle::new([0.0, 0.0], 1.0),
    ];
    assert_eq!(
        advance_pair(angle_overflow, 1.0, 0.0),
        Err(PairFlowError::Unrepresentable)
    );
    assert_eq!(pair[0].pos, [0.0, 0.0]);
}
