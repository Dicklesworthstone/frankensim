//! Source-physics checks for prepared velocity/Jacobian tree queries.
use fs_vpm::filament3d::{
    PreparedWakeTree, WakeKinematics, WakeSegment, WakeTreeAccuracy, segment_kinematics_with_core,
};

const WORK: u64 = 10_000_000;

fn prepare(sources: &[WakeSegment], leaf_capacity: usize) -> PreparedWakeTree {
    let bytes = PreparedWakeTree::required_bytes(sources.len()).unwrap();
    let (tree, report) =
        PreparedWakeTree::prepare(sources, leaf_capacity, bytes, WORK, &mut || false).unwrap();
    assert!(report.work_used > 0);
    assert_eq!(tree.source_count(), sources.len());
    tree
}

fn accuracy(velocity: f64, gradient: f64) -> WakeTreeAccuracy {
    WakeTreeAccuracy {
        absolute_velocity_error_mps: velocity,
        absolute_gradient_error_per_s: gradient,
    }
}

fn sources(n: usize, scale: f64) -> Vec<WakeSegment> {
    (0..n)
        .map(|i| {
            let phase = 2.0 * core::f64::consts::PI * i as f64 / n as f64;
            let z = 0.8 * i as f64 / n as f64 - 0.4;
            WakeSegment {
                start: [scale * phase.cos(), scale * phase.sin(), scale * z],
                end: [
                    scale * (phase + 0.08).cos(),
                    scale * (phase + 0.08).sin(),
                    scale * (z + 0.03),
                ],
                circulation: if i % 3 == 0 { -0.7 } else { 1.0 },
                core2_m2: scale * scale * (0.002 + 0.0005 * (i % 5) as f64),
            }
        })
        .collect()
}

// Original-source-order reference, through the separately tested exact
// segment owner. Compensated summation fixes the intended zero-tolerance order.
fn direct(sources: &[WakeSegment], p: [f64; 3]) -> WakeKinematics {
    let mut out = WakeKinematics::default();
    let mut correction = WakeKinematics::default();
    for s in sources {
        if s.circulation == 0.0 {
            continue;
        }
        let k = segment_kinematics_with_core(p, s.start, s.end, s.core2_m2).unwrap();
        for i in 0..3 {
            add(&mut out.velocity[i], &mut correction.velocity[i], s.circulation * k.velocity[i]);
            for j in 0..3 {
                add(
                    &mut out.gradient[i][j],
                    &mut correction.gradient[i][j],
                    s.circulation * k.gradient[i][j],
                );
            }
        }
    }
    out
}

fn add(sum: &mut f64, correction: &mut f64, value: f64) {
    let increment = value - *correction;
    let next = *sum + increment;
    *correction = (next - *sum) - increment;
    *sum = next;
}

fn velocity_error(a: &WakeKinematics, b: &WakeKinematics) -> f64 {
    (a.velocity[0] - b.velocity[0])
        .hypot(a.velocity[1] - b.velocity[1])
        .hypot(a.velocity[2] - b.velocity[2])
}

fn gradient_error(a: &WakeKinematics, b: &WakeKinematics) -> f64 {
    a.gradient
        .iter()
        .flatten()
        .zip(b.gradient.iter().flatten())
        .fold(0.0_f64, |n, (a, b)| n.hypot(a - b))
}

#[test]
fn zero_tolerance_preserves_original_order_and_snapshot() {
    let mut input = sources(96, 1.0);
    input.reverse();
    input.push(WakeSegment {
        start: [0.0; 3],
        end: [0.0, 0.0, 1.0],
        circulation: 0.0,
        core2_m2: 0.0,
    });
    let before = input.clone();
    let mut tree = prepare(&input, 3);
    let p = [2.1, -0.7, 0.3];
    let reference = direct(&input, p);
    for settings in [accuracy(0.0, 0.0), accuracy(1.0, 0.0), accuracy(0.0, 1.0)] {
        let (value, report) = tree.kinematics(p, settings, WORK, &mut || false).unwrap();
        assert_eq!(value.unwrap(), reference);
        assert_eq!(report.accepted_cells, 0);
        assert_eq!(report.visited_cells, 0);
        assert_eq!(report.velocity_error_bound_mps, 0.0);
        assert_eq!(report.gradient_error_bound_per_s, 0.0);
    }
    assert_eq!(input, before);
}

#[test]
fn mixed_core_expansions_bound_both_fields_across_physical_scales() {
    for scale in [1e-3, 1.0, 1e3] {
        let input = sources(96, scale);
        let mut tree = prepare(&input, 4);
        let p = [18.0 * scale, 5.0 * scale, 7.0 * scale];
        let reference = direct(&input, p);
        let mut coarse_direct = 0;
        for (index, settings) in [
            accuracy(1e-3 / scale, 1e-4 / (scale * scale)),
            accuracy(1e-9 / scale, 1e-10 / (scale * scale)),
        ]
        .into_iter()
        .enumerate()
        {
            let (value, report) = tree.kinematics(p, settings, WORK, &mut || false).unwrap();
            let value = value.unwrap();
            assert!(
                velocity_error(&value, &reference)
                    <= report.velocity_error_bound_mps + 1e-12 / scale
            );
            assert!(
                gradient_error(&value, &reference)
                    <= report.gradient_error_bound_per_s + 1e-12 / (scale * scale)
            );
            assert!(report.velocity_error_bound_mps <= settings.absolute_velocity_error_mps);
            assert!(report.gradient_error_bound_per_s <= settings.absolute_gradient_error_per_s);
            if index == 0 {
                assert!(report.accepted_cells > 0);
                assert!(report.direct_interactions < input.len() as u64);
                coarse_direct = report.direct_interactions;
            } else {
                assert!(report.direct_interactions >= coarse_direct);
            }
        }
    }
}

#[test]
fn analytic_tree_gradient_agrees_with_independent_velocity_differences() {
    let input = sources(96, 1.0);
    let mut tree = prepare(&input, 4);
    let p = [18.0, 5.0, 7.0];
    let (value, report) =
        tree.kinematics(p, accuracy(1e-3, 1e-4), WORK, &mut || false).unwrap();
    let value = value.unwrap();
    assert!(report.accepted_cells > 0);
    let h = 1e-3;
    let mut fd = WakeKinematics::default();
    for j in 0..3 {
        let mut plus = p;
        let mut minus = p;
        plus[j] += h;
        minus[j] -= h;
        let (vp, _) = tree.velocity(plus, 0.0, WORK, &mut || false).unwrap();
        let (vm, _) = tree.velocity(minus, 0.0, WORK, &mut || false).unwrap();
        let vp = vp.unwrap();
        let vm = vm.unwrap();
        for i in 0..3 {
            fd.gradient[i][j] = (vp[i] - vm[i]) / (2.0 * h);
        }
    }
    assert!(gradient_error(&value, &fd) <= report.gradient_error_bound_per_s + 1e-11);
    assert!(value.divergence().abs() < 1e-12);
}

fn direct_core(input: &[WakeSegment], p: [f64; 3], factor: f64) -> bool {
    input.iter().any(|s| {
        if s.circulation == 0.0 {
            return false;
        }
        let ab: [f64; 3] = core::array::from_fn(|i| s.end[i] - s.start[i]);
        let ap: [f64; 3] = core::array::from_fn(|i| p[i] - s.start[i]);
        let l2 = ab[0] * ab[0] + ab[1] * ab[1] + ab[2] * ab[2];
        let t = if l2 == 0.0 {
            0.0
        } else {
            ((ap[0] * ab[0] + ap[1] * ab[1] + ap[2] * ab[2]) / l2).clamp(0.0, 1.0)
        };
        let d: [f64; 3] = core::array::from_fn(|i| p[i] - (s.start[i] + t * ab[i]));
        fs_math::det::sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2])
            < factor * fs_math::det::sqrt(s.core2_m2)
    })
}

#[test]
fn exact_core_lookup_culls_far_nodes_and_keeps_native_mask_semantics() {
    let mut input = sources(96, 1.0);
    input.push(WakeSegment {
        start: [-5.0, 2.0, 3.0],
        end: [-5.0, 2.0, 3.0],
        circulation: 1.0,
        core2_m2: 0.04,
    });
    input.push(WakeSegment {
        start: [-10.0, 0.0, 0.0],
        end: [-10.0, 0.0, 1.0],
        circulation: 0.0,
        core2_m2: 1000.0,
    });
    let mut tree = prepare(&input, 2);
    for p in [
        [-5.0, 2.0, 3.0],
        [-10.0, 0.0, 0.5],
        [1.0, 0.0, -0.4],
        [0.0, 0.0, 0.0],
        [0.7, 0.7, 0.3],
        [3.0, 4.0, 1.0],
    ] {
        let (inside, _) = tree.in_core(p, 4.0, WORK, &mut || false).unwrap();
        assert_eq!(inside, direct_core(&input, p, 4.0));
    }
    let (inside, report) = tree.in_core([1e4, 0.0, 0.0], 4.0, WORK, &mut || false).unwrap();
    assert!(!inside);
    assert_eq!(report.direct_interactions, 0);
    assert_eq!(report.visited_cells, 1);

    let line = [WakeSegment {
        start: [-1.0, 0.0, 0.0],
        end: [1.0, 0.0, 0.0],
        circulation: 1.0,
        core2_m2: 0.0625,
    }];
    let mut line = prepare(&line, 1);
    assert!(!line.in_core([0.0, 1.0, 0.0], 4.0, WORK, &mut || false).unwrap().0);
    assert!(line.in_core([0.0, 0.999, 0.0], 4.0, WORK, &mut || false).unwrap().0);
}

#[test]
fn zero_core_singular_probes_are_masked_without_hiding_regular_axis_probes() {
    let input = [WakeSegment {
        start: [0.0, 0.0, -1.0],
        end: [0.0, 0.0, 1.0],
        circulation: 1.0,
        core2_m2: 0.0,
    }];
    let mut tree = prepare(&input, 1);
    for tolerance in [0.0, 1e-3] {
        assert!(tree.velocity([0.0; 3], tolerance, WORK, &mut || false).unwrap().0.is_none());
        assert!(
            tree.kinematics([0.0; 3], accuracy(tolerance, tolerance), WORK, &mut || false)
                .unwrap().0.is_none()
        );
        assert!(
            tree.kinematics([0.0, 0.0, 2.0], accuracy(tolerance, tolerance), WORK, &mut || false)
                .unwrap().0.is_some()
        );
    }
    let mut absent = input;
    absent[0].circulation = 0.0;
    let mut tree = prepare(&absent, 1);
    let (value, _) = tree.kinematics([0.0; 3], accuracy(1.0, 1.0), WORK, &mut || false).unwrap();
    assert_eq!(value.unwrap(), WakeKinematics::default());
    assert!(!tree.in_core([0.0; 3], 4.0, WORK, &mut || false).unwrap().0);
}

#[test]
fn underflowed_weights_and_moments_fall_back_to_representable_exact_fields() {
    let tiny = WakeSegment {
        start: [-5e-21, 0.0, 0.0],
        end: [5e-21, 0.0, 0.0],
        circulation: 1e-310,
        core2_m2: 1e-40,
    };
    // Root unsigned weight is zero, but the physical field is not.
    let mut tree = prepare(&[tiny, tiny], 1);
    let p = [0.0, 1e-20, 0.0];
    let (value, report) =
        tree.kinematics(p, accuracy(1e-300, 1e-280), WORK, &mut || false).unwrap();
    let value = value.unwrap();
    assert!(value.velocity[2] > 5e-292);
    assert_eq!(value, direct(&[tiny, tiny], p));
    assert_eq!(report.accepted_cells, 0);

    // A zero-weight child inside a positive-weight tree must also open.
    // Both large sources are axial at p and contribute exactly zero velocity.
    let large = WakeSegment {
        start: [0.0, 10.0, 0.0],
        end: [0.0, 11.0, 0.0],
        circulation: 1.0,
        core2_m2: 0.01,
    };
    let mixed = [tiny, tiny, large, large];
    let mut tree = prepare(&mixed, 1);
    let (value, report) = tree.velocity(p, 1e-300, WORK, &mut || false).unwrap();
    assert!(value.unwrap()[2] > 5e-292);
    assert_eq!(report.accepted_cells, 0);

    // Here W is normal, while W*a and the continuous-line quadrupole
    // underflow. An apparently small Taylor remainder cannot certify them.
    let normal_weight = WakeSegment { circulation: 1e-280, ..tiny };
    let input = [normal_weight, normal_weight];
    let p = [0.0, 1e-18, 0.0];
    let mut tree = prepare(&input, 1);
    let (value, report) =
        tree.kinematics(p, accuracy(1e-270, 2e-252), WORK, &mut || false).unwrap();
    assert_eq!(value.unwrap(), direct(&input, p));
    assert_eq!(report.accepted_cells, 0);
}

#[test]
fn nonnormal_derivative_scales_open_nodes_instead_of_claiming_zero_error() {
    for scale in [1e-40, 1e40] {
        let input = sources(16, scale);
        let p = [18.0 * scale, 5.0 * scale, 7.0 * scale];
        let mut tree = prepare(&input, 1);
        let (value, report) =
            tree.kinematics(p, accuracy(1e100, 1e200), WORK, &mut || false).unwrap();
        let value = value.unwrap();
        assert!(value.velocity.iter().any(|v| *v != 0.0));
        assert!(value.gradient.iter().flatten().all(|v| v.is_finite()));
        assert_eq!(report.accepted_cells, 0);
        assert_eq!(report.direct_interactions, input.len() as u64);
    }
}

#[test]
fn memory_work_admission_and_final_cancellation_preserve_reusable_snapshot() {
    let input = sources(64, 1.0);
    let bytes = PreparedWakeTree::required_bytes(input.len()).unwrap();
    assert!(PreparedWakeTree::required_bytes(usize::MAX).is_none());
    assert!(PreparedWakeTree::prepare(&input, 4, bytes - 1, WORK, &mut || false).is_err());
    assert!(PreparedWakeTree::prepare(&input, 0, bytes, WORK, &mut || false).is_err());
    assert!(PreparedWakeTree::prepare(&input, 4, bytes, 0, &mut || false).is_err());
    let mut tree = prepare(&input, 4);
    let p = [3.0, -0.5, 0.4];
    let settings = accuracy(0.0, 0.0);
    let reference = tree.kinematics(p, settings, WORK, &mut || false).unwrap();
    assert_eq!(
        tree.kinematics(p, settings, reference.1.work_used, &mut || false).unwrap(),
        reference
    );
    assert!(
        tree.kinematics(p, settings, reference.1.work_used - 1, &mut || false).is_err()
    );
    // Entry check, first exact batch, then final publication check.
    let mut checkpoints = 0;
    assert!(
        tree.kinematics(p, settings, WORK, &mut || {
            checkpoints += 1;
            checkpoints == 3
        }).is_err()
    );
    assert_eq!(tree.kinematics(p, settings, WORK, &mut || false).unwrap(), reference);
    let mut checkpoints = 0;
    assert!(
        tree.in_core([1e4, 0.0, 0.0], 4.0, WORK, &mut || {
            checkpoints += 1;
            checkpoints == 3
        }).is_err()
    );
    assert!(tree.velocity(p, f64::NAN, WORK, &mut || false).is_err());
    assert!(tree.kinematics([f64::NAN, 0.0, 0.0], settings, WORK, &mut || false).is_err());
    let mut malformed = input;
    malformed[0].core2_m2 = -1.0;
    assert!(PreparedWakeTree::prepare(&malformed, 4, bytes, WORK, &mut || false).is_err());
}
