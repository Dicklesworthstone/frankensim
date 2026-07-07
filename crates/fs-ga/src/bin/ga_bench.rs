//! Transform-benchmark for the wqd.22 acceptance criterion: motors vs the
//! hand-written quaternion/matrix baselines. Emits JSON lines (P10).
//!
//! Doctrine measured here: COMPOSE in motor land (associative, drift-
//! correctable, no gimbal states), APPLY in matrix land (`Mat34` lowered
//! once per motor). The apply-side parity claim is for the lowered path;
//! the raw sandwich is also measured for honesty.

use fs_ga::{Mat34, Motor, Point, Quat, Vec3};
use std::hint::black_box;
use std::time::Instant;

fn lcg(seed: &mut u64) -> f64 {
    *seed = seed
        .wrapping_mul(6364136223846793005)
        .wrapping_add(1442695040888963407);
    ((*seed >> 11) as f64) / (1u64 << 53) as f64
}

fn time_ns<F: FnMut()>(mut f: F, reps: usize) -> f64 {
    let start = Instant::now();
    for _ in 0..reps {
        f();
    }
    start.elapsed().as_secs_f64() * 1e9 / reps as f64
}

#[allow(clippy::too_many_lines)]
fn main() {
    let mut seed = 0xBEBEu64;
    let n_points = 100_000usize;
    let points: Vec<Vec3> = (0..n_points)
        .map(|_| {
            Vec3::new(
                lcg(&mut seed) * 8.0 - 4.0,
                lcg(&mut seed) * 8.0 - 4.0,
                lcg(&mut seed) * 8.0 - 4.0,
            )
        })
        .collect();
    let axis = {
        let v = Vec3::new(0.3, -0.5, 0.81);
        v.scale(1.0 / v.norm())
    };
    let q = Quat::from_axis_angle(axis, 0.73);
    let t = Vec3::new(1.5, -2.0, 0.25);
    let motor = Motor::from_parts(q, t);
    let mat = Mat34::from_motor(&motor).expect("motor lowers");

    // Bulk point transforms.
    let quat_ns = time_ns(
        || {
            let mut acc = 0.0;
            for p in &points {
                let r = q.rotate(*p) + t;
                acc += r.x;
            }
            black_box(acc);
        },
        5,
    ) / n_points as f64;
    let mat_ns = time_ns(
        || {
            let mut acc = 0.0;
            for p in &points {
                acc += mat.apply(*p).x;
            }
            black_box(acc);
        },
        5,
    ) / n_points as f64;
    let sandwich_ns = time_ns(
        || {
            let mut acc = 0.0;
            for p in &points {
                let r = motor
                    .transform_point(Point::new(p.x, p.y, p.z))
                    .expect("Euclidean");
                acc += r.x;
            }
            black_box(acc);
        },
        5,
    ) / n_points as f64;

    // Composition chains (where motors do the real work).
    let incs: Vec<(Quat, Vec3, Motor)> = (0..4096)
        .map(|_| {
            let a = Vec3::new(
                lcg(&mut seed) - 0.5,
                lcg(&mut seed) - 0.5,
                lcg(&mut seed) - 0.5,
            );
            let a = a.scale(1.0 / a.norm());
            let qi = Quat::from_axis_angle(a, lcg(&mut seed) * 0.1);
            let ti = Vec3::new(lcg(&mut seed), lcg(&mut seed), lcg(&mut seed));
            (qi, ti, Motor::from_parts(qi, ti))
        })
        .collect();
    let quat_compose_ns = time_ns(
        || {
            let mut cq = Quat::identity();
            let mut ct = Vec3::new(0.0, 0.0, 0.0);
            for (qi, ti, _) in &incs {
                ct = qi.rotate(ct) + *ti;
                cq = *qi * cq;
            }
            black_box((cq.w, ct.x));
        },
        20,
    ) / incs.len() as f64;
    let motor_compose_ns = time_ns(
        || {
            let mut m = Motor::identity();
            for (_, _, mi) in &incs {
                m = mi.compose(&m);
            }
            black_box(m.0.0[0]);
        },
        20,
    ) / incs.len() as f64;

    println!(
        "{{\"suite\":\"fs-ga/bench\",\"metric\":\"transform_ns_per_point\",\
         \"quat\":{quat_ns:.2},\"mat34_lowered\":{mat_ns:.2},\"motor_sandwich\":{sandwich_ns:.2},\
         \"lowered_over_quat\":{:.3}}}",
        mat_ns / quat_ns
    );
    println!(
        "{{\"suite\":\"fs-ga/bench\",\"metric\":\"compose_ns_per_step\",\
         \"quat_plus_vec\":{quat_compose_ns:.2},\"motor\":{motor_compose_ns:.2},\
         \"motor_over_quat\":{:.3}}}",
        motor_compose_ns / quat_compose_ns
    );
}
