//! G1/G3/G4: accepted traction histories drive staged mechanics without replay drift.
use super::*;
use core::convert::Infallible;
use std::cell::Cell;

fn collect(
    law: &EllipticCoulomb,
    input: EllipticDrivenInput,
    step: f64,
) -> (EllipticDrivenReport, Vec<EllipticDrivenSegment>) {
    let mut segments = Vec::new();
    let report = law
        .advance_driven_segments(
            &law.initial_state(),
            &input,
            controls(step),
            &mut || false,
            &mut |segment| {
                segments.push(*segment);
                Ok::<(), Infallible>(())
            },
        )
        .unwrap();
    (report, segments)
}

#[test]
fn accepted_segments_reconstruct_force_history_motion_work_and_original_driven_report() {
    let law = EllipticCoulomb::new([0.8, 0.9], [0.2, 0.4], [1.0, 0.0]).unwrap();
    let i = input([1.1, 0.7], [[0.8, 0.25], [0.25, 1.2]], [0.4, -0.3], 0.2);
    let (report, segments) = collect(&law, i, 0.025);
    assert_eq!(report, apply(&law, i, 0.025));
    assert!(segments.len() > 1);
    let mut time = 0.0;
    let mut velocity = i.velocity_m_s;
    let mut displacement = [0.0; 2];
    let mut impulse = [0.0; 2];
    let mut heat = 0.0;
    let mut travel = 0.0;
    for segment in &segments {
        close(segment.elapsed_before_s, time, 2e-15);
        assert!(segment.dt_s > 0.0 && segment.dt_s <= 0.025);
        assert_eq!(segment.mode_before, EllipticDrivenMode::Sliding);
        assert_eq!(segment.mode_after, EllipticDrivenMode::Sliding);
        close2(segment.velocity_before_m_s, velocity, 8e-12);
        let recoil = mat(i.mobility_per_kg, segment.traction_n);
        let acceleration: [f64; 2] =
            core::array::from_fn(|j| i.driving_acceleration_m_s2[j] - recoil[j]);
        let h = segment.dt_s;
        // Independent integration of the actual held force, including its
        // first moment; no reported endpoint or displacement drives this state.
        let dx: [f64; 2] =
            core::array::from_fn(|j| h * velocity[j] + 0.5 * h * h * acceleration[j]);
        let next: [f64; 2] = core::array::from_fn(|j| velocity[j] + h * acceleration[j]);
        close2(segment.velocity_after_m_s, next, 8e-12);
        close2(segment.slip_increment_m, dx, 8e-13);
        assert!(dot(segment.traction_n, velocity) >= -1e-12);
        assert!(dot(segment.traction_n, next) >= -1e-12);
        let loss = dot(segment.traction_n, dx);
        close(segment.dissipated_j, loss, 8e-13);
        // Independent fixed-grid Simpson integral of unsigned linear-path speed.
        let n = 256;
        let speed = |k: usize| {
            norm(core::array::from_fn(|j| {
                velocity[j] + h * k as f64 / n as f64 * acceleration[j]
            }))
        };
        let mut integral = speed(0) + speed(n);
        for k in 1..n {
            integral += if k % 2 == 0 { 2.0 } else { 4.0 } * speed(k);
        }
        let path = integral * h / (3.0 * n as f64);
        close(segment.distance_m, path, 2e-12);
        for j in 0..2 {
            displacement[j] += dx[j];
            impulse[j] += h * segment.traction_n[j];
        }
        heat += loss;
        travel += path;
        velocity = next;
        time += h;
    }
    close(time, i.dt_s, 2e-15);
    close2(report.velocity_after_m_s, velocity, 8e-12);
    close2(report.slip_increment_m, displacement, 2e-12);
    close2(report.traction_impulse_ns, impulse, 2e-15);
    close(report.dissipated_j, heat, 3e-12);
    close(report.distance_m, travel, 3e-12);
    close(
        report.drive_work_j,
        dot(report.driving_force_n, displacement),
        3e-12,
    );
    let [[a, b], [_, d]] = i.mobility_per_kg;
    let kinetic = |v: [f64; 2]| {
        0.5 * (d * v[0] * v[0] - 2.0 * b * v[0] * v[1] + a * v[1] * v[1]) / (a * d - b * b)
    };
    close(
        kinetic(velocity) - kinetic(i.velocity_m_s) + heat - report.drive_work_j,
        0.0,
        1e-11,
    );
    let recoil = mat(i.mobility_per_kg, impulse);
    let averaged_path: [f64; 2] = core::array::from_fn(|j| {
        i.velocity_m_s[j] * i.dt_s
            + 0.5 * i.dt_s * (i.driving_acceleration_m_s2[j] * i.dt_s - recoil[j])
    });
    assert!(norm(core::array::from_fn(|j| averaged_path[j] - displacement[j])) > 1e-5);

    // Rejected false-capture/assisting proposals consume numerical budget but
    // never reach a pose visitor or appear as physical intervals.
    let isotropic = EllipticCoulomb::new([0.6; 2], [0.2; 2], [1.0, 0.0]).unwrap();
    let turning = input([-1.0, 0.1], [[1.0, 0.0], [0.0, 1.0]], [3.0, 0.0], 1.0);
    let (report, segments) = collect(&isotropic, turning, 0.2);
    assert!(report.substeps > segments.len());
    assert_eq!(report.first_arrest_time_s, None);
    for segment in segments {
        assert!(dot(segment.traction_n, segment.velocity_before_m_s) >= 0.0);
        assert_eq!(segment.mode_after, EllipticDrivenMode::Sliding);
    }
}

#[test]
fn static_stop_hold_and_reverse_visits_expose_their_physical_event_sides() {
    let law = EllipticCoulomb::new([0.6; 2], [0.2; 2], [1.0, 0.0]).unwrap();
    let w = [[0.75, 0.0], [0.0, 1.5]];
    let mut rest = input([1e-15, -1e-15], w, [2.25, 0.0], 0.4);
    rest.mode = EllipticDrivenMode::Resting;
    let (report, segments) = collect(&law, rest, 0.001);
    assert_eq!(segments.len(), 1);
    let hold = segments[0];
    assert_eq!(hold.mode_before, EllipticDrivenMode::Resting);
    assert_eq!(hold.mode_after, EllipticDrivenMode::Resting);
    assert_eq!(hold.velocity_before_m_s, [0.0; 2]);
    assert_eq!(hold.velocity_after_m_s, [0.0; 2]);
    assert_eq!(hold.dt_s, rest.dt_s);
    assert_eq!(hold.traction_n, [3.0, 0.0]);
    assert_eq!(hold.slip_increment_m, [0.0; 2]);
    assert_eq!(hold.distance_m, 0.0);
    assert_eq!(hold.dissipated_j, 0.0);
    assert_eq!(report, apply(&law, rest, 0.001));

    for q in [0.0, 3.0, 8.0] {
        let i = input([-2.0, 0.0], w, [0.75 * q, 0.0], 2.0);
        let (report, segments) = collect(&law, i, 0.001);
        assert_eq!(segments.len(), 2);
        let stop = 2.0 / (0.75 * (q + 2.0));
        let incoming = segments[0];
        assert_eq!(incoming.mode_before, EllipticDrivenMode::Sliding);
        assert_eq!(incoming.mode_after, EllipticDrivenMode::Resting);
        assert_eq!(incoming.velocity_after_m_s, [0.0; 2]);
        assert!(incoming.arrest_time_is_exact);
        close(incoming.dt_s, stop, 1e-15);
        close(report.first_arrest_time_s.unwrap(), stop, 1e-15);
        let outgoing = segments[1];
        assert_eq!(outgoing.mode_before, EllipticDrivenMode::Resting);
        assert_eq!(outgoing.velocity_before_m_s, [0.0; 2]);
        close(outgoing.elapsed_before_s, stop, 1e-15);
        close(outgoing.dt_s, 2.0 - stop, 1e-15);
        if q <= 6.0 {
            assert_eq!(outgoing.mode_after, EllipticDrivenMode::Resting);
            assert_eq!(outgoing.traction_n, [q, 0.0]);
            assert_eq!(outgoing.dissipated_j, 0.0);
        } else {
            assert_eq!(outgoing.mode_after, EllipticDrivenMode::Sliding);
            assert_eq!(outgoing.traction_n, [2.0, 0.0]);
            assert!(outgoing.velocity_after_m_s[0] > 0.0);
            assert!(outgoing.dissipated_j > 0.0);
        }
    }
}

#[test]
fn segmented_zero_drive_curvature_converges_without_changing_the_exact_aggregate_kernel() {
    let law = EllipticCoulomb::new([0.8, 0.9], [0.2, 0.4], [1.0, 0.0]).unwrap();
    let i = input([1.1, 0.7], [[0.8, 0.25], [0.25, 1.2]], [0.0; 2], 0.2);
    let oracle = reference(i.mobility_per_kg, [0.0; 2], i.velocity_m_s, i.dt_s, 8192);
    let exact = apply(&law, i, 0.05);
    close2(exact.velocity_after_m_s, [oracle[0], oracle[1]], 2e-9);
    close2(exact.slip_increment_m, [oracle[2], oracle[3]], 2e-9);
    close(exact.distance_m, oracle[4], 2e-9);
    close(exact.dissipated_j, oracle[5], 2e-9);
    let mut errors = Vec::new();
    for h in [0.05, 0.025, 0.0125] {
        let (report, segments) = collect(&law, i, h);
        assert!(segments.len() >= 4);
        assert!(segments.iter().all(|segment| segment.dt_s <= h));
        assert_ne!(report.velocity_after_m_s, exact.velocity_after_m_s);
        errors.push([
            norm([
                report.velocity_after_m_s[0] - oracle[0],
                report.velocity_after_m_s[1] - oracle[1],
            ]),
            norm([
                report.slip_increment_m[0] - oracle[2],
                report.slip_increment_m[1] - oracle[3],
            ]),
            (report.distance_m - oracle[4]).abs(),
            (report.dissipated_j - oracle[5]).abs(),
        ]);
        assert_eq!(report.drive_work_j, 0.0);
        close(report.energy_residual_j, 0.0, 2e-11);
        assert_eq!(report.first_arrest_time_s, None);
    }
    for j in 0..4 {
        assert!(
            errors[2][j] < 0.02,
            "segmented zero-drive errors {errors:?}"
        );
        for pair in errors.windows(2) {
            let ratio = pair[0][j] / pair[1][j];
            assert!(
                ratio > 1.4 && ratio < 2.6,
                "first-order zero-drive ladder {errors:?}"
            );
        }
    }
    assert_eq!(exact, apply(&law, i, 0.0125));
}

#[test]
fn callback_native_budget_and_cancellation_failures_leave_published_owners_unchanged() {
    let law = EllipticCoulomb::new([0.8, 0.9], [0.2, 0.4], [1.0, 0.0]).unwrap();
    let i = input([1.1, 0.7], [[0.8, 0.25], [0.25, 1.2]], [0.4, -0.3], 0.2);
    let original = (i.velocity_m_s, [0.0; 2], law.initial_state());
    let mut published = original.clone();
    for failure in 0..3 {
        let mut staged = published.clone();
        let cancelled = Cell::new(false);
        let mut visits = 0;
        let mut c = controls(0.025);
        if failure == 2 {
            c.maximum_substeps = 1;
        }
        let result = law.advance_driven_segments(
            &published.2,
            &i,
            c,
            &mut || cancelled.get(),
            &mut |segment| {
                visits += 1;
                staged.0 = segment.velocity_after_m_s;
                for j in 0..2 {
                    staged.1[j] += segment.slip_increment_m[j];
                }
                if failure == 0 && visits == 2 {
                    return Err("original pose owner refused");
                }
                if failure == 1 {
                    cancelled.set(true);
                }
                Ok(())
            },
        );
        let expected = match failure {
            0 => EllipticDrivenSegmentError::Visitor("original pose owner refused"),
            1 => EllipticDrivenSegmentError::Driven(EllipticDrivenError::Cancelled),
            _ => EllipticDrivenSegmentError::Driven(EllipticDrivenError::SubstepBudget),
        };
        assert_eq!(result.unwrap_err(), expected);
        assert_eq!(visits, if failure == 0 { 2 } else { 1 });
        assert_ne!(staged, published); // The discarded private candidate did move.
        assert_eq!(published, original);
    }
    let mut visits = 0;
    let cancelled =
        law.advance_driven_segments(&published.2, &i, controls(0.025), &mut || true, &mut |_| {
            visits += 1;
            Ok::<(), Infallible>(())
        });
    assert!(matches!(
        cancelled,
        Err(EllipticDrivenSegmentError::Driven(
            EllipticDrivenError::Cancelled
        ))
    ));
    assert_eq!(visits, 0);
    let mut staged = published.clone();
    let retry = law
        .advance_driven_segments(
            &published.2,
            &i,
            controls(0.025),
            &mut || false,
            &mut |segment| {
                staged.0 = segment.velocity_after_m_s;
                for j in 0..2 {
                    staged.1[j] += segment.slip_increment_m[j];
                }
                Ok::<(), Infallible>(())
            },
        )
        .unwrap();
    assert_eq!(retry, apply(&law, i, 0.025));
    staged.2 = retry.state;
    published = staged;
    assert_ne!(published, original);
    close2(published.1, retry.slip_increment_m, 2e-15);
}
