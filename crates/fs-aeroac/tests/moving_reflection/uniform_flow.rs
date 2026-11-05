//! G1/G3 stationary wall images commute with a tangential Galilean transform.

use super::*;
use fs_aeroac::curle3d::moving::flow::UniformFlowFrame;

const WIND: [f64; 3] = [0.0, 0.0, 60.0];

fn frame() -> UniformFlowFrame {
    // The wall points must be translated as well as the two trajectories.
    UniformFlowFrame::new(WIND, 0.5, [100.0, -40.0, 30.0]).unwrap()
}

fn flow_path(
    frame: UniformFlowFrame,
    times: &[f64],
    motion: impl Fn(f64) -> SourceKinematics,
) -> PreparedKinematicTrajectory {
    PreparedKinematicTrajectory::in_uniform_flow(
        medium(),
        &knots(times, motion),
        frame,
        0.8,
        path_budget(),
        &mut || false,
    )
    .unwrap()
}

// Independent SOURCE images and laboratory-frame convected Green potentials.
// A wedge image is a polar rotation, never the production receiver transform.
fn potentials(x: [f64; 3], time: f64, index: usize, wedge: Option<usize>) -> [f64; 4] {
    let transform = |vector| match wedge {
        Some(order) => wedge_vector(vector, index, order),
        None => image(vector, index),
    };
    let motion = |time| {
        let value = source(time);
        if wedge.is_some() {
            wedge_motion(value)
        } else {
            value
        }
    };
    let displacement = |tau: f64| {
        let y = transform(motion(tau).position_m);
        core::array::from_fn::<_, 3, _>(|i| x[i] - y[i] - WIND[i] * (time - tau))
    };
    let length = |d: [f64; 3]| d.iter().map(|v| v * v).sum::<f64>().sqrt();
    let mut lo: f64 = -0.3;
    let mut hi: f64 = 0.4;
    for _ in 0..100 {
        let tau = lo.midpoint(hi);
        if tau + length(displacement(tau)) / 300.0 < time {
            lo = tau;
        } else {
            hi = tau;
        }
    }
    let tau = lo.midpoint(hi);
    let d = displacement(tau);
    let r = length(d);
    let velocity = transform(motion(tau).velocity_m_s);
    let kappa = 1.0
        - (0..3)
            .map(|i| d[i] * (velocity[i] - WIND[i]) / (r * 300.0))
            .sum::<f64>();
    let denominator = 4.0 * PI * r * kappa;
    let original = strength(tau);
    let force = transform(original.force_on_fluid_n);
    [
        original.volume_flow_m3_s / denominator,
        force[0] / denominator,
        force[1] / denominator,
        force[2] / denominator,
    ]
}

fn green_pressure(x: [f64; 3], time: f64, index: usize, wedge: Option<usize>, h: f64) -> f64 {
    // rho*(partial_t + U.grad) Phi_Q - div Phi_F. The microphone's
    // laboratory velocity must never appear in this Eulerian pressure law.
    let mut p = 1.2
        * derivative(
            |dt| {
                potentials(
                    core::array::from_fn(|i| x[i] + WIND[i] * dt),
                    time + dt,
                    index,
                    wedge,
                )[0]
            },
            h,
        );
    for axis in 0..3 {
        p -= derivative(
            |dx| {
                let mut point = x;
                point[axis] += dx;
                potentials(point, time, index, wedge)[axis + 1]
            },
            h,
        );
    }
    p
}

#[test]
fn accelerating_reflections_in_flow_match_independent_green_pressure_and_doppler() {
    for wedge in [None, Some(3)] {
        let motion = |time| {
            let value = source(time);
            if wedge.is_some() {
                wedge_motion(value)
            } else {
                value
            }
        };
        let microphone = |time| {
            let value = receiver(time);
            if wedge.is_some() {
                wedge_motion(value)
            } else {
                value
            }
        };
        let source_path = flow_path(frame(), &[-0.3, 0.0, 0.4], motion);
        let microphone_path = flow_path(frame(), &[0.0, 0.1, 0.2], microphone);
        for parity in [0, 3] {
            let observer = match wedge {
                Some(order) => MovingReflectionObserver::wedge_in_uniform_flow(
                    &microphone_path,
                    wedge_planes(order, parity),
                    order,
                    2 * order,
                    bounds(),
                    &mut || false,
                ),
                None => MovingReflectionObserver::in_uniform_flow(
                    &microphone_path,
                    &walls(parity, 2),
                    4,
                    bounds(),
                    &mut || false,
                ),
            }
            .unwrap();
            let owner = PreparedReflectedMovingReceiverSource::new(
                &source_path,
                &observer,
                config([0.03, 0.15]),
                bounds(),
                pairs(),
                &mut || false,
            )
            .unwrap();
            for time in [0.04, 0.08, 0.14] {
                let mut sum = 0.0;
                for index in 0..owner.path_count() {
                    let sample = owner.prepare_path(index, time, &mut || false).unwrap();
                    let parity_index = if wedge.is_some() {
                        index % 2
                    } else {
                        index.count_ones() as usize
                    };
                    let sign = if parity == 0 || parity_index % 2 == 0 {
                        1.0
                    } else {
                        -1.0
                    };
                    let x = microphone(time).position_m;
                    let coarse = sign * green_pressure(x, time, index, wedge, 2e-4);
                    let fine = sign * green_pressure(x, time, index, wedge, 1e-4);
                    near(fine, coarse, 2e-7, 0.01);
                    let actual = sample
                        .pressure
                        .evaluate(strength(sample.pressure.source_time_s()))
                        .unwrap()
                        .pressure_pa;
                    near(actual, fine, 2e-7, 0.01);
                    sum += fine;
                    let h = 1e-5;
                    let before = owner
                        .prepare_path(index, time - h, &mut || false)
                        .unwrap()
                        .pressure
                        .source_time_s();
                    let after = owner
                        .prepare_path(index, time + h, &mut || false)
                        .unwrap()
                        .pressure
                        .source_time_s();
                    near(
                        sample.observed_source_time_rate,
                        (after - before) / (2.0 * h),
                        3e-8,
                        1.0,
                    );
                }
                near(
                    owner
                        .pressure_at(time, &mut |tau| Ok(strength(tau)), &mut || false)
                        .unwrap(),
                    sum,
                    3e-7,
                    0.01,
                );
            }
        }
    }
}

#[test]
fn uniform_reflection_rejects_normal_flow_foreign_frames_and_incomplete_work() {
    let microphone = flow_path(frame(), &[0.0, 0.2], receiver);
    let source_path = flow_path(frame(), &[-0.3, 0.4], source);
    let planes = walls(0, 2);
    let observer =
        MovingReflectionObserver::in_uniform_flow(&microphone, &planes, 4, bounds(), &mut || false)
            .unwrap();
    let foreign = UniformFlowFrame::new(WIND, 0.0, [0.0; 3]).unwrap();
    let foreign_source = flow_path(foreign, &[-0.3, 0.4], source);
    assert!(
        PreparedReflectedMovingReceiverSource::new(
            &foreign_source,
            &observer,
            config([0.03, 0.15]),
            bounds(),
            pairs(),
            &mut || false,
        )
        .is_err()
    );
    // A flow tangential to the first two planes crosses the third; it cannot
    // be treated as a fixed fluid-frame wall even at a tiny nonzero speed.
    for wind in [WIND, [0.0, 0.0, 1e-100]] {
        let f = UniformFlowFrame::new(wind, 0.0, [0.0; 3]).unwrap();
        let mic = flow_path(f, &[0.0, 0.2], receiver);
        assert!(
            MovingReflectionObserver::in_uniform_flow(&mic, &walls(0, 3), 8, bounds(), &mut || {
                false
            },)
            .is_err()
        );
    }
    assert!(MovingReflectionObserver::in_uniform_flow(
        &microphone, &planes, 3, bounds(), &mut || false,
    ).is_err());
    assert!(
        MovingReflectionObserver::wedge_in_uniform_flow(
            &microphone,
            wedge_planes(3, 1),
            3,
            6,
            bounds(),
            &mut || false,
        )
        .is_err()
    );
    assert!(matches!(
        MovingReflectionObserver::in_uniform_flow(&microphone, &planes, 4, bounds(), &mut || true,),
        Err(AeroacError::Cancelled { .. })
    ));
    let one_path_budget = MovingReceiverBudget {
        maximum_segment_pairs: 1,
        ..pairs()
    };
    assert!(
        PreparedReflectedMovingReceiverSource::new(
            &source_path,
            &observer,
            config([0.03, 0.15]),
            bounds(),
            one_path_budget,
            &mut || false,
        )
        .is_err()
    );
}
