//! The same causal owner observes accepted laboratory motion through flow.
use super::*;
use fs_aeroac::curle3d::moving::flow::UniformFlowFrame;

const WIND: [f64; 3] = [70.0, -20.0, 15.0];

fn frame() -> UniformFlowFrame {
    UniformFlowFrame::new(WIND, 0.0, [0.0; 3]).unwrap()
}

fn flow_path(
    frame: UniformFlowFrame,
    times: &[f64],
    motion: impl Fn(f64) -> SourceKinematics,
) -> PreparedKinematicTrajectory {
    PreparedKinematicTrajectory::in_uniform_flow(
        medium(),
        &samples(times, motion),
        frame,
        0.8,
        path_budget(),
        &mut || false,
    )
    .unwrap()
}

// Independent convected point Green potential in laboratory coordinates.
// Its retarded root uses the explicit polynomial, never a prepared path.
fn potentials_in_flow(x: [f64; 3], time: f64) -> [f64; 4] {
    let displacement = |tau: f64| {
        core::array::from_fn::<_, 3, _>(|i| {
            x[i] - source(tau).position_m[i] - WIND[i] * (time - tau)
        })
    };
    let length = |d: [f64; 3]| d.iter().map(|v| v * v).sum::<f64>().sqrt();
    let mut lo: f64 = -0.3;
    let mut hi = 0.4;
    for _ in 0..90 {
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
    let kappa = 1.0
        - (0..3)
            .map(|i| d[i] * (source(tau).velocity_m_s[i] - WIND[i]) / (r * 300.0))
            .sum::<f64>();
    let denominator = 4.0 * PI * r * kappa;
    let original = strength(tau);
    [
        original.volume_flow_m3_s / denominator,
        original.force_on_fluid_n[0] / denominator,
        original.force_on_fluid_n[1] / denominator,
        original.force_on_fluid_n[2] / denominator,
    ]
}

fn pressure_from_green_derivatives(x: [f64; 3], time: f64, step: f64) -> f64 {
    let derivative = |f: &dyn Fn(f64) -> f64| {
        (-f(2.0 * step) + 8.0 * f(step) - 8.0 * f(-step) + f(-2.0 * step)) / (12.0 * step)
    };
    // rho*(partial_t+U.grad) potential_Q - div potential_F.
    // The moving microphone's velocity MUST NOT enter this pressure derivative.
    let mut pressure = 1.2
        * derivative(&|h| {
            potentials_in_flow(core::array::from_fn(|i| x[i] + WIND[i] * h), time + h)[0]
        });
    for axis in 0..3 {
        pressure -= derivative(&|h| {
            let mut point = x;
            point[axis] += h;
            potentials_in_flow(point, time)[axis + 1]
        });
    }
    pressure
}

#[test]
fn accelerating_source_and_microphone_match_convected_green_pressure_and_clock() {
    let source_path = flow_path(frame(), &[-0.3, 0.0, 0.4], source);
    let microphone_path = flow_path(frame(), &[0.0, 0.1, 0.2], microphone);
    let owner = PreparedMovingReceiverSource::new(
        &source_path,
        &microphone_path,
        config([0.03, 0.15], 0.1),
        observation_budget(),
        &mut || false,
    )
    .unwrap();
    for time in [0.04, 0.08, 0.14] {
        let observed = owner.prepare_receiver(time, &mut || false).unwrap();
        let pressure = observed
            .pressure
            .evaluate(strength(observed.pressure.source_time_s()))
            .unwrap()
            .pressure_pa;
        let x = microphone(time).position_m;
        let coarse = pressure_from_green_derivatives(x, time, 2e-4);
        let fine = pressure_from_green_derivatives(x, time, 1e-4);
        near(fine, coarse, 1e-7, 0.01);
        near(pressure, fine, 1e-7, 0.01);
        let h = 1e-5;
        let before = owner
            .prepare_receiver(time - h, &mut || false)
            .unwrap()
            .pressure
            .source_time_s();
        let after = owner
            .prepare_receiver(time + h, &mut || false)
            .unwrap()
            .pressure
            .source_time_s();
        near(
            observed.observed_source_time_rate,
            (after - before) / (2.0 * h),
            2e-9,
            1.0,
        );
        let tau = observed.pressure.source_time_s();
        let retarded_displacement: [f64; 3] =
            core::array::from_fn(|i| x[i] - source(tau).position_m[i] - WIND[i] * (time - tau));
        let r = retarded_displacement
            .iter()
            .map(|v| v * v)
            .sum::<f64>()
            .sqrt();
        assert!((time - tau - r / 300.0).abs() <= 1.1e-13);
        let laboratory = microphone_path.laboratory_kinematics_at(time).unwrap();
        for axis in 0..3 {
            near(laboratory.position_m[axis], x[axis], 1e-13, 1.0);
            near(
                laboratory.velocity_m_s[axis],
                microphone(time).velocity_m_s[axis],
                1e-13,
                1.0,
            );
            near(
                laboratory.acceleration_m_s2[axis],
                microphone(time).acceleration_m_s2[axis],
                1e-12,
                1.0,
            );
        }
    }
}

#[test]
fn zero_flow_keeps_original_quintics_pressure_and_rate_bitwise() {
    let still = UniformFlowFrame::new([0.0; 3], 0.0, [0.0; 3]).unwrap();
    let source_before = path(&[-0.3, 0.0, 0.4], source);
    let receiver_before = path(&[0.0, 0.1, 0.2], microphone);
    let source_after = flow_path(still, &[-0.3, 0.0, 0.4], source);
    let receiver_after = flow_path(still, &[0.0, 0.1, 0.2], microphone);
    let before = PreparedMovingReceiverSource::new(
        &source_before,
        &receiver_before,
        config([0.03, 0.15], 0.1),
        observation_budget(),
        &mut || false,
    )
    .unwrap();
    let after = PreparedMovingReceiverSource::new(
        &source_after,
        &receiver_after,
        config([0.03, 0.15], 0.1),
        observation_budget(),
        &mut || false,
    )
    .unwrap();
    for time in [0.04, 0.08, 0.14] {
        assert_eq!(
            source_before.kinematics_at(time).unwrap(),
            source_after.kinematics_at(time).unwrap()
        );
        let a = before.prepare_receiver(time, &mut || false).unwrap();
        let b = after.prepare_receiver(time, &mut || false).unwrap();
        assert_eq!(
            a.pressure.source_time_s().to_bits(),
            b.pressure.source_time_s().to_bits()
        );
        assert_eq!(
            a.observed_source_time_rate.to_bits(),
            b.observed_source_time_rate.to_bits()
        );
        let strength = strength(a.pressure.source_time_s());
        assert_eq!(
            a.pressure.evaluate(strength).unwrap().pressure_pa.to_bits(),
            b.pressure.evaluate(strength).unwrap().pressure_pa.to_bits()
        );
    }
}

#[test]
fn large_reference_events_and_fluid_relative_speed_remain_physical() {
    let epoch = 1e12;
    let origin = [1e12, -1e12, 0.0];
    let frame = UniformFlowFrame::new([64.0, -16.0, 8.0], epoch, origin).unwrap();
    let event = SourceKinematics {
        position_m: [origin[0] + 8.0, origin[1] - 2.0, 1.0],
        velocity_m_s: [64.0, -16.0, 8.0],
        acceleration_m_s2: [1.0, -2.0, 3.0],
    };
    let fluid = frame.to_fluid_frame(epoch + 0.125, event).unwrap();
    assert_eq!(fluid.position_m, [0.0; 3]);
    assert_eq!(fluid.velocity_m_s, [0.0; 3]);
    assert_eq!(fluid.acceleration_m_s2, event.acceleration_m_s2);
    assert_eq!(
        frame.to_laboratory_frame(epoch + 0.125, fluid).unwrap(),
        event
    );

    let samples = samples(&[-0.3, 0.4], |t| SourceKinematics {
        position_m: [270.0 * t, 0.0, 0.0],
        velocity_m_s: [270.0, 0.0, 0.0],
        ..SourceKinematics::default()
    });
    assert!(
        PreparedKinematicTrajectory::new(medium(), &samples, 0.2, path_budget(), &mut || false)
            .is_err()
    );
    let frame = UniformFlowFrame::new([240.0, 0.0, 0.0], 0.0, [0.0; 3]).unwrap();
    let path = PreparedKinematicTrajectory::in_uniform_flow(
        medium(),
        &samples,
        frame,
        0.2,
        path_budget(),
        &mut || false,
    )
    .unwrap();
    near(
        path.kinematics_at(0.1).unwrap().velocity_m_s[0],
        30.0,
        1e-13,
        1.0,
    );
    near(
        path.laboratory_kinematics_at(0.1).unwrap().velocity_m_s[0],
        270.0,
        1e-13,
        1.0,
    );
    assert!(
        PreparedKinematicTrajectory::in_uniform_flow(
            medium(),
            &samples,
            frame,
            0.05,
            path_budget(),
            &mut || false
        )
        .is_err()
    );
}

#[test]
fn frames_domains_storage_and_cancellation_cannot_be_mixed_or_bypassed() {
    let source_path = flow_path(frame(), &[-0.3, 0.0, 0.4], source);
    for other in [
        UniformFlowFrame::new([69.0, -20.0, 15.0], 0.0, [0.0; 3]).unwrap(),
        UniformFlowFrame::new(WIND, 1.0, [0.0; 3]).unwrap(),
        UniformFlowFrame::new(WIND, 0.0, [1.0, 0.0, 0.0]).unwrap(),
    ] {
        let receiver = flow_path(other, &[0.0, 0.1, 0.2], microphone);
        assert!(
            PreparedMovingReceiverSource::new(
                &source_path,
                &receiver,
                config([0.03, 0.15], 0.1),
                observation_budget(),
                &mut || false
            )
            .is_err()
        );
    }
    let receiver = path(&[0.0, 0.1, 0.2], microphone);
    assert!(
        PreparedMovingReceiverSource::new(
            &source_path,
            &receiver,
            config([0.03, 0.15], 0.1),
            observation_budget(),
            &mut || false
        )
        .is_err()
    );
    assert!(source_path.laboratory_kinematics_at(0.41).is_err());
    assert!(UniformFlowFrame::new(WIND, f64::NAN, [0.0; 3]).is_err());
    let samples = samples(&[-0.3, 0.0, 0.4], source);
    let mut budget = path_budget();
    budget.maximum_payload_bytes = 1;
    assert!(
        PreparedKinematicTrajectory::in_uniform_flow(
            medium(),
            &samples,
            frame(),
            0.8,
            budget,
            &mut || false
        )
        .is_err()
    );
    assert!(matches!(
        PreparedKinematicTrajectory::in_uniform_flow(
            medium(),
            &samples,
            frame(),
            0.8,
            path_budget(),
            &mut || true
        ),
        Err(AeroacError::Cancelled { .. })
    ));
}
