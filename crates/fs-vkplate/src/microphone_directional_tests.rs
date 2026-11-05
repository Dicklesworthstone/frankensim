//! Ideal pressure/velocity receivers on the actual mechanics clock.
//! Analytic fields are development fixtures, not measured capsule models.

use super::*;

fn source(instrument: &mut Instrument, position: [f64; 3]) {
    let mut body = Body::point("uniform small radiator", 1.0);
    body.radiates = true;
    body.h = 0.01;
    body.origin = position;
    instrument.add(body);
}

#[test]
fn g2_cardioid_proximity_vs_point_source() {
    let mut instrument = Instrument::new(48_000.0);
    source(&mut instrument, [0.0; 3]);
    let distance = 0.06;
    for alpha in [0.0, 0.5, 1.0] {
        for axis in [[0.0, 0.0, -1.0], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0]] {
            instrument
                .mics
                .push(Mic::first_order([0.0, 0.0, distance], alpha, axis).unwrap());
        }
    }
    instrument.prepare();
    let gain = instrument.gas.density * 0.01_f64.powi(2) / (2.0 * PI * distance);
    for mic in &instrument.mics {
        let acceleration = mic.taps.iter().map(|t| t.3 + t.4).sum::<f64>() * instrument.k().powi(2);
        let velocity = mic.vtaps.iter().map(|t| t.3 + t.4).sum::<f64>() * 2.0 * instrument.k();
        let (alpha, axis) = mic.first_order_pattern();
        // With exp(+i*w*t), v_r/p=(1-i/(kr))/(rho*c). The front axis
        // points toward the source, opposite outward radial velocity.
        for kr in [0.05, 0.5, 5.0, 50.0] {
            let omega = kr * instrument.gas.sound_speed / distance;
            let actual = [acceleration / gain, -velocity / (omega * gain)];
            let exact = [
                alpha - (1.0 - alpha) * axis[2],
                (1.0 - alpha) * axis[2] / kr,
            ];
            assert!(
                (actual[0] - exact[0]).hypot(actual[1] - exact[1])
                    < 2e-13 * exact[0].hypot(exact[1]).max(1.0)
            );
        }
    }
}

#[test]
fn g1_directional_pcm_converges_to_retarded_point_field_on_three_clocks() {
    let (frequency, distance, displacement) = (1000.0_f64, 0.06_f64, 1e-7_f64);
    let omega = 2.0 * PI * frequency;
    let mut previous = f64::INFINITY;
    for rate in [24_000.0, 48_000.0, 96_000.0] {
        let mut instrument = Instrument::new(rate);
        source(&mut instrument, [0.0; 3]);
        let phase_step = omega / rate;
        instrument.bodies[0].w[0] = displacement;
        instrument.bodies[0].wp[0] = displacement * phase_step.cos();
        instrument.drives.push(Drive {
            body: 0,
            node: 0,
            knots: vec![(0.0, 0.0)],
        });
        instrument
            .mics
            .push(Mic::first_order([0.0, 0.0, distance], 0.5, [0.0, 0.0, -1.0]).unwrap());
        let gain = instrument.gas.density * 0.01_f64.powi(2) / (2.0 * PI * distance);
        let flight = distance / instrument.gas.sound_speed;
        let force = -4.0 * (0.5 * phase_step).sin().powi(2) * rate * rate * displacement;
        let mut output = [0.0];
        let (mut residual, mut reference) = (0.0, 0.0);
        for n in 0..(12.0 * rate / frequency) as usize {
            instrument.drives[0].knots[0].1 = force * (n as f64 * phase_step).cos();
            instrument.step(&mut output);
            // Observation must leave the independently prescribed trajectory intact.
            let exact_motion = displacement * ((n + 1) as f64 * phase_step).cos();
            assert!((instrument.bodies[0].w[0] - exact_motion).abs() < 2e-10 * displacement);
            if n < (2.0 * rate / frequency) as usize {
                continue;
            }
            let phase = omega * (n as f64 / rate - flight);
            let exact = -gain
                * displacement
                * (omega * omega * phase.cos()
                    + 0.5 * instrument.gas.sound_speed / distance * omega * phase.sin());
            residual += (output[0] - exact).powi(2);
            reference += exact * exact;
        }
        let relative = (residual / reference).sqrt();
        assert!(
            relative < previous / 3.5,
            "rate={rate}, error={relative}, previous={previous}"
        );
        previous = relative;
    }
    assert!(previous < 5e-4);
}

#[test]
fn g2_floor_reflects_particle_velocity_direction_as_well_as_pressure() {
    let mut instrument = Instrument::new(48_000.0);
    let position = [0.02, 0.01, 0.1];
    source(&mut instrument, position);
    let receiver = [0.08, 0.04, 0.2];
    let (alpha, axis) = (0.3, [0.6, 0.8, 0.0]);
    let floor = floor_y(-0.1);
    let mut mic = Mic::first_order(receiver, alpha, axis).unwrap();
    mic.floor = Some(floor);
    instrument.mics.push(mic);
    instrument.prepare();
    // Independent image-SOURCE construction. Production prepares an image
    // RECEIVER and must reflect its axis as well as its position.
    let (mut expected_acceleration, mut expected_velocity) = (0.0, 0.0);
    for origin in [position, floor.image(position)] {
        let delta: [f64; 3] = std::array::from_fn(|i| receiver[i] - origin[i]);
        let radius = delta[0].hypot(delta[1]).hypot(delta[2]);
        let projection: f64 = (0..3).map(|i| axis[i] * delta[i] / radius).sum();
        let gain = instrument.gas.density * 0.01_f64.powi(2) / (2.0 * PI * radius);
        expected_acceleration += gain * (alpha - (1.0 - alpha) * projection);
        expected_velocity -=
            gain * (1.0 - alpha) * instrument.gas.sound_speed / radius * projection;
    }
    let mic = &instrument.mics[0];
    let acceleration = mic.taps.iter().map(|t| t.3 + t.4).sum::<f64>() * instrument.k().powi(2);
    let velocity = mic.vtaps.iter().map(|t| t.3 + t.4).sum::<f64>() * 2.0 * instrument.k();
    assert!((acceleration / expected_acceleration - 1.0).abs() < 2e-13);
    assert!((velocity / expected_velocity - 1.0).abs() < 2e-13);
}

#[test]
fn g0_pattern_admission_and_export_preserve_the_declared_observer() {
    for (alpha, axis) in [
        (f64::NAN, [1.0, 0.0, 0.0]),
        (-0.1, [1.0, 0.0, 0.0]),
        (1.1, [1.0, 0.0, 0.0]),
        (0.5, [0.0; 3]),
        (0.5, [0.0, 0.0, 2.0]),
    ] {
        assert!(Mic::first_order([0.0; 3], alpha, axis).is_err());
    }
    let mut instrument = Instrument::new(48_000.0);
    source(&mut instrument, [0.0; 3]);
    instrument.mics.push(Mic::new([0.0, 0.0, 0.1]));
    let original = instrument.image().mics[0].taps.clone();
    assert!(instrument.set_microphone_pattern(0, 0.5, [0.0; 3]).is_err());
    assert_eq!(instrument.image().mics[0].taps, original);
    instrument
        .set_microphone_pattern(0, 0.5, [0.0, 0.0, -1.0])
        .unwrap();
    assert!(!instrument.prepared);
    assert!(std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| instrument.image())).is_err());
    instrument.step(&mut [0.0]);
    assert!(
        instrument
            .set_microphone_pattern(0, 1.0, [0.0, 0.0, 1.0])
            .is_err()
    );
    assert_eq!(
        instrument.mics[0].first_order_pattern(),
        (0.5, [0.0, 0.0, -1.0])
    );
}
