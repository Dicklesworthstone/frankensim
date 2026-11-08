//! G1/G3: independent retarded spherical fields, coherent cancellation and
//! observation-only readout. No new acoustic source or fitted reference.
use super::*;
use std::f64::consts::TAU;

fn medium() -> RayleighMedium {
    RayleighMedium { density: 1.2, sound_speed: 320.0 }
}

fn patch(x: f64, shape: f64) -> SurfaceSample {
    SurfaceSample { position_m: [x, 0.0, 0.0], area_m2: 0.01, mode_shape: vec![shape] }
}

#[test]
fn complete_field_matches_retarded_spherical_pressure_and_velocity() {
    let rate = 48_000;
    let position = [0.12, -0.16, 0.15];
    let radius = 0.25;
    let direction = position.map(|x| x / radius);
    for ratio in [1, 3, 4] {
        for frequency in [120.0, 7000.0, 16_000.0] {
            let mut observer = BaffledField::new_multirate(
                &[patch(0.0, 1.0)], 1, rate * ratio, rate, position, medium(),
            ).unwrap();
            let delay = radius / medium().sound_speed
                + observer.filter_delay_samples() / f64::from(rate);
            let omega = TAU * frequency;
            let mut error = [0.0; 4];
            let mut signal = [0.0; 4];
            let mut trace = vec![0.0; ratio as usize];
            for frame in 0..2200 {
                for (sub, velocity) in trace.iter_mut().enumerate() {
                    let t = (frame * ratio as usize + sub) as f64 / f64::from(rate * ratio);
                    *velocity = (omega * t).sin();
                }
                let actual = observer.step_trace(&trace).unwrap();
                if frame < 500 { continue; }
                let t = ((frame + 1) * ratio as usize - 1) as f64 / f64::from(rate * ratio) - delay;
                let (flow, acceleration) = (0.01 * (omega * t).sin(), 0.01 * omega * (omega * t).cos());
                let p = medium().density * acceleration / (TAU * radius);
                let radial = (flow / (radius * radius) + acceleration / (medium().sound_speed * radius)) / TAU;
                let expected = [p, radial * direction[0], radial * direction[1], radial * direction[2]];
                let measured = [actual.pressure_pa, actual.particle_velocity_m_s[0],
                    actual.particle_velocity_m_s[1], actual.particle_velocity_m_s[2]];
                for c in 0..4 {
                    error[c] += (measured[c] - expected[c]).powi(2);
                    signal[c] += expected[c] * expected[c];
                }
            }
            for c in 0..4 {
                assert!((error[c] / signal[c]).sqrt() < 7e-4,
                    "ratio={ratio}, frequency={frequency}, channel={c}, relative RMS={}",
                    (error[c] / signal[c]).sqrt());
            }
        }
    }
}

#[test]
fn pressure_null_keeps_velocity_and_can_be_observed_from_new_directions() {
    let surface = [patch(-0.1, 1.0), patch(0.1, -1.0)];
    let position = [0.0, 0.0, 0.2];
    let mut observer = BaffledField::new(&surface, 1, 48_000, position, medium()).unwrap();
    let left = FirstOrder::new(0.0, [-1.0, 0.0, 0.0]).unwrap();
    let right = FirstOrder::new(0.0, [1.0, 0.0, 0.0]).unwrap();
    let mut reference = BaffledPressure::new_first_order(
        &surface, 1, 48_000, position, medium(), left,
    ).unwrap();
    let mut peak_velocity = 0.0_f64;
    for n in 0..1200 {
        let source = [(TAU * 180.0 * n as f64 / 48_000.0).sin()];
        let frame = observer.step(&source).unwrap();
        let expected = reference.step(&source).unwrap();
        assert_eq!(frame.pressure_pa, 0.0);
        assert_eq!(frame.particle_velocity_m_s[1], 0.0);
        assert_eq!(frame.particle_velocity_m_s[2], 0.0);
        peak_velocity = peak_velocity.max(frame.particle_velocity_m_s[0].abs());
        let measured = frame.microphone_pressure_pa(left).unwrap();
        assert!((measured - expected).abs() < 1e-12 * expected.abs().max(1.0));
        assert_eq!(measured, -frame.microphone_pressure_pa(right).unwrap());
        assert_eq!(frame.intensity_w_m2().unwrap(), [0.0; 3]);
        assert_eq!(observer.last_frame(), frame, "orientation queries must not advance history");
    }
    assert!(peak_velocity > 0.01, "a pressure null is not a silent vector field");
}

#[test]
fn intensity_retains_reactive_return_and_correct_cycle_average() {
    let rate = 48_000;
    let frequency = 120.0;
    let distance = 0.1;
    let area = 0.01;
    let mut observer = BaffledField::new(
        &[patch(0.0, 1.0)], 1, rate, [0.0, 0.0, distance], medium(),
    ).unwrap();
    let mut sum = 0.0;
    let mut negative = false;
    for n in 0..2400 {
        let frame = observer.step(&[(TAU * frequency * n as f64 / f64::from(rate)).sin()]).unwrap();
        if n >= 800 {
            let intensity = frame.intensity_w_m2().unwrap();
            assert_eq!([intensity[0], intensity[1]], [0.0; 2]);
            negative |= intensity[2] < 0.0;
            sum += intensity[2];
        }
    }
    // Four complete cycles: the reactive Q*Q' term has zero mean. The
    // remaining mean is rho/(2*c)*(A*omega/(2*pi*R))², independently derived.
    let expected = medium().density / (2.0 * medium().sound_speed)
        * (area * TAU * frequency / (TAU * distance)).powi(2);
    assert!(negative, "near-field stored energy must be allowed to return");
    assert!((sum / 1600.0 / expected - 1.0).abs() < 1e-4);
}

#[test]
fn density_changes_pressure_but_not_prescribed_flow_particle_velocity() {
    let mut a = BaffledField::new(&[patch(0.0, 1.0)], 1, 48_000, [0.1, 0.0, 0.3], medium()).unwrap();
    let dense = RayleighMedium { density: 3.0 * medium().density, ..medium() };
    let mut b = BaffledField::new(&[patch(0.0, 1.0)], 1, 48_000, [0.1, 0.0, 0.3], dense).unwrap();
    for n in 0..1000 {
        let source = [(0.03 * n as f64).sin()];
        let (x, y) = (a.step(&source).unwrap(), b.step(&source).unwrap());
        assert_eq!(y.medium(), dense);
        assert!((y.pressure_pa - 3.0 * x.pressure_pa).abs() < 1e-11);
        for c in 0..3 {
            assert!((y.particle_velocity_m_s[c] - x.particle_velocity_m_s[c]).abs() < 1e-13);
        }
    }
}

#[test]
fn refused_traces_preserve_all_fields_and_the_next_accepted_sample() {
    let mut observer = BaffledField::new_multirate(
        &[patch(0.0, 1.0)], 1, 144_000, 48_000, [0.1, 0.0, 0.3], medium(),
    ).unwrap();
    for n in 0..400 {
        observer.step_trace(&[n as f64 * 1e-5; 3]).unwrap();
    }
    let mut reference = observer.clone();
    let before = observer.last_frame();
    for bad in [&[0.0][..], &[0.0, 0.0, f64::NAN][..], &[f64::MAX; 3][..]] {
        assert!(observer.step_trace(bad).is_err());
        assert_eq!(observer.last_frame(), before);
    }
    for n in 0..400 {
        let trace = [0.01 * (n as f64 * 0.1).cos(); 3];
        assert_eq!(observer.step_trace(&trace).unwrap(), reference.step_trace(&trace).unwrap());
    }
}

#[test]
fn outlet_frontdoor_retains_geometry_clock_and_band_admission() {
    let receiver = CircularOutletReceiver {
        position_m: [0.0, 0.0, 0.2], radial_rings: 4, angular_points: 32,
        maximum_frequency_hz: 1000.0,
    };
    let mut field = BaffledField::circular_outlet(0.007, 48_000, 48_000, receiver, medium()).unwrap();
    let mut pressure = BaffledPressure::circular_outlet_multirate_differentiating(
        0.007,
        48_000,
        48_000,
        receiver,
        medium(),
    )
    .unwrap();
    for n in 0..300 {
        let flow = [1e-5 * (n as f64 * 0.03).sin()];
        assert_eq!(field.step(&flow).unwrap().pressure_pa, pressure.step(&flow).unwrap());
    }
    assert!(BaffledField::circular_outlet(0.007, 48_001, 48_000, receiver, medium()).is_err());
    let high = CircularOutletReceiver { maximum_frequency_hz: 22_000.0, ..receiver };
    assert!(BaffledField::circular_outlet(0.007, 48_000, 48_000, high, medium()).is_err());
    assert!(BaffledField::new(&[patch(0.0, 1.0)], 1, 48_000, [0.0; 3], medium()).is_err());
}
