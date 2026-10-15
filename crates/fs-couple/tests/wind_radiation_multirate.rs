//! Circular wind ports use the same broadband observer as structural sources.
use fs_couple::pcm_wav::baffled::{BaffledPressure, CircularOutletReceiver, RayleighMedium};
use fs_couple::pcm_wav::baffled::openings::CircularOpening;
use std::f64::consts::{PI, TAU};

const OUTPUT: u32 = 48_000;
fn medium() -> RayleighMedium { RayleighMedium { density: 1.2, sound_speed: 320.0 } }
fn receiver() -> CircularOutletReceiver {
    CircularOutletReceiver { position_m: [0.0, 0.0, 0.2], radial_rings: 64,
        angular_points: 64, maximum_frequency_hz: 21_600.0 }
}

#[test]
fn disk_pressure_matches_the_continuous_integral_through_upper_audio_partials() {
    let a = 0.006_f64;
    let area = PI * a * a;
    let cfg = receiver();
    let gas = medium();
    let near = cfg.position_m[2] / gas.sound_speed;
    let far = (cfg.position_m[2].powi(2) + a * a).sqrt() / gas.sound_speed;
    for ratio in [1, 4] {
        for frequency in [440.0, 8000.0, 18000.0, 21600.0] {
            let mechanics = OUTPUT * ratio;
            let mut mic = BaffledPressure::circular_outlet_multirate_differentiating(
                a, mechanics, OUTPUT, cfg, gas).unwrap();
            assert_eq!(mic.filter_delay_samples(), if ratio == 1 { 80.0 } else { 84.0 });
            let latency = mic.filter_delay_samples() / f64::from(OUTPUT);
            let omega = TAU * frequency;
            let amplitude = 1e-5;
            let mut trace = vec![0.0; ratio as usize];
            let (mut error, mut signal) = (0.0, 0.0);
            for frame in 0..1800 {
                for (sub, flow) in trace.iter_mut().enumerate() {
                    let t = (frame * ratio as usize + sub + 1) as f64 / f64::from(mechanics);
                    *flow = amplitude * (omega * t).sin();
                }
                let pressure = mic.step_trace(&trace).unwrap();
                if frame > 500 {
                    let t = (frame + 1) as f64 / f64::from(OUTPUT) - latency;
                    // Exact axial Rayleigh integral: no numerical surface rule,
                    // finite difference or discrete derivative in this reference.
                    let exact = gas.density * gas.sound_speed / area * amplitude
                        * ((omega * (t - near)).sin() - (omega * (t - far)).sin());
                    error += (pressure - exact).powi(2);
                    signal += exact.powi(2);
                }
            }
            assert!((error / signal).sqrt() < 1e-3,
                "ratio={ratio}, f={frequency}, RMS={}", (error / signal).sqrt());
        }
    }
}

#[test]
fn fused_multirate_outlets_keep_every_signed_input_and_the_original_filter_latency() {
    let openings = [CircularOpening { center_m: [0.0, 0.0], radius_m: 0.006 },
        CircularOpening { center_m: [0.1, 0.0], radius_m: 0.004 }];
    let cfg = receiver();
    for ratio in [1, 3, 4, 16] {
        let mechanics = OUTPUT * ratio;
        let mut combined = BaffledPressure::circular_openings_multirate_differentiating(
            &openings, mechanics, OUTPUT, cfg, medium()).unwrap();
        let mut singles: Vec<_> = openings.iter().map(|disk| {
            let local = CircularOutletReceiver { position_m: [cfg.position_m[0] - disk.center_m[0],
                cfg.position_m[1] - disk.center_m[1], cfg.position_m[2]], ..cfg };
            BaffledPressure::circular_outlet_multirate_differentiating(
                disk.radius_m, mechanics, OUTPUT, local, medium()).unwrap()
        }).collect();
        assert_eq!(combined.filter_delay_samples(), singles[0].filter_delay_samples());
        let mut trace = vec![0.0; 2 * ratio as usize];
        let mut one = vec![0.0; ratio as usize];
        let mut copied = None;
        for frame in 0..900 {
            for (sub, q) in trace.chunks_exact_mut(2).enumerate() {
                let t = (frame * ratio as usize + sub + 1) as f64 / f64::from(mechanics);
                q[0] = 1e-5 * (TAU * 18000.0 * t).sin();
                q[1] = -2e-5 * (TAU * 7333.0 * t).sin();
            }
            if frame == 450 { copied = Some(combined.clone()); }
            let pressure = combined.step_trace(&trace).unwrap();
            if let Some(copy) = &mut copied {
                assert_eq!(pressure.to_bits(), copy.step_trace(&trace).unwrap().to_bits());
            }
            let mut components = [0.0; 2];
            combined.mode_pressures(&mut components).unwrap();
            let mut expected = 0.0;
            for i in 0..2 {
                for (sub, q) in one.iter_mut().enumerate() { *q = trace[2 * sub + i]; }
                let p = singles[i].step_trace(&one).unwrap();
                assert!((p - components[i]).abs() < 2e-12 * (1.0 + p.abs()));
                expected += p;
            }
            assert!((pressure - expected).abs() < 2e-12 * (1.0 + expected.abs()));
        }
    }
}

#[test]
fn endpoint_zero_attack_survives_and_refused_traces_do_not_advance_filter_history() {
    let cfg = receiver();
    let disk = [CircularOpening { center_m: [0.0, 0.0], radius_m: 0.006 }];
    let make = || BaffledPressure::circular_openings_multirate_differentiating(
        &disk, 4 * OUTPUT, OUTPUT, cfg, medium()).unwrap();
    let (mut a, mut b) = (make(), make());
    let mut peak = 0.0_f64;
    for frame in 0..600 {
        assert!(a.step_trace(&[0.0; 3]).is_err());
        assert!(a.step_trace(&[0.0, f64::NAN, 0.0, 0.0]).is_err());
        let trace = if frame == 0 { [1e-5, 0.0, 0.0, 0.0] } else { [0.0; 4] };
        let p = a.step_trace(&trace).unwrap();
        assert_eq!(p, b.step_trace(&trace).unwrap());
        if frame < a.delay_samples.0 { assert_eq!(p, 0.0); }
        peak = peak.max(p.abs());
    }
    assert!(peak > 1e-4, "a true substep flow impulse was lost by endpoint sampling");
}

#[test]
fn full_band_is_not_permission_for_incompatible_clocks_or_coarse_physical_geometry() {
    let disk = [CircularOpening { center_m: [0.0, 0.0], radius_m: 0.006 }];
    for (mechanics, output) in [(0, OUTPUT), (OUTPUT, 0), (OUTPUT + 1, OUTPUT),
        (OUTPUT, OUTPUT * 2), (OUTPUT * 17, OUTPUT)] {
        assert!(BaffledPressure::circular_openings_multirate_differentiating(
            &disk, mechanics, output, receiver(), medium()).is_err());
    }
    let mut cfg = receiver(); cfg.maximum_frequency_hz = 21_601.0;
    assert!(BaffledPressure::circular_openings_multirate_differentiating(
        &disk, OUTPUT * 4, OUTPUT, cfg, medium()).is_err());
    cfg = receiver(); cfg.radial_rings = 1; cfg.angular_points = 8;
    assert!(BaffledPressure::circular_openings_multirate_differentiating(
        &disk, OUTPUT * 4, OUTPUT, cfg, medium()).is_err());
    // The old image retains its original temporal admission, not relaxed by
    // extracting the common physical surface quadrature.
    assert!(BaffledPressure::circular_openings(&disk, OUTPUT, receiver(), medium()).is_err());
}
