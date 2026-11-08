//! The output clock observes the same complete field/duct trajectory.
use super::super::super::dg::performance::{
    ThermoviscousDgControl, ThermoviscousDgPerformance, ThermoviscousDgPerformanceConfig,
};
use super::*;
use crate::bernoulli_aperture::performance::ApertureObservation;
use crate::pcm_wav::observation::{PressureRenderer, resample::BandlimitedObservationPolicy};
use std::sync::Arc;

fn policy() -> BandlimitedObservationPolicy {
    BandlimitedObservationPolicy {
        passband_hz: 1000.0,
        stopband_hz: 8000.0,
        phase_count: 128,
        max_preparation_terms: 1 << 22,
        max_history_values: 1 << 16,
        max_output_terms: 1 << 16,
    }
}

#[test]
fn independent_clock_retains_secondary_duct_and_original_molecular_trajectory() {
    let (mesh, cube) = kuhn_cube(1);
    let xyz = points(&cube);
    let gas = gas();
    let interior = ThermoviscousInterior::from_relaxation(
        three_port_field(&mesh, &xyz, &gas),
        source(gas, true, 1.0),
        BAND,
    )
    .unwrap();
    let dt = 7.3e-6; // A physical interval with no integer-Hz reciprocal.
    for rate in [48000, 192000] {
        let output_samples = 48;
        let intervals = (output_samples as f64 / (dt * f64::from(rate))).ceil() as usize;
        let gate = Arc::new(CancelGate::new_clock_free());
        let controls: Vec<_> = (0..intervals)
            .map(|i| ThermoviscousDgControl {
                traction_pa: vec![0.0, 0.0, if i < 8 { 0.1 } else { -0.02 }],
                wall_forces_n_per_sqrt_kg: vec![],
            })
            .collect();
        let make = || {
            let prepared = interior
                .prepare_time(dt, ThermoviscousTimeLimits::default(), &mut || false)
                .unwrap();
            let mut initial = prepared.zero_state();
            initial.molecular_temperature_k[0][0] = 0.002;
            ThermoviscousDg::new(
                &interior,
                prepared,
                time_dg_load(&gas, dt, DgBoundary::Closed),
                1,
                ThermoviscousDgLimits::default(),
                &gate,
            )
            .unwrap()
            .with_additional_duct(0, time_dg_load(&gas, dt, DgBoundary::Matched), &gate)
            .unwrap()
            .with_initial_state(initial)
            .unwrap()
        };
        let config = ThermoviscousDgPerformanceConfig {
            sample_rate_hz: rate,
            maximum_block: 48,
            maximum_control_values: 3 * intervals,
        };
        assert!(
            ThermoviscousDgPerformance::new_at_port(
                make(),
                0,
                controls.clone(),
                ApertureObservation::Inlet,
                config,
                Arc::clone(&gate)
            )
            .is_err(),
            "the original exact-clock constructor must remain strict"
        );
        let mut whole = ThermoviscousDgPerformance::new_bandlimited_at_port(
            make(),
            0,
            controls.clone(),
            ApertureObservation::Inlet,
            config,
            output_samples,
            policy(),
            Arc::clone(&gate),
        )
        .unwrap();
        let mut divided = ThermoviscousDgPerformance::new_bandlimited_at_port(
            make(),
            0,
            controls.clone(),
            ApertureObservation::Inlet,
            ThermoviscousDgPerformanceConfig {
                maximum_block: 5,
                ..config
            },
            output_samples,
            policy(),
            Arc::clone(&gate),
        )
        .unwrap();
        let filter = whole.bandlimited_observation().unwrap();
        assert_eq!(filter.report().input_period_s.to_bits(), dt.to_bits());
        assert_eq!(filter.report().output_rate_hz, rate);
        assert_eq!(filter.report().channels, 1);
        assert_eq!(filter.source_intervals(), intervals);
        let mut a = [0.0; 48];
        let mut b = [0.0; 48];
        whole.block(&mut a).unwrap();
        for chunk in b.chunks_mut(5) {
            divided.block(chunk).unwrap();
        }
        assert_eq!(a, b);
        assert!(a.iter().any(|p| *p != 0.0));

        // Replay the SAME original controls directly through the physical
        // owner, with no observation. No field, duct or molecular state may
        // depend on the output clock, filter latency or callback partition.
        let mut original = make();
        for control in &controls {
            original
                .step(
                    &control.traction_pa,
                    &control.wall_forces_n_per_sqrt_kg,
                    &gate,
                )
                .unwrap();
        }
        for observed in [&whole, &divided] {
            assert_eq!(observed.system().state(), original.state());
            assert_eq!(observed.controls(), controls.as_slice());
            assert_eq!(observed.observation_port(), 0);
            assert_eq!(observed.system().accepted_steps(), intervals as u64);
            assert_eq!(observed.samples_rendered(), output_samples);
            assert_eq!(observed.remaining_samples(), 0);
            for port in [0, 1] {
                assert_eq!(
                    observed.system().load_at_port(port).unwrap().state(),
                    original.load_at_port(port).unwrap().state()
                );
            }
            assert_eq!(
                observed.last_frame().unwrap().balance_residual_j,
                original.last_frame().unwrap().balance_residual_j
            );
        }
    }
}
