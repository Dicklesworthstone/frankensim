//! Original body-force/heat work remains in the same physical trajectory.
use super::*;
use crate::pcm_wav::observation::resample::BandlimitedObservationPolicy;
use crate::render::RenderError;

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
fn bandlimited_volume_playback_preserves_all_source_work_and_both_duct_histories() {
    let (mesh, cube) = kuhn_cube(1);
    let xyz = points(&cube);
    let gas = gas();
    let interior = ThermoviscousInterior::from_relaxation(
        build(&mesh, &xyz, &gas, false),
        source(gas, true, 1.0),
        BAND,
    )
    .unwrap();
    let (force, heat) = profiles(&interior);
    let dt = 1.0 / 96000.0;
    let rate = 44100;
    let count = 64;
    let intervals = (count as f64 / (dt * f64::from(rate))).ceil() as usize;
    let gate = Arc::new(CancelGate::new_clock_free());
    let make = || {
        system(
            &interior,
            &force,
            &heat,
            dt,
            ThermoviscousDgLimits::default(),
            &gate,
        )
        .with_additional_duct(0, time_dg_load(&gas, dt, DgBoundary::Matched), &gate)
        .unwrap()
    };
    let controls: Vec<_> = (0..intervals)
        .map(|i| ThermoviscousDgVolumeControl {
            traction_pa: vec![0.0; 2],
            wall_forces_n_per_sqrt_kg: vec![],
            force_drive_n: vec![if i < 8 {
                1e-6
            } else if i < 16 {
                -5e-7
            } else {
                0.0
            }],
            heat_drive_w: vec![if i < 5 {
                1e-3
            } else if (12..16).contains(&i) {
                -2e-4
            } else {
                0.0
            }],
        })
        .collect();
    let config = ThermoviscousDgPerformanceConfig {
        sample_rate_hz: rate,
        maximum_block: 64,
        maximum_control_values: 4 * intervals,
    };
    let mut whole = ThermoviscousVolumeDgPerformance::new_bandlimited_at_port(
        make(),
        0,
        controls.clone(),
        ApertureObservation::Inlet,
        config,
        count,
        policy(),
        Arc::clone(&gate),
    )
    .unwrap();
    let mut divided = ThermoviscousVolumeDgPerformance::new_bandlimited_at_port(
        make(),
        0,
        controls.clone(),
        ApertureObservation::Inlet,
        ThermoviscousDgPerformanceConfig {
            maximum_block: 5,
            ..config
        },
        count,
        policy(),
        Arc::clone(&gate),
    )
    .unwrap();
    let mut a = [0.0; 64];
    let mut b = [0.0; 64];
    whole.block(&mut a).unwrap();
    for chunk in b.chunks_mut(5) {
        divided.block(chunk).unwrap();
    }
    assert_eq!(a, b);
    assert!(a.iter().any(|p| *p != 0.0));
    let mut original = make();
    for control in &controls {
        let frame = original
            .step_with_volume_sources(
                &control.traction_pa,
                &control.wall_forces_n_per_sqrt_kg,
                &control.force_drive_n,
                &control.heat_drive_w,
                &gate,
            )
            .unwrap();
        source_balance(frame);
    }
    for observed in [&whole, &divided] {
        assert_eq!(observed.system().state(), original.state());
        assert_eq!(
            observed.last_volume_observations(),
            original.last_volume_observations()
        );
        assert_eq!(observed.controls(), controls.as_slice());
        assert_eq!(observed.system().accepted_steps(), intervals as u64);
        assert_eq!(observed.samples_rendered(), count);
        assert_eq!(observed.remaining_samples(), 0);
        for port in [0, 1] {
            assert_eq!(
                observed.system().load_at_port(port).unwrap().state(),
                original.load_at_port(port).unwrap().state()
            );
        }
        source_balance(observed.last_frame().unwrap());
    }
    assert!(
        ThermoviscousVolumeDgPerformance::new_bandlimited(
            make(),
            controls[..intervals - 1].to_vec(),
            ApertureObservation::Inlet,
            config,
            count,
            policy(),
            Arc::clone(&gate)
        )
        .is_err()
    );

    // At 2x upsampling the second output needs NO new physical interval. The
    // retained gate must still cancel it, leaving the previous physical/source
    // work and caller buffer intact before poisoning further continuation.
    let mut stream = ThermoviscousVolumeDgPerformance::new_bandlimited(
        make(),
        controls[..2].to_vec(),
        ApertureObservation::Inlet,
        ThermoviscousDgPerformanceConfig {
            sample_rate_hz: 192000,
            ..config
        },
        4,
        policy(),
        Arc::clone(&gate),
    )
    .unwrap();
    stream.block(&mut [0.0]).unwrap();
    let physical = stream.system().state().clone();
    let source_work = stream.last_volume_observations().cloned();
    gate.request();
    let mut untouched = [f64::MAX];
    assert!(stream.block(&mut untouched).is_err());
    assert_eq!(untouched, [f64::MAX]);
    assert_eq!(stream.system().state(), &physical);
    assert_eq!(stream.last_volume_observations(), source_work.as_ref());
    assert_eq!(stream.system().accepted_steps(), 1);
    assert_eq!(stream.samples_rendered(), 1);
    assert!(matches!(
        stream.block(&mut untouched),
        Err(RenderError::Poisoned)
    ));
}
