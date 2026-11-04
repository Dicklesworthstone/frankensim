//! A translating volume source through full physical fields and a cardioid.
//!
//! Run: `cargo run -p fs-couple --example moving_monopole_pickup -- output.wav`
//!
//! The prescribed 120 Hz volume injection is sampled on its ORIGINAL source
//! clock, including exact time derivatives. The observer retains retarded
//! propagation, source motion and reactive near-field velocity. The ideal
//! cardioid uses p and that velocity, so its proximity response is physical.
//! The explicit 0.01 Pa full-scale is never adjusted to the observed peak.
//! This example assumes linear compact radiation in a homogeneous lossless
//! fluid; it does not model a finite loudspeaker, capsule mechanics or wind.

use fs_aeroac::curle3d::moving::{
    CompactSourceStrength, MovingCompactLimits, MovingCompactObserver, TimedCompactSource,
    UniformSourceMotion,
};
use fs_bem::near_field::FirstOrder;
use fs_couple::pcm_wav::encode_pcm16_wav;
use fs_couple::pcm_wav::observation::moving::field::{CompactFieldSample, MovingMonopoleRenderer};
use fs_couple::pcm_wav::observation::moving::{
    CompactSourceHistory, MovingPressureBudget, MovingPressureConfig,
};
use fs_material::fluid_acoustics::FluidAcousticState;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let output_path = std::env::args()
        .nth(1)
        .unwrap_or_else(|| "moving_monopole_pickup.wav".into());
    let frequency_hz = 120.0;
    let omega = std::f64::consts::TAU * frequency_hz;
    let amplitude_m3_s = 1e-4;
    let source_rate_hz = 16_000.0;
    // Retain explicit prehistory and a finite source horizon. No missing
    // endpoint, initial field or audible tail is synthesized by the observer.
    let history: Vec<_> = (0..=8000)
        .map(|index| {
            let source_time_s = -0.1 + f64::from(index) / source_rate_hz;
            let (sin, cos) = fs_math::det::sin_cos(omega * source_time_s);
            TimedCompactSource {
                source_time_s,
                strength: CompactSourceStrength {
                    volume_flow_m3_s: amplitude_m3_s * sin,
                    volume_flow_rate_m3_s2: amplitude_m3_s * omega * cos,
                    ..Default::default()
                },
            }
        })
        .collect();
    let medium = FluidAcousticState::try_new(1.2, 343.0)?;
    let observer = MovingCompactObserver::new(
        medium,
        [6.0, 0.0, 1.0],
        MovingCompactLimits {
            maximum_mach: 0.1,
            minimum_distance_m: 0.1,
        },
    )?;
    let sources = [CompactSourceHistory {
        motion: UniformSourceMotion::new(0.0, [0.0, 1.5, 0.0], [20.0, 0.0, 0.0])?,
        samples: &history,
        source_bandwidth_hz: frequency_hz,
    }];
    let microphone = FirstOrder::new(0.5, [-1.0, 0.0, 0.0])?;
    let config = MovingPressureConfig {
        first_receiver_time_s: 0.0,
        sample_rate_hz: 8000,
        samples: 2000,
        max_block: 128,
        // Includes Doppler and geometric modulation. Source interpolation is
        // not an antialias filter; refine accepted sampling when changing it.
        receiver_bandwidth_hz: 512.0,
        maximum_time_residual_s: 1e-12,
    };
    let mut renderer = MovingMonopoleRenderer::new(
        observer,
        &sources,
        microphone,
        config,
        MovingPressureBudget {
            maximum_sources: 1,
            maximum_source_samples: history.len(),
            maximum_output_samples: config.samples,
            maximum_payload_bytes: 2_000_000,
            maximum_source_receiver_pairs: config.samples,
        },
        &mut || false,
    )
    .expect("admit the finite source history, receiver clock and field storage");
    let mut fields = [CompactFieldSample::default(); 128];
    let mut pickup = Vec::with_capacity(config.samples as usize);
    let mut pressure_square = 0.0;
    let mut pickup_square = 0.0;
    let mut peak_velocity: f64 = 0.0;
    while renderer.remaining_samples() > 0 {
        let count = renderer.remaining_samples().min(fields.len() as u64) as usize;
        renderer
            .field_block(&mut fields[..count])
            .expect("observe the same physical source clock");
        for field in &fields[..count] {
            pickup.push(field.microphone_pressure_pa);
            pressure_square += field.pressure_pa * field.pressure_pa;
            pickup_square += field.microphone_pressure_pa * field.microphone_pressure_pa;
            let speed = field
                .particle_velocity_m_s
                .iter()
                .map(|u| u * u)
                .sum::<f64>()
                .sqrt();
            peak_velocity = peak_velocity.max(speed);
        }
    }
    let full_scale_pa = 0.01;
    let (wav, clipped_samples) = encode_pcm16_wav(&pickup, config.sample_rate_hz, full_scale_pa)?;
    std::fs::write(&output_path, wav)?;
    println!(
        "Wrote {output_path}: {} samples at {} Hz",
        config.samples, config.sample_rate_hz
    );
    println!("Source: {frequency_hz} Hz, peak Q={amplitude_m3_s:e} m³/s, translation=20 m/s");
    println!(
        "Pressure RMS: {:.8e} Pa",
        (pressure_square / config.samples as f64).sqrt()
    );
    println!(
        "Cardioid RMS: {:.8e} Pa",
        (pickup_square / config.samples as f64).sqrt()
    );
    println!("Peak particle speed: {peak_velocity:.8e} m/s");
    println!("PCM full-scale: {full_scale_pa} Pa; clipped samples: {clipped_samples}");
    Ok(())
}
