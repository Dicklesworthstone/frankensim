//! Receiver-clock reconstruction of finite dry stems on subsonic linear tracks.
//!
//! This is a kinematic observation of a supplied digital waveform. The source
//! clock follows t = tau + |x_receiver(t) - x_source(tau)| / c, hence
//! d tau / dt = (1 - n.v_receiver/c) / (1 - n.v_source/c). No second pitch or
//! amplitude factor is applied. A moving monopole's radiation coefficients
//! belong to fs-aeroac; they cannot be inferred from a dry modal stem.

use super::*;

#[derive(Clone, Copy)]
struct VelocityEnvelope {
    low: [f64; 3],
    high: [f64; 3],
    speed: f64,
}

fn envelope(
    frames: usize,
    position: impl Fn(usize) -> [f64; 3],
    sound_step: f64,
    checkpoint_fn: &mut impl FnMut() -> Result<(), SpatialAudioError>,
) -> Result<VelocityEnvelope, SpatialAudioError> {
    let mut result = VelocityEnvelope {
        low: [f64::INFINITY; 3],
        high: [f64::NEG_INFINITY; 3],
        speed: 0.0,
    };
    for i in 1..frames {
        if i % SPATIAL_AUDIO_CANCELLATION_POLL_FRAMES == 0 {
            checkpoint_fn()?;
        }
        let a = position(i - 1);
        let b = position(i);
        let velocity = core::array::from_fn(|j| (b[j] - a[j]) / sound_step);
        let speed = norm(velocity);
        if !speed.is_finite() || speed >= 1.0 {
            return Err(invalid(
                "strictly subsonic continuous source and receiver tracks",
            ));
        }
        result.speed = result.speed.max(speed);
        for j in 0..3 {
            result.low[j] = result.low[j].min(velocity[j]);
            result.high[j] = result.high[j].max(velocity[j]);
        }
    }
    if frames < 2 {
        result.low = [0.0; 3];
        result.high = [0.0; 3];
    }
    Ok(result)
}

pub(super) fn admit(
    config: SpatialAudioConfig,
    input: SpatialAudioRenderInput<'_>,
    frames: usize,
    checkpoint_fn: &mut impl FnMut() -> Result<(), SpatialAudioError>,
) -> Result<u64, SpatialAudioError> {
    if frames < 2 || frames as u64 > (1_u64 << 53) {
        return Err(invalid(
            "retarded source history with representable sample clocks",
        ));
    }
    let SpatialDelayPolicy::RetardedSinc {
        source_bandwidth_hz,
        receiver_bandwidth_hz,
    } = config.delay_policy
    else {
        return Err(invalid("retarded delay policy"));
    };
    let sound_step = config.speed_of_sound_m_per_s / f64::from(config.sample_rate_hz);
    if !sound_step.is_finite() || sound_step <= 0.0 {
        return Err(invalid("representable sound travel per sample"));
    }
    let receiver_frames = match input.listener {
        ListenerPoseTrack::Static(_) => 1,
        ListenerPoseTrack::PerFrame(poses) => poses.len(),
    };
    let receiver = envelope(
        receiver_frames,
        |i| input.listener.pose(i).position_m,
        sound_step,
        checkpoint_fn,
    )?;
    let mut last_arrival = (frames - 1) as f64;
    for source in input.sources {
        checkpoint_fn()?;
        let source_frames = if matches!(source.positions, SourcePositionTrack::Static(_)) {
            1
        } else {
            frames
        };
        let motion = envelope(
            source_frames,
            |i| source.positions.position(i),
            sound_step,
            checkpoint_fn,
        )?;
        // Both envelopes cover complete linear interpolants, not just endpoint
        // instantaneous Doppler samples. The relative-velocity bound preserves
        // common translation without pretending accelerated simultaneous
        // velocities equal their different emission/reception-time values.
        let difference = norm(core::array::from_fn(|j| {
            (motion.low[j] - receiver.high[j])
                .abs()
                .max((motion.high[j] - receiver.low[j]).abs())
        }));
        let maximum_rate = if difference == 0.0 {
            1.0
        } else {
            let margin = (1.0 - motion.speed.next_up()).next_down();
            if margin <= 0.0 {
                return Err(invalid("representable subsonic speed margin"));
            }
            ((1.0 + receiver.speed.next_up()).next_up() / margin)
                .min((1.0 + difference.next_up() / margin).next_up())
                .next_up()
        };
        if !maximum_rate.is_finite() || source_bandwidth_hz * maximum_rate > receiver_bandwidth_hz {
            return Err(invalid(
                "receiver bandwidth covering both source and receiver Doppler",
            ));
        }
        let arrival = last_reception(
            source.positions.position(frames - 1),
            (frames - 1) as f64,
            input.listener,
            sound_step,
        )?;
        last_arrival = last_arrival.max(arrival);
    }
    if !last_arrival.is_finite() || last_arrival >= (1_u64 << 53) as f64 {
        return Err(invalid("representable retarded output horizon"));
    }
    Ok(last_arrival.ceil() as u64 + 1)
}

fn last_reception(
    source: [f64; 3],
    emission: f64,
    listener: ListenerPoseTrack<'_>,
    sound_step: f64,
) -> Result<f64, SpatialAudioError> {
    let ListenerPoseTrack::PerFrame(poses) = listener else {
        let arrival = emission + distance(source, listener.pose(0).position_m) / sound_step;
        validate_finite(arrival, "final retarded arrival", emission as usize)?;
        return Ok(arrival);
    };
    let residual = |i: usize| {
        let value = i as f64 - emission - distance(source, poses[i].position_m) / sound_step;
        validate_finite(value, "retarded receiver-history residual", i)?;
        Ok::<_, SpatialAudioError>(value)
    };
    let mut low = emission as usize;
    let mut high = poses.len() - 1;
    if residual(high)? < 0.0 {
        return Err(invalid(
            "listener history covering the final retarded arrival",
        ));
    }
    if residual(low)? == 0.0 {
        return Ok(emission);
    }
    while high - low > 1 {
        let mid = low + (high - low) / 2;
        if residual(mid)? < 0.0 {
            low = mid;
        } else {
            high = mid;
        }
    }
    if residual(high)? == 0.0 {
        return Ok(high as f64);
    }
    let q = core::array::from_fn(|j| (poses[low].position_m[j] - source[j]) / sound_step);
    let beta = core::array::from_fn(|j| {
        (poses[high].position_m[j] - poses[low].position_m[j]) / sound_step
    });
    let elapsed = low as f64 - emission;
    let a = 1.0 - dot(beta, beta);
    let b = elapsed - dot(q, beta);
    let radius = norm(q);
    let c = (elapsed - radius) * (elapsed + radius);
    let discriminant = b * b - a * c;
    let root = det::sqrt(discriminant);
    let fraction = if b >= 0.0 {
        -c / (b + root)
    } else {
        (root - b) / a
    };
    let fraction = refine_root(fraction, |x| {
        elapsed + x - norm(core::array::from_fn(|j| q[j] + beta[j] * x))
    })?;
    Ok(low as f64 + fraction)
}

pub(super) fn render(
    config: SpatialAudioConfig,
    input: SpatialAudioRenderInput<'_>,
    output: &mut [StereoAccumulator],
    preflight: &mut Preflight,
    checkpoint_fn: &mut impl FnMut() -> Result<(), SpatialAudioError>,
) -> Result<(), SpatialAudioError> {
    let sound_step = config.speed_of_sound_m_per_s / f64::from(config.sample_rate_hz);
    for source in input.sources.iter().copied() {
        let mut span = 0;
        for (frame, accumulator) in output.iter_mut().enumerate() {
            if frame % SPATIAL_AUDIO_CANCELLATION_POLL_FRAMES == 0 {
                checkpoint_fn()?;
            }
            let listener = input.listener.pose(frame);
            let Some((emission, source_position)) = emission_at(
                source.positions,
                source.signal.len(),
                listener.position_m,
                frame as f64,
                sound_step,
                &mut span,
                checkpoint_fn,
            )?
            else {
                continue;
            };
            let geometry = geometry(
                source_position,
                listener,
                config.minimum_distance_m,
                config.speed_of_sound_m_per_s,
                config.sample_rate_hz,
                config.microphone_directivity,
                frame,
            )?;
            preflight.maximum_distance_m = preflight.maximum_distance_m.max(geometry.distance_m);
            preflight.maximum_delay_frames =
                preflight.maximum_delay_frames.max(geometry.delay_frames);
            if geometry.distance_m < config.minimum_distance_m {
                preflight.minimum_distance_clamp_count += 1;
            }
            let value = sinc_sample(source.signal, emission)?
                * source.gain_linear
                * geometry.attenuation_gain
                * geometry.microphone_gain;
            accumulator.add(
                value * det::sqrt(0.5 * (1.0 - geometry.pan)),
                value * det::sqrt(0.5 * (1.0 + geometry.pan)),
            );
        }
    }
    checkpoint_fn()
}

fn emission_at(
    positions: SourcePositionTrack<'_>,
    frames: usize,
    receiver: [f64; 3],
    reception: f64,
    sound_step: f64,
    span: &mut usize,
    checkpoint_fn: &mut impl FnMut() -> Result<(), SpatialAudioError>,
) -> Result<Option<(f64, [f64; 3])>, SpatialAudioError> {
    if let SourcePositionTrack::Static(position) = positions {
        let emission = reception - distance(position, receiver) / sound_step;
        validate_finite(emission, "retarded emission time", reception as usize)?;
        return Ok(
            (emission >= 0.0 && emission <= (frames - 1) as f64).then_some((emission, position))
        );
    }
    let residual = |i: usize| {
        let value = i as f64 + distance(positions.position(i), receiver) / sound_step - reception;
        validate_finite(value, "retarded source-history residual", i)?;
        Ok::<_, SpatialAudioError>(value)
    };
    if residual(*span)? > 0.0 {
        if *span == 0 {
            return Ok(None);
        }
        return Err(invalid("monotone subsonic retarded source clock"));
    }
    // With both speeds below c, emission time increases with reception time.
    // The cursor crosses each source span once over the complete render.
    while *span + 1 < frames && residual(*span + 1)? <= 0.0 {
        *span += 1;
        if *span % SPATIAL_AUDIO_CANCELLATION_POLL_FRAMES == 0 {
            checkpoint_fn()?;
        }
    }
    if residual(*span)? == 0.0 {
        return Ok(Some((*span as f64, positions.position(*span))));
    }
    if *span + 1 == frames {
        return Ok(None);
    }
    let start = positions.position(*span);
    let end = positions.position(*span + 1);
    let q = core::array::from_fn(|j| (receiver[j] - start[j]) / sound_step);
    let beta = core::array::from_fn(|j| (end[j] - start[j]) / sound_step);
    let elapsed = reception - *span as f64;
    let a = 1.0 - dot(beta, beta);
    let b = elapsed - dot(q, beta);
    let radius = norm(q);
    let c = (elapsed - radius) * (elapsed + radius);
    // Rationalized smaller quadratic root avoids subtracting nearly equal
    // flight times. Source interpolation uses this same segment and fraction.
    let fraction = c / (b + det::sqrt(b * b - a * c));
    let fraction = refine_root(fraction, |x| {
        x + norm(core::array::from_fn(|j| q[j] - beta[j] * x)) - elapsed
    })?;
    Ok(Some((
        *span as f64 + fraction,
        core::array::from_fn(|j| start[j] + (end[j] - start[j]) * fraction),
    )))
}

// Ordinary paths use one analytical root. A rounded root outside its actual
// segment or with a large residual receives bounded bisection, never clamping.
fn refine_root(candidate: f64, residual: impl Fn(f64) -> f64) -> Result<f64, SpatialAudioError> {
    if candidate.is_finite()
        && (0.0..=1.0).contains(&candidate)
        && residual(candidate).abs() <= 1.0e-9
    {
        return Ok(candidate);
    }
    let mut low = 0.0;
    let mut high = 1.0;
    let (lo, hi) = (residual(low), residual(high));
    if !lo.is_finite() || !hi.is_finite() || lo > 0.0 || hi < 0.0 {
        return Err(invalid("representable retarded segment bracket"));
    }
    if lo == 0.0 {
        return Ok(low);
    }
    if hi == 0.0 {
        return Ok(high);
    }
    for _ in 0..48 {
        let mid = 0.5 * (low + high);
        let value = residual(mid);
        if !value.is_finite() {
            return Err(invalid("finite retarded segment residual"));
        }
        if value == 0.0 {
            return Ok(mid);
        }
        if value < 0.0 {
            low = mid;
        } else {
            high = mid;
        }
    }
    let result = 0.5 * (low + high);
    let error = residual(result);
    if !error.is_finite() || error.abs() > 1.0e-9 {
        return Err(invalid("retarded reconstruction time residual"));
    }
    Ok(result)
}

fn sinc_sample(signal: SpatialMonoSignal<'_>, emission: f64) -> Result<f64, SpatialAudioError> {
    let center = emission.floor() as usize;
    let fraction = emission - center as f64;
    if fraction == 0.0 {
        return Ok(signal.sample(center));
    }
    // sin(pi (fraction-j)) alternates sign. The Lanczos window advances with
    // a fixed rotation, so all thirty-two taps need only two trig evaluations.
    let pi = core::f64::consts::PI;
    let mut carrier_sine = -det::sin(pi * fraction);
    let (mut window_sine, mut window_cosine) = det::sin_cos(pi * (fraction + 15.0) / 16.0);
    let (step_sine, step_cosine) = (0.195_090_322_016_128_25, 0.980_785_280_403_230_4);
    let mut weighted = CompensatedSum::ZERO;
    let mut weights = CompensatedSum::ZERO;
    for j in -15_i32..=16 {
        let x = fraction - f64::from(j);
        let weight = if x.abs() < 1.0e-4 {
            let y = pi * x;
            let sinc = 1.0 - y * y / 6.0 + y.powi(4) / 120.0;
            let z = y / 16.0;
            sinc * (1.0 - z * z / 6.0 + z.powi(4) / 120.0)
        } else {
            16.0 * carrier_sine * window_sine / (pi * pi * x * x)
        };
        weights.add(weight);
        if let Some(index) = center.checked_add_signed(j as isize)
            && index < signal.len()
        {
            weighted.add(weight * signal.sample(index));
        }
        let next_sine = window_sine * step_cosine - window_cosine * step_sine;
        window_cosine = window_cosine * step_cosine + window_sine * step_sine;
        window_sine = next_sine;
        carrier_sine = -carrier_sine;
    }
    let result = weighted.finish() / weights.finish();
    validate_finite(result, "retarded interpolated dry source", center)?;
    Ok(result)
}

fn norm(vector: [f64; 3]) -> f64 {
    det::sqrt(dot(vector, vector))
}
fn distance(a: [f64; 3], b: [f64; 3]) -> f64 {
    norm(core::array::from_fn(|j| a[j] - b[j]))
}

#[cfg(test)]
mod tests {
    use super::*;

    /// G1: independent complex-exponential response across the admitted band.
    /// This sampled numerical envelope is not a continuum filter certificate.
    #[test]
    fn lanczos_response_respects_the_admitted_source_sampling_margin() {
        for bin in 0..=128 {
            let frequency = 0.25 * f64::from(bin) / 128.0;
            let cosine: [f64; 128] =
                core::array::from_fn(|i| (core::f64::consts::TAU * frequency * i as f64).cos());
            let sine: [f64; 128] =
                core::array::from_fn(|i| (core::f64::consts::TAU * frequency * i as f64).sin());
            for phase in 0..=32 {
                let time = 64.0 + f64::from(phase) / 32.0;
                let re = sinc_sample(SpatialMonoSignal::Samples(&cosine), time).unwrap();
                let im = sinc_sample(SpatialMonoSignal::Samples(&sine), time).unwrap();
                let angle = core::f64::consts::TAU * frequency * time;
                let error = ((re - angle.cos()).powi(2) + (im - angle.sin()).powi(2)).sqrt();
                assert!(
                    error < 0.001,
                    "frequency={frequency}, phase={phase}, error={error}"
                );
            }
        }
    }

    #[test]
    fn sinc_integer_samples_and_constant_amplitude_are_preserved() {
        let constant = [0.25; 128];
        for fraction in [0.0, f64::EPSILON * 32.0, 0.01, 0.25, 0.5, 0.75, 0.99] {
            let value =
                sinc_sample(SpatialMonoSignal::Samples(&constant), 64.0 + fraction).unwrap();
            assert!((value - 0.25).abs() < 1.0e-14);
        }
        let samples: [f64; 128] = core::array::from_fn(|i| (i as f64 * 0.19).sin());
        for i in 0..samples.len() {
            assert_eq!(
                sinc_sample(SpatialMonoSignal::Samples(&samples), i as f64)
                    .unwrap()
                    .to_bits(),
                samples[i].to_bits()
            );
        }
    }
}
