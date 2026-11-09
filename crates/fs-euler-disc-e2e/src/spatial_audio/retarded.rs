//! Receiver-clock reconstruction of finite dry stems on subsonic linear tracks.
//!
//! This is a kinematic observation of a supplied digital waveform. The source
//! clock follows |x_receiver(t) - x_source(tau) - U*(t-tau)| = c*(t-tau), hence
//! d tau / dt = (1 - n.(v_receiver-U)/c) / (1 - n.(v_source-U)/c). U is the
//! uniform medium velocity; n is the acoustic ray direction in that medium's
//! rest frame. Only event differences are transformed, avoiding subtraction of
//! separately translated absolute positions. No second pitch or amplitude
//! factor is applied. A moving monopole's radiation coefficients belong to
//! fs-aeroac; they cannot be inferred from a dry modal stem.

use super::*;
use fs_ivl::Interval;

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
    flow_step: [f64; 3],
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
        let velocity = core::array::from_fn(|j| (b[j] - a[j] - flow_step[j]) / sound_step);
        let speed = norm(velocity);
        if !speed.is_finite() || speed >= 1.0 {
            return Err(invalid(
                "strictly subsonic source and receiver tracks relative to the medium",
            ));
        }
        result.speed = result.speed.max(speed);
        for j in 0..3 {
            result.low[j] = result.low[j].min(velocity[j]);
            result.high[j] = result.high[j].max(velocity[j]);
        }
    }
    if frames < 2 {
        let velocity = flow_step.map(|v| -v / sound_step);
        result.speed = norm(velocity);
        if !result.speed.is_finite() || result.speed >= 1.0 {
            return Err(invalid(
                "strictly subsonic source and receiver tracks relative to the medium",
            ));
        }
        result.low = velocity;
        result.high = velocity;
    }
    Ok(result)
}

pub(super) fn admit(
    config: SpatialAudioConfig,
    input: SpatialAudioRenderInput<'_>,
    frames: usize,
    checkpoint_fn: &mut impl FnMut() -> Result<(), SpatialAudioError>,
) -> Result<(u64, u64), SpatialAudioError> {
    if frames < 2 || frames as u64 > (1_u64 << 53) {
        return Err(invalid(
            "retarded source history with representable sample clocks",
        ));
    }
    let (source_bandwidth_hz, receiver_bandwidth_hz) = config
        .delay_policy
        .retarded_bands()
        .ok_or_else(|| invalid("retarded delay policy"))?;
    let sound_step = config.speed_of_sound_m_per_s / f64::from(config.sample_rate_hz);
    let flow_step = config
        .delay_policy
        .medium_velocity_m_per_s()
        .map(|v| v / f64::from(config.sample_rate_hz));
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
        flow_step,
        checkpoint_fn,
    )?;
    let mut last_arrival = (frames - 1) as f64;
    let mut refinement_work = 0_u64;
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
            flow_step,
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
        let mut maximum_rate = if difference == 0.0 {
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
        // Refine only a bandwidth refusal, so already-admitted scenes retain
        // their preparation cost. Against a fixed endpoint in resting air,
        // the first radial projection bounds each COMPLETE linear span.
        if flow_step == [0.0; 3]
            && source_bandwidth_hz * maximum_rate > receiver_bandwidth_hz
        {
            let refinement_frames = match (input.listener, source.positions) {
                (ListenerPoseTrack::Static(_), _) => source_frames,
                (_, SourcePositionTrack::Static(_)) => receiver_frames,
                _ => 0,
            };
            if refinement_frames > 1 {
                refinement_work = refinement_work
                    .checked_add((refinement_frames - 1) as u64)
                    .ok_or_else(|| invalid("directional admission work"))?;
                let source_count = input.sources.len() as u64;
                let total_input_frames = source_count
                    .checked_mul(frames as u64)
                    .ok_or_else(|| invalid("spatial source work"))?;
                let listener_frames = match input.listener {
                    ListenerPoseTrack::Static(_) => frames,
                    ListenerPoseTrack::PerFrame(poses) => poses.len(),
                };
                // Before scanning, charge the already-known minimum work.
                // Both output policies retain at least the input horizon;
                // preflight charges the exact propagated/convolved tail later.
                validate_work_budget(
                    total_input_frames,
                    frames as u64,
                    frames as u64,
                    input.room_ir.map_or(0, StereoRoomImpulseResponse::len),
                    retarded_work_units(
                        source_count,
                        frames as u64,
                        listener_frames as u64,
                        refinement_work,
                    )?,
                    config.budget,
                )?;
            }
            if let ListenerPoseTrack::Static(listener) = input.listener {
                maximum_rate = fixed_endpoint_maximum_rate(
                    source_frames,
                    |i| source.positions.position(i),
                    listener.position_m,
                    sound_step,
                    true,
                    maximum_rate,
                    checkpoint_fn,
                )?;
            } else if let SourcePositionTrack::Static(position) = source.positions {
                maximum_rate = fixed_endpoint_maximum_rate(
                    receiver_frames,
                    |i| input.listener.pose(i).position_m,
                    position,
                    sound_step,
                    false,
                    maximum_rate,
                    checkpoint_fn,
                )?;
            }
        }
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
            flow_step,
        )?;
        last_arrival = last_arrival.max(arrival);
    }
    if !last_arrival.is_finite() || last_arrival >= (1_u64 << 53) as f64 {
        return Err(invalid("representable retarded output horizon"));
    }
    Ok((last_arrival.ceil() as u64 + 1, refinement_work))
}

/// Tighten the complete-history bound when one endpoint is fixed in resting
/// air. For either moving endpoint let r = fixed - moving and beta = v/c.
/// On a linear span, d(n.beta)/ds = -(|beta|² - (n.beta)²)/|r| <= 0
/// (s is distance travelled by sound). Thus the source rate 1/(1-n.beta)
/// and the receiver rate 1+n.beta both attain their maximum at span start.
/// A crossing only jumps from approach to recession; every new span is still
/// checked, so velocity changes cannot hide behind the first history sample.
fn fixed_endpoint_maximum_rate(
    frames: usize,
    position: impl Fn(usize) -> [f64; 3],
    fixed: [f64; 3],
    sound_step: f64,
    moving_source: bool,
    fallback: f64,
    checkpoint_fn: &mut impl FnMut() -> Result<(), SpatialAudioError>,
) -> Result<f64, SpatialAudioError> {
    let mut maximum = 0.0_f64;
    for i in 1..frames {
        if i % SPATIAL_AUDIO_CANCELLATION_POLL_FRAMES == 0 {
            checkpoint_fn()?;
        }
        let start = position(i - 1);
        let end = position(i);
        let rate = fixed_endpoint_segment_rate(start, end, fixed, sound_step, moving_source)
            .unwrap_or(fallback);
        maximum = maximum.max(rate);
    }
    Ok(if frames < 2 { 1.0 } else { maximum.min(fallback) })
}

fn fixed_endpoint_segment_rate(
    start: [f64; 3],
    end: [f64; 3],
    fixed: [f64; 3],
    sound_step: f64,
    moving_source: bool,
) -> Option<f64> {
    if start == end {
        return Some(1.0);
    }
    let ray: [Interval; 3] =
        core::array::from_fn(|axis| Interval::point(fixed[axis]) - Interval::point(start[axis]));
    let scale = ray.iter().map(|x| x.abs().hi()).fold(0.0_f64, f64::max);
    if !scale.is_finite() || scale == 0.0 {
        return None;
    }
    // Scale before squaring, and enclose every subtraction/product/division.
    // An unresolved coincident start or near-sonic denominator retains the
    // original speed-only bound instead of guessing a propagation direction.
    let ray = ray.map(|x| x / Interval::point(scale));
    let mut square = Interval::point(0.0);
    let mut radial = Interval::point(0.0);
    for axis in 0..3 {
        square = square + ray[axis] * ray[axis];
        let beta = (Interval::point(end[axis]) - Interval::point(start[axis]))
            / Interval::point(sound_step);
        radial = radial + ray[axis] * beta;
    }
    let radius = square.sqrt();
    if radius.lo() <= 0.0 || !radius.hi().is_finite() {
        return None;
    }
    radial = radial / radius;
    let one = Interval::point(1.0);
    let rate = if moving_source {
        let denominator = one - radial;
        if denominator.lo() <= 0.0 {
            return None;
        }
        (one / denominator).hi()
    } else {
        (one + radial).hi()
    };
    (rate.is_finite() && rate > 0.0).then_some(rate)
}

fn last_reception(
    source: [f64; 3],
    emission: f64,
    listener: ListenerPoseTrack<'_>,
    sound_step: f64,
    flow_step: [f64; 3],
) -> Result<f64, SpatialAudioError> {
    let ListenerPoseTrack::PerFrame(poses) = listener else {
        let arrival = emission
            + static_flight_time(source, listener.pose(0).position_m, sound_step, flow_step)?;
        validate_finite(arrival, "final retarded arrival", emission as usize)?;
        return Ok(arrival);
    };
    let residual = |i: usize| {
        let elapsed = i as f64 - emission;
        let ray = fluid_separation(source, poses[i].position_m, elapsed, flow_step);
        let value = elapsed - norm(ray) / sound_step;
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
    let elapsed = low as f64 - emission;
    let q =
        fluid_separation(source, poses[low].position_m, elapsed, flow_step).map(|v| v / sound_step);
    let beta = core::array::from_fn(|j| {
        (poses[high].position_m[j] - poses[low].position_m[j] - flow_step[j]) / sound_step
    });
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
    let flow_step = config
        .delay_policy
        .medium_velocity_m_per_s()
        .map(|v| v / f64::from(config.sample_rate_hz));
    for source in input.sources.iter().copied() {
        let mut span = 0;
        for (frame, accumulator) in output.iter_mut().enumerate() {
            if frame % SPATIAL_AUDIO_CANCELLATION_POLL_FRAMES == 0 {
                checkpoint_fn()?;
            }
            let listener = input.listener.pose(frame);
            let Some((emission, source_displacement)) = emission_at(
                source.positions,
                source.signal.len(),
                listener.position_m,
                frame as f64,
                sound_step,
                flow_step,
                &mut span,
                checkpoint_fn,
            )?
            else {
                continue;
            };
            let geometry = geometry(
                source_displacement,
                ListenerPose {
                    position_m: [0.0; 3],
                    ..listener
                },
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
    flow_step: [f64; 3],
    span: &mut usize,
    checkpoint_fn: &mut impl FnMut() -> Result<(), SpatialAudioError>,
) -> Result<Option<(f64, [f64; 3])>, SpatialAudioError> {
    if let SourcePositionTrack::Static(position) = positions {
        let delay = static_flight_time(position, receiver, sound_step, flow_step)?;
        let emission = reception - delay;
        validate_finite(emission, "retarded emission time", reception as usize)?;
        let direction = fluid_separation(position, receiver, delay, flow_step).map(|v| -v);
        return Ok(
            (emission >= 0.0 && emission <= (frames - 1) as f64).then_some((emission, direction))
        );
    }
    let residual = |i: usize| {
        let elapsed = reception - i as f64;
        let ray = fluid_separation(positions.position(i), receiver, elapsed, flow_step);
        let value = norm(ray) / sound_step - elapsed;
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
        let direction = fluid_separation(
            positions.position(*span),
            receiver,
            reception - *span as f64,
            flow_step,
        )
        .map(|v| -v);
        return Ok(Some((*span as f64, direction)));
    }
    if *span + 1 == frames {
        return Ok(None);
    }
    let start = positions.position(*span);
    let end = positions.position(*span + 1);
    let elapsed = reception - *span as f64;
    let q = fluid_separation(start, receiver, elapsed, flow_step).map(|v| v / sound_step);
    let beta = core::array::from_fn(|j| (end[j] - start[j] - flow_step[j]) / sound_step);
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
    // Reconstruct the ray from relative segment geometry. Do not add a large
    // common laboratory translation only to subtract it again for directivity.
    let direction = core::array::from_fn(|j| -(q[j] - beta[j] * fraction) * sound_step);
    Ok(Some((*span as f64 + fraction, direction)))
}

/// Propagation time between fixed laboratory events' spatial endpoints. The
/// relative endpoint in the fluid frame moves with velocity -U while the wave
/// travels. Rationalizing the downstream root avoids near-sonic cancellation.
fn static_flight_time(
    source: [f64; 3],
    receiver: [f64; 3],
    sound_step: f64,
    flow_step: [f64; 3],
) -> Result<f64, SpatialAudioError> {
    let separation = core::array::from_fn(|j| receiver[j] - source[j]);
    let radius = norm(separation) / sound_step;
    validate_finite(radius, "uniform-flow source-listener separation", 0)?;
    if radius == 0.0 || flow_step == [0.0; 3] {
        return Ok(radius);
    }
    let direction = separation.map(|v| v / (radius * sound_step));
    let mach = flow_step.map(|v| v / sound_step);
    let margin = 1.0 - dot(mach, mach);
    let along = dot(direction, mach);
    let root = det::sqrt(along * along + margin);
    let delay = if along >= 0.0 {
        radius / (root + along)
    } else {
        radius * (root - along) / margin
    };
    let residual = delay - norm(fluid_separation(source, receiver, delay, flow_step)) / sound_step;
    if !delay.is_finite() || delay < 0.0 || !residual.is_finite() || residual.abs() > 1.0e-9 {
        return Err(invalid("representable uniform-flow propagation delay"));
    }
    Ok(delay)
}

fn fluid_separation(
    source: [f64; 3],
    receiver: [f64; 3],
    elapsed: f64,
    flow_step: [f64; 3],
) -> [f64; 3] {
    core::array::from_fn(|j| receiver[j] - source[j] - flow_step[j] * elapsed)
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

    /// G1/G3: a straight span has an analytic radial extremum, including
    /// flybys and a collinear crossing. Check both independent clock laws.
    #[test]
    fn fixed_endpoint_directional_bound_covers_complete_linear_spans() {
        for (start, end, sound_step) in [
            ([4.0, 0.0, 0.0], [3.25, 0.0, 0.0], 1.0),
            ([4.0, 0.0, 0.0], [4.75, 0.0, 0.0], 1.0),
            ([2.0, 1.0, 0.0], [-2.0, 1.0, 0.0], 8.0),
            ([1.0, 0.0, 0.0], [-1.0, 0.0, 0.0], 4.0),
            ([4.0, 3.0, 2.0], [4.125, 3.25, 2.375], 1.0),
        ] {
            let beta: [f64; 3] =
                core::array::from_fn(|axis| (end[axis] - start[axis]) / sound_step);
            for moving_source in [true, false] {
                let bound =
                    fixed_endpoint_segment_rate(start, end, [0.0; 3], sound_step, moving_source)
                        .unwrap();
                for step in 0..=128 {
                    let fraction = f64::from(step) / 128.0;
                    let ray: [f64; 3] = core::array::from_fn(|axis| {
                        -(start[axis] + fraction * (end[axis] - start[axis]))
                    });
                    let distance = ray.iter().map(|x| x * x).sum::<f64>().sqrt();
                    if distance == 0.0 {
                        continue;
                    }
                    let radial = ray
                        .iter()
                        .zip(beta)
                        .map(|(r, v)| r * v / distance)
                        .sum::<f64>();
                    let exact = if moving_source {
                        1.0 / (1.0 - radial)
                    } else {
                        1.0 + radial
                    };
                    assert!(exact <= bound, "{moving_source}: {exact} > {bound}");
                    if step == 0 {
                        assert!(bound - exact < 1.0e-12, "resolved first-event extremum");
                    }
                }
            }
        }
    }

    #[test]
    fn fixed_endpoint_refinement_retains_coincidence_fallback_and_cancellation() {
        for moving_source in [true, false] {
            assert!(
                fixed_endpoint_segment_rate(
                    [0.0; 3],
                    [0.25, 0.0, 0.0],
                    [0.0; 3],
                    1.0,
                    moving_source,
                )
                .is_none()
            );
            let bound = fixed_endpoint_maximum_rate(
                2,
                |i| [i as f64 * 0.25, 0.0, 0.0],
                [0.0; 3],
                1.0,
                moving_source,
                4.0,
                &mut || Ok(()),
            )
            .unwrap();
            assert_eq!(bound, 4.0);
        }
        let mut polls = 0;
        let result = fixed_endpoint_maximum_rate(
            512,
            |i| [1.0 + i as f64 * 0.25, 0.0, 0.0],
            [0.0; 3],
            1.0,
            true,
            4.0,
            &mut || {
                polls += 1;
                Err(SpatialAudioError::Cancelled)
            },
        );
        assert!(matches!(result, Err(SpatialAudioError::Cancelled)));
        assert_eq!(polls, 1);
    }
}
