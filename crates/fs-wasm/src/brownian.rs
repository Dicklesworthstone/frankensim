//! One-shot and windowed 1-D Brownian frames for BM-01 / BM-05.
//!
//! Step-kernel ids, arithmetic, and quantity bindings are frozen in
//! `docs/FRANKENSIM_BINDING.md` (`kernel-resolution`, `distinct-meaning`).
//! Kernel 2 is the dimensionless unit-Gaussian teaching walk; kernel 3 is
//! the physically scaled Gaussian with per-step variance `2 D dt`. They are
//! not aliases. The TypeScript owner (`walkLaws.ts`) pins fourth-moment
//! factors 1 (coin), 1.8 (uniform), 3 (Gaussian); this export uses those
//! same distributions.
//!
//! Native API: typed [`Refusal`], never `NaN` and never an unexplained empty
//! buffer. The one-shot form is one window from zeros at `start_step = 0`.

use fs_math::det;
use fs_rand::{
    STREAM_CHECKPOINT_CANONICAL_LEN, STREAM_CHECKPOINT_VERSION, STREAM_SEMANTICS_VERSION, Stream,
    StreamCheckpoint, StreamKey, StreamReplayError,
};

/// Declared output budget: `n_particles * (steps + 1)` per call.
pub const BROWNIAN_MAX_OUTPUT_LEN: usize = 2_097_152;
/// Registered Brownian latent stream kernel id (`0x19050001`).
pub const BROWNIAN_STREAM_KERNEL_ID: u32 = 0x1905_0001;
/// Envelope identity for ok-records.
pub const KERNEL_VERSION: &str = "fs-wasm 0.0.1 brownian_frames";

/// Typed refusal or execution-outcome envelope. An empty buffer is never a
/// substitute for this value.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Refusal {
    /// Registered target (`unsupported-kernel`, `invalid-parameter`, …) or
    /// the execution outcome `budget-exhausted`.
    pub code: &'static str,
    /// Readable reason.
    pub message: String,
    /// Ordered repairs a host may offer.
    pub ranked_repairs: Vec<&'static str>,
    /// JSON object, `{}` if none.
    pub details: String,
}

impl Refusal {
    fn new(
        code: &'static str,
        message: String,
        ranked_repairs: Vec<&'static str>,
        details: String,
    ) -> Self {
        Self {
            code,
            message,
            ranked_repairs,
            details,
        }
    }

    /// `{"refusal":{...}}` for the WASM / JS boundary.
    #[must_use]
    pub fn to_json(&self) -> String {
        let repairs: Vec<String> = self
            .ranked_repairs
            .iter()
            .map(|r| format!("\"{}\"", json_escape(r)))
            .collect();
        format!(
            "{{\"refusal\":{{\"code\":\"{}\",\"message\":\"{}\",\"ranked_repairs\":[{}],\"details\":{}}}}}",
            json_escape(self.code),
            json_escape(&self.message),
            repairs.join(","),
            self.details
        )
    }
}

fn json_escape(s: &str) -> String {
    let mut out = String::with_capacity(s.len());
    for ch in s.chars() {
        match ch {
            '\\' => out.push_str("\\\\"),
            '"' => out.push_str("\\\""),
            '\n' => out.push_str("\\n"),
            '\r' => out.push_str("\\r"),
            '\t' => out.push_str("\\t"),
            c if (c as u32) < 0x20 => out.push_str(&format!("\\u{:04x}", c as u32)),
            c => out.push(c),
        }
    }
    out
}

/// Draws consumed per step. `None` for an unknown step kernel.
#[must_use]
pub fn draws_per_step(step_kernel: u32) -> Option<u32> {
    match step_kernel {
        0 | 1 => Some(1),
        2 | 3 => Some(2),
        _ => None,
    }
}

/// Quantity id and unit bound by an admitted kernel.
#[must_use]
pub fn quantity_binding(step_kernel: u32) -> Option<(&'static str, &'static str)> {
    match step_kernel {
        0 | 1 | 3 => Some(("latentPosition1d", "metre")),
        2 => Some(("walkStepCoordinate1d", "step")),
        _ => None,
    }
}

#[derive(Debug, Clone, Copy)]
struct Scales {
    s: f64,
    h: f64,
}

/// Private admitted one-shot / window specification.
#[derive(Debug, Clone)]
pub struct BrownianFramesSpec {
    n_particles: usize,
    start_step: usize,
    steps: usize,
    step_kernel: u32,
    seed: u64,
    scales: Scales,
    start_positions: Vec<f64>,
}

impl BrownianFramesSpec {
    #[must_use]
    pub fn n_particles(&self) -> usize {
        self.n_particles
    }
    #[must_use]
    pub fn steps(&self) -> usize {
        self.steps
    }
    #[must_use]
    pub fn step_kernel(&self) -> u32 {
        self.step_kernel
    }
}

fn stream_for(seed: u64, tile: u32, index: u64) -> Stream {
    Stream::resume(StreamCheckpoint::current(
        StreamKey {
            seed,
            kernel: BROWNIAN_STREAM_KERNEL_ID,
            tile,
        },
        index,
    ))
    .expect("current() checkpoints are this build's versions")
}

fn refuse_unsupported(step_kernel: u32) -> Refusal {
    Refusal::new(
        "unsupported-kernel",
        format!("step_kernel {step_kernel} is not in {{0,1,2,3}}"),
        vec![
            "use step_kernel 0 (coin), 1 (uniform), 2 (unit Gaussian teaching), or 3 (Gaussian with exact D)",
        ],
        format!("{{\"stepKernel\":{step_kernel}}}"),
    )
}

fn admit_common(
    n_particles: usize,
    start_step: usize,
    steps: usize,
    step_kernel: u32,
    seed: u64,
    diffusion: f64,
    dt: f64,
    start_positions: Option<&[f64]>,
) -> Result<BrownianFramesSpec, Refusal> {
    let Some(dps) = draws_per_step(step_kernel) else {
        return Err(refuse_unsupported(step_kernel));
    };
    if !diffusion.is_finite() || !dt.is_finite() {
        let name = if !diffusion.is_finite() {
            "diffusion"
        } else {
            "dt"
        };
        let value = if !diffusion.is_finite() {
            diffusion
        } else {
            dt
        };
        return Err(Refusal::new(
            "nonfinite-input",
            format!("{name} must be finite"),
            vec!["pass a finite diffusion (>= 0) and a finite dt (> 0)"],
            format!("{{\"name\":\"{name}\",\"value\":\"{value:?}\"}}"),
        ));
    }
    if diffusion < 0.0 || dt <= 0.0 {
        let name = if diffusion < 0.0 { "diffusion" } else { "dt" };
        let value = if diffusion < 0.0 { diffusion } else { dt };
        return Err(Refusal::new(
            "invalid-parameter",
            format!("{name} is outside the admitted domain"),
            vec![
                "use diffusion >= 0; use dt > 0",
                "diffusion = 0 is valid and returns zeros",
            ],
            format!("{{\"name\":\"{name}\",\"value\":\"{value:?}\"}}"),
        ));
    }
    if n_particles == 0 {
        return Err(Refusal::new(
            "invalid-parameter",
            "n_particles must be at least 1".into(),
            vec!["use n_particles >= 1 and steps >= 1"],
            "{\"name\":\"n_particles\",\"value\":0}".into(),
        ));
    }
    if n_particles > u32::MAX as usize {
        return Err(Refusal::new(
            "invalid-parameter",
            "n_particles exceeds StreamKey.tile width".into(),
            vec!["n_particles must fit StreamKey.tile (u32)"],
            format!("{{\"name\":\"n_particles\",\"limit\":{}}}", u32::MAX),
        ));
    }
    if steps == 0 {
        let code_msg = if start_positions.is_some() {
            ("invalid-parameter", "window steps must be >= 1")
        } else {
            (
                "invalid-parameter",
                "n_particles must be at least 1 and steps >= 1",
            )
        };
        return Err(Refusal::new(
            code_msg.0,
            code_msg.1.into(),
            if start_positions.is_some() {
                vec!["window steps must be >= 1"]
            } else {
                vec!["use n_particles >= 1 and steps >= 1"]
            },
            "{\"name\":\"steps\",\"value\":0}".into(),
        ));
    }
    let end_step = start_step.checked_add(steps).ok_or_else(|| {
        Refusal::new(
            "budget-exhausted",
            "start_step + steps overflows usize".into(),
            vec![
                "use brownian_frames_window with a smaller steps",
                "reduce n_particles",
            ],
            format!(
                "{{\"requested\":\"overflow\",\"allowed\":{BROWNIAN_MAX_OUTPUT_LEN},\"unit\":\"f64-values\"}}"
            ),
        )
    })?;
    let len = n_particles.checked_mul(steps.checked_add(1).ok_or_else(|| {
        Refusal::new(
            "budget-exhausted",
            "steps + 1 overflows usize".into(),
            vec![
                "use brownian_frames_window with a smaller steps",
                "reduce n_particles",
            ],
            format!(
                "{{\"requested\":\"overflow\",\"allowed\":{BROWNIAN_MAX_OUTPUT_LEN},\"unit\":\"f64-values\"}}"
            ),
        )
    })?)
    .ok_or_else(|| {
        Refusal::new(
            "budget-exhausted",
            "n_particles * (steps + 1) overflows usize".into(),
            vec![
                "use brownian_frames_window with a smaller steps",
                "reduce n_particles",
            ],
            format!(
                "{{\"requested\":\"overflow\",\"allowed\":{BROWNIAN_MAX_OUTPUT_LEN},\"unit\":\"f64-values\"}}"
            ),
        )
    })?;
    if len > BROWNIAN_MAX_OUTPUT_LEN {
        return Err(Refusal::new(
            "budget-exhausted",
            format!("output length {len} exceeds {BROWNIAN_MAX_OUTPUT_LEN}"),
            vec![
                "use brownian_frames_window with a smaller steps",
                "reduce n_particles",
            ],
            format!(
                "{{\"requested\":{len},\"allowed\":{BROWNIAN_MAX_OUTPUT_LEN},\"unit\":\"f64-values\"}}"
            ),
        ));
    }
    let dps = u64::from(dps);
    let start_index = dps
        .checked_mul(start_step as u64)
        .ok_or_else(|| stream_overflow(start_step, steps, dps))?;
    let total_draws = dps
        .checked_mul(end_step as u64)
        .ok_or_else(|| stream_overflow(start_step, steps, dps))?;
    let _ = (start_index, total_draws);
    let starts = match start_positions {
        None => vec![0.0; n_particles],
        Some(buf) => {
            if buf.len() != n_particles {
                return Err(Refusal::new(
                    "invalid-parameter",
                    format!(
                        "start_positions length {} != n_particles {n_particles}",
                        buf.len()
                    ),
                    vec![
                        "pass n_particles finite positions",
                        "start_step 0 requires all zeros",
                    ],
                    "{\"name\":\"start_positions\",\"index\":\"length\"}".into(),
                ));
            }
            for (i, &x) in buf.iter().enumerate() {
                if !x.is_finite() {
                    return Err(Refusal::new(
                        "invalid-parameter",
                        format!("start_positions[{i}] is not finite"),
                        vec![
                            "pass n_particles finite positions",
                            "start_step 0 requires all zeros",
                        ],
                        format!("{{\"name\":\"start_positions\",\"index\":{i}}}"),
                    ));
                }
                if start_step == 0 && x != 0.0 {
                    return Err(Refusal::new(
                        "invalid-parameter",
                        format!("start_step 0 requires start_positions[{i}] == 0.0"),
                        vec![
                            "pass n_particles finite positions",
                            "start_step 0 requires all zeros",
                        ],
                        format!("{{\"name\":\"start_positions\",\"index\":{i}}}"),
                    ));
                }
            }
            buf.to_vec()
        }
    };
    Ok(BrownianFramesSpec {
        n_particles,
        start_step,
        steps,
        step_kernel,
        seed,
        scales: Scales {
            s: det::sqrt(2.0 * diffusion * dt),
            h: det::sqrt(6.0 * diffusion * dt),
        },
        start_positions: starts,
    })
}

fn stream_overflow(start_step: usize, steps: usize, dps: u64) -> Refusal {
    let start_index = dps.saturating_mul(start_step as u64);
    let draws = dps.saturating_mul(steps as u64);
    Refusal::new(
        "stream-index-overflow",
        "draws_per_step * (start_step + steps) exceeds 2^64-1".into(),
        vec![
            "reduce start_step + steps",
            "use a coin or uniform kernel (1 draw/step) if the Gaussian 2-draw counter is the limiter",
        ],
        format!(
            "{{\"startIndex\":\"{start_index}\",\"draws\":\"{draws}\",\"maxIndex\":\"18446744073709551615\"}}"
        ),
    )
}

/// Admit a one-shot call. Does not allocate the position buffer.
pub fn admit_brownian_frames(
    n_particles: usize,
    steps: usize,
    step_kernel: u32,
    seed: u64,
    diffusion: f64,
    dt: f64,
) -> Result<BrownianFramesSpec, Refusal> {
    admit_common(
        n_particles,
        0,
        steps,
        step_kernel,
        seed,
        diffusion,
        dt,
        None,
    )
}

/// Admit a windowed continuation.
pub fn admit_brownian_frames_window(
    n_particles: usize,
    start_step: usize,
    steps: usize,
    step_kernel: u32,
    seed: u64,
    diffusion: f64,
    dt: f64,
    start_positions: &[f64],
) -> Result<BrownianFramesSpec, Refusal> {
    admit_common(
        n_particles,
        start_step,
        steps,
        step_kernel,
        seed,
        diffusion,
        dt,
        Some(start_positions),
    )
}

fn step_from(stream: &mut Stream, kernel: u32, scales: Scales) -> f64 {
    match kernel {
        0 => {
            let u = stream.next_u64();
            if u >> 63 == 1 { scales.s } else { -scales.s }
        }
        1 => {
            let u = stream.next_f64();
            (2.0 * u - 1.0) * scales.h
        }
        2 => stream.next_normal(),
        3 => stream.next_normal() * scales.s,
        _ => unreachable!("admitted kernels only"),
    }
}

fn fill_admitted(spec: &BrownianFramesSpec) -> Vec<f64> {
    let cols = spec.steps + 1;
    let mut out = vec![0.0; spec.n_particles * cols];
    let dps = u64::from(draws_per_step(spec.step_kernel).expect("admitted"));
    let start_index = dps * spec.start_step as u64;
    for p in 0..spec.n_particles {
        let mut stream = stream_for(spec.seed, p as u32, start_index);
        let row = p * cols;
        let mut x = spec.start_positions[p];
        out[row] = x;
        for s in 1..=spec.steps {
            x += step_from(&mut stream, spec.step_kernel, spec.scales);
            out[row + s] = x;
        }
        debug_assert_eq!(stream.index(), start_index + dps * spec.steps as u64);
    }
    out
}

/// Fill an admitted spec. Fields are crate-private, so a forged spec cannot
/// be constructed outside this module.
#[must_use]
pub fn brownian_frames_admitted(spec: &BrownianFramesSpec) -> Vec<f64> {
    fill_admitted(spec)
}

/// One-shot: `admit` then fill. Layout `p * (steps + 1) + s`, start at 0.
pub fn brownian_frames(
    n_particles: usize,
    steps: usize,
    step_kernel: u32,
    seed: u64,
    diffusion: f64,
    dt: f64,
) -> Result<Vec<f64>, Refusal> {
    let spec = admit_brownian_frames(n_particles, steps, step_kernel, seed, diffusion, dt)?;
    Ok(fill_admitted(&spec))
}

/// Windowed continuation. First column equals `start_positions`.
pub fn brownian_frames_window(
    n_particles: usize,
    start_step: usize,
    steps: usize,
    step_kernel: u32,
    seed: u64,
    diffusion: f64,
    dt: f64,
    start_positions: &[f64],
) -> Result<Vec<f64>, Refusal> {
    let spec = admit_brownian_frames_window(
        n_particles,
        start_step,
        steps,
        step_kernel,
        seed,
        diffusion,
        dt,
        start_positions,
    )?;
    Ok(fill_admitted(&spec))
}

/// Validate a retained checkpoint against the derived window identity.
pub fn admit_brownian_checkpoint(
    seed: u64,
    step_kernel: u32,
    start_step: usize,
    particle: usize,
    bytes: &[u8],
) -> Result<(), Refusal> {
    let Some(dps) = draws_per_step(step_kernel) else {
        return Err(refuse_unsupported(step_kernel));
    };
    if particle > u32::MAX as usize {
        return Err(Refusal::new(
            "invalid-parameter",
            "particle exceeds StreamKey.tile width".into(),
            vec!["n_particles must fit StreamKey.tile (u32)"],
            format!("{{\"name\":\"particle\",\"limit\":{}}}", u32::MAX),
        ));
    }
    let expected_index = u64::from(dps)
        .checked_mul(start_step as u64)
        .ok_or_else(|| stream_overflow(start_step, 0, u64::from(dps)))?;
    let expected = StreamKey {
        seed,
        kernel: BROWNIAN_STREAM_KERNEL_ID,
        tile: particle as u32,
    };
    let checkpoint = match StreamCheckpoint::from_canonical_le_bytes(bytes) {
        Ok(c) => c,
        Err(StreamReplayError::InvalidCheckpointLength { actual, expected }) => {
            return Err(Refusal::new(
                "invalid-parameter",
                format!("stream checkpoint has {actual} bytes; expected {expected}"),
                vec!["supply an 83-byte canonical StreamCheckpoint frame"],
                format!(
                    "{{\"reason\":\"checkpoint-disagreement\",\"field\":\"length\",\"declared\":{actual},\"expected\":{expected}}}"
                ),
            ));
        }
        Err(StreamReplayError::InvalidCheckpointMagic) => {
            return Err(Refusal::new(
                "invalid-parameter",
                "stream checkpoint magic is invalid".into(),
                vec!["frame must begin with FSRCKPT\\0"],
                "{\"reason\":\"checkpoint-disagreement\",\"field\":\"magic\"}".into(),
            ));
        }
        Err(StreamReplayError::InvalidCheckpointDomain) => {
            return Err(Refusal::new(
                "invalid-parameter",
                "stream checkpoint belongs to another identity domain".into(),
                vec!["use the fs-rand stream-checkpoint identity domain"],
                "{\"reason\":\"checkpoint-disagreement\",\"field\":\"domain\"}".into(),
            ));
        }
        Err(StreamReplayError::UnknownCheckpointVersion {
            declared,
            supported,
        }) => {
            return Err(Refusal::new(
                "invalid-parameter",
                format!("checkpoint_version {declared} is not {supported}"),
                vec!["rebuild the checkpoint under STREAM_CHECKPOINT_VERSION = 1"],
                format!(
                    "{{\"reason\":\"checkpoint-disagreement\",\"field\":\"checkpoint_version\",\"declared\":{declared},\"expected\":{supported}}}"
                ),
            ));
        }
        Err(StreamReplayError::UnknownStreamSemanticsVersion {
            declared,
            supported,
        }) => {
            return Err(Refusal::new(
                "invalid-parameter",
                format!("stream_semantics_version {declared} is not {supported}"),
                vec!["rebuild the checkpoint under STREAM_SEMANTICS_VERSION = 1"],
                format!(
                    "{{\"reason\":\"checkpoint-disagreement\",\"field\":\"stream_semantics_version\",\"declared\":{declared},\"expected\":{supported}}}"
                ),
            ));
        }
    };
    if checkpoint.key.kernel != expected.kernel {
        return Err(mismatch("kernel", checkpoint.key.kernel, expected.kernel));
    }
    if checkpoint.key.tile != expected.tile {
        return Err(mismatch("tile", checkpoint.key.tile, expected.tile));
    }
    if checkpoint.key.seed != expected.seed {
        return Err(mismatch_u64("seed", checkpoint.key.seed, expected.seed));
    }
    if checkpoint.index != expected_index {
        return Err(mismatch_u64("index", checkpoint.index, expected_index));
    }
    let _ = (
        STREAM_CHECKPOINT_VERSION,
        STREAM_SEMANTICS_VERSION,
        STREAM_CHECKPOINT_CANONICAL_LEN,
    );
    Ok(())
}

fn mismatch(field: &str, declared: u32, expected: u32) -> Refusal {
    Refusal::new(
        "invalid-parameter",
        format!("checkpoint {field} {declared} disagrees with derived {expected}"),
        vec![
            "do not pass a checkpoint; the window derives the draw index from start_step",
            "or pass a checkpoint whose key and index match the derivation",
        ],
        format!(
            "{{\"reason\":\"checkpoint-disagreement\",\"field\":\"{field}\",\"declared\":{declared},\"expected\":{expected}}}"
        ),
    )
}

fn mismatch_u64(field: &str, declared: u64, expected: u64) -> Refusal {
    Refusal::new(
        "invalid-parameter",
        format!("checkpoint {field} {declared} disagrees with derived {expected}"),
        vec![
            "do not pass a checkpoint; the window derives the draw index from start_step",
            "or pass a checkpoint whose key and index match the derivation",
        ],
        format!(
            "{{\"reason\":\"checkpoint-disagreement\",\"field\":\"{field}\",\"declared\":{declared},\"expected\":{expected}}}"
        ),
    )
}

/// Independent reconstruction of one particle for kernels 0 and 1.
#[must_use]
pub fn reconstruct_particle(
    particle: u32,
    steps: usize,
    step_kernel: u32,
    seed: u64,
    diffusion: f64,
    dt: f64,
) -> Option<Vec<f64>> {
    if !matches!(step_kernel, 0 | 1) {
        return None;
    }
    let spec = admit_brownian_frames(1, steps, step_kernel, seed, diffusion, dt).ok()?;
    let mut stream = stream_for(seed, particle, 0);
    let mut x = 0.0;
    let mut out = vec![0.0; steps + 1];
    for s in 1..=steps {
        x += step_from(&mut stream, step_kernel, spec.scales);
        out[s] = x;
    }
    Some(out)
}

/// Stream index after `steps` from `start_step`.
#[must_use]
pub fn derived_stream_index(step_kernel: u32, start_step: usize, steps: usize) -> Option<u64> {
    let dps = u64::from(draws_per_step(step_kernel)?);
    dps.checked_mul((start_step as u64).checked_add(steps as u64)?)
}

/// Ok-envelope JSON without packing the sample buffer as a JSON number list.
#[must_use]
pub fn ok_envelope_json(spec: &BrownianFramesSpec, value_count: usize) -> String {
    let (quantity, unit) = quantity_binding(spec.step_kernel).expect("admitted");
    format!(
        "{{\"ok\":{{\"kernel\":\"{}\",\"export\":\"brownian_frames\",\"layout\":{{\"nParticles\":{},\"steps\":{},\"index\":\"p*(steps+1)+s\"}},\"quantityId\":\"{quantity}\",\"unit\":\"{unit}\",\"stepKernel\":{},\"valueCount\":{value_count}}}}}",
        json_escape(KERNEL_VERSION),
        spec.n_particles,
        spec.steps,
        spec.step_kernel
    )
}
