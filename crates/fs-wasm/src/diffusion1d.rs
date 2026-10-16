//! Explicit FTCS 1-D diffusion for BM-06 (`diffusion1d_frames`).
//!
//! Operator `L` is the mirrored-ghost zero-flux discretisation of `+Δ`
//! without `1/dx²`, assembled through `fs-sparse`. Sign is opposite
//! `laplacian_5pt` (`-Δ`). Time update is unfused so it can agree bit
//! for bit with the TypeScript owner in
//! `src/physics/reference/diffusion/ftcs.ts`.
//!
//! Native API: typed [`Refusal`]. An empty buffer is never a refusal.

use fs_sparse::{Coo, Csr};

/// Declared output budget: `frames * n` per call.
pub const DIFFUSION1D_MAX_OUTPUT_LEN: usize = 2_097_152;
/// Declared step budget: `frames * steps_per_frame` per call.
pub const DIFFUSION1D_MAX_TOTAL_STEPS: usize = 1_048_576;
/// Envelope identity for ok-records.
pub const KERNEL_VERSION: &str = "fs-wasm 0.0.1 diffusion1d_frames";

/// Typed refusal or execution-outcome envelope.
#[derive(Debug, Clone, PartialEq)]
pub struct Refusal {
    pub code: &'static str,
    pub message: String,
    pub ranked_repairs: Vec<String>,
    pub details: String,
}

impl Refusal {
    fn new(
        code: &'static str,
        message: String,
        ranked_repairs: Vec<String>,
        details: String,
    ) -> Self {
        Self {
            code,
            message,
            ranked_repairs,
            details,
        }
    }

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

fn json_f64(x: f64) -> String {
    if x.is_nan() {
        return "null".into();
    }
    if x.is_infinite() {
        return if x.is_sign_positive() {
            "null".into()
        } else {
            "null".into()
        };
    }
    // Debug keeps enough digits to round-trip the planted 0.5000001 case.
    format!("{x:?}")
}

/// Private admitted specification. Fields are crate-private so a forged
/// spec cannot be built outside this module.
#[derive(Debug, Clone)]
pub struct Diffusion1dSpec {
    n: usize,
    frames: usize,
    steps_per_frame: usize,
    #[allow(dead_code)]
    diffusion: f64,
    dx: f64,
    #[allow(dead_code)]
    dt: f64,
    profile: u32,
    r: f64,
    operator: Csr,
}

impl Diffusion1dSpec {
    #[must_use]
    pub fn n(&self) -> usize {
        self.n
    }
    #[must_use]
    pub fn frames(&self) -> usize {
        self.frames
    }
    #[must_use]
    pub fn steps_per_frame(&self) -> usize {
        self.steps_per_frame
    }
    #[must_use]
    pub fn r(&self) -> f64 {
        self.r
    }
    #[must_use]
    pub fn profile(&self) -> u32 {
        self.profile
    }
}

/// Documented ratio: `(diffusion * dt) / (dx * dx)` in that order.
#[must_use]
pub fn stability_ratio(diffusion: f64, dt: f64, dx: f64) -> f64 {
    (diffusion * dt) / (dx * dx)
}

/// Zero-flux `+Δ` stencil without `1/dx²`. Opposite sign of `laplacian_5pt`.
#[must_use]
pub fn assemble_zero_flux_laplacian(n: usize) -> Csr {
    let mut coo = Coo::new(n, n);
    coo.push(0, 0, -1.0);
    coo.push(0, 1, 1.0);
    for i in 1..n - 1 {
        coo.push(i, i - 1, 1.0);
        coo.push(i, i, -2.0);
        coo.push(i, i + 1, 1.0);
    }
    coo.push(n - 1, n - 2, 1.0);
    coo.push(n - 1, n - 1, -1.0);
    coo.assemble()
}

fn adjacent_positive(value: f64, upward: bool) -> f64 {
    let bits = value.to_bits();
    if upward {
        f64::from_bits(bits + 1)
    } else {
        f64::from_bits(bits - 1)
    }
}

fn stable_repair(mut value: f64, upward: bool, ratio: impl Fn(f64) -> f64) -> Option<f64> {
    for _ in 0..4 {
        if !(value > 0.0) || !value.is_finite() {
            return None;
        }
        let r = ratio(value);
        if r > 0.0 && r <= 0.5 {
            return Some(value);
        }
        value = adjacent_positive(value, upward);
    }
    None
}

fn fill_profile(n: usize, dx: f64, profile: u32) -> Result<Vec<f64>, Refusal> {
    let mut field = vec![0.0; n];
    match profile {
        0 => {
            let ic = n / 2;
            field[ic] = 1.0 / dx;
        }
        1 => {
            let mid = n / 2;
            for v in field.iter_mut().take(mid) {
                *v = 1.0;
            }
        }
        2 => {
            let a = n / 4;
            let b = (3 * n) / 4;
            field[a] = 0.5 / dx;
            field[b] = 0.5 / dx;
        }
        _ => unreachable!("admitted profiles only"),
    }
    if !field.iter().all(|v| v.is_finite()) {
        return Err(Refusal::new(
            "nonfinite-input",
            "The initial cell density is outside binary64 range; the grid was not advanced.".into(),
            vec!["choose a representable dx".into()],
            "{\"reason\":\"initial-profile\"}".into(),
        ));
    }
    Ok(field)
}

fn refuse_invalid(parameter: &str, message: &str, details: String) -> Refusal {
    Refusal::new(
        "invalid-parameter",
        message.into(),
        vec![format!("repair {parameter}")],
        details,
    )
}

/// Admit inputs. Never clamps.
pub fn admit_diffusion1d_frames(
    n: usize,
    frames: usize,
    steps_per_frame: usize,
    diffusion: f64,
    dx: f64,
    dt: f64,
    profile: u32,
) -> Result<Diffusion1dSpec, Refusal> {
    if n < 3 {
        return Err(refuse_invalid(
            "n",
            "Use at least three cells.",
            format!("{{\"parameterIds\":[\"n\"],\"n\":{n}}}"),
        ));
    }
    if frames == 0 || steps_per_frame == 0 {
        return Err(refuse_invalid(
            "frames",
            "Use positive whole-number frame and step counts.",
            format!(
                "{{\"parameterIds\":[\"frames\",\"steps_per_frame\"],\"frames\":{frames},\"stepsPerFrame\":{steps_per_frame}}}"
            ),
        ));
    }
    if profile > 2 {
        return Err(Refusal::new(
            "unsupported-kernel",
            format!("profile {profile} is not in {{0,1,2}}"),
            vec!["use profile 0 (spike), 1 (step), or 2 (two spikes)".into()],
            format!("{{\"profile\":{profile}}}"),
        ));
    }
    if !diffusion.is_finite() || !dx.is_finite() || !dt.is_finite() {
        let name = if !diffusion.is_finite() {
            "diffusion"
        } else if !dx.is_finite() {
            "dx"
        } else {
            "dt"
        };
        return Err(Refusal::new(
            "nonfinite-input",
            format!("{name} must be finite"),
            vec!["pass finite diffusion (>= 0) and strictly positive dx and dt".into()],
            format!("{{\"parameterIds\":[\"diffusion\",\"dx\",\"dt\"],\"name\":\"{name}\"}}"),
        ));
    }
    if diffusion < 0.0 || dx <= 0.0 || dt <= 0.0 {
        return Err(refuse_invalid(
            "dt",
            "Use nonnegative diffusion and strictly positive dx and dt.",
            format!(
                "{{\"parameterIds\":[\"diffusion\",\"dx\",\"dt\"],\"diffusion\":{},\"dx\":{},\"dt\":{}}}",
                json_f64(diffusion),
                json_f64(dx),
                json_f64(dt)
            ),
        ));
    }
    let out_len = frames.checked_mul(n).ok_or_else(|| {
        Refusal::new(
            "budget-exhausted",
            "frames * n overflows usize.".into(),
            vec!["reduce frames or n".into()],
            format!("{{\"frames\":{frames},\"n\":{n}}}"),
        )
    })?;
    if out_len > DIFFUSION1D_MAX_OUTPUT_LEN {
        return Err(Refusal::new(
            "budget-exhausted",
            format!("frames * n = {out_len} exceeds DIFFUSION1D_MAX_OUTPUT_LEN"),
            vec!["reduce frames or n".into()],
            format!(
                "{{\"requested\":{out_len},\"allowed\":{DIFFUSION1D_MAX_OUTPUT_LEN},\"frames\":{frames},\"n\":{n}}}"
            ),
        ));
    }
    let total_steps = frames.checked_mul(steps_per_frame).ok_or_else(|| {
        Refusal::new(
            "budget-exhausted",
            "frames * steps_per_frame overflows usize.".into(),
            vec!["reduce frames or steps_per_frame".into()],
            format!("{{\"frames\":{frames},\"stepsPerFrame\":{steps_per_frame}}}"),
        )
    })?;
    if total_steps > DIFFUSION1D_MAX_TOTAL_STEPS {
        return Err(Refusal::new(
            "budget-exhausted",
            format!("frames * steps_per_frame = {total_steps} exceeds DIFFUSION1D_MAX_TOTAL_STEPS"),
            vec!["reduce frames or steps_per_frame".into()],
            format!("{{\"requested\":{total_steps},\"allowed\":{DIFFUSION1D_MAX_TOTAL_STEPS}}}"),
        ));
    }

    let r = if diffusion == 0.0 {
        0.0
    } else {
        stability_ratio(diffusion, dt, dx)
    };
    if !r.is_finite() || (diffusion > 0.0 && r == 0.0) {
        return Err(Refusal::new(
            "nonfinite-input",
            "The stability ratio is outside binary64 range; the grid was not advanced.".into(),
            vec!["rescale dx, dt, or diffusion into the representable range".into()],
            format!(
                "{{\"ratio\":{},\"diffusion\":{},\"dx\":{},\"dt\":{}}}",
                json_f64(r),
                json_f64(diffusion),
                json_f64(dx),
                json_f64(dt)
            ),
        ));
    }
    if r > 0.5 {
        let dt_max = (dx * dx) / (2.0 * diffusion);
        let dx_min = (2.0 * diffusion * dt).sqrt();
        let d_max = (dx * dx) / (2.0 * dt);
        if ![dt_max, dx_min, d_max]
            .iter()
            .all(|v| v.is_finite() && *v > 0.0)
        {
            return Err(Refusal::new(
                "nonfinite-input",
                "The dimensional stability repairs are outside binary64 range.".into(),
                vec!["rescale the grid".into()],
                format!("{{\"ratio\":{},\"limit\":0.5}}", json_f64(r)),
            ));
        }
        let dt_repair = stable_repair(dt_max, false, |v| (diffusion * v) / (dx * dx));
        let dx_repair = stable_repair(dx_min, true, |v| (diffusion * dt) / (v * v));
        let d_repair = stable_repair(d_max, false, |v| (v * dt) / (dx * dx));
        let (Some(dt_r), Some(dx_r), Some(d_r)) = (dt_repair, dx_repair, d_repair) else {
            return Err(Refusal::new(
                "nonfinite-input",
                "No representable repair was found near the dimensional stability boundary.".into(),
                vec!["choose a coarser grid by hand".into()],
                format!(
                    "{{\"ratio\":{},\"limit\":0.5,\"dtMax\":{}}}",
                    json_f64(r),
                    json_f64(dt_max)
                ),
            ));
        };
        return Err(Refusal::new(
            "ftcs-unstable",
            "This time step is too large for the explicit diffusion scheme.".into(),
            vec![
                format!(
                    "Reduce the time step to the explicit scheme's limit (dt = {}).",
                    json_f64(dt_r)
                ),
                format!("Use a coarser spatial grid (dx = {}).", json_f64(dx_r)),
                format!(
                    "Choose a smaller diffusivity; this changes the physical setup (diffusion = {}).",
                    json_f64(d_r)
                ),
            ],
            format!(
                "{{\"ratio\":{},\"limit\":0.5,\"dtMax\":{},\"repairs\":[{{\"parameterId\":\"dt\",\"value\":{}}},{{\"parameterId\":\"dx\",\"value\":{}}},{{\"parameterId\":\"diffusion\",\"value\":{}}}]}}",
                json_f64(r),
                json_f64(dt_max),
                json_f64(dt_r),
                json_f64(dx_r),
                json_f64(d_r)
            ),
        ));
    }

    Ok(Diffusion1dSpec {
        n,
        frames,
        steps_per_frame,
        diffusion,
        dx,
        dt,
        profile,
        r,
        operator: assemble_zero_flux_laplacian(n),
    })
}

fn advance(
    field: &mut [f64],
    y: &mut [f64],
    operator: &Csr,
    r: f64,
    steps: usize,
) -> Result<(), Refusal> {
    if r == 0.0 {
        return Ok(());
    }
    for _ in 0..steps {
        operator.spmv(field, y);
        for i in 0..field.len() {
            let increment = r * y[i];
            let next = field[i] + increment;
            if !next.is_finite() || next < 0.0 {
                return Err(Refusal::new(
                    "nonfinite-input",
                    "A cell became nonfinite or negative; no partial field is returned.".into(),
                    vec!["reduce r or the number of steps".into()],
                    "{\"reason\":\"step-blew-up\"}".into(),
                ));
            }
            field[i] = next;
        }
    }
    Ok(())
}

fn fill_admitted(spec: &Diffusion1dSpec) -> Result<Vec<f64>, Refusal> {
    let mut field = fill_profile(spec.n, spec.dx, spec.profile)?;
    let mut y = vec![0.0; spec.n];
    let mut out = vec![0.0; spec.frames * spec.n];
    out[..spec.n].copy_from_slice(&field);
    for frame in 1..spec.frames {
        advance(
            &mut field,
            &mut y,
            &spec.operator,
            spec.r,
            spec.steps_per_frame,
        )?;
        let start = frame * spec.n;
        out[start..start + spec.n].copy_from_slice(&field);
    }
    Ok(out)
}

#[must_use]
pub fn diffusion1d_frames_admitted(spec: &Diffusion1dSpec) -> Result<Vec<f64>, Refusal> {
    fill_admitted(spec)
}

pub fn diffusion1d_frames(
    n: usize,
    frames: usize,
    steps_per_frame: usize,
    diffusion: f64,
    dx: f64,
    dt: f64,
    profile: u32,
) -> Result<Vec<f64>, Refusal> {
    let spec = admit_diffusion1d_frames(n, frames, steps_per_frame, diffusion, dx, dt, profile)?;
    fill_admitted(&spec)
}

#[must_use]
pub fn ok_envelope_json(spec: &Diffusion1dSpec, value_count: usize) -> String {
    format!(
        "{{\"ok\":{{\"kernel\":\"{}\",\"export\":\"diffusion1d_frames\",\"layout\":{{\"n\":{},\"frames\":{},\"index\":\"f*n+i\"}},\"quantityId\":\"probabilityDensity\",\"profile\":{},\"r\":{},\"valueCount\":{value_count}}}}}",
        json_escape(KERNEL_VERSION),
        spec.n,
        spec.frames,
        spec.profile,
        json_f64(spec.r)
    )
}
