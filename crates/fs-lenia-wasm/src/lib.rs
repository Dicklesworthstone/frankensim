//! fs-lenia-wasm — browser (WASM) surface for the CMA-ES explainer site's
//! continuous Lenia field. Layer: L6.
//!
//! The site's TypeScript fallback steps Lenia with a direct O(N²·R²) spatial
//! convolution, which caps the field at 96². This kernel computes the SAME
//! toroidal convolution through a planned radix-2 FFT — O(N² log N), with
//! precomputed twiddle/bit-reversal tables, blocked transposes, and a
//! transpose-sparing convolution path — so a 256² or 512² display field with
//! a proportionally larger ring kernel steps faster than the 96² fallback,
//! on phones included. The growth mapping runs through a lookup table cached
//! across steps, and the bioluminescent colormap is written straight into an
//! in-memory RGBA buffer the page blits with zero per-pixel JS work.
//!
//! Model (identical to the site's TS fallback, up to floating-point error):
//!
//! - Kernel: concentric ring K(r) = exp(-((r/R - r0)/sigma_k)²/2) for r <= R,
//!   normalized to sum 1, on a toroidal grid.
//! - Growth: G(u) = 2·exp(-((u - mu)/sigma)²/2) - 1.
//! - Update: a <- clamp(a + dt·G(K*a), 0, 1).
//! - Metrics: "interface" = fraction of cells strictly inside (0.08, 0.92);
//!   "mass" = mean activation. The fitness evaluator scores
//!   mean_over_steps(interface - 2·|mass - 0.25|), the site's objective.
//!
//! Contracts every entry inherits (fs-flyer-wasm pattern):
//!
//! - **Typed-refusal JSON envelope** for fallible entries; nothing is
//!   silently clamped and nothing traps across the boundary.
//! - **Determinism.** Pure functions of the field state and scalar inputs:
//!   no wall-clock, no entropy. Same call sequence ⇒ identical states.
//!
//! No-claims: teaching/viz surface, not a general spectral PDE solver.

use core::cell::RefCell;

/// Kernel id baked into envelopes so the page can prove which build is live.
pub const KERNEL_VERSION: &str = "fs-lenia-wasm 0.2.0";

const LUT_SIZE: usize = 2048;

// ---------------------------------------------------------------------------
// Refusal envelope
// ---------------------------------------------------------------------------

/// A typed refusal for the JS boundary.
#[derive(Debug, Clone)]
pub struct Refusal {
    pub code: &'static str,
    pub message: String,
    pub ranked_repairs: Vec<&'static str>,
}

fn escape_json_string(s: &str) -> String {
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

impl Refusal {
    /// Render the refusal as the JSON envelope the JS boundary returns.
    pub fn json(&self) -> String {
        let repairs: Vec<String> = self
            .ranked_repairs
            .iter()
            .map(|r| format!("\"{}\"", escape_json_string(r)))
            .collect();
        format!(
            "{{\"refusal\":{{\"code\":\"{}\",\"message\":\"{}\",\"ranked_repairs\":[{}]}}}}",
            self.code,
            escape_json_string(&self.message),
            repairs.join(",")
        )
    }
}

fn require_finite(name: &'static str, v: f64) -> Result<(), Refusal> {
    if v.is_finite() {
        Ok(())
    } else {
        Err(Refusal {
            code: "input-non-finite",
            message: format!("{name} must be finite, got {v}"),
            ranked_repairs: vec!["pass a finite number"],
        })
    }
}

// ---------------------------------------------------------------------------
// Radix-2 FFT with a precomputed plan (bit-reversal + per-stage twiddle
// tables), blocked square transposes, and a transpose-sparing 2D pass.
// ---------------------------------------------------------------------------

/// Precomputed tables for one transform length. Twiddles are stored with the
/// FORWARD sign convention (w_k = exp(-i·2πk/len)); the inverse conjugates on
/// the fly. Tables replace the per-butterfly rotation recurrence, which is
/// both faster and free of the recurrence's accumulated roundoff.
pub struct FftPlan {
    n: usize,
    bitrev: Vec<u32>,
    // Stages len = 2, 4, ..., n concatenated; stage `len` contributes len/2
    // factors, so the tables hold exactly n - 1 entries.
    tw_re: Vec<f64>,
    tw_im: Vec<f64>,
}

impl FftPlan {
    pub fn new(n: usize) -> Self {
        debug_assert!(n.is_power_of_two() && n >= 2);
        let mut bitrev = vec![0u32; n];
        let mut j = 0usize;
        for i in 1..n {
            let mut bit = n >> 1;
            while j & bit != 0 {
                j ^= bit;
                bit >>= 1;
            }
            j |= bit;
            bitrev[i] = j as u32;
        }
        let mut tw_re = Vec::with_capacity(n - 1);
        let mut tw_im = Vec::with_capacity(n - 1);
        let mut len = 2usize;
        while len <= n {
            let half = len / 2;
            for k in 0..half {
                let angle = -core::f64::consts::TAU * k as f64 / len as f64;
                tw_re.push(angle.cos());
                tw_im.push(angle.sin());
            }
            len <<= 1;
        }
        Self { n, bitrev, tw_re, tw_im }
    }

    pub fn fft(&self, re: &mut [f64], im: &mut [f64], invert: bool) {
        let n = self.n;
        debug_assert!(re.len() == n && im.len() == n);

        for i in 1..n {
            let j = self.bitrev[i] as usize;
            if i < j {
                re.swap(i, j);
                im.swap(i, j);
            }
        }

        let mut tw_base = 0usize;
        let mut len = 2usize;
        while len <= n {
            let half = len / 2;
            let mut block = 0usize;
            while block < n {
                for k in 0..half {
                    let w_re = self.tw_re[tw_base + k];
                    let w_im = if invert { -self.tw_im[tw_base + k] } else { self.tw_im[tw_base + k] };
                    let a = block + k;
                    let b = a + half;
                    let u_re = re[a];
                    let u_im = im[a];
                    let v_re = re[b] * w_re - im[b] * w_im;
                    let v_im = re[b] * w_im + im[b] * w_re;
                    re[a] = u_re + v_re;
                    im[a] = u_im + v_im;
                    re[b] = u_re - v_re;
                    im[b] = u_im - v_im;
                }
                block += len;
            }
            tw_base += half;
            len <<= 1;
        }

        if invert {
            let inv = 1.0 / n as f64;
            for v in re.iter_mut() {
                *v *= inv;
            }
            for v in im.iter_mut() {
                *v *= inv;
            }
        }
    }
}

/// Row padding for the FFT working planes, in f64 slots (one cache line).
/// Power-of-two row strides are a cache-set aliasing disaster during the
/// transpose (at 512² every column access lands in the same set group and
/// the pass runs at memory-conflict speed, ~3× the whole convolution cost
/// under wasm). Padding each row by a line breaks the aliasing; rows stay
/// contiguous for the FFT, only the row-to-row stride changes.
const ROW_PAD: usize = 8;

/// In-place square transpose of the n×n prefix of a stride-row matrix,
/// cache-blocked (each plane is multiple MB at 512²).
fn transpose_square(m: &mut [f64], n: usize, stride: usize) {
    const B: usize = 32;
    let mut ii = 0usize;
    while ii < n {
        let i_end = (ii + B).min(n);
        // Diagonal block: swap only the upper triangle within the block.
        for i in ii..i_end {
            for j in (i + 1)..i_end {
                m.swap(i * stride + j, j * stride + i);
            }
        }
        // Off-diagonal block pairs.
        let mut jj = i_end;
        while jj < n {
            let j_end = (jj + B).min(n);
            for i in ii..i_end {
                for j in jj..j_end {
                    m.swap(i * stride + j, j * stride + i);
                }
            }
            jj += B;
        }
        ii += B;
    }
}

/// Half of a 2D transform: FFT every row, transpose, FFT every row again —
/// leaving the TRANSPOSED 2D spectrum. Composing half(fwd) → pointwise
/// multiply → half(inv) computes an exact circular convolution while
/// skipping two full transposes per convolution, PROVIDED the kernel
/// spectrum is transpose-symmetric. The ring kernel depends only on
/// distance, so k(x, y) = k(y, x) and its spectrum inherits the symmetry.
fn fft2d_half(re: &mut [f64], im: &mut [f64], n: usize, stride: usize, invert: bool, plan: &FftPlan) {
    for row in 0..n {
        let start = row * stride;
        plan.fft(&mut re[start..start + n], &mut im[start..start + n], invert);
    }
    transpose_square(re, n, stride);
    transpose_square(im, n, stride);
    for row in 0..n {
        let start = row * stride;
        plan.fft(&mut re[start..start + n], &mut im[start..start + n], invert);
    }
}

/// Full 2D FFT (spatial orientation preserved). Used for the one-time kernel
/// spectrum build; the per-step path uses `fft2d_half`.
fn fft2d(re: &mut [f64], im: &mut [f64], n: usize, stride: usize, invert: bool, plan: &FftPlan) {
    fft2d_half(re, im, n, stride, invert, plan);
    transpose_square(re, n, stride);
    transpose_square(im, n, stride);
}

// ---------------------------------------------------------------------------
// One toroidal Lenia field with an FFT-diagonalized ring kernel
// ---------------------------------------------------------------------------

pub struct LeniaCore {
    pub n: usize,
    /// Row-to-row stride of the FFT working planes (n + ROW_PAD). The field
    /// itself stays dense; only the spectral scratch is padded.
    pub stride: usize,
    /// Ring kernel radius in cells (fractional radii are honored exactly).
    pub radius: f64,
    pub field: Vec<f64>,
    // FFT scratch (re/im of the working spectrum), stride-padded rows.
    re: Vec<f64>,
    im: Vec<f64>,
    plan: FftPlan,
    // Forward FFT of the normalized kernel (fixed for the core's lifetime),
    // same padded layout as the scratch planes.
    // Transpose-symmetric because the kernel is a pure function of distance,
    // which is what lets convolve() use the transpose-sparing half passes.
    kr: Vec<f64>,
    ki: Vec<f64>,
    // Growth lookup table; rebuilt only when (mu, sigma) actually change
    // (they are live sliders, constant across almost every step).
    lut: Vec<f64>,
    lut_key: Option<(u64, u64)>,
}

impl LeniaCore {
    pub fn new(n: usize, radius: f64, ring_r0: f64, ring_sigma: f64) -> Self {
        let stride = n + ROW_PAD;
        let padded = n * stride;
        let plan = FftPlan::new(n);
        let mut kr = vec![0.0f64; padded];
        let mut ki = vec![0.0f64; padded];

        // Build the ring kernel centered at the origin with toroidal wrap,
        // matching the site's TS fallback: w = exp(-((r/R - r0)/sigma_k)^2/2)
        // for r <= R, then normalize to unit sum so K*a stays in [0, 1].
        let mut total = 0.0f64;
        for y in 0..n {
            let dy = if y > n / 2 { (n - y) as f64 } else { y as f64 };
            for x in 0..n {
                let dx = if x > n / 2 { (n - x) as f64 } else { x as f64 };
                let dist = (dx * dx + dy * dy).sqrt();
                if dist <= radius {
                    let d = (dist / radius - ring_r0) / ring_sigma;
                    let w = (-0.5 * d * d).exp();
                    kr[y * stride + x] = w;
                    total += w;
                }
            }
        }
        let inv_total = 1.0 / total.max(f64::MIN_POSITIVE);
        for w in kr.iter_mut() {
            *w *= inv_total;
        }
        fft2d(&mut kr, &mut ki, n, stride, false, &plan);

        Self {
            n,
            stride,
            radius,
            field: vec![0.0; n * n],
            re: vec![0.0; padded],
            im: vec![0.0; padded],
            plan,
            kr,
            ki,
            lut: vec![0.0; LUT_SIZE + 1],
            lut_key: None,
        }
    }

    pub fn clear(&mut self) {
        self.field.fill(0.0);
    }

    /// Additively seed a hollow gaussian ring (the site's soliton brush):
    /// ring = exp(-((dist - radius*ring_frac)/width)^2), clamped to 1.
    pub fn seed_ring(
        &mut self,
        cx: f64,
        cy: f64,
        radius: f64,
        ring_frac: f64,
        width: f64,
        intensity: f64,
    ) {
        let n = self.n as i64;
        let r = radius.ceil() as i64;
        let cxi = cx.round() as i64;
        let cyi = cy.round() as i64;
        for dy in -r..=r {
            for dx in -r..=r {
                let dist = ((dx * dx + dy * dy) as f64).sqrt();
                if dist <= radius {
                    let x = ((cxi + dx) % n + n) % n;
                    let y = ((cyi + dy) % n + n) % n;
                    let arg = (dist - radius * ring_frac) / width;
                    let ring = (-arg * arg).exp();
                    let idx = (y * n + x) as usize;
                    self.field[idx] = (self.field[idx] + ring * intensity).min(1.0);
                }
            }
        }
    }

    fn rebuild_lut(&mut self, mu: f64, sigma: f64) {
        let key = (mu.to_bits(), sigma.to_bits());
        if self.lut_key == Some(key) {
            return;
        }
        let inv_sigma = 1.0 / sigma.max(1e-4);
        for (i, slot) in self.lut.iter_mut().enumerate() {
            let u = i as f64 / LUT_SIZE as f64;
            let d = (u - mu) * inv_sigma;
            *slot = 2.0 * (-0.5 * d * d).exp() - 1.0;
        }
        self.lut_key = Some(key);
    }

    /// Convolve the current field with the ring kernel into `self.re`
    /// (the neighborhood potential u = K * a, exact toroidal convolution).
    /// The pointwise multiply runs in the TRANSPOSED spectral layout the
    /// half passes leave behind, which is valid because the kernel spectrum
    /// equals its own transpose (distance-only kernel).
    pub fn convolve(&mut self) {
        let n = self.n;
        let stride = self.stride;
        for row in 0..n {
            self.re[row * stride..row * stride + n].copy_from_slice(&self.field[row * n..row * n + n]);
        }
        self.im.fill(0.0);
        fft2d_half(&mut self.re, &mut self.im, n, stride, false, &self.plan);
        for row in 0..n {
            let base = row * stride;
            for i in base..base + n {
                let ar = self.re[i];
                let ai = self.im[i];
                let br = self.kr[i];
                let bi = self.ki[i];
                self.re[i] = ar * br - ai * bi;
                self.im[i] = ar * bi + ai * br;
            }
        }
        fft2d_half(&mut self.re, &mut self.im, n, stride, true, &self.plan);
    }

    /// One Lenia step. Returns (interface fraction, mean mass) of the
    /// updated field — the same two metrics the site's TS fallback reports.
    pub fn step(&mut self, mu: f64, sigma: f64, dt: f64) -> (f64, f64) {
        self.rebuild_lut(mu, sigma);
        self.convolve();

        let n = self.n;
        let stride = self.stride;
        let cells = n * n;
        let scale = LUT_SIZE as f64;
        let mut active = 0usize;
        let mut mass = 0.0f64;
        for y in 0..n {
            let pot_base = y * stride;
            let row_base = y * n;
            for x in 0..n {
                // Kernel is unit-normalized and the field lives in [0, 1], so
                // u only leaves [0, 1] by FFT roundoff; clamp before the LUT.
                let u = self.re[pot_base + x].clamp(0.0, 1.0);
                let t = u * scale;
                let idx = t as usize; // 0..=LUT_SIZE, last slot duplicated below
                let frac = t - idx as f64;
                let g = self.lut[idx] + (self.lut[(idx + 1).min(LUT_SIZE)] - self.lut[idx]) * frac;
                let updated = (self.field[row_base + x] + dt * g).clamp(0.0, 1.0);
                self.field[row_base + x] = updated;
                mass += updated;
                if updated > 0.08 && updated < 0.92 {
                    active += 1;
                }
            }
        }
        (active as f64 / cells as f64, mass / cells as f64)
    }
}

// ---------------------------------------------------------------------------
// The site-facing simulation: a display field, a half-resolution evaluation
// field for CMA-ES rollouts, and an RGBA colormap buffer.
// ---------------------------------------------------------------------------

pub struct LeniaSim {
    pub display: LeniaCore,
    pub eval: LeniaCore,
    eval_seed: Vec<f64>,
    rgba: Vec<u8>,
}

impl LeniaSim {
    pub fn new(size: usize, eval_size: usize, rel_radius: f64) -> Self {
        let display = LeniaCore::new(size, size as f64 * rel_radius, 0.5, 0.18);
        let eval = LeniaCore::new(eval_size, eval_size as f64 * rel_radius, 0.5, 0.18);
        let eval_cells = eval_size * eval_size;
        Self {
            display,
            eval,
            eval_seed: vec![0.0; eval_cells],
            rgba: vec![0u8; size * size * 4],
        }
    }

    /// Box-average the display field down into the evaluation seed
    /// (display size is an integer multiple of the eval size).
    pub fn snapshot_eval_seed(&mut self) {
        let n = self.display.n;
        let m = self.eval.n;
        let factor = n / m;
        let inv = 1.0 / (factor * factor) as f64;
        for ey in 0..m {
            for ex in 0..m {
                let mut sum = 0.0;
                for sy in 0..factor {
                    for sx in 0..factor {
                        sum += self.display.field[(ey * factor + sy) * n + ex * factor + sx];
                    }
                }
                self.eval_seed[ey * m + ex] = sum * inv;
            }
        }
    }

    /// Score a parameter triple from the frozen snapshot: the site's
    /// objective, mean over steps of (interface - 2·|mass - 0.25|).
    pub fn eval_score(&mut self, mu: f64, sigma: f64, dt: f64, steps: usize) -> f64 {
        self.eval.field.copy_from_slice(&self.eval_seed);
        let mut sum = 0.0;
        for _ in 0..steps {
            let (interface, mass) = self.eval.step(mu, sigma, dt);
            sum += interface - (mass - 0.25).abs() * 2.0;
        }
        sum / steps.max(1) as f64
    }

    /// Write the bioluminescent colormap of the display field into the RGBA
    /// buffer (byte-identical ramps to the site's TS fallback).
    pub fn render(&mut self) {
        for (i, &v) in self.display.field.iter().enumerate() {
            let o = i * 4;
            let (r, g, b) = if v < 0.015 {
                (3u8, 7u8, 18u8)
            } else if v < 0.35 {
                let t = v / 0.35;
                ((3.0 + t * 45.0) as u8, (7.0 + t * 180.0) as u8, (18.0 + t * 210.0) as u8)
            } else if v < 0.7 {
                let t = (v - 0.35) / 0.35;
                ((48.0 + t * 155.0) as u8, (187.0 - t * 85.0) as u8, (228.0 + t * 15.0) as u8)
            } else {
                let t = (v - 0.7) / 0.3;
                ((203.0 + t * 52.0) as u8, (102.0 + t * 150.0) as u8, (243.0 + t * 12.0) as u8)
            };
            self.rgba[o] = r;
            self.rgba[o + 1] = g;
            self.rgba[o + 2] = b;
            self.rgba[o + 3] = 255;
        }
    }
}

thread_local! {
    static SIM: RefCell<Option<LeniaSim>> = const { RefCell::new(None) };
}

// ---------------------------------------------------------------------------
// Envelope-returning core entries (shared by native tests and the JS boundary)
// ---------------------------------------------------------------------------

pub fn lenia_init_core(size: u32, eval_size: u32, rel_radius: f64) -> Result<(usize, f64), Refusal> {
    require_finite("rel_radius", rel_radius)?;
    let size = size as usize;
    let eval_size = eval_size as usize;
    if !size.is_power_of_two() || !(64..=512).contains(&size) {
        return Err(Refusal {
            code: "size-out-of-range",
            message: format!("size must be a power of two in 64..=512, got {size}"),
            ranked_repairs: vec!["use 128, 256, or 512"],
        });
    }
    if !eval_size.is_power_of_two() || eval_size < 32 || size % eval_size != 0 {
        return Err(Refusal {
            code: "eval-size-invalid",
            message: format!("eval_size must be a power of two >= 32 dividing size, got {eval_size}"),
            ranked_repairs: vec!["use size / 2"],
        });
    }
    if !(0.01..=0.25).contains(&rel_radius) {
        return Err(Refusal {
            code: "rel-radius-out-of-range",
            message: format!("rel_radius must lie in 0.01..=0.25, got {rel_radius}"),
            ranked_repairs: vec!["use ~0.052 (the site's 5/96 ratio)"],
        });
    }
    let sim = LeniaSim::new(size, eval_size, rel_radius);
    let radius = sim.display.radius;
    SIM.with(|slot| *slot.borrow_mut() = Some(sim));
    Ok((size, radius))
}

fn with_sim<T>(f: impl FnOnce(&mut LeniaSim) -> T) -> Result<T, Refusal> {
    SIM.with(|slot| {
        let mut guard = slot.borrow_mut();
        guard.as_mut().map(f).ok_or(Refusal {
            code: "sim-not-initialized",
            message: "call lenia_init before stepping".to_string(),
            ranked_repairs: vec!["call lenia_init(size, eval_size, rel_radius)"],
        })
    })
}

pub fn lenia_step_core(mu: f64, sigma: f64, dt: f64, steps: u32) -> Result<(f64, f64), Refusal> {
    require_finite("mu", mu)?;
    require_finite("sigma", sigma)?;
    require_finite("dt", dt)?;
    with_sim(|sim| {
        let mut metrics = (0.0, 0.0);
        for _ in 0..steps.max(1) {
            metrics = sim.display.step(mu, sigma, dt);
        }
        metrics
    })
}

pub fn lenia_eval_core(mu: f64, sigma: f64, dt: f64, steps: u32) -> Result<f64, Refusal> {
    require_finite("mu", mu)?;
    require_finite("sigma", sigma)?;
    require_finite("dt", dt)?;
    with_sim(|sim| sim.eval_score(mu, sigma, dt, steps.max(1) as usize))
}

// ---------------------------------------------------------------------------
// wasm32-only JS boundary.
// ---------------------------------------------------------------------------

#[cfg(target_arch = "wasm32")]
mod js {
    use wasm_bindgen::prelude::wasm_bindgen;

    /// Allocate the simulation. size: power of two in 64..=512; eval_size
    /// must divide it (the CMA-ES fitness rollouts run at this resolution);
    /// rel_radius: ring kernel radius as a fraction of size (~0.052 mirrors
    /// the site's 96-cell fallback). Returns `{"ok":{"size","kernelRadius"}}`.
    #[wasm_bindgen]
    #[must_use]
    pub fn lenia_init(size: u32, eval_size: u32, rel_radius: f64) -> String {
        match super::lenia_init_core(size, eval_size, rel_radius) {
            Ok((n, radius)) => format!(
                "{{\"ok\":{{\"kernel\":\"{}\",\"size\":{n},\"kernelRadius\":{radius}}}}}",
                super::KERNEL_VERSION
            ),
            Err(r) => r.json(),
        }
    }

    /// Zero the display field.
    #[wasm_bindgen]
    pub fn lenia_clear() {
        let _ = super::with_sim(|sim| sim.display.clear());
    }

    /// Additively seed a hollow gaussian ring at (cx, cy) in grid cells.
    #[wasm_bindgen]
    pub fn lenia_seed_ring(cx: f64, cy: f64, radius: f64, ring_frac: f64, width: f64, intensity: f64) {
        if ![cx, cy, radius, ring_frac, width, intensity].iter().all(|v| v.is_finite()) {
            return;
        }
        let _ = super::with_sim(|sim| sim.display.seed_ring(cx, cy, radius, ring_frac, width, intensity));
    }

    /// Advance the display field `steps` times. Returns the LAST step's
    /// metrics: `{"ok":{"interface":..,"mass":..}}`.
    #[wasm_bindgen]
    #[must_use]
    pub fn lenia_step(mu: f64, sigma: f64, dt: f64, steps: u32) -> String {
        match super::lenia_step_core(mu, sigma, dt, steps) {
            Ok((interface, mass)) => {
                format!("{{\"ok\":{{\"interface\":{interface},\"mass\":{mass}}}}}")
            }
            Err(r) => r.json(),
        }
    }

    /// Colormap the current display field into the internal RGBA buffer.
    #[wasm_bindgen]
    pub fn lenia_render() {
        let _ = super::with_sim(|sim| sim.render());
    }

    /// Pointer/length of the RGBA buffer inside wasm memory. The buffer is
    /// allocated once at init and never reallocated, so the pointer stays
    /// valid until the next lenia_init; the page rewraps it per frame in
    /// case wasm memory grows.
    #[wasm_bindgen]
    #[must_use]
    pub fn lenia_rgba_ptr() -> u32 {
        super::SIM.with(|slot| {
            slot.borrow()
                .as_ref()
                .map(|sim| sim.rgba.as_ptr() as u32)
                .unwrap_or(0)
        })
    }

    #[wasm_bindgen]
    #[must_use]
    pub fn lenia_rgba_len() -> u32 {
        super::SIM.with(|slot| {
            slot.borrow()
                .as_ref()
                .map(|sim| sim.rgba.len() as u32)
                .unwrap_or(0)
        })
    }

    /// Freeze the current display field (box-averaged to eval resolution) as
    /// the seed every subsequent lenia_eval rollout starts from.
    #[wasm_bindgen]
    pub fn lenia_snapshot_eval() {
        let _ = super::with_sim(|sim| sim.snapshot_eval_seed());
    }

    /// Fitness rollout from the frozen snapshot: `steps` Lenia steps at eval
    /// resolution, scored with the site's objective
    /// mean(interface - 2·|mass - 0.25|). Returns `{"ok":{"score":..}}`.
    #[wasm_bindgen]
    #[must_use]
    pub fn lenia_eval(mu: f64, sigma: f64, dt: f64, steps: u32) -> String {
        match super::lenia_eval_core(mu, sigma, dt, steps) {
            Ok(score) => format!("{{\"ok\":{{\"score\":{score}}}}}"),
            Err(r) => r.json(),
        }
    }

    /// Kernel identity probe (capability check after instantiation).
    #[wasm_bindgen]
    #[must_use]
    pub fn lenia_version() -> String {
        super::KERNEL_VERSION.to_string()
    }
}

// ---------------------------------------------------------------------------
// Tests — correctness of the FFT path against direct convolution, plus
// behavior and determinism checks. Native.
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    /// Deterministic pseudo-random field (LCG; no entropy in tests).
    fn test_field(n: usize, seed: u64) -> Vec<f64> {
        let mut state = seed;
        (0..n * n)
            .map(|_| {
                state = state.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
                ((state >> 33) as f64 / (1u64 << 31) as f64).fract().abs()
            })
            .collect()
    }

    #[test]
    fn fft_round_trips() {
        let n = 64;
        let plan = FftPlan::new(n);
        let mut re = test_field(8, 7)[..n].to_vec();
        let original = re.clone();
        let mut im = vec![0.0; n];
        plan.fft(&mut re, &mut im, false);
        plan.fft(&mut re, &mut im, true);
        for i in 0..n {
            assert!((re[i] - original[i]).abs() < 1e-12, "re[{i}] drifted");
            assert!(im[i].abs() < 1e-12, "im[{i}] nonzero");
        }
    }

    #[test]
    fn planned_fft_matches_direct_dft() {
        // Verify the table-driven FFT against the O(n²) textbook DFT, so the
        // fast path is checked against the definition rather than itself.
        let n = 32usize;
        let plan = FftPlan::new(n);
        let src_re = test_field(8, 11)[..n].to_vec();
        let src_im = test_field(8, 13)[..n].to_vec();

        let mut dft_re = vec![0.0f64; n];
        let mut dft_im = vec![0.0f64; n];
        for (k, (dr, di)) in dft_re.iter_mut().zip(dft_im.iter_mut()).enumerate() {
            for j in 0..n {
                let angle = -core::f64::consts::TAU * (k * j) as f64 / n as f64;
                let (c, s) = (angle.cos(), angle.sin());
                *dr += src_re[j] * c - src_im[j] * s;
                *di += src_re[j] * s + src_im[j] * c;
            }
        }

        let mut re = src_re.clone();
        let mut im = src_im.clone();
        plan.fft(&mut re, &mut im, false);
        for i in 0..n {
            assert!((re[i] - dft_re[i]).abs() < 1e-10, "re[{i}] vs DFT");
            assert!((im[i] - dft_im[i]).abs() < 1e-10, "im[{i}] vs DFT");
        }
    }

    #[test]
    fn blocked_transpose_is_a_transpose() {
        // Non-multiple-of-block size exercises the ragged tail blocks; the
        // strided case exercises the padded-row layout the cores use.
        for pad in [0usize, ROW_PAD] {
            let n = 48usize;
            let stride = n + pad;
            let mut m: Vec<f64> = (0..n * stride).map(|i| i as f64).collect();
            let original = m.clone();
            transpose_square(&mut m, n, stride);
            for i in 0..n {
                for j in 0..n {
                    assert_eq!(m[i * stride + j], original[j * stride + i], "({i},{j}) pad {pad}");
                }
            }
        }
    }

    #[test]
    fn fft_convolution_matches_direct_convolution() {
        // The load-bearing claim of this crate: the FFT path computes the
        // exact toroidal ring convolution the site's TS fallback computes
        // directly. Verify on a 32² field with a fractional radius.
        let n = 32usize;
        let radius = 3.3f64;
        let mut core = LeniaCore::new(n, radius, 0.5, 0.18);
        core.field = test_field(n, 42);

        // Direct convolution with the identical kernel definition.
        let r = radius.ceil() as i64;
        let mut weights = Vec::new();
        let mut total = 0.0;
        for dy in -r..=r {
            for dx in -r..=r {
                let dist = ((dx * dx + dy * dy) as f64).sqrt();
                if dist <= radius {
                    let d = (dist / radius - 0.5) / 0.18;
                    let w = (-0.5 * d * d).exp();
                    weights.push((dx, dy, w));
                    total += w;
                }
            }
        }
        let ni = n as i64;
        let mut direct = vec![0.0f64; n * n];
        for y in 0..ni {
            for x in 0..ni {
                let mut acc = 0.0;
                for &(dx, dy, w) in &weights {
                    let sx = ((x + dx) % ni + ni) % ni;
                    let sy = ((y + dy) % ni + ni) % ni;
                    acc += core.field[(sy * ni + sx) as usize] * w;
                }
                direct[(y * ni + x) as usize] = acc / total;
            }
        }

        core.convolve();
        let mut max_diff = 0.0f64;
        for y in 0..n {
            for x in 0..n {
                max_diff = max_diff.max((core.re[y * core.stride + x] - direct[y * n + x]).abs());
            }
        }
        assert!(max_diff < 1e-10, "FFT vs direct convolution diff {max_diff}");
    }

    #[test]
    fn step_keeps_field_in_unit_interval_and_metrics_sane() {
        let mut core = LeniaCore::new(64, 64.0 * 0.052, 0.5, 0.18);
        core.seed_ring(20.0, 24.0, 7.0, 0.55, 1.5, 0.95);
        core.seed_ring(44.0, 40.0, 7.0, 0.55, 1.5, 0.95);
        for _ in 0..60 {
            let (interface, mass) = core.step(0.152, 0.038, 0.22);
            assert!((0.0..=1.0).contains(&interface));
            assert!((0.0..=1.0).contains(&mass));
        }
        assert!(core.field.iter().all(|v| (0.0..=1.0).contains(v)));
    }

    #[test]
    fn deterministic_replay() {
        let run = || {
            let mut core = LeniaCore::new(64, 64.0 * 0.052, 0.5, 0.18);
            core.seed_ring(30.0, 30.0, 8.0, 0.55, 1.8, 0.9);
            let mut digest = 0.0f64;
            for _ in 0..40 {
                let (interface, mass) = core.step(0.2, 0.04, 0.2);
                digest += interface * 3.0 + mass;
            }
            digest
        };
        assert_eq!(run().to_bits(), run().to_bits(), "same inputs must replay bitwise");
    }

    #[test]
    fn snapshot_downsamples_by_box_average() {
        let mut sim = LeniaSim::new(64, 32, 0.052);
        // Constant field: any box average must reproduce the constant.
        sim.display.field.fill(0.37);
        sim.snapshot_eval_seed();
        assert!(sim.eval_seed.iter().all(|v| (v - 0.37).abs() < 1e-12));
        // A single hot display cell spreads 1/4 of its mass into its block.
        sim.display.field.fill(0.0);
        sim.display.field[0] = 1.0;
        sim.snapshot_eval_seed();
        assert!((sim.eval_seed[0] - 0.25).abs() < 1e-12);
        assert!(sim.eval_seed[1].abs() < 1e-12);
    }

    #[test]
    fn eval_scores_are_finite_and_replayable() {
        let mut sim = LeniaSim::new(64, 32, 0.052);
        sim.display.seed_ring(30.0, 30.0, 8.0, 0.55, 1.8, 0.9);
        sim.snapshot_eval_seed();
        let a = sim.eval_score(0.152, 0.038, 0.22, 18);
        let b = sim.eval_score(0.152, 0.038, 0.22, 18);
        assert!(a.is_finite());
        assert_eq!(a.to_bits(), b.to_bits(), "snapshot evals must be independent and identical");
        assert!((-1.0..=1.0).contains(&a), "score {a} outside plausible range");
    }

    #[test]
    fn init_core_refuses_bad_shapes() {
        assert_eq!(lenia_init_core(100, 50, 0.052).unwrap_err().code, "size-out-of-range");
        assert_eq!(lenia_init_core(256, 96, 0.052).unwrap_err().code, "eval-size-invalid");
        assert_eq!(lenia_init_core(256, 128, 0.5).unwrap_err().code, "rel-radius-out-of-range");
        assert!(lenia_init_core(256, 128, 0.052).is_ok());
        assert!(lenia_step_core(0.15, 0.04, 0.2, 2).is_ok());
    }

    #[test]
    fn render_paints_background_and_core_colors() {
        let mut sim = LeniaSim::new(64, 32, 0.052);
        sim.display.field[5] = 1.0;
        sim.render();
        // Background pixel: deep void #030712.
        assert_eq!(&sim.rgba[0..4], &[3, 7, 18, 255]);
        // Saturated pixel: top of the ramp (255, 252, 255).
        let o = 5 * 4;
        assert_eq!(sim.rgba[o + 3], 255);
        assert!(sim.rgba[o] >= 203);
    }
}

#[cfg(test)]
mod perf_probe {
    use super::*;

    #[test]
    #[ignore = "perf probe, run explicitly"]
    fn probe_step_cost() {
        for size in [256usize, 512] {
            let s = size as f64 / 96.0;
            let mut sim = LeniaSim::new(size, 128, 0.052);
            sim.display.seed_ring(33.6 * s, 43.2 * s, 10.0 * s, 0.55, 2.2 * s, 0.9);
            sim.display.seed_ring(62.4 * s, 52.8 * s, 10.0 * s, 0.55, 2.2 * s, 0.9);
            let t0 = std::time::Instant::now();
            for _ in 0..120 {
                sim.display.step(0.152, 0.038, 0.22);
            }
            let per_display = t0.elapsed().as_secs_f64() / 120.0;
            sim.snapshot_eval_seed();
            let t1 = std::time::Instant::now();
            for _ in 0..12 {
                sim.eval_score(0.152, 0.038, 0.22, 18);
            }
            let per_generation = t1.elapsed().as_secs_f64();
            let t2 = std::time::Instant::now();
            sim.render();
            let render = t2.elapsed().as_secs_f64();
            println!("display step {size}²: {:.3} ms", per_display * 1e3);
            println!("CMA-ES generation (12 cands × 18 steps @128²): {:.1} ms", per_generation * 1e3);
            println!("render {size}²: {:.3} ms", render * 1e3);
        }
    }
}
