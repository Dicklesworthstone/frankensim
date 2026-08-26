//! fs-heatmap-wasm — browser (WASM) surface for the CMA-ES explainer site's
//! objective-landscape heatmaps. Layer: L6.
//!
//! Seven site components rasterize a 2D objective into a background canvas
//! with the same template: pixel → (x, y) → v = field(x, y) → normalization
//! → linear color ramp → RGBA. In JS those loops burn 500k–1M transcendental
//! evaluations on the main thread at mount and again on every landscape /
//! strategy switch — a visible stall on phones. This kernel hosts the field
//! registry and rasterizes the identical template into an in-memory RGBA
//! buffer the page blits with zero per-pixel JS work.
//!
//! Pixel mapping (matches every site builder):
//!   x = xmin + (px / W) · (xmax − xmin)
//!   y = ymax − (py / H) · (ymax − ymin)      (row 0 is the TOP of the domain)
//!
//! Normalization modes (`norm_mode`, with scale `norm_k`):
//!   "log10p1"  clamp01(log10(1 + v) / k)
//!   "tanh"     tanh(v / k)                    (no clamp; mirrors the JS)
//!   "linear"   clamp01(v / k)
//!   "sqrt"     clamp01(sqrt(v) / k)
//!   "log10eps" clamp01(log10(max(1e-4, v + 1e-4)) / k)
//!
//! Ramp (all six coefficients supplied by the caller; `|0`-style truncation):
//!   r = trunc(r0 + rk·n)   g = trunc(g0 + gk·(1−n))   b = trunc(b0 + bk·(1−n))
//!
//! Contracts (fs-flyer-wasm pattern): typed-refusal JSON envelopes, nothing
//! silently clamped, no traps across the boundary, fully deterministic.
//!
//! No-claims: teaching/viz surface. The site keeps its JS loops as the
//! fallback; formulas here must stay identical to them.

use core::cell::RefCell;

/// Kernel id baked into envelopes so the page can prove which build is live.
pub const KERNEL_VERSION: &str = "fs-heatmap-wasm 0.1.0";

const MAX_DIM: u32 = 4096;
const MAX_PIXELS: u64 = 4_194_304;

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

// ---------------------------------------------------------------------------
// Field registry — 1:1 ports of the site's landscape functions. Any change
// on either side must be mirrored (the site's JS loops are the fallback).
// ---------------------------------------------------------------------------

const TAU: f64 = core::f64::consts::TAU;

fn rastrigin(x: f64, y: f64) -> f64 {
    20.0 + x * x + y * y - 10.0 * ((TAU * x).cos() + (TAU * y).cos())
}

fn ackley(x: f64, y: f64) -> f64 {
    let radius = (0.5 * (x * x + y * y)).sqrt();
    let ripple = 0.5 * ((TAU * x).cos() + (TAU * y).cos());
    -20.0 * (-0.2 * radius).exp() - ripple.exp() + 20.0 + core::f64::consts::E
}

/// The constraint demo's box objective: unconstrained minimum at
/// (1.35, 1.35), outside the [0, 1]² feasible box.
fn box_quad(x: f64, y: f64) -> f64 {
    (x - 1.35) * (x - 1.35) + 2.5 * (y - 1.35) * (y - 1.35) + 0.1 * (6.0 * x).sin()
}

fn reflect01(v: f64) -> f64 {
    if (0.0..=1.0).contains(&v) {
        return v;
    }
    let m = v.abs() % 2.0;
    if m > 1.0 { 2.0 - m } else { m }
}

fn logit_sig(v: f64) -> f64 {
    1.0 / (1.0 + (-v * 3.0).exp())
}

/// Resolve a field id to its evaluator. Ids are shared vocabulary with the
/// site (app/lib/frankensimHeatmap.ts).
pub fn field_fn(id: &str) -> Option<fn(f64, f64) -> f64> {
    Some(match id {
        // WasmDemo benchmark suite
        "rosenbrock100" => |x, y| 100.0 * (y - x * x) * (y - x * x) + (1.0 - x) * (1.0 - x),
        "rastrigin" => rastrigin,
        "ackley" => ackley,
        "cigar-y1000" => |x, y| x * x + 1000.0 * y * y,
        "himmelblau" => |x, y| {
            let a = x * x + y - 11.0;
            let b = x + y * y - 7.0;
            a * a + b * b
        },
        "step-ridge" => |x, y| (x * x + y * y).floor() + 0.5 * (x - y).abs(),
        // CmaesIntro sandbox landscapes
        "rosenbrock10" => |x, y| 10.0 * (y - x * x) * (y - x * x) + (1.0 - x) * (1.0 - x),
        "rot-cigar80" => |x, y| {
            let rot = 0.61f64;
            let u = x * rot.cos() + y * rot.sin();
            let v = -x * rot.sin() + y * rot.cos();
            u * u + 80.0 * v * v
        },
        // CovarianceMinimap objectives
        "cigar-x100" => |x, y| 100.0 * x * x + y * y,
        "sphere" => |x, y| x * x + y * y,
        // ActiveCovarianceDemo: sharp banana canyon with deceptive ridges
        "banana-canyon" => |x, y| {
            let c = y - 0.4 * x * x;
            let s = (4.0 * x).sin();
            80.0 * c * c + (x - 1.2) * (x - 1.2) + 5.0 * s * s
        },
        // NoiseExplorer: rippled bowl
        "bowl-ripple" => |x, y| x * x + y * y + 0.3 * (4.0 * x).cos() + 0.3 * (4.0 * y).sin(),
        // ConstraintRepairDemo: the box objective seen THROUGH each repair
        "box-quad-clamp" => |x, y| box_quad(x.clamp(0.0, 1.0), y.clamp(0.0, 1.0)),
        "box-quad-reflect" => |x, y| box_quad(reflect01(x), reflect01(y)),
        "box-quad-logit" => |x, y| box_quad(logit_sig(x - 0.5), logit_sig(y - 0.5)),
        _ => return None,
    })
}

fn apply_norm(mode: &str, v: f64, k: f64) -> Option<f64> {
    Some(match mode {
        "log10p1" => ((1.0 + v).log10() / k).clamp(0.0, 1.0),
        "tanh" => (v / k).tanh(),
        "linear" => (v / k).clamp(0.0, 1.0),
        "sqrt" => (v.sqrt() / k).clamp(0.0, 1.0),
        "log10eps" => ((v + 1e-4).max(1e-4).log10() / k).clamp(0.0, 1.0),
        _ => return None,
    })
}

// ---------------------------------------------------------------------------
// Rasterization core (shared by native tests and the JS boundary)
// ---------------------------------------------------------------------------

pub struct RampSpec {
    pub r0: f64,
    pub rk: f64,
    pub g0: f64,
    pub gk: f64,
    pub b0: f64,
    pub bk: f64,
}

thread_local! {
    static RGBA: RefCell<Vec<u8>> = const { RefCell::new(Vec::new()) };
}

#[allow(clippy::too_many_arguments)]
pub fn heatmap_render_core(
    field: &str,
    width: u32,
    height: u32,
    xmin: f64,
    xmax: f64,
    ymin: f64,
    ymax: f64,
    norm_mode: &str,
    norm_k: f64,
    ramp: &RampSpec,
) -> Result<(), Refusal> {
    for (name, v) in [
        ("xmin", xmin),
        ("xmax", xmax),
        ("ymin", ymin),
        ("ymax", ymax),
        ("norm_k", norm_k),
        ("r0", ramp.r0),
        ("rk", ramp.rk),
        ("g0", ramp.g0),
        ("gk", ramp.gk),
        ("b0", ramp.b0),
        ("bk", ramp.bk),
    ] {
        if !v.is_finite() {
            return Err(Refusal {
                code: "input-non-finite",
                message: format!("{name} must be finite, got {v}"),
                ranked_repairs: vec!["pass a finite number"],
            });
        }
    }
    if width < 8 || height < 8 || width > MAX_DIM || height > MAX_DIM || (width as u64) * (height as u64) > MAX_PIXELS {
        return Err(Refusal {
            code: "size-out-of-range",
            message: format!("width/height must lie in 8..=4096 with at most {MAX_PIXELS} pixels, got {width}x{height}"),
            ranked_repairs: vec!["use the site's native canvas resolution"],
        });
    }
    if !(xmax > xmin) || !(ymax > ymin) {
        return Err(Refusal {
            code: "domain-degenerate",
            message: format!("domain must satisfy xmax > xmin and ymax > ymin, got [{xmin},{xmax}]x[{ymin},{ymax}]"),
            ranked_repairs: vec!["swap or widen the bounds"],
        });
    }
    let Some(f) = field_fn(field) else {
        return Err(Refusal {
            code: "field-unknown",
            message: format!("no field named {field:?} in the registry"),
            ranked_repairs: vec!["see field_fn in fs-heatmap-wasm for the id list"],
        });
    };
    if apply_norm(norm_mode, 0.0, 1.0).is_none() {
        return Err(Refusal {
            code: "norm-unknown",
            message: format!("no normalization mode named {norm_mode:?}"),
            ranked_repairs: vec!["use log10p1, tanh, linear, sqrt, or log10eps"],
        });
    }
    if norm_k == 0.0 {
        return Err(Refusal {
            code: "norm-scale-zero",
            message: "norm_k must be nonzero".to_string(),
            ranked_repairs: vec!["pass the divisor the JS fallback uses"],
        });
    }

    let w = width as usize;
    let h = height as usize;
    RGBA.with(|slot| {
        let mut buf = slot.borrow_mut();
        buf.clear();
        buf.resize(w * h * 4, 0);
        let x_scale = (xmax - xmin) / width as f64;
        let y_scale = (ymax - ymin) / height as f64;
        for py in 0..h {
            let y = ymax - py as f64 * y_scale;
            let row = py * w * 4;
            for px in 0..w {
                let x = xmin + px as f64 * x_scale;
                let v = f(x, y);
                // Registered fields are finite on any finite input; the
                // unwrap-free fallback keeps a hostile NaN from trapping.
                let n = apply_norm(norm_mode, v, norm_k).unwrap_or(0.0);
                let n = if n.is_finite() { n } else { 0.0 };
                // `as u8` truncates toward zero exactly like the JS `| 0`
                // (ramps stay inside 0..=255 by construction on both sides).
                let o = row + px * 4;
                buf[o] = (ramp.r0 + ramp.rk * n) as u8;
                buf[o + 1] = (ramp.g0 + ramp.gk * (1.0 - n)) as u8;
                buf[o + 2] = (ramp.b0 + ramp.bk * (1.0 - n)) as u8;
                buf[o + 3] = 255;
            }
        }
    });
    Ok(())
}

/// Copy of the current buffer for native tests (the wasm boundary hands out
/// ptr/len instead).
pub fn rgba_snapshot() -> Vec<u8> {
    RGBA.with(|slot| slot.borrow().clone())
}

// ---------------------------------------------------------------------------
// wasm32-only JS boundary.
// ---------------------------------------------------------------------------

#[cfg(target_arch = "wasm32")]
mod js {
    use wasm_bindgen::prelude::wasm_bindgen;

    /// Rasterize one landscape heatmap into the internal RGBA buffer.
    /// Returns `{"ok":{"kernel","width","height"}}` or a typed refusal.
    #[wasm_bindgen]
    #[must_use]
    #[allow(clippy::too_many_arguments)]
    pub fn heatmap_render(
        field: String,
        width: u32,
        height: u32,
        xmin: f64,
        xmax: f64,
        ymin: f64,
        ymax: f64,
        norm_mode: String,
        norm_k: f64,
        r0: f64,
        rk: f64,
        g0: f64,
        gk: f64,
        b0: f64,
        bk: f64,
    ) -> String {
        let ramp = super::RampSpec { r0, rk, g0, gk, b0, bk };
        match super::heatmap_render_core(&field, width, height, xmin, xmax, ymin, ymax, &norm_mode, norm_k, &ramp) {
            Ok(()) => format!(
                "{{\"ok\":{{\"kernel\":\"{}\",\"width\":{width},\"height\":{height}}}}}",
                super::KERNEL_VERSION
            ),
            Err(r) => r.json(),
        }
    }

    /// Pointer/length of the RGBA buffer inside wasm memory. Valid until the
    /// next heatmap_render call; the page copies it into an ImageData
    /// immediately (putImageData copies, so no lifetime hazard).
    #[wasm_bindgen]
    #[must_use]
    pub fn heatmap_rgba_ptr() -> u32 {
        super::RGBA.with(|slot| slot.borrow().as_ptr() as u32)
    }

    #[wasm_bindgen]
    #[must_use]
    pub fn heatmap_rgba_len() -> u32 {
        super::RGBA.with(|slot| slot.borrow().len() as u32)
    }

    /// Kernel identity probe (capability check after instantiation).
    #[wasm_bindgen]
    #[must_use]
    pub fn heatmap_version() -> String {
        super::KERNEL_VERSION.to_string()
    }
}

// ---------------------------------------------------------------------------
// Tests — template correctness, orientation, refusals, determinism. The
// cross-implementation parity proof against the site's JS loops lives in the
// site repo's bun test (test-heatmap-parity.mjs), which diffs this kernel's
// output against the exact JS template for all seven component configs.
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    const RAMP: RampSpec = RampSpec { r0: 10.0, rk: 20.0, g0: 20.0, gk: 75.0, b0: 40.0, bk: 115.0 };

    fn render_small(field: &str, norm: &str, k: f64) -> Vec<u8> {
        heatmap_render_core(field, 16, 16, -2.5, 2.5, -2.5, 2.5, norm, k, &RAMP).unwrap();
        rgba_snapshot()
    }

    #[test]
    fn registry_covers_all_site_fields() {
        for id in [
            "rosenbrock100",
            "rastrigin",
            "ackley",
            "cigar-y1000",
            "himmelblau",
            "step-ridge",
            "rosenbrock10",
            "rot-cigar80",
            "cigar-x100",
            "sphere",
            "banana-canyon",
            "bowl-ripple",
            "box-quad-clamp",
            "box-quad-reflect",
            "box-quad-logit",
        ] {
            assert!(field_fn(id).is_some(), "missing field {id}");
            let f = field_fn(id).unwrap();
            assert!(f(0.3, -0.7).is_finite(), "{id} non-finite at probe point");
        }
        assert!(field_fn("nope").is_none());
    }

    #[test]
    fn template_pixel_matches_hand_computation() {
        // sphere, linear norm k=8, 16x16 over [-2.5, 2.5]²: verify one pixel
        // end to end. Pixel (3, 5): x = -2.5 + 3·(5/16), y = 2.5 − 5·(5/16).
        let buf = render_small("sphere", "linear", 8.0);
        let x: f64 = -2.5 + 3.0 * (5.0 / 16.0);
        let y: f64 = 2.5 - 5.0 * (5.0 / 16.0);
        let n = ((x * x + y * y) / 8.0).clamp(0.0, 1.0);
        let o = (5 * 16 + 3) * 4;
        assert_eq!(buf[o], (10.0 + 20.0 * n) as u8);
        assert_eq!(buf[o + 1], (20.0 + 75.0 * (1.0 - n)) as u8);
        assert_eq!(buf[o + 2], (40.0 + 115.0 * (1.0 - n)) as u8);
        assert_eq!(buf[o + 3], 255);
    }

    #[test]
    fn orientation_top_row_uses_ymax() {
        let ramp = RampSpec { r0: 0.0, rk: 200.0, g0: 0.0, gk: 0.0, b0: 0.0, bk: 0.0 };
        // Field increasing in y ⇒ top row (py = 0, y = ymax) must be redder
        // than the bottom row under an increasing red ramp.
        heatmap_render_core("bowl-ripple", 16, 16, 0.0, 1.0, 0.0, 4.0, "linear", 20.0, &ramp).unwrap();
        let buf = rgba_snapshot();
        let top_r = buf[(8 * 4) as usize];
        let bottom_r = buf[((15 * 16 + 8) * 4) as usize];
        assert!(top_r > bottom_r, "top {top_r} vs bottom {bottom_r}");
    }

    #[test]
    fn refusals_are_typed() {
        let ramp = RAMP;
        let bad = |f: &str, w: u32, xmax: f64, norm: &str, k: f64| {
            heatmap_render_core(f, w, 16, -1.0, xmax, -1.0, 1.0, norm, k, &ramp).unwrap_err().code
        };
        assert_eq!(bad("nope", 16, 1.0, "linear", 8.0), "field-unknown");
        assert_eq!(bad("sphere", 4, 1.0, "linear", 8.0), "size-out-of-range");
        assert_eq!(bad("sphere", 16, -1.0, "linear", 8.0), "domain-degenerate");
        assert_eq!(bad("sphere", 16, 1.0, "nope", 8.0), "norm-unknown");
        assert_eq!(bad("sphere", 16, 1.0, "linear", 0.0), "norm-scale-zero");
        assert_eq!(
            heatmap_render_core("sphere", 16, 16, f64::NAN, 1.0, -1.0, 1.0, "linear", 8.0, &ramp)
                .unwrap_err()
                .code,
            "input-non-finite"
        );
    }

    #[test]
    fn deterministic_replay() {
        let a = render_small("rastrigin", "linear", 50.0);
        let b = render_small("rastrigin", "linear", 50.0);
        assert_eq!(a, b, "same inputs must replay byte-identical");
    }

    #[test]
    fn constraint_fields_are_flat_inside_reachable_image() {
        // clamp and reflect map the whole plane into [0, 1]², so the field
        // value at (1.7, 1.7) equals the value at the clamped corner (1, 1).
        let f_clamp = field_fn("box-quad-clamp").unwrap();
        assert_eq!(f_clamp(1.7, 1.7), f_clamp(1.0, 1.0));
        let f_reflect = field_fn("box-quad-reflect").unwrap();
        assert!((f_reflect(1.25, 0.5) - f_reflect(0.75, 0.5)).abs() < 1e-15);
        let f_logit = field_fn("box-quad-logit").unwrap();
        assert!((f_logit(0.5, 0.5) - box_quad(0.5, 0.5)).abs() < 1e-12);
    }
}
