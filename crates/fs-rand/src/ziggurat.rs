//! Ziggurat standard-normal sampler — the PERF path for `next_normal`
//! (bead frankensim-1za9). Box–Muller stays the strict default; this is
//! FAST-MODE-ONLY.
//!
//! # Determinism
//!
//! The 128-layer table is generated at first use from two constants (the tail
//! boundary `R` and the common layer area `V`) using ONLY `fs_math::det`
//! (correctly-rounded pure-IEEE `exp`/`ln`/`sqrt` — NOT platform libm), so the
//! table is bit-identical on every conforming ISA. Sampling consumes stream
//! draws under the same DETERMINISTIC REJECTION contract as `Stream::next_below`
//! (every rejected draw advances the index), so the consumed count is a pure
//! function of the stream content and the whole thing replays bitwise.
//!
//! # No-claim (strict mode)
//!
//! Admitting the ziggurat to STRICT mode needs a proven cross-ISA bitwise-equal
//! run on both reference machines (the `trj` pipeline). Until that lands it is
//! gated FAST-MODE-ONLY; the strict `Stream::next_normal` (Box–Muller) is
//! unchanged. Accuracy IS gated here (moments + a two-sample KS vs Box–Muller).

use crate::Stream;
use fs_math::det;
use std::sync::OnceLock;

/// Layer count (a power of two so the layer index is a bit mask).
const N: usize = 128;
/// `N − 1` as an `i32` bit mask; the const-assert pins the two in sync.
const LAYER_MASK: i32 = 0x7F;
const _: () = assert!(LAYER_MASK as usize + 1 == N);
/// The tail boundary `x[N-1]` and common layer area for `N = 128` (Marsaglia &
/// Tsang) — the "bit-pattern constants" the bead permits; every other table
/// entry is derived from these by `det` transcendentals.
const R: f64 = 3.442_619_855_899;
const V: f64 = 9.912_563_035_262_17e-3;
/// `2³¹` — the magnitude scale (a 32-bit signed draw carries sign+layer+value).
const M1: f64 = 2_147_483_648.0;

struct Tables {
    /// Integer acceptance thresholds `k[i]`.
    k: [u32; N],
    /// Scale `w[i]` mapping the integer magnitude to `x`.
    w: [f64; N],
    /// Layer heights `f[i] = exp(-x[i]²/2)`.
    f: [f64; N],
}

fn nexp_half_sq(x: f64) -> f64 {
    det::exp(-0.5 * x * x)
}

fn build() -> Tables {
    let mut k = [0u32; N];
    let mut w = [0.0f64; N];
    let mut f = [0.0f64; N];
    let mut dn = R;
    let q = V / nexp_half_sq(dn);
    k[0] = (dn / q * M1) as u32;
    k[1] = 0;
    w[0] = q / M1;
    w[N - 1] = dn / M1;
    f[0] = 1.0;
    f[N - 1] = nexp_half_sq(dn);
    let mut tn = dn;
    for i in (1..=(N - 2)).rev() {
        dn = det::sqrt(-2.0 * det::ln(V / dn + nexp_half_sq(dn)));
        k[i + 1] = (dn / tn * M1) as u32;
        tn = dn;
        f[i] = nexp_half_sq(dn);
        w[i] = dn / M1;
    }
    Tables { k, w, f }
}

fn tables() -> &'static Tables {
    static T: OnceLock<Tables> = OnceLock::new();
    T.get_or_init(build)
}

/// A uniform draw in `(0, 1]` (guards the tail's `ln`).
fn uni_pos(s: &mut Stream) -> f64 {
    1.0 - s.next_f64()
}

/// One standard-normal sample by the ziggurat method. Consumes one `u64` per
/// attempt plus, on a wedge/tail fallback, extra uniforms — all deterministic.
#[must_use]
pub fn normal(s: &mut Stream) -> f64 {
    let t = tables();
    loop {
        let hz = s.next_u64() as i32;
        let iz = (hz & LAYER_MASK) as usize;
        if hz.unsigned_abs() < t.k[iz] {
            return f64::from(hz) * t.w[iz];
        }
        // Fallback: exact tail (base layer) or the wedge acceptance test.
        if iz == 0 {
            // Marsaglia's exponential-tail sampler beyond R.
            loop {
                let x = -det::ln(uni_pos(s)) / R;
                let y = -det::ln(uni_pos(s));
                if y + y >= x * x {
                    return if hz > 0 { R + x } else { -(R + x) };
                }
            }
        } else {
            let x = f64::from(hz) * t.w[iz];
            if (t.f[iz] + s.next_f64() * (t.f[iz - 1] - t.f[iz])) < nexp_half_sq(x) {
                return x;
            }
        }
        // Rejected: the index already advanced; draw again.
    }
}
