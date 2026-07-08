//! certified.rs — the certified-numerics showcase:
//!
//! * [`mandelbrot_certified`] — a Mandelbrot render where each pixel's
//!   exterior classification is RIGOROUS: the whole pixel box is iterated with
//!   `fs-ivl` outward-rounded interval arithmetic, so a pixel is only marked
//!   "escaped" when EVERY parameter `c` in the pixel is provably outside the
//!   set (the guaranteed lower bound of `|z|²` exceeds 4). No floating-point
//!   guesswork at the boundary.

use fs_ivl::Interval;
use fs_math::det;

/// Render the Mandelbrot set with rigorous interval arithmetic. For each pixel
/// the parameter `c` ranges over the pixel's box `[cr]×[ci]`; the orbit
/// `z ← z² + c` is iterated in outward-rounded `fs-ivl` intervals. A pixel is
/// certified EXTERIOR the first iteration its `|z|²` interval has a lower bound
/// strictly greater than 4 (every `c` in the box has escaped); otherwise it is
/// reported as in-set / undetermined.
///
/// The view is centred at `(cx, cy)` with half-width `scale` in the real axis
/// (the imaginary extent is aspect-corrected).
///
/// Output: `w*h` values, row-major `index = py*w + px` (`py` increasing
/// downward). Value semantics:
/// - `0.0`  → the pixel did NOT certifiably escape within `maxiter`
///   (rigorously "not proven exterior": interior or boundary).
/// - `> 0.0` → CERTIFIED exterior; the value is the smooth escape iteration
///   `≈ it + 1 − log₂(log|z|)` (in `[0.5, maxiter]`) for continuous colouring.
///
/// `w,h` clamped to `[1,360]`, `maxiter` to `[1,500]`, `scale` to
/// `[1e-6,10]`, `cx,cy` to `[-3,3]`.
pub fn mandelbrot_certified(
    w_in: usize,
    h_in: usize,
    cx: f64,
    cy: f64,
    scale_in: f64,
    maxiter_in: usize,
) -> Vec<f64> {
    let w = w_in.clamp(1, 360);
    let h = h_in.clamp(1, 360);
    let scale = scale_in.clamp(1.0e-6, 10.0);
    let maxiter = maxiter_in.clamp(1, 500);
    let cxc = cx.clamp(-3.0, 3.0);
    let cyc = cy.clamp(-3.0, 3.0);

    let aspect = h as f64 / w as f64;
    let two = Interval::point(2.0);
    let ln2 = std::f64::consts::LN_2;
    let mut out = vec![0.0f64; w * h];

    for py in 0..h {
        let fy0 = cyc - scale * aspect + 2.0 * scale * aspect * (py as f64) / (h as f64);
        let fy1 = cyc - scale * aspect + 2.0 * scale * aspect * ((py as f64) + 1.0) / (h as f64);
        let ci = Interval::new(fy0.min(fy1), fy0.max(fy1));
        for px in 0..w {
            let fx0 = cxc - scale + 2.0 * scale * (px as f64) / (w as f64);
            let fx1 = cxc - scale + 2.0 * scale * ((px as f64) + 1.0) / (w as f64);
            let cr = Interval::new(fx0.min(fx1), fx0.max(fx1));

            let mut zr = Interval::point(0.0);
            let mut zi = Interval::point(0.0);
            let mut escaped = 0.0f64;
            for it in 0..maxiter {
                // z ← z² + c  (all outward-rounded interval arithmetic).
                let new_zr = (zr * zr - zi * zi) + cr;
                let new_zi = (two * zr * zi) + ci;
                zr = new_zr;
                zi = new_zi;
                let r2 = zr * zr + zi * zi;
                if r2.lo() > 4.0 {
                    // Certified escape: smooth iteration from the midpoint |z|.
                    let zn = r2.midpoint().max(4.000_1).sqrt();
                    let nu = (it as f64) + 1.0 - det::ln(det::ln(zn).max(1.0e-9)) / ln2;
                    escaped = nu.max(0.5);
                    break;
                }
            }
            out[py * w + px] = escaped;
        }
    }
    out
}
