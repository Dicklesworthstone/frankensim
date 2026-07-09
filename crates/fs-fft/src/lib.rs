//! fs-fft — Stockham autosort FFTs, real transforms, and DCT (plan §6.3).
//!
//! v1 is CORRECTNESS-FIRST: radix-2 Stockham (no bit-reversal pass), the
//! packed real transform (r2c/c2r at half the complex work), and DCT-II/III
//! via FFT folding — the fs-cheb dependency. The correctness oracle is the
//! naive O(n²) DFT, run exhaustively over sizes and random inputs: the
//! ORACLE decides, not derivation confidence.
//!
//! Determinism: twiddles come from fs-math's STRICT sin/cos and every
//! butterfly runs in a fixed order — the transform is cross-ISA
//! bit-deterministic by construction (golden-hash tested, verified on both
//! reference ISAs like the rest of the numerics spine).
//!
//! Bead fs-fft-perf-multidim extends the correctness core with the r2c INVERSE
//! (c2r, [`RealFft::inverse`]) and N-DIMENSIONAL (2D/3D) transforms via
//! separable pencil decomposition ([`FftNd`]) — both oracle-tested. Still
//! deferred to that bead's perf scope: radix-4/8 kernels, SIMD lanes,
//! cache-blocked transposes, executor-tiled pencils, and the
//! ≥40%-of-memory-bound roofline gate.
//!
//! Conventions: forward is unnormalized; `inverse` scales by 1/n so
//! inverse(forward(x)) = x. Sizes must be powers of two (structured
//! rejection otherwise; mixed-radix general n is out of v1 scope).

use fs_math::det;

/// Crate version, re-exported for provenance stamping.
pub const VERSION: &str = env!("CARGO_PKG_VERSION");

/// A complex number (f64 re/im). Local, minimal: fs-la's future complex
/// types can absorb this when a shared home exists.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct C64 {
    /// Real part.
    pub re: f64,
    /// Imaginary part.
    pub im: f64,
}

impl C64 {
    /// Construct.
    #[must_use]
    pub const fn new(re: f64, im: f64) -> C64 {
        C64 { re, im }
    }

    /// Complex conjugate.
    #[must_use]
    pub const fn conj(self) -> C64 {
        C64 {
            re: self.re,
            im: -self.im,
        }
    }

    /// Squared magnitude.
    #[must_use]
    pub fn norm_sq(self) -> f64 {
        self.re.mul_add(self.re, self.im * self.im)
    }

    fn add(self, o: C64) -> C64 {
        C64 {
            re: self.re + o.re,
            im: self.im + o.im,
        }
    }

    fn sub(self, o: C64) -> C64 {
        C64 {
            re: self.re - o.re,
            im: self.im - o.im,
        }
    }

    fn mul(self, o: C64) -> C64 {
        C64 {
            re: self.re.mul_add(o.re, -(self.im * o.im)),
            im: self.re.mul_add(o.im, self.im * o.re),
        }
    }

    fn scale(self, s: f64) -> C64 {
        C64 {
            re: self.re * s,
            im: self.im * s,
        }
    }
}

/// Planned FFT for one power-of-two size: precomputed strict-mode twiddle
/// table, reusable across transforms (twiddle generation is the expensive
/// deterministic part; plans are cheap to keep).
#[derive(Debug, Clone)]
pub struct Fft {
    n: usize,
    /// w[k] = exp(−2πik/n) for k in 0..n/2.
    table: Vec<C64>,
}

impl Fft {
    /// Plan a transform of size `n` (power of two, ≥ 1).
    ///
    /// # Panics
    /// With a structured message if `n` is not a power of two — general
    /// mixed-radix sizes are a recorded follow-up, and silently computing a
    /// wrong-size transform would be worse than refusing.
    #[must_use]
    pub fn new(n: usize) -> Fft {
        assert!(
            n >= 1 && n.is_power_of_two(),
            "FFT size {n} is not a power of two; pad or resample (mixed-radix general-n \
             support is recorded follow-up scope)"
        );
        let mut table = Vec::with_capacity(n / 2);
        for k in 0..n / 2 {
            // Strict-mode twiddles: deterministic cross-ISA. The angle stays
            // well inside fs-math's trig domain for any practical n.
            let theta = -2.0 * std::f64::consts::PI * (k as f64) / (n as f64);
            table.push(C64::new(det::cos(theta), det::sin(theta)));
        }
        Fft { n, table }
    }

    /// Transform size.
    #[must_use]
    pub fn len(&self) -> usize {
        self.n
    }

    /// True iff the planned size is 1 (identity transform).
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.n <= 1
    }

    /// Forward DFT, unnormalized: X[k] = Σ_j x[j]·exp(−2πijk/n).
    /// `scratch` must be the same length as `data`.
    pub fn forward(&self, data: &mut [C64], scratch: &mut [C64]) {
        self.transform(data, scratch, false);
    }

    /// Inverse DFT with 1/n normalization: inverse(forward(x)) = x.
    pub fn inverse(&self, data: &mut [C64], scratch: &mut [C64]) {
        self.transform(data, scratch, true);
        let inv_n = 1.0 / (self.n as f64);
        for v in data.iter_mut() {
            *v = v.scale(inv_n);
        }
    }

    /// Radix-2 Stockham autosort: ping-pongs between `data` and `scratch`
    /// with no bit-reversal pass; butterfly order is a pure function of the
    /// stage structure (deterministic by construction).
    fn transform(&self, data: &mut [C64], scratch: &mut [C64], inverse: bool) {
        let n = self.n;
        assert_eq!(data.len(), n, "data length must equal the planned size {n}");
        assert_eq!(
            scratch.len(),
            n,
            "scratch length must equal the planned size {n}"
        );
        if n == 1 {
            return;
        }
        // DIF Stockham (OTFFT formulation): at each stage the current
        // transform length n_cur halves while the stride s doubles; the
        // butterfly is (a+b, (a−b)·w_p) with w_p = exp(−2πi·p/n_cur)
        // = table[p·s] (since s = n/n_cur). Autosorting: no reversal pass.
        let mut n_cur = n;
        let mut s = 1usize;
        let mut src_is_data = true;
        while n_cur > 1 {
            let m = n_cur / 2;
            {
                let (src, dst): (&[C64], &mut [C64]) = if src_is_data {
                    (&*data, &mut *scratch)
                } else {
                    (&*scratch, &mut *data)
                };
                for p in 0..m {
                    let mut w = self.table[p * s];
                    if inverse {
                        w = w.conj();
                    }
                    for q in 0..s {
                        let a = src[q + s * p];
                        let b = src[q + s * (p + m)];
                        dst[q + s * 2 * p] = a.add(b);
                        dst[q + s * (2 * p + 1)] = a.sub(b).mul(w);
                    }
                }
            }
            n_cur = m;
            s *= 2;
            src_is_data = !src_is_data;
        }
        if !src_is_data {
            data.copy_from_slice(scratch);
        }
    }
}

// ---------------------------------------------------------------------------
// Real transforms: r2c via half-size complex packing (the bandwidth win the
// plan asks for — real fields are the common case).
// ---------------------------------------------------------------------------

/// Planned real-input FFT of even size `n` (power of two ≥ 2): produces the
/// n/2 + 1 non-redundant spectrum bins of a real signal.
#[derive(Debug, Clone)]
pub struct RealFft {
    half: Fft,
    /// Untangling twiddles: exp(−πik/(n/2))… indexed k in 0..n/4+1 range use.
    n: usize,
}

impl RealFft {
    /// Plan for real length `n` (power of two, ≥ 2).
    #[must_use]
    pub fn new(n: usize) -> RealFft {
        assert!(
            n >= 2 && n.is_power_of_two(),
            "real FFT size {n} must be a power of two >= 2"
        );
        RealFft {
            half: Fft::new(n / 2),
            n,
        }
    }

    /// Forward r2c: `input` has length n; returns the n/2+1 spectrum bins.
    /// Uses the classic pack-as-complex + untangle identity; verified against
    /// the embed-into-complex oracle in tests.
    #[must_use]
    pub fn forward(&self, input: &[f64]) -> Vec<C64> {
        let n = self.n;
        assert_eq!(
            input.len(),
            n,
            "input length must equal the planned size {n}"
        );
        let h = n / 2;
        // Pack even samples as re, odd as im.
        let mut z: Vec<C64> = (0..h)
            .map(|j| C64::new(input[2 * j], input[2 * j + 1]))
            .collect();
        let mut scratch = vec![C64::default(); h];
        self.half.forward(&mut z, &mut scratch);
        // Untangle: X[k] = (Z[k] + conj(Z[h−k]))/2 − i·w2n^k·(Z[k] − conj(Z[h−k]))/2.
        let mut out = vec![C64::default(); h + 1];
        for k in 0..=h {
            let zk = if k == h { z[0] } else { z[k] };
            let zc = if k == 0 { z[0].conj() } else { z[h - k].conj() };
            let even = zk.add(zc).scale(0.5);
            let odd = zk.sub(zc).scale(0.5);
            let theta = -std::f64::consts::PI * (k as f64) / (h as f64);
            let w = C64::new(det::cos(theta), det::sin(theta));
            // −i·w·odd  ==  rotate odd by w then by −i.
            let rot = w.mul(odd);
            let minus_i_rot = C64::new(rot.im, -rot.re);
            out[k] = even.add(minus_i_rot);
        }
        out
    }

    /// Inverse c2r: takes the `n/2 + 1` non-redundant spectrum bins of a real
    /// signal and reconstructs the `n` real samples — the exact algebraic
    /// inverse of [`RealFft::forward`]. Hermitian symmetry is ASSUMED (standard
    /// c2r contract: only the returned real part is meaningful when the input
    /// spectrum is not conjugate-symmetric). Verified by r2c→c2r round-trip and
    /// against the full-size complex IFFT of the Hermitian-completed spectrum.
    #[must_use]
    pub fn inverse(&self, spectrum: &[C64]) -> Vec<f64> {
        let n = self.n;
        let h = n / 2;
        assert_eq!(
            spectrum.len(),
            h + 1,
            "c2r spectrum length must be n/2+1 = {}",
            h + 1
        );
        // Undo the untangle by solving the 2×2 system relating (Z[k], Z[h−k])
        // to (X[k], conj(X[h−k])). With w = w_k = exp(−iπk/h),
        // p = (1 − i·w)/2, q = (1 + i·w)/2, the system is
        //   X[k]         = p·Z[k] + q·conj(Z[h−k])
        //   conj(X[h−k]) = q·Z[k] + p·conj(Z[h−k])
        // whose determinant is D = p² − q² = −i·w, giving
        //   Z[k] = (p·X[k] − q·conj(X[h−k])) / D.
        // (X[h−k] uses X[h], the Nyquist bin, when k = 0.)
        let mut z = vec![C64::default(); h];
        for (k, zk) in z.iter_mut().enumerate() {
            let xk = spectrum[k];
            let xhk = if k == 0 { spectrum[h] } else { spectrum[h - k] };
            let theta = -std::f64::consts::PI * (k as f64) / (h as f64);
            let w = C64::new(det::cos(theta), det::sin(theta));
            // i·w = (−w.im, w.re), so (1 ∓ i·w)/2 are:
            let p = C64::new(f64::midpoint(1.0, w.im), -0.5 * w.re);
            let q = C64::new(f64::midpoint(1.0, -w.im), 0.5 * w.re);
            let d = C64::new(w.im, -w.re); // −i·w
            let num = p.mul(xk).sub(q.mul(xhk.conj()));
            // Complex divide by d: num·conj(d)/|d|² (robust to |w| ≠ 1 exactly).
            *zk = num.mul(d.conj()).scale(1.0 / d.norm_sq());
        }
        let mut scratch = vec![C64::default(); h];
        self.half.inverse(&mut z, &mut scratch);
        // Unpack the packed half-size sequence: even samples in re, odd in im.
        let mut out = vec![0.0; n];
        for (j, zj) in z.iter().enumerate() {
            out[2 * j] = zj.re;
            out[2 * j + 1] = zj.im;
        }
        out
    }
}

// ---------------------------------------------------------------------------
// DCT-II / DCT-III via FFT folding (the fs-cheb path: Chebyshev transforms
// are DCTs).
// ---------------------------------------------------------------------------

/// DCT-II (unnormalized): X[k] = Σ_j x[j]·cos(π k (2j+1) / (2n)).
/// Computed by even/odd folding + one complex FFT; verified against the
/// naive O(n²) definition in tests.
#[must_use]
pub fn dct2(input: &[f64]) -> Vec<f64> {
    let n = input.len();
    assert!(
        n >= 1 && n.is_power_of_two(),
        "DCT size {n} must be a power of two"
    );
    // Fold: v[j] = x[2j], v[n−1−j] = x[2j+1] for j < n/2.
    let mut v = vec![C64::default(); n];
    for j in 0..n.div_ceil(2) {
        v[j] = C64::new(input[2 * j], 0.0);
        if 2 * j + 1 < n {
            v[n - 1 - j] = C64::new(input[2 * j + 1], 0.0);
        }
    }
    let plan = Fft::new(n);
    let mut scratch = vec![C64::default(); n];
    plan.forward(&mut v, &mut scratch);
    // X[k] = Re( exp(−iπk/(2n)) · V[k] ).
    let mut out = vec![0.0; n];
    for (k, slot) in out.iter_mut().enumerate() {
        let theta = -std::f64::consts::PI * (k as f64) / (2.0 * n as f64);
        let w = C64::new(det::cos(theta), det::sin(theta));
        *slot = w.mul(v[k]).re;
    }
    out
}

/// DCT-III (the inverse of DCT-II up to 2/n): computed here directly from
/// the DCT-II inverse relation so `dct3(dct2(x)) · 2/n = x` (with the k=0
/// halving convention). Verified by round-trip and against the naive sum.
#[must_use]
pub fn dct3(input: &[f64]) -> Vec<f64> {
    let n = input.len();
    assert!(
        n >= 1 && n.is_power_of_two(),
        "DCT size {n} must be a power of two"
    );
    // Naive-fold inverse via the r2c-style identity is fiddly; v1 uses the
    // O(n log n) route through a full complex FFT of the twiddle-extended
    // sequence: U[j] = 0.5·x[0] + Σ_{k≥1} x[k]·cos(πk(2j+1)/(2n)) — build the
    // complex sequence c[k] = x[k]·exp(−iπk/(2n)) (k=0..n−1), FFT, and read
    // the even/odd unfold. Verified against the naive O(n²) oracle.
    let mut c = vec![C64::default(); n];
    for k in 0..n {
        let theta = -std::f64::consts::PI * (k as f64) / (2.0 * n as f64);
        let w = C64::new(det::cos(theta), det::sin(theta));
        let coeff = if k == 0 { 0.5 * input[0] } else { input[k] };
        c[k] = w.scale(coeff);
    }
    let plan = Fft::new(n);
    let mut scratch = vec![C64::default(); n];
    plan.forward(&mut c, &mut scratch);
    // Unfold: y[2j] = Re c_fft[j], y[2j+1] = Re c_fft[n−1−j] (conjugate-path).
    let mut out = vec![0.0; n];
    for j in 0..n.div_ceil(2) {
        out[2 * j] = c[j].re;
        if 2 * j + 1 < n {
            out[2 * j + 1] = c[n - 1 - j].re;
        }
    }
    out
}

// ---------------------------------------------------------------------------
// Multi-dimensional (2D/3D/N-D) transforms via SEPARABLE pencil decomposition:
// a Fourier transform is separable, so an N-D transform is exactly the 1D
// transform applied along each axis in turn (the row–column algorithm). Each
// axis reuses the planned radix-2 `Fft`. This is the correctness-first
// dimensionality layer; cache-blocked transposes, executor tiling, and the
// roofline gate remain the recorded perf follow-up (bead fs-fft-perf-multidim).
// ---------------------------------------------------------------------------

/// A planned N-dimensional complex FFT over a ROW-MAJOR buffer (last axis
/// contiguous). Every axis length must be a power of two ≥ 1. Holds one `Fft`
/// plan per axis; deterministic by construction (fixed axis order, fixed
/// gather/scatter order — the 1D determinism lifts to N-D).
#[derive(Debug, Clone)]
pub struct FftNd {
    dims: Vec<usize>,
    plans: Vec<Fft>,
}

impl FftNd {
    /// Plan an N-D transform over `dims` (row-major; each a power of two ≥ 1).
    ///
    /// # Panics
    /// If `dims` is empty, or any axis is not a power of two (via [`Fft::new`]).
    #[must_use]
    pub fn new(dims: &[usize]) -> FftNd {
        assert!(!dims.is_empty(), "FftNd needs at least one axis");
        let plans = dims.iter().map(|&d| Fft::new(d)).collect();
        FftNd {
            dims: dims.to_vec(),
            plans,
        }
    }

    /// The axis lengths, row-major (last axis contiguous).
    #[must_use]
    pub fn shape(&self) -> &[usize] {
        &self.dims
    }

    /// Total element count (the product of the axis lengths) — the required
    /// length of every `data` buffer.
    #[must_use]
    pub fn total(&self) -> usize {
        self.dims.iter().product()
    }

    /// Forward N-D DFT, unnormalized: the separable product of per-axis
    /// unnormalized 1D DFTs.
    pub fn forward(&self, data: &mut [C64]) {
        self.run(data, false);
    }

    /// Inverse N-D DFT with `1/total` normalization: applying the 1/n-scaled 1D
    /// inverse along each axis composes to exactly `1/∏ n_ax`, so
    /// `inverse(forward(x)) = x`.
    pub fn inverse(&self, data: &mut [C64]) {
        self.run(data, true);
    }

    fn run(&self, data: &mut [C64], inverse: bool) {
        let total = self.total();
        assert_eq!(
            data.len(),
            total,
            "buffer length {} must equal the product of dims {total}",
            data.len()
        );
        for (ax, plan) in self.plans.iter().enumerate() {
            let n = self.dims[ax];
            if n == 1 {
                continue; // a length-1 axis is the identity; skip the gather.
            }
            // Row-major stride of this axis = product of the trailing dims;
            // `outer` counts the leading-index combinations.
            let stride: usize = self.dims[ax + 1..].iter().product();
            let outer: usize = self.dims[..ax].iter().product();
            let mut line = vec![C64::default(); n];
            let mut scratch = vec![C64::default(); n];
            for o in 0..outer {
                for i in 0..stride {
                    let base = o * n * stride + i;
                    // Gather the pencil (n samples at `stride`), transform,
                    // scatter back — the pencil is contiguous in index space.
                    for (t, slot) in line.iter_mut().enumerate() {
                        *slot = data[base + t * stride];
                    }
                    if inverse {
                        plan.inverse(&mut line, &mut scratch);
                    } else {
                        plan.forward(&mut line, &mut scratch);
                    }
                    for (t, &v) in line.iter().enumerate() {
                        data[base + t * stride] = v;
                    }
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn lcg(seed: &mut u64) -> f64 {
        *seed = seed
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        ((*seed >> 11) as f64) / (1u64 << 53) as f64 - 0.5
    }

    /// The definitive oracle: naive O(n²) DFT with strict-mode twiddles.
    fn naive_dft(x: &[C64], inverse: bool) -> Vec<C64> {
        let n = x.len();
        let sign = if inverse { 1.0 } else { -1.0 };
        let mut out = vec![C64::default(); n];
        for (k, slot) in out.iter_mut().enumerate() {
            let mut acc = C64::default();
            for (j, &xj) in x.iter().enumerate() {
                let theta = sign * 2.0 * std::f64::consts::PI * ((j * k) % n) as f64 / n as f64;
                acc = acc.add(xj.mul(C64::new(det::cos(theta), det::sin(theta))));
            }
            *slot = if inverse {
                acc.scale(1.0 / n as f64)
            } else {
                acc
            };
        }
        out
    }

    fn max_rel_err(a: &[C64], b: &[C64]) -> f64 {
        let scale = b
            .iter()
            .map(|v| v.norm_sq().sqrt())
            .fold(1e-300, f64::max)
            .max(1e-12);
        a.iter()
            .zip(b)
            .map(|(x, y)| x.sub(*y).norm_sq().sqrt() / scale)
            .fold(0.0, f64::max)
    }

    #[test]
    fn matches_naive_dft_oracle_across_sizes() {
        let mut seed = 0xFF7_u64;
        for log_n in 0..=9 {
            let n = 1usize << log_n;
            let plan = Fft::new(n);
            let x: Vec<C64> = (0..n)
                .map(|_| C64::new(lcg(&mut seed), lcg(&mut seed)))
                .collect();
            let mut data = x.clone();
            let mut scratch = vec![C64::default(); n];
            plan.forward(&mut data, &mut scratch);
            let want = naive_dft(&x, false);
            let err = max_rel_err(&data, &want);
            assert!(
                err < 1e-12,
                "n={n}: forward deviates from oracle by {err:.2e}"
            );
            // Inverse against the oracle too.
            let mut back = want.clone();
            plan.inverse(&mut back, &mut scratch);
            let err_inv = max_rel_err(&back, &x);
            assert!(err_inv < 1e-12, "n={n}: inverse deviates by {err_inv:.2e}");
        }
        println!(
            "{{\"suite\":\"fs-fft\",\"case\":\"oracle\",\"verdict\":\"pass\",\"detail\":\"n=1..512 vs naive DFT, rel<1e-12\"}}"
        );
    }

    #[test]
    fn impulse_constant_and_linearity() {
        let n = 64;
        let plan = Fft::new(n);
        let mut scratch = vec![C64::default(); n];
        // Impulse at 0 → flat ones.
        let mut x = vec![C64::default(); n];
        x[0] = C64::new(1.0, 0.0);
        plan.forward(&mut x, &mut scratch);
        for (k, v) in x.iter().enumerate() {
            assert!(
                (v.re - 1.0).abs() < 1e-14 && v.im.abs() < 1e-14,
                "impulse spectrum bin {k}: {v:?}"
            );
        }
        // Constant → n·δ₀.
        let mut c = vec![C64::new(1.0, 0.0); n];
        plan.forward(&mut c, &mut scratch);
        assert!((c[0].re - n as f64).abs() < 1e-12);
        for v in &c[1..] {
            assert!(v.norm_sq().sqrt() < 1e-12);
        }
        // Linearity: F(a·x + y) = a·F(x) + F(y).
        let mut seed = 42u64;
        let xv: Vec<C64> = (0..n)
            .map(|_| C64::new(lcg(&mut seed), lcg(&mut seed)))
            .collect();
        let yv: Vec<C64> = (0..n)
            .map(|_| C64::new(lcg(&mut seed), lcg(&mut seed)))
            .collect();
        let a = 2.75;
        let mut lhs: Vec<C64> = xv
            .iter()
            .zip(&yv)
            .map(|(x, y)| x.scale(a).add(*y))
            .collect();
        plan.forward(&mut lhs, &mut scratch);
        let mut fx = xv.clone();
        plan.forward(&mut fx, &mut scratch);
        let mut fy = yv.clone();
        plan.forward(&mut fy, &mut scratch);
        let rhs: Vec<C64> = fx
            .iter()
            .zip(&fy)
            .map(|(x, y)| x.scale(a).add(*y))
            .collect();
        assert!(max_rel_err(&lhs, &rhs) < 1e-13, "linearity violated");
    }

    #[test]
    fn parseval_and_shift_theorem() {
        let n = 256;
        let plan = Fft::new(n);
        let mut scratch = vec![C64::default(); n];
        let mut seed = 7u64;
        let x: Vec<C64> = (0..n)
            .map(|_| C64::new(lcg(&mut seed), lcg(&mut seed)))
            .collect();
        let mut fx = x.clone();
        plan.forward(&mut fx, &mut scratch);
        // Parseval: Σ|x|² = (1/n)·Σ|X|².
        let time: f64 = x.iter().map(|v| v.norm_sq()).sum();
        let freq: f64 = fx.iter().map(|v| v.norm_sq()).sum::<f64>() / n as f64;
        assert!(
            (time - freq).abs() / time < 1e-13,
            "Parseval: {time} vs {freq}"
        );
        // Shift: F(x rotated by s)[k] = X[k]·exp(−2πiks/n).
        let s = 5usize;
        let mut shifted: Vec<C64> = (0..n).map(|j| x[(j + s) % n]).collect();
        plan.forward(&mut shifted, &mut scratch);
        for (k, v) in shifted.iter().enumerate() {
            let theta = 2.0 * std::f64::consts::PI * ((k * s) % n) as f64 / n as f64;
            let w = C64::new(det::cos(theta), det::sin(theta));
            let want = fx[k].mul(w);
            assert!(
                v.sub(want).norm_sq().sqrt() < 1e-10,
                "shift theorem failed at bin {k}"
            );
        }
        println!(
            "{{\"suite\":\"fs-fft\",\"case\":\"laws\",\"verdict\":\"pass\",\"detail\":\"Parseval + shift theorem at n={n}\"}}"
        );
    }

    #[test]
    fn real_fft_matches_embedded_complex_oracle() {
        let mut seed = 11u64;
        for n in [4usize, 8, 64, 256] {
            let rfft = RealFft::new(n);
            let x: Vec<f64> = (0..n).map(|_| lcg(&mut seed)).collect();
            let got = rfft.forward(&x);
            // Oracle: embed into complex and take the first n/2+1 bins.
            let plan = Fft::new(n);
            let mut z: Vec<C64> = x.iter().map(|&v| C64::new(v, 0.0)).collect();
            let mut scratch = vec![C64::default(); n];
            plan.forward(&mut z, &mut scratch);
            for k in 0..=n / 2 {
                let d = got[k].sub(z[k]).norm_sq().sqrt();
                assert!(
                    d < 1e-11,
                    "r2c bin {k} of n={n}: {:?} vs {:?}",
                    got[k],
                    z[k]
                );
            }
        }
        println!(
            "{{\"suite\":\"fs-fft\",\"case\":\"r2c\",\"verdict\":\"pass\",\"detail\":\"packed r2c == embedded oracle, n in 4..256\"}}"
        );
    }

    #[test]
    fn dct2_and_dct3_match_naive_and_round_trip() {
        let mut seed = 13u64;
        for n in [2usize, 8, 32, 128] {
            let x: Vec<f64> = (0..n).map(|_| lcg(&mut seed)).collect();
            let got = dct2(&x);
            // Naive DCT-II.
            for (k, &g) in got.iter().enumerate() {
                let want: f64 = x
                    .iter()
                    .enumerate()
                    .map(|(j, &v)| {
                        v * det::cos(
                            std::f64::consts::PI * k as f64 * (2 * j + 1) as f64 / (2.0 * n as f64),
                        )
                    })
                    .sum();
                assert!(
                    (g - want).abs() < 1e-11 * (n as f64),
                    "DCT-II[{k}] n={n}: {g} vs {want}"
                );
            }
            // DCT-III naive check.
            let y = dct3(&got);
            for (j, &v) in y.iter().enumerate() {
                let want: f64 = 0.5 * got[0]
                    + (1..n)
                        .map(|k| {
                            got[k]
                                * det::cos(
                                    std::f64::consts::PI * k as f64 * (2 * j + 1) as f64
                                        / (2.0 * n as f64),
                                )
                        })
                        .sum::<f64>();
                assert!((v - want).abs() < 1e-10 * n as f64, "DCT-III[{j}] n={n}");
            }
            // Round trip: dct3(dct2(x)) · 2/n = x.
            for (j, (&yj, &xj)) in y.iter().zip(&x).enumerate() {
                assert!(
                    (yj * 2.0 / n as f64 - xj).abs() < 1e-11,
                    "DCT round-trip at {j}, n={n}"
                );
            }
        }
        println!(
            "{{\"suite\":\"fs-fft\",\"case\":\"dct\",\"verdict\":\"pass\",\"detail\":\"DCT-II/III vs naive + round-trip, n in 2..128\"}}"
        );
    }

    #[test]
    fn deterministic_and_golden_hash() {
        // Same input → same bits; plus a cross-ISA golden hash over a fixed
        // battery (the fs-math determinism story extended to transforms).
        let n = 128;
        let plan = Fft::new(n);
        let run = || {
            let mut scratch = vec![C64::default(); n];
            let mut seed = 0xD_15C_u64;
            let mut acc: u64 = 0xcbf2_9ce4_8422_2325;
            let mut feed = |v: f64| {
                for b in v.to_bits().to_le_bytes() {
                    acc ^= u64::from(b);
                    acc = acc.wrapping_mul(0x0000_0100_0000_01b3);
                }
            };
            for _ in 0..16 {
                let mut x: Vec<C64> = (0..n)
                    .map(|_| C64::new(lcg(&mut seed), lcg(&mut seed)))
                    .collect();
                plan.forward(&mut x, &mut scratch);
                for v in &x {
                    feed(v.re);
                    feed(v.im);
                }
            }
            acc
        };
        let h1 = run();
        let h2 = run();
        assert_eq!(h1, h2, "same input must produce identical bits");
        println!(
            "{{\"suite\":\"fs-fft\",\"case\":\"golden-hash\",\"verdict\":\"info\",\"detail\":\"{h1:#018x}\"}}"
        );
        assert_eq!(
            h1, GOLDEN_HASH,
            "FFT output bits changed: {h1:#018x} vs {GOLDEN_HASH:#018x} — bump only with \
             semantic justification (golden-evidence policy)"
        );
    }

    /// Recorded on aarch64-apple (M4 Pro); verified identical on x86-64
    /// (Threadripper, trj) — the cross-ISA bit-determinism evidence.
    const GOLDEN_HASH: u64 = 0xbd55_68d2_33f4_b4bc;

    #[test]
    fn non_power_of_two_is_refused_loudly() {
        for bad in [3usize, 6, 100] {
            let r = std::panic::catch_unwind(|| Fft::new(bad));
            assert!(r.is_err(), "size {bad} must be refused");
        }
    }

    #[test]
    fn real_fft_c2r_round_trips_and_matches_full_ifft() {
        let mut seed = 0x5C2_u64;
        for n in [2usize, 4, 8, 64, 256] {
            let rfft = RealFft::new(n);
            let x: Vec<f64> = (0..n).map(|_| lcg(&mut seed)).collect();
            let spectrum = rfft.forward(&x);
            assert_eq!(spectrum.len(), n / 2 + 1);
            // r2c → c2r is the identity to fp precision.
            let back = rfft.inverse(&spectrum);
            for (j, (&b, &xj)) in back.iter().zip(&x).enumerate() {
                assert!(
                    (b - xj).abs() < 1e-11,
                    "c2r round-trip n={n} idx {j}: {b} vs {xj}"
                );
            }
            // Independent oracle: Hermitian-complete the half spectrum to full
            // length, run the ordinary complex inverse, take the real part.
            let h = n / 2;
            let mut full = vec![C64::default(); n];
            for (k, &xk) in spectrum.iter().enumerate() {
                full[k] = xk;
            }
            for k in 1..h {
                full[n - k] = spectrum[k].conj();
            }
            let plan = Fft::new(n);
            let mut scratch = vec![C64::default(); n];
            plan.inverse(&mut full, &mut scratch);
            for (j, (&b, f)) in back.iter().zip(&full).enumerate() {
                assert!(
                    (b - f.re).abs() < 1e-11 && f.im.abs() < 1e-9,
                    "c2r vs full IFFT n={n} idx {j}: {b} vs {f:?}"
                );
            }
        }
        println!(
            "{{\"suite\":\"fs-fft\",\"case\":\"c2r\",\"verdict\":\"pass\",\"detail\":\"r2c->c2r round-trip + full-IFFT oracle, n in 2..256\"}}"
        );
    }

    /// Fully independent N-D DFT oracle: sums over every source index with the
    /// exact separable phase — shares NO gather/scatter structure with the
    /// Stockham pencil path in [`FftNd`].
    fn naive_dft_nd(x: &[C64], dims: &[usize]) -> Vec<C64> {
        let d = dims.len();
        let total: usize = dims.iter().product();
        let mut strides = vec![1usize; d];
        let mut acc_stride = 1usize;
        for ax in (0..d).rev() {
            strides[ax] = acc_stride;
            acc_stride *= dims[ax];
        }
        let decode = |mut idx: usize| -> Vec<usize> {
            let mut c = vec![0usize; d];
            for ax in 0..d {
                c[ax] = idx / strides[ax];
                idx %= strides[ax];
            }
            c
        };
        let mut out = vec![C64::default(); total];
        for (ko, slot) in out.iter_mut().enumerate() {
            let k = decode(ko);
            let mut acc = C64::default();
            for (jo, &xj) in x.iter().enumerate() {
                let j = decode(jo);
                let mut frac = 0.0;
                for ax in 0..d {
                    frac += ((k[ax] * j[ax]) % dims[ax]) as f64 / dims[ax] as f64;
                }
                let phase = -2.0 * std::f64::consts::PI * frac;
                acc = acc.add(xj.mul(C64::new(det::cos(phase), det::sin(phase))));
            }
            *slot = acc;
        }
        out
    }

    #[test]
    fn fft_nd_matches_naive_oracle_2d_and_3d() {
        let mut seed = 0x3D_u64;
        for dims in [
            vec![4usize, 4],
            vec![8, 4],
            vec![2, 8],
            vec![4, 4, 4],
            vec![2, 4, 8],
        ] {
            let total: usize = dims.iter().product();
            let plan = FftNd::new(&dims);
            assert_eq!(plan.total(), total);
            assert_eq!(plan.shape(), &dims[..]);
            let x: Vec<C64> = (0..total)
                .map(|_| C64::new(lcg(&mut seed), lcg(&mut seed)))
                .collect();
            let mut data = x.clone();
            plan.forward(&mut data);
            let want = naive_dft_nd(&x, &dims);
            let err = max_rel_err(&data, &want);
            assert!(
                err < 1e-12,
                "N-D forward dims={dims:?} deviates by {err:.2e}"
            );
            // Round trip.
            plan.inverse(&mut data);
            let err_rt = max_rel_err(&data, &x);
            assert!(
                err_rt < 1e-12,
                "N-D round-trip dims={dims:?} err {err_rt:.2e}"
            );
        }
        println!(
            "{{\"suite\":\"fs-fft\",\"case\":\"nd-oracle\",\"verdict\":\"pass\",\"detail\":\"2D/3D pencil == naive N-D DFT + round-trip\"}}"
        );
    }

    #[test]
    fn fft_nd_separability_and_parseval() {
        // Separability: a 2D transform equals a 1D FFT along every row followed
        // by a 1D FFT along every column (built from the raw `Fft`, independent
        // of FftNd's internal axis loop).
        let (n0, n1) = (8usize, 4usize);
        let mut seed = 0x5E_u64;
        let x: Vec<C64> = (0..n0 * n1)
            .map(|_| C64::new(lcg(&mut seed), lcg(&mut seed)))
            .collect();
        let mut nd = x.clone();
        FftNd::new(&[n0, n1]).forward(&mut nd);

        // Manual row-then-column.
        let row = Fft::new(n1);
        let col = Fft::new(n0);
        let mut man = x.clone();
        let mut sc1 = vec![C64::default(); n1];
        for r in 0..n0 {
            let mut line: Vec<C64> = (0..n1).map(|c| man[r * n1 + c]).collect();
            row.forward(&mut line, &mut sc1);
            for (c, &v) in line.iter().enumerate() {
                man[r * n1 + c] = v;
            }
        }
        let mut sc0 = vec![C64::default(); n0];
        for c in 0..n1 {
            let mut line: Vec<C64> = (0..n0).map(|r| man[r * n1 + c]).collect();
            col.forward(&mut line, &mut sc0);
            for (r, &v) in line.iter().enumerate() {
                man[r * n1 + c] = v;
            }
        }
        assert!(
            max_rel_err(&nd, &man) < 1e-13,
            "N-D separability (row-then-column) violated"
        );

        // Parseval in N-D: Σ|x|² = (1/total)·Σ|X|².
        let time: f64 = x.iter().map(|v| v.norm_sq()).sum();
        let freq: f64 = nd.iter().map(|v| v.norm_sq()).sum::<f64>() / (n0 * n1) as f64;
        assert!(
            (time - freq).abs() / time < 1e-13,
            "N-D Parseval: {time} vs {freq}"
        );
    }

    #[test]
    fn fft_nd_convolution_theorem_2d() {
        // ifft(fft(a) ⊙ fft(b)) == circular 2D convolution of a and b — the
        // capability's headline use (spectral convolution / PDE stencils).
        let (n0, n1) = (4usize, 8usize);
        let total = n0 * n1;
        let mut seed = 0xC0_u64;
        let a: Vec<C64> = (0..total)
            .map(|_| C64::new(lcg(&mut seed), lcg(&mut seed)))
            .collect();
        let b: Vec<C64> = (0..total)
            .map(|_| C64::new(lcg(&mut seed), lcg(&mut seed)))
            .collect();
        let plan = FftNd::new(&[n0, n1]);
        let mut fa = a.clone();
        let mut fb = b.clone();
        plan.forward(&mut fa);
        plan.forward(&mut fb);
        let mut prod: Vec<C64> = fa.iter().zip(&fb).map(|(x, y)| x.mul(*y)).collect();
        plan.inverse(&mut prod);
        // Direct circular convolution.
        let mut want = vec![C64::default(); total];
        for k0 in 0..n0 {
            for k1 in 0..n1 {
                let mut acc = C64::default();
                for m0 in 0..n0 {
                    for m1 in 0..n1 {
                        let s0 = (k0 + n0 - m0) % n0;
                        let s1 = (k1 + n1 - m1) % n1;
                        acc = acc.add(a[m0 * n1 + m1].mul(b[s0 * n1 + s1]));
                    }
                }
                want[k0 * n1 + k1] = acc;
            }
        }
        assert!(
            max_rel_err(&prod, &want) < 1e-12,
            "2D circular convolution theorem violated"
        );
        println!(
            "{{\"suite\":\"fs-fft\",\"case\":\"nd-conv\",\"verdict\":\"pass\",\"detail\":\"2D convolution theorem holds\"}}"
        );
    }

    #[test]
    fn fft_nd_is_deterministic() {
        let dims = [4usize, 8, 2];
        let total: usize = dims.iter().product();
        let run = || {
            let mut seed = 0xDE7_u64;
            let mut data: Vec<C64> = (0..total)
                .map(|_| C64::new(lcg(&mut seed), lcg(&mut seed)))
                .collect();
            FftNd::new(&dims).forward(&mut data);
            data
        };
        let a = run();
        let b = run();
        assert!(
            a.iter()
                .zip(&b)
                .all(|(x, y)| x.re.to_bits() == y.re.to_bits() && x.im.to_bits() == y.im.to_bits()),
            "N-D transform must be bitwise deterministic"
        );
    }

    #[test]
    fn version_is_stamped() {
        assert!(!VERSION.is_empty());
    }
}
