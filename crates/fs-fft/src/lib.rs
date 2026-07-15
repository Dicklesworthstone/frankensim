//! fs-fft — Stockham autosort FFTs, real transforms, and DCT (plan §6.3).
//!
//! The correctness-first core is Stockham autosort (no bit-reversal pass),
//! packed real transforms (r2c/c2r at half the complex work), and DCT-II/III
//! via FFT folding, which is the fs-cheb dependency. The production stage
//! walk is now mixed radix-8/4/2. The correctness oracle remains the naive
//! O(n²) DFT, run exhaustively over sizes and random inputs: the ORACLE
//! decides, not derivation confidence.
//!
//! Determinism: twiddles come from fs-math's STRICT sin/cos and every
//! butterfly runs in a fixed order. The complex stage path has a cross-ISA
//! golden. DCT-II/III have numerical oracles and same-build bit replay here;
//! downstream DCT bit changes are additionally tracked by registered fs-cheb
//! and vessel goldens. The current fs-cheb row is both-ISA verified; the
//! vessel row's post-radix-8 x86 reproduction remains pending.
//!
//! Bead fs-fft-perf-multidim extended the core with the r2c inverse
//! ([`RealFft::inverse`]), N-dimensional separable pencils ([`FftNd`]),
//! mixed radix-8/4/2 stages, and fs-simd capsules. The optional `[F]`
//! six-step path fuses its transpose/copy structure into two full-array
//! gather/scatter passes. It remains default-off after losing to the stage
//! walk on M4; its current x86 verdict, executor-tiled pencils, and
//! quiet-machine certification of the 40% target remain open.
//!
//! Conventions: forward is unnormalized; `inverse` scales by 1/n so
//! inverse(forward(x)) = x. Sizes must be powers of two (structured
//! rejection otherwise; mixed-radix general n is out of v1 scope).

use fs_math::det;

mod simd_view;

/// Crate version, re-exported for provenance stamping.
pub const VERSION: &str = env!("CARGO_PKG_VERSION");

/// Semantic version of the transform bit contract (golden-couplings
/// surface `fs-fft:transform-bits`): the stage decomposition (currently
/// mixed radix-8/4/2), the six-step dispatch predicate and its pass
/// structure, twiddle-application orders, and DCT folding/post-rotation.
/// Any change that moves output bits bumps this and deliberately re-freezes
/// the dependents in golden-couplings.json (docs/GOLDEN_POLICY.md).
pub const TRANSFORM_BIT_SEMANTICS_VERSION: u32 = 1;

/// Full-array memory passes performed by one fused six-step transform.
///
/// Stage A gathers `data` columns and scatters completed first sweeps into
/// `scratch`; stage B scatters completed second-sweep rows back to `data`.
/// Sub-transforms are cache-resident and therefore do not add full-array
/// traffic at the parent size. The roofline lane consumes this constant so
/// its byte model cannot silently retain the pre-fusion six-pass count.
#[doc(hidden)]
pub const SIXSTEP_FULL_ARRAY_PASSES: usize = 2;

/// Evidence identity for the fused two-pass six-step performance model.
///
/// Bump this whenever the six-step implementation or its traffic accounting
/// changes, even when exact-move rewrites leave transform output bits intact.
#[doc(hidden)]
pub const SIXSTEP_PERFORMANCE_MODEL_VERSION: &str = "27d3-6s-fused2";

/// A complex number (f64 re/im). Local, minimal: fs-la's future complex
/// types can absorb this when a shared home exists. `repr(C)` pins the
/// (re, im) interleaved layout the fs-simd stage kernels view as f64.
#[repr(C)]
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
    /// Square sub-plan (size √n) for the six-step path; built exactly when
    /// [`Fft::takes_sixstep`] holds for `n` (recursion terminates: √n < n
    /// for every n ≥ 2).
    sub: Option<Box<Fft>>,
}

/// Six-step engagement threshold (bead 27d3): below this the whole working
/// set is cache-resident and the direct stage walk wins; at and above it
/// the stage walk streams the full array from DRAM every pass, while the
/// fused six-step keeps its sub-transforms cache-resident and makes two
/// full-array gather/scatter passes.
/// Part of the bit contract via dispatch (a pure function of n).
const SIXSTEP_MIN: usize = 1 << 16;

impl Fft {
    /// Does size `n` take the six-step path? A pure function of `n` for a
    /// given build: large enough, an even power of two (n₁ = n₂ = √n,
    /// square transposes), AND the `frontier-sixstep` feature — the path
    /// is correct and golden-frozen but remains SLOWER than the stage walk
    /// on M4 after two-pass fusion and vectorized strip moves (2026-07-11,
    /// n = 2²²: six-step 0.0822-0.0852 s vs stage walk 0.053-0.055 s).
    /// It therefore stays frontier/default-off on that ISA; the current x86
    /// verdict remains pending. Enabling the feature changes large-n output
    /// bits, which is why the six-step golden lives in the gated battery.
    #[doc(hidden)]
    #[must_use]
    pub fn takes_sixstep(n: usize) -> bool {
        cfg!(feature = "frontier-sixstep")
            && n.is_power_of_two()
            && n >= SIXSTEP_MIN
            && n.trailing_zeros().is_multiple_of(2)
    }

    /// Stage-walk entry for the gated conformance battery's cross-path
    /// check (the fs-math `erf_both_paths` pattern): same values as the
    /// dispatched path, different summation order.
    #[doc(hidden)]
    pub fn forward_via_stages(&self, data: &mut [C64], scratch: &mut [C64]) {
        assert_eq!(
            data.len(),
            self.n,
            "data length must equal the planned size"
        );
        assert_eq!(
            scratch.len(),
            self.n,
            "scratch length must equal the planned size"
        );
        self.transform_stages(data, scratch, false);
    }

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
        let sub = if Fft::takes_sixstep(n) {
            Some(Box::new(Fft::new(1 << (n.trailing_zeros() / 2))))
        } else {
            None
        };
        Fft { n, table, sub }
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

    /// Twiddle w^k = exp(−2πik/n) for k < 3n/4: the table stores the
    /// first half-turn; the half-turn symmetry w^(k+n/2) = −w^k extends
    /// it EXACTLY (negation is exact) for the radix-4 stages' 3·p·s
    /// indices.
    /// One radix-8 DIF Stockham stage: even outputs are the 4-point DFT
    /// of the half-sums, odd outputs the 4-point DFT of the half-diffs
    /// twisted by ω₈ʲ (ω₈² is the exact ∓i rotation; ω₈¹/ω₈³ use the
    /// exact `FRAC_1_SQRT_2` literal). One fused operation order per
    /// element — deterministic on every target.
    fn stage_radix8(&self, src: &[C64], dst: &mut [C64], m: usize, s: usize, inverse: bool) {
        let c = std::f64::consts::FRAC_1_SQRT_2;
        let (w8_1, w8_3) = if inverse {
            (C64::new(c, c), C64::new(-c, c))
        } else {
            (C64::new(c, -c), C64::new(-c, -c))
        };
        for p in 0..m {
            let mut w = [C64::default(); 7];
            for (k, slot) in w.iter_mut().enumerate() {
                *slot = self.tw((k + 1) * p * s);
                if inverse {
                    *slot = slot.conj();
                }
            }
            for q in 0..s {
                let a0 = src[q + s * p];
                let a1 = src[q + s * (p + m)];
                let a2 = src[q + s * (p + 2 * m)];
                let a3 = src[q + s * (p + 3 * m)];
                let a4 = src[q + s * (p + 4 * m)];
                let a5 = src[q + s * (p + 5 * m)];
                let a6 = src[q + s * (p + 6 * m)];
                let a7 = src[q + s * (p + 7 * m)];
                // Half split: even outputs = 4-point DFT of u, odd
                // outputs = 4-point DFT of v·ω₈ʲ.
                let u0 = a0.add(a4);
                let u1 = a1.add(a5);
                let u2 = a2.add(a6);
                let u3 = a3.add(a7);
                let v0 = a0.sub(a4);
                let v1 = a1.sub(a5).mul(w8_1);
                // v2·ω₈² = ∓i·v2 — the exact rotation, no multiply.
                let v2r = a2.sub(a6);
                let v2 = if inverse {
                    C64::new(-v2r.im, v2r.re)
                } else {
                    C64::new(v2r.im, -v2r.re)
                };
                let v3 = a3.sub(a7).mul(w8_3);
                // Even 4-point kernel (same ∓i definition as the radix-4
                // stage — one semantic definition).
                let e0 = u0.add(u2);
                let e1 = u0.sub(u2);
                let e2 = u1.add(u3);
                let e3 = u1.sub(u3);
                let e3i = if inverse {
                    C64::new(-e3.im, e3.re)
                } else {
                    C64::new(e3.im, -e3.re)
                };
                // Odd 4-point kernel.
                let f0 = v0.add(v2);
                let f1 = v0.sub(v2);
                let f2 = v1.add(v3);
                let f3 = v1.sub(v3);
                let f3i = if inverse {
                    C64::new(-f3.im, f3.re)
                } else {
                    C64::new(f3.im, -f3.re)
                };
                dst[q + s * 8 * p] = e0.add(e2);
                dst[q + s * (8 * p + 1)] = f0.add(f2).mul(w[0]);
                dst[q + s * (8 * p + 2)] = e1.add(e3i).mul(w[1]);
                dst[q + s * (8 * p + 3)] = f1.add(f3i).mul(w[2]);
                dst[q + s * (8 * p + 4)] = e0.sub(e2).mul(w[3]);
                dst[q + s * (8 * p + 5)] = f0.sub(f2).mul(w[4]);
                dst[q + s * (8 * p + 6)] = e1.sub(e3i).mul(w[5]);
                dst[q + s * (8 * p + 7)] = f1.sub(f3i).mul(w[6]);
            }
        }
    }

    fn tw(&self, k: usize) -> C64 {
        if k < self.table.len() {
            self.table[k]
        } else {
            let t = self.table[k - self.table.len()];
            C64::new(-t.re, -t.im)
        }
    }

    /// Transform dispatch — a PURE FUNCTION of the planned size (part of
    /// the bit contract): sizes passing [`Fft::takes_sixstep`] run the
    /// cache-blocked six-step; everything else runs the direct stage walk.
    fn transform(&self, data: &mut [C64], scratch: &mut [C64], inverse: bool) {
        let n = self.n;
        assert_eq!(data.len(), n, "data length must equal the planned size {n}");
        assert_eq!(
            scratch.len(),
            n,
            "scratch length must equal the planned size {n}"
        );
        if self.sub.is_some() {
            self.transform_sixstep(data, scratch, inverse);
        } else {
            self.transform_stages(data, scratch, inverse);
        }
    }

    /// Cache-blocked six-step FFT (bead 27d3): for n = n₁², view the array
    /// as an n₁×n₁ row-major matrix. Its logical three-transpose/two-sweep
    /// decomposition is fused into two full-array storage passes around
    /// cache-resident row transforms: `data` columns → first sweep + twiddle
    /// → `scratch` columns, then `scratch` rows → second sweep → `data`
    /// columns. Gather/scatter operations are exact element moves; twiddles
    /// come from the same strict table via `tw` (`k·j < n` by construction),
    /// and row sweeps reuse the mixed radix-8/4/2 kernels. Decomposition,
    /// sweep order, and twiddle indices are pure functions of n.
    fn transform_sixstep(&self, data: &mut [C64], scratch: &mut [C64], inverse: bool) {
        const GCOLS: usize = 8;

        let sub = self.sub.as_deref().expect("dispatch guaranteed a sub-plan");
        let n1 = sub.n;
        // FUSED six-step (bead 27d3 copy-back fusion): the three
        // transposes and the final copy-back are folded into
        // line-blocked gather/scatter around the two row-FFT sweeps —
        // TWO full-array memory passes instead of six. The arithmetic
        // is BIT-IDENTICAL to the unfused formulation: the same
        // sub-FFTs see the same values, the same twiddles apply in the
        // same per-element order, and only STORAGE locations moved
        // (the frozen sixstep golden pins exactly this).
        //
        // Line blocking: one gathered/scattered group is GCOLS = 8
        // complex = 128 bytes = one cache line, so the strided side of
        // each pass still reads/writes every line exactly once (a
        // column-at-a-time walk would touch each line 8x).
        debug_assert_eq!(n1 % GCOLS, 0, "sixstep sizes have n1 = 2^e >= 256");
        let mut row_scratch = vec![C64::default(); n1];
        let mut bufs = vec![C64::default(); GCOLS * n1];
        // STAGE A (was T1 + sweep1 + T2): per column group of `data`,
        // gather 8 columns (each line feeds all 8 buffers), sub-FFT +
        // fused twiddle each column j — tw(0) = 1 exactly, and C64::mul
        // by (1, 0) is exact, so no special case — then scatter into
        // the same 8 columns of `scratch`. A column is dead after its
        // own gather, so writing it back cannot alias live input.
        let gather = fs_simd::ops().gath8c64;
        let scatter = fs_simd::ops().scat8c64;
        let mut g = 0;
        while g < n1 {
            // Vectorized column-group gather (bead 27d3 final lever):
            // one q-register move per complex, bitwise vs the scalar
            // twin in fs-simd's tier battery.
            gather(
                simd_view::as_f64(data),
                simd_view::as_f64_mut(&mut bufs),
                n1,
                g,
            );
            for c in 0..GCOLS {
                let j = g + c;
                let buf = &mut bufs[c * n1..(c + 1) * n1];
                // Full dispatch (not transform_stages): a sub-plan large
                // enough to qualify recurses into its own six-step.
                sub.transform(buf, &mut row_scratch, inverse);
                for (k, v) in buf.iter_mut().enumerate() {
                    let mut w = self.tw(k * j);
                    if inverse {
                        w = w.conj();
                    }
                    *v = v.mul(w);
                }
            }
            scatter(
                simd_view::as_f64(&bufs),
                simd_view::as_f64_mut(scratch),
                n1,
                g,
            );
            g += GCOLS;
        }
        // STAGE B (was sweep2 + T3 + copy): rows of `scratch` are the
        // contiguous second-sweep inputs — sub-FFT each group of 8 rows
        // in place, then scatter those finished rows into 8 COLUMNS of
        // `data`. Everything `data` held was consumed by stage A's
        // gathers, so the scatter aliases nothing live, and the result
        // lands in `data` already output-transposed: no third
        // transpose, no copy-back.
        let mut r = 0;
        while r < n1 {
            for c in 0..GCOLS {
                let row = &mut scratch[(r + c) * n1..(r + c + 1) * n1];
                sub.transform(row, &mut row_scratch, inverse);
            }
            // Eight finished contiguous rows ARE a dense 8×n₁ buffer:
            // the same scatter primitive lands them as output columns.
            scatter(
                simd_view::as_f64(&scratch[r * n1..(r + GCOLS) * n1]),
                simd_view::as_f64_mut(data),
                n1,
                r,
            );
            r += GCOLS;
        }
    }

    /// Mixed radix-8/4/2 Stockham autosort (bead 27d3): ping-pongs between
    /// `data` and `scratch` with no bit-reversal pass; butterfly order is
    /// a pure function of the stage structure (deterministic by
    /// construction). Radix-8 stages consume three log₂ bits per
    /// full-array pass, radix-4 two, and one radix-4-or-2 residue absorbs
    /// the modulus. Radix changes legitimately changed twiddle-application
    /// order and hence output bits: the golden hash was bumped with that
    /// justification each time (see the golden test).
    #[allow(clippy::too_many_lines)] // the stage driver IS the decomposition
    fn transform_stages(&self, data: &mut [C64], scratch: &mut [C64], inverse: bool) {
        let n = self.n;
        if n == 1 {
            return;
        }
        let mut n_cur = n;
        let mut s = 1usize;
        let mut src_is_data = true;
        // The stage q-runs go through the fs-simd dispatch table for the
        // large-stride stages (bead 27d3 scope 2): the NEON kernel
        // deinterleaves two complex elements per iteration, the x86 AVX2
        // kernel four, and BOTH are BITWISE-identical to the scalar twin,
        // which is itself the inline loop below — one semantic definition,
        // tier-tested in fs-simd's battery, and this golden did NOT move
        // when either capsule path landed.
        let r4 = fs_simd::ops().r4qrun_f64;
        // Radix-8 stages first (bead 27d3 slice 3): each consumes THREE
        // log₂ bits per full-array pass — ceil(log₂n/3) passes instead of
        // the radix-4 formulation's ceil(log₂n/2). The AVX2/NEON q-run
        // finding showed the transform is BANDWIDTH bound, so pass count
        // is the lever; the ~⅓ traffic cut beats the capsuled radix-4
        // path (measured on the perf lane — see the golden-bump note).
        // The decomposition is a pure function of n (radix-8 while
        // n ≥ 8, then one radix-4 or radix-2 residue), so bits stay a
        // deterministic function of (n, input) on every target. The
        // odd-branch eighth-turn twiddles use the exact FRAC_1_SQRT_2
        // literal — identical on every conforming platform.
        while n_cur >= 8 {
            let m = n_cur / 8;
            {
                let (src, dst): (&[C64], &mut [C64]) = if src_is_data {
                    (&*data, &mut *scratch)
                } else {
                    (&*scratch, &mut *data)
                };
                self.stage_radix8(src, dst, m, s, inverse);
            }
            n_cur = m;
            s *= 8;
            src_is_data = !src_is_data;
        }
        while n_cur >= 4 {
            let m = n_cur / 4;
            {
                let (src, dst): (&[C64], &mut [C64]) = if src_is_data {
                    (&*data, &mut *scratch)
                } else {
                    (&*scratch, &mut *data)
                };
                for p in 0..m {
                    let (mut w1, mut w2, mut w3) =
                        (self.tw(p * s), self.tw(2 * p * s), self.tw(3 * p * s));
                    if inverse {
                        (w1, w2, w3) = (w1.conj(), w2.conj(), w3.conj());
                    }
                    // Threshold measured, not vibed: below s = 64 the
                    // per-p call overhead (slicing + dispatch) costs
                    // more than the NEON win on 8–256-f64 runs; the
                    // inline loop below autovectorizes to ~the same
                    // GB/s there. Bits identical on both paths.
                    if s >= 64 {
                        let wv = [w1.re, w1.im, w2.re, w2.im, w3.re, w3.im];
                        r4(
                            simd_view::as_f64(&src[s * p..s * p + s]),
                            simd_view::as_f64(&src[s * (p + m)..s * (p + m) + s]),
                            simd_view::as_f64(&src[s * (p + 2 * m)..s * (p + 2 * m) + s]),
                            simd_view::as_f64(&src[s * (p + 3 * m)..s * (p + 3 * m) + s]),
                            simd_view::as_f64_mut(&mut dst[s * 4 * p..s * 4 * (p + 1)]),
                            &wv,
                            inverse,
                        );
                        continue;
                    }
                    for q in 0..s {
                        let a = src[q + s * p];
                        let b = src[q + s * (p + m)];
                        let c = src[q + s * (p + 2 * m)];
                        let d = src[q + s * (p + 3 * m)];
                        let t0 = a.add(c);
                        let t1 = a.sub(c);
                        let t2 = b.add(d);
                        let t3 = b.sub(d);
                        // ∓i·t3: forward −i (DIF kernel), inverse +i.
                        let t3i = if inverse {
                            C64::new(-t3.im, t3.re)
                        } else {
                            C64::new(t3.im, -t3.re)
                        };
                        dst[q + s * 4 * p] = t0.add(t2);
                        dst[q + s * (4 * p + 1)] = t1.add(t3i).mul(w1);
                        dst[q + s * (4 * p + 2)] = t0.sub(t2).mul(w2);
                        dst[q + s * (4 * p + 3)] = t1.sub(t3i).mul(w3);
                    }
                }
            }
            n_cur = m;
            s *= 4;
            src_is_data = !src_is_data;
        }
        if n_cur == 2 {
            let (src, dst): (&[C64], &mut [C64]) = if src_is_data {
                (&*data, &mut *scratch)
            } else {
                (&*scratch, &mut *data)
            };
            let mut w = self.table[0];
            if inverse {
                w = w.conj();
            }
            for q in 0..s {
                let a = src[q];
                let b = src[q + s];
                dst[q] = a.add(b);
                dst[q + s] = a.sub(b).mul(w);
            }
            src_is_data = !src_is_data;
        }
        if !src_is_data {
            data.copy_from_slice(scratch);
        }
    }
}

/// Out-of-place cache-blocked transpose of an n₁×n₁ row-major complex
/// matrix: `dst[i·n1 + j] = src[j·n1 + i]`. Pure element moves — NO
/// floating-point arithmetic, so the pass is bit-neutral by construction.
/// 8×8-element tiles (8 C64 = 128 B = one Apple cache line) with
/// SEQUENTIAL writes per destination row chunk (the strided side is the
/// read, which prefetches better than strided read-modify-write swaps).
#[cfg(test)] // production transposes are fused into the sixstep passes (bead 27d3)
fn transpose_into(src: &[C64], dst: &mut [C64], n1: usize) {
    debug_assert_eq!(src.len(), n1 * n1, "transpose needs a square matrix");
    debug_assert_eq!(dst.len(), n1 * n1, "transpose needs a square dst");
    // The tiled move loop is the fs-simd trn1c64 capsule (bead 27d3):
    // one interleaved complex is one 128-bit vector, so the NEON tier
    // moves 16 bytes per instruction with no per-element bounds checks.
    // Pure exact moves — bit-neutral by construction, gated bitwise
    // against the scalar twin in fs-simd's tier battery.
    (fs_simd::ops().trn1c64)(simd_view::as_f64(src), simd_view::as_f64_mut(dst), n1);
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

// ---------------------------------------------------------------------------
// Executor-tiled N-D pencils (bead 27d3): the per-axis pencil loop is
// embarrassingly parallel over OUTER BLOCKS (contiguous n*stride slices,
// disjoint by construction), so those axes run on the TilePool with the
// GEMM-band Mutex-chunk pattern — per-pencil arithmetic and order are
// EXACTLY the serial path's, so output is bitwise identical across
// worker counts (P2 by construction, gated). The FIRST axis has a single
// outer block (outer == 1) and runs gate-checked serial in v1; the
// column-group row-locking design for it is recorded on the bead.
// ---------------------------------------------------------------------------

/// One parallel axis pass: each tile owns one outer block and transforms
/// the `stride` pencils inside it, in the serial path's exact order.
struct PencilBlockKernel<'a> {
    blocks: &'a [std::sync::Mutex<&'a mut [C64]>],
    plan: &'a Fft,
    n: usize,
    stride: usize,
    inverse: bool,
}

impl fs_exec::TileKernel for PencilBlockKernel<'_> {
    type Out = ();

    fn tiles(&self) -> fs_exec::TilePlan {
        fs_exec::TilePlan::new(
            "fs-fft/ndim-pencil-block-v1",
            u64::try_from(self.blocks.len()).expect("block count fits u64"),
        )
    }

    fn run(
        &self,
        tile: u64,
        cx: &fs_exec::Cx<'_>,
    ) -> core::ops::ControlFlow<fs_exec::Cancelled, ()> {
        let tile = usize::try_from(tile).expect("tile index fits usize");
        let mut block = self.blocks[tile].lock().expect("pencil block lock");
        let mut line = vec![C64::default(); self.n];
        let mut scratch = vec![C64::default(); self.n];
        for i in 0..self.stride {
            if cx.checkpoint().is_err() {
                return core::ops::ControlFlow::Break(fs_exec::Cancelled);
            }
            for (t, slot) in line.iter_mut().enumerate() {
                *slot = block[i + t * self.stride];
            }
            if self.inverse {
                self.plan.inverse(&mut line, &mut scratch);
            } else {
                self.plan.forward(&mut line, &mut scratch);
            }
            for (t, &v) in line.iter().enumerate() {
                block[i + t * self.stride] = v;
            }
        }
        core::ops::ControlFlow::Continue(())
    }
}

/// Axis-0 parallelism (single outer block): the (n x stride) matrix is
/// chunked into n contiguous ROW slices behind mutexes; each tile owns a
/// COLUMN GROUP, gathering its pencils row by row, transforming them in
/// cache, and scattering back. Every element is read and written by
/// exactly one tile, so lock order affects timing only — per-pencil
/// arithmetic and order match the serial path bit for bit.
struct PencilColumnKernel<'a> {
    rows: &'a [std::sync::Mutex<&'a mut [C64]>],
    plan: &'a Fft,
    n: usize,
    stride: usize,
    group: usize,
    inverse: bool,
}

impl fs_exec::TileKernel for PencilColumnKernel<'_> {
    type Out = ();

    fn tiles(&self) -> fs_exec::TilePlan {
        fs_exec::TilePlan::new(
            "fs-fft/ndim-pencil-column-v1",
            u64::try_from(self.stride.div_ceil(self.group)).expect("tile count fits u64"),
        )
    }

    fn run(
        &self,
        tile: u64,
        cx: &fs_exec::Cx<'_>,
    ) -> core::ops::ControlFlow<fs_exec::Cancelled, ()> {
        let tile = usize::try_from(tile).expect("tile index fits usize");
        let lo = tile * self.group;
        let g = self.group.min(self.stride - lo);
        // Column-major workspace: line c occupies [c*n, (c+1)*n).
        let mut lines = vec![C64::default(); g * self.n];
        let mut scratch = vec![C64::default(); self.n];
        for (t, row) in self.rows.iter().enumerate() {
            let row = row.lock().expect("pencil row lock");
            for c in 0..g {
                lines[c * self.n + t] = row[lo + c];
            }
        }
        for c in 0..g {
            if cx.checkpoint().is_err() {
                return core::ops::ControlFlow::Break(fs_exec::Cancelled);
            }
            let line = &mut lines[c * self.n..(c + 1) * self.n];
            if self.inverse {
                self.plan.inverse(line, &mut scratch);
            } else {
                self.plan.forward(line, &mut scratch);
            }
        }
        for (t, row) in self.rows.iter().enumerate() {
            let mut row = row.lock().expect("pencil row lock");
            for c in 0..g {
                row[lo + c] = lines[c * self.n + t];
            }
        }
        core::ops::ControlFlow::Continue(())
    }
}

/// Measured facts about one executor-tiled axis pass (bead 3f6c): the
/// pass geometry and its wall time, for small-kernel granularity
/// diagnosis. Measurements only — transform results never depend on
/// them (the P2 law: output is bitwise identical at every worker
/// count), so `wall_ns` is envelope-class exactly like
/// [`fs_exec::RunReport`]'s latency samples.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct NdPassReport {
    /// Axis index in `dims` order.
    pub axis: usize,
    /// `true` for an inverse-direction pass.
    pub inverse: bool,
    /// Executor kernel that ran the pass (the serial degenerate pass
    /// reports `"fs-fft/ndim-axis-serial-v1"`).
    pub kernel: &'static str,
    /// Pencil length (this axis's dimension).
    pub n: usize,
    /// Pencil stride (product of trailing dims).
    pub stride: usize,
    /// Outer block count (product of leading dims).
    pub outer: usize,
    /// Tiles the pass planned (1 for the serial degenerate pass).
    pub tiles: u64,
    /// Tiles completed (equals `tiles` on success).
    pub completed: u64,
    /// Workers the serving pool exposed to the pass.
    pub workers: usize,
    /// Wall time of the pass in ns (saturating; measurement only).
    pub wall_ns: u64,
}

impl FftNd {
    /// Executor-tiled forward N-D DFT: bitwise identical to
    /// [`FftNd::forward`] at every worker count (gated). See
    /// [`FftNd::run_pooled`] for the cancellation contract.
    ///
    /// # Errors
    /// [`fs_exec::RunError`] on cancellation or a contained executor
    /// failure.
    ///
    /// # Panics
    /// If `data.len()` differs from [`FftNd::total`].
    pub fn forward_pooled<P: fs_exec::KernelRunner>(
        &self,
        data: &mut [C64],
        pool: &P,
        gate: &fs_exec::CancelGate,
    ) -> Result<(), fs_exec::RunError> {
        self.run_pooled(data, false, pool, gate)
    }

    /// Executor-tiled inverse (1/total-normalized), bitwise identical to
    /// [`FftNd::inverse`] at every worker count.
    ///
    /// # Errors
    /// As [`FftNd::forward_pooled`].
    ///
    /// # Panics
    /// As [`FftNd::forward_pooled`].
    pub fn inverse_pooled<P: fs_exec::KernelRunner>(
        &self,
        data: &mut [C64],
        pool: &P,
        gate: &fs_exec::CancelGate,
    ) -> Result<(), fs_exec::RunError> {
        self.run_pooled(data, true, pool, gate)
    }

    /// As [`FftNd::forward_pooled`], additionally reporting each axis
    /// pass's geometry and wall time to `observe` (bead 3f6c
    /// diagnostics). The observer sees exactly the passes that ran, in
    /// execution order; it cannot affect results.
    ///
    /// # Errors
    /// As [`FftNd::forward_pooled`]. On cancellation the interrupted
    /// pass is still observed, with `completed < tiles`.
    ///
    /// # Panics
    /// As [`FftNd::forward_pooled`].
    pub fn forward_pooled_observed<P: fs_exec::KernelRunner>(
        &self,
        data: &mut [C64],
        pool: &P,
        gate: &fs_exec::CancelGate,
        observe: &mut dyn FnMut(NdPassReport),
    ) -> Result<(), fs_exec::RunError> {
        self.run_pooled_observed(data, false, pool, gate, observe)
    }

    /// As [`FftNd::inverse_pooled`] with per-pass observation; see
    /// [`FftNd::forward_pooled_observed`].
    ///
    /// # Errors
    /// As [`FftNd::inverse_pooled`].
    ///
    /// # Panics
    /// As [`FftNd::inverse_pooled`].
    pub fn inverse_pooled_observed<P: fs_exec::KernelRunner>(
        &self,
        data: &mut [C64],
        pool: &P,
        gate: &fs_exec::CancelGate,
        observe: &mut dyn FnMut(NdPassReport),
    ) -> Result<(), fs_exec::RunError> {
        self.run_pooled_observed(data, true, pool, gate, observe)
    }

    /// CANCELLATION CONTRACT: this is a scratch-transform API — on `Err`
    /// the buffer contents are UNSPECIFIED (some axes may be transformed,
    /// some pencils of the interrupted axis may be). Callers needing
    /// transactional output stage a copy first; the drained pool and the
    /// structured error are the guarantees, torn-freedom of `data` is not.
    fn run_pooled<P: fs_exec::KernelRunner>(
        &self,
        data: &mut [C64],
        inverse: bool,
        pool: &P,
        gate: &fs_exec::CancelGate,
    ) -> Result<(), fs_exec::RunError> {
        self.run_pooled_observed(data, inverse, pool, gate, &mut |_| {})
    }

    /// The pooled axis loop with per-pass observation. Timing wraps only
    /// the pass itself; the observer runs between passes, outside any
    /// measured window, and results never depend on it.
    #[allow(clippy::too_many_lines)]
    fn run_pooled_observed<P: fs_exec::KernelRunner>(
        &self,
        data: &mut [C64],
        inverse: bool,
        pool: &P,
        gate: &fs_exec::CancelGate,
        observe: &mut dyn FnMut(NdPassReport),
    ) -> Result<(), fs_exec::RunError> {
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
                continue;
            }
            let stride: usize = self.dims[ax + 1..].iter().product();
            let outer: usize = self.dims[..ax].iter().product();
            if outer == 1 && stride == 1 {
                // Degenerate single pencil (1-D or all-leading-1 shapes):
                // serial, gate-checked.
                if gate.is_requested() {
                    return Err(fs_exec::RunError::Cancelled {
                        kernel: "fs-fft/ndim-axis-serial-v1",
                        completed: 0,
                        total: 1,
                    });
                }
                let pass_start = std::time::Instant::now();
                let mut line = vec![C64::default(); n];
                let mut scratch = vec![C64::default(); n];
                line.copy_from_slice(data);
                if inverse {
                    plan.inverse(&mut line, &mut scratch);
                } else {
                    plan.forward(&mut line, &mut scratch);
                }
                data.copy_from_slice(&line);
                observe(NdPassReport {
                    axis: ax,
                    inverse,
                    kernel: "fs-fft/ndim-axis-serial-v1",
                    n,
                    stride,
                    outer,
                    tiles: 1,
                    completed: 1,
                    workers: pool.workers(),
                    wall_ns: saturating_ns(pass_start.elapsed()),
                });
                continue;
            }
            if outer == 1 {
                // Axis-0 column groups (row-locked): full parallelism on
                // the first axis without changing any element's bits.
                let rows: Vec<std::sync::Mutex<&mut [C64]>> =
                    data.chunks_mut(stride).map(std::sync::Mutex::new).collect();
                let group = stride.div_ceil(pool.workers().max(1)).clamp(1, 64);
                let kernel = PencilColumnKernel {
                    rows: &rows,
                    plan,
                    n,
                    stride,
                    group,
                    inverse,
                };
                let pass_start = std::time::Instant::now();
                let (outcome, report) = pool.run_with_gate(&kernel, gate);
                observe(NdPassReport {
                    axis: ax,
                    inverse,
                    kernel: report.kernel,
                    n,
                    stride,
                    outer,
                    tiles: report.total,
                    completed: report.completed,
                    workers: pool.workers(),
                    wall_ns: saturating_ns(pass_start.elapsed()),
                });
                outcome?;
                continue;
            }
            let blocks: Vec<std::sync::Mutex<&mut [C64]>> = data
                .chunks_mut(n * stride)
                .map(std::sync::Mutex::new)
                .collect();
            let kernel = PencilBlockKernel {
                blocks: &blocks,
                plan,
                n,
                stride,
                inverse,
            };
            let pass_start = std::time::Instant::now();
            let (outcome, report) = pool.run_with_gate(&kernel, gate);
            observe(NdPassReport {
                axis: ax,
                inverse,
                kernel: report.kernel,
                n,
                stride,
                outer,
                tiles: report.total,
                completed: report.completed,
                workers: pool.workers(),
                wall_ns: saturating_ns(pass_start.elapsed()),
            });
            outcome?;
        }
        Ok(())
    }
}

/// Envelope-class ns from a duration (saturating; measurement only).
fn saturating_ns(elapsed: std::time::Duration) -> u64 {
    u64::try_from(elapsed.as_nanos()).unwrap_or(u64::MAX)
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
            let got_replay = dct2(&x);
            assert!(
                got.iter()
                    .zip(&got_replay)
                    .all(|(a, b)| a.to_bits() == b.to_bits()),
                "DCT-II same-build replay changed bits at n={n}"
            );
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
            let y_replay = dct3(&got);
            assert!(
                y.iter()
                    .zip(&y_replay)
                    .all(|(a, b)| a.to_bits() == b.to_bits()),
                "DCT-III same-build replay changed bits at n={n}"
            );
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
            "{{\"suite\":\"fs-fft\",\"case\":\"dct\",\"verdict\":\"pass\",\"detail\":\"DCT-II/III vs naive + round-trip + same-build bit replay, n in 2..128\"}}"
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

    /// JUSTIFIED BUMP (bead 27d3, 2026-07-10, second): the transform
    /// moved to mixed radix-8/4/2 Stockham — radix-8 stages consume
    /// three log₂ bits per full-array pass (ceil(log₂n/3) passes), the
    /// bandwidth-bound lever the AVX2 finding pointed at. Twiddle
    /// application order changes with the radix, which legitimately
    /// changes output bits (pre-authorized by the bead). Correctness is
    /// pinned by the unchanged naive-DFT oracle, Parseval, shift,
    /// round-trip, r2c/c2r, DCT, and N-D tests — all green across the
    /// change. HISTORY: radix-2 golden 0xbd55_68d2_33f4_b4bc (M4 + trj);
    /// radix-4/2 golden 0x0506_a4a0_955d_cf8e (M4; x86 row stayed
    /// pending). This value VERIFIED IN ALL FOUR QUADRANTS: aarch64
    /// (M4 Pro) and x86-64 (ts2) in debug and release each
    /// (2026-07-10, ts2:/data/tmp/fsim_xisa).
    const GOLDEN_HASH: u64 = 0x22dd_b617_266e_a792;

    #[test]
    fn non_power_of_two_is_refused_loudly() {
        for bad in [3usize, 6, 100] {
            let r = std::panic::catch_unwind(|| Fft::new(bad));
            assert!(r.is_err(), "size {bad} must be refused");
        }
    }

    // ---- Six-step path (bead 27d3): the dispatched battery lives in
    // tests/sixstep.rs behind `frontier-sixstep`; here we pin the
    // DEFAULT bit contract and the always-compiled transpose helper. ----

    #[test]
    fn sixstep_stays_off_by_default() {
        if cfg!(feature = "frontier-sixstep") {
            return; // the gated battery pins the enabled side
        }
        for n in [1usize << 16, 1 << 20, 1 << 22] {
            assert!(
                !Fft::takes_sixstep(n),
                "default build must keep n=2^{} on the stage walk",
                n.ilog2()
            );
            assert!(
                Fft::new(n).sub.is_none(),
                "default build must not carry a sub-plan"
            );
        }
    }

    #[test]
    fn transpose_into_matches_naive_and_round_trips() {
        let n1 = 24usize; // exercises the tile tails too
        let a0: Vec<C64> = (0..n1 * n1)
            .map(|i| C64::new(i as f64, -(i as f64) - 0.5))
            .collect();
        let mut t = vec![C64::default(); n1 * n1];
        transpose_into(&a0, &mut t, n1);
        for i in 0..n1 {
            for j in 0..n1 {
                assert_eq!(t[i * n1 + j], a0[j * n1 + i], "({i},{j})");
            }
        }
        let mut back = vec![C64::default(); n1 * n1];
        transpose_into(&t, &mut back, n1);
        assert_eq!(back, a0, "transpose twice must be the identity");
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
