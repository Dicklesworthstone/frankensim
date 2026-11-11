//! Joint propagation and integrated quadratic loss without an unstable inverse
//! propagator. For E(t) = exp(A t) and W(t) = integral E(s)^T C E(s) ds,
//! W(2t) = W(t) + E(t)^T W(t) E(t). Only the physical propagator is squared;
//! a rapidly decaying mode never requires forming exp(-A t).

mod affine;
pub(crate) use affine::flow_integrals as affine_flow_integrals;

use std::sync::atomic::{AtomicBool, Ordering};

pub(crate) fn identity(n: usize) -> Vec<f64> {
    let mut a = vec![0.0; n * n];
    for i in 0..n {
        a[i * n + i] = 1.0;
    }
    a
}

pub(crate) fn product(a: &[f64], b: &[f64], n: usize) -> Vec<f64> {
    let mut out = vec![0.0; n * n];
    accumulate(&mut out, a, b, n);
    out
}

/// Side length from which [`accumulate`] uses the declared workers: below it a product takes a few milliseconds and the
/// serial loop beats starting threads.
const PARALLEL_MIN_N: usize = 128;
/// Most threads one product uses: the bound the render's `--workers` option already holds (1 to 256).
const MAX_PRODUCT_WORKERS: usize = 256;

/// Threads one `n` x `n` product uses: one unless `n` reaches [`PARALLEL_MIN_N`] and more than one worker is declared,
/// and never more than [`MAX_PRODUCT_WORKERS`] or the number of four-row blocks.
fn product_threads(declared: usize, n: usize) -> usize {
    if n < PARALLEL_MIN_N { 1 } else { declared.clamp(1, MAX_PRODUCT_WORKERS).min(n.div_ceil(4)) }
}

// Preserve the existing i-k-j addition order, sharing each row of b between
// four output rows. Both the ordinary exponential and integral use this kernel.
// On the workers that `fs_la::eigen_complex::with_lu_workers` declares (one by default) the four-row blocks are shared
// between threads. A block reads its own rows of `a` and all of `b` and writes only its own rows of `out`, and every
// entry keeps its own ascending-k sum of unfused products, so the result is bit-identical for every worker count.
pub(crate) fn accumulate(out: &mut [f64], a: &[f64], b: &[f64], n: usize) {
    accumulate_polled(out, a, b, n, &|| false);
}

/// `out += a * b` for row-major `n` x `n` matrices by the kernel of [`accumulate`], polling `cancelled` once before each
/// four-row block of `out`, on the thread that runs the block, and running the blocks on the workers that
/// `fs_la::eigen_complex::with_lu_workers` declares (one by default, and one below side length 128). Returns `true` when every
/// block ran, and then `out` holds the same bits for every worker count. Returns `false` as soon as a poll answers `true`: `out` is
/// then partial and the caller must discard it. A request is seen by every thread within one block of work.
pub fn accumulate_polled(
    out: &mut [f64], a: &[f64], b: &[f64], n: usize, cancelled: &(impl Fn() -> bool + Sync),
) -> bool {
    if n == 0 {
        return true;
    }
    let stopped = AtomicBool::new(false);
    let threads = product_threads(fs_la::eigen_complex::declared_workers(), n);
    if threads == 1 {
        return accumulate_blocks(out, a, b, n, cancelled, &stopped);
    }
    // Whole four-row blocks per thread: each thread's blocks are exactly the serial blocks.
    let per_thread = n.div_ceil(4).div_ceil(threads) * 4 * n;
    std::thread::scope(|scope| {
        let stopped = &stopped;
        for (a_part, out_part) in a.chunks(per_thread).zip(out.chunks_mut(per_thread)) {
            scope.spawn(move || accumulate_blocks(out_part, a_part, b, n, cancelled, stopped));
        }
    });
    !stopped.load(Ordering::Relaxed)
}

fn accumulate_blocks(
    out: &mut [f64], a: &[f64], b: &[f64], n: usize, cancelled: &(impl Fn() -> bool + Sync), stopped: &AtomicBool,
) -> bool {
    for (a_rows, out_rows) in a.chunks(4 * n).zip(out.chunks_mut(4 * n)) {
        if stopped.load(Ordering::Relaxed) || cancelled() {
            stopped.store(true, Ordering::Relaxed);
            return false;
        }
        for (k, b_row) in b.chunks_exact(n).enumerate() {
            for (a_row, out_row) in a_rows.chunks_exact(n).zip(out_rows.chunks_exact_mut(n)) {
                let scale = a_row[k];
                for (value, &right) in out_row.iter_mut().zip(b_row) {
                    *value += scale * right;
                }
            }
        }
    }
    true
}

pub(crate) struct FlowIntegral {
    pub(crate) transition: Vec<f64>,
    pub(crate) loss_integral: Vec<f64>,
}

fn largest(a: &[f64]) -> f64 {
    a.iter().fold(0.0_f64, |norm, value| norm.max(value.abs()))
}

fn finite(a: &[f64]) -> bool {
    a.iter().all(|value| value.is_finite())
}

/// The caller admits the dense dimension/work budget before entry. The
/// original 32-term/48-doubling bounds remain. C is a quadratic loss form:
/// its symmetric part has the identical scalar value. No eigenvalue clipping,
/// inferred loss, state removal or change to the requested timestep is used.
/// Where the a priori rounding of the f64 doubling chain could exceed the
/// exact step's independent energy balance gate, the same Taylor stage and
/// chain run in double-double instead; everywhere else the f64 chain runs
/// unchanged, bit for bit.
pub(crate) fn flow_integral(
    generator: &[f64],
    loss: &[f64],
    n: usize,
    dt: f64,
) -> Result<FlowIntegral, &'static str> {
    if n == 0
        || n.checked_mul(n) != Some(generator.len())
        || loss.len() != generator.len()
        || !finite(generator)
        || !finite(loss)
        || !dt.is_finite()
        || dt < 0.0
    {
        return Err("invalid quadratic flow integral");
    }
    if dt == 0.0 {
        return Ok(FlowIntegral {
            transition: identity(n),
            loss_integral: vec![0.0; n * n],
        });
    }
    // Bound A and A^T so L(X)=A^T X+X A has norm <= 1/2 on the short
    // interval. C's magnitude does not force unnecessary squarings: its
    // series is linear in C and convergence is measured relatively.
    let mut norm = 0.0_f64;
    for i in 0..n {
        let row: f64 = (0..n).map(|j| (dt * generator[i * n + j]).abs()).sum();
        let column: f64 = (0..n).map(|j| (dt * generator[j * n + i]).abs()).sum();
        norm = norm.max(row).max(column);
    }
    if !norm.is_finite() {
        return Err("nonfinite quadratic flow operator");
    }
    let mut scale = 1.0;
    let mut doublings = 0_u32;
    while norm * scale > 0.25 {
        if doublings == 48 {
            return Err("quadratic flow scaling exceeds 48 doublings");
        }
        scale *= 0.5;
        doublings += 1;
    }
    let h = dt * scale;
    if doubling_rounding_bound(scale, n) > super::EXACT_BALANCE_RELATIVE {
        return compensated_chain(generator, loss, n, h, doublings);
    }
    plain_chain(generator, loss, n, h, doublings)
}

/// A priori relative rounding of the f64 doubling chain where the caller's
/// balance gate measures it. Each doubling forms E <- E E and
/// W <- W + sym(E^T W E) from n-term f64 dot products, each with relative error
/// at most gamma_n = n u / (1 - n u), u = 2^-53 (Higham, Accuracy and Stability
/// of Numerical Algorithms, 2nd ed., section 3.1). For the passive, nearly
/// norm-preserving blocks this protects (a nearly lossless string), the next
/// squaring doubles the error already carried, delta_k <= 2 delta_(k-1) +
/// gamma_n, so s doublings leave delta_s <= (2^s - 1) gamma_n; a quadratic
/// storage H = x^T Q x / 2 carries twice its state's relative error. `scale` is
/// the exact power 2^-s. On the stiff wound-core memory card
/// (frankensim-in05v) the measured balance defect grew as 2^s eps and sat
/// 67-232x under this bound. Double-double rounds at O(u^2) per operation, so
/// the compensated chain's own amplified rounding, about 2^s n u^2, stays far
/// below the gate for every admitted s <= 48 and n <= 12288, the hard
/// `MAX_PREPARED_EXPONENTIAL_DIM` (2^49 12288 u^2 is about 8.5e-14, 12 times under the 1e-12 gate).
fn doubling_rounding_bound(scale: f64, n: usize) -> f64 {
    let nu = n as f64 * (0.5 * f64::EPSILON);
    2.0 * (1.0 / scale - 1.0) * (nu / (1.0 - nu))
}

/// The original f64 Taylor stage and doubling chain, bit for bit.
fn plain_chain(
    generator: &[f64],
    loss: &[f64],
    n: usize,
    h: f64,
    doublings: u32,
) -> Result<FlowIntegral, &'static str> {
    let a: Vec<_> = generator.iter().map(|value| h * value).collect();
    let mut transition = identity(n);
    let mut e_term = transition.clone();
    let mut loss_integral = vec![0.0; n * n];
    for i in 0..n {
        for j in i..n {
            let value = h * f64::midpoint(loss[i * n + j], loss[j * n + i]);
            loss_integral[i * n + j] = value;
            loss_integral[j * n + i] = value;
        }
    }
    let mut w_term = loss_integral.clone();
    let mut converged = false;
    for order in 1..=32 {
        e_term = product(&e_term, &a, n);
        for (sum, term) in transition.iter_mut().zip(&mut e_term) {
            *term /= f64::from(order);
            *sum += *term;
        }
        // W = h C + L(h C)/2! + L^2(h C)/3! + ... . All terms are
        // symmetric, so one product gives both halves of L(term).
        let wa = product(&w_term, &a, n);
        for i in 0..n {
            for j in i..n {
                let value = (wa[i * n + j] + wa[j * n + i]) / f64::from(order + 1);
                w_term[i * n + j] = value;
                w_term[j * n + i] = value;
                loss_integral[i * n + j] += value;
                if i != j {
                    loss_integral[j * n + i] += value;
                }
            }
        }
        if !finite(&transition) || !finite(&loss_integral) {
            return Err("quadratic flow series is unrepresentable");
        }
        let tolerance = 0.25 * f64::EPSILON;
        if largest(&e_term) <= tolerance * largest(&transition)
            && largest(&w_term) <= tolerance * largest(&loss_integral)
        {
            converged = true;
            break;
        }
    }
    if !converged {
        return Err("quadratic flow series did not converge in 32 terms");
    }
    for _ in 0..doublings {
        let we = product(&loss_integral, &transition, n);
        for i in 0..n {
            for j in i..n {
                let mut upper = 0.0;
                let mut lower = 0.0;
                for k in 0..n {
                    upper += transition[k * n + i] * we[k * n + j];
                    lower += transition[k * n + j] * we[k * n + i];
                }
                // Averaging preserves the quadratic form, not a PSD repair.
                loss_integral[i * n + j] += f64::midpoint(upper, lower);
                loss_integral[j * n + i] = loss_integral[i * n + j];
            }
        }
        transition = product(&transition, &transition, n);
        if !finite(&transition) || !finite(&loss_integral) {
            return Err("quadratic flow doubling is unrepresentable");
        }
    }
    Ok(FlowIntegral {
        transition,
        loss_integral,
    })
}

// Double-double building blocks without fused multiply-add: Knuth's two-sum,
// Dekker's split (2^27 + 1) and two-product. Rust never contracts a * b + c.
fn two_sum(a: f64, b: f64) -> (f64, f64) {
    let s = a + b;
    let bb = s - a;
    (s, (a - (s - bb)) + (b - bb))
}

fn quick_two_sum(a: f64, b: f64) -> (f64, f64) {
    let s = a + b;
    (s, b - (s - a))
}

fn split(a: f64) -> (f64, f64) {
    let t = 134_217_729.0 * a;
    let hi = t - (t - a);
    (hi, a - hi)
}

fn two_prod(a: f64, b: f64) -> (f64, f64) {
    let p = a * b;
    let (ah, al) = split(a);
    let (bh, bl) = split(b);
    (p, ((ah * bh - p) + ah * bl + al * bh) + al * bl)
}

fn dd_add((ah, al): (f64, f64), (bh, bl): (f64, f64)) -> (f64, f64) {
    let (s, e) = two_sum(ah, bh);
    let (t, f) = two_sum(al, bl);
    let (s, e) = quick_two_sum(s, e + t);
    quick_two_sum(s, e + f)
}

fn dd_mul((ah, al): (f64, f64), (bh, bl): (f64, f64)) -> (f64, f64) {
    let (p, e) = two_prod(ah, bh);
    quick_two_sum(p, e + (ah * bl + al * bh))
}

fn dd_mul_plain((ah, al): (f64, f64), b: f64) -> (f64, f64) {
    let (p, e) = two_prod(ah, b);
    quick_two_sum(p, e + al * b)
}

fn dd_div((ah, al): (f64, f64), b: f64) -> (f64, f64) {
    let q1 = ah / b;
    let (p, e) = two_prod(q1, b);
    let (s, f) = two_sum(ah, -p);
    quick_two_sum(q1, (s + (f - e + al)) / b)
}

/// An n x n matrix of double-double entries, row major.
#[derive(Clone)]
struct Compensated(Vec<(f64, f64)>);

impl Compensated {
    fn zeros(n: usize) -> Self {
        Self(vec![(0.0, 0.0); n * n])
    }

    fn identity(n: usize) -> Self {
        let mut a = Self::zeros(n);
        for i in 0..n {
            a.0[i * n + i] = (1.0, 0.0);
        }
        a
    }

    fn product(&self, b: &Self, n: usize) -> Self {
        let mut out = Self::zeros(n);
        for i in 0..n {
            for k in 0..n {
                let left = self.0[i * n + k];
                for j in 0..n {
                    out.0[i * n + j] = dd_add(out.0[i * n + j], dd_mul(left, b.0[k * n + j]));
                }
            }
        }
        out
    }

    /// The product with an f64 matrix (the scaled generator, exact as given).
    fn times_plain(&self, b: &[f64], n: usize) -> Self {
        let mut out = Self::zeros(n);
        for i in 0..n {
            for k in 0..n {
                let left = self.0[i * n + k];
                for j in 0..n {
                    out.0[i * n + j] = dd_add(out.0[i * n + j], dd_mul_plain(left, b[k * n + j]));
                }
            }
        }
        out
    }

    fn largest(&self) -> f64 {
        self.0
            .iter()
            .fold(0.0_f64, |norm, value| norm.max(value.0.abs()))
    }

    fn finite(&self) -> bool {
        self.0.iter().all(|value| value.0.is_finite())
    }

    fn high(self) -> Vec<f64> {
        self.0.into_iter().map(|value| value.0).collect()
    }
}

/// `plain_chain` with the Taylor stage and every doubling carried in
/// double-double and rounded to f64 once at the end; the same admissions,
/// term limit and refusals. Its series stops at 0.25 eps^2 so the truncation
/// is not amplified by 2^s either: with ||A h|| <= 1/4 the k-th term is at most
/// 4^-k / k!, under that tolerance by k = 21 of the 32 admitted terms.
fn compensated_chain(
    generator: &[f64],
    loss: &[f64],
    n: usize,
    h: f64,
    doublings: u32,
) -> Result<FlowIntegral, &'static str> {
    let a: Vec<_> = generator.iter().map(|value| h * value).collect();
    let mut transition = Compensated::identity(n);
    let mut e_term = transition.clone();
    let mut loss_integral = Compensated::zeros(n);
    for i in 0..n {
        for j in i..n {
            let value = two_prod(h, f64::midpoint(loss[i * n + j], loss[j * n + i]));
            loss_integral.0[i * n + j] = value;
            loss_integral.0[j * n + i] = value;
        }
    }
    let mut w_term = loss_integral.clone();
    let tolerance = 0.25 * f64::EPSILON * f64::EPSILON;
    let mut converged = false;
    for order in 1..=32 {
        e_term = e_term.times_plain(&a, n);
        for (sum, term) in transition.0.iter_mut().zip(&mut e_term.0) {
            *term = dd_div(*term, f64::from(order));
            *sum = dd_add(*sum, *term);
        }
        let wa = w_term.times_plain(&a, n);
        for i in 0..n {
            for j in i..n {
                let value = dd_div(
                    dd_add(wa.0[i * n + j], wa.0[j * n + i]),
                    f64::from(order + 1),
                );
                w_term.0[i * n + j] = value;
                w_term.0[j * n + i] = value;
                let sum = dd_add(loss_integral.0[i * n + j], value);
                loss_integral.0[i * n + j] = sum;
                loss_integral.0[j * n + i] = sum;
            }
        }
        if !transition.finite() || !loss_integral.finite() {
            return Err("quadratic flow series is unrepresentable");
        }
        if e_term.largest() <= tolerance * transition.largest()
            && w_term.largest() <= tolerance * loss_integral.largest()
        {
            converged = true;
            break;
        }
    }
    if !converged {
        return Err("quadratic flow series did not converge in 32 terms");
    }
    for _ in 0..doublings {
        let we = loss_integral.product(&transition, n);
        for i in 0..n {
            for j in i..n {
                let mut upper = (0.0, 0.0);
                let mut lower = (0.0, 0.0);
                for k in 0..n {
                    upper = dd_add(upper, dd_mul(transition.0[k * n + i], we.0[k * n + j]));
                    lower = dd_add(lower, dd_mul(transition.0[k * n + j], we.0[k * n + i]));
                }
                // The midpoint, as in the f64 chain; halving is exact.
                let (sh, sl) = dd_add(upper, lower);
                let sum = dd_add(loss_integral.0[i * n + j], (0.5 * sh, 0.5 * sl));
                loss_integral.0[i * n + j] = sum;
                loss_integral.0[j * n + i] = sum;
            }
        }
        transition = transition.product(&transition, n);
        if !transition.finite() || !loss_integral.finite() {
            return Err("quadratic flow doubling is unrepresentable");
        }
    }
    Ok(FlowIntegral {
        transition: transition.high(),
        loss_integral: loss_integral.high(),
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::linear_exact::{
        DenseExponentialBudget, UpperBlocks, exponential_with_budget, upper_block_exponential_with_budget,
    };
    use fs_la::eigen_complex::with_lu_workers;

    /// Deterministic entries in [-1, 1) with every bit of the significand in play (xorshift64).
    fn filled(len: usize, seed: u64) -> Vec<f64> {
        let mut state = seed.wrapping_mul(0x9E37_79B9_7F4A_7C15) | 1;
        (0..len)
            .map(|_| {
                state ^= state << 13;
                state ^= state >> 7;
                state ^= state << 17;
                (state >> 11) as f64 / (1_u64 << 53) as f64 * 2.0 - 1.0
            })
            .collect()
    }

    fn bits(values: &[f64]) -> Vec<u64> {
        values.iter().map(|v| v.to_bits()).collect()
    }

    /// Gauntlet G0/G5: the thread count follows the declaration (one when undeclared), the side-length floor and the number of
    /// four-row blocks, and never exceeds the bound of the render's `--workers` option.
    #[test]
    fn product_threads_follow_the_declaration_the_floor_and_the_block_count() {
        assert_eq!(fs_la::eigen_complex::declared_workers(), 1);
        assert_eq!(with_lu_workers(3, fs_la::eigen_complex::declared_workers), 3);
        assert_eq!(fs_la::eigen_complex::declared_workers(), 1);
        assert_eq!(product_threads(1, 5433), 1);
        assert_eq!(product_threads(0, 5433), 1);
        assert_eq!(product_threads(8, PARALLEL_MIN_N - 1), 1);
        assert_eq!(product_threads(8, PARALLEL_MIN_N), 8);
        assert_eq!(product_threads(1000, 129), 33);
        assert_eq!(product_threads(1000, 5433), MAX_PRODUCT_WORKERS);
    }

    /// Gauntlet G5: a product on any declared worker count equals the serial one bit for bit, on sizes that are and are not
    /// multiples of four (the last block is short) and from a nonzero start (the kernel adds into `out`).
    #[test]
    fn threaded_products_equal_the_serial_product_bit_for_bit() {
        for n in [PARALLEL_MIN_N, PARALLEL_MIN_N + 1, PARALLEL_MIN_N + 2, PARALLEL_MIN_N + 3, 150, 203] {
            let (a, b, start) = (filled(n * n, 1), filled(n * n, 2), filled(n * n, 3));
            let mut serial = start.clone();
            accumulate(&mut serial, &a, &b, n);
            assert_ne!(bits(&serial), bits(&start));
            for workers in [2, 3, 5, 8] {
                assert!(product_threads(workers, n) > 1, "n {n}, {workers} workers would run serially");
                let mut threaded = start.clone();
                with_lu_workers(workers, || accumulate(&mut threaded, &a, &b, n));
                assert_eq!(bits(&threaded), bits(&serial), "n {n}, {workers} workers");
            }
        }
    }

    /// Gauntlet G5: the budgeted exponential and the block upper-triangular one, whose products are these, return the serial
    /// bits on every declared worker count. The scale puts the row sums near one, so a squaring follows the Taylor series.
    #[test]
    fn threaded_exponentials_equal_the_serial_ones_through_both_budgeted_entry_points() {
        let d = PARALLEL_MIN_N + 1;
        let scaled = |seed: u64| -> Vec<f64> { filled(d * d, seed).into_iter().map(|v| v * 0.012).collect() };
        let budget = |n: usize| DenseExponentialBudget {
            max_dimension: n,
            max_dense_bytes: 5 * n * n * core::mem::size_of::<f64>(),
            max_multiply_terms: 80 * (n as u64).pow(3),
        };
        let full = scaled(4);
        let blocks = UpperBlocks { first: scaled(5), cross: scaled(6), last: scaled(7) };
        let run = || {
            (
                exponential_with_budget(&full, d, budget(d)).unwrap(),
                upper_block_exponential_with_budget(&blocks, d, budget(2 * d)).unwrap(),
            )
        };
        let (serial_full, serial_blocks) = run();
        for workers in [2, 3, 8] {
            let (threaded_full, threaded_blocks) = with_lu_workers(workers, run);
            assert_eq!(bits(&threaded_full), bits(&serial_full), "exponential, {workers} workers");
            assert_eq!(bits(&threaded_blocks.first), bits(&serial_blocks.first), "first block, {workers} workers");
            assert_eq!(bits(&threaded_blocks.cross), bits(&serial_blocks.cross), "cross block, {workers} workers");
            assert_eq!(bits(&threaded_blocks.last), bits(&serial_blocks.last), "last block, {workers} workers");
        }
    }

    /// Gauntlet G5: the polled kernel polls once before each four-row block (counted, serial and threaded), returns the plain
    /// kernel's bits when no poll answers, and stops within one poll per thread when one does.
    #[test]
    fn polled_accumulate_polls_once_per_block_stops_on_a_request_and_equals_the_plain_kernel() {
        use core::sync::atomic::{AtomicUsize, Ordering};
        assert!(accumulate_polled(&mut [], &[], &[], 0, &|| true));
        let n = PARALLEL_MIN_N + 3;
        let blocks = n.div_ceil(4);
        let (a, b) = (filled(n * n, 21), filled(n * n, 22));
        let mut plain = vec![0.0; n * n];
        accumulate(&mut plain, &a, &b, n);
        for workers in [1, 3, 8] {
            let polls = AtomicUsize::new(0);
            let mut out = vec![0.0; n * n];
            let done = with_lu_workers(workers, || {
                accumulate_polled(&mut out, &a, &b, n, &|| {
                    polls.fetch_add(1, Ordering::Relaxed);
                    false
                })
            });
            assert!(done, "{workers} workers");
            assert_eq!(polls.load(Ordering::Relaxed), blocks, "{workers} workers poll once per block");
            assert_eq!(bits(&out), bits(&plain), "{workers} workers");
            // The first poll answers no and every later one yes: the call returns false after at most one poll per thread
            // and one more.
            let polls = AtomicUsize::new(0);
            let mut out = vec![0.0; n * n];
            let done = with_lu_workers(workers, || {
                accumulate_polled(&mut out, &a, &b, n, &|| polls.fetch_add(1, Ordering::Relaxed) >= 1)
            });
            assert!(!done, "{workers} workers");
            let seen = polls.load(Ordering::Relaxed);
            assert!((2..=workers + 1).contains(&seen), "{workers} workers polled {seen} times");
        }
    }

    /// The D4c film page's shape (d = 5433): one product at one worker and at 8 and 16, equal bit for bit, with each wall time
    /// printed. Run by hand: `cargo test --release -p fs-phs --lib -- --ignored --nocapture --test-threads 1 d4c_shape`.
    #[test]
    #[ignore = "measurement: 1.6e11 multiply-adds per run, a few minutes serial"]
    fn d4c_shape_product_is_bitwise_equal_threaded_and_reports_the_wall_times() {
        let n = 5433;
        let (a, b) = (filled(n * n, 11), filled(n * n, 12));
        let mut serial = vec![0.0; n * n];
        let start = std::time::Instant::now();
        accumulate(&mut serial, &a, &b, n);
        let serial_s = start.elapsed().as_secs_f64();
        let serial_bits = bits(&serial);
        eprintln!("d4c shape n {n}: 1 worker {serial_s:.2} s");
        for workers in [8, 16] {
            let mut threaded = vec![0.0; n * n];
            let start = std::time::Instant::now();
            with_lu_workers(workers, || accumulate(&mut threaded, &a, &b, n));
            let threaded_s = start.elapsed().as_secs_f64();
            eprintln!("d4c shape n {n}: {workers} workers {threaded_s:.2} s, speedup {:.2}", serial_s / threaded_s);
            assert_eq!(bits(&threaded), serial_bits, "{workers} workers");
        }
    }

    fn close(actual: f64, expected: f64, scale: f64) {
        assert!(
            (actual - expected).abs() <= 8.0e-13 * scale,
            "actual={actual:.17e}, expected={expected:.17e}"
        );
    }

    #[test]
    fn passive_decay_resolves_heat_after_a_million_relaxation_times() {
        for rate in [1.0e-8_f64, 0.25, 10.0, 1.0e3, 1.0e6] {
            let flow = flow_integral(&[-rate], &[rate], 1, 1.0).unwrap();
            let decay = (-rate).exp();
            let heat = -0.5 * (-2.0 * rate).exp_m1();
            close(flow.transition[0], decay, decay.max(f64::MIN_POSITIVE));
            close(flow.loss_integral[0], heat, heat);
            close(
                0.5 * flow.transition[0].powi(2) + flow.loss_integral[0],
                0.5,
                0.5,
            );
        }
    }

    #[test]
    fn stiff_driven_storage_closes_independent_source_work_and_heat() {
        // x'=-a*x+u, H=x^2/2. Hold the source and integrate the conjugate
        // output x as additional states: the real exact-ZOH augmentation.
        let rate = 10_000.0;
        let input = 2_500.0;
        let flow = flow_integral(
            &[-rate, input, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0],
            &[rate, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            3,
            1.0,
        )
        .unwrap();
        let initial = [1.0, 1.0, 0.0];
        let end: Vec<f64> = flow
            .transition
            .chunks_exact(3)
            .map(|row| row.iter().zip(initial).map(|(a, b)| a * b).sum())
            .collect();
        let heat: f64 = (0..3)
            .map(|i| {
                (0..3)
                    .map(|j| initial[i] * flow.loss_integral[3 * i + j] * initial[j])
                    .sum::<f64>()
            })
            .sum();
        close(end[0], 0.25, 1.0);
        close(end[2], 0.25 + 0.75 / rate, 1.0);
        close(heat, 625.65625, 625.65625);
        let source_work = input * end[2];
        close(
            0.5 * (end[0].powi(2) - 1.0) + heat,
            source_work,
            source_work,
        );
    }

    #[test]
    fn nilpotent_transport_retains_zero_modes_and_exact_polynomial_work() {
        for t in [0.0, 0.001, 0.5, 25.0] {
            let flow = flow_integral(&[0.0, 1.0, 0.0, 0.0], &[1.0, 0.0, 0.0, 1.0], 2, t).unwrap();
            for (actual, expected) in flow.transition.iter().zip([1.0, t, 0.0, 1.0]) {
                close(*actual, expected, expected.abs().max(1.0));
            }
            for (actual, expected) in
                flow.loss_integral
                    .iter()
                    .zip([t, t * t / 2.0, t * t / 2.0, t + t * t * t / 3.0])
            {
                close(*actual, expected, expected.abs().max(1.0e-30));
            }
        }
    }

    #[test]
    fn noncommuting_damping_and_rotation_match_the_analytic_energy_integral() {
        let (r1, r2, omega, t) = (2.0_f64, 7.0_f64, 3.0_f64, 0.8_f64);
        let generator = [-r1, omega, -omega, -r2];
        let flow = flow_integral(&generator, &[r1, 0.0, 0.0, r2], 2, t).unwrap();
        let first = flow_integral(&generator, &[r1, 0.0, 0.0, 0.0], 2, t).unwrap();
        let second = flow_integral(&generator, &[0.0, 0.0, 0.0, r2], 2, t).unwrap();
        let mu = 0.5 * (r1 + r2);
        let delta = 0.5 * (r2 - r1);
        let nu = (omega * omega - delta * delta).sqrt();
        let decay = (-mu * t).exp();
        let sine = (nu * t).sin() / nu;
        let cosine = (nu * t).cos();
        let e = [
            decay * (cosine + delta * sine),
            decay * omega * sine,
            -decay * omega * sine,
            decay * (cosine - delta * sine),
        ];
        for i in 0..2 {
            for j in 0..2 {
                close(flow.transition[i * 2 + j], e[i * 2 + j], 1.0);
                // H=x^T x/2 implies W=(I-E^T E)/2; independent of the
                // production Taylor recurrence and its semigroup doubling.
                let identity = if i == j { 1.0 } else { 0.0 };
                let expected = 0.5 * (identity - e[i] * e[j] - e[2 + i] * e[2 + j]);
                close(flow.loss_integral[i * 2 + j], expected, 1.0);
                close(
                    first.loss_integral[i * 2 + j] + second.loss_integral[i * 2 + j],
                    expected,
                    1.0,
                );
            }
        }
        assert_eq!(flow.loss_integral[1], flow.loss_integral[2]);
    }

    #[test]
    fn a_long_doubling_chain_holds_the_analytic_flow_to_final_rounding() {
        // A weakly damped fast rotation, x' = A x with H = x^T x / 2, needs
        // 23 doublings; there the f64 chain's bound (about 4e-9) is over the
        // gate.
        let (rate, omega) = (1.0_f64, 1_048_576.0_f64);
        let generator = [-rate, omega, -omega, -rate];
        assert!(
            doubling_rounding_bound(0.5_f64.powi(23), 2) > super::super::EXACT_BALANCE_RELATIVE
        );
        let flow = flow_integral(&generator, &[rate, 0.0, 0.0, rate], 2, 1.0).unwrap();
        let decay = (-rate).exp();
        let (sine, cosine) = omega.sin_cos();
        let e = [decay * cosine, decay * sine, -decay * sine, decay * cosine];
        // E^T E = exp(-2 rate) I, so W = (I - E^T E) / 2 is diagonal.
        let heat = -0.5 * (-2.0 * rate).exp_m1();
        let w = [heat, 0.0, 0.0, heat];
        // Allowance on entries of magnitude <= 1: the result's one rounding to
        // f64 (u), the reference's exp, sin_cos and product (1 ulp, 1 ulp and
        // half an ulp: under 5u relative on e^-1, 1.9u absolute), and the
        // double-double chain's own amplified rounding (2^23 * 2 * u^2, under
        // 1e-8 u). The sum stays under 4u = 2 eps.
        let allowance = 2.0 * f64::EPSILON;
        for k in 0..4 {
            assert!(
                (flow.transition[k] - e[k]).abs() <= allowance,
                "E[{k}] = {:.17e}, analytic {:.17e}",
                flow.transition[k],
                e[k]
            );
            assert!(
                (flow.loss_integral[k] - w[k]).abs() <= allowance,
                "W[{k}] = {:.17e}, analytic {:.17e}",
                flow.loss_integral[k],
                w[k]
            );
        }
    }

    #[test]
    fn under_the_rounding_bound_the_f64_chain_runs_bit_for_bit() {
        // Five exact doublings: norm 0.5 * 10 = 5 needs 2^-5 to reach 1/4. The
        // bound, 62 gamma_2 (about 1.4e-14), is under the gate.
        let generator = [-2.0, 3.0, -3.0, -7.0];
        let loss = [2.0, 0.0, 0.0, 7.0];
        let (t, scale, doublings) = (0.5, 0.5_f64.powi(5), 5);
        assert!(doubling_rounding_bound(scale, 2) <= super::super::EXACT_BALANCE_RELATIVE);
        let flow = flow_integral(&generator, &loss, 2, t).unwrap();
        let plain = plain_chain(&generator, &loss, 2, t * scale, doublings).unwrap();
        let compensated = compensated_chain(&generator, &loss, 2, t * scale, doublings).unwrap();
        let bits = |flow: &FlowIntegral| {
            flow.transition
                .iter()
                .chain(&flow.loss_integral)
                .map(|value| value.to_bits())
                .collect::<Vec<_>>()
        };
        assert_eq!(bits(&flow), bits(&plain));
        // Not vacuous: the compensated chain changes bits on this case.
        assert_ne!(bits(&compensated), bits(&plain));
    }

    #[test]
    fn malformed_or_unrepresentable_integrals_refuse() {
        assert!(flow_integral(&[], &[], 0, 1.0).is_err());
        assert!(flow_integral(&[0.0], &[0.0; 2], 1, 1.0).is_err());
        assert!(flow_integral(&[f64::NAN], &[1.0], 1, 1.0).is_err());
        assert!(flow_integral(&[-1.0], &[f64::INFINITY], 1, 1.0).is_err());
        assert!(flow_integral(&[-1.0], &[1.0], 1, -1.0).is_err());
        assert!(flow_integral(&[-1.0e18], &[1.0], 1, 1.0).is_err());
    }
}
