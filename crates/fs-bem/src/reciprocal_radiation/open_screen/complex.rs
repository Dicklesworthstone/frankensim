//! Holomorphic continuation of the SAME real-frequency open-screen operator.
//!
//! The real owner uses the cosine Green remainder and an angular outgoing Gram
//! form, not the point-quadrature imaginary remainder of solve_laplace. Write
//! c_+(k,d)=-ik integral (n.d)phi exp(-ik d.x) and c_-(k,d)=c_+(k,-d).
//! On the real axis Re(conj(c_i)c_j) equals
//! (c_-,i*c_+,j+c_+,i*c_-,j)/2. Continuing THAT symmetric product, without
//! conjugating complex k, gives an entire finite operator and its derivative.
//!
//! Real-axis matrix arithmetic remains the original owner's operation order.
//! Complex Gram factors use one fixed centered geometry reference; the common
//! phase cancels analytically in each antipodal product. This bounded numerical
//! pole chart is not a positive-real/causal certificate in the whole Laplace
//! half-plane. It does not change solve_laplace or provide a CQ reaction model.
use super::*;

mod kernel;
use kernel::{cosine_remainder, plane_factors};

/// Complete two-matrix and temporary payload/work admission, separate from the
/// retained original geometry cache. No factorization or source elimination.
#[derive(Debug, Clone, Copy)]
pub struct ScreenComplexLimits {
    /// Total complex entries across BOTH H and dH/dk, namely 2*n*n.
    pub maximum_matrix_entries: usize,
    /// Both output matrices, bounded parallel pair buffers, angular-rule and
    /// coefficient arrays. Excludes cached geometry, allocator/thread-stack
    /// overhead and any downstream mechanical pencil or factorization.
    pub maximum_workspace_bytes: usize,
    /// Eight times the original assembly work (its n^3 LU term is excluded).
    /// This conservative derivative admission is not a machine FLOP count.
    /// The original OpenOptions::maximum_work must also admit this work.
    pub maximum_work: u64,
    /// Positive and <=8. Require |Im k|*2R <= this value, where R covers all
    /// vertices and cached quadrature points about the fixed bounding-box center.
    /// Individual centered exponential magnitudes are then at most exp(limit/2).
    /// This bounds numerical growth, not spatial or resonance error.
    pub maximum_growth_exponent: f64,
}

/// The unchanged pressure-test ordering, with analytic derivative in k [1/m].
/// H has units m and dH/dk has units m². For an angular-frequency pencil use
/// k=omega/c and dH/domega=(dH/dk)/c, retaining the pressure unknowns.
#[derive(Debug, Clone)]
pub struct ScreenComplexOperator {
    pub dimension: usize,
    pub matrix: Vec<C64>,
    pub derivative_wavenumber: Vec<C64>,
}

struct Admission {
    origin: Point,
    workers: usize,
    block: usize,
}

#[derive(Clone, Copy)]
struct PairDerivative {
    single: C64,
    shape: [[C64; 3]; 3],
    derivative_single: C64,
    derivative_shape: [[C64; 3]; 3],
}

impl PreparedScreen<'_> {
    /// Admit a COMPLETE complex-frequency domain before sampling its boundary.
    /// The supplied upper bounds may be conservative; both are nonnegative.
    /// Zero belongs to this entire open-screen chart, unlike Burton--Miller's
    /// closed-body i/k coupling. Resources are checked before any allocation.
    pub fn admit_complex_domain(
        &self,
        maximum_abs_k: f64,
        maximum_abs_imag_k: f64,
        limits: ScreenComplexLimits,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<(), HelmholtzError> {
        self.complex_admission(maximum_abs_k, maximum_abs_imag_k, limits, &mut cancelled)?;
        Ok(())
    }

    fn complex_admission(
        &self,
        maximum_abs_k: f64,
        maximum_abs_imag_k: f64,
        limits: ScreenComplexLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Admission, HelmholtzError> {
        poll(cancelled)?;
        if !maximum_abs_k.is_finite()
            || maximum_abs_k < 0.0
            || !maximum_abs_imag_k.is_finite()
            || maximum_abs_imag_k < 0.0
            || maximum_abs_k > self.highest_k
            || !(maximum_abs_k * maximum_abs_k).is_finite()
            || !limits.maximum_growth_exponent.is_finite()
            || limits.maximum_growth_exponent <= 0.0
            || limits.maximum_growth_exponent > 8.0
        {
            return Err(bad("invalid or out-of-band complex open-screen domain"));
        }
        let n = self.n as u128;
        let nf = self.facets.len();
        let workers = self.options.workers.min(nf);
        let block = (ROWS_PER_WORKER * workers).min(nf);
        let entries = 2 * n * n;
        let directions = 64_u128 << (2 * self.options.angular_refinement);
        let pair_bytes = if workers > 1 {
            block as u128 * nf as u128 * core::mem::size_of::<PairDerivative>() as u128
                + 128 * (block + workers) as u128
        } else {
            0
        };
        let bytes = entries * core::mem::size_of::<C64>() as u128
            + pair_bytes
            + 5 * n * core::mem::size_of::<C64>() as u128
            + directions * core::mem::size_of::<(Point, f64)>() as u128
            + 4096;
        let work = 8 * (self.frequency_work() - n * n * n);
        if entries > limits.maximum_matrix_entries as u128
            || bytes > limits.maximum_workspace_bytes as u128
            || work > u128::from(limits.maximum_work)
            || work > u128::from(self.options.maximum_work)
        {
            return Err(bad(
                "complex open-screen matrices or assembly exceed declared resources",
            ));
        }
        let mut lower = self.surface.vertices[0];
        let mut upper = lower;
        for (i, point) in self.surface.vertices.iter().enumerate() {
            if i % 256 == 0 {
                poll(cancelled)?;
            }
            for a in 0..3 {
                lower[a] = lower[a].min(point[a]);
                upper[a] = upper[a].max(point[a]);
            }
        }
        let origin = core::array::from_fn(|a| 0.5 * lower[a] + 0.5 * upper[a]);
        let mut radius = 0.0_f64;
        for (i, &point) in self.surface.vertices.iter().enumerate() {
            if i % 256 == 0 {
                poll(cancelled)?;
            }
            radius = radius.max(norm(sub(point, origin)));
        }
        // Cached quadrature points also participate, including their actual
        // coordinate roundoff rather than only the ideal convex triangle.
        for row in &self.samples {
            poll(cancelled)?;
            for sample in row {
                radius = radius.max(norm(sub(sample.point, origin)));
            }
        }
        let span = 2.0 * radius * (1.0 + 64.0 * f64::EPSILON);
        let growth = maximum_abs_imag_k * span;
        if !span.is_finite()
            || span <= 0.0
            || !growth.is_finite()
            || growth > limits.maximum_growth_exponent
        {
            return Err(bad(
                "complex open-screen domain exceeds its centered exponential-growth chart",
            ));
        }
        poll(cancelled)?;
        Ok(Admission {
            origin,
            workers,
            block,
        })
    }

    /// Assemble H and its analytic dH/dk on a fixed geometric/angular rule.
    /// Both signs of Im(k), either real-frequency sign, and zero are admitted
    /// inside the declared band/growth chart. Real positive k reproduces the
    /// original rigid-screen matrix bit-for-bit. The derivative and complex
    /// values use the analytic cosine/antipodal formula, not conjugated k or
    /// a finite-difference quotient. A failure returns no partial operator.
    pub fn complex_operator(
        &self,
        k: C64,
        limits: ScreenComplexLimits,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<ScreenComplexOperator, HelmholtzError> {
        poll(&mut cancelled)?;
        if !finite(k) {
            return Err(bad("complex open-screen wavenumber is nonfinite"));
        }
        let chart = self.complex_admission(k.abs(), k.im.abs(), limits, &mut cancelled)?;
        let real_axis = k.im == 0.0;
        let k = if real_axis { C64::from_re(k.re) } else { k };
        let k2 = k * k;
        let n = self.n;
        let mut matrix = zeroed(n * n, C64::ZERO)?;
        let mut derivative = zeroed(n * n, C64::ZERO)?;
        for i in 0..n * n {
            if i % n == 0 {
                poll(&mut cancelled)?;
            }
            matrix[i] = C64::from_re(self.static_curl[i]) + k2.scale(self.static_normal[i]);
            derivative[i] = k.scale(2.0 * self.static_normal[i]);
        }
        let mut add = |i, j, pair: PairDerivative| {
            self.add_dynamic(&mut matrix, k2, i, j, pair.single, pair.shape);
            self.add_dynamic(
                &mut derivative,
                k2,
                i,
                j,
                pair.derivative_single,
                pair.derivative_shape,
            );
            self.add_dynamic(&mut derivative, k.scale(2.0), i, j, C64::ZERO, pair.shape);
        };
        let nf = self.facets.len();
        if chart.workers <= 1 {
            for i in 0..nf {
                for j in i..nf {
                    poll(&mut cancelled)?;
                    add(i, j, self.complex_pair(k, real_axis, i, j));
                }
            }
        } else {
            for first in (0..nf).step_by(chart.block) {
                poll(&mut cancelled)?;
                let rows = first..nf.min(first + chart.block);
                let pairs = rows_on_workers(rows.clone(), chart.workers, |i| {
                    (i..nf)
                        .map(|j| self.complex_pair(k, real_axis, i, j))
                        .collect::<Vec<_>>()
                });
                for (i, row) in rows.zip(pairs) {
                    for (j, pair) in (i..nf).zip(row) {
                        add(i, j, pair);
                    }
                }
            }
        }
        // The closure's mutable matrix borrows end before the angular part.
        drop(add);
        for (direction, weight) in directions(self.options.angular_refinement) {
            poll(&mut cancelled)?;
            let coefficients =
                self.complex_coefficients(k, direction, chart.origin, &mut cancelled)?;
            let original = if real_axis {
                self.coefficients(k.re, direction)
            } else {
                Vec::new()
            };
            let [forward, backward, d_forward, d_backward] = coefficients;
            let factor = C64::new(0.0, weight / (8.0 * PI * PI));
            for i in 0..n {
                poll(&mut cancelled)?;
                for j in 0..n {
                    let product = (backward[i] * forward[j] + forward[i] * backward[j]).scale(0.5);
                    let d_product = (d_backward[i] * forward[j]
                        + backward[i] * d_forward[j]
                        + d_forward[i] * backward[j]
                        + forward[i] * d_backward[j])
                        .scale(0.5);
                    let at = i * n + j;
                    if real_axis {
                        // The exact original outgoing Gram arithmetic and sum
                        // order, including its absolute-coordinate phase.
                        let scale = k.re * weight / (8.0 * PI * PI);
                        matrix[at].im += scale * (original[i].conj() * original[j]).re;
                    } else {
                        matrix[at] = matrix[at] + (factor * k) * product;
                    }
                    derivative[at] = derivative[at] + factor * (product + k * d_product);
                }
            }
        }
        if matrix.iter().chain(&derivative).any(|&v| !finite(v)) {
            return Err(bad("complex open-screen operator or derivative overflow"));
        }
        poll(&mut cancelled)?;
        Ok(ScreenComplexOperator {
            dimension: n,
            matrix,
            derivative_wavenumber: derivative,
        })
    }

    fn complex_pair(&self, k: C64, real_axis: bool, i: usize, j: usize) -> PairDerivative {
        let mut out = PairDerivative {
            single: C64::ZERO,
            shape: [[C64::ZERO; 3]; 3],
            derivative_single: C64::ZERO,
            derivative_shape: [[C64::ZERO; 3]; 3],
        };
        for x in &self.samples[i] {
            for y in &self.samples[j] {
                let r = norm(sub(x.point, y.point));
                let (g, dg) = if real_axis {
                    let mut g = green_remainder(k.re, r);
                    g.im = 0.0;
                    (g, C64::from_re(-det::sin(k.re * r) / (4.0 * PI)))
                } else {
                    cosine_remainder(k, r)
                };
                let g = g.scale(x.weight * y.weight);
                let dg = dg.scale(x.weight * y.weight);
                out.single = out.single + g;
                out.derivative_single = out.derivative_single + dg;
                for a in 0..3 {
                    for b in 0..3 {
                        let shape = x.shape[a] * y.shape[b];
                        out.shape[a][b] = out.shape[a][b] + g.scale(shape);
                        out.derivative_shape[a][b] = out.derivative_shape[a][b] + dg.scale(shape);
                    }
                }
            }
        }
        out
    }

    fn complex_coefficients(
        &self,
        k: C64,
        direction: Point,
        origin: Point,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<[Vec<C64>; 4], HelmholtzError> {
        let mut out = [
            zeroed(self.n, C64::ZERO)?,
            zeroed(self.n, C64::ZERO)?,
            zeroed(self.n, C64::ZERO)?,
            zeroed(self.n, C64::ZERO)?,
        ];
        for (f, samples) in self.facets.iter().zip(&self.samples) {
            poll(cancelled)?;
            let normal_dot = dot(f.normal, direction);
            for sample in samples {
                let u = dot(direction, sub(sample.point, origin));
                let terms = plane_factors(k, u, normal_dot);
                for a in 0..3 {
                    if let Some(v) = self.dof[f.ids[a]] {
                        let weight = sample.weight * sample.shape[a];
                        for c in 0..4 {
                            out[c][v] = out[c][v] + terms[c].scale(weight);
                        }
                    }
                }
            }
        }
        Ok(out)
    }
}

#[cfg(test)]
mod tests;
