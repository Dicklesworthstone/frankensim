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
use kernel::{cosine_remainder, cosine_remainder_value, plane_factors, plane_values};

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

/// Values-only harmonic assembly of the identical finite complex operator.
/// No derivative matrix, derivative pair moments or derivative coefficient
/// vectors are formed. Geometry admission is shared with the pole chart.
#[derive(Debug, Clone, Copy)]
pub struct ScreenComplexValueLimits {
    /// Complex entries in the single returned H, namely n*n.
    pub maximum_matrix_entries: usize,
    /// H, bounded parallel pair scratch, angular rule and coefficient vectors.
    /// Excludes the retained geometry, factorization and downstream sources.
    pub maximum_workspace_bytes: usize,
    /// Four times the original assembly work, excluding its n^3 LU term.
    pub maximum_work: u64,
    /// The same centered |Im(k)|*2R chart as ScreenComplexLimits, in (0,8].
    pub maximum_growth_exponent: f64,
}

impl ScreenComplexLimits {
    fn assembly_limits(self) -> ScreenComplexValueLimits {
        ScreenComplexValueLimits {
            maximum_matrix_entries: self.maximum_matrix_entries,
            maximum_workspace_bytes: self.maximum_workspace_bytes,
            maximum_work: self.maximum_work,
            maximum_growth_exponent: self.maximum_growth_exponent,
        }
    }
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
struct Pair<const N: usize> {
    single: [C64; N],
    shape: [[[C64; 3]; 3]; N],
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
        self.complex_admission::<2>(
            maximum_abs_k,
            maximum_abs_imag_k,
            limits.assembly_limits(),
            &mut cancelled,
        )?;
        Ok(())
    }

    fn complex_admission<const N: usize>(
        &self,
        maximum_abs_k: f64,
        maximum_abs_imag_k: f64,
        limits: ScreenComplexValueLimits,
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
        let entries = N as u128 * n * n;
        let directions = 64_u128 << (2 * self.options.angular_refinement);
        let pair_bytes = if workers > 1 {
            block as u128 * nf as u128 * core::mem::size_of::<Pair<N>>() as u128
                + 128 * (block + workers) as u128
        } else {
            0
        };
        let bytes = entries * core::mem::size_of::<C64>() as u128
            + pair_bytes
            + (2 * N as u128 + 1) * n * core::mem::size_of::<C64>() as u128
            + directions * core::mem::size_of::<(Point, f64)>() as u128
            + 4096;
        let work = 4 * N as u128 * (self.frequency_work() - n * n * n);
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
        cancelled: impl FnMut() -> bool,
    ) -> Result<ScreenComplexOperator, HelmholtzError> {
        self.assemble_complex::<2>(k, limits.assembly_limits(), cancelled)
    }

    /// Assemble only H on the same centered cosine/antipodal rule. Positive
    /// real k retains the original real-frequency matrix arithmetic exactly.
    /// This is a finite frequency response, not a causal or passive realization.
    pub fn complex_matrix(
        &self,
        k: C64,
        limits: ScreenComplexValueLimits,
        cancelled: impl FnMut() -> bool,
    ) -> Result<Vec<C64>, HelmholtzError> {
        self.assemble_complex::<1>(k, limits, cancelled)
            .map(|op| op.matrix)
    }

    fn assemble_complex<const N: usize>(
        &self,
        k: C64,
        limits: ScreenComplexValueLimits,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<ScreenComplexOperator, HelmholtzError> {
        poll(&mut cancelled)?;
        if !finite(k) {
            return Err(bad("complex open-screen wavenumber is nonfinite"));
        }
        let chart = self.complex_admission::<N>(k.abs(), k.im.abs(), limits, &mut cancelled)?;
        let real_axis = k.im == 0.0;
        let k = if real_axis { C64::from_re(k.re) } else { k };
        let k2 = k * k;
        let n = self.n;
        let mut matrix = zeroed(n * n, C64::ZERO)?;
        let mut derivative = if N == 2 {
            zeroed(n * n, C64::ZERO)?
        } else {
            Vec::new()
        };
        for i in 0..n * n {
            if i % n == 0 {
                poll(&mut cancelled)?;
            }
            matrix[i] = C64::from_re(self.static_curl[i]) + k2.scale(self.static_normal[i]);
            if N == 2 {
                derivative[i] = k.scale(2.0 * self.static_normal[i]);
            }
        }
        let mut add = |i, j, pair: Pair<N>| {
            self.add_dynamic(&mut matrix, k2, i, j, pair.single[0], pair.shape[0]);
            if N == 2 {
                self.add_dynamic(&mut derivative, k2, i, j, pair.single[1], pair.shape[1]);
                self.add_dynamic(
                    &mut derivative,
                    k.scale(2.0),
                    i,
                    j,
                    C64::ZERO,
                    pair.shape[0],
                );
            }
        };
        let nf = self.facets.len();
        if chart.workers <= 1 {
            for i in 0..nf {
                for j in i..nf {
                    poll(&mut cancelled)?;
                    add(i, j, self.complex_pair::<N>(k, real_axis, i, j));
                }
            }
        } else {
            for first in (0..nf).step_by(chart.block) {
                poll(&mut cancelled)?;
                let rows = first..nf.min(first + chart.block);
                let pairs = rows_on_workers(rows.clone(), chart.workers, |i| {
                    (i..nf)
                        .map(|j| self.complex_pair::<N>(k, real_axis, i, j))
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
            let coefficients = if real_axis && N == 1 {
                core::array::from_fn(|_| Vec::new())
            } else {
                self.complex_coefficients::<N>(k, direction, chart.origin, &mut cancelled)?
            };
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
                    let at = i * n + j;
                    if real_axis {
                        // The exact original outgoing Gram arithmetic and sum
                        // order, including its absolute-coordinate phase.
                        let scale = k.re * weight / (8.0 * PI * PI);
                        matrix[at].im += scale * (original[i].conj() * original[j]).re;
                    }
                    if !real_axis || N == 2 {
                        let product =
                            (backward[i] * forward[j] + forward[i] * backward[j]).scale(0.5);
                        if !real_axis {
                            matrix[at] = matrix[at] + (factor * k) * product;
                        }
                        if N == 2 {
                            let d_product = (d_backward[i] * forward[j]
                                + backward[i] * d_forward[j]
                                + d_forward[i] * backward[j]
                                + forward[i] * d_backward[j])
                                .scale(0.5);
                            derivative[at] = derivative[at] + factor * (product + k * d_product);
                        }
                    }
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

    fn complex_pair<const N: usize>(&self, k: C64, real_axis: bool, i: usize, j: usize) -> Pair<N> {
        let mut out = Pair {
            single: [C64::ZERO; N],
            shape: [[[C64::ZERO; 3]; 3]; N],
        };
        for x in &self.samples[i] {
            for y in &self.samples[j] {
                let r = norm(sub(x.point, y.point));
                let (g, dg) = if real_axis {
                    let mut g = green_remainder(k.re, r);
                    g.im = 0.0;
                    (
                        g,
                        if N == 2 {
                            C64::from_re(-det::sin(k.re * r) / (4.0 * PI))
                        } else {
                            C64::ZERO
                        },
                    )
                } else if N == 1 {
                    (cosine_remainder_value(k, r), C64::ZERO)
                } else {
                    cosine_remainder(k, r)
                };
                let g = g.scale(x.weight * y.weight);
                out.single[0] = out.single[0] + g;
                let dg = if N == 2 {
                    dg.scale(x.weight * y.weight)
                } else {
                    C64::ZERO
                };
                if N == 2 {
                    out.single[1] = out.single[1] + dg;
                }
                for a in 0..3 {
                    for b in 0..3 {
                        let shape = x.shape[a] * y.shape[b];
                        out.shape[0][a][b] = out.shape[0][a][b] + g.scale(shape);
                        if N == 2 {
                            out.shape[1][a][b] = out.shape[1][a][b] + dg.scale(shape);
                        }
                    }
                }
            }
        }
        out
    }

    fn complex_coefficients<const N: usize>(
        &self,
        k: C64,
        direction: Point,
        origin: Point,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<[Vec<C64>; 4], HelmholtzError> {
        let mut out = [
            zeroed(self.n, C64::ZERO)?,
            zeroed(self.n, C64::ZERO)?,
            if N == 2 {
                zeroed(self.n, C64::ZERO)?
            } else {
                Vec::new()
            },
            if N == 2 {
                zeroed(self.n, C64::ZERO)?
            } else {
                Vec::new()
            },
        ];
        for (f, samples) in self.facets.iter().zip(&self.samples) {
            poll(cancelled)?;
            let normal_dot = dot(f.normal, direction);
            for sample in samples {
                let u = dot(direction, sub(sample.point, origin));
                let terms = if N == 2 {
                    plane_factors(k, u, normal_dot)
                } else {
                    let values = plane_values(k, u, normal_dot);
                    [values[0], values[1], C64::ZERO, C64::ZERO]
                };
                for a in 0..3 {
                    if let Some(v) = self.dof[f.ids[a]] {
                        let weight = sample.weight * sample.shape[a];
                        for c in 0..2 * N {
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
