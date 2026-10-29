//! Contour discovery of isolated finite poles, followed by original-equation
//! refinement. Beyn, LAA 436 (2012), 3839–3863, doi:10.1016/j.laa.2011.03.030.
//!
//! The two resolvent moments retain the original coordinates. A thin complex
//! SVD and the existing general QR owner extract seeds, and the shared bordered
//! Newton owner refines them. An independent argument-principle integral must
//! agree with the resolved moment rank. Nested quadrature differences are
//! numerical convergence evidence, NOT rigorous interval/count certificates.

use super::{NonlinearEigError, NonlinearPole, PoleEvaluation, PoleOptions};
use super::{admit, divide, finite, newton_pole_with_cancel, norm2, poll};
use crate::eigen_complex::{
    ComplexSvdError, EigError, LuComplexError, eig_with_cancel, lu_complex_with_cancel,
    svd_complex_with_cancel,
};
use fs_math::{c64::C64, det};

/// Dense contour integration is explicitly bounded, even with unlimited budgets.
pub const MAX_BEYN_NODES: usize = 4096;

/// Counterclockwise ellipse in the caller's spectral units. A circle has equal
/// radii. The operator must be holomorphic on AND throughout its interior.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Ellipse {
    pub center: C64,
    pub real_radius: f64,
    pub imaginary_radius: f64,
}

impl Ellipse {
    fn scale(self) -> f64 {
        self.real_radius.max(self.imaginary_radius)
    }

    pub fn contains(self, value: C64) -> bool {
        let delta = value - self.center;
        let x = delta.re / self.real_radius;
        let y = delta.im / self.imaginary_radius;
        x * x + y * y < 1.0
    }

    fn node(self, k: usize, nodes: usize) -> (C64, C64, C64) {
        let theta = core::f64::consts::TAU * (k as f64 / nodes as f64);
        let (sine, cosine) = det::sin_cos(theta);
        let delta = C64::new(self.real_radius * cosine, self.imaginary_radius * sine);
        // dz/(2*pi*i) for the trapezoid interval 2*pi/nodes.
        let weight = C64::new(self.imaginary_radius * cosine, self.real_radius * sine)
            .scale(1.0 / nodes as f64);
        (self.center + delta, divide(delta, self.scale()), weight)
    }
}

/// Explicit quadrature, resolution, conditioning and computational limits.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BeynOptions {
    /// Even fine-grid node count, 8..=MAX_BEYN_NODES. The even nodes also give
    /// the nested half-grid estimate without another operator evaluation.
    pub nodes: usize,
    /// 1..=n original-coordinate probes. Full n uses the identity; narrower
    /// probes use a fixed integer generator. A saturated narrower probe refuses.
    pub probe_columns: usize,
    /// Numerical rank threshold relative to the largest moment singular value.
    pub rank_relative_tolerance: f64,
    /// Bound on both integrals' distance from the same integer, their imaginary
    /// parts and their mutual difference. Strictly less than one quarter.
    pub count_tolerance: f64,
    /// Bound on nested moment differences relative to accumulated integrand
    /// magnitudes. This is an error indicator, not an enclosure.
    pub quadrature_tolerance: f64,
    /// Refuse when 1/(max(||T||_inf, residual_scale) ||T^-1||_inf) falls below
    /// this floor. The coefficient envelope also resolves scalar near-poles.
    /// This numerical margin is coordinate dependent, not an interval bound.
    pub contour_min_reciprocal_condition: f64,
    /// Inverse-column rowwise backward-residual ceiling:
    /// |A_i*x-b_i|/(|b_i|+||A_i||_1*||x||_inf). This remains meaningful for
    /// structurally zero inverse entries. Up to three corrective solves reuse
    /// the same LU. Physical consumers still check their original pole rows.
    pub linear_residual_tolerance: f64,
    pub svd_tolerance: f64,
    pub svd_max_sweeps: usize,
    pub eigen_max_sweeps: usize,
    pub max_dimension: usize,
    /// Conservative simultaneous numerical payload, excluding callback caches
    /// and source-owned storage, allocator overhead and process RSS.
    pub max_dense_bytes: usize,
    /// Conservative arithmetic bound for integration, SVD, QR and all local
    /// refinements. Callback-owned assembly work is additional.
    pub max_work: u64,
    pub refinement: PoleOptions,
}

#[derive(Debug, Clone)]
pub struct BeynReport {
    pub poles: Vec<NonlinearPole>,
    /// Values before local refinement, useful for a live quadrature ladder.
    pub unrefined_values: Vec<C64>,
    pub argument_integral: C64,
    pub coarse_argument_integral: C64,
    pub argument_count: usize,
    pub rank: usize,
    pub nodes: usize,
    pub moment_relative_change: f64,
    pub minimum_reciprocal_condition: f64,
}

#[derive(Debug, Clone, PartialEq)]
pub enum BeynError {
    Invalid {
        what: &'static str,
    },
    Operator(NonlinearEigError),
    ContourSingular {
        node: usize,
        index: usize,
    },
    ContourIllConditioned {
        node: usize,
        reciprocal_condition: f64,
    },
    LinearResidual {
        node: usize,
        relative: f64,
    },
    CountUnresolved {
        fine: C64,
        coarse: C64,
    },
    QuadratureUnresolved {
        relative_change: f64,
    },
    ProbeRankSaturated {
        probes: usize,
        count: usize,
    },
    CountMismatch {
        count: usize,
        rank: usize,
    },
    Svd(ComplexSvdError),
    Eigen(EigError),
    PoleOutside {
        value: C64,
    },
    UnresolvedCluster {
        first: C64,
        second: C64,
    },
}
impl core::fmt::Display for BeynError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "contour pole discovery refused: {self:?}")
    }
}
impl std::error::Error for BeynError {}
impl From<NonlinearEigError> for BeynError {
    fn from(error: NonlinearEigError) -> Self {
        Self::Operator(error)
    }
}

fn invalid(what: &'static str) -> BeynError {
    BeynError::Invalid { what }
}

fn preflight(n: usize, contour: Ellipse, opts: &BeynOptions) -> Result<(), BeynError> {
    admit(n, &opts.refinement)?;
    let l = opts.probe_columns;
    if n > opts.max_dimension
        || l == 0
        || l > n
        || opts.nodes < 8
        || opts.nodes > MAX_BEYN_NODES
        || opts.nodes % 2 != 0
        || opts.svd_max_sweeps == 0
        || opts.eigen_max_sweeps == 0
        || !finite(contour.center)
        || [contour.real_radius, contour.imaginary_radius]
            .iter()
            .any(|&x| !x.is_finite() || x <= 0.0)
        || !finite(contour.center + C64::new(contour.real_radius, contour.imaginary_radius))
        || !finite(contour.center - C64::new(contour.real_radius, contour.imaginary_radius))
        || contour.center.re + contour.real_radius == contour.center.re
        || contour.center.im + contour.imaginary_radius == contour.center.im
        || [
            opts.rank_relative_tolerance,
            opts.quadrature_tolerance,
            opts.contour_min_reciprocal_condition,
            opts.linear_residual_tolerance,
            opts.svd_tolerance,
        ]
        .iter()
        .any(|&x| !x.is_finite() || x <= 0.0 || x >= 1.0)
        || !opts.count_tolerance.is_finite()
        || opts.count_tolerance <= 0.0
        || opts.count_tolerance >= 0.25
        || opts.svd_tolerance <= f64::EPSILON
    {
        return Err(invalid(
            "invalid contour, resolution, accuracy or dimension budget",
        ));
    }
    // u128 keeps every checked dimension product exact before comparing limits.
    let (nn, ll, nodes) = (n as u128, l as u128, opts.nodes as u128);
    let d = nn + 1;
    let bytes = 16 * (18 * nn * nn + 24 * nn * ll + 12 * ll * ll + 128 * d);
    let work = nodes * (128 * nn * nn * nn + 128 * nn * nn)
        + opts.svd_max_sweeps as u128 * (32 * nn * ll * ll + 32 * ll * ll * ll * ll)
        + opts.eigen_max_sweeps as u128 * 32 * ll * ll * ll * ll
        + ll * opts.refinement.max_iterations as u128 * d * d * (d + 16);
    if bytes > opts.max_dense_bytes as u128 || work > opts.max_work as u128 {
        return Err(invalid(
            "insufficient contour numerical payload or work budget",
        ));
    }
    Ok(())
}

fn add_compensated(sum: &mut C64, correction: &mut C64, term: C64) {
    let adjusted = term - *correction;
    let next = *sum + adjusted;
    *correction = (next - *sum) - adjusted;
    *sum = next;
}

fn matrix_inf(a: &[C64], n: usize) -> f64 {
    a.chunks_exact(n)
        .map(|row| row.iter().map(|z| z.abs()).sum::<f64>())
        .fold(0.0, f64::max)
}

fn make_probes(n: usize, l: usize) -> Vec<C64> {
    let mut probe = vec![C64::ZERO; n * l];
    if n == l {
        for i in 0..n {
            probe[i * l + i] = C64::ONE;
        }
    } else {
        let mut state = 0x5e4d_2c1b_9a87_630f_u64;
        for x in &mut probe {
            let mut next = || {
                state = state
                    .wrapping_mul(6364136223846793005)
                    .wrapping_add(1442695040888963407);
                (state >> 11) as f64 * (1.0 / 9007199254740992.0) - 0.5
            };
            *x = C64::new(next(), next());
        }
    }
    probe
}

fn contour_lu(error: LuComplexError, node: usize) -> BeynError {
    match error {
        LuComplexError::Cancelled => NonlinearEigError::Cancelled.into(),
        LuComplexError::Singular { index } => BeynError::ContourSingular { node, index },
    }
}

/// Find the resolved simple finite poles inside a caller-declared analytic
/// region, without individual eigenvalue/vector seeds.
///
/// The callback contract is the same as `newton_pole_with_cancel`. T must be
/// holomorphic throughout the ellipse; a meromorphic operator's poles would
/// invalidate the eigenvalue count. This assumption is the producer's duty.
/// Two moments cannot resolve linearly dependent eigendirections or more roots
/// than the probe rank. Count/rank mismatch refuses rather than dropping modes.
/// Repeated/defective poles and unresolved near-coincident roots can refuse.
/// Quadrature agreement is sampled evidence and does not prove completeness.
pub fn beyn_with_cancel<F, C>(
    n: usize,
    contour: Ellipse,
    opts: &BeynOptions,
    mut evaluate: F,
    mut cancelled: C,
) -> Result<BeynReport, BeynError>
where
    F: FnMut(C64, &mut dyn FnMut() -> bool) -> Result<PoleEvaluation, NonlinearEigError>,
    C: FnMut() -> bool,
{
    poll(&mut cancelled)?;
    preflight(n, contour, opts)?;
    let l = opts.probe_columns;
    let probe = make_probes(n, l);
    let mut moments = vec![vec![C64::ZERO; n * l]; 4];
    let mut corrections = vec![vec![C64::ZERO; n * l]; 4];
    let mut counts = [C64::ZERO; 2];
    let mut count_corrections = [C64::ZERO; 2];
    let mut magnitude = 0.0;
    let mut minimum_rcond = 1.0_f64;
    let mut fixed_scale = 0.0;
    for k in 0..opts.nodes {
        poll(&mut cancelled)?;
        let (z, coordinate, weight) = contour.node(k, opts.nodes);
        if !finite(z) || !finite(weight) {
            return Err(invalid("contour nodes are not representable"));
        }
        let value = evaluate(z, &mut cancelled)?;
        poll(&mut cancelled)?;
        if value.matrix.len() != n * n
            || value.derivative.len() != n * n
            || value
                .matrix
                .iter()
                .chain(&value.derivative)
                .any(|&a| !finite(a))
            || !value.residual_scale.is_finite()
            || value.residual_scale <= 0.0
        {
            return Err(invalid(
                "finite original operator, derivative and residual scale required",
            ));
        }
        if k == 0 {
            fixed_scale = value.matrix.iter().map(|z| z.abs()).fold(0.0, f64::max);
            if !fixed_scale.is_finite() || fixed_scale == 0.0 {
                return Err(BeynError::ContourSingular { node: k, index: 0 });
            }
        }
        // One CONSTANT normalization over the entire contour preserves both
        // analytic moments. A different unaccounted node scale would not.
        let a: Vec<_> = value
            .matrix
            .iter()
            .map(|&x| divide(x, fixed_scale))
            .collect();
        let derivative: Vec<_> = value
            .derivative
            .iter()
            .map(|&x| divide(x, fixed_scale))
            .collect();
        if a.iter().chain(&derivative).any(|&x| !finite(x))
            || value
                .matrix
                .iter()
                .chain(&value.derivative)
                .zip(a.iter().chain(&derivative))
                .any(|(original, scaled)| {
                    (original.re != 0.0 && scaled.re == 0.0)
                        || (original.im != 0.0 && scaled.im == 0.0)
                })
        {
            return Err(invalid("normalized contour operator is not representable"));
        }
        let lu = lu_complex_with_cancel(&a, n, &mut cancelled).map_err(|e| contour_lu(e, k))?;
        let mut inverse = vec![C64::ZERO; n * n];
        let mut rhs = vec![C64::ZERO; n];
        let mut correction = vec![C64::ZERO; n];
        for j in 0..n {
            poll(&mut cancelled)?;
            rhs.fill(C64::ZERO);
            rhs[j] = C64::ONE;
            lu.solve_with_cancel(&mut rhs, &mut cancelled)
                .map_err(|e| contour_lu(e, k))?;
            for refinement in 0..=3 {
                let mut relative = 0.0_f64;
                let solution_norm = rhs.iter().map(|z| z.abs()).fold(0.0, f64::max);
                for i in 0..n {
                    poll(&mut cancelled)?;
                    let b = if i == j { C64::ONE } else { C64::ZERO };
                    let mut residual = -b;
                    let mut residual_correction = C64::ZERO;
                    let mut scale = b.abs();
                    for t in 0..n {
                        add_compensated(
                            &mut residual,
                            &mut residual_correction,
                            a[i * n + t] * rhs[t],
                        );
                        scale += a[i * n + t].abs() * solution_norm;
                    }
                    correction[i] = -residual;
                    let row = if scale == 0.0 {
                        residual.abs()
                    } else {
                        residual.abs() / scale
                    };
                    if !row.is_finite() || !scale.is_finite() {
                        return Err(BeynError::LinearResidual {
                            node: k,
                            relative: f64::INFINITY,
                        });
                    }
                    relative = relative.max(row);
                }
                if relative <= opts.linear_residual_tolerance {
                    break;
                }
                if refinement == 3 {
                    return Err(BeynError::LinearResidual { node: k, relative });
                }
                // Correct the actual inverse column in the original rows.
                // No small entry is clipped or replaced by a structural guess.
                lu.solve_with_cancel(&mut correction, &mut cancelled)
                    .map_err(|e| contour_lu(e, k))?;
                for i in 0..n {
                    rhs[i] = rhs[i] + correction[i];
                }
            }
            for i in 0..n {
                inverse[i * n + j] = rhs[i];
            }
        }
        if inverse.iter().any(|&x| !finite(x)) {
            return Err(invalid("contour inverse is not representable"));
        }
        let coefficient_norm = value.residual_scale / fixed_scale;
        let rcond = 1.0 / (matrix_inf(&a, n).max(coefficient_norm) * matrix_inf(&inverse, n));
        if !rcond.is_finite() || rcond < opts.contour_min_reciprocal_condition {
            return Err(BeynError::ContourIllConditioned {
                node: k,
                reciprocal_condition: rcond,
            });
        }
        minimum_rcond = minimum_rcond.min(rcond);
        let mut trace = C64::ZERO;
        let mut trace_correction = C64::ZERO;
        for i in 0..n {
            poll(&mut cancelled)?;
            for j in 0..n {
                add_compensated(
                    &mut trace,
                    &mut trace_correction,
                    inverse[i * n + j] * derivative[j * n + i],
                );
            }
        }
        add_compensated(&mut counts[0], &mut count_corrections[0], trace * weight);
        if k % 2 == 0 {
            add_compensated(
                &mut counts[1],
                &mut count_corrections[1],
                (trace * weight).scale(2.0),
            );
        }
        let mut node_moment = vec![C64::ZERO; n * l];
        for i in 0..n {
            poll(&mut cancelled)?;
            for j in 0..l {
                let mut sum = C64::ZERO;
                for t in 0..n {
                    sum = sum + inverse[i * n + t] * probe[t * l + j];
                }
                let at = i * l + j;
                node_moment[at] = sum * weight;
                let terms = [node_moment[at], node_moment[at] * coordinate];
                for p in 0..2 {
                    add_compensated(&mut moments[p][at], &mut corrections[p][at], terms[p]);
                    if k % 2 == 0 {
                        add_compensated(
                            &mut moments[p + 2][at],
                            &mut corrections[p + 2][at],
                            terms[p].scale(2.0),
                        );
                    }
                }
            }
        }
        magnitude += norm2(&node_moment) * (1.0 + coordinate.abs());
    }
    if counts.iter().any(|&x| !finite(x))
        || !magnitude.is_finite()
        || magnitude == 0.0
        || moments.iter().flatten().any(|&x| !finite(x))
    {
        return Err(invalid("contour integrals are not representable"));
    }
    let rounded = counts[0].re.round();
    let count_ok = rounded >= 0.0
        && rounded <= n as f64
        && counts.iter().all(|z| {
            (z.re - rounded).abs() <= opts.count_tolerance && z.im.abs() <= opts.count_tolerance
        })
        && (counts[0] - counts[1]).abs() <= opts.count_tolerance;
    if !count_ok {
        return Err(BeynError::CountUnresolved {
            fine: counts[0],
            coarse: counts[1],
        });
    }
    let count = rounded as usize;
    let difference: Vec<_> = (0..2)
        .flat_map(|p| moments[p].iter().zip(&moments[p + 2]).map(|(&a, &b)| a - b))
        .collect();
    let change = norm2(&difference) / magnitude;
    if !change.is_finite() || change > opts.quadrature_tolerance {
        return Err(BeynError::QuadratureUnresolved {
            relative_change: change,
        });
    }
    let mut report = BeynReport {
        poles: Vec::new(),
        unrefined_values: Vec::new(),
        argument_integral: counts[0],
        coarse_argument_integral: counts[1],
        argument_count: count,
        rank: 0,
        nodes: opts.nodes,
        moment_relative_change: change,
        minimum_reciprocal_condition: minimum_rcond,
    };
    if count == 0 {
        let zero_error = det::hypot(norm2(&moments[0]), norm2(&moments[1])) / magnitude;
        if zero_error > opts.quadrature_tolerance {
            return Err(BeynError::QuadratureUnresolved {
                relative_change: zero_error,
            });
        }
        poll(&mut cancelled)?;
        return Ok(report);
    }
    if count > l || (count == l && l < n) {
        return Err(BeynError::ProbeRankSaturated { probes: l, count });
    }
    let svd = svd_complex_with_cancel(
        &moments[0],
        n,
        l,
        opts.svd_tolerance,
        opts.svd_max_sweeps,
        &mut cancelled,
    )
    .map_err(BeynError::Svd)?;
    let threshold = svd.singular_values[0] * opts.rank_relative_tolerance;
    let rank = svd
        .singular_values
        .iter()
        .take_while(|&&s| s > threshold)
        .count();
    if rank != count {
        return Err(BeynError::CountMismatch { count, rank });
    }
    report.rank = rank;
    let mut projected = vec![C64::ZERO; n * rank];
    for i in 0..n {
        poll(&mut cancelled)?;
        for j in 0..rank {
            for t in 0..l {
                projected[i * rank + j] =
                    projected[i * rank + j] + moments[1][i * l + t] * svd.v[t * l + j];
            }
            projected[i * rank + j] = divide(projected[i * rank + j], svd.singular_values[j]);
        }
    }
    let mut small = vec![C64::ZERO; rank * rank];
    for i in 0..rank {
        poll(&mut cancelled)?;
        for j in 0..rank {
            for t in 0..n {
                small[i * rank + j] =
                    small[i * rank + j] + svd.u[t * l + i].conj() * projected[t * rank + j];
            }
        }
    }
    let values = eig_with_cancel(&small, rank, opts.eigen_max_sweeps, &mut cancelled)
        .map_err(BeynError::Eigen)?;
    for coordinate in values {
        poll(&mut cancelled)?;
        let value = contour.center + coordinate.scale(contour.scale());
        if !finite(value) || !contour.contains(value) {
            return Err(BeynError::PoleOutside { value });
        }
        let mut shifted = small.clone();
        for i in 0..rank {
            shifted[i * rank + i] = shifted[i * rank + i] - coordinate;
        }
        let null = svd_complex_with_cancel(
            &shifted,
            rank,
            rank,
            opts.svd_tolerance,
            opts.svd_max_sweeps,
            &mut cancelled,
        )
        .map_err(BeynError::Svd)?;
        let mut seed = vec![C64::ZERO; n];
        for i in 0..n {
            poll(&mut cancelled)?;
            for j in 0..rank {
                seed[i] = seed[i] + svd.u[i * l + j] * null.v[j * rank + rank - 1];
            }
        }
        let pole = newton_pole_with_cancel(
            n,
            value,
            &seed,
            &opts.refinement,
            &mut evaluate,
            &mut cancelled,
        )?;
        if !contour.contains(pole.value) {
            return Err(BeynError::PoleOutside { value: pole.value });
        }
        for old in &report.poles {
            let separation = (old.value - pole.value).abs();
            let tolerance = 8.0
                * opts.refinement.step_tolerance
                * opts
                    .refinement
                    .spectral_scale
                    .max(old.value.abs())
                    .max(pole.value.abs());
            if separation <= tolerance {
                return Err(BeynError::UnresolvedCluster {
                    first: old.value,
                    second: pole.value,
                });
            }
        }
        report.unrefined_values.push(value);
        report.poles.push(pole);
    }
    report.poles.sort_by(|a, b| {
        a.value
            .re
            .total_cmp(&b.value.re)
            .then_with(|| a.value.im.total_cmp(&b.value.im))
    });
    report
        .unrefined_values
        .sort_by(|a, b| a.re.total_cmp(&b.re).then_with(|| a.im.total_cmp(&b.im)));
    poll(&mut cancelled)?;
    Ok(report)
}

#[cfg(test)]
#[path = "beyn_tests.rs"]
mod tests;
