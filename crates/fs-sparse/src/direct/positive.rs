//! Positive scalar-pivot whitening of the existing sparse LDLᵀ factor.
//!
//! With P A Pᵀ = L D Lᵀ, S = Pᵀ L⁻ᵀ D⁻¹/² satisfies Sᵀ A S = I.
//! The two actions are transposes in the Euclidean metric. They preserve
//! symmetry of Sᵀ H S for an indefinite symmetric H; A⁻¹H generally does not.
//! No factor perturbation or implicit treatment of a 2×2 D block is allowed.

use super::{DBlock, LdltFactor};

/// A whitening admission or action failed.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PositiveFactorError {
    /// The sparse factor has negative or missing positive pivots.
    NotPositive,
    /// This bounded lane accepts scalar D only; prepare with 2×2 disabled.
    BlockPivotsUnsupported,
    /// A factor, action, or inverse square root is not finite.
    NonFinite,
    /// Input, output, or caller-owned workspace has the wrong length.
    InvalidShape,
    /// Cancellation observed before publishing an action.
    Cancelled,
}
impl core::fmt::Display for PositiveFactorError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "positive sparse whitening: {self:?}")
    }
}
impl std::error::Error for PositiveFactorError {}

/// An immutable admitted positive factor and its scalar inverse square roots.
/// Preparation uses O(n) extra storage. Each action uses caller-owned O(n)
/// scratch and O(nnz(L)) arithmetic, with cancellation between sparse columns.
pub struct PositiveWhitening<'a> {
    factor: &'a LdltFactor,
    inverse_sqrt: Vec<f64>,
}

impl<'a> PositiveWhitening<'a> {
    /// Admit the complete factor, not only sampled positive Rayleigh quotients.
    /// An SPD factor with 2×2 pivots is explicitly unsupported in this lane.
    pub fn new(
        factor: &'a LdltFactor,
        cancelled: &impl Fn() -> bool,
    ) -> Result<Self, PositiveFactorError> {
        poll(cancelled)?;
        if factor.inertia.negative != 0 || factor.inertia.positive != factor.n {
            return Err(PositiveFactorError::NotPositive);
        }
        if factor.stats.pivots_2x2 != 0 {
            return Err(PositiveFactorError::BlockPivotsUnsupported);
        }
        let mut inverse_sqrt = vec![0.0; factor.n];
        for node in &factor.snodes {
            poll(cancelled)?;
            if node.l.iter().any(|v| !v.is_finite()) {
                return Err(PositiveFactorError::NonFinite);
            }
            for (k, block) in node.d.iter().enumerate() {
                let DBlock::Single(diagonal) = *block else {
                    return Err(PositiveFactorError::BlockPivotsUnsupported);
                };
                if !diagonal.is_finite() || diagonal <= 0.0 {
                    return Err(PositiveFactorError::NotPositive);
                }
                inverse_sqrt[node.block_vars[k]] = 1.0 / diagonal.sqrt();
            }
        }
        if inverse_sqrt.iter().any(|v| !v.is_finite() || *v <= 0.0) {
            return Err(PositiveFactorError::NonFinite);
        }
        poll(cancelled)?;
        Ok(Self {
            factor,
            inverse_sqrt,
        })
    }

    /// Matrix dimension.
    #[must_use]
    pub fn dim(&self) -> usize {
        self.factor.n
    }

    /// S x: whitened coordinates to original variable order.
    /// On failure the output is unspecified; callers must not publish it.
    pub fn apply(
        &self,
        x: &[f64],
        y: &mut [f64],
        scratch: &mut [f64],
        cancelled: &impl Fn() -> bool,
    ) -> Result<(), PositiveFactorError> {
        crate::fma::positive_apply_dispatch(self, x, y, scratch, cancelled)
    }

    // The same safe arithmetic is inlined into the registered FMA capsule.
    #[inline(always)]
    pub(crate) fn apply_body(
        &self,
        x: &[f64],
        y: &mut [f64],
        scratch: &mut [f64],
        cancelled: &dyn Fn() -> bool,
    ) -> Result<(), PositiveFactorError> {
        self.shape(x, y, scratch)?;
        poll(cancelled)?;
        for i in 0..self.dim() {
            scratch[i] = x[i] * self.inverse_sqrt[i];
        }
        for node in self.factor.snodes.iter().rev() {
            let w = node.block_vars.len();
            let m = w + node.panel_vars.len();
            for k in (0..w).rev() {
                poll(cancelled)?;
                let mut acc = scratch[node.block_vars[k]];
                for i in k + 1..w {
                    acc = node.l[i + k * m].mul_add(-scratch[node.block_vars[i]], acc);
                }
                for (i, &variable) in node.panel_vars.iter().enumerate() {
                    acc = node.l[w + i + k * m].mul_add(-scratch[variable], acc);
                }
                scratch[node.block_vars[k]] = acc;
            }
        }
        for (k, &original) in self.factor.perm.iter().enumerate() {
            y[original] = scratch[k];
        }
        self.finish(y, cancelled)
    }

    /// Sᵀ x: original variable order to whitened coordinates.
    /// On failure the output is unspecified; callers must not publish it.
    pub fn apply_transpose(
        &self,
        x: &[f64],
        y: &mut [f64],
        scratch: &mut [f64],
        cancelled: &impl Fn() -> bool,
    ) -> Result<(), PositiveFactorError> {
        crate::fma::positive_transpose_dispatch(self, x, y, scratch, cancelled)
    }

    #[inline(always)]
    pub(crate) fn transpose_body(
        &self,
        x: &[f64],
        y: &mut [f64],
        scratch: &mut [f64],
        cancelled: &dyn Fn() -> bool,
    ) -> Result<(), PositiveFactorError> {
        self.shape(x, y, scratch)?;
        poll(cancelled)?;
        for (k, &original) in self.factor.perm.iter().enumerate() {
            scratch[k] = x[original];
        }
        for node in &self.factor.snodes {
            let w = node.block_vars.len();
            let m = w + node.panel_vars.len();
            for k in 0..w {
                poll(cancelled)?;
                let vk = scratch[node.block_vars[k]];
                for i in k + 1..w {
                    let value = &mut scratch[node.block_vars[i]];
                    *value = node.l[i + k * m].mul_add(-vk, *value);
                }
                for (i, &variable) in node.panel_vars.iter().enumerate() {
                    let value = &mut scratch[variable];
                    *value = node.l[w + i + k * m].mul_add(-vk, *value);
                }
            }
        }
        for i in 0..self.dim() {
            y[i] = scratch[i] * self.inverse_sqrt[i];
        }
        self.finish(y, cancelled)
    }

    fn shape(&self, x: &[f64], y: &[f64], scratch: &[f64]) -> Result<(), PositiveFactorError> {
        if x.len() != self.dim() || y.len() != self.dim() || scratch.len() != self.dim() {
            return Err(PositiveFactorError::InvalidShape);
        }
        if x.iter().any(|v| !v.is_finite()) {
            return Err(PositiveFactorError::NonFinite);
        }
        Ok(())
    }

    fn finish(&self, y: &[f64], cancelled: &dyn Fn() -> bool) -> Result<(), PositiveFactorError> {
        poll(cancelled)?;
        if y.iter().any(|v| !v.is_finite()) {
            return Err(PositiveFactorError::NonFinite);
        }
        Ok(())
    }
}

fn poll(cancelled: &dyn Fn() -> bool) -> Result<(), PositiveFactorError> {
    if cancelled() {
        Err(PositiveFactorError::Cancelled)
    } else {
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::super::{DirectLimits, LdltError};
    use super::*;
    use crate::{Coo, Csr, DirectOrdering, LdltOptions, SymbolicLdlt};
    use std::cell::Cell;

    fn matrix() -> Csr {
        let n = 9;
        let mut c = Coo::new(n, n);
        for i in 0..n {
            c.push(i, i, 3.0 + 0.1 * i as f64);
        }
        for i in 1..n {
            c.push(0, i, -0.3);
            c.push(i, 0, -0.3);
        }
        c.assemble()
    }
    fn limits() -> DirectLimits {
        DirectLimits {
            max_symbolic_entries: 1000,
            max_analysis_work: 1_000_000,
            max_numeric_work: 1_000_000,
            max_front_dim: 8,
        }
    }
    fn close(a: &[f64], b: &[f64]) {
        for (&x, &y) in a.iter().zip(b) {
            assert!((x - y).abs() < 2e-12, "{x} != {y}");
        }
    }

    #[test]
    fn positive_whitening_is_a_transpose_pair_in_original_permuted_coordinates() {
        let a = matrix();
        let n = a.nrows();
        let symbolic =
            SymbolicLdlt::analyze_bounded(&a, DirectOrdering::Amd, &limits(), &|| false).unwrap();
        assert_ne!(symbolic.permutation(), &(0..n).collect::<Vec<_>>());
        let options = LdltOptions {
            allow_2x2: false,
            max_factor_bytes: Some(100_000),
            ..LdltOptions::default()
        };
        let factor = symbolic
            .factor_bounded(&a, &options, &limits(), &|| false)
            .unwrap();
        let old = SymbolicLdlt::analyze(&a, DirectOrdering::Amd)
            .unwrap()
            .factor(&a, &options)
            .unwrap();
        let rhs: Vec<_> = (0..n).map(|i| 1.0 + i as f64 * 0.3).collect();
        assert_eq!(factor.solve(&rhs), old.solve(&rhs));
        let s = PositiveWhitening::new(&factor, &|| false).unwrap();
        let mut scratch = vec![0.; n];
        let mut physical = vec![0.; n];
        let mut force = vec![0.; n];
        let mut whitened = vec![0.; n];
        for j in 0..n {
            let mut unit = vec![0.; n];
            unit[j] = 1.;
            s.apply(&unit, &mut physical, &mut scratch, &|| false)
                .unwrap();
            a.spmv(&physical, &mut force);
            s.apply_transpose(&force, &mut whitened, &mut scratch, &|| false)
                .unwrap();
            close(&whitened, &unit);
        }
        s.apply_transpose(&rhs, &mut whitened, &mut scratch, &|| false)
            .unwrap();
        s.apply(&whitened, &mut physical, &mut scratch, &|| false)
            .unwrap();
        close(&physical, &factor.solve(&rhs));
    }

    #[test]
    fn cancellation_and_all_sparse_preparation_caps_refuse() {
        let a = matrix();
        assert!(matches!(
            SymbolicLdlt::analyze_bounded(
                &a,
                DirectOrdering::NestedDissection(crate::NdOptions {
                    leaf_size: 64,
                    refinement_passes: 4
                }),
                &DirectLimits::UNLIMITED,
                &|| false
            ),
            Err(LdltError::InvalidOptions { .. })
        ));
        let symbolic =
            SymbolicLdlt::analyze_bounded(&a, DirectOrdering::Amd, &limits(), &|| false).unwrap();
        let options = LdltOptions {
            allow_2x2: false,
            ..LdltOptions::default()
        };
        for stop in [1, 3, 12] {
            let polls = Cell::new(0);
            let answer = SymbolicLdlt::analyze_bounded(&a, DirectOrdering::Amd, &limits(), &|| {
                polls.set(polls.get() + 1);
                polls.get() >= stop
            });
            assert!(matches!(answer, Err(LdltError::Cancelled)));
        }
        for stop in [1, 3, 12] {
            let polls = Cell::new(0);
            let answer = symbolic.factor_bounded(&a, &options, &limits(), &|| {
                polls.set(polls.get() + 1);
                polls.get() >= stop
            });
            assert!(matches!(answer, Err(LdltError::Cancelled)));
        }
        let mut cap = limits();
        cap.max_symbolic_entries = 1;
        assert!(matches!(
            SymbolicLdlt::analyze_bounded(&a, DirectOrdering::Amd, &cap, &|| false),
            Err(LdltError::SymbolicBudget { .. })
        ));
        let mut cap = limits();
        cap.max_analysis_work = 1;
        assert!(matches!(
            SymbolicLdlt::analyze_bounded(&a, DirectOrdering::Amd, &cap, &|| false),
            Err(LdltError::WorkLimit { .. })
        ));
        let mut cap = limits();
        cap.max_numeric_work = 1;
        assert!(matches!(
            symbolic.factor_bounded(&a, &options, &cap, &|| false),
            Err(LdltError::WorkLimit { .. })
        ));
        let mut cap = limits();
        cap.max_front_dim = 1;
        assert!(matches!(
            symbolic.factor_bounded(&a, &options, &cap, &|| false),
            Err(LdltError::FrontBudget { .. })
        ));
        let factor = symbolic.factor(&a, &options).unwrap();
        let s = PositiveWhitening::new(&factor, &|| false).unwrap();
        let mut out = vec![0.; 9];
        let mut scratch = vec![0.; 9];
        assert_eq!(
            s.apply(&[1.; 9], &mut out, &mut scratch, &|| true),
            Err(PositiveFactorError::Cancelled)
        );
    }

    #[test]
    fn nonpositive_and_positive_block_factors_are_explicit_refusals() {
        let mut c = Coo::new(2, 2);
        c.push(0, 0, -1.);
        c.push(1, 1, 2.);
        let a = c.assemble();
        let factor = SymbolicLdlt::analyze(&a, DirectOrdering::Natural)
            .unwrap()
            .factor(&a, &LdltOptions::default())
            .unwrap();
        assert!(matches!(
            PositiveWhitening::new(&factor, &|| false),
            Err(PositiveFactorError::NotPositive)
        ));
        // Exact hand-known LDL factor of [[2,1],[1,2]], L=I and one SPD D
        // block. It must refuse; treating block diagonals as scalar D would
        // produce a mathematically incorrect whitening.
        let mut c = Coo::new(2, 2);
        for (i, j, v) in [(0, 0, 2.), (0, 1, 1.), (1, 0, 1.), (1, 1, 2.)] {
            c.push(i, j, v);
        }
        let a = c.assemble();
        let mut factor = SymbolicLdlt::analyze(&a, DirectOrdering::Natural)
            .unwrap()
            .factor(&a, &LdltOptions::default())
            .unwrap();
        assert_eq!(factor.snodes.len(), 1);
        factor.snodes[0].l = vec![1., 0., 0., 1.];
        factor.snodes[0].d = vec![super::super::DBlock::Pair(2., 1., 2.)];
        factor.stats.pivots_1x1 = 0;
        factor.stats.pivots_2x2 = 1;
        assert!(matches!(
            PositiveWhitening::new(&factor, &|| false),
            Err(PositiveFactorError::BlockPivotsUnsupported)
        ));
    }
    #[inline(never)]
    fn portable(s: &PositiveWhitening<'_>, x: &[f64], transpose: bool) -> Vec<f64> {
        let mut out = vec![0.; s.dim()];
        let mut scratch = vec![0.; s.dim()];
        if transpose {
            s.transpose_body(x, &mut out, &mut scratch, &|| false)
                .unwrap();
        } else {
            s.apply_body(x, &mut out, &mut scratch, &|| false).unwrap();
        }
        out
    }

    #[test]
    fn whitening_fma_dispatch_preserves_portable_bits_at_finite_rhs_scales() {
        let a = matrix();
        let factor = SymbolicLdlt::analyze(&a, DirectOrdering::Amd)
            .unwrap()
            .factor(
                &a,
                &LdltOptions {
                    allow_2x2: false,
                    ..LdltOptions::default()
                },
            )
            .unwrap();
        let s = PositiveWhitening::new(&factor, &|| false).unwrap();
        for scale in [1e-100, 1., 1e100] {
            let rhs: Vec<f64> = (0..s.dim()).map(|i| scale * (i as f64 - 3.5)).collect();
            for transpose in [false, true] {
                let mut dispatched = vec![0.; s.dim()];
                let mut scratch = vec![0.; s.dim()];
                if transpose {
                    s.apply_transpose(&rhs, &mut dispatched, &mut scratch, &|| false)
                        .unwrap();
                } else {
                    s.apply(&rhs, &mut dispatched, &mut scratch, &|| false)
                        .unwrap();
                }
                let expected = portable(&s, &rhs, transpose);
                assert_eq!(
                    dispatched.iter().map(|x| x.to_bits()).collect::<Vec<_>>(),
                    expected.iter().map(|x| x.to_bits()).collect::<Vec<_>>()
                );
            }
        }
    }
}
