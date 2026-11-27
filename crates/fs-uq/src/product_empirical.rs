//! Content-addressed empirical marginals for Gaussian-copula propagation.

use std::sync::Arc;

use fs_blake3::{ContentHash, DomainHasher};

/// Identifier of the right-closed inverse-ECDF step rule.
pub const INVERSE_ECDF_RULE_ID: &str = "inverse-ecdf-step-v1";

/// Maximum number of observations in one empirical marginal (8 MiB of binary64).
pub const MAX_EMPIRICAL_SAMPLES: usize = 1 << 20;

/// An equally weighted finite sample. Repetitions retain their empirical mass.
/// Physical units belong to the enclosing parameter declaration.
#[derive(Clone, Debug, PartialEq)]
pub struct EmpiricalMarginal {
    samples: Arc<[f64]>,
    sample_hash: ContentHash,
}

impl EmpiricalMarginal {
    /// Admit a nonempty sample in canonical ascending binary64 total order.
    /// All observations must be finite, with negative zero before positive zero.
    /// The declared quantile rule must equal [`INVERSE_ECDF_RULE_ID`].
    pub fn new(samples: Vec<f64>, rule_id: &str) -> Result<Self, &'static str> {
        if rule_id != INVERSE_ECDF_RULE_ID {
            return Err("unknown empirical quantile rule");
        }
        if samples.is_empty() || samples.len() > MAX_EMPIRICAL_SAMPLES {
            return Err("empirical sample length must be in 1..=1048576");
        }
        if !samples.iter().all(|x| x.is_finite()) {
            return Err("empirical observations must be finite");
        }
        if !samples.is_sorted_by(|a, b| a.total_cmp(b).is_le()) {
            return Err("empirical observations must be in canonical ascending order");
        }

        let mut hash = DomainHasher::new("org.frankensim.uq.empirical-sample.v1");
        for sample in &samples {
            hash.update(&sample.to_le_bytes());
        }
        Ok(Self {
            samples: samples.into(),
            sample_hash: hash.finalize(),
        })
    }

    /// Canonically ordered observations, including ties.
    #[must_use]
    pub fn samples(&self) -> &[f64] {
        &self.samples
    }

    /// Number of source observations, including ties.
    #[must_use]
    pub fn sample_count(&self) -> usize {
        self.samples.len()
    }

    /// Observed minimum and maximum. The empirical law represents no tails
    /// beyond these endpoints and assigns mass only to observed values.
    #[must_use]
    pub fn observed_support(&self) -> [f64; 2] {
        [self.samples[0], self.samples[self.samples.len() - 1]]
    }

    /// BLAKE3 domain hash of the ordered, little-endian binary64 sample bytes.
    #[must_use]
    pub const fn sample_hash(&self) -> ContentHash {
        self.sample_hash
    }

    /// Declared inverse-CDF rule bound into a copula checkpoint.
    #[must_use]
    pub const fn rule_id(&self) -> &'static str {
        INVERSE_ECDF_RULE_ID
    }

    /// Order statistic at one-based rank max(1, ceil(n * u)) for admitted u.
    pub(crate) fn quantile(&self, u: f64) -> f64 {
        debug_assert!((0.0..=1.0).contains(&u));
        let n = self.sample_count() as f64;
        let whole = (u * n) as usize;
        // The rounded product can be exactly k even when u is just above
        // k/n. A fused subtraction retains the sign of the exact remainder.
        let above = u.mul_add(n, -(whole as f64)) > 0.0;
        let rank = whole + usize::from(above);
        self.samples[rank.saturating_sub(1)]
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn step_rule_includes_the_upper_boundary_and_the_observed_endpoints() {
        let law = EmpiricalMarginal::new(vec![1.0, 4.0, 4.0, 9.0], INVERSE_ECDF_RULE_ID).unwrap();
        for (u, value) in [
            (0.0, 1.0_f64),
            (0.25_f64.next_down(), 1.0),
            (0.25, 1.0),
            (0.25_f64.next_up(), 4.0),
            (0.5, 4.0),
            (0.75, 4.0),
            (0.75_f64.next_up(), 9.0),
            (1.0, 9.0),
        ] {
            assert_eq!(law.quantile(u).to_bits(), value.to_bits());
        }
        let three = EmpiricalMarginal::new(vec![1.0, 2.0, 3.0], INVERSE_ECDF_RULE_ID).unwrap();
        let below = 1.0_f64 / 3.0;
        let above = below.next_up();
        assert_eq!((above * 3.0).to_bits(), 1.0_f64.to_bits());
        assert_eq!(three.quantile(below).to_bits(), 1.0_f64.to_bits());
        assert_eq!(three.quantile(above).to_bits(), 2.0_f64.to_bits());

        let one = EmpiricalMarginal::new(vec![-0.0], INVERSE_ECDF_RULE_ID).unwrap();
        for u in [0.0, 0.5, 1.0] {
            assert_eq!(one.quantile(u).to_bits(), (-0.0_f64).to_bits());
        }
    }

    #[test]
    fn sample_identity_preserves_ties_and_exact_bits() {
        let samples = vec![-0.0_f64, 0.0, 5.0, 5.0];
        let bytes: Vec<_> = samples.iter().flat_map(|x| x.to_le_bytes()).collect();
        let law = EmpiricalMarginal::new(samples, INVERSE_ECDF_RULE_ID).unwrap();
        assert_eq!(
            law.sample_hash(),
            fs_blake3::hash_domain("org.frankensim.uq.empirical-sample.v1", &bytes)
        );
        assert_eq!(law.sample_count(), 4);
        assert_eq!(
            law.observed_support().map(f64::to_bits),
            [(-0.0_f64).to_bits(), 5.0_f64.to_bits()]
        );
        let changed =
            EmpiricalMarginal::new(vec![0.0, 0.0, 5.0, 5.0], INVERSE_ECDF_RULE_ID).unwrap();
        assert_ne!(law.sample_hash(), changed.sample_hash());
    }

    #[test]
    fn refuses_unknown_rule_empty_unordered_nonfinite_and_oversized_samples() {
        assert!(EmpiricalMarginal::new(vec![0.0], "interpolate").is_err());
        for samples in [
            vec![],
            vec![2.0, 1.0],
            vec![0.0, -0.0],
            vec![f64::NAN],
            vec![f64::NEG_INFINITY],
            vec![f64::INFINITY],
            vec![0.0; MAX_EMPIRICAL_SAMPLES + 1],
        ] {
            assert!(EmpiricalMarginal::new(samples, INVERSE_ECDF_RULE_ID).is_err());
        }
    }
}
