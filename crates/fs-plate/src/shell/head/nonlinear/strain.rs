//! Complete energy metric on the existing condensed quadratic strain field.
use super::{MembraneReduction, PlateError, bad, dot, mul};

impl MembraneReduction {
    /// Original symmetric modal-pair order: strain is `sum S_ij q_i q_j`.
    /// Diagonal pairs occur once; mixed coefficients already contain both
    /// product-rule contributions. This is geometry, not evolving memory.
    #[must_use]
    pub fn strain_pairs(&self) -> &[(usize, usize)] {
        &self.pairs
    }

    /// Borrow the original area-weighted plane-stress metrics and condensed
    /// engineering-strain pair columns, without copying geometry or history.
    /// A cold consumer can use its existing positive factor/QR owner to change
    /// energy coordinates; the original matrices remain the authority.
    pub fn strain_facets(&self) -> impl ExactSizeIterator<Item = (&[f64; 9], &[[f64; 3]])> {
        self.facets
            .iter()
            .map(|facet| (&facet.a, facet.strains.as_slice()))
    }

    /// Complete Gram `G_ab = sum_facets S_a^T A S_b` of the SAME statically
    /// equilibrated membrane field. For `p_a=q_i*q_j`, stretching energy is
    /// exactly `p^T G p/2`, including every mixed-mode term. Installed tension
    /// and linear bending are excluded. No coefficient or small direction is
    /// discarded. Entries are row-major in [`Self::strain_pairs`] order.
    ///
    /// A uniform scalar proportional Maxwell spectrum may factor this metric
    /// once and retain histories in its range. Because every S_a is already
    /// equilibrated against the in-plane trial space, that common spectrum and
    /// compatible initial history preserve the original static condensation.
    /// Heterogeneous spectra, independently evolving tensor directions, or an
    /// arbitrary uncondensed initial strain do NOT inherit that property.
    ///
    /// `maximum_terms` bounds facets*pair_count^2*12 before allocation. The
    /// caller owns factorization/admission budgets and the constitutive clock.
    pub fn stretching_strain_gram(&self, maximum_terms: usize) -> Result<Vec<f64>, PlateError> {
        let n = self.pairs.len();
        let entries = n
            .checked_mul(n)
            .ok_or_else(|| bad("membrane strain Gram size overflow"))?;
        if entries
            .checked_mul(self.facets.len())
            .and_then(|v| v.checked_mul(12))
            .is_none_or(|v| v > maximum_terms)
        {
            return Err(bad(
                "complete membrane strain Gram exceeds its cold projection budget",
            ));
        }
        let mut gram = vec![0.0; entries];
        for facet in &self.facets {
            for i in 0..n {
                let stress = mul(&facet.a, facet.strains[i]);
                for j in 0..=i {
                    let value = dot(facet.strains[j], stress);
                    gram[i * n + j] += value;
                }
            }
        }
        for i in 0..n {
            for j in 0..i {
                gram[j * n + i] = gram[i * n + j];
            }
        }
        if gram.iter().any(|v| !v.is_finite()) {
            return Err(bad("membrane strain Gram is unrepresentable"));
        }
        Ok(gram)
    }
}

#[cfg(test)]
mod tests {
    use super::super::Facet;
    use super::*;

    #[test]
    fn complete_pair_metric_equals_original_positive_facet_energy() {
        let pairs = vec![(0, 0), (0, 1), (1, 1)];
        let law = MembraneReduction {
            omegas: vec![2.0, 3.0],
            linear: vec![4.0, 0.0, 0.0, 9.0],
            pairs,
            solve_residual: 0.0,
            facets: vec![
                Facet {
                    a: [4.0, 1.0, 0.2, 1.0, 3.0, -0.1, 0.2, -0.1, 2.0],
                    gradients: vec![[0.0; 2]; 2],
                    strains: vec![[1.0, 0.1, 0.0], [0.2, 0.3, 0.4], [-0.2, 0.5, 0.1]],
                },
                Facet {
                    a: [2.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 3.0],
                    gradients: vec![[0.0; 2]; 2],
                    strains: vec![[0.0, 0.3, 0.2], [0.1, -0.4, 0.2], [0.6, 0.1, 0.5]],
                },
            ],
        };
        let gram = law.stretching_strain_gram(216).unwrap();
        assert!(law.stretching_strain_gram(215).is_err());
        for q in [[0.1, 0.2], [-0.7, 0.3], [0.0, 0.9]] {
            let product: Vec<_> = law
                .strain_pairs()
                .iter()
                .map(|&(i, j)| q[i] * q[j])
                .collect();
            let mut energy = 0.0;
            for i in 0..3 {
                for j in 0..3 {
                    assert_eq!(gram[i * 3 + j], gram[j * 3 + i]);
                    energy += 0.5 * product[i] * gram[i * 3 + j] * product[j];
                }
            }
            let original = law.stretching_energy(&q);
            assert!((energy - original).abs() < 2e-14 * original);
        }
    }
}
