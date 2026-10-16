//! Coherent-draw uncertainty (bead wf-root-guzez.5.1.3, E4.0c). Plan
//! §5.2.1: uncertainty is a COHERENT draw of spline coefficients — a
//! low-rank function-space realization — never an independent interval
//! re-sampled per query. That law is STRUCTURAL here: the only way to get
//! an uncertain value is realize-then-query, so two queries of one
//! realization are automatically consistent, and there is no per-query
//! interval API to misuse.
//!
//! Draw weights derive deterministically from a realization id string via
//! `fs-blake3` domain hashing (the id is minted upstream by the identity
//! machinery — ModelUncertaintyRealizationId per replay-identity-schema-v1;
//! this module is a pure function of it).

use crate::Refusal;
use crate::fit::ResidualSurface;

/// Domain separator for realization-weight derivation.
pub const REALIZATION_DOMAIN: &str = "org.frankensim.fs-airfoil.coef-realization.v1";

/// Cap on uncertainty modes (low-rank law: a mode count near the
/// coefficient count is a smell, not a model).
pub const MAX_MODES: usize = 16;

/// A low-rank uncertain surface: mean coefficients plus scaled coefficient
/// modes. A realization is `mean + Σ w_k · mode_k` with weights derived
/// coherently from ONE realization id.
#[derive(Clone, Debug, PartialEq)]
pub struct UncertainSurface {
    /// The mean (fitted) surface.
    pub mean: ResidualSurface,
    /// Coefficient modes, each the same length as `mean.coef`, PRE-SCALED
    /// (a mode's norm IS its one-sigma-equivalent influence).
    pub modes: Vec<Vec<f64>>,
}

impl UncertainSurface {
    /// Validate mode shapes against the mean surface.
    ///
    /// # Errors
    /// `uncertainty-modes-invalid` (count above [`MAX_MODES`], length
    /// mismatch, or non-finite entries).
    pub fn admit(&self) -> Result<(), Refusal> {
        if self.modes.len() > MAX_MODES {
            return Err(Refusal {
                code: "uncertainty-modes-invalid",
                message: format!(
                    "{} modes exceed the low-rank cap {MAX_MODES}",
                    self.modes.len()
                ),
                ranked_repairs: vec![
                    "truncate to the dominant modes; this is a LOW-RANK law".into(),
                ],
            });
        }
        let n = self.mean.coef.len();
        for (k, mode) in self.modes.iter().enumerate() {
            if mode.len() != n {
                return Err(Refusal {
                    code: "uncertainty-modes-invalid",
                    message: format!("mode {k} has {} coefficients, mean has {n}", mode.len()),
                    ranked_repairs: vec![
                        "modes live in the SAME coefficient space as the mean".into(),
                    ],
                });
            }
            if !mode.iter().all(|v| v.is_finite()) {
                return Err(Refusal {
                    code: "uncertainty-modes-invalid",
                    message: format!("mode {k} contains a non-finite entry"),
                    ranked_repairs: vec!["repair the mode construction upstream".into()],
                });
            }
        }
        Ok(())
    }

    /// Deterministic coherent weights in [−1, 1] for `realization_id`:
    /// weight k is derived from the fs-blake3 domain hash of
    /// (REALIZATION_DOMAIN, id, k). A pure function — the same id yields
    /// the same weights forever, on every platform.
    #[must_use]
    pub fn weights(realization_id: &str, n_modes: usize) -> Vec<f64> {
        (0..n_modes)
            .map(|k| {
                let mut payload = realization_id.as_bytes().to_vec();
                payload.extend_from_slice(&(k as u64).to_le_bytes());
                let digest = fs_blake3::hash_domain(REALIZATION_DOMAIN, &payload);
                let word = u64::from_le_bytes(digest.as_bytes()[0..8].try_into().unwrap());
                // Map to [-1, 1] with 53-bit uniformity.
                ((word >> 11) as f64) / ((1u64 << 52) as f64) - 1.0
            })
            .collect()
    }

    /// Realize ONE coherent surface for `realization_id`. Every subsequent
    /// query of the returned surface shares the same draw — per-query
    /// independent intervals are impossible by construction.
    ///
    /// # Errors
    /// Mode-shape refusals from [`Self::admit`];
    /// `realization-id-empty` (an anonymous draw is not a realization).
    pub fn realize(&self, realization_id: &str) -> Result<ResidualSurface, Refusal> {
        if realization_id.is_empty() {
            return Err(Refusal {
                code: "realization-id-empty",
                message: "a coherent draw must bind a realization id".into(),
                ranked_repairs: vec![
                    "pass the ModelUncertaintyRealizationId (replay-identity-schema-v1)".into(),
                ],
            });
        }
        self.admit()?;
        let w = Self::weights(realization_id, self.modes.len());
        let mut surface = self.mean.clone();
        for (weight, mode) in w.iter().zip(&self.modes) {
            for (c, m) in surface.coef.iter_mut().zip(mode) {
                *c += weight * m;
            }
        }
        Ok(surface)
    }
}
