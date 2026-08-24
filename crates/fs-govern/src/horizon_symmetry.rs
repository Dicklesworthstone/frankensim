//! Horizon trigger 13b (bead `frankensim-epic-addendum-xpck.5.5`): the
//! real-workload symmetry prevalence and isotypic solver gate for Proposal 13b
//! ("Symmetry harvesting", [`crate::proposals`]).
//!
//! Evaluates representative workload populations for exploitable exact or
//! certified-approximate symmetry. The isotypic solver activates ONLY when
//! at least 15% of the representative workload denominator qualifies AND
//! an independent full-solve falsifier remains available. Below 15%, low-cost
//! detection and abstraction-coarsening signals are preserved while solver
//! promotion is deferred ([`SymmetryDisposition::DetectionOnly`]).

use fs_blake3::hash_bytes;

/// Minimum qualifying prevalence threshold: at least 15% of representative workloads.
pub const SYMMETRY_PREVALENCE_MIN: f64 = 0.15;

/// Maximum admitted asymmetry residual $\epsilon_{\text{asym}}$ for certified approximate symmetry.
pub const MAX_ASYMMETRY_RESIDUAL: f64 = 1e-4;

/// One representative workload evaluated for exploitable symmetry.
#[derive(Debug, Clone, PartialEq)]
pub struct WorkloadSymmetryAssessment {
    /// Workload identifier.
    pub workload_id: String,
    /// Detected symmetry group name (e.g. "C2", "C4v", "D2h", "C1" for asymmetric).
    pub group_name: String,
    /// Group order $|G| \ge 1$ ($|G| \ge 2$ represents non-trivial symmetry).
    pub group_order: usize,
    /// Certified asymmetry residual $\epsilon_{\text{asym}} \ge 0$.
    pub asymmetry_residual: f64,
    /// Whether an independent unreduced full-solve oracle is available to falsify block-diagonalization.
    pub full_solve_falsifier_available: bool,
}

impl WorkloadSymmetryAssessment {
    /// True if this workload qualifies as having exploitable certified symmetry.
    #[must_use]
    pub fn is_qualifying(&self) -> bool {
        self.group_order >= 2
            && self.asymmetry_residual.is_finite()
            && self.asymmetry_residual <= MAX_ASYMMETRY_RESIDUAL
            && self.full_solve_falsifier_available
    }
}

/// Typed refusals for malformed symmetry population assessments.
#[derive(Debug, Clone, PartialEq)]
pub enum Trigger13bRefusal {
    /// Empty workload population.
    EmptyPopulation,
    /// Non-finite or negative asymmetry residual.
    InvalidResidual { workload_id: String, val: f64 },
    /// Group order is zero.
    ZeroGroupOrder { workload_id: String },
}

/// Activation verdict for Proposal 13b.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Trigger13bVerdict {
    /// Prevalence >= 15% with falsifiers: promote and activate isotypic solver.
    ActivateIsotypicSolver,
    /// Prevalence < 15%: preserve low-cost detection, defer solver promotion.
    DetectionOnly,
}

/// Overall population disposition.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SymmetryDisposition {
    /// Promoted and activated.
    Activate,
    /// Evaluated and deferred (detection-only active).
    DetectionOnly,
    /// No representative workload population evaluated yet.
    NoData,
}

/// Immutable receipt of a Trigger 13b evaluation.
#[derive(Debug, Clone, PartialEq)]
pub struct Trigger13bReceipt {
    pub proposal: &'static str,
    pub disposition: SymmetryDisposition,
    pub verdict: Trigger13bVerdict,
    pub qualifying_count: usize,
    pub total_count: usize,
    pub prevalence: f64,
    pub receipt_hash: String,
    pub reason: String,
}

/// Evaluate a population of representative workloads against the 15% symmetry prevalence gate.
///
/// # Errors
/// Returns [`Trigger13bRefusal`] if any workload has malformed quantities.
pub fn evaluate_trigger_13b(
    population: &[WorkloadSymmetryAssessment],
) -> Result<Trigger13bVerdict, Trigger13bRefusal> {
    if population.is_empty() {
        return Err(Trigger13bRefusal::EmptyPopulation);
    }
    for w in population {
        if w.group_order == 0 {
            return Err(Trigger13bRefusal::ZeroGroupOrder { workload_id: w.workload_id.clone() });
        }
        if !w.asymmetry_residual.is_finite() || w.asymmetry_residual < 0.0 {
            return Err(Trigger13bRefusal::InvalidResidual {
                workload_id: w.workload_id.clone(),
                val: w.asymmetry_residual,
            });
        }
    }

    let qualifying = population.iter().filter(|w| w.is_qualifying()).count();
    let prevalence = qualifying as f64 / population.len() as f64;

    if prevalence >= SYMMETRY_PREVALENCE_MIN {
        Ok(Trigger13bVerdict::ActivateIsotypicSolver)
    } else {
        Ok(Trigger13bVerdict::DetectionOnly)
    }
}

/// Mint an immutable decision receipt for Proposal 13b.
#[must_use]
pub fn mint_trigger_13b_receipt(
    population_opt: Option<&[WorkloadSymmetryAssessment]>,
) -> Trigger13bReceipt {
    let Some(population) = population_opt else {
        let hash = hash_bytes(b"org.frankensim.horizon-trigger-13b.nodata.v1").to_hex();
        return Trigger13bReceipt {
            proposal: "13b",
            disposition: SymmetryDisposition::NoData,
            verdict: Trigger13bVerdict::DetectionOnly,
            qualifying_count: 0,
            total_count: 0,
            prevalence: 0.0,
            receipt_hash: hash,
            reason: "no representative real-workload symmetry census exists in the program yet".into(),
        };
    };

    match evaluate_trigger_13b(population) {
        Ok(Trigger13bVerdict::ActivateIsotypicSolver) => {
            let qualifying = population.iter().filter(|w| w.is_qualifying()).count();
            let total = population.len();
            let prevalence = qualifying as f64 / total as f64;
            let mut payload = Vec::new();
            payload.extend_from_slice(b"org.frankensim.horizon-trigger-13b.activate.v1");
            payload.extend_from_slice(prevalence.to_le_bytes().as_slice());
            let hash = hash_bytes(&payload).to_hex();
            Trigger13bReceipt {
                proposal: "13b",
                disposition: SymmetryDisposition::Activate,
                verdict: Trigger13bVerdict::ActivateIsotypicSolver,
                qualifying_count: qualifying,
                total_count: total,
                prevalence,
                receipt_hash: hash,
                reason: format!("symmetry prevalence ({:.1}%, {}/{}) satisfies the 15% activation threshold with full-solve falsifiers", prevalence * 100.0, qualifying, total),
            }
        }
        Ok(Trigger13bVerdict::DetectionOnly) => {
            let qualifying = population.iter().filter(|w| w.is_qualifying()).count();
            let total = population.len();
            let prevalence = qualifying as f64 / total as f64;
            let hash = hash_bytes(b"org.frankensim.horizon-trigger-13b.defer.v1").to_hex();
            Trigger13bReceipt {
                proposal: "13b",
                disposition: SymmetryDisposition::DetectionOnly,
                verdict: Trigger13bVerdict::DetectionOnly,
                qualifying_count: qualifying,
                total_count: total,
                prevalence,
                receipt_hash: hash,
                reason: format!("symmetry prevalence ({:.1}%, {}/{}) sits below the 15% activation threshold; preserving low-cost detection and deferring solver promotion", prevalence * 100.0, qualifying, total),
            }
        }
        Err(refusal) => {
            let hash = hash_bytes(format!("{:?}", refusal).as_bytes()).to_hex();
            Trigger13bReceipt {
                proposal: "13b",
                disposition: SymmetryDisposition::DetectionOnly,
                verdict: Trigger13bVerdict::DetectionOnly,
                qualifying_count: 0,
                total_count: population.len(),
                prevalence: 0.0,
                receipt_hash: hash,
                reason: format!("inadmissible symmetry population ({refusal:?}); preserving detection only"),
            }
        }
    }
}
