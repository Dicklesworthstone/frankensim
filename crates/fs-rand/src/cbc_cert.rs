//! Per-prefix minimality certificates for exact component-by-component
//! lattice construction, with an independent checker (bead 6ys.20,
//! certificate tranche).
//!
//! A [`CbcPrefixCertificate`] binds, for one scanned component: the exact
//! reduction root (the committed prefix itself — the products a scan scored
//! against are a pure function of `(n, prefix)`, so declaring the prefix
//! declares them exactly, with no hash indirection), the chosen candidate,
//! the exact winning numerator, the exact runner-up numerator or its
//! absence, the winning equality class, the tie rule, the admissible-set
//! rule, and the common-denominator derivation the numerators drop.
//!
//! The checker is INDEPENDENT — it recomputes with its own arithmetic from
//! the declared inputs and never trusts executor state — and has two
//! honestly named modes:
//!
//! - [`verify_consistency`] recomputes the winning, tie-class, and
//!   runner-up scores from the declared prefix in `O(n · |claims|)` exact
//!   work: cheap, but it proves the DECLARED candidates score as claimed,
//!   not that no other candidate scores lower.
//! - [`audit_minimality`] rescans every admissible candidate at the
//!   declared prefix (`O(n²)` exact work per certificate): the full
//!   minimality proof by exhaustion.
//!
//! NO-CLAIM: a compact sub-quadratic minimality proof (branch-and-bound /
//! sheaf-glued sections) is the bead's [M] ratchet, not this tranche; the
//! first component's unit-residue-permutation theorem certificate is the
//! [F] ratchet, so certificates here start at the first SCANNED component.
//! Certificates are plain data with no fs-blake3 identity minting; identity
//! governance for durable certificate stores belongs to consumers.

use crate::qmc::{ExactNat, exact_kernel_numerator, gcd, lattice_residue};

/// Version of the certificate schema and checker semantics.
pub const CBC_CERTIFICATE_SCHEMA_VERSION: u32 = 1;

/// The declared tie rule token (the only rule this schema admits).
pub const TIE_RULE_LOWEST_CANDIDATE: &str = "lowest-candidate-wins";

/// The declared admissible-set rule token (units modulo `n`).
pub const ADMISSIBLE_RULE_UNITS: &str = "units-modulo-n";

/// One scanned component's exact selection evidence. Construction is
/// executor-side ([`crate::cbc_exec::CbcExecutor`]); checking is here.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CbcPrefixCertificate {
    /// Number of lattice points `n`.
    pub point_count: u32,
    /// The committed prefix INCLUDING the certified component (which is the
    /// last element). This is the exact reduction root: the scan's products
    /// are a pure function of `(n, prefix[..len-1])`.
    pub prefix: Vec<u32>,
    /// Exact winning-score numerator, normalized little-endian base-2³²
    /// limbs, over the dropped denominator `(6n²)^(prefix.len())`.
    pub winning_score_limbs: Vec<u32>,
    /// Candidates whose exact score equals the winning score, ascending.
    /// Always contains the chosen candidate as its minimum. Mirror symmetry
    /// (`c` and `n−c` share a residue multiset) makes real ties structural,
    /// not exotic.
    pub tie_class: Vec<u32>,
    /// The smallest exact score strictly above the winning score, with its
    /// lowest achieving candidate; `None` when every admissible candidate
    /// ties the winner.
    pub runner_up: Option<(Vec<u32>, u32)>,
    /// Exponent `e` in the dropped common denominator `(6n²)^e` (the number
    /// of kernel factors in every compared numerator at this step).
    pub denominator_exponent: u32,
    /// Tie rule token; must equal [`TIE_RULE_LOWEST_CANDIDATE`].
    pub tie_rule: &'static str,
    /// Admissible-set rule token; must equal [`ADMISSIBLE_RULE_UNITS`].
    pub admissible_rule: &'static str,
}

impl CbcPrefixCertificate {
    /// The certified (chosen) component.
    #[must_use]
    pub fn chosen(&self) -> u32 {
        *self
            .prefix
            .last()
            .expect("a certificate binds a non-empty prefix")
    }
}

/// Checker refusals: every variant names what disagreed. Fail-closed — any
/// tampering with scores, limbs, ties, candidates, or derivation tokens
/// lands in exactly one of these.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum CbcCertError {
    /// The prefix was empty or the certified component is the theorem-fixed
    /// first component (certificates start at the first scanned component).
    MalformedPrefix,
    /// `point_count < 3`, a component out of `1..n`, or a non-coprime
    /// prefix component.
    InadmissiblePrefix,
    /// A rule token differs from the schema's declared rule.
    UnknownRule,
    /// The denominator exponent does not equal the prefix length.
    DenominatorMismatch,
    /// The tie class is empty, unsorted, out of range, non-coprime, or its
    /// minimum is not the chosen candidate.
    MalformedTieClass,
    /// A declared tie-class member's recomputed score differs from the
    /// declared winning score.
    TieClassScoreMismatch {
        /// The disagreeing candidate.
        candidate: u32,
    },
    /// The recomputed runner-up score or candidate differs from the
    /// declaration, or the declared runner-up does not score strictly above
    /// the winner.
    RunnerUpMismatch,
    /// Full audit: some admissible candidate scores strictly below the
    /// declared winning score.
    NotMinimal {
        /// The candidate that beats the declared winner.
        candidate: u32,
    },
    /// Full audit: the true equality class differs from the declaration.
    TieClassIncomplete,
}

/// Recompute the prefix products (the reduction root) from `(n, prefix)`.
fn products_for(n: u32, prefix: &[u32]) -> Vec<ExactNat> {
    let point_count = usize::try_from(n).expect("checker point count fits usize");
    let mut products = vec![ExactNat::one(); point_count];
    for &component in prefix {
        for (point, product) in products.iter_mut().enumerate() {
            let residue = lattice_residue(point, component, n);
            product.mul_assign_factor(exact_kernel_numerator(n, residue));
        }
    }
    products
}

/// Exact score of `candidate` against `products` (normalized).
fn score_for(n: u32, products: &[ExactNat], candidate: u32) -> ExactNat {
    let mut score = ExactNat::zero();
    for (point, product) in products.iter().enumerate() {
        let residue = lattice_residue(point, candidate, n);
        score.add_mul_factor(product, exact_kernel_numerator(n, residue));
    }
    score.normalize();
    score
}

/// Structural checks shared by both modes. Returns the scan prefix (the
/// declared prefix without its certified last component).
fn structural<'a>(certificate: &'a CbcPrefixCertificate) -> Result<&'a [u32], CbcCertError> {
    let n = certificate.point_count;
    if n < 3 {
        return Err(CbcCertError::InadmissiblePrefix);
    }
    if certificate.prefix.len() < 2 {
        return Err(CbcCertError::MalformedPrefix);
    }
    if certificate.tie_rule != TIE_RULE_LOWEST_CANDIDATE
        || certificate.admissible_rule != ADMISSIBLE_RULE_UNITS
    {
        return Err(CbcCertError::UnknownRule);
    }
    let expected_exponent =
        u32::try_from(certificate.prefix.len()).map_err(|_| CbcCertError::MalformedPrefix)?;
    if certificate.denominator_exponent != expected_exponent {
        return Err(CbcCertError::DenominatorMismatch);
    }
    for &component in &certificate.prefix {
        if component == 0 || component >= n || gcd(component, n) != 1 {
            return Err(CbcCertError::InadmissiblePrefix);
        }
    }
    let chosen = certificate.chosen();
    if certificate.tie_class.is_empty()
        || certificate.tie_class.first() != Some(&chosen)
        || !certificate
            .tie_class
            .windows(2)
            .all(|pair| pair[0] < pair[1])
        || certificate
            .tie_class
            .iter()
            .any(|&candidate| candidate == 0 || candidate >= n || gcd(candidate, n) != 1)
    {
        return Err(CbcCertError::MalformedTieClass);
    }
    if let Some((_, runner)) = &certificate.runner_up {
        if *runner == 0 || *runner >= n || gcd(*runner, n) != 1 {
            return Err(CbcCertError::RunnerUpMismatch);
        }
        if certificate.tie_class.contains(runner) {
            return Err(CbcCertError::RunnerUpMismatch);
        }
    }
    Ok(&certificate.prefix[..certificate.prefix.len() - 1])
}

/// Compact consistency check: recompute the declared candidates' exact
/// scores from the declared reduction root. Proves the declaration is
/// internally exact; does NOT prove global minimality (see
/// [`audit_minimality`]).
///
/// # Errors
/// The first [`CbcCertError`] the certificate fails.
pub fn verify_consistency(certificate: &CbcPrefixCertificate) -> Result<(), CbcCertError> {
    let scan_prefix = structural(certificate)?;
    let n = certificate.point_count;
    let products = products_for(n, scan_prefix);
    for &candidate in &certificate.tie_class {
        let score = score_for(n, &products, candidate);
        if score.limbs() != certificate.winning_score_limbs.as_slice() {
            return Err(CbcCertError::TieClassScoreMismatch { candidate });
        }
    }
    if let Some((declared_limbs, runner)) = &certificate.runner_up {
        let score = score_for(n, &products, *runner);
        if score.limbs() != declared_limbs.as_slice() {
            return Err(CbcCertError::RunnerUpMismatch);
        }
        let mut winning = ExactNat::zero();
        winning.add_mul_factor(&exact_from_limbs(&certificate.winning_score_limbs), 1);
        winning.normalize();
        if score.magnitude_cmp(&winning) != core::cmp::Ordering::Greater {
            return Err(CbcCertError::RunnerUpMismatch);
        }
    }
    Ok(())
}

/// Full minimality audit by exhaustion: rescan every admissible candidate
/// at the declared reduction root and require the declared winner, equality
/// class, and runner-up to be exactly what the scan finds.
///
/// # Errors
/// The first [`CbcCertError`] the certificate fails.
pub fn audit_minimality(certificate: &CbcPrefixCertificate) -> Result<(), CbcCertError> {
    let scan_prefix = structural(certificate)?;
    let n = certificate.point_count;
    let products = products_for(n, scan_prefix);
    let winning = exact_from_limbs(&certificate.winning_score_limbs);
    let mut true_ties = Vec::new();
    let mut true_runner: Option<(ExactNat, u32)> = None;
    for candidate in 1..n {
        if gcd(candidate, n) != 1 {
            continue;
        }
        let score = score_for(n, &products, candidate);
        match score.magnitude_cmp(&winning) {
            core::cmp::Ordering::Less => {
                return Err(CbcCertError::NotMinimal { candidate });
            }
            core::cmp::Ordering::Equal => true_ties.push(candidate),
            core::cmp::Ordering::Greater => {
                let replace = match &true_runner {
                    None => true,
                    Some((best_above, _)) => {
                        score.magnitude_cmp(best_above) == core::cmp::Ordering::Less
                    }
                };
                if replace {
                    true_runner = Some((score, candidate));
                }
            }
        }
    }
    if true_ties != certificate.tie_class {
        return Err(CbcCertError::TieClassIncomplete);
    }
    match (&certificate.runner_up, true_runner) {
        (None, None) => Ok(()),
        (Some((declared_limbs, declared_candidate)), Some((true_score, true_candidate))) => {
            if true_score.limbs() == declared_limbs.as_slice()
                && *declared_candidate == true_candidate
            {
                Ok(())
            } else {
                Err(CbcCertError::RunnerUpMismatch)
            }
        }
        _ => Err(CbcCertError::RunnerUpMismatch),
    }
}

/// Rebuild an [`ExactNat`] from declared limbs (normalizing defensively).
fn exact_from_limbs(limbs: &[u32]) -> ExactNat {
    let mut value = ExactNat::zero();
    for (index, &limb) in limbs.iter().enumerate() {
        if limb == 0 {
            continue;
        }
        let mut term = ExactNat::one();
        // Shift the term to limb position `index` by repeated base moves,
        // then scale by the limb. Certificate limbs are short (they grow
        // with prefix length only), so this stays cheap and dependency-free.
        for _ in 0..index {
            term.mul_assign_factor(1_u128 << 32);
        }
        term.mul_assign_factor(u128::from(limb));
        value.add_mul_factor(&term, 1);
    }
    value.normalize();
    value
}
