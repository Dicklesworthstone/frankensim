//! fs-grammar-e2e — GrammarForge: certified-fabricable geometric program
//! discovery. Layer: L4 (ASCENT).
//!
//! # The campaign
//!
//! A CAD model is one hand-built artifact with no guarantees. This instead
//! ILLUMINATES the diverse family of CSG PROGRAMS that approximate a target
//! shape and are fabricable, composing crates never designed to meet:
//!
//! - **Programs as data** ([`fs_shapeprog`]): a candidate is a CSG program —
//!   `sphere(r₁)@−d ∪ (sphere(r₂)⊕o)@+d`. Its fidelity is the worst-case SDF
//!   discrepancy from the target over a sample grid.
//! - **Certificate-preserving simplification** ([`fs_shapeprog::simplify`]): the
//!   rewrite engine drops redundant tiny offsets and applies geometric
//!   identities, each with a fidelity certificate (`Exact` or `Approximate{bound}`),
//!   so the simplified program is provably within `max_error` of the original —
//!   and the campaign INDEPENDENTLY re-measures the discrepancy to confirm the
//!   certificate holds (certifying the certifier). The local offset-radius
//!   admission threshold and returned global envelope remain distinct.
//! - **Manufacturability** ([`fs_fab`]): a minimum-feature-size constraint scores
//!   each program's smallest feature — the fabrication margin.
//! - **Illumination** ([`fs_archive`]): MAP-Elites over (program size × fab
//!   margin) keeps the best-matching program in every complexity/fabricability
//!   niche — the diverse atlas, not one model.
//! - **Honest colors** ([`fs_evidence`]): a program that matches within tolerance,
//!   is fab-satisfied, and simplifies soundly is `Verified`.
//!
//! Deterministic; no dependencies beyond the composed crates.

use fs_archive::MapElites;
use fs_evidence::Color;
use fs_fab::min_feature_size;
use fs_shapeprog::{Geom, SimplifyRefusal, max_sdf_discrepancy, simplify};

/// The target shape: a "peanut" — two unit spheres at `x = ±0.8`.
#[must_use]
pub fn target() -> Geom {
    Geom::sphere(1.0)
        .translate([-0.8, 0.0, 0.0])
        .union(Geom::sphere(1.0).translate([0.8, 0.0, 0.0]))
}

/// Build a candidate program from parameters `[r1, r2, d, o]`.
#[must_use]
pub fn build_program(r1: f64, r2: f64, d: f64, o: f64) -> Geom {
    let left = Geom::sphere(r1).translate([-d, 0.0, 0.0]);
    let right = Geom::sphere(r2).offset(o).translate([d, 0.0, 0.0]);
    left.union(right)
}

/// A deterministic 3-D sample grid over `[-2, 2]³` for SDF discrepancy.
#[must_use]
fn sample_points() -> Vec<[f64; 3]> {
    let mut pts = Vec::new();
    let n = 7;
    for i in 0..n {
        for j in 0..n {
            for k in 0..n {
                let c = |t: usize| -2.0 + 4.0 * t as f64 / (n - 1) as f64;
                pts.push([c(i), c(j), c(k)]);
            }
        }
    }
    pts
}

/// Typed outcome of independently checking one ShapeProg simplification or
/// aggregating multiple checks.
///
/// Stable wire codes are consumed by `fs-wasm`; do not reorder variants or
/// change their discriminants without migrating that wire contract. Codes are
/// identifiers, not severity ranks.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[repr(u8)]
pub enum SimplificationCheckStatus {
    /// A finite conservative outward sample check is no larger than the
    /// finite certificate.
    Certified = 0,
    /// Both programs have the core-recognized structural-empty `+∞` semantics.
    StructuralEmptyAgreement = 1,
    /// The simplifier transactionally refused and returned its original input.
    SimplifierRefused = 2,
    /// A nominally successful simplification exposed a non-finite certificate.
    /// This is defensive consumer validation of the core contract.
    NonFiniteCertificate = 3,
    /// A nominally successful simplification exposed a finite negative
    /// certificate, violating the nonnegative-envelope contract.
    NegativeCertificate = 4,
    /// Independent discrepancy evidence was empty, non-finite, or otherwise
    /// refused by the fail-closed ShapeProg checker.
    DiscrepancyEvidenceRefused = 5,
    /// A conservative outward finite-sample check exceeded the certificate.
    CertificateCheckExceeded = 6,
    /// An assessment from a different local radius threshold was mixed into a
    /// campaign summary.
    ThresholdMismatch = 7,
}

impl SimplificationCheckStatus {
    /// Stable integer representation for report/wire serialization.
    #[must_use]
    pub const fn wire_code(self) -> u8 {
        self as u8
    }

    /// Whether this individual/aggregate status authorizes the sampled
    /// simplification-soundness claim.
    #[must_use]
    pub const fn is_sound(self) -> bool {
        matches!(self, Self::Certified | Self::StructuralEmptyAgreement)
    }

    const fn severity(self) -> u8 {
        match self {
            Self::Certified | Self::StructuralEmptyAgreement => 0,
            Self::SimplifierRefused => 1,
            Self::NegativeCertificate => 2,
            Self::NonFiniteCertificate => 3,
            Self::DiscrepancyEvidenceRefused => 4,
            Self::ThresholdMismatch => 5,
            Self::CertificateCheckExceeded => 6,
        }
    }

    fn combine(self, next: Self) -> Self {
        if self.is_sound() && next.is_sound() {
            if matches!(self, Self::StructuralEmptyAgreement)
                && matches!(next, Self::StructuralEmptyAgreement)
            {
                Self::StructuralEmptyAgreement
            } else {
                Self::Certified
            }
        } else if next.severity() > self.severity() {
            next
        } else {
            self
        }
    }
}

/// One sealed simplification/certifier result record. Equality is exact over
/// IEEE-754 bits, including signed zero and NaN payloads in refused inputs.
#[derive(Debug, Clone)]
#[non_exhaustive]
pub struct SimplificationAssessment {
    /// Strict local `|offset radius|` admission threshold supplied to ShapeProg.
    radius_threshold: f64,
    /// Program nodes before the simplification transaction.
    size_before: usize,
    /// Program nodes after success, or the unchanged size after refusal.
    size_after: usize,
    /// Finite compositional uniform error envelope, when the core produced one.
    certified_error: Option<f64>,
    /// Exact bits of a rejected non-finite or negative nominal certificate.
    rejected_certificate_bits: Option<u64>,
    /// Finite outward sampled discrepancy witness, when evidence was admitted.
    sampled_discrepancy: Option<f64>,
    /// Typed consumer-side interpretation of the check.
    status: SimplificationCheckStatus,
    /// Exact structured ShapeProg refusal, when status is `SimplifierRefused`.
    refusal: Option<SimplifyRefusal>,
}

impl PartialEq for SimplificationAssessment {
    fn eq(&self, other: &Self) -> bool {
        self.radius_threshold.to_bits() == other.radius_threshold.to_bits()
            && self.size_before == other.size_before
            && self.size_after == other.size_after
            && self.certified_error.map(f64::to_bits) == other.certified_error.map(f64::to_bits)
            && self.rejected_certificate_bits == other.rejected_certificate_bits
            && self.sampled_discrepancy.map(f64::to_bits)
                == other.sampled_discrepancy.map(f64::to_bits)
            && self.status == other.status
            && self.refusal == other.refusal
    }
}

impl Eq for SimplificationAssessment {}

impl SimplificationAssessment {
    /// Exact local offset-radius admission threshold.
    #[must_use]
    pub const fn radius_threshold(&self) -> f64 {
        self.radius_threshold
    }

    /// Program size before the transaction.
    #[must_use]
    pub const fn size_before(&self) -> usize {
        self.size_before
    }

    /// Program size after success or rollback.
    #[must_use]
    pub const fn size_after(&self) -> usize {
        self.size_after
    }

    /// Finite compositional certificate, if one was admitted.
    #[must_use]
    pub const fn certified_error(&self) -> Option<f64> {
        self.certified_error
    }

    /// Exact bits of a rejected invalid nominal certificate, if any.
    #[must_use]
    pub const fn rejected_certificate_bits(&self) -> Option<u64> {
        self.rejected_certificate_bits
    }

    /// Finite outward sampled discrepancy, if evidence was admitted.
    #[must_use]
    pub const fn sampled_discrepancy(&self) -> Option<f64> {
        self.sampled_discrepancy
    }

    /// Typed certifier status.
    #[must_use]
    pub const fn status(&self) -> SimplificationCheckStatus {
        self.status
    }

    /// Exact structured ShapeProg refusal, if any.
    #[must_use]
    pub fn refusal(&self) -> Option<&SimplifyRefusal> {
        self.refusal.as_ref()
    }
}

/// Simplify one program and independently check the returned certificate.
///
/// The local radius threshold and global certified discrepancy are deliberately
/// separate fields. The latter can exceed the former because a context-free
/// dropped-offset envelope is `2*|radius|`, sequential losses add, and retained
/// rounded parents can contribute additional envelopes.
#[must_use]
pub fn assess_simplification(
    program: &Geom,
    radius_threshold: f64,
    samples: &[[f64; 3]],
) -> SimplificationAssessment {
    let size_before = program.size();
    let simplified = simplify(program, radius_threshold);
    let size_after = simplified.program.size();

    if let Some(refusal) = simplified.refusal.clone() {
        return SimplificationAssessment {
            radius_threshold,
            size_before,
            size_after,
            certified_error: None,
            rejected_certificate_bits: None,
            sampled_discrepancy: None,
            status: SimplificationCheckStatus::SimplifierRefused,
            refusal: Some(refusal),
        };
    }

    let bound = simplified.max_error;
    if !bound.is_finite() {
        return SimplificationAssessment {
            radius_threshold,
            size_before,
            size_after,
            certified_error: None,
            rejected_certificate_bits: Some(bound.to_bits()),
            sampled_discrepancy: None,
            status: SimplificationCheckStatus::NonFiniteCertificate,
            refusal: None,
        };
    }
    if bound < 0.0 {
        return SimplificationAssessment {
            radius_threshold,
            size_before,
            size_after,
            certified_error: None,
            rejected_certificate_bits: Some(bound.to_bits()),
            sampled_discrepancy: None,
            status: SimplificationCheckStatus::NegativeCertificate,
            refusal: None,
        };
    }

    let actual = max_sdf_discrepancy(program, &simplified.program, samples);
    if !actual.is_finite() {
        return SimplificationAssessment {
            radius_threshold,
            size_before,
            size_after,
            certified_error: Some(bound),
            rejected_certificate_bits: None,
            sampled_discrepancy: None,
            status: SimplificationCheckStatus::DiscrepancyEvidenceRefused,
            refusal: None,
        };
    }

    let structural_empty_agreement = actual == 0.0
        && samples.iter().all(|&point| {
            program.sdf(point) == f64::INFINITY && simplified.program.sdf(point) == f64::INFINITY
        });
    let status = if structural_empty_agreement {
        SimplificationCheckStatus::StructuralEmptyAgreement
    } else if actual <= bound {
        SimplificationCheckStatus::Certified
    } else {
        SimplificationCheckStatus::CertificateCheckExceeded
    };
    SimplificationAssessment {
        radius_threshold,
        size_before,
        size_after,
        certified_error: Some(bound),
        rejected_certificate_bits: None,
        sampled_discrepancy: Some(actual),
        status,
        refusal: None,
    }
}

/// Aggregate certificate accounting shared by native and WASM GrammarForge.
/// Equality uses the same exact IEEE-754 identity as threshold aggregation.
#[derive(Debug, Clone, Copy)]
#[non_exhaustive]
pub struct SimplificationSummary {
    /// Strict local offset-radius threshold used by every assessment.
    radius_threshold: f64,
    /// Number of incorporated assessment records.
    assessments: usize,
    /// Assessments whose returned program contains fewer AST nodes. This
    /// syntactic count does not itself authorize a certificate claim.
    simplified_count: usize,
    /// Total input AST nodes.
    size_before: usize,
    /// Total successful/rollback output AST nodes.
    size_after: usize,
    /// Largest finite compositional certificate among all assessments.
    max_certified_error: f64,
    /// Largest admitted conservative outward finite-sample check.
    max_sampled_discrepancy: f64,
    /// Aggregate typed status. An all-structural aggregate preserves
    /// `StructuralEmptyAgreement`; mixed finite/structural success is
    /// `Certified`, with the separate structural count visible below.
    status: SimplificationCheckStatus,
    /// Transactional simplifier refusals.
    simplifier_refusals: usize,
    /// Non-finite nominal certificates.
    non_finite_certificates: usize,
    /// Invalid finite negative nominal certificates.
    negative_certificates: usize,
    /// Independent discrepancy-evidence refusals.
    discrepancy_evidence_refusals: usize,
    /// Admitted structural-empty agreements.
    structural_empty_agreements: usize,
    /// Conservative outward finite-sample checks exceeding their certificate.
    certificate_check_exceedances: usize,
    /// Records created under another radius threshold.
    threshold_mismatches: usize,
}

impl PartialEq for SimplificationSummary {
    fn eq(&self, other: &Self) -> bool {
        self.radius_threshold.to_bits() == other.radius_threshold.to_bits()
            && self.assessments == other.assessments
            && self.simplified_count == other.simplified_count
            && self.size_before == other.size_before
            && self.size_after == other.size_after
            && self.max_certified_error.to_bits() == other.max_certified_error.to_bits()
            && self.max_sampled_discrepancy.to_bits() == other.max_sampled_discrepancy.to_bits()
            && self.status == other.status
            && self.simplifier_refusals == other.simplifier_refusals
            && self.non_finite_certificates == other.non_finite_certificates
            && self.negative_certificates == other.negative_certificates
            && self.discrepancy_evidence_refusals == other.discrepancy_evidence_refusals
            && self.structural_empty_agreements == other.structural_empty_agreements
            && self.certificate_check_exceedances == other.certificate_check_exceedances
            && self.threshold_mismatches == other.threshold_mismatches
    }
}

impl Eq for SimplificationSummary {}

impl SimplificationSummary {
    /// Start empty accounting for one exact threshold identity.
    #[must_use]
    pub const fn new(radius_threshold: f64) -> Self {
        Self {
            radius_threshold,
            assessments: 0,
            simplified_count: 0,
            size_before: 0,
            size_after: 0,
            max_certified_error: 0.0,
            max_sampled_discrepancy: 0.0,
            status: SimplificationCheckStatus::Certified,
            simplifier_refusals: 0,
            non_finite_certificates: 0,
            negative_certificates: 0,
            discrepancy_evidence_refusals: 0,
            structural_empty_agreements: 0,
            certificate_check_exceedances: 0,
            threshold_mismatches: 0,
        }
    }

    /// Exact local offset-radius admission threshold shared by all records.
    #[must_use]
    pub const fn radius_threshold(&self) -> f64 {
        self.radius_threshold
    }

    /// Number of incorporated assessment records.
    #[must_use]
    pub const fn assessments(&self) -> usize {
        self.assessments
    }

    /// Number of returned programs that reduced AST size, independent of status.
    #[must_use]
    pub const fn simplified_count(&self) -> usize {
        self.simplified_count
    }

    /// Total input AST nodes.
    #[must_use]
    pub const fn size_before(&self) -> usize {
        self.size_before
    }

    /// Total successful/rollback output AST nodes.
    #[must_use]
    pub const fn size_after(&self) -> usize {
        self.size_after
    }

    /// Largest finite compositional certificate among all records.
    #[must_use]
    pub const fn max_certified_error(&self) -> f64 {
        self.max_certified_error
    }

    /// Largest admitted outward finite-sample discrepancy check.
    #[must_use]
    pub const fn max_sampled_discrepancy(&self) -> f64 {
        self.max_sampled_discrepancy
    }

    /// Aggregate typed status.
    #[must_use]
    pub const fn status(&self) -> SimplificationCheckStatus {
        self.status
    }

    /// Number of transactional simplifier refusals.
    #[must_use]
    pub const fn simplifier_refusals(&self) -> usize {
        self.simplifier_refusals
    }

    /// Number of non-finite nominal certificates.
    #[must_use]
    pub const fn non_finite_certificates(&self) -> usize {
        self.non_finite_certificates
    }

    /// Number of invalid finite negative nominal certificates.
    #[must_use]
    pub const fn negative_certificates(&self) -> usize {
        self.negative_certificates
    }

    /// Number of independent discrepancy-evidence refusals.
    #[must_use]
    pub const fn discrepancy_evidence_refusals(&self) -> usize {
        self.discrepancy_evidence_refusals
    }

    /// Number of admitted structural-empty agreements.
    #[must_use]
    pub const fn structural_empty_agreements(&self) -> usize {
        self.structural_empty_agreements
    }

    /// Number of conservative outward finite-sample checks exceeding their
    /// certificate.
    #[must_use]
    pub const fn certificate_check_exceedances(&self) -> usize {
        self.certificate_check_exceedances
    }

    /// Number of records created under another radius threshold.
    #[must_use]
    pub const fn threshold_mismatches(&self) -> usize {
        self.threshold_mismatches
    }

    /// Incorporate one assessment without laundering any non-certified state.
    pub fn observe(&mut self, assessment: &SimplificationAssessment) {
        let is_first_assessment = self.assessments == 0;
        self.assessments += 1;
        self.size_before += assessment.size_before;
        self.size_after += assessment.size_after;
        if assessment.size_after < assessment.size_before {
            self.simplified_count += 1;
        }
        if let Some(bound) = assessment.certified_error {
            self.max_certified_error = self.max_certified_error.max(bound);
        }
        if let Some(actual) = assessment.sampled_discrepancy {
            self.max_sampled_discrepancy = self.max_sampled_discrepancy.max(actual);
        }

        let threshold_mismatch =
            assessment.radius_threshold.to_bits() != self.radius_threshold.to_bits();
        let mut incoming_status = assessment.status;
        if threshold_mismatch {
            self.threshold_mismatches += 1;
            incoming_status = incoming_status.combine(SimplificationCheckStatus::ThresholdMismatch);
        }
        match assessment.status {
            SimplificationCheckStatus::Certified => {}
            SimplificationCheckStatus::StructuralEmptyAgreement => {
                self.structural_empty_agreements += 1;
            }
            SimplificationCheckStatus::SimplifierRefused => {
                self.simplifier_refusals += 1;
            }
            SimplificationCheckStatus::NonFiniteCertificate => {
                self.non_finite_certificates += 1;
            }
            SimplificationCheckStatus::NegativeCertificate => {
                self.negative_certificates += 1;
            }
            SimplificationCheckStatus::DiscrepancyEvidenceRefused => {
                self.discrepancy_evidence_refusals += 1;
            }
            SimplificationCheckStatus::CertificateCheckExceeded => {
                self.certificate_check_exceedances += 1;
            }
            SimplificationCheckStatus::ThresholdMismatch => {
                if !threshold_mismatch {
                    self.threshold_mismatches += 1;
                }
            }
        }
        self.status = if is_first_assessment {
            incoming_status
        } else {
            self.status.combine(incoming_status)
        };
    }

    /// Whether every observed result is admitted and certificate-sound.
    #[must_use]
    pub const fn is_sound(&self) -> bool {
        self.assessments > 0 && self.status.is_sound() && self.threshold_mismatches == 0
    }

    /// Whether the aggregate is both sound and complete for an expected set.
    /// This prevents a sound strict subset from authorizing an “all checked”
    /// campaign claim.
    #[must_use]
    pub const fn is_complete_and_sound(&self, expected_assessments: usize) -> bool {
        self.is_sound() && self.assessments == expected_assessments
    }
}

/// The campaign report.
#[derive(Debug, Clone)]
pub struct GrammarReport {
    /// Fraction of (size × fab-margin) niches filled.
    pub coverage: f64,
    /// Quality-diversity score (Σ elite fitness, where fitness is `1/(1+discrepancy)`).
    pub qd_score: f64,
    /// Number of filled niches.
    pub num_elites: usize,
    /// The best (lowest) SDF discrepancy from the target.
    pub best_discrepancy: f64,
    /// The best program's parameters `[r1, r2, d, o]`.
    pub best_params: [f64; 4],
    /// Shared, typed simplification/certificate accounting.
    pub simplification: SimplificationSummary,
    /// Compatibility mirror of `simplification.simplified_count()`.
    pub simplified_count: usize,
    /// Compatibility mirror of `simplification.size_before()`.
    pub size_before: usize,
    /// Compatibility mirror of `simplification.size_after()`.
    pub size_after: usize,
    /// Compatibility mirror of `simplification.max_certified_error()`.
    pub max_certified_error: f64,
    /// Compatibility mirror of
    /// `simplification.is_complete_and_sound(num_elites)`.
    pub simplification_sound: bool,
    /// Elites that satisfy the minimum-feature-size fabrication constraint.
    pub fab_satisfied: usize,
    /// The headline color: `Verified` iff the best design matches within
    /// tolerance, is fab-satisfied, and simplifies soundly.
    pub headline_color: Color,
}

/// Run the GrammarForge campaign. `match_tol` is the target-match discrepancy
/// threshold; `simplify_radius_threshold` is only the strict local offset-radius
/// admission threshold, not the returned compositional error envelope.
#[must_use]
pub fn run_campaign(match_tol: f64, simplify_radius_threshold: f64) -> GrammarReport {
    let target = target();
    let samples = sample_points();
    // A minimum-feature-size rule that actually discriminates: the thinnest
    // (r=0.7) spheres fail it, so not every program is fabricable.
    let fab = min_feature_size(0.8);

    let r_vals = [0.7, 0.9, 1.0, 1.1];
    let d_vals = [0.6, 0.8, 1.0];
    let o_vals = [0.0, 0.02, 0.05];

    // MAP-Elites over (total material `r1+r2`, dipole separation `d`) — the
    // behavioral axes of the peanut family.
    let mut archive = MapElites::new(vec![1.3, 0.5], vec![2.3, 1.1], vec![6, 4]);
    for &r1 in &r_vals {
        for &r2 in &r_vals {
            for &d in &d_vals {
                for &o in &o_vals {
                    let prog = build_program(r1, r2, d, o);
                    let disc = max_sdf_discrepancy(&prog, &target, &samples);
                    // Closeness score in (0, 1] — higher is a better match; the
                    // archive requires a non-negative fitness.
                    let fitness = 1.0 / (1.0 + disc);
                    let descriptor = vec![r1 + r2, d];
                    archive.add(vec![r1, r2, d, o], descriptor, fitness);
                }
            }
        }
    }

    // Post-process the elites: simplification soundness + fabrication tally.
    let mut simplification = SimplificationSummary::new(simplify_radius_threshold);
    let mut fab_satisfied = 0;
    for e in archive.elites() {
        let prog = build_program(e.solution[0], e.solution[1], e.solution[2], e.solution[3]);
        let assessment = assess_simplification(&prog, simplify_radius_threshold, &samples);
        simplification.observe(&assessment);
        if fab.satisfied(e.solution[0].min(e.solution[1])) {
            fab_satisfied += 1;
        }
    }

    let best = archive.best().expect("archive has at least one elite");
    let best_discrepancy = 1.0 / best.fitness - 1.0;
    let best_params = [
        best.solution[0],
        best.solution[1],
        best.solution[2],
        best.solution[3],
    ];
    let best_fab_ok = fab.satisfied(best.solution[0].min(best.solution[1]));
    let simplification_sound = simplification.is_complete_and_sound(archive.num_elites());
    let headline_color = if best_discrepancy <= match_tol && best_fab_ok && simplification_sound {
        // declared-color-ok: demo headline candidate from local discrepancy/fabricability checks; admitted only at a consumer's authority boundary (6pf9)
        Color::Verified {
            lo: 0.0,
            hi: best_discrepancy,
        }
    } else {
        Color::Estimated {
            estimator: "grammar-open".to_string(),
            dispersion: best_discrepancy,
        }
    };

    GrammarReport {
        coverage: archive.coverage(),
        qd_score: archive.qd_score(),
        num_elites: archive.num_elites(),
        best_discrepancy,
        best_params,
        simplification,
        simplified_count: simplification.simplified_count,
        size_before: simplification.size_before,
        size_after: simplification.size_after,
        max_certified_error: simplification.max_certified_error,
        simplification_sound,
        fab_satisfied,
        headline_color,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn defensive_aggregate_states_are_counted_and_fail_closed() {
        let non_finite_certificate = SimplificationAssessment {
            radius_threshold: 0.03,
            size_before: 2,
            size_after: 1,
            certified_error: None,
            rejected_certificate_bits: Some(f64::NAN.to_bits()),
            sampled_discrepancy: None,
            status: SimplificationCheckStatus::NonFiniteCertificate,
            refusal: None,
        };
        let exceeded_certificate_check = SimplificationAssessment {
            radius_threshold: 0.03,
            size_before: 3,
            size_after: 2,
            certified_error: Some(0.01),
            rejected_certificate_bits: None,
            sampled_discrepancy: Some(0.02),
            status: SimplificationCheckStatus::CertificateCheckExceeded,
            refusal: None,
        };
        let negative_certificate = SimplificationAssessment {
            radius_threshold: 0.03,
            size_before: 2,
            size_after: 2,
            certified_error: None,
            rejected_certificate_bits: Some((-0.01_f64).to_bits()),
            sampled_discrepancy: None,
            status: SimplificationCheckStatus::NegativeCertificate,
            refusal: None,
        };

        let mut summary = SimplificationSummary::new(0.03);
        summary.observe(&non_finite_certificate);
        summary.observe(&negative_certificate);
        summary.observe(&exceeded_certificate_check);

        assert_eq!(summary.assessments(), 3);
        assert_eq!(summary.simplified_count(), 2);
        assert_eq!(summary.non_finite_certificates(), 1);
        assert_eq!(summary.negative_certificates(), 1);
        assert_eq!(
            non_finite_certificate.rejected_certificate_bits(),
            Some(f64::NAN.to_bits())
        );
        assert_eq!(
            negative_certificate.rejected_certificate_bits(),
            Some((-0.01_f64).to_bits())
        );
        assert_eq!(summary.certificate_check_exceedances(), 1);
        assert_eq!(summary.max_certified_error().to_bits(), 0.01_f64.to_bits());
        assert_eq!(
            summary.max_sampled_discrepancy().to_bits(),
            0.02_f64.to_bits()
        );
        assert_eq!(
            summary.status(),
            SimplificationCheckStatus::CertificateCheckExceeded,
            "the highest-severity state must remain visible"
        );
        assert!(!summary.is_sound());
        assert!(!summary.is_complete_and_sound(3));
    }
}
