//! fs-conform — the restriction-map plugin conformance SDK (plan addendum,
//! Proposal 7). Layer: L2.
//!
//! The restriction-map layer is where the hard engineering hides: the sheaf
//! organizes bookkeeping GIVEN the trace/conversion operators, and it will
//! faithfully propagate GARBAGE WITH CERTIFICATES ATTACHED if those operators
//! are bad (risk R6). This crate turns that weakest point into an ecosystem
//! play: third parties ship [`Converter`]s (chart-to-chart operators — the Rep
//! Router edges), and a CONFORMANCE SUITE auto-generated from the sheaf axioms
//! certifies each converter into a [`Tier`]. A converter reaches a tier ONLY by
//! passing all three axioms:
//!
//! 1. **Functoriality** — composition agrees (`f∘g == direct`) and identities
//!    act as identities.
//! 2. **Adjoint consistency** — the declared transpose really is the adjoint
//!    the ledger uses: `⟨A x, y⟩ == ⟨x, Aᵀ y⟩`.
//! 3. **Tolerance honesty** — against MANUFACTURED solutions with known
//!    interface traces, the measured error must not exceed the converter's
//!    DECLARED error model. A converter that understates its error FAILS.
//!
//! R6 mitigation: [`certify`] is applied to FIRST-PARTY converters with the
//! same severity as third-party ones. The certified tier is meant to be stamped
//! on every ledger entry the converter touches. Deterministic; no dependencies.

/// A chart-to-chart converter (a Rep Router edge / restriction map). Kept
/// object-safe so heterogeneous third-party converters share one SDK surface.
pub trait Converter {
    /// A stable id (stamped alongside the tier on ledger entries).
    fn id(&self) -> &str;
    /// The source chart dimension.
    fn source_dim(&self) -> usize;
    /// The target chart dimension.
    fn target_dim(&self) -> usize;
    /// Apply the conversion (source → target).
    fn apply(&self, x: &[f64]) -> Vec<f64>;
    /// The DECLARED adjoint/transpose (target → source).
    fn adjoint(&self, y: &[f64]) -> Vec<f64>;
    /// The DECLARED error bound of the converter's error model.
    fn declared_error(&self) -> f64;
}

/// A manufactured solution: an input with its KNOWN exact converted output.
#[derive(Debug, Clone, PartialEq)]
pub struct ManufacturedCase {
    /// The source-chart input.
    pub input: Vec<f64>,
    /// The known-exact target-chart output.
    pub exact_output: Vec<f64>,
}

/// A functoriality witness: `after ∘ self` must equal `direct` on `probes`.
pub struct Composition<'a> {
    /// The converter applied AFTER `self` (target → C).
    pub after: &'a dyn Converter,
    /// The claimed direct converter (source → C).
    pub direct: &'a dyn Converter,
    /// Source-chart probe vectors.
    pub probes: Vec<Vec<f64>>,
}

/// The conformance suite for one converter.
pub struct ConformanceSuite<'a> {
    /// `(x, y)` pairs (source, target) for the adjoint identity.
    pub adjoint_pairs: Vec<(Vec<f64>, Vec<f64>)>,
    /// Manufactured tolerance-honesty cases.
    pub manufactured: Vec<ManufacturedCase>,
    /// An optional functoriality witness (composition).
    pub composition: Option<Composition<'a>>,
    /// Numerical tolerance for the axiom checks.
    pub tolerance: f64,
}

impl ConformanceSuite<'_> {
    /// An empty suite with the given numerical tolerance.
    #[must_use]
    pub fn new(tolerance: f64) -> ConformanceSuite<'static> {
        ConformanceSuite {
            adjoint_pairs: Vec::new(),
            manufactured: Vec::new(),
            composition: None,
            tolerance,
        }
    }
}

/// The certified conformance tier.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Tier {
    /// Failed a hard axiom — NOT certified (do not trust its certificates).
    Rejected,
    /// Certified, coarse declared error.
    Bronze,
    /// Certified, tight declared error.
    Silver,
    /// Certified, very tight declared error.
    Gold,
}

/// The conformance report for a converter.
#[derive(Debug, Clone, PartialEq)]
pub struct ConformanceReport {
    /// The converter id.
    pub converter: String,
    /// Did composition/identity hold? (`true` if no witness supplied.)
    pub functoriality: bool,
    /// Did the adjoint identity hold?
    pub adjoint_consistent: bool,
    /// Did the declared error model contain the measured error?
    pub tolerance_honest: bool,
    /// The worst measured manufactured-solution error.
    pub measured_error: f64,
    /// The awarded tier.
    pub tier: Tier,
    /// Human-readable findings (reasons for any failure).
    pub findings: Vec<String>,
}

impl ConformanceReport {
    /// Was the converter certified (any tier above `Rejected`)?
    #[must_use]
    pub fn certified(&self) -> bool {
        self.tier != Tier::Rejected
    }
}

fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(x, y)| x * y).sum()
}

fn dist(a: &[f64], b: &[f64]) -> f64 {
    a.iter()
        .zip(b)
        .map(|(x, y)| (x - y) * (x - y))
        .sum::<f64>()
        .sqrt()
}

/// Check adjoint consistency `⟨A x, y⟩ == ⟨x, Aᵀ y⟩` over the pairs.
#[must_use]
pub fn check_adjoint(c: &dyn Converter, pairs: &[(Vec<f64>, Vec<f64>)], tol: f64) -> bool {
    pairs.iter().all(|(x, y)| {
        let lhs = dot(&c.apply(x), y);
        let rhs = dot(x, &c.adjoint(y));
        (lhs - rhs).abs() <= tol
    })
}

/// Check tolerance honesty; returns `(honest, worst_measured_error)`.
#[must_use]
pub fn check_tolerance_honesty(
    c: &dyn Converter,
    cases: &[ManufacturedCase],
    tol: f64,
) -> (bool, f64) {
    let measured = cases
        .iter()
        .map(|m| dist(&c.apply(&m.input), &m.exact_output))
        .fold(0.0_f64, f64::max);
    (measured <= c.declared_error() + tol, measured)
}

/// Check functoriality: `after(self(x)) == direct(x)` on the probes.
#[must_use]
pub fn check_functoriality(c: &dyn Converter, comp: &Composition, tol: f64) -> bool {
    comp.probes.iter().all(|x| {
        let composed = comp.after.apply(&c.apply(x));
        let direct = comp.direct.apply(x);
        composed.len() == direct.len() && dist(&composed, &direct) <= tol
    })
}

/// Check that a converter claiming to be an identity acts as one.
#[must_use]
pub fn check_identity(c: &dyn Converter, probes: &[Vec<f64>], tol: f64) -> bool {
    c.source_dim() == c.target_dim() && probes.iter().all(|x| dist(&c.apply(x), x) <= tol)
}

/// Certify a converter against its suite. It reaches a tier ABOVE `Rejected`
/// only by passing every supplied axiom; the tier level reflects how tight an
/// (honestly met) error model it declares.
#[must_use]
pub fn certify(c: &dyn Converter, suite: &ConformanceSuite) -> ConformanceReport {
    let mut findings = Vec::new();

    let functoriality = match &suite.composition {
        Some(comp) => {
            let ok = check_functoriality(c, comp, suite.tolerance);
            if !ok {
                findings.push(format!(
                    "functoriality: {} ∘ {} != direct",
                    comp.after.id(),
                    c.id()
                ));
            }
            ok
        }
        None => true,
    };

    let adjoint_consistent = check_adjoint(c, &suite.adjoint_pairs, suite.tolerance);
    if !adjoint_consistent {
        findings.push(
            "adjoint consistency: <Ax,y> != <x,Aᵀy> (declared transpose is not the adjoint)"
                .to_string(),
        );
    }

    let (tolerance_honest, measured_error) =
        check_tolerance_honesty(c, &suite.manufactured, suite.tolerance);
    if !tolerance_honest {
        findings.push(format!(
            "tolerance honesty: measured error {measured_error:.3e} exceeds declared {:.3e}",
            c.declared_error()
        ));
    }

    let tier = if functoriality && adjoint_consistent && tolerance_honest {
        tier_for_error(c.declared_error())
    } else {
        Tier::Rejected
    };

    ConformanceReport {
        converter: c.id().to_string(),
        functoriality,
        adjoint_consistent,
        tolerance_honest,
        measured_error,
        tier,
        findings,
    }
}

/// The tier awarded to a converter that passed every axiom, by declared error.
fn tier_for_error(declared: f64) -> Tier {
    if declared <= 1e-6 {
        Tier::Gold
    } else if declared <= 1e-3 {
        Tier::Silver
    } else {
        Tier::Bronze
    }
}
