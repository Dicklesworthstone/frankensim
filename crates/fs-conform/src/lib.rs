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
    /// An optional identity witness: probes on which a converter CLAIMED to be
    /// the identity map must act as one (`source_dim == target_dim`,
    /// `apply(x) == x`). `None` for converters that are not identities.
    pub identity: Option<Vec<Vec<f64>>>,
    /// Numerical tolerance for the axiom checks.
    pub tolerance: f64,
}

impl ConformanceSuite<'_> {
    /// An incomplete empty suite with the given numerical tolerance. Populate
    /// adjoint and manufactured evidence before calling [`certify`].
    #[must_use]
    pub fn new(tolerance: f64) -> ConformanceSuite<'static> {
        ConformanceSuite {
            adjoint_pairs: Vec::new(),
            manufactured: Vec::new(),
            composition: None,
            identity: None,
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

fn valid_tolerance(tol: f64) -> bool {
    tol.is_finite() && tol >= 0.0
}

fn finite_vector(values: &[f64]) -> bool {
    values.iter().all(|value| value.is_finite())
}

fn dot(a: &[f64], b: &[f64]) -> Option<f64> {
    if a.len() != b.len() || !finite_vector(a) || !finite_vector(b) {
        return None;
    }
    let mut total = 0.0;
    for (&x, &y) in a.iter().zip(b) {
        let product = x * y;
        if !product.is_finite() || (x != 0.0 && y != 0.0 && product == 0.0) {
            return None;
        }
        let next = total + product;
        if !next.is_finite()
            || (product != 0.0 && next == total)
            || (total != 0.0 && next == product)
        {
            return None;
        }
        total = next;
    }
    Some(total)
}

fn dist(a: &[f64], b: &[f64]) -> Option<f64> {
    if a.len() != b.len() || !finite_vector(a) || !finite_vector(b) {
        return None;
    }
    let mut squared = 0.0;
    for (&x, &y) in a.iter().zip(b) {
        let delta = x - y;
        let term = delta * delta;
        if !delta.is_finite() || !term.is_finite() || (delta != 0.0 && term == 0.0) {
            return None;
        }
        let next = squared + term;
        if !next.is_finite() || (term != 0.0 && next == squared) || (squared != 0.0 && next == term)
        {
            return None;
        }
        squared = next;
    }
    let distance = squared.sqrt();
    distance.is_finite().then_some(distance)
}

/// Check adjoint consistency `⟨A x, y⟩ == ⟨x, Aᵀ y⟩` over the pairs.
#[must_use]
pub fn check_adjoint(c: &dyn Converter, pairs: &[(Vec<f64>, Vec<f64>)], tol: f64) -> bool {
    if pairs.is_empty() || !valid_tolerance(tol) {
        return false;
    }
    let (source_dim, target_dim) = (c.source_dim(), c.target_dim());
    pairs.iter().all(|(x, y)| {
        if x.len() != source_dim || y.len() != target_dim || !finite_vector(x) || !finite_vector(y)
        {
            return false;
        }
        let applied = c.apply(x);
        let adjoint = c.adjoint(y);
        if applied.len() != target_dim || adjoint.len() != source_dim {
            return false;
        }
        let (Some(lhs), Some(rhs)) = (dot(&applied, y), dot(x, &adjoint)) else {
            return false;
        };
        let delta = lhs - rhs;
        delta.is_finite() && delta.abs() <= tol
    })
}

fn check_tolerance_honesty_with_declared(
    c: &dyn Converter,
    cases: &[ManufacturedCase],
    tol: f64,
    declared: f64,
) -> (bool, f64) {
    if cases.is_empty() || !valid_tolerance(tol) || !declared.is_finite() || declared < 0.0 {
        return (false, f64::INFINITY);
    }
    let admitted_bound = declared + tol;
    if !admitted_bound.is_finite() {
        return (false, f64::INFINITY);
    }
    let (source_dim, target_dim) = (c.source_dim(), c.target_dim());
    let mut measured = 0.0_f64;
    for case in cases {
        if case.input.len() != source_dim
            || case.exact_output.len() != target_dim
            || !finite_vector(&case.input)
            || !finite_vector(&case.exact_output)
        {
            return (false, f64::INFINITY);
        }
        let applied = c.apply(&case.input);
        if applied.len() != target_dim {
            return (false, f64::INFINITY);
        }
        let Some(error) = dist(&applied, &case.exact_output) else {
            return (false, f64::INFINITY);
        };
        measured = measured.max(error);
    }
    (measured <= admitted_bound, measured)
}

/// Check tolerance honesty; returns `(honest, worst_measured_error)`.
#[must_use]
pub fn check_tolerance_honesty(
    c: &dyn Converter,
    cases: &[ManufacturedCase],
    tol: f64,
) -> (bool, f64) {
    check_tolerance_honesty_with_declared(c, cases, tol, c.declared_error())
}

/// Check functoriality: `after(self(x)) == direct(x)` on the probes.
#[must_use]
pub fn check_functoriality(c: &dyn Converter, comp: &Composition, tol: f64) -> bool {
    if comp.probes.is_empty()
        || !valid_tolerance(tol)
        || c.target_dim() != comp.after.source_dim()
        || c.source_dim() != comp.direct.source_dim()
        || comp.after.target_dim() != comp.direct.target_dim()
    {
        return false;
    }
    let (source_dim, middle_dim, target_dim) =
        (c.source_dim(), c.target_dim(), comp.after.target_dim());
    comp.probes.iter().all(|x| {
        if x.len() != source_dim || !finite_vector(x) {
            return false;
        }
        let middle = c.apply(x);
        if middle.len() != middle_dim || !finite_vector(&middle) {
            return false;
        }
        let composed = comp.after.apply(&middle);
        let direct = comp.direct.apply(x);
        if composed.len() != target_dim || direct.len() != target_dim {
            return false;
        }
        dist(&composed, &direct).is_some_and(|distance| distance <= tol)
    })
}

/// Check that a converter claiming to be an identity acts as one.
#[must_use]
pub fn check_identity(c: &dyn Converter, probes: &[Vec<f64>], tol: f64) -> bool {
    if probes.is_empty() || !valid_tolerance(tol) || c.source_dim() != c.target_dim() {
        return false;
    }
    let dim = c.source_dim();
    probes.iter().all(|x| {
        if x.len() != dim || !finite_vector(x) {
            return false;
        }
        let applied = c.apply(x);
        applied.len() == dim && dist(&applied, x).is_some_and(|distance| distance <= tol)
    })
}

/// Certify a converter against its suite. It reaches a tier ABOVE `Rejected`
/// only by passing every supplied axiom; the tier level reflects how tight an
/// (honestly met) error model it declares. Adjoint and manufactured evidence
/// must be non-empty; any supplied composition or identity witness must carry
/// at least one probe.
#[must_use]
pub fn certify(c: &dyn Converter, suite: &ConformanceSuite) -> ConformanceReport {
    let mut findings = Vec::new();
    let declared_error = c.declared_error();
    if !valid_tolerance(suite.tolerance) || !declared_error.is_finite() || declared_error < 0.0 {
        findings.push(
            "admission: tolerance and declared error must be finite and non-negative".to_string(),
        );
        return ConformanceReport {
            converter: c.id().to_string(),
            functoriality: false,
            adjoint_consistent: false,
            tolerance_honest: false,
            measured_error: f64::INFINITY,
            tier: Tier::Rejected,
            findings,
        };
    }

    // Functoriality: composition agrees AND (if the converter claims to be an
    // identity) it acts as the identity.
    let composition_ok = match &suite.composition {
        Some(comp) if comp.probes.is_empty() => {
            findings.push("functoriality: supplied composition has no probes".to_string());
            false
        }
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
    let identity_ok = match &suite.identity {
        Some(probes) if probes.is_empty() => {
            findings.push("identity: supplied identity witness has no probes".to_string());
            false
        }
        Some(probes) => {
            let ok = check_identity(c, probes, suite.tolerance);
            if !ok {
                findings.push(format!(
                    "identity: {} claims to be an identity but apply(x) != x",
                    c.id()
                ));
            }
            ok
        }
        None => true,
    };
    let functoriality = composition_ok && identity_ok;

    let adjoint_consistent =
        !suite.adjoint_pairs.is_empty() && check_adjoint(c, &suite.adjoint_pairs, suite.tolerance);
    if !adjoint_consistent {
        findings.push(if suite.adjoint_pairs.is_empty() {
            "adjoint consistency: no witness pairs supplied".to_string()
        } else {
            "adjoint consistency: <Ax,y> != <x,Aᵀy> (declared transpose is not the adjoint)"
                .to_string()
        });
    }

    let (tolerance_honest, measured_error) = if suite.manufactured.is_empty() {
        (false, f64::INFINITY)
    } else {
        check_tolerance_honesty_with_declared(
            c,
            &suite.manufactured,
            suite.tolerance,
            declared_error,
        )
    };
    if !tolerance_honest {
        findings.push(if suite.manufactured.is_empty() {
            "tolerance honesty: no manufactured cases supplied".to_string()
        } else {
            format!(
                "tolerance honesty: measured error {measured_error:.3e} exceeds declared {:.3e}",
                declared_error
            )
        });
    }

    let tier = if functoriality && adjoint_consistent && tolerance_honest {
        tier_for_error(declared_error)
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
    if !declared.is_finite() || declared < 0.0 {
        Tier::Rejected
    } else if declared <= 1e-6 {
        Tier::Gold
    } else if declared <= 1e-3 {
        Tier::Silver
    } else {
        Tier::Bronze
    }
}
