//! Operating envelopes, duty cycles, and load-combination EVALUATION.
//!
//! [`scenario`](crate::scenario) *declares* load cases and factored
//! combinations; it does not run them. This module is the evaluation half:
//!
//! - [`OperatingEnvelope`] — named axes (ambient temperature, total power,
//!   pressure, discrete fan states including failure modes) with declared
//!   co-occurrence limits between axis pairs;
//! - [`OperatingEnvelope::enumerate_corners`] — the deterministic worst-case
//!   corner set, which for a coupled envelope is **not** the naive hypercube;
//! - [`EnvelopeDutyCycle`] — dwell- or fraction-weighted profiles over
//!   envelope points, with steady-weighted aggregation labelled as the
//!   approximation it is;
//! - [`compile_case_set`] — factored [`Combination`]s compiled into a
//!   [`CaseSet`] carrying per-case provenance, plus [`governing_case`], which
//!   answers the question engineers actually ask: *which condition sizes the
//!   design?*
//!
//! # What an envelope is, stated narrowly
//!
//! An envelope is a declared **support**: the set of operating points the
//! author asserts the system may be required to occupy. It is not a
//! probability distribution, and a corner is not a likelihood. Corner
//! enumeration answers "what are the extreme admissible points", which is a
//! worst-case question; it says nothing about how often any of them occur.
//!
//! # Co-occurrence limits are not correlation coefficients
//!
//! A correlation coefficient constrains *moments*; it does not constrain the
//! *support*. Two axes can be strongly correlated and still jointly attain
//! every corner of their box. So a declared `rho` cannot license dropping a
//! corner, and this module refuses to accept one: the only thing that shrinks
//! a corner set here is a [`CouplingRelation::CoOccurrenceLimit`], which is an
//! explicit statement that a region of the box is *unreachable*, carrying the
//! author's rationale for why.
//!
//! Declaring the co-occurrence *unknown* is a first-class option
//! ([`CouplingRelation::Unknown`]) and is the honest default when nobody has
//! characterised the pair. It never shrinks the corner set — ignorance is not
//! evidence of infeasibility — and it emits a caveat that travels into the
//! error budget.

use core::cmp::Ordering;
use core::fmt;

use fs_la::{LuWorkspace, LuWorkspaceError};
use fs_math::dd::Dd;
use fs_qty::{Dims, QtyAny};

use crate::scenario::{Combination, LoadCase, Violation};

/// Temperature dimensions (K).
const TEMP_DIMS: Dims = Dims([0, 0, 0, 1, 0, 0]);
/// Power dimensions (W = m^2 kg s^-3).
const POWER_DIMS: Dims = Dims([2, 1, -3, 0, 0, 0]);
/// Time dimensions (s).
const TIME_DIMS: Dims = Dims([0, 0, 1, 0, 0, 0]);

/// Relative tolerance for accepting a clipped vertex as lying on a box edge.
///
/// A co-occurrence line meeting a box edge produces a coordinate that is
/// mathematically inside the edge span but can land a rounding step outside
/// it. Admitting the point and clamping it to the span keeps a real corner;
/// rejecting it would silently *drop* an extreme operating point, which is the
/// unsafe direction.
pub const EDGE_SPAN_TOLERANCE: f64 = 1e-12;

/// Maximum active constraint sets tried by the default shared-axis enumerator.
///
/// This bounds numerical work independently of the number of retained corners.
/// Call [`OperatingEnvelope::enumerate_corners_with_candidate_limit`] to choose
/// a different limit. Independent axes and single-limit pairs keep their
/// existing analytic enumeration and do not consume this allowance.
pub const DEFAULT_ENVELOPE_CANDIDATE_LIMIT: usize = 262_144;

/// Relative tolerance on the sum of declared duty fractions.
pub const DUTY_FRACTION_TOLERANCE: f64 = 1e-12;

/// Resource admission for the combinatorial products in this module.
///
/// Independent corner blocks and case-set crossing form Cartesian products,
/// whose size is bounded before the full product is allocated. An envelope with ten
/// five-corner blocks describes ten million corners; the refusal names the
/// count so the author can see what they asked for. Shared-axis components
/// additionally have a separate active-set work allowance.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct EnvelopeBudget {
    /// Maximum enumerated corners.
    pub max_corners: usize,
    /// Maximum compiled cases.
    pub max_cases: usize,
    /// Maximum axes in one envelope.
    pub max_axes: usize,
}

/// Default admission: generous for hand-authored envelopes, far below the
/// point where a product becomes an accidental denial of service.
pub const DEFAULT_ENVELOPE_BUDGET: EnvelopeBudget = EnvelopeBudget {
    max_corners: 65_536,
    max_cases: 262_144,
    max_axes: 64,
};

impl Default for EnvelopeBudget {
    fn default() -> Self {
        DEFAULT_ENVELOPE_BUDGET
    }
}

/// How a discrete axis state relates to the design intent.
///
/// The distinction is load-bearing for reporting: a QoI governed by a `Failed`
/// state is a fault-tolerance result, and a report that cannot tell it apart
/// from a nominal result is misleading.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum StateKind {
    /// The intended operating configuration.
    Nominal,
    /// Reduced capability, still within intended operation.
    Degraded,
    /// A declared failure mode.
    Failed,
}

impl fmt::Display for StateKind {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let text = match self {
            StateKind::Nominal => "nominal",
            StateKind::Degraded => "degraded",
            StateKind::Failed => "failed",
        };
        f.write_str(text)
    }
}

/// One declared state of a discrete axis.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DiscreteState {
    /// State name, unique within its axis.
    pub name: String,
    /// Relation to design intent.
    pub kind: StateKind,
}

impl DiscreteState {
    /// A nominal state.
    #[must_use]
    pub fn nominal(name: impl Into<String>) -> Self {
        DiscreteState {
            name: name.into(),
            kind: StateKind::Nominal,
        }
    }

    /// A declared failure state.
    #[must_use]
    pub fn failed(name: impl Into<String>) -> Self {
        DiscreteState {
            name: name.into(),
            kind: StateKind::Failed,
        }
    }
}

/// The admissible values of one envelope axis.
#[derive(Debug, Clone, PartialEq)]
pub enum AxisDomain {
    /// A closed interval in the axis's own units.
    Continuous {
        /// Lower bound.
        low: QtyAny,
        /// Upper bound, strictly above `low`.
        high: QtyAny,
    },
    /// An enumerated set of named states.
    Discrete {
        /// Declared states; exactly one must be [`StateKind::Nominal`].
        states: Vec<DiscreteState>,
    },
}

/// A named axis of an operating envelope.
#[derive(Debug, Clone, PartialEq)]
pub struct EnvelopeAxis {
    /// Axis name, unique within the envelope.
    pub name: String,
    /// Admissible values.
    pub domain: AxisDomain,
}

impl EnvelopeAxis {
    /// A continuous axis over `[low, high]` in the given dimensions.
    #[must_use]
    pub fn continuous(name: impl Into<String>, low: QtyAny, high: QtyAny) -> Self {
        EnvelopeAxis {
            name: name.into(),
            domain: AxisDomain::Continuous { low, high },
        }
    }

    /// A discrete axis over declared states.
    #[must_use]
    pub fn discrete(name: impl Into<String>, states: Vec<DiscreteState>) -> Self {
        EnvelopeAxis {
            name: name.into(),
            domain: AxisDomain::Discrete { states },
        }
    }

    /// The axis's dimensions, or `None` for a discrete axis.
    #[must_use]
    pub fn dims(&self) -> Option<Dims> {
        match &self.domain {
            AxisDomain::Continuous { low, .. } => Some(low.dims),
            AxisDomain::Discrete { .. } => None,
        }
    }
}

/// What is declared about how two axes co-occur.
#[derive(Debug, Clone, PartialEq)]
pub enum CouplingRelation {
    /// The two axes are declared independent: every combination of their
    /// extremes is admissible. This is a *declaration*, not a default — an
    /// undeclared pair is [`CouplingRelation::Unknown`].
    DeclaredIndependent,
    /// Co-occurrence is not characterised.
    ///
    /// No corner may be dropped: the corner set stays the independence
    /// product, which is a superset of the true reachable set. That is
    /// conservative for a worst-case question and inadmissible as the basis
    /// for any probability-weighted statement.
    Unknown {
        /// Why it is unknown, and what would characterise it.
        rationale: String,
    },
    /// A declared linear limit on simultaneous attainment:
    /// `a_coeff * a + b_coeff * b <= bound`.
    ///
    /// Coefficients carry dimensions, so `a_coeff * a`, `b_coeff * b` and
    /// `bound` must all share one dimension. Point admission evaluates the
    /// declared closed inequality without a numerical slack allowance.
    CoOccurrenceLimit {
        /// Coefficient on the first axis.
        a_coeff: QtyAny,
        /// Coefficient on the second axis.
        b_coeff: QtyAny,
        /// Right-hand side.
        bound: QtyAny,
        /// The physical reason the region beyond the limit is unreachable.
        rationale: String,
    },
}

/// A declared relation between two named axes.
#[derive(Debug, Clone, PartialEq)]
pub struct AxisCoupling {
    /// First axis name.
    pub a: String,
    /// Second axis name.
    pub b: String,
    /// What is declared about the pair.
    pub relation: CouplingRelation,
}

/// A caveat that must travel with any result derived from the envelope.
///
/// Distinct from [`Violation`]: a violation means the declaration is wrong and
/// nothing may proceed; a caveat means the result is admissible but carries a
/// stated limitation into the error budget.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct EnvelopeCaveat {
    /// Stable machine-readable code.
    pub code: &'static str,
    /// What the limitation is.
    pub what: String,
    /// What it means for anything computed from this result.
    pub consequence: String,
}

/// A declared operating envelope.
#[derive(Debug, Clone, PartialEq)]
pub struct OperatingEnvelope {
    /// Envelope name.
    pub name: String,
    /// Declared axes.
    pub axes: Vec<EnvelopeAxis>,
    /// Declared pairwise couplings. Shared axes intersect all declared limits.
    pub couplings: Vec<AxisCoupling>,
}

/// One coordinate of an envelope point.
#[derive(Debug, Clone, PartialEq)]
pub enum AxisPoint {
    /// A value on a continuous axis.
    Continuous(QtyAny),
    /// A named state of a discrete axis.
    Discrete(String),
}

/// A point in the envelope: one coordinate per declared axis.
#[derive(Debug, Clone, PartialEq)]
pub struct EnvelopePoint {
    /// `(axis name, value)` in the envelope's declared axis order.
    pub coordinates: Vec<(String, AxisPoint)>,
}

impl EnvelopePoint {
    /// Look up one coordinate by axis name.
    #[must_use]
    pub fn get(&self, axis: &str) -> Option<&AxisPoint> {
        self.coordinates
            .iter()
            .find(|(name, _)| name == axis)
            .map(|(_, value)| value)
    }

    /// The continuous value on `axis`, or `None` if absent or discrete.
    #[must_use]
    pub fn continuous(&self, axis: &str) -> Option<QtyAny> {
        match self.get(axis) {
            Some(AxisPoint::Continuous(q)) => Some(*q),
            _ => None,
        }
    }

    /// The discrete state name on `axis`, or `None` if absent or continuous.
    #[must_use]
    pub fn discrete(&self, axis: &str) -> Option<&str> {
        match self.get(axis) {
            Some(AxisPoint::Discrete(state)) => Some(state.as_str()),
            _ => None,
        }
    }
}

/// The deterministic worst-case corner set of an envelope.
#[derive(Debug, Clone, PartialEq)]
pub struct CornerSet {
    /// The corners, in deterministic enumeration order.
    pub corners: Vec<EnvelopePoint>,
    /// Per-block explanation of how the count arose. Logged so a reader can
    /// see *why* a coupled envelope has the corners it has rather than taking
    /// the enumerator's word for it.
    pub rationale: Vec<String>,
    /// Limitations that travel into the budget.
    pub caveats: Vec<EnvelopeCaveat>,
}

impl OperatingEnvelope {
    /// Validate the declaration.
    ///
    /// Returns every finding, each with a fix hint, in a deterministic order.
    /// An empty vector establishes declaration validity; enumeration separately
    /// checks joint feasibility, numerical resolution and caller budgets.
    #[must_use]
    pub fn validate(&self) -> Vec<Violation> {
        self.validate_with_axis_limit(DEFAULT_ENVELOPE_BUDGET.max_axes)
    }

    #[allow(clippy::too_many_lines)]
    fn validate_with_axis_limit(&self, max_axes: usize) -> Vec<Violation> {
        let mut findings = Vec::new();

        if self.name.trim().is_empty() {
            findings.push(Violation {
                code: "envelope-name-empty",
                what: "the envelope has no name".to_string(),
                fix: "give the envelope a name; it is cited by every case the envelope generates"
                    .to_string(),
            });
        }

        if self.axes.is_empty() {
            findings.push(Violation {
                code: "envelope-axes-empty",
                what: format!("envelope `{}` declares no axes", self.name),
                fix: "declare at least one axis, or drop the envelope and state the single \
                      operating point directly"
                    .to_string(),
            });
        }

        if self.axes.len() > max_axes {
            findings.push(Violation {
                code: "envelope-axes-overflow",
                what: format!(
                    "envelope `{}` declares {} axes, above the admitted maximum of {}",
                    self.name,
                    self.axes.len(),
                    max_axes
                ),
                fix: "split the envelope, or raise max_axes on an explicit EnvelopeBudget"
                    .to_string(),
            });
        }

        for (index, axis) in self.axes.iter().enumerate() {
            if axis.name.trim().is_empty() {
                findings.push(Violation {
                    code: "envelope-axis-name-empty",
                    what: format!("axis {index} has no name"),
                    fix: "name the axis; couplings and points reference axes by name".to_string(),
                });
            }
            if self.axes[..index]
                .iter()
                .any(|other| other.name == axis.name)
            {
                findings.push(Violation {
                    code: "envelope-axis-name-duplicate",
                    what: format!("axis name `{}` is declared more than once", axis.name),
                    fix: "rename one of them; axis names must be unique within an envelope"
                        .to_string(),
                });
            }
            findings.extend(validate_domain(&axis.name, &axis.domain));
        }

        for coupling in &self.couplings {
            findings.extend(self.validate_coupling(coupling));
        }

        findings
    }

    #[allow(clippy::too_many_lines)]
    fn validate_coupling(&self, coupling: &AxisCoupling) -> Vec<Violation> {
        let mut findings = Vec::new();

        if coupling.a == coupling.b {
            findings.push(Violation {
                code: "envelope-coupling-self",
                what: format!("axis `{}` is coupled to itself", coupling.a),
                fix: "couple two distinct axes, or drop the coupling".to_string(),
            });
            return findings;
        }

        let (Some(axis_a), Some(axis_b)) = (self.axis(&coupling.a), self.axis(&coupling.b)) else {
            for name in [&coupling.a, &coupling.b] {
                if self.axis(name).is_none() {
                    findings.push(Violation {
                        code: "envelope-coupling-axis-missing",
                        what: format!("coupling references undeclared axis `{name}`"),
                        fix: "declare the axis, or correct the name in the coupling".to_string(),
                    });
                }
            }
            return findings;
        };

        let CouplingRelation::CoOccurrenceLimit {
            a_coeff,
            b_coeff,
            bound,
            rationale,
        } = &coupling.relation
        else {
            if let CouplingRelation::Unknown { rationale } = &coupling.relation
                && rationale.trim().is_empty()
            {
                findings.push(Violation {
                    code: "envelope-coupling-unknown-unexplained",
                    what: format!(
                        "coupling `{}`/`{}` is declared unknown with no rationale",
                        coupling.a, coupling.b
                    ),
                    fix: "state why it is uncharacterised and what would characterise it; the \
                          caveat is reported to the reader and an empty one tells them nothing"
                        .to_string(),
                });
            }
            return findings;
        };

        if rationale.trim().is_empty() {
            findings.push(Violation {
                code: "envelope-limit-unexplained",
                what: format!(
                    "the co-occurrence limit on `{}`/`{}` has no rationale",
                    coupling.a, coupling.b
                ),
                fix: "state the physical reason the excluded region is unreachable. A limit that \
                      drops worst-case corners without a stated reason is an unaudited assumption"
                    .to_string(),
            });
        }

        let (Some(dims_a), Some(dims_b)) = (axis_a.dims(), axis_b.dims()) else {
            findings.push(Violation {
                code: "envelope-limit-discrete-axis",
                what: format!(
                    "the co-occurrence limit on `{}`/`{}` references a discrete axis; a linear \
                     limit is not expressible over named states",
                    coupling.a, coupling.b
                ),
                fix: "express the restriction by declaring the state set that is actually \
                      reachable, or split into one envelope per discrete state"
                    .to_string(),
            });
            return findings;
        };

        for (label, coeff, axis_dims) in [
            (coupling.a.as_str(), a_coeff, dims_a),
            (coupling.b.as_str(), b_coeff, dims_b),
        ] {
            match coeff.dims.checked_plus(axis_dims) {
                Some(product) if product == bound.dims => {}
                Some(product) => findings.push(Violation {
                    code: "envelope-limit-dims",
                    what: format!(
                        "coefficient on `{label}` gives term dimensions {:?} but the bound has \
                         {:?}",
                        product.0, bound.dims.0
                    ),
                    fix: "every term and the bound must share one dimension; scale the \
                          coefficient so coefficient x axis matches the bound"
                        .to_string(),
                }),
                None => findings.push(Violation {
                    code: "envelope-limit-dims-overflow",
                    what: format!(
                        "coefficient on `{label}` overflows the dimension exponents when \
                         multiplied by the axis"
                    ),
                    fix: "use a coefficient whose exponents stay in range".to_string(),
                }),
            }
        }

        for (label, value) in [
            ("a_coeff", a_coeff.value),
            ("b_coeff", b_coeff.value),
            ("bound", bound.value),
        ] {
            if !value.is_finite() {
                findings.push(Violation {
                    code: "envelope-limit-nonfinite",
                    what: format!(
                        "{label} of the `{}`/`{}` limit is non-finite ({value})",
                        coupling.a, coupling.b
                    ),
                    fix: "supply finite coefficients and bound".to_string(),
                });
            }
        }

        if a_coeff.value == 0.0 && b_coeff.value == 0.0 {
            findings.push(Violation {
                code: "envelope-limit-degenerate",
                what: format!(
                    "the `{}`/`{}` limit has both coefficients zero, so it constrains nothing \
                     about either axis",
                    coupling.a, coupling.b
                ),
                fix: "give at least one axis a nonzero coefficient, or declare the pair \
                      independent"
                    .to_string(),
            });
        }

        if findings.is_empty()
            && let (
                AxisDomain::Continuous {
                    low: a_lo,
                    high: a_hi,
                },
                AxisDomain::Continuous {
                    low: b_lo,
                    high: b_hi,
                },
            ) = (&axis_a.domain, &axis_b.domain)
        {
            let corners = [
                (a_lo.value, b_lo.value),
                (a_hi.value, b_lo.value),
                (a_hi.value, b_hi.value),
                (a_lo.value, b_hi.value),
            ];
            let admitted = corners
                .iter()
                .filter(|(a, b)| a_coeff.value * a + b_coeff.value * b <= bound.value)
                .count();
            if admitted == 0 {
                findings.push(Violation {
                    code: "envelope-limit-empty",
                    what: format!(
                        "the `{}`/`{}` limit excludes every corner of the declared box, so the \
                         envelope has no admissible operating points",
                        coupling.a, coupling.b
                    ),
                    fix: "loosen the bound or widen the axis ranges; an envelope that admits \
                          nothing cannot be run"
                        .to_string(),
                });
            }
        }

        findings
    }

    /// The axis with this name, if declared.
    #[must_use]
    pub fn axis(&self, name: &str) -> Option<&EnvelopeAxis> {
        self.axes.iter().find(|axis| axis.name == name)
    }

    /// Check that a point lies in the envelope.
    ///
    /// # Errors
    /// Returns every reason the point is inadmissible: a missing or extra
    /// coordinate, a kind or dimension mismatch, an out-of-range value, an
    /// undeclared discrete state, or a violated co-occurrence limit.
    pub fn admits(&self, point: &EnvelopePoint) -> Result<(), Vec<Violation>> {
        let mut findings = Vec::new();

        for (name, _) in &point.coordinates {
            if self.axis(name).is_none() {
                findings.push(Violation {
                    code: "envelope-point-axis-unknown",
                    what: format!("point names axis `{name}`, which the envelope does not declare"),
                    fix: "remove the coordinate, or declare the axis".to_string(),
                });
            }
        }

        for axis in &self.axes {
            let matches: Vec<&AxisPoint> = point
                .coordinates
                .iter()
                .filter(|(name, _)| name == &axis.name)
                .map(|(_, value)| value)
                .collect();
            match matches.as_slice() {
                [] => findings.push(Violation {
                    code: "envelope-point-axis-missing",
                    what: format!("point has no coordinate on axis `{}`", axis.name),
                    fix: "supply a value for every declared axis; a missing coordinate is an \
                          undeclared default"
                        .to_string(),
                }),
                [value] => findings.extend(admit_coordinate(axis, value)),
                _ => findings.push(Violation {
                    code: "envelope-point-axis-repeated",
                    what: format!(
                        "point supplies {} coordinates on axis `{}`",
                        matches.len(),
                        axis.name
                    ),
                    fix: "supply exactly one coordinate per axis".to_string(),
                }),
            }
        }

        if findings.is_empty() {
            for coupling in &self.couplings {
                if let CouplingRelation::CoOccurrenceLimit {
                    a_coeff,
                    b_coeff,
                    bound,
                    rationale,
                } = &coupling.relation
                    && let (Some(a), Some(b)) =
                        (point.continuous(&coupling.a), point.continuous(&coupling.b))
                {
                    let lhs = a_coeff.value * a.value + b_coeff.value * b.value;
                    if let Some(admitted) =
                        limit_contains(a_coeff.value, a.value, b_coeff.value, b.value, bound.value)
                    {
                        if admitted {
                            continue;
                        }
                        findings.push(Violation {
                            code: "envelope-point-limit-violated",
                            what: format!(
                                "point violates the `{}`/`{}` co-occurrence limit: {lhs} exceeds \
                                 the bound {} ({rationale})",
                                coupling.a, coupling.b, bound.value
                            ),
                            fix: "move the point inside the declared limit, or revise the limit \
                                  if the point is genuinely reachable"
                                .to_string(),
                        });
                    } else {
                        findings.push(Violation {
                            code: "envelope-point-limit-nonfinite",
                            what: format!("the `{}`/`{}` limit cannot be evaluated finitely at this point", coupling.a, coupling.b),
                            fix: "rescale the coefficients and bound; nonfinite arithmetic cannot establish membership".to_string(),
                        });
                    }
                }
            }
        }

        if findings.is_empty() {
            Ok(())
        } else {
            Err(findings)
        }
    }

    /// Enumerate the deterministic worst-case corner set.
    ///
    /// The envelope decomposes into independent blocks — an unpaired
    /// continuous axis, a connected set of limited axes, or a discrete axis — and
    /// the corner set is the Cartesian product of the blocks' own extreme
    /// points. A coupled pair contributes the vertices of its box clipped by
    /// the declared limit, which is where a coupled envelope stops agreeing
    /// with the naive hypercube: clipping a corner off a rectangle can *add* a
    /// vertex (four becomes five) as easily as remove one.
    ///
    /// # Errors
    /// Returns declaration, output-budget, candidate-work or numerical findings.
    /// Shared-axis work uses [`DEFAULT_ENVELOPE_CANDIDATE_LIMIT`]; callers can
    /// choose a different cap with [`Self::enumerate_corners_with_candidate_limit`].
    pub fn enumerate_corners(&self, budget: EnvelopeBudget) -> Result<CornerSet, Vec<Violation>> {
        self.enumerate_corners_with_candidate_limit(budget, DEFAULT_ENVELOPE_CANDIDATE_LIMIT)
    }

    /// Enumerate with an explicit total active-set allowance for shared limits.
    ///
    /// A connected n-axis component with m distinct inequalities (including
    /// its box bounds) tries each of the `m choose n` active sets. The complete
    /// count is admitted before solving; exact singular sets are not isolated
    /// vertices. Nonfinite or poorly conditioned solves refuse rather than
    /// silently discarding a possible extreme point. Returned shared-axis
    /// vertices are numerical estimates with an explicit caveat, not interval
    /// certificates. Independent blocks retain their original arithmetic.
    ///
    /// # Errors
    /// Invalid declarations, joint infeasibility, work/output limits, allocation
    /// failure, or ambiguous numerical solves publish no partial corner set.
    #[allow(clippy::too_many_lines)]
    pub fn enumerate_corners_with_candidate_limit(
        &self,
        budget: EnvelopeBudget,
        max_candidates: usize,
    ) -> Result<CornerSet, Vec<Violation>> {
        let findings = self.validate_with_axis_limit(budget.max_axes);
        if !findings.is_empty() {
            return Err(findings);
        }

        // Only explicit inequalities join blocks. An Unknown or independence
        // declaration must neither hide another limit nor consume either axis.
        let mut components: Vec<usize> = (0..self.axes.len()).collect();
        for coupling in &self.couplings {
            if matches!(
                coupling.relation,
                CouplingRelation::CoOccurrenceLimit { .. }
            ) {
                let a = self
                    .axes
                    .iter()
                    .position(|axis| axis.name == coupling.a)
                    .expect("validated axis");
                let b = self
                    .axes
                    .iter()
                    .position(|axis| axis.name == coupling.b)
                    .expect("validated axis");
                let (keep, merge) = (components[a], components[b]);
                for component in &mut components {
                    if *component == merge {
                        *component = keep;
                    }
                }
            }
        }

        let mut blocks: Vec<Vec<Vec<(String, AxisPoint)>>> = Vec::new();
        let mut rationale = Vec::new();
        let mut caveats = Vec::new();
        let mut consumed: Vec<&str> = Vec::new();
        let mut remaining_candidates = max_candidates;
        // Keyed on the exact axis pair, not on a substring of the caveat text:
        // axis names nest (`power` inside `total-power`), so a `contains`
        // check would silently suppress a second pair's caveat.
        let mut caveated: Vec<(&str, &str)> = Vec::new();

        for axis in &self.axes {
            for coupling in &self.couplings {
                if coupling.a != axis.name && coupling.b != axis.name {
                    continue;
                }
                if let CouplingRelation::Unknown { rationale: why } = &coupling.relation
                    && !caveated.contains(&(coupling.a.as_str(), coupling.b.as_str()))
                    && !caveated.contains(&(coupling.b.as_str(), coupling.a.as_str()))
                {
                    caveated.push((coupling.a.as_str(), coupling.b.as_str()));
                    let shares_limit = self.couplings.iter().any(|other| {
                        matches!(other.relation, CouplingRelation::CoOccurrenceLimit { .. })
                            && [other.a.as_str(), other.b.as_str()]
                                .iter()
                                .any(|name| *name == coupling.a || *name == coupling.b)
                    });
                    caveats.push(EnvelopeCaveat {
                        code: "envelope-coupling-unknown",
                        what: format!(
                            "co-occurrence of `{}` and `{}` is declared unknown: {why}",
                            coupling.a, coupling.b
                        ),
                        consequence: if shares_limit {
                            "this declaration drops no corners; only explicit limits \
                                      restrict the support, which remains a SUPERSET of the true \
                                      reachable set. That is conservative for a worst-case claim \
                                      and inadmissible as the basis for any probability-weighted \
                                      statement over the pair."
                        } else {
                            "no corner was dropped for this pair, so the corner set \
                             is the independence product — a SUPERSET of the true \
                             reachable set. That is conservative for a worst-case \
                             claim and inadmissible as the basis for any \
                             probability-weighted statement over the pair."
                        }
                        .to_string(),
                    });
                }
            }
        }

        for (axis_index, axis) in self.axes.iter().enumerate() {
            if consumed.contains(&axis.name.as_str()) {
                continue;
            }
            let component: Vec<&EnvelopeAxis> = self
                .axes
                .iter()
                .enumerate()
                .filter(|(index, _)| components[*index] == components[axis_index])
                .map(|(_, member)| member)
                .collect();
            let limits: Vec<&AxisCoupling> = self
                .couplings
                .iter()
                .filter(|coupling| {
                    matches!(
                        coupling.relation,
                        CouplingRelation::CoOccurrenceLimit { .. }
                    ) && component.iter().any(|member| member.name == coupling.a)
                })
                .collect();
            match limits.as_slice() {
                [coupling] => {
                    let CouplingRelation::CoOccurrenceLimit {
                        a_coeff,
                        b_coeff,
                        bound,
                        rationale: why,
                    } = &coupling.relation
                    else {
                        unreachable!("guarded by the match arm above")
                    };
                    // `coupling.a` and `coupling.b` fix the coefficient order;
                    // the iteration order of `self.axes` does not.
                    let axis_a = self
                        .axis(&coupling.a)
                        .expect("validated: coupling axes are declared");
                    let axis_b = self
                        .axis(&coupling.b)
                        .expect("validated: coupling axes are declared");
                    let (
                        AxisDomain::Continuous {
                            low: a_lo,
                            high: a_hi,
                        },
                        AxisDomain::Continuous {
                            low: b_lo,
                            high: b_hi,
                        },
                    ) = (&axis_a.domain, &axis_b.domain)
                    else {
                        unreachable!("validated: a limit's axes are continuous")
                    };

                    let vertices = clip_box_to_limit(
                        (a_lo.value, a_hi.value),
                        (b_lo.value, b_hi.value),
                        a_coeff.value,
                        b_coeff.value,
                        bound.value,
                    );
                    // Exact equality is the right test, not a tolerance: a
                    // surviving box corner is COPIED verbatim from the axis
                    // bound, never computed, so it is bit-identical or it is a
                    // different point. A tolerance here would count a vertex
                    // that merely landed near a corner as being that corner,
                    // and understate how much the limit changed.
                    #[allow(clippy::float_cmp)]
                    let dropped = 4 - vertices
                        .iter()
                        .filter(|(a, b)| {
                            (*a == a_lo.value || *a == a_hi.value)
                                && (*b == b_lo.value || *b == b_hi.value)
                        })
                        .count();
                    rationale.push(format!(
                        "axes `{}`x`{}`: box corners 4 -> {} vertices after the declared \
                         co-occurrence limit ({dropped} box corner(s) unreachable, {} vertex/ices \
                         introduced on the limit boundary). Rationale: {why}",
                        coupling.a,
                        coupling.b,
                        vertices.len(),
                        vertices.len() + dropped - 4
                    ));

                    let block = vertices
                        .into_iter()
                        .map(|(a, b)| {
                            vec![
                                (
                                    coupling.a.clone(),
                                    AxisPoint::Continuous(QtyAny::new(a, a_lo.dims)),
                                ),
                                (
                                    coupling.b.clone(),
                                    AxisPoint::Continuous(QtyAny::new(b, b_lo.dims)),
                                ),
                            ]
                        })
                        .collect();
                    blocks.push(block);
                    consumed.push(coupling.a.as_str());
                    consumed.push(coupling.b.as_str());
                }
                [] => {
                    let block = axis_block(axis);
                    rationale.push(format!(
                        "axis `{}`: {} extreme point(s), uncoupled",
                        axis.name,
                        block.len()
                    ));
                    blocks.push(block);
                    consumed.push(axis.name.as_str());
                }
                _ => {
                    let (vertices, tried) = clip_coupled_box(
                        &component,
                        &limits,
                        budget.max_corners,
                        remaining_candidates,
                    )
                    .map_err(|finding| vec![finding])?;
                    remaining_candidates -= tried;
                    let names = component
                        .iter()
                        .map(|member| member.name.as_str())
                        .collect::<Vec<_>>()
                        .join(" x ");
                    rationale.push(format!(
                        "axes `{names}`: {} joint-limit vertices from {tried} active sets; \
                         {} declared co-occurrence limits intersected",
                        vertices.len(),
                        limits.len(),
                    ));
                    for coupling in &limits {
                        if let CouplingRelation::CoOccurrenceLimit { rationale: why, .. } =
                            &coupling.relation
                        {
                            rationale
                                .push(format!("limit `{}`/`{}`: {why}", coupling.a, coupling.b));
                        }
                    }
                    caveats.push(EnvelopeCaveat {
                        code: "envelope-corners-numerical",
                        what: format!("shared-axis vertices on `{names}` use normalized floating-point active-set solves"),
                        consequence: "finite residual and conditioning checks are numerical admission, \
                                      not a certified enclosure or a proof of a nonlinear QoI extremum"
                            .to_string(),
                    });
                    let block = vertices
                        .into_iter()
                        .map(|values| {
                            component
                                .iter()
                                .zip(values)
                                .map(|(member, value)| {
                                    (
                                        member.name.clone(),
                                        AxisPoint::Continuous(QtyAny::new(
                                            value,
                                            member.dims().expect("limited axis is continuous"),
                                        )),
                                    )
                                })
                                .collect()
                        })
                        .collect();
                    blocks.push(block);
                    consumed.extend(component.iter().map(|member| member.name.as_str()));
                }
            }
        }

        let mut total: usize = 1;
        for block in &blocks {
            total = match total.checked_mul(block.len()) {
                Some(product) => product,
                None => usize::MAX,
            };
            if total > budget.max_corners {
                return Err(vec![Violation {
                    code: "envelope-corner-overflow",
                    what: format!(
                        "the corner product of envelope `{}` exceeds the admitted maximum of {} \
                         corners",
                        self.name, budget.max_corners
                    ),
                    fix: "reduce the axis count, split the envelope, or raise max_corners on an \
                          explicit EnvelopeBudget. The product is refused before allocation, so \
                          nothing was enumerated"
                        .to_string(),
                }]);
            }
        }

        let mut corners: Vec<EnvelopePoint> = vec![EnvelopePoint {
            coordinates: Vec::new(),
        }];
        for block in blocks {
            let mut next = Vec::with_capacity(corners.len() * block.len());
            for prefix in &corners {
                for choice in &block {
                    let mut coordinates = prefix.coordinates.clone();
                    coordinates.extend(choice.iter().cloned());
                    next.push(EnvelopePoint { coordinates });
                }
            }
            corners = next;
        }

        // Present coordinates in declared axis order regardless of the block
        // order the couplings imposed, so two envelopes that differ only in
        // coupling declaration order produce byte-comparable points.
        for corner in &mut corners {
            corner.coordinates.sort_by_key(|(name, _)| {
                self.axes
                    .iter()
                    .position(|axis| &axis.name == name)
                    .unwrap_or(usize::MAX)
            });
        }

        // Order the corners themselves canonically too, for the same reason
        // one step up: a corner's INDEX is what a compiled case cites as its
        // provenance, so if swapping the two names in a coupling declaration
        // permuted the list, it would silently renumber every case in every
        // report while changing nothing physical. Continuous axes order by
        // value; discrete axes order by declared state position, so a fan
        // ordering stays as the author wrote it rather than going alphabetical.
        corners.sort_by(|left, right| {
            for axis in &self.axes {
                let ordering = match (left.get(&axis.name), right.get(&axis.name)) {
                    (Some(AxisPoint::Continuous(a)), Some(AxisPoint::Continuous(b))) => {
                        a.value.total_cmp(&b.value)
                    }
                    (Some(AxisPoint::Discrete(a)), Some(AxisPoint::Discrete(b))) => {
                        let position = |name: &String| match &axis.domain {
                            AxisDomain::Discrete { states } => states
                                .iter()
                                .position(|state| &state.name == name)
                                .unwrap_or(usize::MAX),
                            AxisDomain::Continuous { .. } => usize::MAX,
                        };
                        position(a).cmp(&position(b))
                    }
                    _ => Ordering::Equal,
                };
                if ordering != Ordering::Equal {
                    return ordering;
                }
            }
            Ordering::Equal
        });

        rationale.push(format!("corner product: {} point(s)", corners.len()));

        Ok(CornerSet {
            corners,
            rationale,
            caveats,
        })
    }
}

/// The extreme points of one uncoupled axis.
fn axis_block(axis: &EnvelopeAxis) -> Vec<Vec<(String, AxisPoint)>> {
    match &axis.domain {
        AxisDomain::Continuous { low, high } => vec![
            vec![(axis.name.clone(), AxisPoint::Continuous(*low))],
            vec![(axis.name.clone(), AxisPoint::Continuous(*high))],
        ],
        AxisDomain::Discrete { states } => states
            .iter()
            .map(|state| vec![(axis.name.clone(), AxisPoint::Discrete(state.name.clone()))])
            .collect(),
    }
}

fn validate_domain(axis_name: &str, domain: &AxisDomain) -> Vec<Violation> {
    let mut findings = Vec::new();
    match domain {
        AxisDomain::Continuous { low, high } => {
            if low.dims != high.dims {
                findings.push(Violation {
                    code: "envelope-axis-dims",
                    what: format!(
                        "axis `{axis_name}` has bounds with different dimensions ({:?} and {:?})",
                        low.dims.0, high.dims.0
                    ),
                    fix: "give both bounds the same dimensions".to_string(),
                });
            }
            for (label, value) in [("low", low.value), ("high", high.value)] {
                if !value.is_finite() {
                    findings.push(Violation {
                        code: "envelope-axis-nonfinite",
                        what: format!(
                            "axis `{axis_name}` has a non-finite {label} bound ({value})"
                        ),
                        fix: "supply finite bounds; an unbounded axis has no corners".to_string(),
                    });
                }
            }
            if low.value.is_finite() && high.value.is_finite() && low.value >= high.value {
                findings.push(Violation {
                    code: "envelope-axis-empty",
                    what: format!(
                        "axis `{axis_name}` has low {} not below high {}",
                        low.value, high.value
                    ),
                    fix: "give the axis a nonempty range. A pinned value is not an axis: fold it \
                          into the scenario as a declared constant so it is visible where it acts"
                        .to_string(),
                });
            }
        }
        AxisDomain::Discrete { states } => {
            if states.is_empty() {
                findings.push(Violation {
                    code: "envelope-axis-no-states",
                    what: format!("discrete axis `{axis_name}` declares no states"),
                    fix: "declare the states the axis may take".to_string(),
                });
            }
            for (index, state) in states.iter().enumerate() {
                if state.name.trim().is_empty() {
                    findings.push(Violation {
                        code: "envelope-state-name-empty",
                        what: format!("state {index} of axis `{axis_name}` has no name"),
                        fix: "name every state; reports cite states by name".to_string(),
                    });
                }
                if states[..index].iter().any(|other| other.name == state.name) {
                    findings.push(Violation {
                        code: "envelope-state-duplicate",
                        what: format!(
                            "state `{}` of axis `{axis_name}` is declared more than once",
                            state.name
                        ),
                        fix: "state names must be unique within an axis".to_string(),
                    });
                }
            }
            let nominal = states
                .iter()
                .filter(|state| state.kind == StateKind::Nominal)
                .count();
            if !states.is_empty() && nominal != 1 {
                findings.push(Violation {
                    code: "envelope-state-nominal-count",
                    what: format!(
                        "discrete axis `{axis_name}` declares {nominal} nominal states; exactly \
                         one is required"
                    ),
                    fix: "mark exactly one state Nominal. Without it a report cannot say whether \
                          a governing case is a design condition or a fault condition"
                        .to_string(),
                });
            }
        }
    }
    findings
}

fn admit_coordinate(axis: &EnvelopeAxis, value: &AxisPoint) -> Vec<Violation> {
    let mut findings = Vec::new();
    match (&axis.domain, value) {
        (AxisDomain::Continuous { low, high }, AxisPoint::Continuous(q)) => {
            if q.dims != low.dims {
                findings.push(Violation {
                    code: "envelope-point-dims",
                    what: format!(
                        "coordinate on axis `{}` has dimensions {:?} but the axis is {:?}",
                        axis.name, q.dims.0, low.dims.0
                    ),
                    fix: "supply the value in the axis's dimensions".to_string(),
                });
            } else if !q.value.is_finite() {
                findings.push(Violation {
                    code: "envelope-point-nonfinite",
                    what: format!(
                        "coordinate on axis `{}` is non-finite ({})",
                        axis.name, q.value
                    ),
                    fix: "supply a finite value".to_string(),
                });
            } else if q.value < low.value || q.value > high.value {
                findings.push(Violation {
                    code: "envelope-point-out-of-range",
                    what: format!(
                        "coordinate {} on axis `{}` is outside the declared range [{}, {}]",
                        q.value, axis.name, low.value, high.value
                    ),
                    fix: "move the point into the declared range, or widen the axis if the point \
                          is genuinely required"
                        .to_string(),
                });
            }
        }
        (AxisDomain::Discrete { states }, AxisPoint::Discrete(name)) => {
            if !states.iter().any(|state| &state.name == name) {
                findings.push(Violation {
                    code: "envelope-point-state-unknown",
                    what: format!(
                        "coordinate on axis `{}` names state `{name}`, which the axis does not \
                         declare",
                        axis.name
                    ),
                    fix: "use a declared state, or declare the state on the axis".to_string(),
                });
            }
        }
        (AxisDomain::Continuous { .. }, AxisPoint::Discrete(_)) => findings.push(Violation {
            code: "envelope-point-kind",
            what: format!(
                "axis `{}` is continuous but the point supplies a discrete state",
                axis.name
            ),
            fix: "supply a dimensioned value".to_string(),
        }),
        (AxisDomain::Discrete { .. }, AxisPoint::Continuous(_)) => findings.push(Violation {
            code: "envelope-point-kind",
            what: format!(
                "axis `{}` is discrete but the point supplies a continuous value",
                axis.name
            ),
            fix: "supply a declared state name".to_string(),
        }),
    }
    findings
}

/// Vertices of the box `[a_lo, a_hi] x [b_lo, b_hi]` intersected with the
/// closed half-plane `a_coeff * a + b_coeff * b <= bound`.
///
/// Candidate-and-filter rather than sequential polygon clipping: the four box
/// corners that satisfy the limit, plus the limit line's crossings of the four
/// box edges. Duplicates — a crossing that lands exactly on a corner, which is
/// what a limit passing *through* a corner produces — are merged on exact
/// equality after a canonical sort.
fn clip_box_to_limit(
    (a_lo, a_hi): (f64, f64),
    (b_lo, b_hi): (f64, f64),
    a_coeff: f64,
    b_coeff: f64,
    bound: f64,
) -> Vec<(f64, f64)> {
    let admits = |a: f64, b: f64| a_coeff * a + b_coeff * b <= bound;
    let mut vertices: Vec<(f64, f64)> = Vec::new();

    for (a, b) in [(a_lo, b_lo), (a_hi, b_lo), (a_hi, b_hi), (a_lo, b_hi)] {
        if admits(a, b) {
            vertices.push((a, b));
        }
    }

    // Crossings of the two edges of constant `a`: solve for b.
    if b_coeff != 0.0 {
        for a in [a_lo, a_hi] {
            let b = (bound - a_coeff * a) / b_coeff;
            if let Some(b) = clamp_to_span(b, b_lo, b_hi) {
                vertices.push((a, b));
            }
        }
    }
    // Crossings of the two edges of constant `b`: solve for a.
    if a_coeff != 0.0 {
        for b in [b_lo, b_hi] {
            let a = (bound - b_coeff * b) / a_coeff;
            if let Some(a) = clamp_to_span(a, a_lo, a_hi) {
                vertices.push((a, b));
            }
        }
    }

    vertices.sort_by(|left, right| {
        left.0
            .total_cmp(&right.0)
            .then_with(|| left.1.total_cmp(&right.1))
    });
    vertices.dedup();
    vertices
}

/// Admit `value` as lying on `[lo, hi]`, clamping a rounding-step overshoot.
///
/// Returns `None` when the value is outside the span by more than
/// [`EDGE_SPAN_TOLERANCE`] relative to the span width.
fn clamp_to_span(value: f64, lo: f64, hi: f64) -> Option<f64> {
    if !value.is_finite() {
        return None;
    }
    let slack = (hi - lo).abs() * EDGE_SPAN_TOLERANCE;
    if value < lo - slack || value > hi + slack {
        return None;
    }
    Some(value.clamp(lo, hi))
}

/// A row of the shared-axis problem in unit-box coordinates. Rows are sorted
/// canonically and columns follow declared axis order, so reordering the
/// coupling declarations cannot change pivots.
#[derive(PartialEq)]
struct EnvelopeHalfspace {
    coefficients: Vec<f64>,
    bound: f64,
}

fn corner_numerical(what: impl Into<String>) -> Violation {
    Violation {
        code: "envelope-corner-numerical",
        what: what.into(),
        fix: "rescale or simplify the declared limits; no partial corner set is returned \
              when a possible vertex cannot be resolved reliably"
            .to_string(),
    }
}

fn candidate_overflow(limit: usize) -> Violation {
    Violation {
        code: "envelope-candidate-overflow",
        what: format!(
            "shared-axis enumeration exceeds the remaining allowance of {limit} active sets"
        ),
        fix: "split independent operating regions or explicitly raise max_candidates; \
              the complete active-set count is checked before numerical solves"
            .to_string(),
    }
}

fn corner_capacity() -> Violation {
    Violation {
        code: "envelope-corner-capacity",
        what: "corner workspace size overflowed or its allocation failed".to_string(),
        fix: "use fewer coupled axes or a smaller candidate/output budget".to_string(),
    }
}

fn envelope_zeros(count: usize) -> Result<Vec<f64>, Violation> {
    let mut values = Vec::new();
    values
        .try_reserve_exact(count)
        .map_err(|_| corner_capacity())?;
    values.resize(count, 0.0);
    Ok(values)
}

fn coordinate_order(left: &[f64], right: &[f64]) -> Ordering {
    left.iter()
        .zip(right)
        .map(|(a, b)| a.total_cmp(b))
        .find(|order| *order != Ordering::Equal)
        .unwrap_or(Ordering::Equal)
}

/// Finite closed-limit membership, without relaxing the physical declaration.
/// Numerical vertices that cannot retain this after mapping back to physical
/// units refuse publication; a small solve residual is not membership proof.
fn limit_contains(a_coeff: f64, a: f64, b_coeff: f64, b: f64, bound: f64) -> Option<bool> {
    let first = a_coeff * a;
    let second = b_coeff * b;
    let lhs = first + second;
    if !first.is_finite() || !second.is_finite() || !lhs.is_finite() || !bound.is_finite() {
        return None;
    }
    Some(lhs <= bound)
}

/// Count combinations without iterating them or overflowing a machine word.
fn active_set_count(rows: usize, dimension: usize, limit: usize) -> Result<usize, Violation> {
    let mut count = 1_u128;
    for i in 0..dimension.min(rows - dimension) {
        count = count
            .checked_mul((rows - i) as u128)
            .ok_or_else(|| candidate_overflow(limit))?
            / (i + 1) as u128;
        if count > limit as u128 {
            return Err(candidate_overflow(limit));
        }
    }
    if count > limit as u128 {
        return Err(candidate_overflow(limit));
    }
    Ok(count as usize)
}

/// Account for the final rounding of a compensated, dimensionless row value.
/// A small physical bound must not disappear behind large translated origins.
fn rounded_envelope_scalar(value: Dd, error: &mut f64) -> Result<f64, Violation> {
    if !value.hi.is_finite() || !value.lo.is_finite() {
        return Err(corner_numerical("co-occurrence normalization is nonfinite"));
    }
    let rounded = value.to_f64();
    *error += (value - Dd::from_f64(rounded)).abs().to_f64();
    if !error.is_finite() || *error > EDGE_SPAN_TOLERANCE {
        return Err(corner_numerical(
            "co-occurrence coordinate transformation exceeds the normalized span tolerance",
        ));
    }
    Ok(rounded)
}

/// Normalize the physical box before doing linear algebra. Each inequality is
/// scaled by its largest coefficient, never by an arbitrary unit convention.
/// Overflow/underflow in this change of coordinates is a refusal, not a zero
/// coefficient that could silently discard a physical constraint.
fn coupled_halfspaces(
    axes: &[&EnvelopeAxis],
    limits: &[&AxisCoupling],
) -> Result<(Vec<EnvelopeHalfspace>, Vec<(f64, f64)>), Violation> {
    let n = axes.len();
    let mut spans = Vec::new();
    spans.try_reserve_exact(n).map_err(|_| corner_capacity())?;
    for axis in axes {
        let AxisDomain::Continuous { low, high } = axis.domain else {
            unreachable!("validated: limits contain only continuous axes")
        };
        let width = high.value - low.value;
        if !width.is_finite() || width <= 0.0 {
            return Err(corner_numerical(format!(
                "axis `{}` has an unrepresentable span",
                axis.name
            )));
        }
        spans.push((low.value, high.value));
    }
    let count = n
        .checked_mul(2)
        .and_then(|count| count.checked_add(limits.len()))
        .ok_or_else(corner_capacity)?;
    let mut rows = Vec::new();
    rows.try_reserve_exact(count)
        .map_err(|_| corner_capacity())?;
    for column in 0..n {
        for (coefficient, bound) in [(-1.0, 0.0), (1.0, 1.0)] {
            let mut coefficients = envelope_zeros(n)?;
            coefficients[column] = coefficient;
            rows.push(EnvelopeHalfspace {
                coefficients,
                bound,
            });
        }
    }
    for limit in limits {
        let CouplingRelation::CoOccurrenceLimit {
            a_coeff,
            b_coeff,
            bound,
            ..
        } = limit.relation
        else {
            unreachable!("component contains only limits")
        };
        let mut coefficients = envelope_zeros(n)?;
        coefficients[axes
            .iter()
            .position(|axis| axis.name == limit.a)
            .expect("component axis")] = a_coeff.value;
        coefficients[axes
            .iter()
            .position(|axis| axis.name == limit.b)
            .expect("component axis")] = b_coeff.value;
        let mut rhs = Dd::from_f64(bound.value);
        let mut scale = 0.0_f64;
        let mut magnitude = bound.value.abs();
        for (&coefficient, &(low, high)) in coefficients.iter().zip(&spans) {
            let width = Dd::from_f64(high) - Dd::from_f64(low);
            let offset = Dd::from_f64(coefficient) * Dd::from_f64(low);
            let transformed = Dd::from_f64(coefficient) * width;
            rhs = rhs - offset;
            if !rhs.hi.is_finite()
                || !rhs.lo.is_finite()
                || !transformed.hi.is_finite()
                || !transformed.lo.is_finite()
                || (coefficient != 0.0 && transformed.hi == 0.0)
            {
                return Err(corner_numerical(
                    "co-occurrence coordinate transformation is not finite and representable",
                ));
            }
            scale = scale.max(transformed.abs().to_f64());
            magnitude = magnitude
                .max(offset.abs().to_f64())
                .max(transformed.abs().to_f64());
        }
        // Dd's documented per-operation relative error is at most 2^-103.
        // A row has only two nonzero terms. Include a conservative aggregate
        // allowance and a subnormal floor, then charge the f64 row rounding.
        let mut transform_error =
            (32.0 * f64::EPSILON * f64::EPSILON * magnitude + f64::from_bits(64)) / scale;
        for (coefficient, &(low, high)) in coefficients.iter_mut().zip(&spans) {
            let old = *coefficient;
            let width = Dd::from_f64(high) - Dd::from_f64(low);
            *coefficient = rounded_envelope_scalar(
                Dd::from_f64(old) * width / Dd::from_f64(scale),
                &mut transform_error,
            )?;
            if old != 0.0 && *coefficient == 0.0 {
                return Err(corner_numerical(
                    "co-occurrence row scaling underflows a coefficient",
                ));
            }
            if *coefficient == 0.0 {
                *coefficient = 0.0;
            } // Canonical positive zero.
        }
        let scaled_rhs = rounded_envelope_scalar(rhs / Dd::from_f64(scale), &mut transform_error)?;
        if (rhs.hi != 0.0 || rhs.lo != 0.0) && scaled_rhs == 0.0 {
            return Err(corner_numerical(
                "co-occurrence row scaling cannot represent the bound",
            ));
        }
        // A row separated from the whole unit box by more than the arithmetic
        // allowance adds no vertex. Retain touching or unresolved rows, and
        // always keep the original physical check at publication.
        let maximum = coefficients.iter().map(|value| value.max(0.0)).sum::<f64>();
        if scaled_rhs - maximum <= transform_error {
            rows.push(EnvelopeHalfspace {
                coefficients,
                bound: if scaled_rhs == 0.0 { 0.0 } else { scaled_rhs },
            });
        }
    }
    rows.sort_by(|left, right| {
        coordinate_order(&left.coefficients, &right.coefficients)
            .then_with(|| left.bound.total_cmp(&right.bound))
    });
    rows.dedup();
    Ok((rows, spans))
}

/// The active constraints of any vertex of a bounded polytope contain a
/// full-rank n-row subset. Enumerating those subsets therefore also finds
/// corners at intersections of several limits, which pairwise clipping misses.
#[allow(clippy::too_many_lines)]
fn clip_coupled_box(
    axes: &[&EnvelopeAxis],
    limits: &[&AxisCoupling],
    max_corners: usize,
    max_candidates: usize,
) -> Result<(Vec<Vec<f64>>, usize), Violation> {
    let n = axes.len();
    let (rows, spans) = coupled_halfspaces(axes, limits)?;
    let count = active_set_count(rows.len(), n, max_candidates)?;
    let mut lu = LuWorkspace::new(n).map_err(|_| corner_capacity())?;
    let mut matrix = envelope_zeros(n.checked_mul(n).ok_or_else(corner_capacity)?)?;
    let mut rhs = envelope_zeros(n)?;
    let mut solution = envelope_zeros(n)?;
    let mut inverse_column = envelope_zeros(n)?;
    let mut inverse_rows = envelope_zeros(n)?;
    let mut active: Vec<usize> = (0..n).collect();
    let mut vertices: Vec<Vec<f64>> = Vec::new();
    loop {
        for (i, &row) in active.iter().enumerate() {
            matrix[i * n..(i + 1) * n].copy_from_slice(&rows[row].coefficients);
            rhs[i] = rows[row].bound;
        }
        match lu.factor(&matrix) {
            Err(LuWorkspaceError::Singular { .. }) => {}
            Err(error) => {
                return Err(corner_numerical(format!(
                    "active-set factorization refused: {error}"
                )));
            }
            Ok(()) => {
                lu.solve_factored_into(&rhs, &mut solution)
                    .map_err(|error| {
                        corner_numerical(format!("active-set solution refused: {error}"))
                    })?;
                // Compute an inverse norm with the same prepared factor, rather
                // than treating a small residual as a conditioning guarantee.
                inverse_rows.fill(0.0);
                for column in 0..n {
                    rhs.fill(0.0);
                    rhs[column] = 1.0;
                    lu.solve_factored_into(&rhs, &mut inverse_column)
                        .map_err(|error| {
                            corner_numerical(format!("active-set condition solve refused: {error}"))
                        })?;
                    for (sum, value) in inverse_rows.iter_mut().zip(&inverse_column) {
                        *sum += value.abs();
                    }
                }
                let inverse_norm = inverse_rows.iter().copied().fold(0.0_f64, f64::max);
                let matrix_norm = matrix
                    .chunks_exact(n)
                    .map(|row| row.iter().map(|x| x.abs()).sum::<f64>())
                    .fold(0.0_f64, f64::max);
                let condition = matrix_norm * inverse_norm;
                if !condition.is_finite()
                    || condition * f64::EPSILON * n as f64 > EDGE_SPAN_TOLERANCE
                {
                    return Err(corner_numerical(format!(
                        "active-set condition estimate {condition} cannot resolve a vertex at the declared span tolerance"
                    )));
                }
                let residual = active
                    .iter()
                    .enumerate()
                    .map(|(i, &row)| {
                        (matrix[i * n..(i + 1) * n]
                            .iter()
                            .zip(&solution)
                            .fold(0.0, |sum, (a, x)| a.mul_add(*x, sum))
                            - rows[row].bound)
                            .abs()
                    })
                    .fold(0.0_f64, f64::max);
                if !residual.is_finite() || residual * inverse_norm > EDGE_SPAN_TOLERANCE {
                    return Err(corner_numerical(
                        "active-set residual cannot resolve the vertex inside the axis-span tolerance",
                    ));
                }
                let inside_box = solution.iter_mut().all(|coordinate| {
                    if let Some(clamped) = clamp_to_span(*coordinate, 0.0, 1.0) {
                        *coordinate = if clamped == 0.0 { 0.0 } else { clamped };
                        true
                    } else {
                        false
                    }
                });
                let feasible = inside_box
                    && rows.iter().all(|row| {
                        let lhs = row
                            .coefficients
                            .iter()
                            .zip(&solution)
                            .fold(0.0, |sum, (a, x)| a.mul_add(*x, sum));
                        let scale = row
                            .coefficients
                            .iter()
                            .zip(&solution)
                            .map(|(a, x)| (a * x).abs())
                            .sum::<f64>()
                            + row.bound.abs();
                        lhs.is_finite() && lhs - row.bound <= 32.0 * f64::EPSILON * scale
                    });
                if feasible {
                    let mut values = envelope_zeros(n)?;
                    for ((value, &coordinate), &(low, high)) in
                        values.iter_mut().zip(&solution).zip(&spans)
                    {
                        let coordinate = clamp_to_span(coordinate, 0.0, 1.0).ok_or_else(|| {
                            corner_numerical(
                                "feasible active-set solution lies outside the normalized box",
                            )
                        })?;
                        *value = if coordinate == 0.0 {
                            low
                        } else if coordinate == 1.0 {
                            high
                        } else {
                            (Dd::from_f64(low)
                                + (Dd::from_f64(high) - Dd::from_f64(low))
                                    * Dd::from_f64(coordinate))
                            .to_f64()
                        };
                        let restored = (Dd::from_f64(*value) - Dd::from_f64(low))
                            / (Dd::from_f64(high) - Dd::from_f64(low));
                        let round_trip_error = (restored - Dd::from_f64(coordinate)).abs();
                        if !round_trip_error.hi.is_finite()
                            || !round_trip_error.lo.is_finite()
                            || round_trip_error.to_f64() > EDGE_SPAN_TOLERANCE
                        {
                            return Err(corner_numerical(
                                "mapping the vertex to physical units loses its normalized round-trip coordinate",
                            ));
                        }
                    }
                    for limit in limits {
                        let CouplingRelation::CoOccurrenceLimit {
                            a_coeff,
                            b_coeff,
                            bound,
                            ..
                        } = limit.relation
                        else {
                            unreachable!("component contains only limits")
                        };
                        let a = values[axes
                            .iter()
                            .position(|axis| axis.name == limit.a)
                            .expect("component axis")];
                        let b = values[axes
                            .iter()
                            .position(|axis| axis.name == limit.b)
                            .expect("component axis")];
                        if limit_contains(a_coeff.value, a, b_coeff.value, b, bound.value)
                            != Some(true)
                        {
                            return Err(corner_numerical(
                                "mapping the vertex back to physical units did not preserve finite limit membership",
                            ));
                        }
                    }
                    if let Err(index) =
                        vertices.binary_search_by(|other| coordinate_order(other, &values))
                    {
                        if vertices.len() >= max_corners {
                            return Err(Violation {
                                code: "envelope-corner-overflow",
                                what: format!("shared-axis component exceeds {max_corners} retained corners"),
                                fix: "split the envelope or raise max_corners explicitly; no partial corner set is returned".to_string(),
                            });
                        }
                        vertices.try_reserve(1).map_err(|_| corner_capacity())?;
                        vertices.insert(index, values);
                    }
                }
            }
        }
        let Some(next) = (0..n)
            .rev()
            .find(|&index| active[index] < rows.len() - n + index)
        else {
            break;
        };
        active[next] += 1;
        for index in next + 1..n {
            active[index] = active[index - 1] + 1;
        }
    }
    if vertices.is_empty() {
        return Err(Violation {
            code: "envelope-limits-empty",
            what: "the co-occurrence limits jointly leave no admissible box vertex".to_string(),
            fix: "check all limits together; individually nonempty limits may have an empty intersection".to_string(),
        });
    }
    Ok((vertices, count))
}

/// How a duty cycle's per-point weights are declared.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DutyWeighting {
    /// Weights are dwell durations (s). Fractions are DERIVED from them, so
    /// the two can never disagree, and the quasi-steady condition can be
    /// checked because absolute time is known.
    Dwells,
    /// Weights are dimensionless fractions that must sum to one. No absolute
    /// time is declared, so the quasi-steady condition cannot be checked at
    /// this layer.
    Fractions,
}

/// One dwell of a duty cycle.
#[derive(Debug, Clone, PartialEq)]
pub struct DutyPoint {
    /// Point name, unique within the cycle; cited in reports.
    pub name: String,
    /// Where in the envelope the system sits during this dwell.
    pub point: EnvelopePoint,
    /// Dwell duration (s) or dimensionless fraction, per the cycle's
    /// [`DutyWeighting`]. Declared once — never both.
    pub weight: QtyAny,
}

/// A duty cycle over envelope points.
#[derive(Debug, Clone, PartialEq)]
pub struct EnvelopeDutyCycle {
    /// Cycle name.
    pub name: String,
    /// How `weight` is to be read.
    pub weighting: DutyWeighting,
    /// The dwells, in declaration order.
    pub points: Vec<DutyPoint>,
}

/// What a weighted aggregate actually is.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AggregationBasis {
    /// A time-weighted mean of per-point STEADY results. It equals the true
    /// time average only where every dwell is long compared with the system's
    /// thermal time constant.
    SteadyApproximation,
}

/// A duty-weighted quantity of interest.
#[derive(Debug, Clone, PartialEq)]
pub struct WeightedQoi {
    /// The weighted value.
    pub value: f64,
    /// Derived time fractions, in point order; these always sum to one.
    pub fractions: Vec<f64>,
    /// What the number is.
    pub basis: AggregationBasis,
    /// Limitations that travel into the budget.
    pub caveats: Vec<EnvelopeCaveat>,
}

/// Whether a duty cycle's dwells justify treating each point as steady.
#[derive(Debug, Clone, PartialEq)]
pub struct QuasiSteadyVerdict {
    /// The smallest `dwell / time_constant` across the cycle.
    pub min_dwell_ratio: f64,
    /// The point that attains it.
    pub shortest_point: String,
    /// The ratio the caller required.
    pub required_ratio: f64,
    /// Whether every dwell met the requirement.
    pub satisfied: bool,
}

impl EnvelopeDutyCycle {
    /// Validate the cycle against the envelope it runs over.
    ///
    /// # Errors
    /// Returns every finding: an unnamed or duplicated point, a weight with
    /// the wrong dimensions for the declared weighting, a non-positive weight,
    /// declared fractions that do not sum to one, or a point the envelope does
    /// not admit.
    pub fn validate_against(&self, envelope: &OperatingEnvelope) -> Result<(), Vec<Violation>> {
        let mut findings = Vec::new();

        if self.name.trim().is_empty() {
            findings.push(Violation {
                code: "duty-name-empty",
                what: "the duty cycle has no name".to_string(),
                fix: "name the cycle; aggregated results cite it".to_string(),
            });
        }
        if self.points.is_empty() {
            findings.push(Violation {
                code: "duty-points-empty",
                what: format!("duty cycle `{}` declares no points", self.name),
                fix: "declare at least one dwell".to_string(),
            });
        }

        let expected_dims = match self.weighting {
            DutyWeighting::Dwells => TIME_DIMS,
            DutyWeighting::Fractions => Dims::NONE,
        };

        for (index, dwell) in self.points.iter().enumerate() {
            if dwell.name.trim().is_empty() {
                findings.push(Violation {
                    code: "duty-point-name-empty",
                    what: format!("point {index} of duty cycle `{}` has no name", self.name),
                    fix: "name every dwell".to_string(),
                });
            }
            if self.points[..index]
                .iter()
                .any(|other| other.name == dwell.name)
            {
                findings.push(Violation {
                    code: "duty-point-name-duplicate",
                    what: format!(
                        "point name `{}` appears more than once in duty cycle `{}`",
                        dwell.name, self.name
                    ),
                    fix: "point names must be unique within a cycle".to_string(),
                });
            }
            if dwell.weight.dims != expected_dims {
                findings.push(Violation {
                    code: "duty-weight-dims",
                    what: format!(
                        "point `{}` has weight dimensions {:?} but the cycle declares {:?} \
                         weighting, which requires {:?}",
                        dwell.name, dwell.weight.dims.0, self.weighting, expected_dims.0
                    ),
                    fix: "supply dwell durations in seconds for Dwells weighting, or \
                          dimensionless fractions for Fractions weighting"
                        .to_string(),
                });
            }
            if !dwell.weight.value.is_finite() || dwell.weight.value <= 0.0 {
                findings.push(Violation {
                    code: "duty-weight-range",
                    what: format!(
                        "point `{}` has weight {}, which is not strictly positive and finite",
                        dwell.name, dwell.weight.value
                    ),
                    fix: "give every dwell a positive weight. A zero-weight dwell contributes \
                          nothing and should be removed, not carried"
                        .to_string(),
                });
            }
            if let Err(mut reasons) = envelope.admits(&dwell.point) {
                for reason in &mut reasons {
                    reason.what = format!("duty point `{}`: {}", dwell.name, reason.what);
                }
                findings.extend(reasons);
            }
        }

        if self.weighting == DutyWeighting::Fractions && findings.is_empty() {
            let sum: f64 = self.points.iter().map(|dwell| dwell.weight.value).sum();
            if (sum - 1.0).abs() > DUTY_FRACTION_TOLERANCE {
                findings.push(Violation {
                    code: "duty-fractions-sum",
                    what: format!(
                        "declared fractions of duty cycle `{}` sum to {sum}, not one",
                        self.name
                    ),
                    fix: "normalise the fractions, or declare Dwells weighting and let the \
                          fractions be derived"
                        .to_string(),
                });
            }
        }

        if findings.is_empty() {
            Ok(())
        } else {
            Err(findings)
        }
    }

    /// The derived time fractions, which always sum to one.
    ///
    /// # Errors
    /// Returns a finding when the weights do not admit normalisation.
    pub fn fractions(&self) -> Result<Vec<f64>, Vec<Violation>> {
        let total: f64 = self.points.iter().map(|dwell| dwell.weight.value).sum();
        if !total.is_finite() || total <= 0.0 {
            return Err(vec![Violation {
                code: "duty-weight-total",
                what: format!(
                    "the weights of duty cycle `{}` total {total}, which cannot be normalised",
                    self.name
                ),
                fix: "give every dwell a positive finite weight".to_string(),
            }]);
        }
        Ok(self
            .points
            .iter()
            .map(|dwell| dwell.weight.value / total)
            .collect())
    }

    /// Aggregate one QoI over the cycle, weighted by dwell.
    ///
    /// `per_point` supplies the steady QoI value at each dwell, in the cycle's
    /// point order.
    ///
    /// # Errors
    /// Returns findings when `per_point` has the wrong length, contains a
    /// non-finite value, or the weights cannot be normalised.
    pub fn weighted_aggregate(&self, per_point: &[f64]) -> Result<WeightedQoi, Vec<Violation>> {
        if per_point.len() != self.points.len() {
            return Err(vec![Violation {
                code: "duty-qoi-length",
                what: format!(
                    "duty cycle `{}` has {} points but {} QoI values were supplied",
                    self.name,
                    self.points.len(),
                    per_point.len()
                ),
                fix: "supply exactly one value per dwell, in declaration order".to_string(),
            }]);
        }
        for (dwell, value) in self.points.iter().zip(per_point) {
            if !value.is_finite() {
                return Err(vec![Violation {
                    code: "duty-qoi-nonfinite",
                    what: format!("QoI at point `{}` is non-finite ({value})", dwell.name),
                    fix: "supply finite QoI values; a failed solve must be reported as a refusal, \
                          not aggregated"
                        .to_string(),
                }]);
            }
        }

        let fractions = self.fractions()?;
        let value = fractions
            .iter()
            .zip(per_point)
            .map(|(fraction, qoi)| fraction * qoi)
            .sum();

        let mut caveats = vec![EnvelopeCaveat {
            code: "duty-steady-approximation",
            what: format!(
                "the aggregate over duty cycle `{}` is a dwell-weighted mean of per-point STEADY \
                 results",
                self.name
            ),
            consequence: "it equals the true time average only where every dwell is long \
                          compared with the system's thermal time constant. Transient evaluation \
                          of a duty cycle is fs_conduction::duty; this layer weights steady \
                          results and does not integrate anything."
                .to_string(),
        }];

        if self.weighting == DutyWeighting::Fractions {
            caveats.push(EnvelopeCaveat {
                code: "duty-no-absolute-time",
                what: format!(
                    "duty cycle `{}` declares dimensionless fractions, so no absolute dwell \
                     duration is known",
                    self.name
                ),
                consequence: "the quasi-steady condition behind the aggregate CANNOT be checked \
                              at all — not merely unmet, but unmeasurable. Declare Dwells \
                              weighting to make it checkable."
                    .to_string(),
            });
        }

        Ok(WeightedQoi {
            value,
            fractions,
            basis: AggregationBasis::SteadyApproximation,
            caveats,
        })
    }

    /// Measure the quasi-steady condition the weighted aggregate assumes.
    ///
    /// This is the check the [`AggregationBasis::SteadyApproximation`] caveat
    /// refers to: it quantifies how much dwell each point gets in units of the
    /// system's thermal time constant, rather than asserting that the
    /// approximation holds.
    ///
    /// # Errors
    /// Returns a finding when the cycle declares fractions rather than dwells
    /// (there is then no absolute time to compare), when `time_constant` is
    /// not a positive finite duration, or when the cycle has no points.
    pub fn quasi_steady_check(
        &self,
        time_constant: QtyAny,
        required_ratio: f64,
    ) -> Result<QuasiSteadyVerdict, Vec<Violation>> {
        if self.weighting != DutyWeighting::Dwells {
            return Err(vec![Violation {
                code: "duty-quasi-steady-untimed",
                what: format!(
                    "duty cycle `{}` declares dimensionless fractions, so dwell duration is \
                     unknown and the quasi-steady condition cannot be evaluated",
                    self.name
                ),
                fix: "declare Dwells weighting with durations in seconds. Fractions carry no \
                      absolute time, so no ratio to a time constant exists"
                    .to_string(),
            }]);
        }
        if time_constant.dims != TIME_DIMS
            || !time_constant.value.is_finite()
            || time_constant.value <= 0.0
        {
            return Err(vec![Violation {
                code: "duty-time-constant",
                what: format!(
                    "the supplied time constant ({} with dimensions {:?}) is not a positive \
                     finite duration",
                    time_constant.value, time_constant.dims.0
                ),
                fix: "supply a positive time constant in seconds".to_string(),
            }]);
        }
        if !required_ratio.is_finite() || required_ratio <= 0.0 {
            return Err(vec![Violation {
                code: "duty-required-ratio",
                what: format!("the required dwell ratio {required_ratio} is not positive finite"),
                fix: "supply a positive ratio; five time constants is a common choice for \
                      settling to under one percent"
                    .to_string(),
            }]);
        }
        let Some(first) = self.points.first() else {
            return Err(vec![Violation {
                code: "duty-points-empty",
                what: format!("duty cycle `{}` declares no points", self.name),
                fix: "declare at least one dwell".to_string(),
            }]);
        };

        let mut min_ratio = first.weight.value / time_constant.value;
        let mut shortest = first.name.clone();
        for dwell in &self.points[1..] {
            let ratio = dwell.weight.value / time_constant.value;
            if ratio < min_ratio {
                min_ratio = ratio;
                shortest.clone_from(&dwell.name);
            }
        }

        Ok(QuasiSteadyVerdict {
            min_dwell_ratio: min_ratio,
            shortest_point: shortest,
            required_ratio,
            satisfied: min_ratio >= required_ratio,
        })
    }
}

/// One runnable case compiled from a factored combination.
#[derive(Debug, Clone, PartialEq)]
pub struct CompiledCase {
    /// The combination this case came from.
    pub combination: String,
    /// The resolved `(load case name, factor)` terms — the provenance that
    /// lets a report say *why* this case exists.
    pub terms: Vec<(String, f64)>,
    /// The envelope point the case is evaluated at, when the case set was
    /// crossed with one.
    pub envelope_point: Option<String>,
}

/// A compiled set of runnable cases.
#[derive(Debug, Clone, PartialEq)]
pub struct CaseSet {
    /// The cases, in deterministic order.
    pub cases: Vec<CompiledCase>,
}

/// Which direction of a QoI is worse.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum QoiSense {
    /// A larger value is worse (junction temperature, pressure drop).
    LargerIsWorse,
    /// A smaller value is worse (margin, flow rate).
    SmallerIsWorse,
}

/// Which case sizes the design for one QoI.
#[derive(Debug, Clone, PartialEq)]
pub struct GoverningCase {
    /// The QoI's name.
    pub qoi: String,
    /// Which direction is worse.
    pub sense: QoiSense,
    /// The worst value attained.
    pub value: f64,
    /// Indices into the case set attaining it. More than one means a genuine
    /// tie, reported in full rather than resolved by declaration order.
    pub governing: Vec<usize>,
}

/// Compile factored combinations into runnable cases with provenance.
///
/// [`crate::scenario::Scenario::validate`] remains the authority on whether a
/// declaration is *valid*; this function re-resolves the case references it
/// must resolve in order to build provenance, and reports what it finds using
/// the same violation codes so a caller sees one vocabulary.
///
/// # Errors
/// Returns findings for an unnamed or duplicated combination, a combination
/// with no terms, a term naming an undeclared load case, a repeated case
/// within one combination, or a non-finite factor.
pub fn compile_case_set(
    cases: &[LoadCase],
    combinations: &[Combination],
) -> Result<CaseSet, Vec<Violation>> {
    let mut findings = Vec::new();

    for (index, combination) in combinations.iter().enumerate() {
        if combination.name.trim().is_empty() {
            findings.push(Violation {
                code: "combo-name-empty",
                what: format!("combination {index} has no name"),
                fix: "name every combination; compiled cases cite it as provenance".to_string(),
            });
        }
        if combinations[..index]
            .iter()
            .any(|other| other.name == combination.name)
        {
            findings.push(Violation {
                code: "combo-name-duplicate",
                what: format!(
                    "combination name `{}` is declared more than once",
                    combination.name
                ),
                fix: "combination names must be unique; two cases with one name cannot be told \
                      apart in a governing-case report"
                    .to_string(),
            });
        }
        if combination.terms.is_empty() {
            findings.push(Violation {
                code: "combo-case-empty",
                what: format!("combination `{}` has no terms", combination.name),
                fix: "reference at least one load case, or remove the combination".to_string(),
            });
        }
        for (term_index, (case_name, factor)) in combination.terms.iter().enumerate() {
            if !cases.iter().any(|case| &case.name == case_name) {
                findings.push(Violation {
                    code: "combo-case-missing",
                    what: format!(
                        "combination `{}` references load case `{case_name}`, which is not \
                         declared",
                        combination.name
                    ),
                    fix: "declare the load case, or correct the reference".to_string(),
                });
            }
            if combination.terms[..term_index]
                .iter()
                .any(|(earlier, _)| earlier == case_name)
            {
                findings.push(Violation {
                    code: "combo-term-duplicate",
                    what: format!(
                        "combination `{}` references load case `{case_name}` more than once",
                        combination.name
                    ),
                    fix: "reference each case once with its total factor. Summing repeated \
                          references silently would make the effective factor invisible at the \
                          declaration"
                        .to_string(),
                });
            }
            if !factor.is_finite() {
                findings.push(Violation {
                    code: "combo-factor",
                    what: format!(
                        "combination `{}` gives load case `{case_name}` a non-finite factor \
                         ({factor})",
                        combination.name
                    ),
                    fix: "supply a finite factor".to_string(),
                });
            }
        }
    }

    if !findings.is_empty() {
        return Err(findings);
    }

    Ok(CaseSet {
        cases: combinations
            .iter()
            .map(|combination| CompiledCase {
                combination: combination.name.clone(),
                terms: combination.terms.clone(),
                envelope_point: None,
            })
            .collect(),
    })
}

impl CaseSet {
    /// Cross every case with every named envelope point.
    ///
    /// # Errors
    /// Returns an `envelope-case-overflow` finding when the product exceeds
    /// `budget.max_cases`, checked before allocation, and a
    /// `envelope-case-point-duplicate` finding when two points share a name.
    pub fn cross_with_points(
        &self,
        points: &[(String, EnvelopePoint)],
        budget: EnvelopeBudget,
    ) -> Result<CaseSet, Vec<Violation>> {
        let mut findings = Vec::new();
        for (index, (name, _)) in points.iter().enumerate() {
            if points[..index].iter().any(|(earlier, _)| earlier == name) {
                findings.push(Violation {
                    code: "envelope-case-point-duplicate",
                    what: format!("envelope point name `{name}` appears more than once"),
                    fix: "give each point a unique name; the name is a compiled case's only \
                          record of where it was evaluated"
                        .to_string(),
                });
            }
        }
        if !findings.is_empty() {
            return Err(findings);
        }

        let total = self.cases.len().checked_mul(points.len());
        match total {
            Some(count) if count <= budget.max_cases => {}
            other => {
                return Err(vec![Violation {
                    code: "envelope-case-overflow",
                    what: format!(
                        "{} combination(s) x {} envelope point(s) is {} case(s), above the \
                         admitted maximum of {}",
                        self.cases.len(),
                        points.len(),
                        other.map_or_else(
                            || "an overflowing number of".to_string(),
                            |count| count.to_string()
                        ),
                        budget.max_cases
                    ),
                    fix: "reduce the point set (corner enumeration rather than a full sweep), \
                          reduce the combinations, or raise max_cases on an explicit \
                          EnvelopeBudget. The product is refused before allocation"
                        .to_string(),
                }]);
            }
        }

        let mut cases = Vec::with_capacity(self.cases.len() * points.len());
        for case in &self.cases {
            for (name, _) in points {
                cases.push(CompiledCase {
                    combination: case.combination.clone(),
                    terms: case.terms.clone(),
                    envelope_point: Some(name.clone()),
                });
            }
        }
        Ok(CaseSet { cases })
    }
}

/// Identify which case governs one QoI.
///
/// `values` supplies the QoI at each case, in case-set order. A tie reports
/// every case that attains the extremum: a governing-case report that silently
/// picked the first would hide that two conditions size the design equally.
///
/// # Errors
/// Returns findings when `values` is empty or contains a non-finite entry.
pub fn governing_case(
    qoi: &str,
    sense: QoiSense,
    values: &[f64],
) -> Result<GoverningCase, Vec<Violation>> {
    if values.is_empty() {
        return Err(vec![Violation {
            code: "governing-no-cases",
            what: format!("no case values were supplied for QoI `{qoi}`"),
            fix: "evaluate at least one case before asking which governs".to_string(),
        }]);
    }
    for (index, value) in values.iter().enumerate() {
        if !value.is_finite() {
            return Err(vec![Violation {
                code: "governing-nonfinite",
                what: format!("QoI `{qoi}` at case {index} is non-finite ({value})"),
                fix: "a case that failed to solve must be reported as a refusal; treating it as a \
                      value would let a NaN silently never govern"
                    .to_string(),
            }]);
        }
    }

    let worse = |candidate: f64, incumbent: f64| match sense {
        QoiSense::LargerIsWorse => candidate > incumbent,
        QoiSense::SmallerIsWorse => candidate < incumbent,
    };

    let mut extremum = values[0];
    for value in &values[1..] {
        if worse(*value, extremum) {
            extremum = *value;
        }
    }
    // Exact equality again, deliberately: `extremum` was SELECTED from this
    // same slice, so every case that attains it is bit-identical to it. A
    // tolerance would report near-governing cases as governing, which is a
    // different — and unasked — engineering question.
    #[allow(clippy::float_cmp)]
    let governing = values
        .iter()
        .enumerate()
        .filter(|(_, value)| **value == extremum)
        .map(|(index, _)| index)
        .collect();

    Ok(GoverningCase {
        qoi: qoi.to_string(),
        sense,
        value: extremum,
        governing,
    })
}

/// Build a bounded fan-state axis: nominal plus one single-unit failure each.
///
/// `n` units give `n + 1` states, not `2^n`. Exhaustive multi-unit failure
/// enumeration is deliberately not offered: the count is exponential in a
/// quantity nobody bounds at declaration time, and the failure combinations
/// worth analysing are a judgement about redundancy, not a product. Declare
/// the multi-failure states you mean, explicitly.
///
/// # Errors
/// Returns findings for an empty unit list, an unnamed unit, or a duplicated
/// unit name.
pub fn single_failure_axis(
    axis_name: impl Into<String>,
    nominal_state: impl Into<String>,
    units: &[&str],
) -> Result<EnvelopeAxis, Vec<Violation>> {
    let axis_name = axis_name.into();
    let mut findings = Vec::new();

    if units.is_empty() {
        findings.push(Violation {
            code: "failure-units-empty",
            what: format!("axis `{axis_name}` was asked for single-failure states over no units"),
            fix: "name the units that can fail".to_string(),
        });
    }
    for (index, unit) in units.iter().enumerate() {
        if unit.trim().is_empty() {
            findings.push(Violation {
                code: "failure-unit-name-empty",
                what: format!("unit {index} of axis `{axis_name}` has no name"),
                fix: "name every unit; the failure state is named after it".to_string(),
            });
        }
        if units[..index].contains(unit) {
            findings.push(Violation {
                code: "failure-unit-duplicate",
                what: format!("unit `{unit}` of axis `{axis_name}` is named more than once"),
                fix: "unit names must be unique; two failure states would otherwise collide"
                    .to_string(),
            });
        }
    }
    if !findings.is_empty() {
        return Err(findings);
    }

    let mut states = vec![DiscreteState::nominal(nominal_state)];
    for unit in units {
        states.push(DiscreteState::failed(format!("{unit}-failed")));
    }
    Ok(EnvelopeAxis::discrete(axis_name, states))
}

/// A worked reference envelope: ambient temperature, total dissipated power
/// with a declared thermal-throttle limit against ambient, and a two-fan axis
/// with single-failure states.
///
/// Provided as a named constructor rather than a test fixture because it is
/// the envelope the CONTRACT's worked example refers to. It is *representative,
/// not canonical*: no repo-wide "reference cooling scenario" artifact exists.
#[must_use]
pub fn reference_cooling_envelope() -> OperatingEnvelope {
    let fans = single_failure_axis("fan-state", "both-running", &["fan-a", "fan-b"])
        .expect("two distinct named units");
    OperatingEnvelope {
        name: "reference-cooling".to_string(),
        axes: vec![
            EnvelopeAxis::continuous(
                "ambient-temperature",
                QtyAny::new(300.0, TEMP_DIMS),
                QtyAny::new(320.0, TEMP_DIMS),
            ),
            EnvelopeAxis::continuous(
                "total-power",
                QtyAny::new(100.0, POWER_DIMS),
                QtyAny::new(200.0, POWER_DIMS),
            ),
            fans,
        ],
        couplings: vec![AxisCoupling {
            a: "ambient-temperature".to_string(),
            b: "total-power".to_string(),
            relation: CouplingRelation::CoOccurrenceLimit {
                a_coeff: QtyAny::new(2.5, POWER_DIMS.checked_minus(TEMP_DIMS).expect("in range")),
                b_coeff: QtyAny::dimensionless(1.0),
                bound: QtyAny::new(975.0, POWER_DIMS),
                rationale: "the controller throttles dissipation as ambient rises, so full power \
                            at the top of the ambient range is not attainable"
                    .to_string(),
            },
        }],
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_box_with_no_limit_keeps_all_four_corners() {
        let vertices = clip_box_to_limit((0.0, 1.0), (0.0, 1.0), 1.0, 1.0, 10.0);
        assert_eq!(vertices.len(), 4);
    }

    #[test]
    fn clipping_a_corner_off_a_box_adds_a_vertex() {
        // 2.5 a + b <= 975 over [300,320] x [100,200] cuts the (320,200)
        // corner and introduces two: 4 - 1 + 2 = 5.
        let vertices = clip_box_to_limit((300.0, 320.0), (100.0, 200.0), 2.5, 1.0, 975.0);
        assert_eq!(
            vertices,
            vec![
                (300.0, 100.0),
                (300.0, 200.0),
                (310.0, 200.0),
                (320.0, 100.0),
                (320.0, 175.0),
            ]
        );
    }

    #[test]
    fn a_limit_through_a_corner_does_not_duplicate_it() {
        // 2.5 a + b <= 950 passes exactly through (300, 200).
        let vertices = clip_box_to_limit((300.0, 320.0), (100.0, 200.0), 2.5, 1.0, 950.0);
        assert_eq!(
            vertices,
            vec![
                (300.0, 100.0),
                (300.0, 200.0),
                (320.0, 100.0),
                (320.0, 150.0),
            ]
        );
    }

    #[test]
    fn the_reference_envelope_validates() {
        assert_eq!(reference_cooling_envelope().validate(), Vec::new());
    }
}
