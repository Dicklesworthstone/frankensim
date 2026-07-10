//! The THREE-COLOR epistemic schema (Proposal 3): every ledger quantity
//! carries a COLOR — **verified** (interval-certified numerics),
//! **validated** (anchored to experimental data in a stated REGIME), or
//! **estimated** (cross-model probes, surrogates) — and the composition
//! algebra is CONSERVATIVE BY TYPE: an estimate can never be laundered
//! into a certificate, and validation is a REGIONAL property that
//! auto-demotes the moment execution leaves its regime.
//!
//! Why: a pipeline can be interval-certified end-to-end and still be
//! precisely WRONG, because model-form error does not compose
//! algebraically. The colors keep the boundary between numerical
//! certainty and modeling uncertainty visible in the type system.
//!
//! Layer discipline (bead qmao.1 polish notes): the enum, payloads, and
//! pairwise algebra live HERE in fs-evidence so any layer can color a
//! `Certified<T>` without touching HELM; write-time enforcement lives
//! in fs-ledger over already-colored values.

use crate::{ModelEvidence, NumericalCertificate, NumericalKind, ValidityDomain};
use std::collections::BTreeMap;

/// The epistemic color, with its color-specific payload.
#[derive(Debug, Clone, PartialEq)]
pub enum Color {
    /// Interval-certified numerics: the payload IS the bound.
    Verified {
        /// Certified lower bound.
        lo: f64,
        /// Certified upper bound.
        hi: f64,
    },
    /// Anchored to experimental data INSIDE a stated regime.
    Validated {
        /// The region of feature space where the anchoring holds
        /// (Reynolds range, strain range, …).
        regime: ValidityDomain,
        /// Identity of the anchoring dataset.
        dataset: String,
    },
    /// Cross-model discrepancy probes, surrogates, heuristics.
    Estimated {
        /// The estimator's identity.
        estimator: String,
        /// The estimator's own dispersion (∞ = no spread claim).
        dispersion: f64,
    },
}

/// Lattice rank (higher = stronger claim). The composition result can
/// never OUTRANK the weakest operand — that is the no-laundering law.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum ColorRank {
    /// Estimated.
    Estimated,
    /// Validated.
    Validated,
    /// Verified.
    Verified,
}

impl Color {
    /// The lattice rank.
    #[must_use]
    pub fn rank(&self) -> ColorRank {
        match self {
            Color::Verified { .. } => ColorRank::Verified,
            Color::Validated { .. } => ColorRank::Validated,
            Color::Estimated { .. } => ColorRank::Estimated,
        }
    }

    /// Stable name for ledger rows.
    #[must_use]
    pub fn name(&self) -> &'static str {
        match self {
            Color::Verified { .. } => "verified",
            Color::Validated { .. } => "validated",
            Color::Estimated { .. } => "estimated",
        }
    }

    /// Canonical JSON payload for ledger rows.
    #[must_use]
    pub fn payload_json(&self) -> String {
        match self {
            Color::Verified { lo, hi } => {
                format!("{{\"interval\":[{},{}]}}", json_f64(*lo), json_f64(*hi))
            }
            Color::Validated { regime, dataset } => {
                use core::fmt::Write as _;
                let mut axes = String::new();
                for (k, (lo, hi)) in regime.bounds() {
                    if !axes.is_empty() {
                        axes.push(',');
                    }
                    let _ = write!(
                        axes,
                        "{}:[{},{}]",
                        json_string(k),
                        json_f64(*lo),
                        json_f64(*hi)
                    );
                }
                format!(
                    "{{\"dataset\":{},\"regime\":{{{axes}}}}}",
                    json_string(dataset)
                )
            }
            Color::Estimated {
                estimator,
                dispersion,
            } => format!(
                "{{\"estimator\":{},\"dispersion\":{}}}",
                json_string(estimator),
                json_f64(*dispersion)
            ),
        }
    }
}

fn json_f64(value: f64) -> String {
    if value.is_finite() {
        format!("{value:.6e}")
    } else {
        format!("\"non-finite:{value}\"")
    }
}

fn json_string(value: &str) -> String {
    use core::fmt::Write as _;
    let mut out = String::with_capacity(value.len() + 2);
    out.push('"');
    for ch in value.chars() {
        match ch {
            '"' => out.push_str("\\\""),
            '\\' => out.push_str("\\\\"),
            '\n' => out.push_str("\\n"),
            '\r' => out.push_str("\\r"),
            '\t' => out.push_str("\\t"),
            c if u32::from(c) < 0x20 => {
                let _ = write!(out, "\\u{:04x}", u32::from(c));
            }
            c => out.push(c),
        }
    }
    out.push('"');
    out
}

/// How verified intervals combine under the ledger operation.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum IntervalOp {
    /// Sum of quantities: bounds add.
    Add,
    /// Product: the four-corner hull.
    Mul,
    /// Anything else: the conservative hull of both bounds.
    Hull,
}

/// A demotion record (validated → estimated on regime exit).
#[derive(Debug, Clone, PartialEq)]
pub struct Demotion {
    /// The dataset whose regime was exited.
    pub dataset: String,
    /// The axis that left its range.
    pub axis: String,
    /// The offending state value.
    pub value: f64,
}

/// Teaching errors for the color algebra.
#[derive(Debug, Clone, PartialEq)]
pub enum ColorError {
    /// A `Verified` claim was attempted from a non-enclosure
    /// certificate — the laundering refusal.
    LaunderingRefused {
        /// What the certificate actually was.
        actual: &'static str,
    },
    /// A regime axis referenced by the state is missing from the tag.
    IncompleteState {
        /// The axis the regime declares but the state omits.
        axis: String,
    },
}

impl core::fmt::Display for ColorError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            ColorError::LaunderingRefused { actual } => write!(
                f,
                "refusing to mark this quantity `verified`: its certificate is \
                 `{actual}`, not an enclosure — estimates cannot be laundered into \
                 certificates (attach a signed waiver to the LEDGER write if a \
                 human accepts responsibility, and the waiver will travel in \
                 provenance)"
            ),
            ColorError::IncompleteState { axis } => write!(
                f,
                "the execution state does not report regime axis `{axis}`; a \
                 validated claim cannot be checked against an unreported axis \
                 (report it, or the value demotes)"
            ),
        }
    }
}

impl std::error::Error for ColorError {}

/// Color a quantity `verified` FROM an interval-certified numerical
/// certificate — the only door in. Anything weaker refuses.
///
/// # Errors
/// [`ColorError::LaunderingRefused`] for estimate/no-claim kinds.
pub fn verified_from(cert: &NumericalCertificate) -> Result<Color, ColorError> {
    match cert.kind {
        NumericalKind::Exact | NumericalKind::Enclosure => Ok(Color::Verified {
            lo: cert.lo,
            hi: cert.hi,
        }),
        NumericalKind::Estimate => Err(ColorError::LaunderingRefused { actual: "estimate" }),
        NumericalKind::NoClaim => Err(ColorError::LaunderingRefused { actual: "no-claim" }),
    }
}

/// Derive the honest color of an existing [`crate::Evidence`] receipt:
/// enclosure-grade numerics → verified; model evidence with a bounded
/// validity domain → validated (regime = that domain); otherwise →
/// estimated with the model/card identity.
#[must_use]
pub fn color_of(numerical: &NumericalCertificate, model: &ModelEvidence) -> Color {
    if model.cards.is_empty()
        && matches!(
            cert_kind(numerical),
            NumericalKind::Exact | NumericalKind::Enclosure
        )
    {
        return Color::Verified {
            lo: numerical.lo,
            hi: numerical.hi,
        };
    }
    if !model.cards.is_empty() && !model.validity.bounds().is_empty() {
        return Color::Validated {
            regime: model.validity.clone(),
            dataset: model.cards.join("+"),
        };
    }
    Color::Estimated {
        estimator: if model.cards.is_empty() {
            "uncarded-numerics".to_string()
        } else {
            model.cards.join("+")
        },
        dispersion: model.discrepancy_rel.abs(),
    }
}

fn cert_kind(c: &NumericalCertificate) -> NumericalKind {
    c.kind
}

/// Outward-round an arithmetic interval endpoint pair (one ulp each way) so a
/// composed enclosure stays a TRUE enclosure under round-to-nearest — otherwise
/// the lower bound could round UP (or the upper bound DOWN) and EXCLUDE the
/// exact result, laundering the very precision a `Verified` color certifies. A
/// NaN endpoint (∞ + −∞) degrades to the whole real line, the only sound
/// enclosure. Mirrors `NumericalCertificate::compose` in `lib.rs`.
fn outward_round(lo: f64, hi: f64) -> (f64, f64) {
    if lo.is_nan() || hi.is_nan() {
        (f64::NEG_INFINITY, f64::INFINITY)
    } else {
        (lo.next_down(), hi.next_up())
    }
}

/// The TOTAL, conservative pairwise composition: the result's rank is
/// the MINIMUM of the operands' ranks (no laundering), verified bounds
/// combine per `op` (outward-rounded to preserve enclosure), validated
/// regimes INTERSECT (both anchors must hold), and estimated dispersions
/// add (conservative).
#[must_use]
pub fn compose(a: &Color, b: &Color, op: IntervalOp) -> Color {
    match (a, b) {
        (Color::Verified { lo: a0, hi: a1 }, Color::Verified { lo: b0, hi: b1 }) => {
            let (lo, hi) = match op {
                IntervalOp::Add => outward_round(a0 + b0, a1 + b1),
                IntervalOp::Mul => {
                    let c = [a0 * b0, a0 * b1, a1 * b0, a1 * b1];
                    if c.iter().any(|x| x.is_nan()) {
                        // A NaN corner (0×∞) is indeterminate AND the min/max
                        // fold below would SILENTLY DROP it (f64::min/max ignore
                        // NaN), reporting a bogus tight interval. The whole real
                        // line is the only sound enclosure.
                        (f64::NEG_INFINITY, f64::INFINITY)
                    } else {
                        outward_round(
                            c.iter().copied().fold(f64::INFINITY, f64::min),
                            c.iter().copied().fold(f64::NEG_INFINITY, f64::max),
                        )
                    }
                }
                // Hull is exact endpoint selection (no arithmetic → no rounding).
                IntervalOp::Hull => (a0.min(*b0), a1.max(*b1)),
            };
            Color::Verified { lo, hi }
        }
        // Validated ⊕ verified stays validated (the weaker anchor);
        // validated ⊕ validated intersects regimes.
        (Color::Validated { regime, dataset }, Color::Verified { .. })
        | (Color::Verified { .. }, Color::Validated { regime, dataset }) => Color::Validated {
            regime: regime.clone(),
            dataset: dataset.clone(),
        },
        (
            Color::Validated {
                regime: r1,
                dataset: d1,
            },
            Color::Validated {
                regime: r2,
                dataset: d2,
            },
        ) => {
            let regime = intersect_domains(r1, r2);
            if regime.is_empty() {
                // Mutually unsatisfiable regimes: there is NO state where both
                // anchors hold, so the composition cannot honestly stay
                // Validated. There is also no defensible spread claim, so the
                // demotion carries infinite dispersion just like a regime exit.
                Color::Estimated {
                    estimator: format!("disjoint-regimes:{d1}&{d2}"),
                    dispersion: f64::INFINITY,
                }
            } else {
                Color::Validated {
                    regime,
                    dataset: if d1 == d2 {
                        d1.clone()
                    } else {
                        format!("{d1}&{d2}")
                    },
                }
            }
        }
        // Anything ⊕ estimated → estimated. No exceptions here; the
        // waiver door lives at the LEDGER, in provenance.
        (
            Color::Estimated {
                estimator,
                dispersion,
            },
            other,
        )
        | (
            other,
            Color::Estimated {
                estimator,
                dispersion,
            },
        ) => {
            let (other_disp, other_id) = match other {
                Color::Estimated {
                    estimator: e2,
                    dispersion: v,
                } => (*v, e2.clone()),
                Color::Verified { .. } => (0.0, "verified".to_string()),
                Color::Validated { dataset, .. } => (0.0, format!("validated:{dataset}")),
            };
            Color::Estimated {
                estimator: if other_id == "verified" {
                    estimator.clone()
                } else {
                    format!("{estimator}+{other_id}")
                },
                dispersion: dispersion + other_disp,
            }
        }
    }
}

/// Intersection of two validity domains (axis-wise; an axis present in
/// either constrains the result — both anchors must hold). Delegates to
/// [`ValidityDomain::intersect`], which PRESERVES emptiness (`lo > hi`) on a
/// disjoint axis so callers can detect an unsatisfiable regime via
/// [`ValidityDomain::is_empty`]. (A previous `hi.max(lo)` clamp here collapsed a
/// disjoint intersection into a phantom single point, silently claiming
/// validity at a state where neither anchor holds.)
#[must_use]
pub fn intersect_domains(a: &ValidityDomain, b: &ValidityDomain) -> ValidityDomain {
    a.intersect(b)
}

/// Check a validated color against the CURRENT execution state:
/// inside the regime → unchanged; outside (or unreported axis) →
/// AUTOMATIC DEMOTION to estimated, with the flag returned. Verified
/// and estimated colors pass through untouched.
#[must_use]
pub fn check_regime(color: &Color, state: &BTreeMap<String, f64>) -> (Color, Option<Demotion>) {
    let Color::Validated { regime, dataset } = color else {
        return (color.clone(), None);
    };
    if regime.bounds().is_empty() {
        return demote(dataset, "<undeclared-regime>", f64::NAN);
    }
    for (axis, &(lo, hi)) in regime.bounds() {
        let Some(&v) = state.get(axis) else {
            return demote(dataset, axis, f64::NAN);
        };
        if !lo.is_finite() || !hi.is_finite() || lo > hi || !v.is_finite() || v < lo || v > hi {
            return demote(dataset, axis, v);
        }
    }
    (color.clone(), None)
}

fn demote(dataset: &str, axis: &str, value: f64) -> (Color, Option<Demotion>) {
    (
        Color::Estimated {
            estimator: format!("regime-exit:{dataset}@{axis}"),
            dispersion: f64::INFINITY,
        },
        Some(Demotion {
            dataset: dataset.to_string(),
            axis: axis.to_string(),
            value,
        }),
    )
}
