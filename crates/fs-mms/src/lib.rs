//! G1 harness (bead frankensim-epic-gauntlet-6nb.2): manufactured-solution
//! refinement ladders, convergence-order fitting, and THE 0.2 slope gate —
//! the tier that catches the quiet death of accuracy.
//!
//! A discretization runs a refinement ladder against a manufactured exact
//! solution and reports (h, error-norm) points; [`fit_order`] takes the
//! least-squares slope of log(error) versus log(h) through
//! `fs_math::det::ln` (cross-ISA deterministic bits), and [`OrderGate`]
//! fails the build when the observed order deviates from the theoretical
//! order by more than 0.2 — primal and ADJOINT ladders get the identical
//! treatment (dual consistency is verified, not assumed; the gate does
//! not care which side the errors came from, the case record does).
//!
//! The battery matrix ([`MmsMatrix`]) is DECLARED IN DATA: every
//! (frontend, element family, boundary condition) row is either covered
//! by a named test or an explicit gap with a reason, so coverage holes
//! are visible and lintable instead of silently absent.
//!
//! No-claims: the harness fits and gates orders it is HANDED — it does
//! not discretize, does not generate forcing terms from symbolic
//! solutions (that is the fs-opdsl integration, recorded as a matrix gap
//! until fs-opdsl lands), and does not ledger results (fingerprint-bound
//! persistence is fs-obs/fs-ledger scope; records here are JSON lines).

use core::fmt::{self, Write as _};

/// The bead's gate half-width: observed order may deviate from the
/// theoretical order by at most this much.
pub const ORDER_GATE_TOLERANCE: f64 = 0.2;

fn escape_json_string(value: &str) -> String {
    let mut escaped = String::with_capacity(value.len());
    for ch in value.chars() {
        match ch {
            '"' => escaped.push_str("\\\""),
            '\\' => escaped.push_str("\\\\"),
            '\u{0008}' => escaped.push_str("\\b"),
            '\u{000c}' => escaped.push_str("\\f"),
            '\n' => escaped.push_str("\\n"),
            '\r' => escaped.push_str("\\r"),
            '\t' => escaped.push_str("\\t"),
            '\u{0000}'..='\u{001f}' => {
                write!(&mut escaped, "\\u{:04x}", u32::from(ch))
                    .expect("writing to a String cannot fail");
            }
            _ => escaped.push(ch),
        }
    }
    escaped
}

fn json_fixed_number(value: f64) -> String {
    if value.is_finite() {
        format!("{value:.6}")
    } else {
        "null".to_string()
    }
}

fn json_scientific_number(value: f64) -> String {
    if value.is_finite() {
        format!("{value:.3e}")
    } else {
        "null".to_string()
    }
}

/// One refinement ladder: mesh parameters and matching error norms.
#[derive(Debug, Clone, PartialEq)]
pub struct RefinementLadder {
    hs: Vec<f64>,
    errors: Vec<f64>,
}

/// A typed refusal from ladder admission or order gating.
#[derive(Debug, Clone, PartialEq)]
pub struct MmsError {
    rule: &'static str,
    detail: String,
}

impl MmsError {
    fn new(rule: &'static str, detail: impl Into<String>) -> MmsError {
        MmsError {
            rule,
            detail: detail.into(),
        }
    }

    /// Stable machine-readable rule slug.
    #[must_use]
    pub const fn rule(&self) -> &'static str {
        self.rule
    }

    /// Human-readable detail.
    #[must_use]
    pub fn detail(&self) -> &str {
        &self.detail
    }
}

impl fmt::Display for MmsError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}: {}", self.rule, self.detail)
    }
}

impl std::error::Error for MmsError {}

impl RefinementLadder {
    /// Admit a ladder: at least three rungs (a two-point slope hides
    /// pre-asymptotic bends), strictly decreasing positive `h`, strictly
    /// positive finite errors.
    pub fn new(hs: Vec<f64>, errors: Vec<f64>) -> Result<RefinementLadder, MmsError> {
        if hs.len() != errors.len() {
            return Err(MmsError::new(
                "mms-ladder-shape",
                format!("{} h values vs {} errors", hs.len(), errors.len()),
            ));
        }
        if hs.len() < 3 {
            return Err(MmsError::new(
                "mms-ladder-shape",
                "a ladder needs at least three rungs for an honest slope",
            ));
        }
        for window in hs.windows(2) {
            if !(window[1] < window[0]) {
                return Err(MmsError::new(
                    "mms-ladder-order",
                    "h must be strictly decreasing",
                ));
            }
        }
        for &value in hs.iter().chain(errors.iter()) {
            if !(value.is_finite() && value > 0.0) {
                return Err(MmsError::new(
                    "mms-ladder-domain",
                    format!("ladder values must be finite and positive, got {value:e}"),
                ));
            }
        }
        Ok(RefinementLadder { hs, errors })
    }

    /// The mesh parameters.
    #[must_use]
    pub fn hs(&self) -> &[f64] {
        &self.hs
    }

    /// The error norms.
    #[must_use]
    pub fn errors(&self) -> &[f64] {
        &self.errors
    }
}

/// A least-squares convergence-order fit.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct OrderFit {
    /// The observed order: the log-log slope.
    pub observed: f64,
    /// The log-space intercept (log error at h = 1).
    pub intercept: f64,
    /// Root-mean-square log-space residual of the fit (envelope of how
    /// straight the ladder actually is).
    pub rms_residual: f64,
}

/// Least-squares slope of log(error) versus log(h). Logs route through
/// `fs_math::det::ln`, so the fit is bit-identical across ISAs.
#[must_use]
pub fn fit_order(ladder: &RefinementLadder) -> OrderFit {
    let n = ladder.hs.len() as f64;
    let xs: Vec<f64> = ladder.hs.iter().map(|&h| fs_math::det::ln(h)).collect();
    let ys: Vec<f64> = ladder.errors.iter().map(|&e| fs_math::det::ln(e)).collect();
    let sx: f64 = xs.iter().sum();
    let sy: f64 = ys.iter().sum();
    let sxx: f64 = xs.iter().map(|x| x * x).sum();
    let sxy: f64 = xs.iter().zip(&ys).map(|(x, y)| x * y).sum();
    let denom = n * sxx - sx * sx;
    let observed = (n * sxy - sx * sy) / denom;
    let intercept = (sy - observed * sx) / n;
    let mut ss = 0.0f64;
    for (x, y) in xs.iter().zip(&ys) {
        let r = y - (observed * x + intercept);
        ss += r * r;
    }
    OrderFit {
        observed,
        intercept,
        rms_residual: (ss / n).sqrt(),
    }
}

/// Which side of the primal/dual pair a ladder came from (the gate is
/// identical; the record keeps the distinction visible).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum LadderSide {
    /// Primal solution error.
    Primal,
    /// Adjoint/dual solution error (dual consistency verified, not
    /// assumed).
    Adjoint,
}

impl LadderSide {
    const fn name(self) -> &'static str {
        match self {
            Self::Primal => "primal",
            Self::Adjoint => "adjoint",
        }
    }
}

/// THE gate: observed order must sit within [`ORDER_GATE_TOLERANCE`] of
/// the theoretical order.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct OrderGate {
    /// The scheme's theoretical order.
    pub theoretical: f64,
}

/// One gated convergence verdict (a JSON-lines record).
#[derive(Debug, Clone, PartialEq)]
pub struct OrderVerdict {
    /// The case name.
    pub case: String,
    /// Primal or adjoint.
    pub side: LadderSide,
    /// The fit.
    pub fit: OrderFit,
    /// The theoretical order gated against.
    pub theoretical: f64,
    /// |observed − theoretical|.
    pub deviation: f64,
}

impl OrderVerdict {
    fn numeric_fields_are_finite(&self) -> bool {
        self.fit.observed.is_finite()
            && self.fit.intercept.is_finite()
            && self.fit.rms_residual.is_finite()
            && self.theoretical.is_finite()
            && self.deviation.is_finite()
    }

    /// The verdict as one JSON line.
    #[must_use]
    pub fn json_line(&self, pass: bool) -> String {
        let pass = pass && self.numeric_fields_are_finite();
        format!(
            "{{\"mms\":\"order\",\"case\":\"{}\",\"side\":\"{}\",\"observed\":{},\
             \"theoretical\":{},\"deviation\":{},\"rms_residual\":{},\"pass\":{pass}}}",
            escape_json_string(&self.case),
            self.side.name(),
            json_fixed_number(self.fit.observed),
            json_fixed_number(self.theoretical),
            json_fixed_number(self.deviation),
            json_scientific_number(self.fit.rms_residual)
        )
    }
}

impl OrderGate {
    /// Gate a fitted ladder: build-fails (returns the refusal) when the
    /// observed order deviates from the theoretical order by more than
    /// [`ORDER_GATE_TOLERANCE`]. The passing verdict is returned for the
    /// caller's structured log.
    pub fn check(
        &self,
        case: &str,
        side: LadderSide,
        ladder: &RefinementLadder,
    ) -> Result<OrderVerdict, MmsError> {
        if !(self.theoretical.is_finite() && self.theoretical > 0.0) {
            return Err(MmsError::new(
                "mms-gate-domain",
                "theoretical order must be finite and positive",
            ));
        }
        let fit = fit_order(ladder);
        let deviation = (fit.observed - self.theoretical).abs();
        let verdict = OrderVerdict {
            case: case.to_owned(),
            side,
            fit,
            theoretical: self.theoretical,
            deviation,
        };
        if !verdict.numeric_fields_are_finite() {
            return Err(MmsError::new("mms-gate-domain", verdict.json_line(false)));
        }
        if deviation > ORDER_GATE_TOLERANCE {
            return Err(MmsError::new("mms-order-gate", verdict.json_line(false)));
        }
        Ok(verdict)
    }
}

/// Coverage status of one battery-matrix row.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Coverage {
    /// Covered by a named test.
    Covered {
        /// The test path/name carrying the MMS battery.
        test: String,
    },
    /// An explicit, reasoned gap — visible, never silent.
    Gap {
        /// Why the row is not covered yet.
        reason: String,
    },
}

/// One declared battery-matrix row: frontend × element family × BC type.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MmsMatrixRow {
    /// The PDE frontend (e.g. "feec-body-fitted", "cutfem-sdf",
    /// "iga-patch").
    pub frontend: String,
    /// The element family (e.g. "p1-simplicial", "cut-p1", "nurbs-p2").
    pub family: String,
    /// The boundary-condition type (e.g. "dirichlet", "neumann",
    /// "mortar-seam", "sliver-cut").
    pub bc: String,
    /// Covered-by-test or an explicit gap.
    pub coverage: Coverage,
}

/// The declared MMS battery matrix: coverage in data, gaps lintable.
#[derive(Debug, Clone, Default)]
pub struct MmsMatrix {
    /// The declared rows.
    pub rows: Vec<MmsMatrixRow>,
}

impl MmsMatrix {
    /// The explicit gaps (the lint output: every hole has a reason).
    #[must_use]
    pub fn gaps(&self) -> Vec<&MmsMatrixRow> {
        self.rows
            .iter()
            .filter(|r| matches!(r.coverage, Coverage::Gap { .. }))
            .collect()
    }

    /// One JSON line per row, coverage visible.
    #[must_use]
    pub fn json_lines(&self) -> Vec<String> {
        self.rows
            .iter()
            .map(|r| {
                let (status, detail) = match &r.coverage {
                    Coverage::Covered { test } => ("covered", test.as_str()),
                    Coverage::Gap { reason } => ("gap", reason.as_str()),
                };
                format!(
                    "{{\"mms\":\"matrix\",\"frontend\":\"{}\",\"family\":\"{}\",\"bc\":\"{}\",\
                     \"status\":\"{status}\",\"detail\":\"{}\"}}",
                    escape_json_string(&r.frontend),
                    escape_json_string(&r.family),
                    escape_json_string(&r.bc),
                    escape_json_string(detail)
                )
            })
            .collect()
    }
}
