//! The Error Ledger and the Time Ledger (plan §11.4, P4): end-to-end
//! attribution trees that make "how accurate is this number and where did
//! the error come from" — and "where did the seconds go" — QUERIES.
//!
//! Honesty model: the Error Ledger is ATTRIBUTION BOOKKEEPING over
//! contributions that are themselves estimates or certificates (their
//! rigor class is carried per entry); rigorous enclosure composition lives
//! in fs-evidence/fs-ivl. Composition here is first-order additive —
//! conservative for the error sources the plan names (they add in the
//! worst case), with the conservativeness LAW checked on fixtures where
//! stage errors are known bounds. Any unattributed mass must be declared
//! as residual: `total() = Σ contributions + declared_residual`, and the
//! completeness lint refuses trees whose declared residual is negative or
//! NaN — no silent error mass.

use core::fmt;

/// Where a piece of QoI error came from (the plan's canonical sources).
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum ErrorSource {
    /// Geometry/chart tolerance (conversion receipts).
    Geometry,
    /// Discretization (mesh/grid/order truncation).
    Discretization,
    /// Algebraic residual (solver stopping tolerance).
    Algebraic,
    /// Surrogate band (conformal/ROM error).
    Surrogate,
    /// Statistical noise (MC/MLMC half-widths).
    Statistical,
    /// Model-form discrepancy (closures, constitutive laws).
    ModelForm,
}

impl ErrorSource {
    /// Stable lowercase name.
    #[must_use]
    pub fn name(self) -> &'static str {
        match self {
            ErrorSource::Geometry => "geometry",
            ErrorSource::Discretization => "discretization",
            ErrorSource::Algebraic => "algebraic",
            ErrorSource::Surrogate => "surrogate",
            ErrorSource::Statistical => "statistical",
            ErrorSource::ModelForm => "model_form",
        }
    }
}

/// How trustworthy a contribution's magnitude is.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Rigor {
    /// Certificate-backed bound (interval/conformal/e-process).
    Certified,
    /// A-posteriori estimate (DWR, residual-based).
    Estimated,
    /// A-priori rate model (order-based extrapolation).
    RateModel,
}

/// One attributed error contribution.
#[derive(Debug, Clone, PartialEq)]
pub struct Contribution {
    /// The source class.
    pub source: ErrorSource,
    /// Which operator/stage produced it.
    pub label: String,
    /// Absolute contribution to the QoI error (≥ 0).
    pub abs: f64,
    /// The magnitude's trust class.
    pub rigor: Rigor,
}

/// A study's error attribution tree.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct ErrorLedger {
    /// Every attributed contribution, in pipeline order.
    pub contributions: Vec<Contribution>,
    /// Declared unattributed mass (≥ 0; "we know we don't know this much").
    pub declared_residual: f64,
}

/// A structured attribution defect.
#[derive(Debug, Clone, PartialEq)]
pub enum LedgerDefect {
    /// A contribution is negative or non-finite.
    BadContribution {
        /// The offending label.
        label: String,
    },
    /// The declared residual is negative or non-finite.
    BadResidual,
}

impl fmt::Display for LedgerDefect {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            LedgerDefect::BadContribution { label } => write!(
                f,
                "contribution {label:?} is negative or non-finite — error mass must be \
                 a nonnegative magnitude"
            ),
            LedgerDefect::BadResidual => {
                write!(f, "declared residual must be nonnegative and finite")
            }
        }
    }
}

impl ErrorLedger {
    /// An empty ledger.
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }

    /// Append a contribution.
    pub fn attribute(&mut self, c: Contribution) {
        self.contributions.push(c);
    }

    /// Merge another ledger (pipeline concatenation; first-order additive
    /// composition — conservative for independent worst cases).
    pub fn compose(&mut self, other: &ErrorLedger) {
        self.contributions.extend(other.contributions.iter().cloned());
        self.declared_residual += other.declared_residual;
    }

    /// The completeness lint: every magnitude nonnegative and finite —
    /// no silent error mass, ever.
    ///
    /// # Errors
    /// The first [`LedgerDefect`] found.
    pub fn lint(&self) -> Result<(), LedgerDefect> {
        for c in &self.contributions {
            if !(c.abs >= 0.0) || !c.abs.is_finite() {
                return Err(LedgerDefect::BadContribution { label: c.label.clone() });
            }
        }
        if !(self.declared_residual >= 0.0) || !self.declared_residual.is_finite() {
            return Err(LedgerDefect::BadResidual);
        }
        Ok(())
    }

    /// Total attributed error bound (Σ contributions + declared residual).
    #[must_use]
    pub fn total(&self) -> f64 {
        self.contributions.iter().map(|c| c.abs).sum::<f64>() + self.declared_residual
    }

    /// Per-source subtotals, deterministic order.
    #[must_use]
    pub fn by_source(&self) -> Vec<(ErrorSource, f64)> {
        let mut sources = [
            ErrorSource::Geometry,
            ErrorSource::Discretization,
            ErrorSource::Algebraic,
            ErrorSource::Surrogate,
            ErrorSource::Statistical,
            ErrorSource::ModelForm,
        ]
        .map(|s| (s, 0.0f64));
        for c in &self.contributions {
            for slot in &mut sources {
                if slot.0 == c.source {
                    slot.1 += c.abs;
                }
            }
        }
        sources.into_iter().filter(|(_, v)| *v > 0.0).collect()
    }

    /// The dominant source and its mass (what escalation should attack).
    #[must_use]
    pub fn dominant(&self) -> Option<(ErrorSource, f64)> {
        self.by_source().into_iter().max_by(|a, b| a.1.total_cmp(&b.1))
    }

    /// The `explain_error` payload: one JSON object with per-source
    /// subtotals, entries, residual, and the dominant source.
    #[must_use]
    pub fn explain(&self) -> String {
        use std::fmt::Write as _;
        let mut out = String::from("{\"error_ledger\":{\"by_source\":{");
        for (i, (s, v)) in self.by_source().iter().enumerate() {
            if i > 0 {
                out.push(',');
            }
            let _ = write!(out, "\"{}\":{v:?}", s.name());
        }
        out.push_str("},\"entries\":[");
        for (i, c) in self.contributions.iter().enumerate() {
            if i > 0 {
                out.push(',');
            }
            let _ = write!(
                out,
                "{{\"label\":{:?},\"source\":\"{}\",\"abs\":{:?},\"rigor\":\"{:?}\"}}",
                c.label,
                c.source.name(),
                c.abs,
                c.rigor
            );
        }
        let _ = write!(
            out,
            "],\"residual\":{:?},\"total\":{:?},\"dominant\":{:?}}}}}",
            self.declared_residual,
            self.total(),
            self.dominant().map(|(s, _)| s.name())
        );
        out
    }
}

/// One pipeline stage's time accounting.
#[derive(Debug, Clone, PartialEq)]
pub struct TimeStage {
    /// Operator/stage name.
    pub op: String,
    /// Predicted (p10, p50, p90) seconds, when a model existed.
    pub predicted: Option<(f64, f64, f64)>,
    /// Measured wall seconds, when the stage ran.
    pub measured_s: Option<f64>,
}

/// A study's wall-clock attribution.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct TimeLedger {
    /// Stages in execution order.
    pub stages: Vec<TimeStage>,
}

impl TimeLedger {
    /// An empty ledger.
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }

    /// Append a stage.
    pub fn record(&mut self, stage: TimeStage) {
        self.stages.push(stage);
    }

    /// Total measured seconds (only stages that ran).
    #[must_use]
    pub fn total_measured_s(&self) -> f64 {
        self.stages.iter().filter_map(|s| s.measured_s).sum()
    }

    /// Total predicted median seconds (only stages with models).
    #[must_use]
    pub fn total_p50_s(&self) -> f64 {
        self.stages.iter().filter_map(|s| s.predicted.map(|p| p.1)).sum()
    }

    /// Fraction of measured stages whose actual landed inside [p10, p90]
    /// (the calibration audit; `None` when nothing is comparable).
    #[must_use]
    pub fn calibration(&self) -> Option<f64> {
        let comparable: Vec<&TimeStage> =
            self.stages.iter().filter(|s| s.predicted.is_some() && s.measured_s.is_some()).collect();
        if comparable.is_empty() {
            return None;
        }
        let inside = comparable
            .iter()
            .filter(|s| {
                let (p10, _, p90) = s.predicted.expect("filtered");
                let m = s.measured_s.expect("filtered");
                m >= p10 && m <= p90
            })
            .count();
        Some(inside as f64 / comparable.len() as f64)
    }

    /// The `explain_time` payload (JSON).
    #[must_use]
    pub fn explain(&self) -> String {
        use std::fmt::Write as _;
        let mut out = String::from("{\"time_ledger\":{\"stages\":[");
        for (i, s) in self.stages.iter().enumerate() {
            if i > 0 {
                out.push(',');
            }
            let _ = write!(out, "{{\"op\":{:?}", s.op);
            if let Some((p10, p50, p90)) = s.predicted {
                let _ = write!(out, ",\"p10\":{p10:?},\"p50\":{p50:?},\"p90\":{p90:?}");
            }
            if let Some(m) = s.measured_s {
                let _ = write!(out, ",\"measured\":{m:?}");
            }
            out.push('}');
        }
        let _ = write!(
            out,
            "],\"total_measured\":{:?},\"total_p50\":{:?},\"calibration\":{:?}}}}}",
            self.total_measured_s(),
            self.total_p50_s(),
            self.calibration()
        );
        out
    }
}
