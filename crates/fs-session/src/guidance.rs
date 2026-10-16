//! ERRORS AS GUIDANCE (plan §11.3): every failure a structured value —
//! machine-readable code, diagnosis, and candidate fixes RANKED BY THE
//! COST MODEL. Admission findings (the canonical `BudgetInfeasible`
//! worked example) surface through this same channel, as do constraint
//! infeasibility diagnoses upstream. "A refusal that teaches is worth
//! ten silent successes."

use fs_ir::admission::{Finding, RankedFix};
use std::fmt::Write as _;

/// A teaching failure.
#[derive(Debug, Clone, PartialEq)]
pub struct Guidance {
    /// Stable machine-readable code (e.g. `budget-infeasible`).
    pub code: String,
    /// What went wrong, with context.
    pub diagnosis: String,
    /// Candidate fixes, best first (cost-model-ranked where estimable).
    pub fixes: Vec<RankedFix>,
}

impl Guidance {
    /// Lift an admission finding into guidance (the enforcement channel).
    #[must_use]
    pub fn from_finding(f: &Finding) -> Guidance {
        Guidance {
            code: format!("{}-rejection", f.check),
            diagnosis: f.what.clone(),
            fixes: f.fixes.clone(),
        }
    }

    /// Canonical rendering (log payloads; deterministic).
    #[must_use]
    pub fn render(&self) -> String {
        let mut out = format!("[{}] {}\n", self.code, self.diagnosis);
        for (i, fix) in self.fixes.iter().enumerate() {
            let wall = fix
                .predicted_wall_s
                .map_or(String::new(), |w| format!(" (predicted wall {w:.1}s)"));
            let _ = writeln!(out, "  fix#{i}: {}{wall} — {}", fix.action, fix.qoi_impact);
        }
        out
    }
}
