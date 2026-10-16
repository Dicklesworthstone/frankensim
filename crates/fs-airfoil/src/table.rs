//! Provenance-bound coefficient tables (bead wf-root-guzez.5.1.2, E4.0b).
//! Every table binds a source-dossier record AND the frozen convention
//! block (frame-conventions-v1 reexpression-v1 rule) — a table without
//! them is a typed refusal, structurally. Wing / canard / rudder / prop
//! tables are SEPARATE (plan §5.2.1): a `SurfaceKind` tag is part of the
//! table identity and cross-kind reuse refuses.

use crate::Refusal;
use crate::fit::{ResidualSurface, verify_regime_continuity};

/// Frozen convention ids this crate binds against (E1.4 artifacts).
pub const AXES_ID: &str = "frd-body-v1";
/// Frozen moment-sign convention id.
pub const MOMENT_SIGNS_ID: &str = "moment-signs-v1";
/// Frozen angle convention id.
pub const ANGLES_ID: &str = "angles-v1";

/// Which surface a table describes. Tables are surface-specific; the 1901
/// tunnel trends may inform several, but each table carries its own record.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SurfaceKind {
    /// Main-wing sections.
    Wing,
    /// Canard (front elevator) sections.
    Canard,
    /// Vertical rudder sections.
    Rudder,
    /// Propeller blade sections.
    Prop,
}

/// The dataset re-expression block (reexpression-v1). Every field must
/// resolve to the frozen ids; `wind_reference` names the instrument
/// lineage of any wind-referred quantity ("not-applicable" for tunnel
/// tables with their own q∞).
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ConventionBlock {
    /// Body-axes id (must equal [`AXES_ID`]).
    pub axes_id: String,
    /// Moment-sign id (must equal [`MOMENT_SIGNS_ID`]).
    pub moment_signs_id: String,
    /// Angle-convention id (must equal [`ANGLES_ID`]).
    pub angles_id: String,
    /// Reference-area convention (geometry-conventions-v1 language).
    pub reference_area_convention: String,
    /// WindReference lineage or "not-applicable".
    pub wind_reference: String,
}

impl ConventionBlock {
    /// Validate against the frozen ids (fail-closed, field named).
    ///
    /// # Errors
    /// `convention-block-mismatch` naming the offending field;
    /// `convention-block-missing` for empty free-text fields.
    pub fn validate(&self) -> Result<(), Refusal> {
        let mismatch = |field: &str, got: &str, want: &str| {
            Refusal {
            code: "convention-block-mismatch",
            message: format!("{field} = {got:?}, frozen id is {want:?}"),
            ranked_repairs: vec![
                "re-express the dataset at ingest under the frozen conventions — never patch at query time".into(),
            ],
        }
        };
        if self.axes_id != AXES_ID {
            return Err(mismatch("axes_id", &self.axes_id, AXES_ID));
        }
        if self.moment_signs_id != MOMENT_SIGNS_ID {
            return Err(mismatch(
                "moment_signs_id",
                &self.moment_signs_id,
                MOMENT_SIGNS_ID,
            ));
        }
        if self.angles_id != ANGLES_ID {
            return Err(mismatch("angles_id", &self.angles_id, ANGLES_ID));
        }
        if self.reference_area_convention.is_empty() || self.wind_reference.is_empty() {
            return Err(Refusal {
                code: "convention-block-missing",
                message: "reference_area_convention and wind_reference must be declared".into(),
                ranked_repairs: vec![
                    "state the geometry-conventions-v1 area convention and the WindReference lineage (or 'not-applicable')".into(),
                ],
            });
        }
        Ok(())
    }
}

/// One regime patch: a fitted residual surface tagged with its regime.
#[derive(Clone, Debug, PartialEq)]
pub struct RegimePatch {
    /// Regime label (attached / transitional / separated / post-stall).
    pub regime: &'static str,
    /// The fitted residual surface over (α, log Re, δ).
    pub surface: ResidualSurface,
}

/// A provenance-bound coefficient table: ordered regime patches tiling an
/// α interval, for ONE surface kind and ONE coefficient channel.
#[derive(Clone, Debug, PartialEq)]
pub struct CoefficientTable {
    /// Which surface this table belongs to.
    pub kind: SurfaceKind,
    /// Coefficient channel ("cl-residual", "cd-residual", "cm-residual").
    pub channel: &'static str,
    /// Source-dossier record id (source-dossier-v1).
    pub dossier_record: String,
    /// The re-expression block (validated fail-closed).
    pub conventions: ConventionBlock,
    /// Regime patches in ascending α order, abutting exactly.
    pub patches: Vec<RegimePatch>,
}

impl CoefficientTable {
    /// Validate the whole table: provenance, conventions, per-patch
    /// constraints, and C⁰ continuity across every regime boundary.
    ///
    /// # Errors
    /// `provenance-missing`, convention refusals, `table-empty`, patch
    /// constraint refusals, and regime-continuity refusals (all typed).
    pub fn validate(&self, continuity_tol: f64) -> Result<(), Refusal> {
        if self.dossier_record.is_empty() {
            return Err(Refusal {
                code: "provenance-missing",
                message: format!(
                    "{:?}/{} table has no dossier record",
                    self.kind, self.channel
                ),
                ranked_repairs: vec!["bind the source-dossier-v1 record id".into()],
            });
        }
        self.conventions.validate()?;
        if self.patches.is_empty() {
            return Err(Refusal {
                code: "table-empty",
                message: "a coefficient table must carry at least one regime patch".into(),
                ranked_repairs: vec!["fit at least the attached regime".into()],
            });
        }
        for patch in &self.patches {
            patch.surface.verify_constraints()?;
        }
        for pair in self.patches.windows(2) {
            verify_regime_continuity(&pair[0].surface, &pair[1].surface, continuity_tol)?;
        }
        Ok(())
    }

    /// Evaluate the residual at (α, log Re, δ), selecting the regime patch
    /// by α. The table must have been validated first.
    ///
    /// # Errors
    /// `alpha-outside-table` when α falls outside every patch (the caller's
    /// `admit_query` bounds the global domain; this bounds the FITTED one —
    /// the applicability-domain law distinguishes them).
    pub fn eval(&self, x: [f64; 3]) -> Result<f64, Refusal> {
        for patch in &self.patches {
            let ax = &patch.surface.axes[0];
            if x[0] >= ax.lo && x[0] <= ax.hi {
                return Ok(patch.surface.eval(x));
            }
        }
        Err(Refusal {
            code: "alpha-outside-table",
            message: format!(
                "alpha {} outside the fitted regime partition [{}, {}]",
                x[0],
                self.patches
                    .first()
                    .map_or(f64::NAN, |p| p.surface.axes[0].lo),
                self.patches
                    .last()
                    .map_or(f64::NAN, |p| p.surface.axes[0].hi)
            ),
            ranked_repairs: vec![
                "the table does not extrapolate; widen the fit only with sourced data".into(),
            ],
        })
    }

    /// Strict fitted-domain evaluation (E4.0c): refuses when ANY axis of
    /// the selected patch is exceeded, instead of the spline's silent
    /// clamp on the secondary axes. This is the applicability-domain law
    /// for queries: out-of-domain gets a refusal that STATES the fitted
    /// box, never an unconstrained extrapolation.
    ///
    /// # Errors
    /// `alpha-outside-table` (no patch covers α);
    /// `query-outside-fitted-domain` (the covering patch's log Re or δ
    /// box is exceeded — the box is stated).
    pub fn eval_strict(&self, x: [f64; 3]) -> Result<f64, Refusal> {
        for patch in &self.patches {
            let ax = &patch.surface.axes;
            if x[0] >= ax[0].lo && x[0] <= ax[0].hi {
                for (axis, value) in ax.iter().zip(x.iter()).skip(1) {
                    if axis.n_coef != 1 && (*value < axis.lo || *value > axis.hi) {
                        return Err(Refusal {
                            code: "query-outside-fitted-domain",
                            message: format!(
                                "{} = {} outside the fitted box [{}, {}] of the covering regime patch",
                                axis.name, value, axis.lo, axis.hi
                            ),
                            ranked_repairs: vec![
                                "the table does not extrapolate; the analytic baseline alone applies out here".into(),
                                "widen the fitted box only with sourced data".into(),
                            ],
                        });
                    }
                }
                return Ok(patch.surface.eval(x));
            }
        }
        self.eval(x) // reuse the alpha-outside-table refusal path
    }
}
