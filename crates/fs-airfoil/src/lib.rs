//! fs-airfoil — generic airfoil section machinery (L2). Bead
//! frankensim-wf-root-guzez.5.1.1 (E4.0a, Wright Flyer program).
//!
//! Spec: COMPREHENSIVE_PLAN_FOR_REAL_TIME_WRIGHT_FLYER_SIM_WITH_FRANKENSIM.md
//! §5.2.1 (ROUND 6 steady state). Consumed by wing, canard, rudder AND
//! propeller crates — sections live at L2 precisely so no L3→L3 layering
//! violation is needed for props.
//!
//! E4.0a scope: section geometry under the frozen convention artifacts
//! (frame-conventions-v1 / geometry-conventions-v1), the two ANALYTIC
//! baselines (thin-airfoil with parabolic camber; fully-separated flat
//! plate), exact wind↔body decomposition and moment-reference transfer,
//! and the admission domain with typed refusals. Coefficient tables,
//! constrained B-spline residuals, coherent uncertainty, and indicial
//! kernels are E4.0b/E4.0c.
//!
//! Conventions (frozen, E1.4): angles in RADIANS; α positive nose-up;
//! pitch moment positive nose-up (about +y of frd-body-v1); section
//! moment reference lengths are the chord. cl/cd are wind-axis; cn/ca are
//! body-normal/axial with cn positive toward the suction side.

use fs_math::det;

/// Admitted angle-of-attack domain [rad], inclusive: the full circle
/// expressed as [−π, π]. Queries outside refuse (they indicate an
/// un-normalized angle upstream, not a physical state).
pub const MAX_ABS_ALPHA_RAD: f64 = core::f64::consts::PI;

/// Admitted chord domain [m] (Wright-era sections; refuse outside).
pub const MIN_CHORD_M: f64 = 0.01;
/// Upper chord bound [m] of the admitted domain.
pub const MAX_CHORD_M: f64 = 10.0;

/// Admitted camber-ratio domain (f/c). 0.15 is far beyond any Wright
/// section (1903 flew ~1/20 = 0.05) but inside thin-airfoil plausibility.
pub const MAX_CAMBER_RATIO: f64 = 0.15;

/// Admitted log10(Reynolds) domain for section queries. The Wright flight
/// envelope sits near 1e6; the domain brackets tunnel (1e4) through
/// full-scale (1e8). Outside → applicability-domain refusal.
pub const MIN_LOG10_RE: f64 = 4.0;
/// Upper log10(Reynolds) bound of the admitted domain.
pub const MAX_LOG10_RE: f64 = 8.0;

/// Drag coefficient of a 2-D flat plate normal to the stream (α = 90°),
/// the anchor of the separated baseline. Declared model constant
/// (classical bluff-body value); provenance-bound tables refine it.
pub const FLAT_PLATE_CD90: f64 = 1.98;

// The separated small-angle slope (CD90) must undercut the attached 2π slope
// or the baselines would be conflatable — enforced at compile time.
const _: () = assert!(FLAT_PLATE_CD90 < core::f64::consts::TAU);

/// A typed refusal (workspace law: refusals are data, never panics).
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Refusal {
    /// Stable machine-readable code.
    pub code: &'static str,
    /// Human-readable diagnosis.
    pub message: String,
    /// Ranked repairs, most likely fix first.
    pub ranked_repairs: Vec<String>,
}

/// Section geometry under the frozen conventions, provenance-bound.
#[derive(Clone, Debug, PartialEq)]
pub struct SectionGeometry {
    /// Chord [m].
    pub chord_m: f64,
    /// Max camber as a fraction of chord (f/c); parabolic camber line
    /// assumed by the analytic baseline.
    pub camber_ratio: f64,
    /// Provenance: the source-dossier record id this geometry derives from
    /// (dataset re-expression rule, frame-conventions-v1).
    pub dossier_record: String,
    /// Digitization class per uncertainty-v1 (e.g. "drawings-pm-2-5-mm").
    pub digitization_class: String,
}

/// Wind-axis section coefficients at one query point.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SectionCoefficients {
    /// Lift coefficient (wind axis, positive up).
    pub cl: f64,
    /// Drag coefficient (wind axis, positive downstream).
    pub cd: f64,
    /// Pitching moment about the quarter chord, positive nose-up.
    pub cm_quarter: f64,
}

/// Body-axis (normal/axial) representation with an explicit moment
/// reference station — the representation that stays meaningful through
/// deep stall where cl/cd lose their small-angle intuition.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct NormalAxial {
    /// Normal-force coefficient (perpendicular to chord, positive suction side).
    pub cn: f64,
    /// Axial (chordwise) force coefficient, positive toward the trailing edge.
    pub ca: f64,
    /// Pitching moment about `x_ref_over_c`, positive nose-up.
    pub cm_ref: f64,
    /// Moment reference station as x/c from the leading edge.
    pub x_ref_over_c: f64,
}

fn require_finite(name: &str, value: f64) -> Result<(), Refusal> {
    if value.is_finite() {
        Ok(())
    } else {
        Err(Refusal {
            code: "non-finite-input",
            message: format!("{name} = {value:?} must be finite"),
            ranked_repairs: vec!["check upstream unit conversions for NaN/inf".into()],
        })
    }
}

impl SectionGeometry {
    /// Validate the geometry against the admitted domains.
    ///
    /// # Errors
    /// `non-finite-input`, `chord-outside-domain`, `camber-outside-domain`,
    /// `provenance-missing`.
    pub fn admit(&self) -> Result<(), Refusal> {
        require_finite("chord_m", self.chord_m)?;
        require_finite("camber_ratio", self.camber_ratio)?;
        if !(MIN_CHORD_M..=MAX_CHORD_M).contains(&self.chord_m) {
            return Err(Refusal {
                code: "chord-outside-domain",
                message: format!(
                    "chord {} m outside admitted [{MIN_CHORD_M}, {MAX_CHORD_M}]",
                    self.chord_m
                ),
                ranked_repairs: vec!["check the metres conversion (1903 chord is 1.981 m)".into()],
            });
        }
        if !(0.0..=MAX_CAMBER_RATIO).contains(&self.camber_ratio) {
            return Err(Refusal {
                code: "camber-outside-domain",
                message: format!(
                    "camber ratio {} outside admitted [0, {MAX_CAMBER_RATIO}]",
                    self.camber_ratio
                ),
                ranked_repairs: vec![
                    "camber is a FRACTION of chord (1903: 0.05), not a percentage".into(),
                ],
            });
        }
        if self.dossier_record.is_empty() || self.digitization_class.is_empty() {
            return Err(Refusal {
                code: "provenance-missing",
                message: "section geometry must bind a dossier record and digitization class"
                    .into(),
                ranked_repairs: vec![
                    "cite the source-dossier-v1 record id this geometry derives from".into(),
                ],
            });
        }
        Ok(())
    }
}

/// Validate a query point (α, log10 Re) against the admitted domain.
///
/// # Errors
/// `non-finite-input`, `alpha-outside-domain`, `reynolds-outside-domain`
/// — the refusal STATES the admitted domain (applicability-domain law).
pub fn admit_query(alpha_rad: f64, log10_re: f64) -> Result<(), Refusal> {
    require_finite("alpha_rad", alpha_rad)?;
    require_finite("log10_re", log10_re)?;
    if alpha_rad.abs() > MAX_ABS_ALPHA_RAD {
        return Err(Refusal {
            code: "alpha-outside-domain",
            message: format!(
                "alpha {alpha_rad} rad outside admitted [-{MAX_ABS_ALPHA_RAD}, {MAX_ABS_ALPHA_RAD}]"
            ),
            ranked_repairs: vec!["normalize the angle to (-pi, pi] upstream".into()],
        });
    }
    if !(MIN_LOG10_RE..=MAX_LOG10_RE).contains(&log10_re) {
        return Err(Refusal {
            code: "reynolds-outside-domain",
            message: format!(
                "log10(Re) {log10_re} outside admitted [{MIN_LOG10_RE}, {MAX_LOG10_RE}]"
            ),
            ranked_repairs: vec![
                "the section database does not extrapolate; widen the domain only with sourced data"
                    .into(),
            ],
        });
    }
    Ok(())
}

/// Thin-airfoil baseline for a PARABOLIC camber line with max camber
/// `camber_ratio` at mid-chord (classical exact results):
///
/// - zero-lift angle  α₀ = −2·(f/c)
/// - lift             cl = 2π·(α + 2·f/c)
/// - quarter-chord moment cm_c/4 = −π·(f/c)  (aerodynamic centre at c/4)
/// - cd = 0 (inviscid; the viscous residual is E4.0b table territory)
///
/// Valid as a BASELINE in the attached regime only; callers needing
/// post-stall use [`flat_plate_separated`] or (later) blended tables.
///
/// # Errors
/// Query refusals from [`admit_query`]; camber refusals as in
/// [`SectionGeometry::admit`].
pub fn thin_airfoil(
    alpha_rad: f64,
    camber_ratio: f64,
    log10_re: f64,
) -> Result<SectionCoefficients, Refusal> {
    admit_query(alpha_rad, log10_re)?;
    require_finite("camber_ratio", camber_ratio)?;
    if !(0.0..=MAX_CAMBER_RATIO).contains(&camber_ratio) {
        return Err(Refusal {
            code: "camber-outside-domain",
            message: format!(
                "camber ratio {camber_ratio} outside admitted [0, {MAX_CAMBER_RATIO}]"
            ),
            ranked_repairs: vec!["camber is a fraction of chord".into()],
        });
    }
    let two_pi = 2.0 * core::f64::consts::PI;
    Ok(SectionCoefficients {
        cl: two_pi * (alpha_rad + 2.0 * camber_ratio),
        cd: 0.0,
        cm_quarter: -core::f64::consts::PI * camber_ratio,
    })
}

/// Fully-separated flat-plate baseline, valid over the WHOLE α circle:
///
/// - normal force  cn = CD90·sin α   (odd in α; |cn| = CD90 at ±90°)
/// - axial force   ca = 0            (no leading-edge suction when separated)
/// - centre of pressure x_cp/c = 0.25 + 0.25·(1 − cos α)  (c/4 at small α
///   → mid-chord at 90°, symmetric continuation beyond)
/// - cm about c/4:  cm = −cn·(x_cp − 0.25)
///
/// This is a DECLARED analytic model (documented shape, classical CD90
/// anchor), the low envelope of the post-stall regime; measured tables
/// refine it in E4.0b. It makes no Wright-specific claim (the synthesized
/// post-stall source record a2-synthesized-stall owns that boundary).
///
/// # Errors
/// Query refusals from [`admit_query`].
pub fn flat_plate_separated(alpha_rad: f64, log10_re: f64) -> Result<NormalAxial, Refusal> {
    admit_query(alpha_rad, log10_re)?;
    let s = det::sin(alpha_rad);
    let c = det::cos(alpha_rad);
    let cn = FLAT_PLATE_CD90 * s;
    let x_cp = 0.25 + 0.25 * (1.0 - c);
    Ok(NormalAxial {
        cn,
        ca: 0.0,
        cm_ref: -cn * (x_cp - 0.25),
        x_ref_over_c: 0.25,
    })
}

/// Wind-axis (cl, cd) → body-axis (cn, ca) at angle of attack α.
/// Exact rotation: cn = cl·cos α + cd·sin α; ca = cd·cos α − cl·sin α.
#[must_use]
pub fn wind_to_body(cl: f64, cd: f64, alpha_rad: f64) -> (f64, f64) {
    let (s, c) = (det::sin(alpha_rad), det::cos(alpha_rad));
    (cl * c + cd * s, cd * c - cl * s)
}

/// Body-axis (cn, ca) → wind-axis (cl, cd): the inverse rotation.
#[must_use]
pub fn body_to_wind(cn: f64, ca: f64, alpha_rad: f64) -> (f64, f64) {
    let (s, c) = (det::sin(alpha_rad), det::cos(alpha_rad));
    (cn * c - ca * s, cn * s + ca * c)
}

/// Transfer a pitching moment between reference stations (x/c from the
/// leading edge): cm_B = cm_A + cn·(x_B − x_A). Exact for any cn;
/// axial force contributes no pitch at zero camber-line offset (section
/// convention: moments taken about the chord line).
#[must_use]
pub fn transfer_moment(cm_a: f64, cn: f64, x_a_over_c: f64, x_b_over_c: f64) -> f64 {
    cm_a + cn * (x_b_over_c - x_a_over_c)
}
