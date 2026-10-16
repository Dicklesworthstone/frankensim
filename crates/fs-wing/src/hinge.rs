//! Whole-canard generalized hinge-load interface (bead
//! wf-root-guzez.5.4, E4.2b). Plan §5.1.3: fs-wing integrates the
//! COUPLED-solution hinge load about the modeled hinge axis — fs-airfoil
//! supplies section data only and never sees the aircraft-level hinge
//! load. Acceleration-dependent apparent-mass terms are EXCLUDED BY
//! CONSTRUCTION: this interface consumes only the circulatory solution
//! (Γ, freestream) plus quasi-steady section couples; apparent mass
//! lives in the added-mass blocks (fs-flyer, AddedMassOnly owner), and
//! the battery's hostile twin shows the quadrature oracle catching any
//! smuggled non-circulatory term at steady state.

use crate::{Panel, Refusal, SurfaceId};

/// Axis-unit admission tolerance: | |a| − 1 | ≤ this.
pub const AXIS_UNIT_TOL: f64 = 1.0e-9;

/// A hinge axis: point + unit direction (FRD body frame).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct HingeAxis {
    /// A point on the axis [m].
    pub point_m: [f64; 3],
    /// Unit direction (canard pivot: spanwise, +y).
    pub axis_unit: [f64; 3],
}

/// One strip's quasi-steady section couple (from the caller's section
/// closure — a pure couple, transported freely to the axis).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SectionCouple {
    /// Couple magnitude about the strip's spanwise direction [N·m]
    /// (positive nose-up in the strip frame).
    pub moment_nm: f64,
    /// Strip spanwise unit direction.
    pub span_unit: [f64; 3],
}

/// Per-panel hinge contribution (the per-item oracle surface).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct HingeItem {
    /// Panel index in the full system.
    pub panel: usize,
    /// This panel's moment about the axis [N·m].
    pub moment_nm: f64,
}

/// The integrated hinge load.
#[derive(Clone, Debug, PartialEq)]
pub struct HingeLoadReport {
    /// Total circulatory hinge moment about the axis [N·m].
    pub circulatory_nm: f64,
    /// Section-couple contribution [N·m].
    pub section_nm: f64,
    /// Grand total [N·m].
    pub total_nm: f64,
    /// Per-panel circulatory items (selected surfaces only).
    pub items: Vec<HingeItem>,
    /// Net force on the selected surfaces [N] (for axis-shift identities).
    pub force_n: [f64; 3],
}

impl HingeAxis {
    /// Admit the axis.
    ///
    /// # Errors
    /// `hinge-axis-invalid` (non-finite point, non-unit direction —
    /// tested at the tolerance cap AND one ulp past it).
    pub fn admit(&self) -> Result<(), Refusal> {
        let finite = self
            .point_m
            .iter()
            .chain(self.axis_unit.iter())
            .all(|v| v.is_finite());
        let n = (self.axis_unit[0] * self.axis_unit[0]
            + self.axis_unit[1] * self.axis_unit[1]
            + self.axis_unit[2] * self.axis_unit[2])
            .sqrt();
        if !finite || (n - 1.0).abs() > AXIS_UNIT_TOL {
            return Err(Refusal {
                code: "hinge-axis-invalid",
                message: format!("point {:?}, |axis| = {n:?}", self.point_m),
                ranked_repairs: vec!["normalize the axis direction".into()],
            });
        }
        Ok(())
    }
}

/// Integrate the whole-surface hinge load about the axis from the
/// coupled solution: per selected panel, F = ρΓ(V×seg) at the bound
/// midpoint, m = ((r − p) × F)·â; section couples project as free
/// vectors (m·(ŝ·â)). NO acceleration inputs exist on this interface —
/// apparent mass cannot enter by construction.
///
/// # Errors
/// `hinge-axis-invalid`; `hinge-selection-empty` (no surfaces named or
/// none matched); `gamma-length-mismatch`; `freestream-invalid`.
#[allow(clippy::too_many_arguments)]
pub fn hinge_load(
    panels: &[Panel],
    gamma: &[f64],
    freestream_mps: [f64; 3],
    rho_kg_m3: f64,
    surfaces: &[SurfaceId],
    axis: &HingeAxis,
    section_couples: &[SectionCouple],
) -> Result<HingeLoadReport, Refusal> {
    axis.admit()?;
    if surfaces.is_empty() {
        return Err(Refusal {
            code: "hinge-selection-empty",
            message: "no surfaces named".into(),
            ranked_repairs: vec!["name the hinged surfaces (e.g. both canard planes)".into()],
        });
    }
    if gamma.len() != panels.len() {
        return Err(Refusal {
            code: "gamma-length-mismatch",
            message: format!("{} panels vs {} gamma", panels.len(), gamma.len()),
            ranked_repairs: vec!["pass the coupled solution's gamma unmodified".into()],
        });
    }
    let finite = freestream_mps.iter().all(|v| v.is_finite());
    if !finite || rho_kg_m3 <= 0.0 || !rho_kg_m3.is_finite() {
        return Err(Refusal {
            code: "freestream-invalid",
            message: format!("V {freestream_mps:?}, rho {rho_kg_m3:?}"),
            ranked_repairs: vec!["evaluate at a physical air state".into()],
        });
    }
    let a = axis.axis_unit;
    let p = axis.point_m;
    let mut items = Vec::new();
    let mut total = 0.0f64;
    let mut force = [0.0f64; 3];
    let mut matched = false;
    for (j, panel) in panels.iter().enumerate() {
        if !surfaces.contains(&panel.surface) {
            continue;
        }
        matched = true;
        let seg = [
            panel.b[0] - panel.a[0],
            panel.b[1] - panel.a[1],
            panel.b[2] - panel.a[2],
        ];
        let s = rho_kg_m3 * gamma[j];
        let f = [
            s * (freestream_mps[1] * seg[2] - freestream_mps[2] * seg[1]),
            s * (freestream_mps[2] * seg[0] - freestream_mps[0] * seg[2]),
            s * (freestream_mps[0] * seg[1] - freestream_mps[1] * seg[0]),
        ];
        let r = [
            0.5 * (panel.a[0] + panel.b[0]) - p[0],
            0.5 * (panel.a[1] + panel.b[1]) - p[1],
            0.5 * (panel.a[2] + panel.b[2]) - p[2],
        ];
        let rxf = [
            r[1] * f[2] - r[2] * f[1],
            r[2] * f[0] - r[0] * f[2],
            r[0] * f[1] - r[1] * f[0],
        ];
        let m = rxf[0] * a[0] + rxf[1] * a[1] + rxf[2] * a[2];
        items.push(HingeItem {
            panel: j,
            moment_nm: m,
        });
        total += m;
        force[0] += f[0];
        force[1] += f[1];
        force[2] += f[2];
    }
    if !matched {
        return Err(Refusal {
            code: "hinge-selection-empty",
            message: format!("{surfaces:?} matched no panels"),
            ranked_repairs: vec!["check the surface ids against the layout".into()],
        });
    }
    let mut section = 0.0f64;
    for c in section_couples {
        let dot = c.span_unit[0] * a[0] + c.span_unit[1] * a[1] + c.span_unit[2] * a[2];
        section += c.moment_nm * dot;
    }
    Ok(HingeLoadReport {
        circulatory_nm: total,
        section_nm: section,
        total_nm: total + section,
        items,
        force_n: force,
    })
}
