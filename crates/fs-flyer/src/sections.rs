//! Wright Flyer section datasets on fs-airfoil (bead wf-root-guzez.5.2,
//! E4.1). The four SurfaceKind datasets built from what the FROZEN
//! registry admits today: analytic baselines (fs-airfoil E4.0a) + a
//! declared zero attached-regime residual (Estimated) + the flat-plate
//! post-stall continuation under the a2-synthesized-stall record's
//! prior/trend role. NOTHING here is fitted to any measurement: the 1901
//! anchors are holdout by doctrine (the battery checks trends against
//! them; the builder never sees them).
//!
//! Full-scale table ingestion (Ames/LFST) remains DATA-BLOCKED — the
//! reductions were never openly published (E1.2/E1.6 research); when
//! tables arrive they ingest under reexpression-v1 with partition
//! assignment, replacing the zero residuals dataset-by-dataset.

use crate::Refusal;
use fs_airfoil::fit::{BsplineAxis, ResidualSurface};
use fs_airfoil::table::{CoefficientTable, ConventionBlock, RegimePatch, SurfaceKind};
use fs_airfoil::{flat_plate_separated, thin_airfoil};
use fs_math::det;
use fs_wing::nonlinear::StripRegime;

/// Shared quasi-steady section regime for the aircraft's existing Estimated
/// thin-airfoil/flat-plate model. Both force channels must use the same weight
/// at the same solved effective incidence, including twist and induction.
fn separated_weight(alpha_rad: f64) -> (f64, StripRegime) {
    let magnitude = alpha_rad.abs();
    if magnitude <= 0.30 {
        (0.0, StripRegime::Attached)
    } else if magnitude < 0.45 {
        let t = (magnitude - 0.30) / 0.15;
        (t * t * (3.0 - 2.0 * t), StripRegime::Blended)
    } else {
        (1.0, StripRegime::Separated)
    }
}

pub(crate) fn section_lift(alpha_rad: f64, camber: f64) -> (f64, StripRegime) {
    let (weight, regime) = separated_weight(alpha_rad);
    let attached = core::f64::consts::TAU * (alpha_rad + 2.0 * camber);
    if regime == StripRegime::Attached {
        return (attached, regime);
    }
    let separated = fs_airfoil::FLAT_PLATE_CD90 * det::sin(alpha_rad) * det::cos(alpha_rad);
    if regime == StripRegime::Separated {
        (separated, regime)
    } else {
        (attached * (1.0 - weight) + separated * weight, regime)
    }
}

/// Non-lifting profile plus separated pressure drag, without induced drag.
/// The existing flat-plate normal force Cn=Cd90*sin(alpha), Ca=0 transforms
/// to Cd=Cd90*sin²(alpha). This must accompany its separated lift channel:
/// leaving the attached Cd in place discards the plate's pressure work.
pub(crate) fn section_profile_drag(alpha_rad: f64, attached_cd: f64) -> f64 {
    let (weight, regime) = separated_weight(alpha_rad);
    if regime == StripRegime::Attached {
        return attached_cd;
    }
    let sin_alpha = det::sin(alpha_rad);
    let separated = fs_airfoil::FLAT_PLATE_CD90 * sin_alpha * sin_alpha;
    if regime == StripRegime::Separated {
        separated
    } else {
        attached_cd * (1.0 - weight) + separated * weight
    }
}

/// Evidence-lineage metadata every dataset carries (E4.1 DONE-WHEN).
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct SectionLineage {
    /// Source-dossier record id.
    pub dossier_record: &'static str,
    /// Independence group.
    pub independence_group: &'static str,
    /// Evidence-color ceiling for values from this dataset.
    pub ceiling: &'static str,
    /// Role restriction (the synthesized-stall law).
    pub role: &'static str,
}

/// One surface's section dataset: baseline camber + residual table +
/// lineage.
#[derive(Clone, Debug, PartialEq)]
pub struct SectionDataset {
    /// Which surface.
    pub kind: SurfaceKind,
    /// Baseline camber ratio (f/c) fed to the thin-airfoil baseline.
    pub camber_ratio: f64,
    /// The residual table (attached regime; zero-residual declared v1).
    pub residual: CoefficientTable,
    /// Lineage metadata.
    pub lineage: SectionLineage,
}

fn conventions() -> ConventionBlock {
    ConventionBlock {
        axes_id: "frd-body-v1".into(),
        moment_signs_id: "moment-signs-v1".into(),
        angles_id: "angles-v1".into(),
        reference_area_convention: "per-surface projected (geometry-conventions-v1)".into(),
        wind_reference: "not-applicable".into(),
    }
}

fn zero_residual(kind: SurfaceKind, record: &'static str) -> CoefficientTable {
    // A declared ZERO residual over the attached regime: the baseline IS
    // the v1 model (Estimated). Zero coefficients trivially satisfy every
    // shape constraint and C0 continuity.
    let surface = ResidualSurface {
        axes: [
            BsplineAxis {
                name: "alpha_rad",
                lo: -0.35,
                hi: 0.35,
                n_coef: 4,
            },
            BsplineAxis {
                name: "log10_re",
                lo: 4.0,
                hi: 8.0,
                n_coef: 1,
            },
            BsplineAxis {
                name: "delta_rad",
                lo: -0.6,
                hi: 0.6,
                n_coef: 1,
            },
        ],
        coef: vec![0.0; 4],
        constraints: vec![],
    };
    CoefficientTable {
        kind,
        channel: "cl-residual",
        dossier_record: record.into(),
        conventions: conventions(),
        patches: vec![RegimePatch {
            regime: "attached",
            surface,
        }],
    }
}

/// Build the four v1 datasets (wing, canard, rudder, prop).
///
/// # Errors
/// Propagates table validation refusals (structurally none for the zero
/// residual, but the gate RUNS — a corrupted convention id would refuse).
pub fn build_v1_datasets() -> Result<Vec<SectionDataset>, Refusal> {
    let specs: [(SurfaceKind, f64, &'static str, &'static str, &'static str); 4] = [
        (
            SurfaceKind::Wing,
            0.05,
            "a1-wright-1901-tunnel",
            "wright-1901-tunnel",
            "flown 1/20 camber (flyer-reference camber_ratio, verified)",
        ),
        (
            SurfaceKind::Canard,
            0.05,
            "a2-simmodels-deters",
            "ames-aiaa-1999",
            "canard camber assumed = wing class (variable-camber mechanism is E4.6b's model axis)",
        ),
        (
            SurfaceKind::Rudder,
            0.0,
            "a2-simmodels-deters",
            "ames-aiaa-1999",
            "symmetric vertical surfaces",
        ),
        (
            SurfaceKind::Prop,
            0.05,
            "a2-props-bentend",
            "lfst-wright-experience",
            "prop sections are an Estimated reconstruction (prop-geometry-v1 1903-absence rule)",
        ),
    ];
    let mut out = Vec::with_capacity(4);
    for (kind, camber, record, group, _why) in specs {
        let residual = zero_residual(kind, record);
        residual.validate(1e-9).map_err(|e| Refusal {
            code: e.code,
            message: format!("{kind:?}: {}", e.message),
            ranked_repairs: e.ranked_repairs,
        })?;
        out.push(SectionDataset {
            kind,
            camber_ratio: camber,
            residual,
            lineage: SectionLineage {
                dossier_record: record,
                independence_group: group,
                ceiling: "Estimated",
                role: if kind == SurfaceKind::Prop {
                    "reconstruction; BEMT method-validation uses the 1911 calibration curves"
                } else {
                    "baseline+zero-residual v1; post-stall continuation under a2-synthesized-stall prior/trend role ONLY (never Wright-specific deep-stall validation)"
                },
            },
        });
    }
    Ok(out)
}

/// Evaluate a dataset's lift coefficient at (alpha, log10 Re): 2-D
/// thin-airfoil baseline + residual in the attached regime, blending to
/// the flat-plate separated branch beyond it; finite-wing corrected for
/// aspect ratio `ar` (lifting-line: CL3D = CL2D / (1 + 2/AR) — a
/// DOCUMENTED simple model; the real multisurface solve is E4.2).
///
/// # Errors
/// fs-airfoil admission refusals.
pub fn cl_3d(ds: &SectionDataset, alpha_rad: f64, log10_re: f64, ar: f64) -> Result<f64, Refusal> {
    let attached_limit = 0.30; // blend start [rad]
    let full_sep = 0.45;
    let map_err = |e: fs_airfoil::Refusal| Refusal {
        code: e.code,
        message: e.message,
        ranked_repairs: e.ranked_repairs,
    };
    let thin =
        thin_airfoil(alpha_rad.clamp(-0.35, 0.35), ds.camber_ratio, log10_re).map_err(map_err)?;
    let resid = ds
        .residual
        .eval([alpha_rad.clamp(-0.35, 0.35), log10_re, 0.0])
        .map_err(map_err)?;
    let cl_attached = (thin.cl + resid) / (1.0 + 2.0 / ar);
    let sep = flat_plate_separated(alpha_rad, log10_re).map_err(map_err)?;
    let (cl_sep, _) = fs_airfoil::body_to_wind(sep.cn, sep.ca, alpha_rad);
    let a = alpha_rad.abs();
    if a <= attached_limit {
        Ok(cl_attached)
    } else if a >= full_sep {
        Ok(cl_sep)
    } else {
        let t = (a - attached_limit) / (full_sep - attached_limit);
        let s = t * t * (3.0 - 2.0 * t);
        Ok(cl_attached * (1.0 - s) + cl_sep * s)
    }
}

#[cfg(test)]
mod drag_tests {
    use super::{section_lift, section_profile_drag};
    use fs_wing::nonlinear::StripRegime;

    #[test]
    fn g0_separated_lift_and_drag_recover_the_same_plate_normal_force() {
        for alpha in [-1.2_f64, -0.7, -0.45, 0.45, 0.7, 1.2] {
            let (cl, regime) = section_lift(alpha, 0.05);
            let cd = section_profile_drag(alpha, 0.015);
            let (cn, ca) = fs_airfoil::wind_to_body(cl, cd, alpha);
            assert_eq!(regime, StripRegime::Separated);
            assert!((cn - fs_airfoil::FLAT_PLATE_CD90 * alpha.sin()).abs() < 1e-12);
            assert!(
                ca.abs() < 1e-12,
                "a separated plate carries no suction force"
            );
            assert!(
                cd > 0.0,
                "profile work must dissipate in both incidence directions"
            );
        }
    }

    #[test]
    fn g0_drag_blend_is_even_continuous_and_keeps_the_attached_baseline() {
        let cd0 = 0.015;
        for alpha in [0.0, 0.15, 0.30] {
            assert_eq!(section_profile_drag(alpha, cd0).to_bits(), cd0.to_bits());
        }
        let mut last = cd0;
        for k in 0..=300 {
            let alpha = 0.30 + 0.15 * f64::from(k) / 300.0;
            let cd = section_profile_drag(alpha, cd0);
            assert_eq!(cd, section_profile_drag(-alpha, cd0));
            assert!(cd >= last);
            last = cd;
        }
        for boundary in [0.30, 0.45] {
            assert!(
                (section_profile_drag(boundary - 1e-8, cd0)
                    - section_profile_drag(boundary + 1e-8, cd0))
                .abs()
                    < 5e-8
            );
        }
    }
}
