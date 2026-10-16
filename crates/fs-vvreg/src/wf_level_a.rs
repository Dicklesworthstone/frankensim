//! Wright Flyer Level-A/B corpus rows (bead
//! `frankensim-wf-root-guzez.11.6.1`, E10.2-i).
//!
//! These rows AGGREGATE the Wright Flyer program's executed V/H
//! evidence into the corpus: each row's reference value is the
//! measure-then-pinned number its battery executed, and the retained
//! manifest carries the receipt digest that binds it (V-06a images,
//! V-08b1 dense referee, V-10 hybrid-wake shape, V-20 blade cover,
//! the E10.1 harness aggregate, and the H-07 synthetic recovery).
//! This is aggregation, not first execution — every cited receipt
//! was executed by its own battery; the corpus records identity and
//! applicability, and every query remains non-certifying.

use fs_qty::Dims;

/// Retained tab-separated manifest backing every Wright Flyer row.
pub(crate) const WF_LEVEL_A_MANIFEST: &[u8] =
    include_bytes!("../../../data/vv-corpus/wright-flyer/wf-level-a-v1.tsv");

/// Dimensionless.
pub const DIMENSIONLESS: Dims = Dims([0, 0, 0, 0, 0, 0]);
/// Coherent-SI dimensions of length.
pub const WF_LENGTH_DIMS: Dims = Dims([1, 0, 0, 0, 0, 0]);
/// Coherent-SI dimensions of speed.
pub const WF_SPEED_DIMS: Dims = Dims([1, 0, -1, 0, 0, 0]);

/// Wright Flyer coverage family.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum WfFamily {
    /// Flat-plane image system (V-06a class).
    GroundImages,
    /// Dense unsteady prescribed-wake referee (V-08b1 class).
    UnsteadyReferee,
    /// fs-vpm hybrid wake vs the dense reference (V-10 class).
    HybridWake,
    /// Swept-contact blade capsules (V-20 class).
    SweptContact,
    /// The E10.1 referee-harness aggregate.
    RefereeHarness,
    /// The H-07 historical-inference machinery.
    HistoricalInference,
}

impl WfFamily {
    /// Stable manifest spelling.
    #[must_use]
    pub const fn name(self) -> &'static str {
        match self {
            Self::GroundImages => "ground-images",
            Self::UnsteadyReferee => "unsteady-referee",
            Self::HybridWake => "hybrid-wake",
            Self::SweptContact => "swept-contact",
            Self::RefereeHarness => "referee-harness",
            Self::HistoricalInference => "historical-inference",
        }
    }
}

/// Whether a row's referee is analytic or an independent in-repo code.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum WfKind {
    /// Analytic/closed-form or structurally exact reference (level A).
    AnalyticReference,
    /// Independent in-repo formulation comparison (level B).
    CrossCodeReference,
}

impl WfKind {
    /// Stable manifest spelling.
    #[must_use]
    pub const fn name(self) -> &'static str {
        match self {
            Self::AnalyticReference => "analytic-reference",
            Self::CrossCodeReference => "cross-code-reference",
        }
    }
}

/// One dimensioned context-of-use coordinate.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct WfContext {
    /// Stable query-axis name.
    pub name: &'static str,
    /// Coherent-SI dimensions.
    pub dims: Dims,
    /// Inclusive lower bound.
    pub lo: f64,
    /// Inclusive upper bound.
    pub hi: f64,
}

/// One receipt-bound Wright Flyer corpus row.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct WfLevelACase {
    /// Stable corpus dataset id.
    pub id: &'static str,
    /// Human-readable title.
    pub title: &'static str,
    /// Coverage family.
    pub family: WfFamily,
    /// Analytic or cross-code referee.
    pub kind: WfKind,
    /// Stable scalar metric name.
    pub metric: &'static str,
    /// Metric dimensions.
    pub metric_dims: Dims,
    /// The measure-then-pinned executed value.
    pub reference_value_si: f64,
    /// Absolute tolerance.
    pub atol: f64,
    /// Relative tolerance.
    pub rtol: f64,
    /// The receipt digest (hex) binding this row to its executed
    /// battery evidence.
    pub receipt_digest: &'static str,
    /// Complete query context.
    pub context: &'static [WfContext],
}

const fn ctx(name: &'static str, dims: Dims, lo: f64, hi: f64) -> WfContext {
    WfContext { name, dims, lo, hi }
}

/// The registered Wright Flyer rows.
#[must_use]
pub fn wf_level_a_cases() -> &'static [WfLevelACase] {
    const CASES: &[WfLevelACase] = &[
        WfLevelACase {
            id: "wf-a-image-onplane-residual",
            title: "V-06a: hybrid-wake on-plane normal-velocity cancellation",
            family: WfFamily::GroundImages,
            kind: WfKind::AnalyticReference,
            metric: "wall-normal-residual-rel",
            metric_dims: DIMENSIONLESS,
            reference_value_si: 0.0,
            atol: 1e-9,
            rtol: 0.0,
            receipt_digest: "cebf414b1ba1b5086b71afb372ab0b3f8bebf39f056e066d34f334b3d827f503",
            context: &[ctx("ground-clearance-m", WF_LENGTH_DIMS, 0.5, 3.0)],
        },
        WfLevelACase {
            id: "wf-a-wakeref-wagner-start",
            title: "V-08b1: dense referee Wagner-class starting ratio",
            family: WfFamily::UnsteadyReferee,
            kind: WfKind::CrossCodeReference,
            metric: "wagner-start-ratio",
            metric_dims: DIMENSIONLESS,
            reference_value_si: 0.918_041_760_939_216,
            atol: 0.0,
            rtol: 1e-9,
            receipt_digest: "289bbe393d8b79dfddcbc92becfe4e42b0eb1a501f9561601a165926a21ce1f1",
            context: &[ctx("freestream-mps", WF_SPEED_DIMS, 13.0, 13.0)],
        },
        WfLevelACase {
            id: "wf-a-farfield-v10-shape",
            title: "V-10: hybrid-wake buildup shape vs the dense reference",
            family: WfFamily::HybridWake,
            kind: WfKind::CrossCodeReference,
            metric: "v10-shape-rms",
            metric_dims: DIMENSIONLESS,
            reference_value_si: 0.128_472_169_394_044_4,
            atol: 0.05,
            rtol: 0.0,
            receipt_digest: "bc6c175b04fe249f5747bd286aa11a994f02c696bbfd524a9730e6d6a073b448",
            context: &[ctx("overlap-ticks", DIMENSIONLESS, 237.0, 237.0)],
        },
        WfLevelACase {
            id: "wf-a-blade-cover-capsules",
            title: "V-20: blade-collision cover certificate at the registered rotor",
            family: WfFamily::SweptContact,
            kind: WfKind::AnalyticReference,
            metric: "capsules-per-blade",
            metric_dims: DIMENSIONLESS,
            reference_value_si: 16.0,
            atol: 0.0,
            rtol: 0.0,
            receipt_digest: "6a7422ddc89749610ebc58a0aad41fb75522dd58362e0be2efc73ea4d1d7c28d",
            context: &[ctx("geometry-uncertainty-m", WF_LENGTH_DIMS, 0.01, 0.01)],
        },
        WfLevelACase {
            id: "wf-a-harness-worst-rel",
            title: "E10.1: pinned referee-harness worst discrepancy",
            family: WfFamily::RefereeHarness,
            kind: WfKind::AnalyticReference,
            metric: "harness-worst-abs-rel",
            metric_dims: DIMENSIONLESS,
            reference_value_si: 0.387_084_355_634_945_5,
            atol: 0.0,
            rtol: 1e-9,
            receipt_digest: "289bbe393d8b79dfddcbc92becfe4e42b0eb1a501f9561601a165926a21ce1f1",
            context: &[ctx("pinned-alpha-rad", DIMENSIONLESS, 0.03, 0.07)],
        },
        WfLevelACase {
            id: "wf-a-h07-slope-recovery",
            title: "H-07: synthetic-truth slope recovery through the signed pipeline",
            family: WfFamily::HistoricalInference,
            kind: WfKind::AnalyticReference,
            metric: "posterior-slope-recovery",
            metric_dims: DIMENSIONLESS,
            reference_value_si: -0.7,
            atol: 0.1,
            rtol: 0.0,
            receipt_digest: "b17cb8e3620e8fc8c7134cec9a4176d5603ad3a3a8708edae987c1ed18bc8d46",
            context: &[ctx("lofo-folds", DIMENSIONLESS, 4.0, 4.0)],
        },
    ];
    CASES
}
