//! Level-A thermal analytic references and manufactured-solution targets.
//!
//! These rows define exact or textbook-limit reference values and G1 order
//! targets. They are deliberately separate from solver observations: until a
//! consuming thermal kernel retains a refinement ladder or comparison receipt,
//! every row remains reference-only and every corpus query remains a numerical
//! no-claim.

use fs_qty::Dims;

/// Retained tab-separated manifest backing every Level-A thermal row.
pub(crate) const THERMAL_LEVEL_A_MANIFEST: &[u8] =
    include_bytes!("../../../data/vv-corpus/thermal-level-a/thermal-level-a-v1.tsv");

/// Coherent-SI dimensions of length.
pub const LENGTH_DIMS: Dims = Dims([1, 0, 0, 0, 0, 0]);
/// Coherent-SI dimensions of temperature.
pub const TEMPERATURE_DIMS: Dims = Dims([0, 0, 0, 1, 0, 0]);
/// Coherent-SI dimensions of heat flux, W/m².
pub const HEAT_FLUX_DIMS: Dims = Dims([0, 1, -3, 0, 0, 0]);
/// Coherent-SI dimensions of thermal conductance, W/K.
pub const THERMAL_CONDUCTANCE_DIMS: Dims = Dims([2, 1, -3, -1, 0, 0]);
/// Coherent-SI dimensions of thermal resistance, K/W.
pub const THERMAL_RESISTANCE_DIMS: Dims = Dims([-2, -1, 3, 1, 0, 0]);

/// Thermal reference family used for coverage accounting.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum ThermalLevelAFamily {
    /// Closed-form steady conduction, including mixed boundaries.
    SteadyConduction,
    /// Straight-fin efficiency.
    Fin,
    /// Lumped-capacitance transient limit.
    LumpedTransient,
    /// Fully developed laminar duct limiting values.
    ConvectionLimit,
    /// Closed-form surface-radiation geometry.
    Radiation,
    /// Series thermal-contact resistance.
    Contact,
    /// Manufactured primal convergence target.
    ManufacturedPrimal,
    /// Manufactured adjoint convergence target.
    ManufacturedAdjoint,
}

impl ThermalLevelAFamily {
    /// Stable manifest spelling.
    #[must_use]
    pub const fn name(self) -> &'static str {
        match self {
            Self::SteadyConduction => "steady-conduction",
            Self::Fin => "fin",
            Self::LumpedTransient => "lumped-transient",
            Self::ConvectionLimit => "convection-limit",
            Self::Radiation => "radiation",
            Self::Contact => "contact",
            Self::ManufacturedPrimal => "mms-primal",
            Self::ManufacturedAdjoint => "mms-adjoint",
        }
    }
}

/// Whether a row is an analytic value or an unexecuted G1 order target.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ThermalLevelAKind {
    /// Closed-form or canonical limiting value.
    AnalyticReference,
    /// Theoretical convergence-order target awaiting a retained ladder.
    ManufacturedTarget,
}

impl ThermalLevelAKind {
    /// Stable manifest spelling.
    #[must_use]
    pub const fn name(self) -> &'static str {
        match self {
            Self::AnalyticReference => "analytic-reference",
            Self::ManufacturedTarget => "manufactured-target",
        }
    }
}

/// Acceptance rule attached to a Level-A row.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum ThermalLevelAAcceptance {
    /// Absolute/relative comparison envelope for a scalar reference.
    Tolerance {
        /// Absolute tolerance in the metric's coherent SI units.
        atol: f64,
        /// Dimensionless relative tolerance.
        rtol: f64,
    },
    /// Two-sided G1 order gate around a theoretical slope.
    OrderGate {
        /// Expected log-error/log-mesh-size slope.
        theoretical: f64,
        /// Maximum absolute slope deviation.
        tolerance: f64,
    },
}

/// One dimensioned context-of-use coordinate.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ThermalLevelAContext {
    /// Stable query-axis name.
    pub name: &'static str,
    /// Coherent-SI dimensions.
    pub dims: Dims,
    /// Inclusive lower bound.
    pub lo: f64,
    /// Inclusive upper bound.
    pub hi: f64,
}

/// One reference-only Level-A thermal case.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ThermalLevelACase {
    /// Stable corpus dataset id.
    pub id: &'static str,
    /// Human-readable title.
    pub title: &'static str,
    /// Coverage family.
    pub family: ThermalLevelAFamily,
    /// Analytic value or G1 target.
    pub kind: ThermalLevelAKind,
    /// Stable scalar metric name.
    pub metric: &'static str,
    /// Metric dimensions.
    pub metric_dims: Dims,
    /// Frozen analytic value or theoretical order.
    pub reference_value_si: f64,
    /// Exact formula/target semantics retained in the manifest.
    pub formula: &'static str,
    /// Comparison or order envelope.
    pub acceptance: ThermalLevelAAcceptance,
    /// Complete query context.
    pub context: &'static [ThermalLevelAContext],
    /// Explicit reason this definition does not yet carry solver evidence.
    pub no_claim_reason: &'static str,
}

const REFERENCE_ONLY: &str =
    "reference definition only; no tracked thermal-kernel comparison receipt is bound";
const TARGET_ONLY: &str =
    "G1 target only; no retained refinement ladder from a thermal kernel is bound";

const fn context(name: &'static str, dims: Dims, lo: f64, hi: f64) -> ThermalLevelAContext {
    ThermalLevelAContext { name, dims, lo, hi }
}

const fn analytic(
    id: &'static str,
    title: &'static str,
    family: ThermalLevelAFamily,
    metric: &'static str,
    metric_dims: Dims,
    reference_value_si: f64,
    formula: &'static str,
    context: &'static [ThermalLevelAContext],
) -> ThermalLevelACase {
    ThermalLevelACase {
        id,
        title,
        family,
        kind: ThermalLevelAKind::AnalyticReference,
        metric,
        metric_dims,
        reference_value_si,
        formula,
        acceptance: ThermalLevelAAcceptance::Tolerance {
            atol: 1.0e-12,
            rtol: 1.0e-12,
        },
        context,
        no_claim_reason: REFERENCE_ONLY,
    }
}

const fn mms(
    id: &'static str,
    title: &'static str,
    family: ThermalLevelAFamily,
    theoretical: f64,
    context: &'static [ThermalLevelAContext],
) -> ThermalLevelACase {
    ThermalLevelACase {
        id,
        title,
        family,
        kind: ThermalLevelAKind::ManufacturedTarget,
        metric: "observed-l2-order",
        metric_dims: Dims::NONE,
        reference_value_si: theoretical,
        formula: "G1 L2 slope target with absolute deviation at most 0.2",
        acceptance: ThermalLevelAAcceptance::OrderGate {
            theoretical,
            tolerance: 0.2,
        },
        context,
        no_claim_reason: TARGET_ONLY,
    }
}

static P1_MMS_CONTEXT: [ThermalLevelAContext; 2] = [
    context("element-degree", Dims::NONE, 1.0, 1.0),
    context("mesh-size-m", LENGTH_DIMS, 0.0625, 0.25),
];

static P2_MMS_CONTEXT: [ThermalLevelAContext; 2] = [
    context("element-degree", Dims::NONE, 2.0, 2.0),
    context("mesh-size-m", LENGTH_DIMS, 0.0625, 0.25),
];

static THERMAL_LEVEL_A_CASES: [ThermalLevelACase; 19] = [
    analytic(
        "thermal-a-slab-dirichlet",
        "Planar slab with two prescribed temperatures",
        ThermalLevelAFamily::SteadyConduction,
        "outward-heat-flux",
        HEAT_FLUX_DIMS,
        4_000.0,
        "q = k * (T_hot - T_cold) / L; k=20 W/(m K), deltaT=40 K, L=0.2 m",
        &[context("slab-thickness-m", LENGTH_DIMS, 0.2, 0.2)],
    ),
    analytic(
        "thermal-a-slab-robin",
        "Planar slab with prescribed-hot and convective-cold boundaries",
        ThermalLevelAFamily::SteadyConduction,
        "outward-heat-flux",
        HEAT_FLUX_DIMS,
        5_000.0,
        "q = (T_hot - T_inf) / (L/k + 1/h); deltaT=100 K, L=0.1 m, k=10 W/(m K), h=100 W/(m2 K)",
        &[context("biot-number", Dims::NONE, 1.0, 1.0)],
    ),
    analytic(
        "thermal-a-slab-uniform-source",
        "Symmetric slab with uniform volumetric heating",
        ThermalLevelAFamily::SteadyConduction,
        "center-temperature-rise",
        TEMPERATURE_DIMS,
        12.5,
        "deltaT_center = qdot * L^2 / (8 k); qdot=100000 W/m3, L=0.1 m, k=10 W/(m K)",
        &[context("slab-thickness-m", LENGTH_DIMS, 0.1, 0.1)],
    ),
    analytic(
        "thermal-a-rectangle-linear",
        "Two-dimensional rectangular affine-temperature patch",
        ThermalLevelAFamily::SteadyConduction,
        "probe-temperature",
        TEMPERATURE_DIMS,
        320.0,
        "T(x,y) = 300 + 20 x + 40 y K at x=0.5 m, y=0.25 m",
        &[context("rectangle-aspect-ratio", Dims::NONE, 2.0, 2.0)],
    ),
    analytic(
        "thermal-a-cylinder-shell",
        "Axisymmetric cylindrical-shell conductance",
        ThermalLevelAFamily::SteadyConduction,
        "thermal-conductance",
        THERMAL_CONDUCTANCE_DIMS,
        135.970_804_254_815_8,
        "G = 2*pi*k*Lz/ln(ro/ri); k=15 W/(m K), Lz=1 m, ri=0.05 m, ro=0.1 m",
        &[context("radius-ratio", Dims::NONE, 2.0, 2.0)],
    ),
    analytic(
        "thermal-a-sphere-shell",
        "Spherical-shell conductance",
        ThermalLevelAFamily::SteadyConduction,
        "thermal-conductance",
        THERMAL_CONDUCTANCE_DIMS,
        18.849_555_921_538_76,
        "G = 4*pi*k/(1/ri - 1/ro); k=15 W/(m K), ri=0.05 m, ro=0.1 m",
        &[context("radius-ratio", Dims::NONE, 2.0, 2.0)],
    ),
    analytic(
        "thermal-a-fin-efficiency",
        "Adiabatic-tip straight-fin efficiency",
        ThermalLevelAFamily::Fin,
        "fin-efficiency",
        Dims::NONE,
        0.761_594_155_955_764_9,
        "eta = tanh(mL)/(mL) at mL=1",
        &[context("m-times-l", Dims::NONE, 1.0, 1.0)],
    ),
    analytic(
        "thermal-a-lumped-transient",
        "Lumped-capacitance one-time-constant decay",
        ThermalLevelAFamily::LumpedTransient,
        "normalized-temperature-excess",
        Dims::NONE,
        0.367_879_441_171_442_33,
        "theta/theta0 = exp(-t/tau) at t/tau=1; valid only in the declared small-Biot regime",
        &[
            context("biot-number", Dims::NONE, 0.0, 0.1),
            context("normalized-time", Dims::NONE, 1.0, 1.0),
        ],
    ),
    analytic(
        "thermal-a-duct-nu-cwt",
        "Fully developed circular-duct constant-wall-temperature limit",
        ThermalLevelAFamily::ConvectionLimit,
        "nusselt-number",
        Dims::NONE,
        3.66,
        "Nu = 3.66 for hydrodynamically and thermally fully developed laminar circular-duct flow with constant wall temperature",
        &[context("reynolds-number", Dims::NONE, 1.0, 2_300.0)],
    ),
    analytic(
        "thermal-a-duct-nu-chf",
        "Fully developed circular-duct constant-heat-flux limit",
        ThermalLevelAFamily::ConvectionLimit,
        "nusselt-number",
        Dims::NONE,
        4.36,
        "Nu = 4.36 for hydrodynamically and thermally fully developed laminar circular-duct flow with uniform wall heat flux",
        &[context("reynolds-number", Dims::NONE, 1.0, 2_300.0)],
    ),
    analytic(
        "thermal-a-parallel-plate-view-factor",
        "Infinite parallel-plate view factor",
        ThermalLevelAFamily::Radiation,
        "view-factor-12",
        Dims::NONE,
        1.0,
        "F12 = 1 for the infinite parallel-plate limiting geometry",
        &[context("gap-to-extent-ratio", Dims::NONE, 0.0, 0.0)],
    ),
    analytic(
        "thermal-a-contact-series",
        "Two-layer plus interface thermal resistance in series",
        ThermalLevelAFamily::Contact,
        "thermal-resistance",
        THERMAL_RESISTANCE_DIMS,
        0.3,
        "R = L1/(k1*A) + Rc + L2/(k2*A); terms are 0.1, 0.1, and 0.1 K/W",
        &[context(
            "interface-resistance-k-per-w",
            THERMAL_RESISTANCE_DIMS,
            0.1,
            0.1,
        )],
    ),
    mms(
        "thermal-a-mms-p1-dirichlet",
        "P1 isotropic Dirichlet thermal MMS target",
        ThermalLevelAFamily::ManufacturedPrimal,
        2.0,
        &P1_MMS_CONTEXT,
    ),
    mms(
        "thermal-a-mms-p2-dirichlet",
        "P2 isotropic Dirichlet thermal MMS target",
        ThermalLevelAFamily::ManufacturedPrimal,
        3.0,
        &P2_MMS_CONTEXT,
    ),
    mms(
        "thermal-a-mms-p1-anisotropic-nonlinear",
        "P1 anisotropic temperature-dependent conductivity MMS target",
        ThermalLevelAFamily::ManufacturedPrimal,
        2.0,
        &P1_MMS_CONTEXT,
    ),
    mms(
        "thermal-a-mms-p1-neumann",
        "P1 mixed-Neumann thermal MMS target",
        ThermalLevelAFamily::ManufacturedPrimal,
        2.0,
        &P1_MMS_CONTEXT,
    ),
    mms(
        "thermal-a-mms-p1-robin",
        "P1 Robin thermal MMS target",
        ThermalLevelAFamily::ManufacturedPrimal,
        2.0,
        &P1_MMS_CONTEXT,
    ),
    mms(
        "thermal-a-mms-p1-adjoint",
        "P1 heat-adjoint consistency-order target",
        ThermalLevelAFamily::ManufacturedAdjoint,
        2.0,
        &P1_MMS_CONTEXT,
    ),
    mms(
        "thermal-a-mms-p2-adjoint",
        "P2 heat-adjoint consistency-order target",
        ThermalLevelAFamily::ManufacturedAdjoint,
        3.0,
        &P2_MMS_CONTEXT,
    ),
];

/// Complete, stable Level-A thermal reference and target catalog.
#[must_use]
pub fn thermal_level_a_cases() -> &'static [ThermalLevelACase] {
    &THERMAL_LEVEL_A_CASES
}
