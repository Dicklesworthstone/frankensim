//! fs-dfo — derivative-free optimization engines (plan §9.3): CMA-ES in
//! its information-geometric form (natural-gradient flow on the Gaussian
//! family — the framing that DICTATES the step-size/covariance couplings
//! and buys invariance properties by construction), BIPOP restarts, and
//! a Nelder–Mead polish baseline.
//!
//! Layer: L4 ASCENT. Engines are IR-agnostic (closure objectives) by
//! design — routing from the fs-opt problem IR is a small wiring bead
//! once that crate stabilizes (deliberate collision avoidance; see bead
//! 7tv.4's comment trail).
//!
//! DETERMINISM: all sampling flows from keyed Philox streams; ordering
//! uses `total_cmp` or explicitly versioned deterministic tie policies.
//! In particular, NSGA-III equal-rank mating keeps the first ordered
//! uniform draw instead of treating the lower live population index as
//! fitter. The whole evolution is a pure function of the seed (bitwise
//! rerun-tested and golden-hashed).

pub mod cma;
pub mod dro;
pub mod moo;
pub mod neldermead;
pub mod ot;
pub mod steer;
pub mod wfg;

pub use fs_robust::{EmpiricalCvarReport, RobustError, empirical_cvar};

pub use cma::{
    BIPOP_ADMISSION_SCHEMA_VERSION, BIPOP_EVALUATION_SCHEMA_VERSION, BIPOP_REPORT_SCHEMA_VERSION,
    BIPOP_RESTART_SCHEMA_VERSION, BIPOP_ROOT_IDENTITY_KIND, BIPOP_ROOT_IDENTITY_SCHEMA_DECLARATION,
    BIPOP_ROOT_IDENTITY_SCHEMA_VERSION, BIPOP_STUDY_IDENTITY_DOMAIN,
    BIPOP_STUDY_IDENTITY_SCHEMA_DECLARATION, BIPOP_STUDY_IDENTITY_SCHEMA_VERSION,
    BIPOP_TRACE_IDENTITY_DOMAIN, BIPOP_TRACE_IDENTITY_SCHEMA_DECLARATION,
    BIPOP_TRACE_IDENTITY_SCHEMA_VERSION, BipopAdmission, BipopError, BipopEvaluationRecord,
    BipopLane, BipopLedgerError, BipopReplayAdmissionError, BipopReport, BipopRestartRecord,
    BipopRootInputs, BipopStudyAdmissionError, BipopStudyIdentity, BipopTraceIdentity, CmaParams,
    CmaReport, CmaStopReason, admit_bipop, bipop_cmaes, cmaes, try_bipop_cmaes,
};
pub use dro::{DroReport, wasserstein_worst_case};
pub use moo::{
    HvArchive, Individual, MoeadParams, NsgaParams, crowding_distance, das_dennis, dominates,
    hypervolume, knee_point, mc_hypervolume, moead, non_dominated_sort, nsga2, nsga3,
};
pub use neldermead::nelder_mead;
pub use ot::{OtReport, cost_sq_1d, monotone_cost_1d, sinkhorn};
pub use steer::{SteerEvent, SteeredStudy, StudyState};

/// Crate version, re-exported for provenance stamping.
pub const VERSION: &str = env!("CARGO_PKG_VERSION");
