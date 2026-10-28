//! fs-uq — uncertainty quantification (plan §8.8). Layer: L4 ASCENT
//! (propagation WRAPS solvers the way ASCENT optimizers do, and its
//! risk outputs feed ASCENT's robust formulations; sitting at L4 also
//! reuses fs-bo's deterministic Φ⁻¹ for QMC normal germs instead of
//! duplicating the polynomial).
//!
//! Slice 1: Karhunen–Loève random fields with CAPTURED-VARIANCE
//! evidence, polynomial chaos by regression with known-answer gates,
//! the QMC workhorse with its MC advantage MEASURED, and multilevel
//! Monte Carlo with the telescoping sum AUDITED. Seismic machinery
//! (Kanai–Tajimi, CQC, IDA fragility) and e-process anytime-valid
//! stopping are the bead's split lanes.

pub mod adaptive;
pub mod anytime;
pub mod chance;
/// Bounded selection of complete experiments from local observation Jacobians.
pub mod experimental_design;
pub mod kl;
pub mod mlmc;
pub mod pce;
pub mod product_batch;
pub mod product_checkpoint;
pub mod product_compliance;
pub mod product_control_variate;
pub mod product_copula;
pub mod product_copula_qmc;
pub mod product_epistemic;
pub mod product_execution;
pub mod product_importance;
pub mod product_lhs;
pub mod product_mlmc;
pub mod product_mlmc_allocation;
pub mod product_mlmc_pilot;
pub mod product_multi;
pub mod product_pce;
pub mod product_pce_control;
pub mod product_pce_model;
pub mod product_plan;
pub mod product_qmc;
pub mod product_sensitivity;
pub mod product_shapley;
pub mod seismic;

pub use fs_robust::{EmpiricalCvarReport, RobustError, cvar, empirical_cvar};
pub use product_batch::{UqBatchExecution, UqBatchInput, UqBatchReport};
pub use product_checkpoint::UqCheckpointError;
pub use product_compliance::UqComplianceError;
pub use product_control_variate::{LinearControlEstimate, LinearControlVariate, UqControlError};
pub use product_copula::GaussianCopulaExecution;
pub use product_copula_qmc::GaussianCopulaQmcExecution;
pub use product_epistemic::{
    EpistemicCompliance, EpistemicConfig, EpistemicExecution, EpistemicReport, EpistemicStatus,
    EpistemicWitness,
};
pub use product_execution::UqExecution;
pub use product_importance::proposal::execution::{
    ImportancePilotRound, PilotImportanceConfig, PilotImportanceExecution, PilotImportanceReport,
    PilotImportanceStage,
};
pub use product_importance::proposal::{
    ImportanceFitConfig, ImportanceFitError, ImportanceProposalFit,
};
pub use product_importance::{
    ImportanceConfig, ImportanceEstimate, ImportanceExecution, ImportanceObservation,
    ImportanceReport,
};
pub use product_lhs::mean_control::{LhsLinearControlEstimate, LhsLinearControlVariate};
pub use product_lhs::{LhsEstimate, LhsExecution, LhsPlan, LhsReport};
pub use product_mlmc::{MlmcEstimate, MlmcExecution, MlmcExecutionReport, MlmcLevelReport};
pub use product_mlmc_allocation::{
    MlmcAllocation, MlmcAllocationError, MlmcAllocationPolicy, allocate_mlmc,
};
pub use product_mlmc_pilot::{
    PilotMlmcConfig, PilotMlmcExecution, PilotMlmcReport, PilotMlmcStage,
};
pub use product_multi::{MultiQoiAssessmentError, MultiQoiCovariance};
pub use product_multi::{MultiQoiExecution, MultiQoiReport, QoiDefinition, QoiReport};
pub use product_pce::{PceConfig, PceExecution, PceReport, PceSummary};
pub use product_pce_control::{
    PceControlError, PceMeanConfig, PceMeanControl, PceMeanEstimate, PceMeanExecution,
    PceMeanReport, PceMeanStage,
};
pub use product_pce_model::{PceBasis, ProductPceModel};
pub use product_plan::{
    CorrelationModel, ParameterUncertainty, PropagationMethod, UncertaintyKind, UqPlan,
    UqPropagator, UqResult, UqStatus,
};
pub use product_qmc::{QmcConfig, QmcEstimate, QmcExecution, QmcReport};
pub use product_sensitivity::{SobolEffect, SobolEstimate, SobolExecution, SobolReport};
pub use product_shapley::{ShapleyEffect, ShapleyEstimate, ShapleyExecution, ShapleyReport};

pub use adaptive::{AdaptiveReport, adaptive_mlmc};
pub use anytime::{AnytimeEstimate, estimate_probability_anytime};
pub use chance::{ChanceSolution, NoFeasibleChanceCandidate, chance_constrained_min};
pub use kl::{CovarianceKind, KlExpansion};
pub use mlmc::{MlmcReport, mlmc_estimate};
pub use pce::{PceModel, fit_pce};
pub use seismic::{
    AccelerationUnit, BilinearIdaModel, FragilityPoint, GroundMotionRecord,
    GroundMotionRecordError, GroundMotionSampling, IdaIntensity, KanaiTajimi, RecordedIda,
    bilinear_peak_ductility, cqc, ida_fragility, ida_fragility_from_records,
    ida_fragility_with_sampling, sdof_peak, srss,
};

/// Crate version, re-exported for provenance stamping.
pub const VERSION: &str = env!("CARGO_PKG_VERSION");
