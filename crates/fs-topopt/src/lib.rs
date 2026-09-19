//! fs-topopt — density-based topology optimization (plan §9.5 [S]).
//! Layer: L4 ASCENT.
//!
//! SIMP with the modern hygiene stack: Helmholtz PDE filtering for
//! mesh-independent length-scale control (REUSING the Poisson
//! machinery — one solver, two jobs), Heaviside projection with β
//! continuation (crisp designs without premature lock-in),
//! penalization continuation, EXACT chain-rule sensitivities through
//! the whole density pipeline (SIMP ∘ projection ∘ filter — every
//! stage linear or with a closed-form derivative, verified against
//! finite differences at multiple continuation stages), and the
//! classic optimality-criteria driver for compliance/volume
//! (fs-ascent's augmented Lagrangian is the general constrained
//! path; OC is the documented default for this problem class).
//!
//! NAMING: the plan's atlas used "fs-topo" for this stack; that crate
//! name now carries the L2 topology-CERTIFICATE machinery
//! (persistence, cubical homology), so the optimization stack lives
//! here as fs-topopt. The feature-gated `sdf3` path connects 3-D raw
//! implicit cuts to density optimization on Cartesian and adaptive backgrounds.
//! `sdf3_goal` drives refinement with enriched goal-weighted residuals.
//! Composed continuum certification remains a separate obligation.

pub mod continuation;
pub mod control;
pub mod eigenfreq;
pub mod elasticity;
pub mod filter;
pub mod gradient_check;
#[cfg(feature = "cutfem-marquee")]
pub mod marquee;
pub mod multi_load;
pub mod oc;
pub mod pipeline;
pub mod robust;
#[cfg(feature = "cutfem-marquee")]
pub mod sdf3;
#[cfg(feature = "cutfem-marquee")]
pub mod sdf3_goal;
pub mod stress;

pub use continuation::{
    ContinuationStageReport, ContinuationTermination, MultiLoadContinuationReport,
    controlled_gradient_checked_multi_load_continuation, controlled_multi_load_continuation,
    multi_load_continuation,
};
pub use control::{EvaluationStop, SolveBudget, SolveControl, SolveProgress, SolveWork};
pub use eigenfreq::{
    eigenfrequency_objective, eigenvalue_gradient, lowest_eigenpairs, mass_interp, smooth_min,
};
pub use elasticity::DensityElasticity;
pub use filter::{DensityFilter, heaviside, heaviside_derivative};
pub use gradient_check::{
    GradientCheckOptions, GradientDirection, GradientProbe, MultiLoadGradientCheck,
    controlled_multi_load_gradient_check,
};
pub use multi_load::{
    MultiLoadOcIteration, MultiLoadOcOptions, MultiLoadOcReport, MultiLoadOcTermination,
    controlled_multi_load_optimality_criteria, multi_load_optimality_criteria,
};
pub use oc::{OcReport, optimality_criteria};
pub use pipeline::{DesignPipeline, SimpParams};
pub use robust::{RobustPipeline, RobustReport, ThreeField, robust_optimality_criteria};
pub use stress::{StressReport, stress_aggregate, von_mises, von_mises_derivative};

/// Crate version, re-exported for provenance stamping.
pub const VERSION: &str = env!("CARGO_PKG_VERSION");
