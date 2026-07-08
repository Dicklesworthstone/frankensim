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
//! here as fs-topopt. CutFEM-octree execution (topology evolving with
//! ZERO remeshing) is the marquee follow-up lane recorded on the
//! bead.

pub mod eigenfreq;
pub mod elasticity;
pub mod filter;
pub mod oc;
pub mod pipeline;
pub mod robust;
pub mod stress;

pub use eigenfreq::{
    eigenfrequency_objective, eigenvalue_gradient, lowest_eigenpairs, mass_interp, smooth_min,
};
pub use elasticity::DensityElasticity;
pub use filter::{DensityFilter, heaviside, heaviside_derivative};
pub use oc::{OcReport, optimality_criteria};
pub use pipeline::{DesignPipeline, SimpParams};
pub use robust::{RobustPipeline, RobustReport, ThreeField, robust_optimality_criteria};
pub use stress::{StressReport, stress_aggregate, von_mises, von_mises_derivative};

/// Crate version, re-exported for provenance stamping.
pub const VERSION: &str = env!("CARGO_PKG_VERSION");
