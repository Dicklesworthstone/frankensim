//! fs-lattice — lattice/infill optimization (plan §9.5 [F], bead
//! 7tv.14). Layer: L4 (ASCENT).
//!
//! Three stages, 2D smoke tier and honest about it: periodic
//! UNIT-CELL HOMOGENIZATION (effective elasticity from cell problems,
//! audited against Voigt bounds, dilute-limit analytics, symmetry and
//! positive-definiteness — physics gates, not vibes), GRADED
//! macro-optimization through the fitted property manifold with
//! adjoint gradients, and DE-HOMOGENIZATION back to explicit
//! micro-geometry re-analyzed at full resolution — with the
//! separation-of-scales validity flag doing real work when gradation
//! is too sharp. 3D TPMS families (gyroid stiffness curves vs
//! literature), fs-fab manufacturability audits, and FrankenNetworkx
//! conforming lattice generation are recorded successors.

pub mod dehomog;
pub mod graded;
pub mod homogenize;

pub use dehomog::{HoleArray, fullres_compliance};
pub use graded::{GradedDesign, PropertyFit, graded_compliance_opt};
pub use homogenize::{EffectiveTensor, Homogenizer, UnitCell, voigt_bound};

/// Crate version, re-exported for provenance stamping.
pub const VERSION: &str = env!("CARGO_PKG_VERSION");
