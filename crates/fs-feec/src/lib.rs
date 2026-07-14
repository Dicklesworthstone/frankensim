//! fs-feec — exterior-calculus core (plan §8.1, Bet 3). Layer: L3
//! FLUX.
//!
//! Fields are COCHAINS on cell complexes; the exterior derivative d is
//! the EXACT integer incidence operator from fs-rep-mesh (dd = 0
//! bitwise is its contract invariant — inherited, not re-proven here);
//! ALL metric information and ALL approximation lives in the Hodge
//! stars and Whitney mass matrices. The payoff: grad→curl→div
//! identities hold to machine precision BY CONSTRUCTION, eliminating
//! algebraic-complex defects. That identity alone does not eliminate
//! spurious pressure/EM modes or checkerboarding: such formulation-level
//! claims additionally require a conforming subcomplex, a bounded commuting
//! projection, correct boundary and gauge treatment, admissible quadrature
//! and coefficients, and coercivity or inf-sup evidence.
//!
//! Element geometry (Jacobians, determinants, inverses, grams) runs
//! through fs-la's batched small-dense kernels — the layout consumer
//! those kernels were built for.

pub mod assembly;
pub mod betti;
pub mod cochain;
pub mod cohomology;
pub mod fixtures;
pub mod highorder;
pub mod hodge;
pub mod whitney;

pub use assembly::{incidence_to_csr, stiffness};
pub use betti::{betti_numbers, integer_rank};
pub use cochain::{Cochain, cell_count};
pub use cohomology::{HodgeParts, circulation, deflate_harmonics, harmonic_basis, hodge_decompose};
pub use fixtures::{kuhn_cube, masked_cube_grid, on_unit_cube_boundary, single_tet, two_tets};
pub use highorder::derham::TensorDeRham;
pub use highorder::hex::{TensorSpace, pcg_matfree};
pub use highorder::quad1d::{element_matrices, gauss_legendre, legendre, lobatto_shapes};
pub use highorder::simplex::SimplexSpace;
pub use highorder::vecfam::{
    DgSpace, Family, VecSpace, build_element, curl_matrix, dg_cell_dofs, div_matrix, grad_matrix,
    nedelec_entity_dofs, rt_entity_dofs, tri_quad3d,
};
pub use hodge::{galerkin_star, hodge_diagonal_barycentric};
pub use whitney::{
    ElementGeometry, deram0, deram1, deram2, deram3, element_geometry, mass_matrix, sort_parity,
};

/// Crate version, re-exported for provenance stamping.
pub const VERSION: &str = env!("CARGO_PKG_VERSION");
