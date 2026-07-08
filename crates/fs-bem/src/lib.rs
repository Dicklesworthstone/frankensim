//! fs-bem — Laplace BEM panel methods (plan §8.3 [F], bead tfz.20):
//! potential-flow screening for exterior aerodynamics, O(N)-class via
//! fs-fmm. Inviscid honesty labels apply everywhere: this is the
//! ornithoid flagship's WIDE-SEARCH stage, not a viscous truth source.
//!
//! Layer: L3.
//! - [`panel3d`]: 3D exterior flow — constant source panels on
//!   fs-rep-mesh surfaces, collocation Neumann conditions, GMRES with
//!   the FMM-accelerated gradient matvec (three Chebyshev passes for
//!   the vector kernel, dotted with target normals) and the dense
//!   direct path as the oracle; the sphere's analytic surface speed is
//!   the G2 gate, single-layer reciprocity the G0 identity.
//! - [`panel2d`]: 2D Hess–Smith airfoils — constant sources per panel
//!   plus one bound vortex, the KUTTA condition closing the system;
//!   thin-airfoil lift slope as the reference band; the ADJOINT
//!   (one transposed dense solve) is the committed gradient path,
//!   FD-gated.
//! - [`wake2d`]: unsteady free wakes — point-vortex sheets shed at the
//!   trailing edge with Kutta-determined strength (Kelvin circulation
//!   conservation), convected by the regularized induced flow; the
//!   impulsive-start fixture shows the Wagner-like circulation
//!   transient and stable roll-up.

pub mod panel2d;
pub mod panel3d;
pub mod wake2d;

pub use panel2d::{Airfoil2d, PanelSolution2d, naca4_symmetric};
pub use panel3d::{SpherePanels, solve_exterior};
pub use wake2d::{WakeSim, WakeStep};

/// Crate version, re-exported for provenance stamping.
pub const VERSION: &str = env!("CARGO_PKG_VERSION");

#[cfg(test)]
mod tests {
    #[test]
    fn version_is_stamped() {
        assert!(!super::VERSION.is_empty());
    }
}
