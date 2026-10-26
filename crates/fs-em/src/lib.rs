//! Material-resolved electromagnetic field problems on the shared numerical owners.
//!
//! The first rung is linear, isotropic, fixed-state planar magnetostatics:
//! -div(mu^-1 grad(A_z)) = J_z, B = (d_y A_z, -d_x A_z), H = B/mu.
//! The mesh, weighted Hodge, sparse assembly and checked Krylov recurrence remain
//! owned by fs-rep-mesh, fs-feec, fs-sparse and fs-solver. Material admission is
//! shared with the existing source-derived circuit inductance.
//! Explicit permanent-magnet cells add source-backed reversible recoil
//! B=mu_recoil H+Br with specimen orientation and checked operating limits.
//! Their exact weak source terms reuse this field solve and sampling path;
//! stored energy uses the declared fixed-remanence reference.
//!
//! Distributed signed unit-current profiles extract reciprocal self/mutual
//! inductance from the actual field/material geometry. The full magnetic
//! storage block drives the existing Kirchhoff circuit DAE, with explicit
//! end/leakage data and circuit-state field reconstruction in the same frame.
//! Fixed-remanence windings additionally retain their zero-current field and
//! total-flux offset. Native circuit coordinates remain incremental linkage;
//! actual reconstructed fields and accepted Ohmic/circuit endpoints must stay
//! within every original recoil envelope.
//!
//! Explicit geometry/frame correspondence transfers the computed field into
//! existing moving-conductor/current/heat/flow owners. That coupling is one-way
//! in the low-induction regime: their induced currents do not update this field.
//! `planar::eddy` additionally retains stationary planar induction history and
//! spatial Ohmic heat, including floating conductor current-return constraints.
//! The stationary induction-history owner still refuses permanent-magnet bias
//! until its source/storage trajectory accounts for it explicitly. There is
//! no hysteresis, irreversible demagnetization, moving
//! conductor induction feedback or distributed magnetostriction claim.
//! `planar::nonlinear` supplies a separate fixed-state reversible H(B) solve,
//! with source-integrated energy and the same Hodge, gauge and transfer owners.
//! Its nonlinear fields can drive existing imposed-field flow consumers; they
//! do not acquire the linear winding/induction superposition assumptions.
//! `planar::eddy::nonlinear` evolves those same H(B) fields with the original
//! floating-conductor projection and spatial Joule heating. Its backward-Euler
//! numerical dissipation is explicit and is never relabeled as material heat.
//! [`force`] transfers a static field's closed exterior Maxwell-stress integral
//! to the existing mechanical wrench owner, including admitted enclosed recoil
//! magnets. Its contour remains in an induced-linear exterior medium; it does
//! not evolve geometry or magnetic work.
//! The dielectric field uses the same contour integral for electric force and
//! torque. Native circuit electrode potentials can reconstruct its original
//! spatial field and mechanical load without another field solve.

/// Spatial source-bound dielectric fields and geometry-derived capacitance.
pub mod electrostatic;
pub mod force;
pub mod planar;
pub mod transfer;
pub mod winding;
