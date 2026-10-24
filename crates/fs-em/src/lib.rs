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
//!
//! Explicit geometry/frame correspondence transfers the computed field into
//! existing moving-conductor/current/heat/flow owners. That coupling is one-way
//! in the low-induction regime: their induced currents do not update this field.
//! `planar::eddy` additionally retains stationary planar induction history and
//! spatial Ohmic heat, including floating conductor current-return constraints.
//! The homogeneous winding and induction owners refuse permanent-magnet bias
//! until their source/storage trajectories account for it explicitly. There is
//! no nonlinear ferromagnetism, hysteresis, irreversible demagnetization, moving
//! conductor induction feedback or magnetization-force claim.

pub mod planar;
pub mod transfer;
pub mod winding;
