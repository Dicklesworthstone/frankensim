//! Material-resolved electromagnetic field problems on the shared numerical owners.
//!
//! The first rung is linear, isotropic, fixed-state planar magnetostatics:
//! -div(mu^-1 grad(A_z)) = J_z, B = (d_y A_z, -d_x A_z), H = B/mu.
//! The mesh, weighted Hodge, sparse assembly and checked Krylov recurrence remain
//! owned by fs-rep-mesh, fs-feec, fs-sparse and fs-solver. Material admission is
//! shared with the existing source-derived circuit inductance.
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
//! There is no nonlinear ferromagnetism, permanent magnet, hysteresis, moving
//! conductor induction feedback, magnetization-force or magnetic parameter-work
//! claim.

pub mod planar;
pub mod transfer;
pub mod winding;
