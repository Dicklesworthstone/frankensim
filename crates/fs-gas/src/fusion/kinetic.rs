//! Isotropic alpha-particle energy distributions on the original fusion sources.
//!
//! The bounded energy grid retains particle number and kinetic energy at birth.
//! Its collision kernel evolves energy diffusion as well as slowing, with
//! separate reciprocal electron and heavy-particle heat accounts. Particles
//! remain kinetic at low energy: there is no terminal energy cutoff, implicit
//! thermalization, or invented removal of their pressure and heat capacity.
//!
//! Both energy boundaries are reflecting. Finite energy support and grid
//! resolution are explicit numerical approximations; detailed balance on a
//! fixed grid is not a claim of exact continuum Maxwellian quadrature.

mod collision;
mod grid;

pub use grid::AlphaEnergyGrid;

pub use super::network::kinetic::{
    KineticAlphaBalance, KineticAlphaTransport, KineticFusionResult, KineticFusionState,
    advance_kinetic_fusion,
};

pub(crate) use collision::AlphaCollisionOperator;

/// Fixed upper bound on success-path collision workspace.
pub(crate) const MAX_ALPHA_BINS: usize = 256;
