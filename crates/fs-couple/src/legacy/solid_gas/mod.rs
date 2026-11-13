//! Finite solid/vapor exchange with reciprocal gas and solid recoil.
//!
//! The original material law supplies the molecular mass, momentum and energy
//! streams. Both finite momentum stores and both original caloric stores consume
//! that same exchange. Latent energy and molecular translational energy are
//! already included in the streams and are never added as extra heat.
//!
//! These are uniform-body, fixed-volume SOURCE operators. The supplied surface
//! is a reference patch: spatial gas transport, translation of the geometry,
//! rotation, pressure-volume work and surface remeshing retain separate owners.
//! Molecular reflection already supplies the pressure load on this patch; an
//! additional pressure boundary for the same molecules would double-count it.

/// One inventory-bounded, energy-conserving molecular source transaction.
pub mod kinetic;

/// Error-controlled repeated exchange on the same original physical owners.
pub mod trajectory;
