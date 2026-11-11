//! Interior acoustics on original spatial complexes (music bead M2.4): see [`helmholtz`] and [`axisymmetric`].

pub mod helmholtz;

/// Resolved velocity, pressure and temperature in viscous acoustic gaps.
pub mod thermoviscous;

/// Source-trace Galerkin coupling between pressure and resolved thermoviscous FE.
pub mod hybrid;

/// Full m=0 meridional pressure fields with exact swept volume and surface forms.
pub mod axisymmetric;

/// The mixed componentwise/normwise row gate of the refined solves.
pub(crate) mod backward_error;
