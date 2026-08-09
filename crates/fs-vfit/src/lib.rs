//! # fs-vfit — passive rational approximation of frequency responses
//!
//! The bridge between offline frequency-domain results (BEM radiation
//! loads, TMM bore impedances, coupled-body mobilities) and
//! runtime-realizable filters: identify a compact stable pole-residue
//! model from tabulated `H(i*omega)` samples, CERTIFY passivity (grid +
//! Hamiltonian eigenvalue test) with convex residue repair when the raw
//! fit is active, and discretize by prewarped bilinear transform to
//! state-space and biquad-cascade forms.
//!
//! Two independent identification front ends — relaxed vector fitting
//! ([`vf`]) and the Loewner matrix pencil ([`loewner`]) — catch each
//! other's artifacts; the cross-check is part of the conformance
//! battery, not an afterthought.
//!
//! Sign conventions: time dependence `e^{+i*omega*t}` is NOT assumed
//! anywhere; everything is phrased on the Laplace axis `s = i*omega`
//! with real impulse responses enforced by conjugate-closed storage
//! ([`model::PoleTerm`]). Passivity is impedance-form positive
//! realness `Re H(i*omega) >= 0`.

pub mod discretize;
pub mod loewner;
pub mod model;
pub mod passivity;
pub mod vf;

pub use model::{PoleTerm, RationalModel, StateSpace};
pub use vf::{FitOptions, FitOutcome, FitReport, VfError, WeightPreset, vector_fit};
