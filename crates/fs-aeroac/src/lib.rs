//! # fs-aeroac — low-Mach aeroacoustic source models
//!
//! Bead frankensim-fsim-aeroacoustic-sources-9ok02 (IN PROGRESS —
//! this crate is the first slice): the acoustic-analogy half of the
//! honest hybrid. fs-lbm computes incompressible unsteady base flows
//! (its honest contract: NO acoustic radiation); this crate turns
//! extracted source terms into radiated spectra via Curle's analogy
//! in the FREQUENCY domain, where the 2D free-space Green's function
//! is the outgoing Hankel function — no time-tail convolution, no
//! 3D-formula-over-2D-data trap.
//!
//! HONEST SCOPE (the load-bearing law, pinned as data and by test):
//! outputs are RELATIVE spectral shapes and scaling laws, never
//! absolute SPL — see [`SCOPE_STATEMENT`].
//!
//! Modules:
//! - [`bessel`]: self-contained deterministic J0/J1/Y0/Y1 and the
//!   outgoing Hankel functions, built exclusively from `fs_math::det`
//!   primitives (series with double-double accumulation below the
//!   crossover, programmatically generated asymptotic sums above it —
//!   NO transcribed coefficient tables to mis-copy), certified by the
//!   exact Wronskian identity and a cross-implementation oracle.
//! - [`curle2d`]: 2D frequency-domain Curle dipole radiation
//!   (e^{-i omega t} convention, matching the workspace acoustics
//!   doctrine).
//! - [`bickley`]: Rayleigh-equation instability oracle for the
//!   Bickley jet `U = sech^2(y)` — the analytic reference the fs-lbm
//!   jet runs are validated against. Its own pins are SELF-VERIFIED
//!   exact eigenmodes (machine-zero ODE residuals re-proven per run).
//! - [`regime`]: recorded 2D CentralMoment tonal-lock landscape (bead
//!   l011o), the inverse-cascade refusal to catalog those spectra as
//!   broadband flute noise, a spectral-flatness measurement, and a
//!   D3Q19 operator smoke that does not mint a 3D broadband table.

pub mod bessel;
pub mod bickley;
pub mod curle2d;
pub mod jetlab;
pub mod noisetable;
pub mod regime;

/// The honest-scope statement every exported artifact embeds (the
/// bead's marketing-mutation guard asserts its presence): 2D line
/// sources radiate with cylindrical spreading and the analogy here is
/// a SHAPE and SCALING authority only.
pub const SCOPE_STATEMENT: &str = "fs-aeroac outputs are relative spectral shapes and scaling laws \
from a 2D incompressible base flow via Curle's analogy; they are NOT absolute SPL predictions. \
2D line-source (cylindrical-spreading) results additionally need a 2D-to-3D span correction \
before any comparison to measured levels.";

/// Typed refusal.
#[derive(Debug, Clone, PartialEq)]
pub enum AeroacError {
    /// Non-finite input.
    NonFinite {
        /// Where.
        what: &'static str,
    },
    /// Physically invalid parameter (non-positive wavenumber, zero
    /// radius, empty grid...).
    InvalidParameter {
        /// What.
        what: &'static str,
    },
    /// An iterative solve did not converge; the partial state is
    /// refused, not returned.
    NotConverged {
        /// Which solve.
        what: &'static str,
        /// Final residual magnitude.
        residual: f64,
    },
}

impl core::fmt::Display for AeroacError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            AeroacError::NonFinite { what } => write!(f, "non-finite input: {what}"),
            AeroacError::InvalidParameter { what } => write!(f, "invalid parameter: {what}"),
            AeroacError::NotConverged { what, residual } => {
                write!(f, "not converged: {what} (residual {residual:e})")
            }
        }
    }
}

impl std::error::Error for AeroacError {}
