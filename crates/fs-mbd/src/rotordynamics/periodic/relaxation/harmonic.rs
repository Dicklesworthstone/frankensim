//! Globally phased harmonic and Fourier loads on every original memory state.

use super::super::{
    CycleForcing, RotorError, RotorFloquetBudget, RotorFourierForce, RotorHarmonicForce,
    RotorPeriodicSegment,
};
use super::{PreparedRotorMaxwellCycle, prepare};
use crate::rotordynamics::relaxation::MaxwellDampingBranch;

/// Prepare a complete hereditary rotor cycle under continuous sinusoidal loads.
///
/// Each segment supplies original generalized force/torque cosine and sine
/// amplitudes relative to the COMMON cycle origin. The positive integer order
/// is in 1..=MAX_ROTOR_FORCE_HARMONIC and sets omega=2*pi*order/period.
/// Phase advances through every coefficient or authored amplitude switch.
/// For exp(+i omega t) phasors, cosine is the real part and sine the negative
/// imaginary part. Source rotation direction is independent of the fixed signed
/// shaft speed; synchronous forcing requires the caller's period/order to match.
///
/// All original q, q_dot and individual Maxwell stresses remain in the same
/// generator, affine map and homogeneous spectrum. Direct mechanical forcing
/// is extended with zero stress-source entries; the reciprocal original
/// generator transfers motion into material memory. No equivalent damping,
/// strain-reset assumption, material-pole cancellation, or held-sinusoid
/// approximation replaces that law. Supplied stiffness is the relaxed K.
///
/// The existing normalized harmonic exponential owns source rotation and
/// coupling. Two auxiliary coordinates carry no physical energy and never
/// enter the returned Floquet spectrum. The bounded largest exponential has
/// 2*n+branches+2 coordinates (at most 130); full-state byte/work admission and
/// original force/memory residuals apply before publication. Warm advances and
/// particular-orbit queries reuse the same prepared maps as held forcing.
///
/// This is a finite segmented model, not a speed ramp, smooth-modulation error
/// certificate, nonlinear resonance analysis or full physical-work receipt.
/// Cancellation or numerical refusal returns no partially prepared cycle.
#[allow(clippy::too_many_arguments)]
pub fn prepare_maxwell_harmonic_cycle(
    spin_rad_s: f64,
    segments: &[RotorPeriodicSegment<'_>],
    branches: &[MaxwellDampingBranch],
    forces: &[RotorHarmonicForce<'_>],
    harmonic_order: usize,
    relative_generator_residual_tolerance: f64,
    budget: RotorFloquetBudget,
    cancelled: impl FnMut() -> bool,
) -> Result<PreparedRotorMaxwellCycle, RotorError> {
    prepare(
        spin_rad_s,
        segments,
        branches,
        None,
        CycleForcing::Harmonic {
            rows: forces,
            order: harmonic_order,
        },
        relative_generator_residual_tolerance,
        budget,
        true,
        cancelled,
    )
}

/// Prepare every original mechanical/material state under DC plus Fourier loads.
///
/// Existing RotorFourierForce rows retain original force/torque units and the
/// same common cycle origin as prepare_maxwell_harmonic_cycle. Orders are
/// strictly increasing in 1..=1024, with at most 32 tones. Empty orders give the
/// held-load cycle. Amplitude and coefficient switches never reset source
/// phase, displacement, velocity, or any retained branch stress.
///
/// One full-mass factor and one augmented physical generator per segment serve
/// all tones. The existing source owner integrates each tone and sums offsets
/// in double-double arithmetic, retaining only one homogeneous transition.
/// All 2*n+branches physical multipliers survive; source oscillators add none.
/// Cold work is linear in tone count and is completely admitted, including
/// silent declared tones. Retained maps and warm advance/orbit work do not grow
/// with tone count. No dynamic-stiffness inverse is needed at resonance.
///
/// Coefficients, branch laws, M/G and signed spin obey prepare_maxwell_cycle's
/// original restrictions. This API does not infer a truncation error for an
/// arbitrary measured waveform or replace a separate physical-work audit.
#[allow(clippy::too_many_arguments)]
pub fn prepare_maxwell_fourier_cycle(
    spin_rad_s: f64,
    segments: &[RotorPeriodicSegment<'_>],
    branches: &[MaxwellDampingBranch],
    forces: &[RotorFourierForce<'_>],
    harmonic_orders: &[usize],
    relative_generator_residual_tolerance: f64,
    budget: RotorFloquetBudget,
    cancelled: impl FnMut() -> bool,
) -> Result<PreparedRotorMaxwellCycle, RotorError> {
    prepare(
        spin_rad_s,
        segments,
        branches,
        None,
        CycleForcing::Fourier {
            rows: forces,
            orders: harmonic_orders,
        },
        relative_generator_residual_tolerance,
        budget,
        true,
        cancelled,
    )
}
