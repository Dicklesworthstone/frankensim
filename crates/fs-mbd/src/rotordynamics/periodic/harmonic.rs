//! Continuous rotating loads in the same complete periodic rotor preparation.
use super::{
    AffineSegmentMap, CycleForcing, DenseExponentialBudget, MassSolve, PreparedRotorCycle,
    RotorError, RotorFloquetBudget, RotorPeriodicSegment, affine_segment_flow,
    exponential_with_budget, finite, poll, prepare_cycle, product_three,
};
use fs_math::{dd::Dd, det};

/// Bound on complete forcing revolutions in one declared coefficient cycle.
pub const MAX_ROTOR_FORCE_HARMONIC: usize = 1024;

/// Original force/torque amplitudes relative to the COMMON cycle time origin.
///
/// The physical load is `cosine*cos(omega*t) + sine*sin(omega*t)`, where
/// `omega = 2*pi*order/period`. For an `exp(+i omega t)` force phasor, use
/// its real part as `cosine` and minus its imaginary part as `sine`.
/// Entries have the model's conjugate force/torque units. Each segment may
/// declare different amplitudes; those explicit jumps introduce no impulse.
#[derive(Clone, Copy, Debug)]
pub struct RotorHarmonicForce<'a> {
    /// Original generalized force/torque amplitudes multiplying cos(omega*t).
    pub cosine: &'a [f64],
    /// Original generalized force/torque amplitudes multiplying sin(omega*t).
    pub sine: &'a [f64],
}

/// Prepare a full periodic rotor under a continuous harmonic physical drive.
///
/// `harmonic_order` is a positive integer bounded by MAX_ROTOR_FORCE_HARMONIC.
/// The common phase advances across every segment; it never restarts at a
/// bearing/stiffness switch. The source completes exactly that many revolutions
/// per declared period, so the retained affine cycle can be repeated. A signed
/// rotating unbalance is expressed by the sine amplitudes, independently of the
/// signed shaft speed entering the original gyroscopic operator. Synchronous
/// operation requires the caller's period/order to match its shaft speed.
///
/// Two dimensionless source coordinates evolve as a skew rotation, and the
/// existing native exponential integrates their coupling into all original
/// physical coordinates. No inverse dynamic stiffness or held approximation
/// of a sinusoid is used, including at resonance. The auxiliary coordinates
/// carry neither physical energy nor Floquet multipliers. Their forcing scale
/// is normalized before exponentiation; the largest exponential is 66 by 66.
/// The original full mass residuals are checked independently for BOTH loads.
///
/// Numerical exponential/residual error is distinct from model validity. Full
/// force-work receipts, smooth coefficient modulation, speed ramps and harmonic
/// superposition are outside this single-harmonic segmented-coefficient API.
/// Failure/cancellation publishes no partial plan. Warm advance/orbit queries
/// reuse the same maps and solve owner as held forcing.
#[allow(clippy::too_many_arguments)] // explicit physical cycle, source and caller budgets
pub fn prepare_harmonic_cycle(
    spin_rad_s: f64,
    segments: &[RotorPeriodicSegment<'_>],
    forces: &[RotorHarmonicForce<'_>],
    harmonic_order: usize,
    relative_mass_residual_tolerance: f64,
    budget: RotorFloquetBudget,
    cancelled: impl FnMut() -> bool,
) -> Result<PreparedRotorCycle, RotorError> {
    prepare_cycle(
        spin_rad_s,
        segments,
        CycleForcing::Harmonic {
            rows: forces,
            order: harmonic_order,
        },
        relative_mass_residual_tolerance,
        budget,
        cancelled,
    )
}

#[allow(clippy::too_many_arguments)] // shared original mass and physical phase
pub(super) fn segment_flow(
    mass: &MassSolve<'_>,
    generator: &[f64],
    force: RotorHarmonicForce<'_>,
    order: usize,
    start: f64,
    duration: f64,
    period: f64,
    tolerance: f64,
    budget: DenseExponentialBudget,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(AffineSegmentMap, f64), RotorError> {
    poll(cancelled)?;
    let d = 2 * mass.model.dofs();
    let (cosine, cosine_residual) =
        mass.force_column(force.cosine, duration, period, tolerance, cancelled)?;
    let (sine, sine_residual) =
        mass.force_column(force.sine, duration, period, tolerance, cancelled)?;
    let residual = cosine_residual.max(sine_residual);
    let scale = cosine
        .iter()
        .chain(&sine)
        .fold(0.0_f64, |s, x| s.max(x.abs()));
    if scale == 0.0 {
        return Ok((
            affine_segment_flow(generator, &cosine, d, budget, cancelled)?,
            residual,
        ));
    }
    let e = d + 2;
    let mut augmented = vec![0.0; e * e];
    for i in 0..d {
        poll(cancelled)?;
        augmented[i * e..i * e + d].copy_from_slice(&generator[i * d..(i + 1) * d]);
        for (j, column) in [&cosine, &sine].iter().enumerate() {
            let normalized = finite(column[i] / scale)?;
            if column[i] != 0.0 && normalized == 0.0 {
                return Err(RotorError::Invalid("harmonic force direction underflow"));
            }
            augmented[i * e + d + j] = normalized;
        }
    }
    // The shared generator already includes duration. In normalized local
    // time, [cos(theta),sin(theta)]' = omega*h * [-sin(theta),cos(theta)].
    let angle = finite((core::f64::consts::TAU * order as f64) * (duration / period))?;
    if angle == 0.0 {
        return Err(RotorError::Invalid("harmonic segment phase underflow"));
    }
    augmented[d * e + d + 1] = -angle;
    augmented[(d + 1) * e + d] = angle;
    poll(cancelled)?;
    let result =
        exponential_with_budget(&augmented, e, budget).map_err(RotorError::PeriodicExponential)?;
    poll(cancelled)?;
    // Reduce cycles before converting to radians, keeping argument reduction
    // independent of the number of prior revolutions and of segment count.
    let turns = (Dd::from_f64(start) * Dd::from_f64(order as f64) / Dd::from_f64(period)).to_f64();
    let phase = core::f64::consts::TAU * (turns - turns.floor());
    let (sin, cos) = det::sin_cos(phase);
    let mut transition = vec![0.0; d * d];
    let mut offset = Vec::with_capacity(d);
    for i in 0..d {
        poll(cancelled)?;
        transition[i * d..(i + 1) * d].copy_from_slice(&result[i * e..i * e + d]);
        let value = (Dd::from_f64(result[i * e + d]) * Dd::from_f64(cos)
            + Dd::from_f64(result[i * e + d + 1]) * Dd::from_f64(sin))
        .to_f64();
        offset.push(product_three(finite(value)?, scale, 1.0)?);
    }
    Ok((AffineSegmentMap { transition, offset }, residual))
}
