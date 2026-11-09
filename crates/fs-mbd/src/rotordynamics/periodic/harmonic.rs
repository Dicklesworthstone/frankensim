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
/// Numerical exponential/residual error is distinct from model validity.
/// Smooth coefficient modulation and speed ramps are outside this segmented
/// API. The companion power::prepare_harmonic_power_cycle adds independently
/// integrated physical work. Use prepare_fourier_cycle for simultaneous steady
/// loading and multiple harmonics of the same period.
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

/// One normalized native source and its phase at the common cycle origin.
/// Motion and quadratic work use these same original-mass solves and rotation.
pub(super) struct HarmonicSource {
    pub(super) columns: [Vec<f64>; 2],
    pub(super) scale: f64,
    pub(super) angle: f64,
    pub(super) phase: [f64; 2],
    pub(super) maximum_mass_residual: f64,
}

#[derive(Clone, Copy)]
pub(super) enum SourceClock {
    Periodic(usize),
    Physical(f64),
}

impl SourceClock {
    pub(super) fn phase(
        self,
        start: Dd,
        duration: f64,
        period: f64,
    ) -> Result<(f64, [f64; 2]), RotorError> {
        let (angle, phase) = match self {
            Self::Periodic(order) => {
                let angle = finite((core::f64::consts::TAU * order as f64) * (duration / period))?;
                let turns = (start * Dd::from_f64(order as f64) / Dd::from_f64(period)).to_f64();
                (angle, core::f64::consts::TAU * (turns - turns.floor()))
            }
            Self::Physical(omega) => {
                // Treat omega and t as the original binary inputs. Split 2*pi
                // retains fractional turns before reducing the DD product;
                // converting omega*t to f64 first would silently lose phase.
                let product = Dd::from_f64(omega) * start;
                let angle = finite((Dd::from_f64(omega) * Dd::from_f64(duration)).to_f64())?;
                if !omega.is_finite()
                    || omega <= 0.0
                    || !start.to_f64().is_finite()
                    || !product.to_f64().is_finite()
                    || (start != Dd::ZERO && product.to_f64() == 0.0)
                    || product.to_f64().abs() > 1_099_511_627_776.0
                    || angle > core::f64::consts::TAU * MAX_ROTOR_FORCE_HARMONIC as f64
                {
                    return Err(RotorError::Invalid(
                        "bounded representable physical source frequency and global phase required",
                    ));
                }
                let tau = Dd::from_f64(core::f64::consts::TAU)
                    + Dd::from_f64(2.449_293_598_294_706_4e-16);
                let turns = (product / tau).to_f64().floor();
                (
                    angle,
                    finite((product - Dd::from_f64(turns) * tau).to_f64())?,
                )
            }
        };
        if angle == 0.0 {
            return Err(RotorError::Invalid("harmonic segment phase underflow"));
        }
        let (sin, cos) = det::sin_cos(phase);
        Ok((angle, [cos, sin]))
    }
}

#[allow(clippy::too_many_arguments)] // original physical source, time and mass admission
pub(super) fn prepare_source(
    mass: &MassSolve<'_>,
    force: RotorHarmonicForce<'_>,
    order: usize,
    start: f64,
    duration: f64,
    period: f64,
    tolerance: f64,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<HarmonicSource, RotorError> {
    prepare_source_with_clock(
        mass,
        force,
        SourceClock::Periodic(order),
        Dd::from_f64(start),
        duration,
        period,
        tolerance,
        cancelled,
    )
}

#[allow(clippy::too_many_arguments)]
pub(super) fn prepare_source_with_clock(
    mass: &MassSolve<'_>,
    force: RotorHarmonicForce<'_>,
    clock: SourceClock,
    start: Dd,
    duration: f64,
    period: f64,
    tolerance: f64,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<HarmonicSource, RotorError> {
    poll(cancelled)?;
    let (cosine, cosine_residual) =
        mass.force_column(force.cosine, duration, period, tolerance, cancelled)?;
    let (sine, sine_residual) =
        mass.force_column(force.sine, duration, period, tolerance, cancelled)?;
    let maximum_mass_residual = cosine_residual.max(sine_residual);
    let scale = cosine
        .iter()
        .chain(&sine)
        .fold(0.0_f64, |s, x| s.max(x.abs()));
    let mut columns = [cosine, sine];
    if scale == 0.0 {
        return Ok(HarmonicSource {
            columns,
            scale,
            angle: 0.0,
            phase: [1.0, 0.0],
            maximum_mass_residual,
        });
    }
    for column in &mut columns {
        for value in column {
            let normalized = finite(*value / scale)?;
            if *value != 0.0 && normalized == 0.0 {
                return Err(RotorError::Invalid("harmonic force direction underflow"));
            }
            *value = normalized;
        }
    }
    let (angle, phase) = clock.phase(start, duration, period)?;
    Ok(HarmonicSource {
        columns,
        scale,
        angle,
        phase,
        maximum_mass_residual,
    })
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
    segment_flow_for_states(
        mass,
        generator,
        2 * mass.model.dofs(),
        force,
        order,
        start,
        duration,
        period,
        tolerance,
        budget,
        cancelled,
    )
}

// The same original force/phase owner advances every admitted physical state.
// Hereditary callers retain their full generator and append zero DIRECT source
// entries for memory; the generator still transfers motion into every stress.
#[allow(clippy::too_many_arguments)]
pub(super) fn segment_flow_for_states(
    mass: &MassSolve<'_>,
    generator: &[f64],
    d: usize,
    force: RotorHarmonicForce<'_>,
    order: usize,
    start: f64,
    duration: f64,
    period: f64,
    tolerance: f64,
    budget: DenseExponentialBudget,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(AffineSegmentMap, f64), RotorError> {
    segment_flow_for_states_with_clock(
        mass,
        generator,
        d,
        force,
        SourceClock::Periodic(order),
        Dd::from_f64(start),
        duration,
        period,
        tolerance,
        budget,
        cancelled,
    )
}

#[allow(clippy::too_many_arguments)]
pub(super) fn segment_flow_for_states_with_clock(
    mass: &MassSolve<'_>,
    generator: &[f64],
    d: usize,
    force: RotorHarmonicForce<'_>,
    clock: SourceClock,
    start: Dd,
    duration: f64,
    period: f64,
    tolerance: f64,
    budget: DenseExponentialBudget,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(AffineSegmentMap, f64), RotorError> {
    if d < 2 * mass.model.dofs() || d.checked_mul(d) != Some(generator.len()) {
        return Err(RotorError::Invalid(
            "complete periodic state generator and original mechanical prefix required",
        ));
    }
    let mut source = prepare_source_with_clock(
        mass, force, clock, start, duration, period, tolerance, cancelled,
    )?;
    for column in &mut source.columns {
        column.resize(d, 0.0);
    }
    if source.scale == 0.0 {
        return Ok((
            affine_segment_flow(generator, &source.columns[0], d, budget, cancelled)?,
            source.maximum_mass_residual,
        ));
    }
    let e = d + 2;
    let mut augmented = vec![0.0; e * e];
    for i in 0..d {
        poll(cancelled)?;
        augmented[i * e..i * e + d].copy_from_slice(&generator[i * d..(i + 1) * d]);
        for j in 0..2 {
            augmented[i * e + d + j] = source.columns[j][i];
        }
    }
    augmented[d * e + d + 1] = -source.angle;
    augmented[(d + 1) * e + d] = source.angle;
    poll(cancelled)?;
    let result =
        exponential_with_budget(&augmented, e, budget).map_err(RotorError::PeriodicExponential)?;
    poll(cancelled)?;
    let mut transition = vec![0.0; d * d];
    let mut offset = Vec::with_capacity(d);
    for i in 0..d {
        poll(cancelled)?;
        transition[i * d..(i + 1) * d].copy_from_slice(&result[i * e..i * e + d]);
        let value = (Dd::from_f64(result[i * e + d]) * Dd::from_f64(source.phase[0])
            + Dd::from_f64(result[i * e + d + 1]) * Dd::from_f64(source.phase[1]))
        .to_f64();
        offset.push(product_three(finite(value)?, source.scale, 1.0)?);
    }
    Ok((
        AffineSegmentMap { transition, offset },
        source.maximum_mass_residual,
    ))
}
