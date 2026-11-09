//! Periodic physical forcing with a steady component and several harmonics.
use super::{
    AffineSegmentMap, CycleForcing, DenseExponentialBudget, MassSolve, PreparedRotorCycle,
    RotorError, RotorFloquetBudget, RotorHarmonicForce, RotorPeriodicSegment, affine_segment_flow,
    finite, harmonic, poll, prepare_cycle,
};
use fs_math::dd::Dd;

/// Maximum simultaneous nonzero-frequency components in one declared cycle.
pub const MAX_ROTOR_FORCE_TONES: usize = 32;

/// Original generalized forces/torques during one coefficient segment.
///
/// At common cycle time t, the force is
/// `steady + sum_j (harmonics[j].cosine*cos(omega_j*t)
///                 + harmonics[j].sine*sin(omega_j*t))`,
/// where `omega_j = 2*pi*harmonic_orders[j]/period`. The arrays retain the
/// physical model's coordinate order and power-conjugate load units.
#[derive(Clone, Copy, Debug)]
pub struct RotorFourierForce<'a> {
    /// The physical DC force/torque, held through this segment.
    pub steady: &'a [f64],
    /// One cosine/sine pair for each common declared harmonic order.
    pub harmonics: &'a [RotorHarmonicForce<'a>],
}

/// Prepare the complete forced rotor cycle under a finite Fourier waveform.
///
/// `harmonic_orders` is strictly increasing, contains at most 32 entries,
/// and every order is in 1..=MAX_ROTOR_FORCE_HARMONIC. Each force row must
/// contain exactly those tones. Empty orders give the existing held-load
/// result, including a valid zero-force plan when every steady row is zero.
///
/// Every tone uses the same physical cycle origin; coefficient or amplitude
/// switches do not reset its phase or the original rotor motion. All tones
/// repeat after the complete period. This is a linear superposition of actual
/// generalized loads; it does not infer nonlinear bearing forces, gear-mesh
/// spectra, speed ramps, or a truncation error for a measured waveform.
///
/// Preparation reuses one mass factor and generator per segment, integrates
/// each tone through the existing normalized harmonic exponential, and sums
/// their particular responses in double-double arithmetic. The homogeneous
/// map is retained ONCE from the steady-source flow; physical eigenvalues are
/// never duplicated or contaminated by auxiliary forcing frequencies. At an
/// undamped resonance the transient secular response remains represented.
///
/// One exponential is at most 66 by 66 regardless of tone count. Cold work is
/// linear in the number of tones, O(segments*(1+tones)*(2n+2)^3); retained
/// maps and warm trajectory/orbit costs are independent of tone count. Caller
/// budgets reserve the complete cold work before assembly, including zero
/// tones that are subsequently skipped. Each nonzero force column receives
/// the existing independent original-mass residual check. Cancellation or
/// numerical refusal publishes no partial prepared cycle.
///
/// Residuals describe the declared finite model, not a complete physical
/// force-work audit or an error bound on an underlying continuous waveform.
/// The companion power::prepare_fourier_power_cycle integrates that physical
/// work with all cross-tone terms retained.
#[allow(clippy::too_many_arguments)] // declared physical cycle and explicit caller budgets
pub fn prepare_fourier_cycle(
    spin_rad_s: f64,
    segments: &[RotorPeriodicSegment<'_>],
    forces: &[RotorFourierForce<'_>],
    harmonic_orders: &[usize],
    relative_mass_residual_tolerance: f64,
    budget: RotorFloquetBudget,
    cancelled: impl FnMut() -> bool,
) -> Result<PreparedRotorCycle, RotorError> {
    prepare_cycle(
        spin_rad_s,
        segments,
        CycleForcing::Fourier {
            rows: forces,
            orders: harmonic_orders,
        },
        relative_mass_residual_tolerance,
        budget,
        cancelled,
    )
}

#[allow(clippy::too_many_arguments)] // common original mass, phase and bounded exponential
pub(super) fn segment_flow(
    mass: &MassSolve<'_>,
    generator: &[f64],
    force: RotorFourierForce<'_>,
    orders: &[usize],
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
        orders,
        start,
        duration,
        period,
        tolerance,
        budget,
        cancelled,
    )
}

// One compensated source sum and homogeneous transition serve mechanical and
// hereditary cycles alike; only the admitted physical state count differs.
#[allow(clippy::too_many_arguments)]
pub(super) fn segment_flow_for_states(
    mass: &MassSolve<'_>,
    generator: &[f64],
    d: usize,
    force: RotorFourierForce<'_>,
    orders: &[usize],
    start: f64,
    duration: f64,
    period: f64,
    tolerance: f64,
    budget: DenseExponentialBudget,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(AffineSegmentMap, f64), RotorError> {
    segment_flow_for_states_with_clocks(
        mass,
        generator,
        d,
        force,
        orders.iter().copied().map(harmonic::SourceClock::Periodic),
        Dd::from_f64(start),
        duration,
        period,
        tolerance,
        budget,
        cancelled,
    )
}

#[allow(clippy::too_many_arguments)]
pub(super) fn segment_flow_for_states_with_clocks(
    mass: &MassSolve<'_>,
    generator: &[f64],
    d: usize,
    force: RotorFourierForce<'_>,
    clocks: impl Iterator<Item = harmonic::SourceClock>,
    start: Dd,
    duration: f64,
    period: f64,
    tolerance: f64,
    budget: DenseExponentialBudget,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(AffineSegmentMap, f64), RotorError> {
    poll(cancelled)?;
    if d < 2 * mass.model.dofs() || d.checked_mul(d) != Some(generator.len()) {
        return Err(RotorError::Invalid(
            "complete periodic state generator and original mechanical prefix required",
        ));
    }
    let (mut steady, mut residual) =
        mass.force_column(force.steady, duration, period, tolerance, cancelled)?;
    steady.resize(d, 0.0);
    let mut flow = affine_segment_flow(generator, &steady, d, budget, cancelled)?;
    let mut total: Vec<_> = flow
        .offset
        .iter()
        .map(|&value| Dd::from_f64(value))
        .collect();
    for (clock, &tone) in clocks.zip(force.harmonics) {
        poll(cancelled)?;
        // The admitted zero load has exactly zero particular response. Avoid
        // repeating a homogeneous exponential solely for a silent component.
        if tone
            .cosine
            .iter()
            .chain(tone.sine)
            .all(|&value| value == 0.0)
        {
            continue;
        }
        let (contribution, checked) = harmonic::segment_flow_for_states_with_clock(
            mass, generator, d, tone, clock, start, duration, period, tolerance, budget, cancelled,
        )?;
        residual = residual.max(checked);
        for (sum, value) in total.iter_mut().zip(contribution.offset) {
            *sum = *sum + Dd::from_f64(value);
        }
    }
    for (value, sum) in flow.offset.iter_mut().zip(total) {
        *value = finite(sum.to_f64())?;
    }
    poll(cancelled)?;
    Ok((flow, residual))
}
