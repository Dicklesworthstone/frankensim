//! Period maps of declared piecewise-constant linear rotor-bearing cycles.
//!
//! Every segment retains the original equation
//! M q'' + (C + Omega G) q' + K q = 0. Signed spin, full mass and skew G
//! stay fixed; full C/K may change, including cross-coupled or active entries.
//! Coordinates and velocities are continuous at each coefficient switch.
//! Changing stiffness can exchange energy with its modulator; no change is
//! silently assigned heat.
//!
//! Native bounded matrix exponentials advance each constant segment. Their
//! chronological product E_last ... E_first is the complete monodromy map.
//! Its numerical multipliers can expose parametric resonance even when every
//! frozen operating point has stable poles. No modal truncation, artificial
//! damping or eigenvalue clipping is applied.
//!
//! This is the declared segmented model. Approximating smooth modulation by
//! such segments additionally requires coefficient/time refinement. Neither
//! eigenvalue error bounds nor a stability certificate are inferred. In
//! particular, rounding around unit spectral radius is not classified.

use super::{
    MAX_ROTOR_DOFS, ROTOR_QR_SWEEPS_PER_WINDOW, RotorBearingModel, RotorError, finite, poll,
};
use fs_la::eigen_complex::{EigError, eig_with_cancel};
use fs_la::factor::{Cholesky, cholesky};
use fs_math::{c64::C64, dd::Dd, det};
use fs_phs::linear_exact::{DenseExponentialBudget, exponential_with_budget};

mod forcing;
pub use forcing::{RotorPeriodicOrbit, RotorPeriodicTrajectory, prepare_forced_cycle};
mod harmonic;
pub use harmonic::{MAX_ROTOR_FORCE_HARMONIC, RotorHarmonicForce, prepare_harmonic_cycle};
mod fourier;
pub use fourier::{MAX_ROTOR_FORCE_TONES, RotorFourierForce, prepare_fourier_cycle};
/// Exact signed work and stiffness-modulation power on original forced cycles.
pub mod power;
mod relaxation;
pub use relaxation::power::{
    PreparedRotorMaxwellPowerCycle, PreparedRotorMaxwellPowerInterval, RotorMaxwellCyclePower,
    RotorMaxwellIntervalPower, RotorMaxwellSegmentWork, prepare_maxwell_forced_power_cycle,
    prepare_maxwell_fourier_power_cycle, prepare_maxwell_harmonic_power_cycle,
    prepare_maxwell_power_interval,
};
pub use relaxation::{
    PreparedRotorMaxwellCycle, PreparedRotorMaxwellInterval, RotorIntervalForcing,
    RotorMaxwellFloquetResponse, RotorMaxwellIntervalResponse, RotorMaxwellIntervalTrajectory,
    RotorMaxwellKinetics, RotorMaxwellPeriodicOrbit, RotorMaxwellPeriodicTrajectory,
    RotorMaxwellState, periodic_maxwell_floquet, prepare_maxwell_cycle,
    prepare_maxwell_forced_cycle, prepare_maxwell_fourier_cycle, prepare_maxwell_harmonic_cycle,
    prepare_maxwell_interval,
};

/// Prescribed physical loading on the original coordinates of each cycle segment.
///
/// Rows follow the supplied segment order. Harmonic phase is continuous in the
/// full cycle time: order `h` means angular frequency `2*pi*h/period`, with time
/// zero at the first segment boundary. Coefficients may switch at boundaries;
/// source phase never restarts there. Preparation checks all dimensions, finite
/// values and harmonic bounds through the same admission as the named forcing
/// entry points.
#[derive(Clone, Copy, Debug)]
pub enum RotorPeriodicForcing<'a> {
    /// Zero applied generalized force throughout the cycle.
    Homogeneous,
    /// One constant generalized-force vector per segment.
    Held(&'a [&'a [f64]]),
    /// One cosine/sine pair per segment at a common positive harmonic order.
    Harmonic {
        /// Physical force coefficients in original coordinate order.
        rows: &'a [RotorHarmonicForce<'a>],
        /// Positive integer multiple of the full-cycle fundamental frequency.
        order: usize,
    },
    /// A steady force and a bounded set of commensurate tones per segment.
    Fourier {
        /// Physical force coefficients, with harmonics in `orders` order.
        rows: &'a [RotorFourierForce<'a>],
        /// Strictly increasing positive integer harmonic orders.
        orders: &'a [usize],
    },
}

use RotorPeriodicForcing as CycleForcing;

impl CycleForcing<'_> {
    fn source_coordinates(self) -> usize {
        match self {
            Self::Homogeneous => 0,
            Self::Held(_) => 1,
            Self::Harmonic { .. } => 2,
            Self::Fourier { orders, .. } => {
                if orders.is_empty() {
                    1
                } else {
                    2
                }
            }
        }
    }

    fn exponentials_per_segment(self) -> usize {
        match self {
            Self::Fourier { orders, .. } => 1 + orders.len(),
            _ => 1,
        }
    }

    fn admit(self, n: usize, segments: usize) -> Result<(), RotorError> {
        let valid_row = |row: &[f64]| row.len() == n && row.iter().all(|x| x.is_finite());
        let valid = match self {
            Self::Homogeneous => true,
            Self::Held(rows) => rows.len() == segments && rows.iter().all(|row| valid_row(row)),
            Self::Harmonic { rows, order } => {
                (1..=MAX_ROTOR_FORCE_HARMONIC).contains(&order)
                    && rows.len() == segments
                    && rows
                        .iter()
                        .all(|row| valid_row(row.cosine) && valid_row(row.sine))
            }
            Self::Fourier { rows, orders } => {
                orders.len() <= MAX_ROTOR_FORCE_TONES
                    && orders
                        .iter()
                        .all(|order| (1..=MAX_ROTOR_FORCE_HARMONIC).contains(order))
                    && orders.windows(2).all(|pair| pair[0] < pair[1])
                    && rows.len() == segments
                    && rows.iter().all(|row| {
                        valid_row(row.steady)
                            && row.harmonics.len() == orders.len()
                            && row
                                .harmonics
                                .iter()
                                .all(|tone| valid_row(tone.cosine) && valid_row(tone.sine))
                    })
            }
        };
        if valid {
            Ok(())
        } else {
            Err(RotorError::Invalid(
                "matching finite segment loads and bounded strictly increasing harmonic orders required",
            ))
        }
    }
}

/// Hard bound on coefficient switches in a single declared physical cycle.
pub const MAX_ROTOR_PERIODIC_SEGMENTS: usize = 256;

/// One positive-duration constant-coefficient portion of a periodic rotor cycle.
///
/// Models must share the same coordinate order, units, frame, full mass and G.
/// Their C/K are actual coefficients for this portion of the cycle. The final
/// segment is followed by the first, without an impulse or state reset.
#[derive(Clone, Copy, Debug)]
pub struct RotorPeriodicSegment<'a> {
    /// Original full physical matrices; borrowed unchanged.
    pub model: &'a RotorBearingModel,
    /// Physical duration [s], finite and strictly positive.
    pub duration_seconds: f64,
}

/// Explicit admission for the complete period map and its numerical spectrum.
///
/// The existing 32-coordinate rotor bound still applies. Dense bytes include
/// all owner workspace and returned arrays, excluding the borrowed models.
/// Maxwell cycles additionally count every original material state, up to
/// 64 branches; larger valid models may require a larger explicit work budget.
/// Work conservatively reserves 80 dense products per segment exponential,
/// and additionally mass solves, chronological products and bounded complex QR.
/// Fourier forcing reserves one held-source exponential plus one per tone.
/// Open-interval preparation uses the same limits without QR work or sweep
/// admission; it does not compute a Floquet spectrum.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct RotorFloquetBudget {
    /// Supplied segment count ceiling, in 1..=MAX_ROTOR_PERIODIC_SEGMENTS.
    pub max_segments: usize,
    /// Native QR sweeps per active window, in 1..=ROTOR_QR_SWEEPS_PER_WINDOW.
    /// Ignored by open-interval preparation, which performs no spectral solve.
    pub max_qr_sweeps_per_window: usize,
    /// Aggregate live dense payload ceiling in bytes.
    pub max_dense_bytes: usize,
    /// Conservative scalar multiply/add term ceiling for the whole query.
    pub max_multiply_terms: u64,
}

impl Default for RotorFloquetBudget {
    fn default() -> Self {
        Self {
            max_segments: 64,
            max_qr_sweeps_per_window: ROTOR_QR_SWEEPS_PER_WINDOW,
            max_dense_bytes: 4 * 1024 * 1024,
            max_multiply_terms: 4_000_000_000,
        }
    }
}

impl RotorFloquetBudget {
    fn admit(
        self,
        n: usize,
        segments: usize,
        source_coordinates: usize,
        exponentials_per_segment: usize,
    ) -> Result<DenseExponentialBudget, RotorError> {
        self.admit_states(
            2 * n,
            segments,
            source_coordinates,
            exponentials_per_segment,
            source_coordinates != 0,
        )
    }

    fn admit_states(
        self,
        d: usize,
        segments: usize,
        source_coordinates: usize,
        exponentials_per_segment: usize,
        retain_segments: bool,
    ) -> Result<DenseExponentialBudget, RotorError> {
        self.admit_states_with_spectrum(
            d,
            segments,
            source_coordinates,
            exponentials_per_segment,
            retain_segments,
            true,
        )
    }

    fn admit_interval_states(
        self,
        d: usize,
        segments: usize,
        source_coordinates: usize,
        exponentials_per_segment: usize,
        retain_segments: bool,
    ) -> Result<DenseExponentialBudget, RotorError> {
        self.admit_states_with_spectrum(
            d,
            segments,
            source_coordinates,
            exponentials_per_segment,
            retain_segments,
            false,
        )
    }

    #[allow(clippy::too_many_arguments)] // one shared physical-work admission
    fn admit_states_with_spectrum(
        self,
        d: usize,
        segments: usize,
        source_coordinates: usize,
        exponentials_per_segment: usize,
        retain_segments: bool,
        spectrum: bool,
    ) -> Result<DenseExponentialBudget, RotorError> {
        if self.max_segments == 0
            || self.max_segments > MAX_ROTOR_PERIODIC_SEGMENTS
            || segments > self.max_segments
            || (spectrum
                && (self.max_qr_sweeps_per_window == 0
                    || self.max_qr_sweeps_per_window > ROTOR_QR_SWEEPS_PER_WINDOW))
        {
            return Err(RotorError::PeriodicBudget("segment or QR sweep admission"));
        }
        // Original states <= 128 and segments <= 256 were checked before
        // these products. Material states consume the same full-state budget.
        // Held/harmonic sources add one/two normalized coordinates to the shared
        // exponential, never physical modes. Forced plans retain each
        // segment map and reserve one complete boundary trajectory and solve.
        let e = d + source_coordinates;
        let retained = if retain_segments {
            segments * (d * d + 3 * d + 16) + d * d + 16 * d
        } else {
            0
        };
        // A Fourier sum retains the held-source map and compensated offset
        // while evaluating one tone at a time. No harmonic state bank is kept.
        let superposition = if exponentials_per_segment > 1 {
            d * d + 4 * d
        } else {
            0
        };
        let bytes = (16 * e * e + 32 * e + retained + superposition) * core::mem::size_of::<f64>();
        let cube = (d as u64).pow(3);
        let exponential_work =
            (e as u64).pow(3) * 84 * segments as u64 * exponentials_per_segment as u64;
        let spectrum_work = if spectrum {
            64 * self.max_qr_sweeps_per_window as u64
        } else {
            0
        };
        // Keep the original fixed mass/assembly allowance in either path;
        // interval propagation omits only the unused spectral iteration work.
        let work =
            exponential_work + cube * (64 + spectrum_work + if retain_segments { 32 } else { 0 });
        if bytes > self.max_dense_bytes || work > self.max_multiply_terms {
            return Err(RotorError::PeriodicBudget(
                "dense byte or whole-cycle work admission",
            ));
        }
        Ok(DenseExponentialBudget {
            max_dimension: e,
            max_dense_bytes: self.max_dense_bytes,
            max_multiply_terms: self.max_multiply_terms,
        })
    }
}

/// Complete numerical Floquet analysis of the original segmented model.
#[derive(Clone, Debug, PartialEq)]
pub struct RotorFloquetResponse {
    /// Signed constant shaft speed used in every C + Omega G [rad/s].
    pub spin_rad_s: f64,
    /// Sum of the declared physical segment durations [s].
    pub period_seconds: f64,
    /// Number of original generalized displacement coordinates.
    pub dofs: usize,
    /// Original-coordinate row-major 2n-by-2n map, in [q, q_dot] state order.
    /// Acting on an initial state gives its homogeneous one-period endpoint.
    pub monodromy: Vec<f64>,
    /// All 2n multipliers, in native canonical (real, imaginary) order.
    /// These are numerical estimates without eigenvectors or error enclosures.
    pub multipliers: Vec<C64>,
    /// Maximum computed multiplier magnitude. Unit-radius rounding is retained.
    pub spectral_radius: f64,
    /// ln(spectral_radius) / period [1/s], without a stability classification.
    pub growth_rate_per_s: f64,
    /// Largest independently reconstructed, mass-scaled original generator
    /// or forcing-column residual. Not an eigenvalue or motion bound.
    pub maximum_mass_residual: f64,
}

/// Cached full affine rotor cycle with a fixed physical force schedule.
///
/// Construct with prepare_forced_cycle, prepare_harmonic_cycle or
/// prepare_fourier_cycle. Held loads and continuous harmonic sources share the
/// same complete coupled operators.
/// Coefficient and force changes introduce no state reset.
/// Warm trajectory calls use the retained maps without new exponentials.
/// The homogeneous Floquet report does not classify this forced orbit as
/// attracting: unstable perturbations remain unstable.
#[derive(Clone, Debug)]
pub struct PreparedRotorCycle {
    report: RotorFloquetResponse,
    state_scale: Vec<f64>,
    monodromy: Vec<f64>,
    offset: Vec<f64>,
    physical_offset: Vec<f64>,
    segment_maps: Vec<AffineSegmentMap>,
    boundary_times: Vec<f64>,
}

#[derive(Clone, Debug)]
struct AffineSegmentMap {
    transition: Vec<f64>,
    offset: Vec<f64>,
}

/// Analyze a homogeneous cycle of original constant-coefficient rotor segments.
///
/// M and G must agree coefficient by coefficient across all segments. Spin
/// remains constant, so no spin-acceleration, changing-inertia or impulse law
/// is inferred. All C/K entries survive, including nonconservative ones.
///
/// An invertible mass/time coordinate scaling conditions the full state.
/// Native Cholesky solves form the generator; original separate C and Omega G
/// enter the independent column-residual checks. The mass residual tolerance
/// must be finite and in (0,1). No inverse, mode reduction or damping floor is
/// substituted for the full physical matrices.
///
/// Cancellation is polled during assembly, products and native QR, and before
/// and after each shared exponential. One exponential is an indivisible call
/// bounded to 64 states, 32 Taylor products and 48 squarings; this is not an
/// operating-system latency guarantee. Any refusal publishes no partial map.
pub fn periodic_floquet(
    spin_rad_s: f64,
    segments: &[RotorPeriodicSegment<'_>],
    relative_mass_residual_tolerance: f64,
    budget: RotorFloquetBudget,
    cancelled: impl FnMut() -> bool,
) -> Result<RotorFloquetResponse, RotorError> {
    Ok(prepare_cycle(
        spin_rad_s,
        segments,
        CycleForcing::Homogeneous,
        relative_mass_residual_tolerance,
        budget,
        cancelled,
    )?
    .report)
}

fn prepare_cycle(
    spin_rad_s: f64,
    segments: &[RotorPeriodicSegment<'_>],
    forces: CycleForcing<'_>,
    relative_mass_residual_tolerance: f64,
    budget: RotorFloquetBudget,
    mut cancelled: impl FnMut() -> bool,
) -> Result<PreparedRotorCycle, RotorError> {
    poll(&mut cancelled)?;
    if segments.is_empty()
        || segments.len() > MAX_ROTOR_PERIODIC_SEGMENTS
        || !spin_rad_s.is_finite()
        || !relative_mass_residual_tolerance.is_finite()
        || relative_mass_residual_tolerance <= 0.0
        || relative_mass_residual_tolerance >= 1.0
    {
        return Err(RotorError::Invalid(
            "finite spin, nonempty bounded cycle and mass tolerance in (0,1) required",
        ));
    }
    let reference = segments[0].model;
    let n = reference.dofs();
    forces.admit(n, segments.len())?;
    let forced = forces.source_coordinates() != 0;
    let dense_budget = budget.admit(
        n,
        segments.len(),
        forces.source_coordinates(),
        forces.exponentials_per_segment(),
    )?;
    let mut boundary_times = if forced {
        Vec::with_capacity(segments.len() + 1)
    } else {
        Vec::new()
    };
    if forced {
        boundary_times.push(0.0);
    }
    let mut period = Dd::ZERO;
    for segment in segments {
        poll(&mut cancelled)?;
        if segment.model.dofs() != n
            || segment.model.mass() != reference.mass()
            || segment.model.gyro_per_rad_s() != reference.gyro_per_rad_s()
            || !segment.duration_seconds.is_finite()
            || segment.duration_seconds <= 0.0
        {
            return Err(RotorError::Invalid(
                "periodic segments require identical M/G and positive finite durations",
            ));
        }
        period = period + Dd::from_f64(segment.duration_seconds);
        let boundary = finite(period.to_f64())?;
        if forced {
            boundary_times.push(boundary);
        }
    }
    let period_seconds = finite(period.to_f64())?;
    let mass = MassSolve::new(reference, &mut cancelled)?;
    // z = [S q, T S q_dot], S_ii=sqrt(M_ii), T=physical period.
    // This retains free coordinates and makes the time scale common to all
    // segments; no per-segment basis changes contaminate the matrix product.
    let d = 2 * n;
    let mut state_scale = mass.scale.clone();
    for &scale in &mass.scale {
        let velocity_scale = finite(period_seconds * scale)?;
        if velocity_scale == 0.0 {
            return Err(RotorError::Invalid("periodic coordinate scale underflow"));
        }
        state_scale.push(velocity_scale);
    }
    let mut monodromy = vec![0.0; d * d];
    for i in 0..d {
        monodromy[i * d + i] = 1.0;
    }
    let mut maximum_mass_residual = 0.0_f64;
    let mut offset = vec![0.0; d];
    let mut segment_maps = if forced {
        Vec::with_capacity(segments.len())
    } else {
        Vec::new()
    };
    for (index, segment) in segments.iter().enumerate() {
        poll(&mut cancelled)?;
        let (generator, residual) = mass.segment_generator(
            segment,
            spin_rad_s,
            period_seconds,
            relative_mass_residual_tolerance,
            &mut cancelled,
        )?;
        maximum_mass_residual = maximum_mass_residual.max(residual);
        let flow = match forces {
            CycleForcing::Homogeneous => {
                affine_segment_flow(&generator, &vec![0.0; d], d, dense_budget, &mut cancelled)?
            }
            CycleForcing::Held(rows) => {
                let (column, residual) = mass.force_column(
                    rows[index],
                    segment.duration_seconds,
                    period_seconds,
                    relative_mass_residual_tolerance,
                    &mut cancelled,
                )?;
                maximum_mass_residual = maximum_mass_residual.max(residual);
                affine_segment_flow(&generator, &column, d, dense_budget, &mut cancelled)?
            }
            CycleForcing::Harmonic { rows, order } => {
                let (flow, residual) = harmonic::segment_flow(
                    &mass,
                    &generator,
                    rows[index],
                    order,
                    boundary_times[index],
                    segment.duration_seconds,
                    period_seconds,
                    relative_mass_residual_tolerance,
                    dense_budget,
                    &mut cancelled,
                )?;
                maximum_mass_residual = maximum_mass_residual.max(residual);
                flow
            }
            CycleForcing::Fourier { rows, orders } => {
                let (flow, residual) = fourier::segment_flow(
                    &mass,
                    &generator,
                    rows[index],
                    orders,
                    boundary_times[index],
                    segment.duration_seconds,
                    period_seconds,
                    relative_mass_residual_tolerance,
                    dense_budget,
                    &mut cancelled,
                )?;
                maximum_mass_residual = maximum_mass_residual.max(residual);
                flow
            }
        };
        // Left multiplication is physical chronological order. The source
        // offset has the same order: b_next = E_segment b + b_segment.
        if forced {
            offset =
                forcing::apply_affine(&flow.transition, &flow.offset, &offset, d, &mut cancelled)?;
        }
        monodromy = multiply(&flow.transition, &monodromy, d, &mut cancelled)?;
        if forced {
            segment_maps.push(flow);
        }
    }
    let complex: Vec<_> = monodromy
        .iter()
        .map(|&value| C64::new(value, 0.0))
        .collect();
    let multipliers = eig_with_cancel(&complex, d, budget.max_qr_sweeps_per_window, &mut cancelled)
        .map_err(|error| match error {
            EigError::Cancelled => RotorError::Cancelled,
            other => RotorError::PeriodicSpectrum(other),
        })?;
    let mut spectral_radius = 0.0_f64;
    for value in &multipliers {
        finite(value.re)?;
        finite(value.im)?;
        let magnitude = finite(det::hypot(value.re, value.im))?;
        if magnitude == 0.0 {
            return Err(RotorError::Invalid(
                "periodic multiplier is unresolved or underflowed",
            ));
        }
        spectral_radius = spectral_radius.max(magnitude);
    }
    // A finite-time fundamental map is invertible. No vanished direction
    // can justify a zero multiplier or an infinite decay claim.
    if spectral_radius == 0.0 {
        return Err(RotorError::Invalid("periodic multiplier radius underflow"));
    }
    let growth_rate_per_s = finite(det::ln(spectral_radius) / period_seconds)?;
    let balanced_monodromy = monodromy.clone();
    for i in 0..d {
        poll(&mut cancelled)?;
        for j in 0..d {
            monodromy[i * d + j] =
                restore_entry(monodromy[i * d + j], state_scale[i], state_scale[j])?;
        }
    }
    poll(&mut cancelled)?;
    let report = RotorFloquetResponse {
        spin_rad_s,
        period_seconds,
        dofs: n,
        monodromy,
        multipliers,
        spectral_radius,
        growth_rate_per_s,
        maximum_mass_residual,
    };
    let mut physical_offset = Vec::with_capacity(d);
    for i in 0..d {
        poll(&mut cancelled)?;
        physical_offset.push(restore_entry(offset[i], state_scale[i], 1.0)?);
    }
    poll(&mut cancelled)?;
    Ok(PreparedRotorCycle {
        report,
        state_scale,
        monodromy: balanced_monodromy,
        offset,
        physical_offset,
        segment_maps,
        boundary_times,
    })
}

struct MassSolve<'a> {
    model: &'a RotorBearingModel,
    scale: Vec<f64>,
    factor: Cholesky,
}

impl<'a> MassSolve<'a> {
    fn new(
        model: &'a RotorBearingModel,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, RotorError> {
        let n = model.dofs();
        let scale: Vec<_> = (0..n).map(|i| det::sqrt(model.mass()[i * n + i])).collect();
        let mut scaled = vec![0.0; n * n];
        for i in 0..n {
            poll(cancelled)?;
            for j in 0..n {
                scaled[i * n + j] = finite((model.mass()[i * n + j] / scale[i]) / scale[j])?;
            }
        }
        let factor = cholesky(&scaled, n).map_err(RotorError::TransientLinear)?;
        poll(cancelled)?;
        Ok(Self {
            model,
            scale,
            factor,
        })
    }

    fn segment_generator(
        &self,
        segment: &RotorPeriodicSegment<'_>,
        spin: f64,
        period: f64,
        tolerance: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(Vec<f64>, f64), RotorError> {
        let n = self.model.dofs();
        let d = 2 * n;
        let h = segment.duration_seconds;
        let fraction = finite(h / period)?;
        if fraction == 0.0 {
            return Err(RotorError::Invalid("periodic segment duration underflow"));
        }
        let mut generator = vec![0.0; d * d];
        let mut maximum = 0.0_f64;
        for j in 0..n {
            poll(cancelled)?;
            generator[j * d + n + j] = fraction;
            for stiffness in [true, false] {
                let (column, residual) =
                    self.solve_column(segment.model, spin, stiffness, j, tolerance, cancelled)?;
                maximum = maximum.max(residual);
                for i in 0..n {
                    if stiffness {
                        generator[(n + i) * d + j] = product_three(-column[i], h, period)?;
                    } else {
                        generator[(n + i) * d + n + j] = product_three(-column[i], h, 1.0)?;
                    }
                }
            }
        }
        Ok((generator, maximum))
    }

    #[allow(clippy::too_many_arguments)] // explicit source coefficient and original column audit
    fn solve_column(
        &self,
        model: &RotorBearingModel,
        spin: f64,
        stiffness: bool,
        column: usize,
        tolerance: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<([f64; MAX_ROTOR_DOFS], f64), RotorError> {
        let n = self.model.dofs();
        let mut rounded = [0.0; MAX_ROTOR_DOFS];
        let mut original = [Dd::ZERO; MAX_ROTOR_DOFS];
        for i in 0..n {
            let index = i * n + column;
            if stiffness {
                rounded[i] = model.stiffness()[index];
                original[i] = Dd::from_f64(rounded[i]);
            } else {
                rounded[i] =
                    finite(spin.mul_add(model.gyro_per_rad_s()[index], model.damping()[index]))?;
                original[i] = Dd::from_f64(model.damping()[index])
                    + Dd::from_f64(spin) * Dd::from_f64(model.gyro_per_rad_s()[index]);
            }
        }
        self.solve_rhs(
            &rounded[..n],
            &original[..n],
            self.scale[column],
            tolerance,
            cancelled,
        )
    }

    fn force_column(
        &self,
        force: &[f64],
        duration: f64,
        period: f64,
        tolerance: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(Vec<f64>, f64), RotorError> {
        let n = self.model.dofs();
        let mut original = [Dd::ZERO; MAX_ROTOR_DOFS];
        for i in 0..n {
            original[i] = Dd::from_f64(force[i]);
        }
        let (acceleration, residual) =
            self.solve_rhs(force, &original[..n], 1.0, tolerance, cancelled)?;
        let mut column = vec![0.0; 2 * n];
        for i in 0..n {
            column[n + i] = product_three(acceleration[i], duration, period)?;
        }
        Ok((column, residual))
    }

    fn solve_rhs(
        &self,
        rounded_rhs: &[f64],
        original_rhs: &[Dd],
        column_scale: f64,
        tolerance: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<([f64; MAX_ROTOR_DOFS], f64), RotorError> {
        let n = self.model.dofs();
        let mut solution = [0.0; MAX_ROTOR_DOFS];
        for (i, value) in solution.iter_mut().enumerate().take(n) {
            *value = finite((rounded_rhs[i] / self.scale[i]) / column_scale)?;
            if rounded_rhs[i] != 0.0 && *value == 0.0 {
                return Err(RotorError::Invalid(
                    "periodic coefficient or force scaling underflow",
                ));
            }
        }
        poll(cancelled)?;
        self.factor.solve(&mut solution[..n]);
        poll(cancelled)?;
        let mut residual_norm = 0.0_f64;
        let mut action_norm = 0.0_f64;
        let mut rhs_norm = 0.0_f64;
        for i in 0..n {
            poll(cancelled)?;
            let row_scale = Dd::from_f64(self.scale[i]);
            let mut action = Dd::ZERO;
            for (j, &value) in solution.iter().enumerate().take(n) {
                finite(value)?;
                // Reconstruct directly from the original mass, independently
                // of the rounded Cholesky matrix and its triangular factors.
                let coefficient = Dd::from_f64(self.model.mass()[i * n + j])
                    / row_scale
                    / Dd::from_f64(self.scale[j]);
                action = action + coefficient * Dd::from_f64(value);
            }
            let rhs = original_rhs[i] / row_scale / Dd::from_f64(column_scale);
            residual_norm = residual_norm.max(finite((action - rhs).to_f64().abs())?);
            action_norm = action_norm.max(finite(action.to_f64().abs())?);
            rhs_norm = rhs_norm.max(finite(rhs.to_f64().abs())?);
        }
        let scale = action_norm.max(rhs_norm);
        let relative = if scale == 0.0 {
            residual_norm
        } else {
            finite((residual_norm / scale) / (action_norm / scale + rhs_norm / scale))?
        };
        if relative > tolerance {
            return Err(RotorError::Residual {
                relative,
                limit: tolerance,
            });
        }
        Ok((solution, relative))
    }
}

fn affine_segment_flow(
    generator: &[f64],
    source: &[f64],
    d: usize,
    budget: DenseExponentialBudget,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<AffineSegmentMap, RotorError> {
    poll(cancelled)?;
    let source_scale = source
        .iter()
        .fold(0.0_f64, |scale, value| scale.max(value.abs()));
    if source_scale == 0.0 {
        let transition = exponential_with_budget(generator, d, budget)
            .map_err(RotorError::PeriodicExponential)?;
        poll(cancelled)?;
        return Ok(AffineSegmentMap {
            transition,
            offset: vec![0.0; d],
        });
    }
    // Normalizing only the source coordinate prevents a large/small load
    // amplitude from dictating the exponential's scaling and work. It adds
    // no physical mode or fictitious stored energy.
    let e = d + 1;
    let mut augmented = vec![0.0; e * e];
    for i in 0..d {
        poll(cancelled)?;
        augmented[i * e..i * e + d].copy_from_slice(&generator[i * d..(i + 1) * d]);
        let normalized = finite(source[i] / source_scale)?;
        if source[i] != 0.0 && normalized == 0.0 {
            return Err(RotorError::Invalid("periodic force direction underflow"));
        }
        augmented[i * e + d] = normalized;
    }
    poll(cancelled)?;
    let result =
        exponential_with_budget(&augmented, e, budget).map_err(RotorError::PeriodicExponential)?;
    poll(cancelled)?;
    let mut transition = vec![0.0; d * d];
    let mut offset = Vec::with_capacity(d);
    for i in 0..d {
        poll(cancelled)?;
        transition[i * d..(i + 1) * d].copy_from_slice(&result[i * e..i * e + d]);
        offset.push(product_three(result[i * e + d], source_scale, 1.0)?);
    }
    Ok(AffineSegmentMap { transition, offset })
}

fn multiply(
    left: &[f64],
    right: &[f64],
    n: usize,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Vec<f64>, RotorError> {
    let mut result = vec![0.0; n * n];
    for i in 0..n {
        poll(cancelled)?;
        for j in 0..n {
            let mut value = 0.0_f64;
            for k in 0..n {
                value = left[i * n + k].mul_add(right[k * n + j], value);
            }
            result[i * n + j] = finite(value)?;
        }
    }
    Ok(result)
}

// Avoid avoidable intermediate overflow/underflow in h*T*coefficient.
// This changes only arithmetic association, never the original coefficient.
fn product_three(a: f64, b: f64, c: f64) -> Result<f64, RotorError> {
    if a == 0.0 {
        return Ok(a);
    }
    for (x, y, z) in [(a, b, c), (a, c, b), (b, c, a)] {
        let pair = x * y;
        if pair.is_finite() && pair != 0.0 {
            let value = pair * z;
            if value.is_finite() && value != 0.0 {
                return Ok(value);
            }
        }
    }
    Err(RotorError::Invalid(
        "periodic generator scaling is unrepresentable",
    ))
}

// S^-1 P S for the original [q,v] state; zero entries stay exactly zero.
fn restore_entry(value: f64, row_scale: f64, column_scale: f64) -> Result<f64, RotorError> {
    if value == 0.0 {
        return Ok(value);
    }
    let ratio = column_scale / row_scale;
    if ratio.is_finite() && ratio > 0.0 {
        let restored = value * ratio;
        if restored.is_finite() && restored != 0.0 {
            return Ok(restored);
        }
    }
    let product = value * column_scale;
    if product.is_finite() && product != 0.0 {
        let restored = product / row_scale;
        if restored.is_finite() && restored != 0.0 {
            return Ok(restored);
        }
    }
    let quotient = value / row_scale;
    if quotient.is_finite() && quotient != 0.0 {
        let restored = quotient * column_scale;
        if restored.is_finite() && restored != 0.0 {
            return Ok(restored);
        }
    }
    Err(RotorError::Invalid(
        "original periodic state map is unrepresentable",
    ))
}
