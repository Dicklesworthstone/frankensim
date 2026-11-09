//! Open physical-time windows with uninterrupted prescribed source phase.

use super::super::{
    CycleForcing, MAX_ROTOR_FORCE_TONES, MAX_ROTOR_PERIODIC_SEGMENTS, RotorError,
    RotorFloquetBudget, RotorFourierForce, RotorHarmonicForce, RotorPeriodicSegment,
    harmonic::SourceClock, poll,
};
use super::{
    MaxwellFlow, RotorMaxwellKinetics, RotorMaxwellPeriodicTrajectory, RotorMaxwellState,
    prepare_flow,
};
use crate::rotordynamics::{MAX_ROTOR_DOFS, relaxation::MaxwellDampingBranch};
use fs_math::dd::Dd;

/// Physical loads `DC + sum(c*cos(omega*t) + s*sin(omega*t))` at absolute t.
/// Frequencies are positive rad/s, independent of window length and shaft spin.
/// Sine signs express the direction of a rotating load. The supplied amplitudes
/// may change between coefficient segments without resetting source phase.
#[derive(Clone, Copy, Debug)]
pub enum RotorIntervalForcing<'a> {
    Homogeneous,
    Held(&'a [&'a [f64]]),
    Harmonic {
        rows: &'a [RotorHarmonicForce<'a>],
        angular_frequency_rad_s: f64,
    },
    Fourier {
        rows: &'a [RotorFourierForce<'a>],
        angular_frequencies_rad_s: &'a [f64],
    },
}

impl RotorIntervalForcing<'_> {
    /// Allocation-free bounded source admission, also for callers with no
    /// interval to advance. Physical phase/time admission requires an actual
    /// window and remains in its preparation owner.
    pub fn validate(
        &self,
        dofs: usize,
        segment_count: usize,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<(), RotorError> {
        poll(&mut cancelled)?;
        if !(1..=MAX_ROTOR_DOFS).contains(&dofs)
            || segment_count == 0
            || segment_count > MAX_ROTOR_PERIODIC_SEGMENTS
        {
            return Err(RotorError::Invalid(
                "bounded physical source coordinates and segments required",
            ));
        }
        Forcing::interval(*self)?.shape.admit(dofs, segment_count)?;
        poll(&mut cancelled)
    }
}

/// Original physical interval map, with no periodicity or stability claim.
#[derive(Clone, Debug, PartialEq)]
pub struct RotorMaxwellIntervalResponse {
    pub spin_rad_s: f64,
    pub start_seconds: f64,
    pub duration_seconds: f64,
    pub dofs: usize,
    pub memory_states: usize,
    /// Row-major original [q,v,sigma] homogeneous transition over this window.
    pub transition: Vec<f64>,
    pub maximum_mass_residual: f64,
    pub maximum_memory_residual: f64,
}

/// The shared physical trajectory, with ABSOLUTE boundary times for a window.
pub type RotorMaxwellIntervalTrajectory = RotorMaxwellPeriodicTrajectory;

/// Immutable open window. Advancement reuses its original full-memory maps;
/// it does not repeat the forcing phase, solve an orbit, or compute a spectrum.
#[derive(Clone, Debug)]
pub struct PreparedRotorMaxwellInterval {
    pub(super) report: RotorMaxwellIntervalResponse,
    pub(super) flow: MaxwellFlow,
}

impl PreparedRotorMaxwellInterval {
    pub fn interval_report(&self) -> &RotorMaxwellIntervalResponse {
        &self.report
    }
    pub fn forcing_offset(&self) -> &[f64] {
        &self.flow.physical_offset
    }
    pub fn boundary_times_seconds(&self) -> &[f64] {
        &self.flow.boundary_times
    }

    pub fn advance(
        &self,
        initial: &RotorMaxwellState,
        relative_residual_tolerance: f64,
        cancelled: impl FnMut() -> bool,
    ) -> Result<RotorMaxwellIntervalTrajectory, RotorError> {
        self.flow.advance(
            &self.report.transition,
            self.report.dofs,
            self.report.memory_states,
            initial,
            relative_residual_tolerance,
            cancelled,
        )
    }
}

/// Prepare arbitrary positive-duration coefficient intervals at absolute time.
///
/// Physical frequencies remain fixed when the interval is subdivided. All q,
/// velocity and original branch stresses remain continuous. Optional kinetics
/// are ratios of the ORIGINAL branch times, with unchanged storage laws.
/// M/G and signed spin remain fixed; C/K and force amplitudes may be held per
/// segment. There is no implicit switch from the final model back to the first.
///
/// Each tone permits at most 1024 rotations per segment and absolute phase
/// |omega*t| <= 2^40 radians. DD multiplication and split-2*pi reduction retain
/// the phase of the original binary inputs within that explicit range. Time
/// boundaries must be finite and representably distinct. Their absolute
/// rounding error divided by the adjacent duration, and each tone's phase
/// handoff error in radians, must not exceed the supplied generator tolerance.
/// Internal source phases retain the DD start+elapsed time. Refusal publishes no
/// partial window. The inherited dense budget performs no QR work here.
#[allow(clippy::too_many_arguments)]
pub fn prepare_maxwell_interval(
    spin_rad_s: f64,
    start_seconds: f64,
    segments: &[RotorPeriodicSegment<'_>],
    branches: &[MaxwellDampingBranch],
    kinetics: Option<RotorMaxwellKinetics<'_>>,
    forcing: RotorIntervalForcing<'_>,
    relative_generator_residual_tolerance: f64,
    budget: RotorFloquetBudget,
    cancelled: impl FnMut() -> bool,
) -> Result<PreparedRotorMaxwellInterval, RotorError> {
    prepare_flow(
        spin_rad_s,
        start_seconds,
        segments,
        branches,
        kinetics,
        Forcing::interval(forcing)?,
        relative_generator_residual_tolerance,
        budget,
        true,
        false,
        cancelled,
    )
}

// Reuse the established amplitude/layout admission and source-count budgets.
// Synthetic integer orders describe SHAPE only; physical clocks below own phase.
#[derive(Clone, Copy)]
pub(super) struct Forcing<'a> {
    pub(super) shape: CycleForcing<'a>,
    physical: Option<PhysicalFrequencies<'a>>,
}

#[derive(Clone, Copy)]
enum PhysicalFrequencies<'a> {
    One(f64),
    Many(&'a [f64]),
}

impl<'a> Forcing<'a> {
    pub(super) fn periodic(shape: CycleForcing<'a>) -> Self {
        Self {
            shape,
            physical: None,
        }
    }

    pub(super) fn interval(forcing: RotorIntervalForcing<'a>) -> Result<Self, RotorError> {
        const ORDERS: [usize; MAX_ROTOR_FORCE_TONES] = [
            1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24,
            25, 26, 27, 28, 29, 30, 31, 32,
        ];
        let (shape, physical) = match forcing {
            RotorIntervalForcing::Homogeneous => (CycleForcing::Homogeneous, None),
            RotorIntervalForcing::Held(rows) => (CycleForcing::Held(rows), None),
            RotorIntervalForcing::Harmonic {
                rows,
                angular_frequency_rad_s: omega,
            } => {
                if !omega.is_finite() || omega <= 0.0 {
                    return Err(RotorError::Invalid(
                        "positive finite physical frequency required",
                    ));
                }
                (
                    CycleForcing::Harmonic { rows, order: 1 },
                    Some(PhysicalFrequencies::One(omega)),
                )
            }
            RotorIntervalForcing::Fourier {
                rows,
                angular_frequencies_rad_s: frequencies,
            } => {
                if frequencies.len() > MAX_ROTOR_FORCE_TONES
                    || frequencies.iter().any(|w| !w.is_finite() || *w <= 0.0)
                    || frequencies.windows(2).any(|w| w[0] >= w[1])
                {
                    return Err(RotorError::Invalid(
                        "bounded strictly increasing positive physical frequencies required",
                    ));
                }
                (
                    CycleForcing::Fourier {
                        rows,
                        orders: &ORDERS[..frequencies.len()],
                    },
                    Some(PhysicalFrequencies::Many(frequencies)),
                )
            }
        };
        Ok(Self { shape, physical })
    }

    pub(super) fn clock(self, index: usize, order: usize) -> SourceClock {
        match self.physical {
            Some(PhysicalFrequencies::One(omega)) => SourceClock::Physical(omega),
            Some(PhysicalFrequencies::Many(values)) => SourceClock::Physical(values[index]),
            None => SourceClock::Periodic(order),
        }
    }

    pub(super) fn admit_phases(
        self,
        times: &[f64],
        phase_times: &[Dd],
        segments: &[RotorPeriodicSegment<'_>],
        tolerance: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), RotorError> {
        let Some(frequencies) = self.physical else {
            return Ok(());
        };
        let single;
        let values = match frequencies {
            PhysicalFrequencies::One(omega) => {
                single = [omega];
                &single[..]
            }
            PhysicalFrequencies::Many(values) => values,
        };
        for (index, segment) in segments.iter().enumerate() {
            poll(cancelled)?;
            for &omega in values {
                poll(cancelled)?;
                for boundary in [index, index + 1] {
                    SourceClock::Physical(omega).phase(
                        phase_times[boundary],
                        segment.duration_seconds,
                        1.0,
                    )?;
                    let error = (Dd::from_f64(omega)
                        * (phase_times[boundary] - Dd::from_f64(times[boundary])))
                    .to_f64();
                    if !error.is_finite() || error.abs() > tolerance {
                        return Err(RotorError::Invalid(
                            "physical source phase handoff exceeds the supplied tolerance",
                        ));
                    }
                }
            }
        }
        Ok(())
    }
}
