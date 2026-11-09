//! Complete hereditary rotor period maps, without material-state elimination.
//!
//! M q'' + (C + Omega G) q' + K_relaxed q + sum(b_j sigma_j) = f,
//! sigma_j' + sigma_j/tau_j = k_j b_j^T q'.
//!
//! All original [q, q', sigma] states remain continuous at C/K/clock switches.
//! Branch stiffness and strain directions, mass, gyroscopic matrix and signed
//! spin stay fixed. Optional held kinetics change tau from its original reference
//! on each segment; changing only this rate transfers no stored energy.
//! No equivalent viscous damping, material-pole cancellation or state reset
//! replaces the hereditary law. Changing K can exchange modulation work;
//! these affine maps do not label that exchange as material heat.

use super::{
    AffineSegmentMap, CycleForcing, MAX_ROTOR_PERIODIC_SEGMENTS, MassSolve, RotorBearingModel,
    RotorError, RotorFloquetBudget, RotorPeriodicSegment, affine_segment_flow, finite,
    forcing as affine, multiply, poll, product_three, restore_entry,
};
use crate::rotordynamics::relaxation::{MAX_ROTOR_MAXWELL_BRANCHES, MaxwellDampingBranch};
use fs_la::eigen_complex::{EigError, eig_with_cancel};
use fs_math::{c64::C64, dd::Dd, det};

mod forcing;
pub use forcing::{RotorMaxwellPeriodicOrbit, prepare_maxwell_forced_cycle};
mod harmonic;
pub use harmonic::{prepare_maxwell_fourier_cycle, prepare_maxwell_harmonic_cycle};
mod kinetics;
pub use kinetics::RotorMaxwellKinetics;
mod interval;
use interval::Forcing;
pub use interval::{
    PreparedRotorMaxwellInterval, RotorIntervalForcing, RotorMaxwellIntervalResponse,
    RotorMaxwellIntervalTrajectory, prepare_maxwell_interval,
};
pub mod power;

#[cfg(test)]
mod tests;

/// Physical homogeneous period map of every mechanical and material state.
#[derive(Clone, Debug, PartialEq)]
pub struct RotorMaxwellFloquetResponse {
    /// Constant signed shaft speed used in every C + Omega G [rad/s].
    pub spin_rad_s: f64,
    /// Sum of the original segment durations [s].
    pub period_seconds: f64,
    /// Number n of original generalized displacement coordinates.
    pub dofs: usize,
    /// Number b of original scalar stresses, in supplied branch order.
    pub memory_states: usize,
    /// Original-coordinate row-major (2n+b)-square map in [q, q_dot, sigma]
    /// order. Stresses are physical conjugates to b_j^T q, not strain or
    /// energy-scaled coordinates. Dark and equal-rate branches remain here.
    pub monodromy: Vec<f64>,
    /// All 2n+b numerical multipliers in the native canonical order.
    /// They have no eigenvalue error enclosure or stability certificate.
    pub multipliers: Vec<C64>,
    /// Largest computed multiplier magnitude, without clipping at unity.
    pub spectral_radius: f64,
    /// ln(spectral_radius) / period [1/s], retaining positive growth.
    pub growth_rate_per_s: f64,
    /// Maximum independent original full-mass solve residual, including
    /// every original branch and applied-force column. Not a trajectory error bound.
    pub maximum_mass_residual: f64,
    /// Maximum independent original memory-force/constitutive generator
    /// residual after coordinate/time scaling. Not an exponential bound.
    pub maximum_memory_residual: f64,
}

/// Actual displacement, velocity and stress of each individually retained arm.
#[derive(Clone, Debug, PartialEq)]
pub struct RotorMaxwellState {
    /// Original generalized displacement q.
    pub displacement: Vec<f64>,
    /// Original generalized velocity q_dot [q/s].
    pub velocity: Vec<f64>,
    /// Original conjugate stress sigma_j, in supplied branch order.
    /// An arbitrary initial memory is admissible; no relaxed initial state
    /// or strain-proportional stress is silently imposed.
    pub branch_stress: Vec<f64>,
}

/// Original physical state at the start and every coefficient boundary.
#[derive(Clone, Debug, PartialEq)]
pub struct RotorMaxwellPeriodicTrajectory {
    /// Physical boundary times. Cycles begin at zero; open intervals retain
    /// their authored absolute start and end times.
    pub boundary_times_seconds: Vec<f64>,
    /// Continuous [q, q_dot, sigma] values; a switch does not reset stresses.
    pub boundary_states: Vec<RotorMaxwellState>,
    /// Independently checked original increment equation
    /// (x_end-x_start) - (P-I)x_start - b = 0, normalized row by row.
    pub relative_map_residual: f64,
}

/// Immutable full-state cycle, ready for repeated arbitrary-memory advances.
///
/// Every warm call uses cached segment maps. Native exponentials and full-mass
/// solves occur only during preparation. This plan owns its complete maps and
/// cannot be changed by subsequently editing the source models or branches.
#[derive(Clone, Debug)]
pub struct PreparedRotorMaxwellCycle {
    report: RotorMaxwellFloquetResponse,
    flow: MaxwellFlow,
}

// One retained flow owner serves periodic analysis and open physical windows.
// Only the periodic wrapper adds a spectrum or fixed-point solve.
#[derive(Clone, Debug)]
struct MaxwellFlow {
    state_scale: Vec<f64>,
    monodromy: Vec<f64>,
    offset: Vec<f64>,
    physical_offset: Vec<f64>,
    segment_maps: Vec<AffineSegmentMap>,
    boundary_times: Vec<f64>,
    phase_times: Vec<Dd>,
}

/// Analyze a homogeneous periodic rotor with its original Maxwell memories.
///
/// The shared segment type supplies RELAXED K and complete C. Do not also add
/// the unrelaxed branch stiffness to K. M/G and signed spin must be constant;
/// branch strain vectors, stiffnesses and relaxation times are immutable for
/// the entire cycle. Every one of the original 2n+b states is retained, even
/// when its strain vanishes or multiple arms have the same relaxation time.
///
/// The same bounded native exponential and QR owners used by periodic_floquet
/// operate on the augmented state. No frozen-pole stability claim replaces
/// the chronologically composed map. Cancellation is polled around each native
/// exponential and during assembly, products and QR. A single exponential is
/// an indivisible, admitted call of at most 128 states, 32 Taylor products and
/// 48 squarings; this is not an operating-system latency bound.
///
/// The finite tolerance in (0,1) gates both original mass and memory rows.
/// A zero/unrepresentable multiplier refuses; near-unit rounding is retained
/// without classifying stability. Smooth modulation requires independent
/// segment refinement. This API gives numerical estimates, not enclosures.
pub fn periodic_maxwell_floquet(
    spin_rad_s: f64,
    segments: &[RotorPeriodicSegment<'_>],
    branches: &[MaxwellDampingBranch],
    relative_generator_residual_tolerance: f64,
    budget: RotorFloquetBudget,
    cancelled: impl FnMut() -> bool,
) -> Result<RotorMaxwellFloquetResponse, RotorError> {
    Ok(prepare(
        spin_rad_s,
        segments,
        branches,
        None,
        CycleForcing::Homogeneous,
        relative_generator_residual_tolerance,
        budget,
        false,
        cancelled,
    )?
    .report)
}

/// Prepare the same complete Floquet map and retain each boundary transition.
///
/// Physical conventions and numerical limits are those of
/// [`periodic_maxwell_floquet`]. The byte/work budget additionally reserves
/// every retained full-state segment map and one complete boundary trajectory.
pub fn prepare_maxwell_cycle(
    spin_rad_s: f64,
    segments: &[RotorPeriodicSegment<'_>],
    branches: &[MaxwellDampingBranch],
    relative_generator_residual_tolerance: f64,
    budget: RotorFloquetBudget,
    cancelled: impl FnMut() -> bool,
) -> Result<PreparedRotorMaxwellCycle, RotorError> {
    prepare(
        spin_rad_s,
        segments,
        branches,
        None,
        CycleForcing::Homogeneous,
        relative_generator_residual_tolerance,
        budget,
        true,
        cancelled,
    )
}

#[allow(clippy::too_many_arguments)] // physical inputs and explicit retention admission
fn prepare(
    spin_rad_s: f64,
    segments: &[RotorPeriodicSegment<'_>],
    branches: &[MaxwellDampingBranch],
    kinetics: Option<RotorMaxwellKinetics<'_>>,
    forces: CycleForcing<'_>,
    tolerance: f64,
    budget: RotorFloquetBudget,
    retain_segments: bool,
    mut cancelled: impl FnMut() -> bool,
) -> Result<PreparedRotorMaxwellCycle, RotorError> {
    let interval = prepare_flow(
        spin_rad_s,
        0.0,
        segments,
        branches,
        kinetics,
        Forcing::periodic(forces),
        tolerance,
        budget,
        retain_segments,
        true,
        &mut cancelled,
    )?;
    interval.into_cycle(budget, retain_segments, &mut cancelled)
}

#[allow(clippy::too_many_arguments)]
fn prepare_flow(
    spin_rad_s: f64,
    start_seconds: f64,
    segments: &[RotorPeriodicSegment<'_>],
    branches: &[MaxwellDampingBranch],
    kinetics: Option<RotorMaxwellKinetics<'_>>,
    forces: Forcing<'_>,
    tolerance: f64,
    budget: RotorFloquetBudget,
    retain_segments: bool,
    spectrum: bool,
    mut cancelled: impl FnMut() -> bool,
) -> Result<PreparedRotorMaxwellInterval, RotorError> {
    poll(&mut cancelled)?;
    if segments.is_empty()
        || segments.len() > MAX_ROTOR_PERIODIC_SEGMENTS
        || branches.len() > MAX_ROTOR_MAXWELL_BRANCHES
        || !spin_rad_s.is_finite()
        || !start_seconds.is_finite()
        || !tolerance.is_finite()
        || tolerance <= 0.0
        || tolerance >= 1.0
    {
        return Err(RotorError::Invalid(
            "bounded nonempty Maxwell cycle, finite spin and tolerance in (0,1) required",
        ));
    }
    let reference = segments[0].model;
    let n = reference.dofs();
    let d = 2 * n + branches.len();
    forces.shape.admit(n, segments.len())?;
    if let Some(kinetics) = kinetics {
        kinetics.validate(segments.len(), branches, &mut cancelled)?;
    }
    let forced = forces.shape.source_coordinates() != 0;
    let admission = if spectrum {
        RotorFloquetBudget::admit_states
    } else {
        RotorFloquetBudget::admit_interval_states
    };
    let dense_budget = admission(
        budget,
        d,
        segments.len(),
        forces.shape.source_coordinates(),
        forces.shape.exponentials_per_segment(),
        retain_segments,
    )?;
    // Reuse the fixed-speed owner's physical branch admission. This does not
    // eliminate or merge any branch and performs no spectral computation.
    reference.at_speed_with_maxwell_memory(spin_rad_s, branches, &mut cancelled)?;
    let mut boundary_times = Vec::with_capacity(segments.len() + 1);
    boundary_times.push(start_seconds);
    let mut phase_times = Vec::with_capacity(segments.len() + 1);
    phase_times.push(Dd::from_f64(start_seconds));
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
                "Maxwell segments require identical M/G and positive finite durations",
            ));
        }
        period = period + Dd::from_f64(segment.duration_seconds);
        let exact_time = Dd::from_f64(start_seconds) + period;
        let time = finite(exact_time.to_f64())?;
        if !spectrum {
            let rounding = exact_time - Dd::from_f64(time);
            let relative =
                finite((rounding / Dd::from_f64(segment.duration_seconds)).to_f64())?.abs();
            if time <= boundary_times[boundary_times.len() - 1] || relative > tolerance {
                return Err(RotorError::Invalid(
                    "absolute interval time resolution exceeds the supplied tolerance",
                ));
            }
        }
        boundary_times.push(time);
        // Legacy cycles retain their original rounded relative-time phase.
        // Open windows keep exact authored start+elapsed phase internally.
        phase_times.push(if spectrum {
            Dd::from_f64(time)
        } else {
            exact_time
        });
    }
    let period_seconds = finite(period.to_f64())?;
    forces.admit_phases(
        &boundary_times,
        &phase_times,
        segments,
        tolerance,
        &mut cancelled,
    )?;
    let assembly = MaxwellGenerator::new(
        reference,
        branches,
        period_seconds,
        tolerance,
        &mut cancelled,
    )?;
    let mass = &assembly.mass;
    let state_scale = &assembly.state_scale;
    let mut maximum_mass_residual = assembly.maximum_mass_residual;
    let mut monodromy = vec![0.0; d * d];
    for i in 0..d {
        monodromy[i * d + i] = 1.0;
    }
    let zero = vec![0.0; d];
    let mut offset = vec![0.0; d];
    let mut maximum_memory_residual = 0.0_f64;
    let mut segment_maps = Vec::with_capacity(if retain_segments { segments.len() } else { 0 });
    for (index, segment) in segments.iter().enumerate() {
        poll(&mut cancelled)?;
        let time_scales = kinetics.map(|k| k.relaxation_time_scales[index].as_slice());
        let (generator, residual, memory_residual) =
            assembly.segment(segment, time_scales, spin_rad_s, tolerance, &mut cancelled)?;
        maximum_mass_residual = maximum_mass_residual.max(residual);
        maximum_memory_residual = maximum_memory_residual.max(memory_residual);
        let h = segment.duration_seconds;
        let (flow, force_residual) = match forces.shape {
            CycleForcing::Homogeneous => (
                affine_segment_flow(&generator, &zero, d, dense_budget, &mut cancelled)?,
                0.0,
            ),
            CycleForcing::Held(rows) => {
                let (mut column, residual) =
                    mass.force_column(rows[index], h, period_seconds, tolerance, &mut cancelled)?;
                // Applied generalized force acts on the original momentum
                // rows. Memory responds through the complete coupled flow.
                column.resize(d, 0.0);
                (
                    affine_segment_flow(&generator, &column, d, dense_budget, &mut cancelled)?,
                    residual,
                )
            }
            CycleForcing::Harmonic { rows, order } => {
                super::harmonic::segment_flow_for_states_with_clock(
                    &mass,
                    &generator,
                    d,
                    rows[index],
                    forces.clock(0, order),
                    phase_times[index],
                    h,
                    period_seconds,
                    tolerance,
                    dense_budget,
                    &mut cancelled,
                )?
            }
            CycleForcing::Fourier { rows, orders } => {
                super::fourier::segment_flow_for_states_with_clocks(
                    &mass,
                    &generator,
                    d,
                    rows[index],
                    orders
                        .iter()
                        .enumerate()
                        .map(|(index, &order)| forces.clock(index, order)),
                    phase_times[index],
                    h,
                    period_seconds,
                    tolerance,
                    dense_budget,
                    &mut cancelled,
                )?
            }
        };
        maximum_mass_residual = maximum_mass_residual.max(force_residual);
        // Chronological affine composition preserves all original stresses:
        // b_next = E_segment b + b_segment, P_next = E_segment P.
        if forced {
            offset =
                affine::apply_affine(&flow.transition, &flow.offset, &offset, d, &mut cancelled)?;
        }
        monodromy = multiply(&flow.transition, &monodromy, d, &mut cancelled)?;
        if retain_segments {
            segment_maps.push(flow);
        }
    }
    let balanced_monodromy = monodromy.clone();
    let mut physical_offset = Vec::with_capacity(d);
    for i in 0..d {
        poll(&mut cancelled)?;
        physical_offset.push(restore_entry(offset[i], state_scale[i], 1.0)?);
    }
    for i in 0..d {
        poll(&mut cancelled)?;
        for j in 0..d {
            monodromy[i * d + j] =
                restore_entry(monodromy[i * d + j], state_scale[i], state_scale[j])?;
        }
    }
    poll(&mut cancelled)?;
    Ok(PreparedRotorMaxwellInterval {
        report: RotorMaxwellIntervalResponse {
            spin_rad_s,
            start_seconds,
            duration_seconds: period_seconds,
            dofs: n,
            memory_states: branches.len(),
            transition: monodromy,
            maximum_mass_residual,
            maximum_memory_residual,
        },
        flow: MaxwellFlow {
            state_scale: assembly.state_scale,
            monodromy: balanced_monodromy,
            offset,
            physical_offset,
            segment_maps,
            boundary_times,
            phase_times,
        },
    })
}

impl PreparedRotorMaxwellInterval {
    fn into_cycle(
        self,
        budget: RotorFloquetBudget,
        retain_segments: bool,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<PreparedRotorMaxwellCycle, RotorError> {
        let Self { report, mut flow } = self;
        let d = flow.state_scale.len();
        let complex: Vec<_> = flow
            .monodromy
            .iter()
            .map(|&value| C64::new(value, 0.0))
            .collect();
        let multipliers = eig_with_cancel(
            &complex,
            d,
            budget.max_qr_sweeps_per_window,
            &mut *cancelled,
        )
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
                    "periodic memory multiplier is unresolved or underflowed",
                ));
            }
            spectral_radius = spectral_radius.max(magnitude);
        }
        let growth_rate_per_s = finite(det::ln(spectral_radius) / report.duration_seconds)?;
        if !retain_segments {
            flow.monodromy = Vec::new();
        }
        poll(cancelled)?;
        Ok(PreparedRotorMaxwellCycle {
            report: RotorMaxwellFloquetResponse {
                spin_rad_s: report.spin_rad_s,
                period_seconds: report.duration_seconds,
                dofs: report.dofs,
                memory_states: report.memory_states,
                monodromy: report.transition,
                multipliers,
                spectral_radius,
                growth_rate_per_s,
                maximum_mass_residual: report.maximum_mass_residual,
                maximum_memory_residual: report.maximum_memory_residual,
            },
            flow,
        })
    }
}

// Shared original hereditary assembly for both native motion and physical work.
// The same balanced generator prevents an accounting-only model from drifting
// away from the stresses actually carried through the cycle.
struct MaxwellGenerator<'a> {
    mass: MassSolve<'a>,
    branches: &'a [MaxwellDampingBranch],
    period: f64,
    state_scale: Vec<f64>,
    coupling: Vec<(Vec<f64>, Vec<f64>)>,
    maximum_mass_residual: f64,
}

impl<'a> MaxwellGenerator<'a> {
    fn new(
        model: &'a RotorBearingModel,
        branches: &'a [MaxwellDampingBranch],
        period: f64,
        tolerance: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, RotorError> {
        let mass = MassSolve::new(model, cancelled)?;
        let n = model.dofs();
        // z = [S q, T S v, T sigma_j/sqrt(k_j)], S_ii = sqrt(M_ii).
        // Branch energy stays sigma_j^2/(2 k_j), not equivalent damping.
        let mut state_scale = mass.scale.clone();
        for &scale in &mass.scale {
            state_scale.push(product_three(period, scale, 1.0)?);
        }
        let mut coupling = Vec::with_capacity(branches.len());
        let mut maximum_mass_residual = 0.0_f64;
        for branch in branches {
            poll(cancelled)?;
            let root_k = det::sqrt(branch.stiffness);
            let stress_scale = finite(period / root_k)?;
            if stress_scale == 0.0 {
                return Err(RotorError::Invalid(
                    "periodic memory coordinate scale underflow",
                ));
            }
            state_scale.push(stress_scale);
            let original: Vec<_> = branch.strain.iter().copied().map(Dd::from_f64).collect();
            let (acceleration, residual) =
                mass.solve_rhs(&branch.strain, &original, 1.0, tolerance, cancelled)?;
            maximum_mass_residual = maximum_mass_residual.max(residual);
            let mut force = Vec::with_capacity(n);
            let mut strain = Vec::with_capacity(n);
            for i in 0..n {
                force.push(product_three(-acceleration[i], root_k, 1.0)?);
                strain.push(restore_entry(branch.strain[i], mass.scale[i], root_k)?);
            }
            coupling.push((force, strain));
        }
        Ok(Self {
            mass,
            branches,
            period,
            state_scale,
            coupling,
            maximum_mass_residual,
        })
    }

    fn segment(
        &self,
        segment: &RotorPeriodicSegment<'_>,
        time_scales: Option<&[f64]>,
        spin_rad_s: f64,
        tolerance: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(Vec<f64>, f64, f64), RotorError> {
        let n = self.mass.model.dofs();
        let d = self.state_scale.len();
        let (mechanical, mass_residual) =
            self.mass
                .segment_generator(segment, spin_rad_s, self.period, tolerance, cancelled)?;
        let mut generator = vec![0.0; d * d];
        for i in 0..2 * n {
            generator[i * d..i * d + 2 * n]
                .copy_from_slice(&mechanical[i * 2 * n..(i + 1) * 2 * n]);
        }
        let h = segment.duration_seconds;
        for (j, (branch, (force, strain))) in self.branches.iter().zip(&self.coupling).enumerate() {
            poll(cancelled)?;
            let row = 2 * n + j;
            let time =
                kinetics::relaxation_time(branch.relaxation_time_s, time_scales.map(|row| row[j]))?;
            generator[row * d + row] = finite(-h / time)?;
            if generator[row * d + row] == 0.0 {
                return Err(RotorError::Invalid("periodic memory relaxation underflow"));
            }
            for i in 0..n {
                generator[(n + i) * d + row] = product_three(force[i], h, 1.0)?;
                generator[row * d + n + i] = product_three(strain[i], h, 1.0)?;
            }
        }
        let memory_residual = memory_residual(
            &self.mass,
            self.branches,
            time_scales,
            &generator,
            &self.state_scale,
            h,
            self.period,
            tolerance,
            cancelled,
        )?;
        Ok((generator, mass_residual, memory_residual))
    }
}

// Reconstruct from original M, b, k, tau and the actual coordinate scales.
// Thus neither a rounded Cholesky factor nor the assembled reciprocal pair
// can validate itself. The branch-force row uses a full-mass normwise residual;
// each constitutive entry is independently scaled, including uncoupled decay.
#[allow(clippy::too_many_arguments)]
fn memory_residual(
    mass: &MassSolve<'_>,
    branches: &[MaxwellDampingBranch],
    time_scales: Option<&[f64]>,
    generator: &[f64],
    scale: &[f64],
    h: f64,
    period: f64,
    tolerance: f64,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<f64, RotorError> {
    let n = mass.model.dofs();
    let d = scale.len();
    let mut maximum = 0.0_f64;
    for (j, branch) in branches.iter().enumerate() {
        poll(cancelled)?;
        let row = 2 * n + j;
        let force_factor = Dd::from_f64(period) / Dd::from_f64(scale[row]);
        let strain_factor =
            (Dd::from_f64(scale[row]) / Dd::from_f64(period)) * Dd::from_f64(branch.stiffness);
        let mut residual_norm = 0.0_f64;
        let mut action_norm = 0.0_f64;
        let mut rhs_norm = 0.0_f64;
        for i in 0..n {
            poll(cancelled)?;
            let mut action = Dd::ZERO;
            for k in 0..n {
                let coefficient = Dd::from_f64(mass.model.mass()[i * n + k])
                    / Dd::from_f64(mass.scale[i])
                    / Dd::from_f64(mass.scale[k]);
                action = action + coefficient * Dd::from_f64(generator[(n + k) * d + row]);
            }
            let original_force =
                scaled_original(branch.strain[i], force_factor, mass.scale[i])? * Dd::from_f64(h);
            residual_norm = residual_norm.max(finite((action + original_force).to_f64().abs())?);
            action_norm = action_norm.max(finite(action.to_f64().abs())?);
            rhs_norm = rhs_norm.max(finite(original_force.to_f64().abs())?);
            let original_strain =
                scaled_original(branch.strain[i], strain_factor, mass.scale[i])? * Dd::from_f64(h);
            maximum = maximum.max(relative_difference(
                Dd::from_f64(generator[row * d + n + i]),
                original_strain,
            )?);
        }
        let peak = action_norm.max(rhs_norm);
        maximum = maximum.max(if peak == 0.0 {
            residual_norm
        } else {
            finite((residual_norm / peak) / (action_norm / peak + rhs_norm / peak))?
        });
        let time =
            kinetics::relaxation_time(branch.relaxation_time_s, time_scales.map(|row| row[j]))?;
        maximum = maximum.max(relative_difference(
            Dd::from_f64(generator[row * d + row]),
            -Dd::from_f64(h) / Dd::from_f64(time),
        )?);
    }
    if maximum > tolerance {
        return Err(RotorError::Residual {
            relative: maximum,
            limit: tolerance,
        });
    }
    Ok(maximum)
}

fn scaled_original(value: f64, factor: Dd, divisor: f64) -> Result<Dd, RotorError> {
    if value == 0.0 {
        return Ok(Dd::ZERO);
    }
    let value = Dd::from_f64(value);
    let divisor = Dd::from_f64(divisor);
    for result in [
        (value * factor) / divisor,
        (value / divisor) * factor,
        (factor / divisor) * value,
    ] {
        if result.to_f64().is_finite() && result.to_f64() != 0.0 {
            return Ok(result);
        }
    }
    Err(RotorError::Invalid(
        "original periodic memory coefficient is unrepresentable",
    ))
}

fn relative_difference(actual: Dd, expected: Dd) -> Result<f64, RotorError> {
    let actual_norm = finite(actual.to_f64().abs())?;
    let expected_norm = finite(expected.to_f64().abs())?;
    let residual = finite((actual - expected).to_f64().abs())?;
    let peak = actual_norm.max(expected_norm);
    if peak == 0.0 {
        Ok(residual)
    } else {
        finite((residual / peak) / (actual_norm / peak + expected_norm / peak))
    }
}

impl PreparedRotorMaxwellCycle {
    /// Complete homogeneous map and all mechanical/material multipliers.
    pub fn floquet_report(&self) -> &RotorMaxwellFloquetResponse {
        &self.report
    }

    /// Advance arbitrary physical displacement, velocity and original memory.
    ///
    /// Each returned boundary is round-tripped into the next cached map, so
    /// the public stress is exactly the state carried across a C/K switch.
    /// The final independent DD check uses original physical increments;
    /// large unchanged coordinates do not set its residual denominator.
    /// A rejection or cancellation publishes no partial trajectory.
    pub fn advance(
        &self,
        initial: &RotorMaxwellState,
        relative_residual_tolerance: f64,
        cancelled: impl FnMut() -> bool,
    ) -> Result<RotorMaxwellPeriodicTrajectory, RotorError> {
        self.flow.advance(
            &self.report.monodromy,
            self.report.dofs,
            self.report.memory_states,
            initial,
            relative_residual_tolerance,
            cancelled,
        )
    }
}

impl MaxwellFlow {
    fn advance(
        &self,
        physical_transition: &[f64],
        n: usize,
        memories: usize,
        initial: &RotorMaxwellState,
        relative_residual_tolerance: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<RotorMaxwellPeriodicTrajectory, RotorError> {
        poll(&mut cancelled)?;
        if !relative_residual_tolerance.is_finite()
            || relative_residual_tolerance <= 0.0
            || relative_residual_tolerance >= 1.0
        {
            return Err(RotorError::Invalid(
                "periodic memory trajectory tolerance must be in (0,1)",
            ));
        }
        let values = state_values(initial, n, memories)?;
        let d = self.state_scale.len();
        let mut balanced = self.encode(&values, &mut cancelled)?;
        let mut states = Vec::with_capacity(self.segment_maps.len() + 1);
        states.push(initial.clone());
        for segment in &self.segment_maps {
            let candidate = affine::apply_affine(
                &segment.transition,
                &segment.offset,
                &balanced,
                d,
                &mut cancelled,
            )?;
            let mut physical = Vec::with_capacity(d);
            for (value, scale) in candidate.iter().zip(&self.state_scale) {
                poll(&mut cancelled)?;
                physical.push(restore_entry(*value, *scale, 1.0)?);
            }
            balanced = self.encode(&physical, &mut cancelled)?;
            states.push(RotorMaxwellState {
                displacement: physical[..n].to_vec(),
                velocity: physical[n..2 * n].to_vec(),
                branch_stress: physical[2 * n..].to_vec(),
            });
        }
        let final_values = state_values(&states[states.len() - 1], n, memories)?;
        let relative_map_residual = affine::affine_residual(
            physical_transition,
            &self.physical_offset,
            &values,
            &final_values,
            |row| {
                affine::row_activity(states.iter().map(|state| {
                    if row < n {
                        state.displacement[row]
                    } else if row < 2 * n {
                        state.velocity[row - n]
                    } else {
                        state.branch_stress[row - 2 * n]
                    }
                }))
            },
            relative_residual_tolerance,
            &mut cancelled,
        )?;
        poll(&mut cancelled)?;
        Ok(RotorMaxwellPeriodicTrajectory {
            boundary_times_seconds: self.boundary_times.clone(),
            boundary_states: states,
            relative_map_residual,
        })
    }

    fn encode(
        &self,
        values: &[f64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Vec<f64>, RotorError> {
        let mut result = Vec::with_capacity(values.len());
        for (&value, &scale) in values.iter().zip(&self.state_scale) {
            poll(cancelled)?;
            result.push(restore_entry(value, 1.0, scale)?);
        }
        Ok(result)
    }
}

fn state_values(
    state: &RotorMaxwellState,
    n: usize,
    memories: usize,
) -> Result<Vec<f64>, RotorError> {
    if state.displacement.len() != n
        || state.velocity.len() != n
        || state.branch_stress.len() != memories
        || state
            .displacement
            .iter()
            .chain(&state.velocity)
            .chain(&state.branch_stress)
            .any(|value| !value.is_finite())
    {
        return Err(RotorError::Invalid(
            "matching finite original rotor and memory state required",
        ));
    }
    Ok(state
        .displacement
        .iter()
        .chain(&state.velocity)
        .chain(&state.branch_stress)
        .copied()
        .collect())
}
