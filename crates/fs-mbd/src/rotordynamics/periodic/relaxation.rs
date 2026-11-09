//! Complete hereditary rotor period maps, without material-state elimination.
//!
//! M q'' + (C + Omega G) q' + K_relaxed q + sum(b_j sigma_j) = f,
//! sigma_j' + sigma_j/tau_j = k_j b_j^T q'.
//!
//! All original [q, q', sigma] states remain continuous at C/K switches.
//! Branch definitions, mass, gyroscopic matrix and signed spin stay fixed.
//! No equivalent viscous damping, material-pole cancellation or state reset
//! replaces the hereditary law. Changing K can exchange modulation work;
//! these affine maps do not label that exchange as material heat.

use super::{
    AffineSegmentMap, CycleForcing, MAX_ROTOR_PERIODIC_SEGMENTS, MassSolve, RotorError,
    RotorFloquetBudget, RotorPeriodicSegment, affine_segment_flow, finite, forcing as affine,
    multiply, poll, product_three, restore_entry,
};
use crate::rotordynamics::relaxation::{MAX_ROTOR_MAXWELL_BRANCHES, MaxwellDampingBranch};
use fs_la::eigen_complex::{EigError, eig_with_cancel};
use fs_math::{c64::C64, dd::Dd, det};

mod forcing;
pub use forcing::{RotorMaxwellPeriodicOrbit, prepare_maxwell_forced_cycle};
mod harmonic;
pub use harmonic::{prepare_maxwell_fourier_cycle, prepare_maxwell_harmonic_cycle};

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
    /// Physical boundary times, beginning at zero and ending at one period.
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
    state_scale: Vec<f64>,
    monodromy: Vec<f64>,
    offset: Vec<f64>,
    physical_offset: Vec<f64>,
    segment_maps: Vec<AffineSegmentMap>,
    boundary_times: Vec<f64>,
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
    forces: CycleForcing<'_>,
    tolerance: f64,
    budget: RotorFloquetBudget,
    retain_segments: bool,
    mut cancelled: impl FnMut() -> bool,
) -> Result<PreparedRotorMaxwellCycle, RotorError> {
    poll(&mut cancelled)?;
    if segments.is_empty()
        || segments.len() > MAX_ROTOR_PERIODIC_SEGMENTS
        || branches.len() > MAX_ROTOR_MAXWELL_BRANCHES
        || !spin_rad_s.is_finite()
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
    forces.admit(n, segments.len())?;
    let forced = forces.source_coordinates() != 0;
    let dense_budget = budget.admit_states(
        d,
        segments.len(),
        forces.source_coordinates(),
        forces.exponentials_per_segment(),
        retain_segments,
    )?;
    // Reuse the fixed-speed owner's physical branch admission. This does not
    // eliminate or merge any branch and performs no spectral computation.
    reference.at_speed_with_maxwell_memory(spin_rad_s, branches, &mut cancelled)?;
    let mut boundary_times = Vec::with_capacity(segments.len() + 1);
    boundary_times.push(0.0);
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
        boundary_times.push(finite(period.to_f64())?);
    }
    let period_seconds = finite(period.to_f64())?;
    let mass = MassSolve::new(reference, &mut cancelled)?;
    // One common invertible balance for every segment:
    // z = [S q, T S v, T sigma_j/sqrt(k_j)], S_ii = sqrt(M_ii).
    // The reciprocal memory/velocity couplings then share sqrt(k_j) b_j.
    // Branch energy is sigma_j^2/(2 k_j); no effective damping is introduced.
    let mut state_scale = mass.scale.clone();
    for &scale in &mass.scale {
        state_scale.push(product_three(period_seconds, scale, 1.0)?);
    }
    let mut coupling = Vec::with_capacity(branches.len());
    let mut maximum_mass_residual = 0.0_f64;
    for branch in branches {
        poll(&mut cancelled)?;
        let root_k = det::sqrt(branch.stiffness);
        let stress_scale = finite(period_seconds / root_k)?;
        if stress_scale == 0.0 {
            return Err(RotorError::Invalid(
                "periodic memory coordinate scale underflow",
            ));
        }
        state_scale.push(stress_scale);
        let original: Vec<_> = branch.strain.iter().copied().map(Dd::from_f64).collect();
        let (acceleration, residual) =
            mass.solve_rhs(&branch.strain, &original, 1.0, tolerance, &mut cancelled)?;
        maximum_mass_residual = maximum_mass_residual.max(residual);
        let mut force = Vec::with_capacity(n);
        let mut strain = Vec::with_capacity(n);
        for i in 0..n {
            force.push(product_three(-acceleration[i], root_k, 1.0)?);
            strain.push(restore_entry(branch.strain[i], mass.scale[i], root_k)?);
        }
        coupling.push((force, strain));
    }
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
        let (mechanical, residual) = mass.segment_generator(
            segment,
            spin_rad_s,
            period_seconds,
            tolerance,
            &mut cancelled,
        )?;
        maximum_mass_residual = maximum_mass_residual.max(residual);
        let mut generator = vec![0.0; d * d];
        for i in 0..2 * n {
            generator[i * d..i * d + 2 * n]
                .copy_from_slice(&mechanical[i * 2 * n..(i + 1) * 2 * n]);
        }
        let h = segment.duration_seconds;
        for (j, (branch, (force, strain))) in branches.iter().zip(&coupling).enumerate() {
            poll(&mut cancelled)?;
            let row = 2 * n + j;
            generator[row * d + row] = finite(-h / branch.relaxation_time_s)?;
            if generator[row * d + row] == 0.0 {
                return Err(RotorError::Invalid("periodic memory relaxation underflow"));
            }
            for i in 0..n {
                generator[(n + i) * d + row] = product_three(force[i], h, 1.0)?;
                generator[row * d + n + i] = product_three(strain[i], h, 1.0)?;
            }
        }
        maximum_memory_residual = maximum_memory_residual.max(memory_residual(
            &mass,
            branches,
            &generator,
            &state_scale,
            h,
            period_seconds,
            tolerance,
            &mut cancelled,
        )?);
        let (flow, force_residual) = match forces {
            CycleForcing::Homogeneous => (
                affine_segment_flow(&generator, &zero, d, dense_budget, &mut cancelled)?,
                0.0,
            ),
            CycleForcing::Held(rows) => {
                let (mut column, residual) = mass.force_column(
                    rows[index], h, period_seconds, tolerance, &mut cancelled,
                )?;
                // Applied generalized force acts on the original momentum
                // rows. Memory responds through the complete coupled flow.
                column.resize(d, 0.0);
                (
                    affine_segment_flow(&generator, &column, d, dense_budget, &mut cancelled)?,
                    residual,
                )
            }
            CycleForcing::Harmonic { rows, order } => harmonic::segment_flow(
                &mass, &generator, d, rows[index], order, boundary_times[index],
                h, period_seconds, tolerance, dense_budget, &mut cancelled,
            )?,
            CycleForcing::Fourier { rows, orders } => harmonic::fourier_segment_flow(
                &mass, &generator, d, rows[index], orders, boundary_times[index],
                h, period_seconds, tolerance, dense_budget, &mut cancelled,
            )?,
        };
        maximum_mass_residual = maximum_mass_residual.max(force_residual);
        // Chronological affine composition preserves all original stresses:
        // b_next = E_segment b + b_segment, P_next = E_segment P.
        if forced {
            offset = affine::apply_affine(
                &flow.transition, &flow.offset, &offset, d, &mut cancelled,
            )?;
        }
        monodromy = multiply(&flow.transition, &monodromy, d, &mut cancelled)?;
        if retain_segments {
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
                "periodic memory multiplier is unresolved or underflowed",
            ));
        }
        spectral_radius = spectral_radius.max(magnitude);
    }
    let growth_rate_per_s = finite(det::ln(spectral_radius) / period_seconds)?;
    let balanced_monodromy = if retain_segments { monodromy.clone() } else { Vec::new() };
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
    Ok(PreparedRotorMaxwellCycle {
        report: RotorMaxwellFloquetResponse {
            spin_rad_s,
            period_seconds,
            dofs: n,
            memory_states: branches.len(),
            monodromy,
            multipliers,
            spectral_radius,
            growth_rate_per_s,
            maximum_mass_residual,
            maximum_memory_residual,
        },
        state_scale,
        monodromy: balanced_monodromy,
        offset,
        physical_offset,
        segment_maps,
        boundary_times,
    })
}

// Reconstruct from original M, b, k, tau and the actual coordinate scales.
// Thus neither a rounded Cholesky factor nor the assembled reciprocal pair
// can validate itself. The branch-force row uses a full-mass normwise residual;
// each constitutive entry is independently scaled, including uncoupled decay.
#[allow(clippy::too_many_arguments)]
fn memory_residual(
    mass: &MassSolve<'_>,
    branches: &[MaxwellDampingBranch],
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
        maximum = maximum.max(relative_difference(
            Dd::from_f64(generator[row * d + row]),
            -Dd::from_f64(h) / Dd::from_f64(branch.relaxation_time_s),
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
        let values = state_values(initial, self.report.dofs, self.report.memory_states)?;
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
            let n = self.report.dofs;
            states.push(RotorMaxwellState {
                displacement: physical[..n].to_vec(),
                velocity: physical[n..2 * n].to_vec(),
                branch_stress: physical[2 * n..].to_vec(),
            });
        }
        let final_values = state_values(
            &states[states.len() - 1],
            self.report.dofs,
            self.report.memory_states,
        )?;
        let relative_map_residual = affine::affine_residual(
            &self.report.monodromy,
            &self.physical_offset,
            &values,
            &final_values,
            |row| {
                let n = self.report.dofs;
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
