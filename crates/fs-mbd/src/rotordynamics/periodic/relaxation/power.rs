//! Original source work, reversible storage and individual Maxwell heats.
//!
//! The full native [q, v, sigma] trajectory owns motion. Independently cached
//! quadratic integrals observe f.v, v.C.v, v.K_skew.q and sigma_j²/(k_j tau_j)
//! on that same generator, retaining all Fourier interference. Endpoint storage
//! is kinetic + relaxed elastic + sum sigma_j²/(2 k_j). A coefficient switch
//! exchanges q.(K_next-K_old).q/2; neither it nor signed active-bearing work is
//! relabeled as heat. Gyroscopic work is identically zero, with G still in motion.

pub use super::super::power::RotorCyclePowerError;
use super::super::power::{
    admit_energy_tolerance, check_energy, condense_sources, insert_source, represented, rounded,
    storage_change, switch_work,
};
use super::super::{
    CycleForcing, MAX_ROTOR_PERIODIC_SEGMENTS, RotorBearingModel, RotorError, RotorFloquetBudget,
    RotorFourierForce, RotorHarmonicForce, RotorPeriodicSegment, harmonic, poll, restore_entry,
};
use super::{
    MaxwellFlow, MaxwellGenerator, PreparedRotorMaxwellCycle, PreparedRotorMaxwellInterval,
    RotorIntervalForcing, RotorMaxwellKinetics, RotorMaxwellPeriodicOrbit,
    RotorMaxwellPeriodicTrajectory, RotorMaxwellState, interval::Forcing, kinetics, prepare_flow,
};
use crate::rotordynamics::relaxation::{MAX_ROTOR_MAXWELL_BRANCHES, MaxwellDampingBranch};
use fs_math::dd::Dd;
use fs_phs::linear_exact::{DenseExponentialBudget, quadratic_flow_integral_with_budget};

/// One original interval, its independent work channels and outgoing K switch.
#[derive(Clone, Debug, PartialEq)]
pub struct RotorMaxwellSegmentWork {
    pub applied_work_j: f64,
    /// Signed integral v.C.v: active bearing supply remains negative work.
    pub velocity_force_work_j: f64,
    /// Signed integral v.K_skew.q, excluded from relaxed elastic storage.
    pub nonconservative_stiffness_work_j: f64,
    pub kinetic_energy_change_j: f64,
    pub relaxed_elastic_energy_change_j: f64,
    /// Original branch order, including dark and equal-rate branches.
    pub branch_elastic_energy_change_j: Vec<f64>,
    /// Individual positive-law integrals sigma_j²/(k_j tau_j), not a loss
    /// inferred from total storage change. A negative numerical heat refuses.
    pub branch_heat_j: Vec<f64>,
    pub stiffness_switch_work_j: f64,
    /// Delta all storage + signed C/K_skew work + branch heat - force work,
    /// before the outgoing stiffness switch.
    pub balance_residual_j: f64,
}

/// Complete original work/heat account over one declared time window.
#[derive(Clone, Debug, PartialEq)]
pub struct RotorMaxwellCyclePower {
    pub segments: Vec<RotorMaxwellSegmentWork>,
    pub applied_work_j: f64,
    pub velocity_force_work_j: f64,
    pub nonconservative_stiffness_work_j: f64,
    pub stiffness_switch_work_j: f64,
    /// Cycles use the first model on both endpoints, including their seam work.
    /// Open intervals use the first/last models and only internal switch work.
    pub kinetic_energy_change_j: f64,
    pub relaxed_elastic_energy_change_j: f64,
    pub branch_elastic_energy_change_j: Vec<f64>,
    pub branch_heat_j: Vec<f64>,
    /// Delta all storage + C/K_skew work + branch heat - force - switch work.
    pub balance_residual_j: f64,
    /// Sum of absolute segment defects; cancellation cannot hide a bad interval.
    pub absolute_segment_residual_j: f64,
    pub average_applied_power_w: f64,
    pub average_velocity_force_power_w: f64,
    pub average_nonconservative_stiffness_power_w: f64,
    pub average_stiffness_switch_power_w: f64,
    pub average_branch_heat_power_w: Vec<f64>,
}

struct SegmentIntegrals<'a> {
    model: &'a RotorBearingModel,
    source_scale: f64,
    // Applied, velocity, circulatory stiffness, then each original branch heat.
    forms: Vec<Vec<f64>>,
}

/// Native full-state motion with independently cached physical work channels.
/// Borrowed source laws cannot change while this immutable plan is in use.
pub struct PreparedRotorMaxwellPowerCycle<'a> {
    motion: PreparedRotorMaxwellCycle,
    accounting: MaxwellPowerIntegrals<'a>,
}

struct MaxwellPowerIntegrals<'a> {
    branches: &'a [MaxwellDampingBranch],
    integrals: Vec<SegmentIntegrals<'a>>,
}

/// Open-window work/heat, with storage evaluated at its actual endpoint models.
pub type RotorMaxwellIntervalPower = RotorMaxwellCyclePower;

/// Cached physical-time motion and all original power channels, without a
/// periodic seam, Floquet spectrum or periodic-orbit operation.
pub struct PreparedRotorMaxwellPowerInterval<'a> {
    motion: PreparedRotorMaxwellInterval,
    accounting: MaxwellPowerIntegrals<'a>,
}

/// Prepare actual source work and branch heat on an arbitrary physical window.
/// Frequencies/global phase follow prepare_maxwell_interval. Internal K switches
/// retain their work; the last model remains in effect at the end. Composing
/// separate windows with a changed boundary K requires that external switch's
/// explicit storage-work accounting. Pure clock changes have no switch work.
#[allow(clippy::too_many_arguments)]
pub fn prepare_maxwell_power_interval<'a>(
    spin_rad_s: f64,
    start_seconds: f64,
    segments: &[RotorPeriodicSegment<'a>],
    branches: &'a [MaxwellDampingBranch],
    kinetics: Option<RotorMaxwellKinetics<'_>>,
    forcing: RotorIntervalForcing<'_>,
    relative_generator_residual_tolerance: f64,
    budget: RotorFloquetBudget,
    mut cancelled: impl FnMut() -> bool,
) -> Result<PreparedRotorMaxwellPowerInterval<'a>, RotorError> {
    let (motion, accounting) = prepare_power_window(
        spin_rad_s,
        start_seconds,
        segments,
        branches,
        kinetics,
        Forcing::interval(forcing)?,
        relative_generator_residual_tolerance,
        budget,
        false,
        &mut cancelled,
    )?;
    Ok(PreparedRotorMaxwellPowerInterval { motion, accounting })
}

/// Prepare held physical loads with full hereditary storage and heat accounting.
///
/// The cycle must retain the original fixed M/G, signed spin and branch laws.
/// Each supplied K is relaxed. C and K may switch. This entry point uses original
/// branch times; RotorMaxwellKinetics supplies held dissipative-clock changes
/// without storage work. Changes in branch k/strain or mass require their own
/// storage/history-transfer and parameter-work law.
/// The caller's complete budget includes native motion, retained forms and
/// exact quadratic-flow preparation. No extra trajectory integrator is added.
#[allow(clippy::too_many_arguments)]
pub fn prepare_maxwell_forced_power_cycle<'a>(
    spin_rad_s: f64,
    segments: &[RotorPeriodicSegment<'a>],
    branches: &'a [MaxwellDampingBranch],
    held_forces: &[&[f64]],
    relative_generator_residual_tolerance: f64,
    budget: RotorFloquetBudget,
    cancelled: impl FnMut() -> bool,
) -> Result<PreparedRotorMaxwellPowerCycle<'a>, RotorError> {
    prepare_power_cycle(
        spin_rad_s,
        segments,
        branches,
        None,
        CycleForcing::Held(held_forces),
        relative_generator_residual_tolerance,
        budget,
        cancelled,
    )
}

/// Integrate actual continuously harmonic force work on the full transient.
/// Global source phase and every branch stress continue through segment edges.
/// Auxiliary source rotations carry no physical storage or Floquet multipliers.
#[allow(clippy::too_many_arguments)]
pub fn prepare_maxwell_harmonic_power_cycle<'a>(
    spin_rad_s: f64,
    segments: &[RotorPeriodicSegment<'a>],
    branches: &'a [MaxwellDampingBranch],
    forces: &[RotorHarmonicForce<'_>],
    harmonic_order: usize,
    relative_generator_residual_tolerance: f64,
    budget: RotorFloquetBudget,
    cancelled: impl FnMut() -> bool,
) -> Result<PreparedRotorMaxwellPowerCycle<'a>, RotorError> {
    prepare_power_cycle(
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
        cancelled,
    )
}

/// Integrate simultaneous DC/tone source work and each material heat together.
///
/// DC-tone and tone-tone interference survive coefficient/amplitude switches
/// and arbitrary initial memory. Summing separate tone powers would lose them.
/// The existing stable quadratic owner observes all source coordinates jointly,
/// then the existing source-phase congruence condenses them to one affine slot.
/// Cold work is O(segments*(3+branches)*(2n+branches+1+2tones)^3); retained forms
/// and warm work are O(segments*(3+branches)*(2n+branches+1)^2). All dense work
/// and payload are admitted. No timing, positivity of C/K, thermal destination,
/// smooth-modulation error, continuous peak power or attraction claim follows.
#[allow(clippy::too_many_arguments)]
pub fn prepare_maxwell_fourier_power_cycle<'a>(
    spin_rad_s: f64,
    segments: &[RotorPeriodicSegment<'a>],
    branches: &'a [MaxwellDampingBranch],
    forces: &[RotorFourierForce<'_>],
    harmonic_orders: &[usize],
    relative_generator_residual_tolerance: f64,
    budget: RotorFloquetBudget,
    cancelled: impl FnMut() -> bool,
) -> Result<PreparedRotorMaxwellPowerCycle<'a>, RotorError> {
    prepare_power_cycle(
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
        cancelled,
    )
}

#[allow(clippy::too_many_arguments)]
pub(super) fn prepare_power_cycle<'a>(
    spin_rad_s: f64,
    segments: &[RotorPeriodicSegment<'a>],
    branches: &'a [MaxwellDampingBranch],
    kinetics: Option<RotorMaxwellKinetics<'_>>,
    forces: CycleForcing<'_>,
    tolerance: f64,
    budget: RotorFloquetBudget,
    mut cancelled: impl FnMut() -> bool,
) -> Result<PreparedRotorMaxwellPowerCycle<'a>, RotorError> {
    let (motion, accounting) = prepare_power_window(
        spin_rad_s,
        0.0,
        segments,
        branches,
        kinetics,
        Forcing::periodic(forces),
        tolerance,
        budget,
        true,
        &mut cancelled,
    )?;
    Ok(PreparedRotorMaxwellPowerCycle {
        motion: motion.into_cycle(budget, true, &mut cancelled)?,
        accounting,
    })
}

#[allow(clippy::too_many_arguments)]
fn prepare_power_window<'a>(
    spin_rad_s: f64,
    start_seconds: f64,
    segments: &[RotorPeriodicSegment<'a>],
    branches: &'a [MaxwellDampingBranch],
    kinetics: Option<RotorMaxwellKinetics<'_>>,
    forces: Forcing<'_>,
    tolerance: f64,
    budget: RotorFloquetBudget,
    spectrum: bool,
    mut cancelled: impl FnMut() -> bool,
) -> Result<(PreparedRotorMaxwellInterval, MaxwellPowerIntegrals<'a>), RotorError> {
    poll(&mut cancelled)?;
    let Some(first) = segments.first() else {
        return Err(RotorError::Invalid("nonempty Maxwell power cycle required"));
    };
    if segments.len() > MAX_ROTOR_PERIODIC_SEGMENTS || branches.len() > MAX_ROTOR_MAXWELL_BRANCHES {
        return Err(RotorError::PeriodicBudget(
            "Maxwell power segment/branch limit",
        ));
    }
    let n = first.model.dofs();
    forces.shape.admit(n, segments.len())?;
    if let Some(kinetics) = kinetics {
        kinetics.validate(segments.len(), branches, &mut cancelled)?;
    }
    let sources = match forces.shape {
        CycleForcing::Held(_) => 1,
        CycleForcing::Harmonic { .. } => 2,
        CycleForcing::Fourier { orders, .. } => 1 + 2 * orders.len(),
        CycleForcing::Homogeneous => 0,
    };
    let d = 2 * n + branches.len();
    let w = d + 1;
    let e = d + sources;
    let channels = 3 + branches.len();
    let count = segments.len();
    // All public model/branch/segment/source caps precede these products.
    // Reserve concurrent Gramian workspace, all retained forms, full results,
    // reassembly and the complete native motion/query admission.
    let extra_bytes =
        ((32 + channels) * e * e + 8 * e + count * (channels * w * w + 128 + 8 * branches.len()))
            * 8;
    let extra_work =
        (e as u64).pow(3) * (channels as u64 * 256 + 16) * count as u64 + 4 * (d as u64).pow(3);
    let native_budget = RotorFloquetBudget {
        max_dense_bytes: budget
            .max_dense_bytes
            .checked_sub(extra_bytes)
            .ok_or(RotorError::PeriodicBudget("Maxwell power integral bytes"))?,
        max_multiply_terms: budget
            .max_multiply_terms
            .checked_sub(extra_work)
            .ok_or(RotorError::PeriodicBudget("Maxwell power integral work"))?,
        ..budget
    };
    let motion = prepare_flow(
        spin_rad_s,
        start_seconds,
        segments,
        branches,
        kinetics,
        forces,
        tolerance,
        native_budget,
        true,
        spectrum,
        &mut cancelled,
    )?;
    let period = motion.report.duration_seconds;
    let assembly = MaxwellGenerator::new(first.model, branches, period, tolerance, &mut cancelled)?;
    let gramian_budget = DenseExponentialBudget {
        max_dimension: e,
        max_dense_bytes: 24 * e * e * 8,
        max_multiply_terms: (e as u64).pow(3) * 256,
    };
    let single_order = match forces.shape {
        CycleForcing::Harmonic { order, .. } => [order],
        _ => [0],
    };
    let mut integrals = Vec::with_capacity(count);
    for (index, segment) in segments.iter().enumerate() {
        poll(&mut cancelled)?;
        let time_scales = kinetics.map(|k| k.relaxation_time_scales[index].as_slice());
        let (generator, _, _) =
            assembly.segment(segment, time_scales, spin_rad_s, tolerance, &mut cancelled)?;
        let mut augmented = vec![0.0; e * e];
        for i in 0..d {
            augmented[i * e..i * e + d].copy_from_slice(&generator[i * d..(i + 1) * d]);
        }
        let mut forms = vec![vec![0.0; e * e]; channels];
        let mut scales = vec![0.0; sources];
        let mut phases = vec![0.0; sources];
        let (steady, tones, orders): (Option<&[f64]>, &[RotorHarmonicForce<'_>], &[usize]) =
            match forces.shape {
                CycleForcing::Held(rows) => (Some(rows[index]), &[], &[]),
                CycleForcing::Harmonic { rows, .. } => {
                    (None, core::slice::from_ref(&rows[index]), &single_order)
                }
                CycleForcing::Fourier { rows, orders } => {
                    (Some(rows[index].steady), rows[index].harmonics, orders)
                }
                CycleForcing::Homogeneous => (None, &[], &[]),
            };
        let mut slot = 0;
        if let Some(force) = steady {
            let (mut column, _) = assembly.mass.force_column(
                force,
                segment.duration_seconds,
                period,
                tolerance,
                &mut cancelled,
            )?;
            let scale = column.iter().fold(0.0_f64, |a, b| a.max(b.abs()));
            let normalization = if scale == 0.0 { 1.0 } else { scale };
            for value in &mut column {
                *value = represented(Dd::from_f64(*value) / Dd::from_f64(normalization))?;
            }
            insert_source(
                &mut augmented,
                &mut forms[0],
                e,
                &motion.flow.state_scale,
                d,
                &column,
                force,
                normalization,
                segment.duration_seconds,
            )?;
            scales[0] = scale;
            phases[0] = 1.0;
            slot = 1;
        }
        for (tone, (&order, &force)) in orders.iter().zip(tones).enumerate() {
            poll(&mut cancelled)?;
            let source = harmonic::prepare_source_with_clock(
                &assembly.mass,
                force,
                forces.clock(tone, order),
                motion.flow.phase_times[index],
                segment.duration_seconds,
                period,
                tolerance,
                &mut cancelled,
            )?;
            let normalization = if source.scale == 0.0 {
                1.0
            } else {
                source.scale
            };
            for (j, amplitude) in [force.cosine, force.sine].into_iter().enumerate() {
                insert_source(
                    &mut augmented,
                    &mut forms[0],
                    e,
                    &motion.flow.state_scale,
                    d + slot + j,
                    &source.columns[j],
                    amplitude,
                    normalization,
                    segment.duration_seconds,
                )?;
                scales[slot + j] = source.scale;
                phases[slot + j] = source.phase[j];
            }
            augmented[(d + slot) * e + d + slot + 1] = -source.angle;
            augmented[(d + slot + 1) * e + d + slot] = source.angle;
            slot += 2;
        }
        let h = Dd::from_f64(segment.duration_seconds);
        let half = Dd::from_f64(0.5);
        for i in 0..n {
            poll(&mut cancelled)?;
            let vi = Dd::from_f64(motion.flow.state_scale[n + i]);
            for j in 0..n {
                let vj = Dd::from_f64(motion.flow.state_scale[n + j]);
                forms[1][(n + i) * e + n + j] =
                    represented(h * Dd::from_f64(segment.model.damping()[i * n + j]) / vi / vj)?;
                let skew = half
                    * (Dd::from_f64(segment.model.stiffness()[i * n + j])
                        - Dd::from_f64(segment.model.stiffness()[j * n + i]));
                let entry =
                    represented(half * h * skew / vi / Dd::from_f64(motion.flow.state_scale[j]))?;
                forms[2][(n + i) * e + j] = entry;
                forms[2][j * e + n + i] = entry;
            }
        }
        for (j, branch) in branches.iter().enumerate() {
            poll(&mut cancelled)?;
            let row = 2 * n + j;
            let scale = Dd::from_f64(motion.flow.state_scale[row]);
            let time =
                kinetics::relaxation_time(branch.relaxation_time_s, time_scales.map(|row| row[j]))?;
            forms[3 + j][row * e + row] = represented(
                h / scale / Dd::from_f64(branch.stiffness) / scale / Dd::from_f64(time),
            )?;
        }
        let maximum = scales.iter().fold(0.0_f64, |a, b| a.max(b.abs()));
        let source_scale = if maximum == 0.0 { 1.0 } else { maximum };
        for (phase, scale) in phases.iter_mut().zip(scales) {
            let normalized = represented(
                Dd::from_f64(scale) / Dd::from_f64(source_scale) * Dd::from_f64(*phase),
            )?;
            if scale != 0.0 && *phase != 0.0 && normalized == 0.0 {
                return Err(RotorError::Invalid(
                    "Maxwell power source-phase scaling underflow",
                ));
            }
            *phase = normalized;
        }
        for form in &mut forms {
            poll(&mut cancelled)?;
            if form.iter().any(|value| *value != 0.0) {
                let flow =
                    quadratic_flow_integral_with_budget(&augmented, form, e, 1.0, gramian_budget)
                        .map_err(RotorError::PeriodicExponential)?;
                poll(&mut cancelled)?;
                *form = flow.integral;
            }
            *form = condense_sources(form, d, &phases, &mut cancelled)?;
        }
        integrals.push(SegmentIntegrals {
            model: segment.model,
            source_scale,
            forms,
        });
    }
    poll(&mut cancelled)?;
    Ok((
        motion,
        MaxwellPowerIntegrals {
            branches,
            integrals,
        },
    ))
}

impl PreparedRotorMaxwellPowerCycle<'_> {
    pub fn motion_cycle(&self) -> &PreparedRotorMaxwellCycle {
        &self.motion
    }

    /// Advance arbitrary original memory and check independent physical energy.
    /// Every segment, sum of absolute segment defects, and full-cycle defect
    /// must satisfy the same caller-authored absolute energy allowance.
    pub fn advance(
        &self,
        initial: &RotorMaxwellState,
        relative_residual_tolerance: f64,
        energy_tolerance_j: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<(RotorMaxwellPeriodicTrajectory, RotorMaxwellCyclePower), RotorCyclePowerError>
    {
        admit_energy_tolerance(energy_tolerance_j)?;
        let trajectory =
            self.motion
                .advance(initial, relative_residual_tolerance, &mut cancelled)?;
        let power = self.accounting.account(
            &self.motion.flow,
            self.motion.report.period_seconds,
            &trajectory,
            true,
            energy_tolerance_j,
            &mut cancelled,
        )?;
        Ok((trajectory, power))
    }

    /// Account the unique native periodic particular orbit. Its original
    /// Floquet report remains authoritative about any computed instability;
    /// a solved orbit and energy closure do not prove attraction or passivity.
    pub fn periodic_orbit(
        &self,
        relative_residual_tolerance: f64,
        maximum_resolvent_gain: f64,
        energy_tolerance_j: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<(RotorMaxwellPeriodicOrbit, RotorMaxwellCyclePower), RotorCyclePowerError> {
        admit_energy_tolerance(energy_tolerance_j)?;
        let orbit = self.motion.periodic_orbit(
            relative_residual_tolerance,
            maximum_resolvent_gain,
            &mut cancelled,
        )?;
        let power = self.accounting.account(
            &self.motion.flow,
            self.motion.report.period_seconds,
            &orbit.trajectory,
            true,
            energy_tolerance_j,
            &mut cancelled,
        )?;
        Ok((orbit, power))
    }
}

impl PreparedRotorMaxwellPowerInterval<'_> {
    pub fn motion_interval(&self) -> &PreparedRotorMaxwellInterval {
        &self.motion
    }

    pub fn advance(
        &self,
        initial: &RotorMaxwellState,
        relative_residual_tolerance: f64,
        energy_tolerance_j: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<(RotorMaxwellPeriodicTrajectory, RotorMaxwellIntervalPower), RotorCyclePowerError>
    {
        admit_energy_tolerance(energy_tolerance_j)?;
        let trajectory =
            self.motion
                .advance(initial, relative_residual_tolerance, &mut cancelled)?;
        let power = self.accounting.account(
            &self.motion.flow,
            self.motion.report.duration_seconds,
            &trajectory,
            false,
            energy_tolerance_j,
            &mut cancelled,
        )?;
        Ok((trajectory, power))
    }
}

impl MaxwellPowerIntegrals<'_> {
    fn account(
        &self,
        motion: &MaxwellFlow,
        duration: f64,
        trajectory: &RotorMaxwellPeriodicTrajectory,
        closed: bool,
        tolerance: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<RotorMaxwellCyclePower, RotorCyclePowerError> {
        let d = motion.state_scale.len();
        let w = d + 1;
        let mut totals = vec![Dd::ZERO; 3 + self.branches.len()];
        let mut total_switch = Dd::ZERO;
        let mut absolute = Dd::ZERO;
        let mut segments = Vec::with_capacity(self.integrals.len());
        let mut z = vec![0.0; w];
        for (index, cached) in self.integrals.iter().enumerate() {
            poll(cancelled)?;
            let before = &trajectory.boundary_states[index];
            let after = &trajectory.boundary_states[index + 1];
            for (i, (&value, &scale)) in before
                .displacement
                .iter()
                .chain(&before.velocity)
                .chain(&before.branch_stress)
                .zip(&motion.state_scale)
                .enumerate()
            {
                z[i] = restore_entry(value, 1.0, scale)?;
            }
            z[d] = cached.source_scale;
            let mut work = vec![Dd::ZERO; totals.len()];
            for (result, form) in work.iter_mut().zip(&cached.forms) {
                for i in 0..w {
                    poll(cancelled)?;
                    for j in 0..w {
                        *result = *result
                            + Dd::from_f64(z[i])
                                * Dd::from_f64(form[i * w + j])
                                * Dd::from_f64(z[j]);
                    }
                }
            }
            let kinetic = storage_change(cached.model.mass(), &before.velocity, &after.velocity);
            let elastic = storage_change(
                cached.model.stiffness(),
                &before.displacement,
                &after.displacement,
            );
            let mut branch_changes = Vec::with_capacity(self.branches.len());
            let mut branch_heat = Vec::with_capacity(self.branches.len());
            let mut memory_change = Dd::ZERO;
            let mut memory_heat = Dd::ZERO;
            for (j, branch) in self.branches.iter().enumerate() {
                poll(cancelled)?;
                let change =
                    branch_storage_change(branch, before.branch_stress[j], after.branch_stress[j])?;
                let heat = rounded(work[3 + j])?;
                if heat < 0.0 {
                    return Err(RotorError::Invalid(
                        "negative Maxwell heat is numerically unresolved",
                    )
                    .into());
                }
                branch_changes.push(rounded(change)?);
                branch_heat.push(heat);
                memory_change = memory_change + change;
                memory_heat = memory_heat + work[3 + j];
            }
            let next = if index + 1 < self.integrals.len() {
                self.integrals[index + 1].model
            } else if closed {
                self.integrals[0].model
            } else {
                cached.model
            };
            let switch = switch_work(
                cached.model.stiffness(),
                next.stiffness(),
                &after.displacement,
            );
            let residual = rounded(
                kinetic + elastic + memory_change + work[1] + work[2] + memory_heat - work[0],
            )?;
            check_energy(residual, tolerance)?;
            absolute = absolute + Dd::from_f64(residual.abs());
            for (total, &value) in totals.iter_mut().zip(&work) {
                *total = *total + value;
            }
            total_switch = total_switch + switch;
            segments.push(RotorMaxwellSegmentWork {
                applied_work_j: rounded(work[0])?,
                velocity_force_work_j: rounded(work[1])?,
                nonconservative_stiffness_work_j: rounded(work[2])?,
                kinetic_energy_change_j: rounded(kinetic)?,
                relaxed_elastic_energy_change_j: rounded(elastic)?,
                branch_elastic_energy_change_j: branch_changes,
                branch_heat_j: branch_heat,
                stiffness_switch_work_j: rounded(switch)?,
                balance_residual_j: residual,
            });
        }
        let first = &trajectory.boundary_states[0];
        let last = &trajectory.boundary_states[trajectory.boundary_states.len() - 1];
        let model = self.integrals[0].model;
        let kinetic = storage_change(model.mass(), &first.velocity, &last.velocity);
        let mut elastic =
            storage_change(model.stiffness(), &first.displacement, &last.displacement);
        if !closed {
            elastic = elastic
                + switch_work(
                    model.stiffness(),
                    self.integrals[self.integrals.len() - 1].model.stiffness(),
                    &last.displacement,
                );
        }
        let mut branch_changes = Vec::with_capacity(self.branches.len());
        let mut branch_heat = Vec::with_capacity(self.branches.len());
        let mut average_branch_heat = Vec::with_capacity(self.branches.len());
        let mut memory_change = Dd::ZERO;
        let mut memory_heat = Dd::ZERO;
        let period = Dd::from_f64(duration);
        for (j, branch) in self.branches.iter().enumerate() {
            poll(cancelled)?;
            let change =
                branch_storage_change(branch, first.branch_stress[j], last.branch_stress[j])?;
            branch_changes.push(rounded(change)?);
            branch_heat.push(rounded(totals[3 + j])?);
            average_branch_heat.push(rounded(totals[3 + j] / period)?);
            memory_change = memory_change + change;
            memory_heat = memory_heat + totals[3 + j];
        }
        let residual = rounded(
            kinetic + elastic + memory_change + totals[1] + totals[2] + memory_heat
                - totals[0]
                - total_switch,
        )?;
        let absolute_segment_residual_j = rounded(absolute)?;
        check_energy(residual, tolerance)?;
        check_energy(absolute_segment_residual_j, tolerance)?;
        poll(cancelled)?;
        Ok(RotorMaxwellCyclePower {
            segments,
            applied_work_j: rounded(totals[0])?,
            velocity_force_work_j: rounded(totals[1])?,
            nonconservative_stiffness_work_j: rounded(totals[2])?,
            stiffness_switch_work_j: rounded(total_switch)?,
            kinetic_energy_change_j: rounded(kinetic)?,
            relaxed_elastic_energy_change_j: rounded(elastic)?,
            branch_elastic_energy_change_j: branch_changes,
            branch_heat_j: branch_heat,
            balance_residual_j: residual,
            absolute_segment_residual_j,
            average_applied_power_w: rounded(totals[0] / period)?,
            average_velocity_force_power_w: rounded(totals[1] / period)?,
            average_nonconservative_stiffness_power_w: rounded(totals[2] / period)?,
            average_stiffness_switch_power_w: rounded(total_switch / period)?,
            average_branch_heat_power_w: average_branch_heat,
        })
    }
}

fn branch_storage_change(
    branch: &MaxwellDampingBranch,
    before: f64,
    after: f64,
) -> Result<Dd, RotorError> {
    let delta = Dd::from_f64(after) - Dd::from_f64(before);
    let sum = Dd::from_f64(after) + Dd::from_f64(before);
    if delta == Dd::ZERO || sum == Dd::ZERO {
        return Ok(Dd::ZERO);
    }
    let k = Dd::from_f64(branch.stiffness);
    let half = Dd::from_f64(0.5);
    // Keep the original k, independently of the motion's rounded sqrt(k).
    // Polarization retains small changes; ordering avoids premature squaring.
    for value in [
        half * (delta / k) * sum,
        half * (delta * sum) / k,
        half * delta * (sum / k),
    ] {
        if let Ok(result) = represented(value) {
            if result != 0.0 {
                return Ok(value);
            }
        }
    }
    Err(RotorError::Invalid(
        "Maxwell storage change is unrepresentable",
    ))
}
