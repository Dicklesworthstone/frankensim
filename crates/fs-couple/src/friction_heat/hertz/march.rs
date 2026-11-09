//! Time-resolved Hertz friction, heat release and finite-cell cooling.
//!
//! One full phase and two half phases start from the same accepted checkpoint.
//! Only the two-half-phase native/thermal states advance. Their physical-unit
//! differences estimate local error in endpoint source sampling and backward
//! Euler transport; they are not certified enclosures or global error bounds.
//! The original laws, finite cells, source receipts and arrest events remain
//! authoritative. In particular, rest still conducts without creating heat.

use super::*;
use fs_contact::interface_binding::thermal_friction::rate_dependent::SourceHertzFrictionError;
use fs_math::dd::Dd;
use fs_mbd::tangential_friction::hertz::driven::{RigidHertzDrivenError, RigidHertzDrivenReport};
use fs_mbd::tangential_friction::hertz::rate_dependent::driven::{
    RigidHertzRateDrivenError, RigidHertzRateDrivenReport,
};
use fs_mbd::tangential_friction::hertz::rate_dependent::{
    RigidHertzRateTransientError, RigidHertzRateTransientReport,
};
use fs_mbd::tangential_friction::hertz::trajectory::{
    RigidHertzTrajectoryError, RigidHertzTrajectoryReport,
};
use fs_mbd::tangential_friction::hertz::transient::{
    HertzTransientControls, RigidHertzTransientError, RigidHertzTransientReport,
};
use fs_mbd::{DynamicsError, RigidKineticEnergyChange, Vec3};

#[derive(Clone, Copy, Debug)]
pub enum HertzHeatMarchModel {
    InvariantDirection,
    Curved(HertzTransientControls),
}

/// Original native reports, including any unresolved curved extinction tail.
pub enum HertzHeatMarchMotion {
    InvariantDirection(RigidHertzTrajectoryReport),
    Curved(RigidHertzTransientReport),
    RateDependent(RigidHertzRateTransientReport),
    Driven(RigidHertzDrivenReport),
    RateDriven(RigidHertzRateDrivenReport),
}
impl HertzHeatMarchMotion {
    fn travel(&self, radius: f64) -> [f64; 3] {
        let (slip, spin) = match self {
            Self::InvariantDirection(m) => (
                m.relaxation.slip_increment_m,
                m.relaxation.spin_increment_rad,
            ),
            Self::Curved(m)
            | Self::RateDependent(RigidHertzRateTransientReport { mechanical: m, .. }) => {
                (m.transient.slip_increment_m, m.transient.spin_increment_rad)
            }
            Self::Driven(m)
            | Self::RateDriven(RigidHertzRateDrivenReport { mechanical: m, .. }) => {
                (m.driven.slip_increment_m, m.driven.spin_increment_rad)
            }
        };
        [slip[0], slip[1], radius * spin]
    }
    // Travel, friction heat and supplied mechanical work, respectively.
    fn uncertainty(&self, radius: f64) -> [f64; 3] {
        match self {
            // Driven held-traction time resolution has a separate native
            // refinement control, not an embedded trajectory error estimate.
            Self::InvariantDirection(_) | Self::Driven(_) => [0.0; 3],
            Self::RateDriven(m) => {
                let travel = m.estimated_scaled_travel_error_m;
                // The original native supply is analytic free-body work plus
                // g dot (travel-free_travel), with g=W^-1*a. This estimate
                // cannot disappear through coarse/fine work cancellation.
                let work = if travel == 0.0 {
                    0.0
                } else {
                    norm(m.mechanical.driven.driving_scaled_wrench_n) * travel
                };
                [travel, m.estimated_heat_error_j, work]
            }
            Self::Curved(m)
            | Self::RateDependent(RigidHertzRateTransientReport { mechanical: m, .. }) => [
                m.transient.estimated_scaled_travel_error_m
                    + norm([
                        m.transient.tail_slip_bound_m[0],
                        m.transient.tail_slip_bound_m[1],
                        radius * m.transient.tail_spin_bound_rad,
                    ]),
                m.transient.estimated_heat_error_j,
                0.0,
            ],
        }
    }
    fn has_driven_arrest(&self) -> bool {
        match self {
            Self::Driven(m)
            | Self::RateDriven(RigidHertzRateDrivenReport { mechanical: m, .. }) => {
                m.driven.first_arrest_time_s.is_some()
            }
            _ => false,
        }
    }
    fn external_work_j(&self) -> f64 {
        match self {
            Self::Driven(m)
            | Self::RateDriven(RigidHertzRateDrivenReport { mechanical: m, .. }) => {
                m.external_work_j
            }
            _ => 0.0,
        }
    }
}

/// Absolute physical scales. Relative scales use phase increments, never
/// background body momentum or the cells' absolute formation energies.
#[derive(Clone, Copy, Debug)]
pub struct HertzHeatMarchTolerance {
    pub linear_momentum_ns: f64,
    pub angular_momentum_nms: f64,
    pub cell_energy_j: f64,
    /// Euclidean norm of [net slip x, net slip y, radius * relative rotation].
    pub scaled_travel_m: f64,
    /// Separate absolute error gate for friction heat AND external mechanical work.
    pub dissipated_j: f64,
    pub relative: f64,
}
#[derive(Clone, Copy, Debug)]
pub struct HertzHeatMarchControls {
    pub initial_step_s: f64,
    /// Smallest comparison interval, except a shorter final remainder.
    pub minimum_step_s: f64,
    pub maximum_step_s: f64,
    pub maximum_attempts: usize,
    /// Every accepted comparison retains exactly two original phase receipts.
    pub maximum_accepted_phases: usize,
    /// Includes discarded coarse phases and rejected comparisons.
    pub maximum_thermal_trials: usize,
    pub tolerance: HertzHeatMarchTolerance,
    /// Independent whole-span gate, not multiplied by the number of phases.
    pub energy_tolerance_j: f64,
}
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct HertzHeatMarchDefect {
    pub linear_momentum_ns: f64,
    pub angular_momentum_nms: f64,
    pub cell_energy_j: f64,
    pub scaled_travel_m: f64,
    pub dissipated_j: f64,
    pub external_mechanical_work_j: f64,
    /// Different retained modes imply different future constitutive branches.
    pub mode_disagreement: bool,
    /// Coarse/fine driven histories disagree whether an arrest occurred, even
    /// when a subsequent reversal leaves both endpoints in the sliding mode.
    pub arrest_disagreement: bool,
    pub normalized: f64,
}
pub struct HertzHeatMarchPhase<P, F = SourceHertzWrench> {
    pub start_s: f64,
    pub end_s: f64,
    pub mechanical: HertzHeatMarchMotion,
    pub heat: FrictionHeatNetworkReport,
    pub friction_source: F,
    pub source: P,
    pub nonlinear: NewtonReport,
    pub trials: usize,
}
/// Publish these original A/B momenta and the final thermal owner together.
/// No intermediate or discarded trial can modify the borrowed input owners.
pub struct HeatedRigidHertzMarchCandidate<'m, P, F = SourceHertzWrench> {
    pub binding: RigidContactBinding,
    pub bodies_after: [RigidBodyState; 2],
    pub thermal: DepositionSimulation<'m>,
    pub phases: Vec<HertzHeatMarchPhase<P, F>>,
    pub attempts: usize,
    pub rejected_attempts: usize,
    pub thermal_trials: usize,
    pub maximum_accepted_normalized_defect: f64,
    pub slip_increment_m: [f64; 2],
    pub spin_increment_rad: f64,
    pub dissipated_j: f64,
    pub body_heat_j: [f64; 2],
    pub external_heat_j: f64,
    /// Original continuous drive supply, including complementary body motion.
    pub external_mechanical_work_j: f64,
    pub kinetic_energy_change_j: f64,
    pub thermal_energy_change_j: f64,
    /// Mechanical, thermal, and combined physical accounts, respectively.
    pub energy_residual_j: [f64; 3],
}
#[derive(Debug)]
pub enum HertzHeatMarchError<E> {
    Invalid(&'static str),
    InvariantDirection(RigidHertzHeatError<E, RigidHertzTrajectoryError>),
    Curved(RigidHertzHeatError<E, RigidHertzTransientError>),
    RateDependent(RigidHertzHeatError<E, RigidHertzRateTransientError, SourceHertzFrictionError>),
    Driven(RigidHertzHeatError<E, RigidHertzDrivenError>),
    RateDriven(RigidHertzHeatError<E, RigidHertzRateDrivenError, SourceHertzFrictionError>),
    Dynamics(DynamicsError),
    Cancelled,
    Budget(&'static str),
    Unresolved {
        step_s: f64,
        defect: HertzHeatMarchDefect,
    },
    Energy {
        residual_j: [f64; 3],
    },
}
impl<E: core::fmt::Debug> core::fmt::Display for HertzHeatMarchError<E> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "Hertz friction/heat march: {self:?}")
    }
}
impl<E: std::error::Error + 'static> std::error::Error for HertzHeatMarchError<E> {}

pub(super) struct Stage<'m, P, F, H> {
    bodies: [RigidBodyState; 2],
    thermal: DepositionSimulation<'m>,
    receipt: HertzHeatMarchPhase<P, F>,
    history: H,
}
pub(super) fn stage<'m, P, M, F, H>(
    c: HeatedRigidHertzCandidate<'m, P, M, F>,
    start_s: f64,
    motion: fn(M) -> HertzHeatMarchMotion,
    history: H,
) -> Stage<'m, P, F, H> {
    Stage {
        bodies: c.bodies_after,
        thermal: c.thermal,
        history,
        receipt: HertzHeatMarchPhase {
            start_s,
            end_s: c.target_time_s,
            mechanical: motion(c.mechanical),
            heat: c.heat,
            friction_source: c.friction_source,
            source: c.source,
            nonlinear: c.nonlinear,
            trials: c.trials,
        },
    }
}

/// Refine one held-geometry span with original source/thermal feedback.
///
/// Coefficients refresh at each solved half-phase temperature. Physical heat
/// release and subsequent cooling therefore resolve in time instead of being
/// averaged across the complete requested horizon. Callbacks must be pure
/// candidate evaluations; their payloads survive only for accepted phases.
/// Any source, event, numerical, budget or cancellation refusal leaves the
/// entire caller state untouched, including after successful private substeps.
///
/// Poses, patch geometry, normal load and radius remain fixed. Pending thermal
/// events must be cut by their original owner. This adds no source history,
/// moving footprint, normal evolution, drive, or radial flash-temperature law.
/// Native mechanical/quadrature tolerances remain separate accuracy controls;
/// reported curved travel/tail and heat estimates also enter the comparison.
#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
pub fn propose_rigid_hertz_heat_march<'m, P, E>(
    bodies: &[RigidBodyState],
    properties: &[MassProperties],
    contact: RigidHertzHeatContact<'_>,
    thermal: &DepositionSimulation<'m>,
    route: &FrictionHeatRoute,
    duration_s: f64,
    phase_controls: RigidHertzHeatControls,
    model: HertzHeatMarchModel,
    controls: HertzHeatMarchControls,
    temperature_guess_k: f64,
    resolve: impl Fn(&DepositionSimulation<'m>, f64) -> Result<RigidHertzHeatSource<P>, E>,
    gate: &CancelGate,
) -> Result<HeatedRigidHertzMarchCandidate<'m, P>, HertzHeatMarchError<E>> {
    propose_hertz_heat_phases(
        bodies,
        properties,
        contact,
        thermal,
        (),
        duration_s,
        phase_controls,
        controls,
        temperature_guess_k,
        |b, props, pair, h, _, dt, local, energy, guess| {
            let start = h.time_s();
            Ok(match model {
                HertzHeatMarchModel::InvariantDirection => stage(
                    trajectory::propose_rigid_hertz_trajectory_heat(
                        b, props, pair, h, route, start, dt, local, energy, guess, &resolve, gate,
                    )
                    .map_err(HertzHeatMarchError::InvariantDirection)?,
                    start,
                    HertzHeatMarchMotion::InvariantDirection,
                    (),
                ),
                HertzHeatMarchModel::Curved(path) => stage(
                    transient::propose_rigid_hertz_transient_heat(
                        b, props, pair, h, route, start, dt, local, path, energy, guess, &resolve,
                        gate,
                    )
                    .map_err(HertzHeatMarchError::Curved)?,
                    start,
                    HertzHeatMarchMotion::Curved,
                    (),
                ),
            })
        },
        gate,
    )
    .map(|(candidate, ())| candidate)
}

/// One checkpoint/refinement/physical-acceptance owner for every Hertz curve.
#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
pub(super) fn propose_hertz_heat_phases<'m, P, E, C, F, H: PartialEq>(
    bodies: &[RigidBodyState],
    properties: &[MassProperties],
    contact: RigidHertzHeatContact<'_, C>,
    thermal: &DepositionSimulation<'m>,
    initial_history: H,
    duration_s: f64,
    phase_controls: RigidHertzHeatControls,
    controls: HertzHeatMarchControls,
    temperature_guess_k: f64,
    advance_phase: impl Fn(
        &[RigidBodyState; 2],
        &[MassProperties; 2],
        RigidHertzHeatContact<'_, C>,
        &DepositionSimulation<'m>,
        &H,
        f64,
        RigidHertzHeatControls,
        &[f64],
        f64,
    ) -> Result<Stage<'m, P, F, H>, HertzHeatMarchError<E>>,
    gate: &CancelGate,
) -> Result<(HeatedRigidHertzMarchCandidate<'m, P, F>, H), HertzHeatMarchError<E>> {
    let poll = || {
        if gate.is_requested() {
            Err(HertzHeatMarchError::Cancelled)
        } else {
            Ok(())
        }
    };
    poll()?;
    let binding = contact.binding;
    let target = thermal.time_s() + duration_s;
    let t = controls.tolerance;
    if bodies.len() != properties.len()
        || bodies.len() > phase_controls.mechanical.maximum_bodies
        || binding.body_a == binding.body_b
        || binding.body_a >= bodies.len()
        || binding.body_b >= bodies.len()
        || thermal.cells().is_empty()
        || thermal.cells().len() > phase_controls.maximum_cells
        || !target.is_finite()
        || target <= thermal.time_s()
        || controls.maximum_attempts == 0
        || controls.maximum_accepted_phases < 2
        || controls.maximum_thermal_trials == 0
        || [
            duration_s,
            controls.initial_step_s,
            controls.minimum_step_s,
            controls.maximum_step_s,
            t.linear_momentum_ns,
            t.angular_momentum_nms,
            t.cell_energy_j,
            t.scaled_travel_m,
            t.dissipated_j,
            temperature_guess_k,
        ]
        .into_iter()
        .any(|x| !x.is_finite() || x <= 0.0)
        || [t.relative, controls.energy_tolerance_j]
            .into_iter()
            .any(|x| !x.is_finite() || x < 0.0)
        || controls.minimum_step_s > controls.initial_step_s
        || controls.initial_step_s > controls.maximum_step_s
    {
        return Err(HertzHeatMarchError::Invalid(
            "bounded original owners, steps and physical error scales required",
        ));
    }
    let initial = [bodies[binding.body_a], bodies[binding.body_b]];
    let props = [properties[binding.body_a], properties[binding.body_b]];
    let mut pair = contact;
    pair.binding.body_a = 0;
    pair.binding.body_b = 1;
    let (mut accepted, mut cells) = (initial, thermal.clone());
    let mut history = initial_history;
    let zero_energy = vec![0.0; thermal.cells().len()];
    let (mut ti, mut step_s) = (temperature_guess_k, controls.initial_step_s);
    let (mut attempts, mut rejected, mut trials) = (0usize, 0usize, 0usize);
    let mut phases = Vec::new();
    let mut maximum_defect = 0.0_f64;
    // Net travel (3), friction heat, external heat, signed body heat (2), drive work.
    let mut totals = [Dd::from_f64(0.0); 8];
    while cells.time_s() < target {
        poll()?;
        if attempts >= controls.maximum_attempts {
            return Err(HertzHeatMarchError::Budget("comparisons"));
        }
        if phases.len() > controls.maximum_accepted_phases - 2 {
            return Err(HertzHeatMarchError::Budget("accepted phases"));
        }
        attempts += 1;
        let start = cells.time_s();
        let remaining = target - start;
        let planned = step_s.min(remaining);
        let planned = if remaining > planned && remaining - planned < controls.minimum_step_s {
            0.5 * remaining
        } else {
            planned
        };
        let end = if planned == remaining {
            target
        } else {
            start + planned
        };
        let half = start + 0.5 * (end - start);
        if half <= start || half >= end {
            return Err(HertzHeatMarchError::Invalid(
                "comparison clock is not representable",
            ));
        }
        let mut phase = |b: &[RigidBodyState; 2],
                         h: &DepositionSimulation<'m>,
                         state: &H,
                         stop: f64,
                         guess: f64| {
            poll()?;
            let remaining = controls.maximum_thermal_trials - trials;
            if remaining == 0 {
                return Err(HertzHeatMarchError::Budget("thermal trials"));
            }
            let mut local = phase_controls;
            local.maximum_trials = local.maximum_trials.min(remaining);
            let start = h.time_s();
            let out = advance_phase(
                b,
                &props,
                pair,
                h,
                state,
                stop - start,
                local,
                &zero_energy,
                guess,
            )?;
            if out.thermal.time_s() != stop || out.receipt.trials > remaining {
                return Err(HertzHeatMarchError::Invalid(
                    "phase clock or trial budget disagrees",
                ));
            }
            trials += out.receipt.trials;
            Ok(out)
        };
        let coarse = phase(&accepted, &cells, &history, end, ti)?;
        let first = phase(&accepted, &cells, &history, half, ti)?;
        let second = phase(
            &first.bodies,
            &first.thermal,
            &first.history,
            end,
            first.receipt.heat.interfaces[0].interface_temperature_k,
        )?;
        let defect = compare(
            accepted,
            &cells,
            &coarse,
            &first,
            &second,
            contact.contact_radius_m,
            t,
            gate,
        )?;
        let dt = end - start;
        if defect.normalized > 1.0 {
            rejected += 1;
            if planned <= controls.minimum_step_s {
                return Err(HertzHeatMarchError::Unresolved { step_s: dt, defect });
            }
            step_s = (dt * (0.8 / det::sqrt(defect.normalized)).clamp(0.1, 0.5))
                .max(controls.minimum_step_s);
            continue;
        }
        poll()?;
        maximum_defect = maximum_defect.max(defect.normalized);
        for phase in [&first.receipt, &second.receipt] {
            let travel = phase.mechanical.travel(contact.contact_radius_m);
            for (sum, value) in totals.iter_mut().zip([
                travel[0],
                travel[1],
                travel[2],
                phase.heat.dissipated_j,
                phase.heat.external_heat_j,
                phase.heat.interfaces[0].body_heat_j[0],
                phase.heat.interfaces[0].body_heat_j[1],
                phase.mechanical.external_work_j(),
            ]) {
                *sum = *sum + Dd::from_f64(value);
            }
        }
        ti = second.receipt.heat.interfaces[0].interface_temperature_k;
        accepted = second.bodies;
        cells = second.thermal;
        history = second.history;
        phases
            .try_reserve(2)
            .map_err(|_| HertzHeatMarchError::Invalid("accepted phase allocation refused"))?;
        phases.extend([first.receipt, second.receipt]);
        let factor = if defect.normalized == 0.0 {
            2.0
        } else {
            (0.8 / det::sqrt(defect.normalized)).clamp(0.5, 2.0)
        };
        step_s = (dt * factor).clamp(controls.minimum_step_s, controls.maximum_step_s);
    }
    let mut kinetic = RigidKineticEnergyChange::new();
    for i in 0..2 {
        kinetic
            .add_body(props[i], initial[i], accepted[i])
            .map_err(HertzHeatMarchError::Dynamics)?;
    }
    let kinetic = kinetic
        .finish_retained()
        .map_err(HertzHeatMarchError::Dynamics)?;
    let mut stored = Dd::from_f64(0.0);
    for (index, (id, old)) in thermal.cells().iter().enumerate() {
        if index % 256 == 0 {
            poll()?;
        }
        stored = stored
            + Dd::from_f64(old.mass_kg)
                * (Dd::from_f64(cells.cells()[id].specific_h_j_kg)
                    - Dd::from_f64(old.specific_h_j_kg));
    }
    let residual_j = [
        kinetic + totals[3] - totals[7],
        stored - totals[3] - totals[4],
        kinetic + stored - totals[4] - totals[7],
    ]
    .map(Dd::to_f64);
    if residual_j
        .iter()
        .any(|x| !x.is_finite() || x.abs() > controls.energy_tolerance_j)
    {
        return Err(HertzHeatMarchError::Energy { residual_j });
    }
    let sums = totals.map(Dd::to_f64);
    let spin = sums[2] / contact.contact_radius_m;
    if sums.iter().chain([&spin]).any(|x| !x.is_finite()) {
        return Err(HertzHeatMarchError::Invalid(
            "unrepresentable accumulated physical account",
        ));
    }
    poll()?;
    Ok((
        HeatedRigidHertzMarchCandidate {
            binding,
            bodies_after: accepted,
            thermal: cells,
            phases,
            attempts,
            rejected_attempts: rejected,
            thermal_trials: trials,
            maximum_accepted_normalized_defect: maximum_defect,
            slip_increment_m: [sums[0], sums[1]],
            spin_increment_rad: spin,
            dissipated_j: sums[3],
            body_heat_j: [sums[5], sums[6]],
            external_heat_j: sums[4],
            external_mechanical_work_j: sums[7],
            kinetic_energy_change_j: kinetic.to_f64(),
            thermal_energy_change_j: stored.to_f64(),
            energy_residual_j: residual_j,
        },
        history,
    ))
}

fn norm(v: [f64; 3]) -> f64 {
    det::hypot(det::hypot(v[0], v[1]), v[2])
}
fn vector(v: Vec3) -> [f64; 3] {
    [v.x, v.y, v.z]
}
#[allow(clippy::too_many_arguments)]
fn compare<P, E, F, H: PartialEq>(
    old: [RigidBodyState; 2],
    thermal: &DepositionSimulation<'_>,
    coarse: &Stage<'_, P, F, H>,
    first: &Stage<'_, P, F, H>,
    fine: &Stage<'_, P, F, H>,
    radius: f64,
    tolerance: HertzHeatMarchTolerance,
    gate: &CancelGate,
) -> Result<HertzHeatMarchDefect, HertzHeatMarchError<E>> {
    let mut out = HertzHeatMarchDefect {
        mode_disagreement: coarse.history != fine.history,
        arrest_disagreement: coarse.receipt.mechanical.has_driven_arrest()
            != (first.receipt.mechanical.has_driven_arrest()
                || fine.receipt.mechanical.has_driven_arrest()),
        ..HertzHeatMarchDefect::default()
    };
    if out.mode_disagreement || out.arrest_disagreement {
        // Neither a velocity cutoff nor a matching final mode can erase a
        // distinct driven history. Use the existing bounded refinement path.
        out.normalized = 2.0;
    }
    let mut scale = |error: f64, absolute: f64, change: f64| {
        let denominator = absolute + tolerance.relative * change;
        if !error.is_finite() || !denominator.is_finite() || denominator <= 0.0 {
            return Err(HertzHeatMarchError::Invalid(
                "unrepresentable physical error estimate",
            ));
        }
        let normalized = error / denominator;
        if !normalized.is_finite() {
            return Err(HertzHeatMarchError::Invalid(
                "unrepresentable normalized error",
            ));
        }
        out.normalized = out.normalized.max(normalized);
        Ok(())
    };
    for i in 0..2 {
        for (getter, absolute, maximum) in [
            (
                RigidBodyState::linear_momentum_world as fn(RigidBodyState) -> Vec3,
                tolerance.linear_momentum_ns,
                &mut out.linear_momentum_ns,
            ),
            (
                RigidBodyState::angular_momentum_body as fn(RigidBodyState) -> Vec3,
                tolerance.angular_momentum_nms,
                &mut out.angular_momentum_nms,
            ),
        ] {
            let c = getter(coarse.bodies[i]);
            let f = getter(fine.bodies[i]);
            let o = getter(old[i]);
            let error = norm(vector(c.sub(f)));
            scale(
                error,
                absolute,
                norm(vector(c.sub(o))).max(norm(vector(f.sub(o)))),
            )?;
            *maximum = maximum.max(error);
        }
    }
    for (index, (id, old)) in thermal.cells().iter().enumerate() {
        if index % 256 == 0 && gate.is_requested() {
            return Err(HertzHeatMarchError::Cancelled);
        }
        let c = coarse.thermal.cells()[id].specific_h_j_kg;
        let f = fine.thermal.cells()[id].specific_h_j_kg;
        let error = old.mass_kg * (c - f).abs();
        scale(
            error,
            tolerance.cell_energy_j,
            old.mass_kg
                * (c - old.specific_h_j_kg)
                    .abs()
                    .max((f - old.specific_h_j_kg).abs()),
        )?;
        out.cell_energy_j = out.cell_energy_j.max(error);
    }
    let c = coarse.receipt.mechanical.travel(radius);
    let a = first.receipt.mechanical.travel(radius);
    let b = fine.receipt.mechanical.travel(radius);
    let f = core::array::from_fn(|i| a[i] + b[i]);
    let mut uncertainty = [0.0; 3];
    for phase in [&coarse.receipt, &first.receipt, &fine.receipt] {
        let u = phase.mechanical.uncertainty(radius);
        for i in 0..3 {
            uncertainty[i] += u[i];
        }
    }
    out.scaled_travel_m = norm(core::array::from_fn(|i| c[i] - f[i])) + uncertainty[0];
    scale(
        out.scaled_travel_m,
        tolerance.scaled_travel_m,
        norm(c).max(norm(f)),
    )?;
    let c = coarse.receipt.heat.dissipated_j;
    let f = first.receipt.heat.dissipated_j + fine.receipt.heat.dissipated_j;
    out.dissipated_j = (c - f).abs() + uncertainty[1];
    scale(out.dissipated_j, tolerance.dissipated_j, c.max(f))?;
    let c = coarse.receipt.mechanical.external_work_j();
    let f = first.receipt.mechanical.external_work_j() + fine.receipt.mechanical.external_work_j();
    out.external_mechanical_work_j = (c - f).abs() + uncertainty[2];
    scale(
        out.external_mechanical_work_j,
        tolerance.dissipated_j,
        c.abs().max(f.abs()),
    )?;
    Ok(out)
}
