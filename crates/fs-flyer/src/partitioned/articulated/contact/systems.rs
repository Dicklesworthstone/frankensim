//! Simultaneous sustained contact between independently assembled fluid-bearing trees.
//!
//! One common midpoint iterate supplies every dry/fluid mass, moving-fluid bias,
//! implicit joint impedance and material-point map. Only incident tree blocks
//! enter a contact row; no dense global mass or pairwise sequential solve is used.

use super::*;

/// One complete independently assembled system, in stable caller index order.
#[derive(Clone, Copy, Debug)]
pub struct ContactSystem<'a> {
    pub model: &'a ArticulatedModel,
    pub state: &'a ArticulatedState,
    pub impedance: &'a [JointImpedance],
    pub memory_inputs: &'a [f64],
    /// Per-tree physical tolerances. The common iteration budget is the minimum
    /// of all supplied `max_iterations`, so no tree exceeds its declared budget.
    pub smooth: ArticulatedStepConfig<'a>,
}

/// Every accepted numerical configuration path, for mandatory whole-path admission.
/// Witness identities and frozen world contact frames must remain applicable;
/// endpoint or finitely sampled gaps alone do not establish feature coverage.
pub struct ContactSystemsPath<'a> {
    pub systems: &'a [ContactSystem<'a>],
    pub midpoints: &'a [ArticulatedState],
    pub endpoints: &'a [&'a ArticulatedState],
    /// Ascending contact identity order; prescribed planes are at window start.
    pub contacts: &'a [ArticulatedFrictionalContact],
    pub start_time_s: f64,
    pub duration_s: f64,
}
impl ContactSystemsPath<'_> {
    /// Sample the original Lie configuration interpolant of every tree.
    /// Twists differentiate that interpolant; this is not a force reintegration
    /// or a certified enclosure of the true dynamical flow.
    pub fn sample(
        &self,
        fraction: f64,
        mut is_cancelled: impl FnMut() -> bool,
    ) -> Result<Vec<Kinematics>, Refusal> {
        self.systems
            .iter()
            .enumerate()
            .map(|(index, system)| {
                ContactPath {
                    model: system.model,
                    initial: system.state,
                    midpoint: &self.midpoints[index],
                    endpoint: self.endpoints[index],
                    contacts: self.contacts,
                    start_time_s: self.start_time_s,
                    duration_s: self.duration_s,
                }
                .sample(fraction, &mut is_cancelled)
            })
            .collect()
    }
}

/// One accepted simultaneous step; reaction receipts use stable contact order.
#[derive(Clone, Debug)]
pub struct ContactSystemsStep<P> {
    /// Accepted states in original system order, evaluated at a common stage clock.
    pub steps: Vec<ArticulatedStep>,
    /// Original auxiliary half-step solve, aligned with the velocity increment
    /// actually accepted by the smooth owner. Its `dissipated_energy` describes
    /// that auxiliary projection; it is NOT physical friction work or heat.
    pub reaction_solve: FrictionalContactResult,
    pub contact_ids: Vec<u64>,
    /// Actual native material-point velocities, in normal/tangent/tangent order.
    pub midpoint_relative_velocities: Vec<[f64; 3]>,
    pub endpoint_relative_velocities: Vec<[f64; 3]>,
    /// Actual material gaps, in the same stable identity order.
    pub initial_gaps: Vec<f64>,
    pub endpoint_gaps: Vec<f64>,
    /// Original normal speeds are retained, including arithmetic-sized negative
    /// values. No input state or gap is projected to erase them.
    pub initial_normal_speeds: Vec<f64>,
    /// Per-contact arithmetic allowance at the accepted midpoint, m/s. It scales
    /// uncancelled original material/support motion and the ORIGINAL free velocity
    /// increment before contact. User solver tolerances do not enter this value.
    pub normal_speed_roundoff_mps: Vec<f64>,
    /// Residuals reconstructed from accepted midpoint velocities and the aligned
    /// reaction impulses, independently of the auxiliary solver's receipt.
    pub midpoint_law_residuals: FrictionalContactResiduals,
    /// Actual midpoint velocity minus its aligned auxiliary prediction, m/s.
    pub midpoint_velocity_residual_mps: f64,
    /// Reaction impulse paired with actual midpoint material velocities, J.
    pub reaction_work_j: f64,
    /// Normal reaction paired with actual relative midpoint normal speed, J.
    pub normal_work_j: f64,
    /// SIGNED tangential impulse work at the actual midpoint, J; dissipative
    /// within the declared residual budget. No heat/deformation split is inferred.
    pub friction_work_j: f64,
    /// Work delivered by prescribed translating supports, J.
    pub support_work_j: f64,
    /// External contact impulse from prescribed supports; internal pairs cancel.
    pub external_linear_impulse_world: Vec3,
    pub external_angular_impulse_world: Vec3,
    pub path_admission: P,
}

/// Advance fixed-feature sustained contact with complete independent-tree inertia.
///
/// The pure load callback sees ALL synchronized midpoint states and motion maps
/// before any reaction is applied. Return one load per tree in input order; this
/// permits cross-tree force coupling while retaining separate assembled inertia.
/// Hydrodynamic mass coupling between distinct trees is not inferred.
///
/// Supports one to 128 trees, at most 4096 total links, 256 links/eight joints per
/// tree, and 21 contacts. Each contact references these original system indices.
/// Strict actual gaps, arithmetic-only incoming-speed admission, unresolved-mode
/// refusals and actual reaction work retain the single-tree contract. An incoming
/// impact requires the physical impact owner; no configuration is projected.
///
/// The shared midpoint owner publishes only a fully converged simultaneous
/// candidate, with the previous response aligned to that candidate. The required
/// callback admits the complete collection of configuration paths and omitted
/// feature coverage before any result is returned. General nonlinear smooth
/// momentum/energy accuracy remains refinement-based, not exact conservation.
#[allow(clippy::too_many_arguments)]
pub fn step_with_contacts<F, P>(
    systems: &[ContactSystem<'_>],
    t_s: f64,
    dt_s: f64,
    contacts: &[ArticulatedFrictionalContact],
    config: FrictionalContactConfig,
    evaluate: F,
    admit_path: impl FnOnce(&ContactSystemsPath<'_>) -> Result<P, Refusal>,
    mut is_cancelled: impl FnMut() -> bool,
) -> Result<ContactSystemsStep<P>, Refusal>
where
    F: Fn(f64, &[ArticulatedStage<'_>]) -> Result<Vec<ArticulatedLoads>, Refusal>,
{
    poll(&mut is_cancelled)?;
    if systems.is_empty()
        || systems.len() > 128
        || systems.iter().any(|s| {
            s.model.link_count() > super::super::MAX_ARTICULATED_LINKS
                || s.model.dof_count() > super::super::MAX_CONTROL_COORDS
        })
        || systems
            .iter()
            .try_fold(0usize, |n, s| n.checked_add(s.model.link_count()))
            .is_none_or(|n| n > 4096)
        || contacts.len() > 21
        || contacts.len() > config.maximum_contacts
        || !t_s.is_finite()
        || !dt_s.is_finite()
        || dt_s <= 0.0
    {
        return Err(invalid(
            "invalid bounded systems, contacts or step duration",
        ));
    }
    let mut ordered = contacts.to_vec();
    ordered.sort_unstable_by_key(|c| c.id);
    if ordered.windows(2).any(|pair| pair[0].id == pair[1].id) {
        return Err(invalid("duplicate contact identity"));
    }
    for c in &ordered {
        if c.point.system_index >= systems.len()
            || matches!(c.other, ContactCounterpart::Body(p) if p.system_index >= systems.len())
            || !c.activation_tolerance.is_finite()
            || c.activation_tolerance < 0.0
            || !c.friction_coefficient.is_finite()
            || c.friction_coefficient < 0.0
            || !matches!(
                c.normal_law,
                NormalVelocityLaw::Target { speed: 0.0 }
                    | NormalVelocityLaw::NewtonRestitution { coefficient: 0.0 }
            )
        {
            return Err(invalid(
                "sustained contacts require admitted trees, finite friction and zero normal target/restitution",
            ));
        }
    }
    let mut initial_native = Vec::with_capacity(systems.len());
    let mut initial_velocities = Vec::with_capacity(systems.len());
    for system in systems {
        poll(&mut is_cancelled)?;
        admit_state(system.model, &system.state.motion)?;
        initial_native.push(kinematics(system.model, system.state)?);
        initial_velocities.push(velocity(system.state));
    }
    let native_refs: Vec<_> = initial_native.iter().collect();
    let before = frames(&native_refs, &ordered, 0.0, &mut is_cancelled)?;
    let initial_gaps = nonnegative_gaps(&before)?;
    let inputs: Vec<_> = systems
        .iter()
        .map(|s| StepSystem {
            model: s.model,
            impedance: s.impedance,
            state: s.state,
            memory_inputs: s.memory_inputs,
            config: s.smooth,
        })
        .collect();
    let mut previous = None;
    let mut accepted = None;
    let steps = step_systems_with_velocity_response(
        &inputs,
        t_s,
        dt_s,
        evaluate,
        |contexts, increments, cancel| {
            let normal_speed_roundoff =
                normal_speed_admission(contexts, increments, &ordered, &before, &initial_native)?;
            let native: Vec<_> = contexts
                .iter()
                .map(|c| c.kinematics.native_kinematics())
                .collect();
            let geometry = frames(&native, &ordered, 0.5 * dt_s, cancel)?;
            let next = solve_midpoint(
                contexts,
                increments,
                &initial_velocities,
                &ordered,
                &initial_gaps,
                &geometry,
                config,
                cancel,
            )?;
            // All trees publish the evaluated old-delta candidate together.
            // Its impulse is the PREVIOUS response, never the unapplied next one.
            let aligned = previous.replace(next);
            if let Some(solve) = aligned {
                let (residuals, velocity_residual, ready) =
                    actual_law(&solve, &geometry, &ordered, config)?;
                accepted = ready.then(|| AcceptedReaction {
                    solve,
                    midpoints: contexts.iter().map(|c| c.midpoint.clone()).collect(),
                    geometry,
                    residuals,
                    velocity_residual,
                    normal_speed_roundoff,
                });
                Ok(ready)
            } else {
                Ok(false)
            }
        },
        &mut is_cancelled,
    )?;
    let accepted = accepted.ok_or_else(|| unresolved("no aligned accepted reaction"))?;
    nonnegative_gaps(&accepted.geometry)?;
    let endpoint_native: Vec<_> = systems
        .iter()
        .zip(&steps)
        .map(|(s, step)| kinematics(s.model, &step.state))
        .collect::<Result<_, _>>()?;
    let endpoint_refs: Vec<_> = endpoint_native.iter().collect();
    let after = frames(&endpoint_refs, &ordered, dt_s, &mut is_cancelled)?;
    let endpoint_gaps = nonnegative_gaps(&after)?;
    for &i in &accepted.solve.active_indices {
        let p = accepted.solve.impulses[i];
        let first = before[i].rows.map(|r| r.relative_normal_speed);
        let last = after[i].rows.map(|r| r.relative_normal_speed);
        let c = ordered[i];
        if p[0] > config.impulse_tolerance
            && (accepted.geometry[i].rows[0].gap > c.activation_tolerance
                || endpoint_gaps[i] > c.activation_tolerance
                || first[0] > config.velocity_tolerance
                || last[0] < -config.velocity_tolerance)
        {
            return Err(unresolved(
                "resolve separation or recontact before applying a sustained normal reaction",
            ));
        }
        let first_slip = finite(first[1].hypot(first[2]))?;
        if p[1].hypot(p[2]) > config.impulse_tolerance
            && first_slip > config.velocity_tolerance
            && finite(first[1] * last[1] + first[2] * last[2])?
                < -config.velocity_tolerance * first_slip
        {
            return Err(unresolved(
                "the step crosses an unresolved friction stopping or turning event",
            ));
        }
    }
    if matches!(
        accepted.solve.determinacy,
        FrictionalContactDeterminacy::UniqueVelocity(_)
    ) {
        let mut common = None;
        for &i in &accepted.solve.active_indices {
            if let ContactCounterpart::PrescribedPlane { velocity_world, .. } = ordered[i].other {
                if common.is_some_and(|v| v != velocity_world) {
                    return Err(unresolved(
                        "redundant reactions do not establish work for independently moving supports",
                    ));
                }
                common = Some(velocity_world);
            }
        }
    }
    let midpoint_native: Vec<_> = systems
        .iter()
        .zip(&accepted.midpoints)
        .map(|(s, state)| kinematics(s.model, state))
        .collect::<Result<_, _>>()?;
    let midpoint_refs: Vec<_> = midpoint_native.iter().collect();
    let (external_linear_impulse_world, external_angular_impulse_world) =
        spatial_impulses(&accepted.geometry, &midpoint_refs, &accepted.solve.impulses)?;
    let (mut reaction_work, mut normal_work, mut friction_work, mut support_work) =
        (0.0, 0.0, 0.0, 0.0);
    for (g, p) in accepted.geometry.iter().zip(&accepted.solve.impulses) {
        normal_work = finite(normal_work + p[0] * g.rows[0].relative_normal_speed)?;
        for a in 0..3 {
            let relative = finite(p[a] * g.rows[a].relative_normal_speed)?;
            let support = finite(p[a] * g.rows[a].plane_speed)?;
            reaction_work = finite(reaction_work + relative + support)?;
            support_work = finite(support_work + support)?;
            if a > 0 {
                friction_work = finite(friction_work + relative)?;
            }
        }
    }
    if friction_work > config.complementarity_tolerance
        || normal_work.abs() > config.complementarity_tolerance * (ordered.len().max(1) as f64)
    {
        return Err(unresolved(
            "accepted contact work exceeds the declared physical residual budget",
        ));
    }
    poll(&mut is_cancelled)?;
    let endpoint_states: Vec<_> = steps.iter().map(|step| &step.state).collect();
    let path_admission = admit_path(&ContactSystemsPath {
        systems,
        midpoints: &accepted.midpoints,
        endpoints: &endpoint_states,
        contacts: &ordered,
        start_time_s: t_s,
        duration_s: dt_s,
    })?;
    poll(&mut is_cancelled)?;
    Ok(ContactSystemsStep {
        steps,
        contact_ids: ordered.iter().map(|c| c.id).collect(),
        midpoint_relative_velocities: accepted
            .geometry
            .iter()
            .map(|g| g.rows.map(|r| r.relative_normal_speed))
            .collect(),
        endpoint_relative_velocities: after
            .iter()
            .map(|g| g.rows.map(|r| r.relative_normal_speed))
            .collect(),
        initial_gaps,
        endpoint_gaps,
        initial_normal_speeds: before
            .iter()
            .map(|g| g.rows[0].relative_normal_speed)
            .collect(),
        normal_speed_roundoff_mps: accepted.normal_speed_roundoff,
        midpoint_law_residuals: accepted.residuals,
        midpoint_velocity_residual_mps: accepted.velocity_residual,
        reaction_work_j: reaction_work,
        normal_work_j: normal_work,
        friction_work_j: friction_work,
        support_work_j: support_work,
        external_linear_impulse_world,
        external_angular_impulse_world,
        reaction_solve: accepted.solve,
        path_admission,
    })
}
