//! One conservative activity solve for original slab or radial inventories.
//!
//! Geometry supplies only its physical face paths. Source admission, nonlinear
//! storage, the tridiagonal solve, physical water accounts and publication
//! limits are shared so partitioned shells do not acquire a parallel model.

use super::*;
use fs_la::tridiagonal::TridiagonalWorkspace;

/// Borrowed original material context; no surrogate material or refitted law.
pub(in crate::moisture) struct LayeredSource<'a> {
    pub card: &'a MaterialCard,
    pub transport: &'a ResolvedMaterialStatePoint,
    pub property: &'a str,
    pub sorption: &'a MonotoneMoistureSorption,
}

/// Actual accepted cell inventories, independent of the spatial reduction.
#[derive(Clone, Copy)]
pub(in crate::moisture) struct LayeredStorage<'a> {
    pub sources: &'a [LayeredSource<'a>],
    pub layer_for_cell: &'a [usize],
    pub activity: &'a [f64],
    pub original: &'a [f64],
    pub dry_masses: &'a [f64],
    pub activity_bounds: [f64; 2],
    pub time_s: f64,
    pub dry_mass_kg: f64,
    pub water_mass_kg: f64,
}

pub(in crate::moisture) struct LayeredSystem {
    pub operator: DiffusionStep,
    pub rhs: Vec<f64>,
    pub boundary: [f64; 2],
    pub surfaces: [Option<f64>; 2],
    pub interfaces: Vec<MoistureLayerInterfaceReport>,
}

impl LayeredSystem {
    fn storage(
        state: &LayeredStorage<'_>,
        activity: &[f64],
        reference: f64,
        cx: &Cx<'_>,
    ) -> Result<Self, MoistureError> {
        let n = activity.len();
        let mut result = Self {
            operator: DiffusionStep {
                diagonal: Vec::with_capacity(n),
                storage_boundary: Vec::with_capacity(n),
                off: vec![0.0; n - 1],
            },
            rhs: Vec::with_capacity(n),
            boundary: [0.0; 2],
            surfaces: [None; 2],
            interfaces: Vec::new(),
        };
        for i in 0..n {
            poll(cx)?;
            let source = &state.sources[state.layer_for_cell[i]];
            let mass = state.dry_masses[i];
            let capacity = mass
                * source
                    .sorption
                    .mean_capacity(state.activity[i], activity[i])?;
            if !positive(capacity) {
                return Err(MoistureError::Invalid(
                    "layered water capacity is unrepresentable",
                ));
            }
            let offset = state.original[i] - source.sorption.moisture(state.activity[i])?;
            result.operator.diagonal.push(capacity);
            result.operator.storage_boundary.push(capacity);
            result
                .rhs
                .push(capacity.mul_add(state.activity[i] - reference, mass * offset));
        }
        Ok(result)
    }

    fn validate(&self) -> Result<(), MoistureError> {
        if self.operator.diagonal.iter().any(|&x| !positive(x))
            || self.rhs.iter().any(|x| !x.is_finite())
            || self
                .interfaces
                .iter()
                .any(|face| !face.water_left_to_right_kg.is_finite())
        {
            return Err(MoistureError::Invalid("layered diffusion system overflow"));
        }
        Ok(())
    }
}

pub(in crate::moisture) struct LayeredCandidate {
    pub activity: Vec<f64>,
    pub moisture: Vec<f64>,
    pub report: MoistureLayeredReport,
    first_load: Option<StageLoad>,
}

/// Physical first-stage face transfers, already weighted for stage two.
/// Nonlinear storage is never replaced by an activity increment or a frozen
/// capacity. Each internal face has exactly one debit and one credit.
pub(in crate::moisture) struct StageLoad {
    water: Vec<f64>,
    boundary: [f64; 2],
    interfaces: Vec<([usize; 2], f64)>,
    bounds: [f64; 2],
}

impl StageLoad {
    fn capture(
        state: &LayeredStorage<'_>,
        endpoint: &LayeredSystem,
        activity: &[f64],
        reservoir: [Option<f64>; 2],
        boundary: [f64; 2],
        weight: f64,
        cx: &Cx<'_>,
    ) -> Result<Self, MoistureError> {
        let mut water = vec![0.0; activity.len()];
        for (i, &conductance) in endpoint.operator.off.iter().enumerate() {
            if i % 128 == 0 {
                poll(cx)?;
            }
            let flux = weight * (conductance * (activity[i + 1] - activity[i]));
            water[i] += flux;
            water[i + 1] -= flux;
        }
        let mut bounds = [f64::INFINITY, f64::NEG_INFINITY];
        for (i, &a) in state.activity.iter().enumerate() {
            if i % 128 == 0 {
                poll(cx)?;
            }
            bounds[0] = bounds[0].min(a);
            bounds[1] = bounds[1].max(a);
        }
        let boundary = boundary.map(|q| weight * q);
        for (side, cell) in [0, activity.len() - 1].into_iter().enumerate() {
            water[cell] += boundary[side];
            if endpoint.boundary[side] != 0.0 {
                let a = reservoir[side].expect("active layered boundary");
                bounds[0] = bounds[0].min(a);
                bounds[1] = bounds[1].max(a);
            }
        }
        let mut interfaces = Vec::with_capacity(endpoint.interfaces.len());
        for face in &endpoint.interfaces {
            poll(cx)?;
            let transfer = weight * face.water_left_to_right_kg;
            if !transfer.is_finite() {
                return Err(MoistureError::Invalid(
                    "layered SDIRK interface load overflow",
                ));
            }
            interfaces.push((face.layers, transfer));
        }
        if water.iter().chain(&boundary).any(|q| !q.is_finite()) {
            return Err(MoistureError::Invalid("layered SDIRK face load overflow"));
        }
        Ok(Self {
            water,
            boundary,
            interfaces,
            bounds,
        })
    }

    fn add_to(&self, system: &mut LayeredSystem, cx: &Cx<'_>) -> Result<(), MoistureError> {
        debug_assert_eq!(self.water.len(), system.rhs.len());
        for (i, (rhs, &water)) in system.rhs.iter_mut().zip(&self.water).enumerate() {
            if i % 128 == 0 {
                poll(cx)?;
            }
            *rhs += water;
        }
        if self.interfaces.len() != system.interfaces.len() {
            return Err(MoistureError::Invalid(
                "layered SDIRK interface topology changed",
            ));
        }
        for (face, &(layers, water)) in system.interfaces.iter_mut().zip(&self.interfaces) {
            poll(cx)?;
            if face.layers != layers {
                return Err(MoistureError::Invalid(
                    "layered SDIRK interface order changed",
                ));
            }
            face.water_left_to_right_kg += water;
        }
        Ok(())
    }

    fn check_range(
        &self,
        activity: &[f64],
        reservoir: [Option<f64>; 2],
        active: [f64; 2],
        cx: &Cx<'_>,
    ) -> Result<(), MoistureError> {
        let mut bounds = self.bounds;
        for side in 0..2 {
            if active[side] != 0.0 {
                let a = reservoir[side].expect("active layered boundary");
                bounds[0] = bounds[0].min(a);
                bounds[1] = bounds[1].max(a);
            }
        }
        let roundoff = 128.0 * f64::EPSILON;
        for (cell, &a) in activity.iter().enumerate() {
            if cell % 128 == 0 {
                poll(cx)?;
            }
            if !a.is_finite()
                || !(0.0..=1.0).contains(&a)
                || a < bounds[0] - roundoff
                || a > bounds[1] + roundoff
            {
                return Err(MoistureError::ActivityStageRange {
                    cell,
                    activity: a,
                    bounds,
                });
            }
        }
        Ok(())
    }
}

/// Native geometry adds faces to the same implicit storage operator. A caller
/// may resolve a physical deformation metric at each proposed activity before
/// assembling paths; that endpoint remains subject to the original residual.
#[allow(clippy::too_many_arguments)]
pub(in crate::moisture) fn advance(
    state: LayeredStorage<'_>,
    dt: f64,
    reservoir: [Option<f64>; 2],
    step: MoistureStepConfig,
    controls: MoistureNonlinearConfig,
    paths: impl FnMut(
        &[KirchhoffMoistureDiffusivity],
        &[f64],
        f64,
        &mut LayeredSystem,
        &Cx<'_>,
    ) -> Result<(), MoistureError>,
    cx: &Cx<'_>,
) -> Result<LayeredCandidate, MoistureError> {
    advance_stage(state, dt, reservoir, step, controls, None, None, paths, cx)
}

/// Alexander's two-stage L-stable method on physical original water storage.
/// The accepted first face rates contribute (1-gamma)*dt; the final implicit
/// rates contribute gamma*dt. Both stages use the original nonlinear kernel.
#[allow(clippy::too_many_arguments)]
pub(in crate::moisture) fn sdirk2(
    state: LayeredStorage<'_>,
    dt: f64,
    reservoir: [Option<f64>; 2],
    step: MoistureStepConfig,
    controls: MoistureNonlinearConfig,
    mut paths: impl FnMut(
        f64,
        &[KirchhoffMoistureDiffusivity],
        &[f64],
        f64,
        &mut LayeredSystem,
        &Cx<'_>,
    ) -> Result<(), MoistureError>,
    cx: &Cx<'_>,
) -> Result<(LayeredCandidate, MoistureLayeredReport), MoistureError> {
    poll(cx)?;
    let gamma = crate::moisture::sdirk::GAMMA;
    let implicit_dt = sdirk_interval(state.time_s, dt)?;
    let mut first = advance_stage(
        state,
        implicit_dt,
        reservoir,
        step,
        controls,
        None,
        Some((1.0 - gamma) / gamma),
        |profiles, activity, reference, system, cx| {
            paths(implicit_dt, profiles, activity, reference, system, cx)
        },
        cx,
    )?;
    let known = first
        .first_load
        .take()
        .expect("requested first-stage face load");
    let next = advance_stage(
        state,
        dt,
        reservoir,
        step,
        controls,
        Some(&known),
        None,
        |profiles, activity, reference, system, cx| {
            paths(implicit_dt, profiles, activity, reference, system, cx)
        },
        cx,
    )?;
    poll(cx)?;
    Ok((next, first.report))
}

/// Shared representability gate for both fixed and finite-reservoir stages.
pub(in crate::moisture) fn sdirk_interval(time_s: f64, dt: f64) -> Result<f64, MoistureError> {
    let implicit_dt = crate::moisture::sdirk::GAMMA * dt;
    let first_time = time_s + implicit_dt;
    let end = time_s + dt;
    if !positive(dt)
        || !positive(implicit_dt)
        || !end.is_finite()
        || first_time <= time_s
        || first_time >= end
    {
        return Err(MoistureError::Invalid(
            "layered SDIRK needs two representable forward stage times",
        ));
    }
    Ok(implicit_dt)
}

/// Reassemble only the already accepted first-stage physical faces, once after
/// the finite reservoir root has converged. No extra nonlinear solve is run,
/// and no activity/storage difference is substituted for a physical face rate.
#[allow(clippy::too_many_arguments)]
pub(in crate::moisture) fn capture_sdirk_load(
    state: LayeredStorage<'_>,
    dt: f64,
    first_activity: &[f64],
    first_time_s: f64,
    reservoir: [Option<f64>; 2],
    paths: impl FnOnce(
        f64,
        &[KirchhoffMoistureDiffusivity],
        &[f64],
        f64,
        &mut LayeredSystem,
        &Cx<'_>,
    ) -> Result<(), MoistureError>,
    cx: &Cx<'_>,
) -> Result<StageLoad, MoistureError> {
    poll(cx)?;
    let implicit_dt = sdirk_interval(state.time_s, dt)?;
    if first_activity.len() != state.activity.len()
        || first_time_s.to_bits() != (state.time_s + implicit_dt).to_bits()
    {
        return Err(MoistureError::Invalid(
            "layered SDIRK load needs the admitted first-stage layout and clock",
        ));
    }
    let (reference, _, profiles) = prepare_profiles(&state, reservoir, None, cx)?;
    let mut system = LayeredSystem::storage(&state, first_activity, reference, cx)?;
    paths(
        implicit_dt,
        &profiles,
        first_activity,
        reference,
        &mut system,
        cx,
    )?;
    system.validate()?;
    let boundary = std::array::from_fn(|side| {
        let cell = if side == 0 {
            0
        } else {
            first_activity.len() - 1
        };
        reservoir[side].map_or(0.0, |a| system.boundary[side] * (a - first_activity[cell]))
    });
    let gamma = crate::moisture::sdirk::GAMMA;
    StageLoad::capture(
        &state,
        &system,
        first_activity,
        reservoir,
        boundary,
        (1.0 - gamma) / gamma,
        cx,
    )
}

/// Final native stage from ORIGINAL storage with the retained physical load.
/// The finite-reservoir owner is free to supply a different endpoint activity.
#[allow(clippy::too_many_arguments)]
pub(in crate::moisture) fn finish_sdirk(
    state: LayeredStorage<'_>,
    dt: f64,
    reservoir: [Option<f64>; 2],
    step: MoistureStepConfig,
    controls: MoistureNonlinearConfig,
    known: &StageLoad,
    mut paths: impl FnMut(
        f64,
        &[KirchhoffMoistureDiffusivity],
        &[f64],
        f64,
        &mut LayeredSystem,
        &Cx<'_>,
    ) -> Result<(), MoistureError>,
    cx: &Cx<'_>,
) -> Result<LayeredCandidate, MoistureError> {
    let implicit_dt = sdirk_interval(state.time_s, dt)?;
    advance_stage(
        state,
        dt,
        reservoir,
        step,
        controls,
        Some(known),
        None,
        |profiles, activity, reference, system, cx| {
            paths(implicit_dt, profiles, activity, reference, system, cx)
        },
        cx,
    )
}

fn prepare_profiles(
    state: &LayeredStorage<'_>,
    reservoir: [Option<f64>; 2],
    known: Option<&StageLoad>,
    cx: &Cx<'_>,
) -> Result<(f64, [f64; 2], Vec<KirchhoffMoistureDiffusivity>), MoistureError> {
    if reservoir
        .into_iter()
        .flatten()
        .any(|a| !a.is_finite() || a < state.activity_bounds[0] || a > state.activity_bounds[1])
    {
        return Err(MoistureError::Invalid(
            "boundary activity lies outside the layers' common source support",
        ));
    }
    let mut span = [f64::INFINITY, f64::NEG_INFINITY];
    for (i, a) in state
        .activity
        .iter()
        .copied()
        .chain(reservoir.into_iter().flatten())
        .enumerate()
    {
        if i % 128 == 0 {
            poll(cx)?;
        }
        span[0] = span[0].min(a);
        span[1] = span[1].max(a);
    }
    // Keep the existing old/current algebraic reference even when retained
    // first-stage source support needs a wider preparation interval.
    let reference = span[0];
    if let Some(known) = known {
        // Retain first-stage source support even if a future finite-reservoir
        // second stage supplies a different boundary. This is source support,
        // not permission for a physical activity overshoot.
        span[0] = span[0].min(known.bounds[0]);
        span[1] = span[1].max(known.bounds[1]);
    }
    let mut profiles = Vec::with_capacity(state.sources.len());
    for source in state.sources {
        poll(cx)?;
        profiles.push(KirchhoffMoistureDiffusivity::from_binding(
            source.card,
            source.transport,
            source.property,
            [
                source.sorption.moisture(span[0])?,
                source.sorption.moisture(span[1])?,
            ],
        )?);
    }
    Ok((reference, span, profiles))
}

#[allow(clippy::too_many_arguments)]
fn advance_stage(
    state: LayeredStorage<'_>,
    dt: f64,
    reservoir: [Option<f64>; 2],
    step: MoistureStepConfig,
    controls: MoistureNonlinearConfig,
    known: Option<&StageLoad>,
    capture_weight: Option<f64>,
    mut paths: impl FnMut(
        &[KirchhoffMoistureDiffusivity],
        &[f64],
        f64,
        &mut LayeredSystem,
        &Cx<'_>,
    ) -> Result<(), MoistureError>,
    cx: &Cx<'_>,
) -> Result<LayeredCandidate, MoistureError> {
    poll(cx)?;
    controls.validate()?;
    let end = state.time_s + dt;
    if !positive(dt)
        || !end.is_finite()
        || end <= state.time_s
        || !positive(step.water_balance_tolerance_kg)
        || !positive(step.maximum_content_increment)
    {
        return Err(MoistureError::Invalid(
            "representable forward layered interval and positive limits required",
        ));
    }
    debug_assert_eq!(state.activity.len(), state.original.len());
    debug_assert_eq!(state.activity.len(), state.dry_masses.len());
    debug_assert_eq!(state.activity.len(), state.layer_for_cell.len());
    let (reference, span, profiles) = prepare_profiles(&state, reservoir, known, cx)?;
    // Exact uniform activity has a zero RHS even in an arbitrarily stiff solve.
    let mut assemble = |activity: &[f64]| -> Result<LayeredSystem, MoistureError> {
        let mut system = LayeredSystem::storage(&state, activity, reference, cx)?;
        paths(&profiles, activity, reference, &mut system, cx)?;
        if let Some(known) = known {
            known.add_to(&mut system, cx)?;
        }
        system.validate()?;
        Ok(system)
    };
    let mut iterate = state.activity.to_vec();
    let mut system = assemble(&iterate)?;
    let mut workspace =
        TridiagonalWorkspace::new(state.activity.len()).map_err(MoistureError::Tridiagonal)?;
    let mut solves = 0usize;
    for iteration in 1..=controls.maximum_iterations {
        poll(cx)?;
        let solution =
            system
                .operator
                .solve_tridiagonal(&system.rhs, step.linear, &mut workspace, cx)?;
        solves = solves
            .checked_add(solution.report.iters)
            .ok_or(MoistureError::Invalid("layered solve work count overflow"))?;
        if !solution.report.converged_euclidean() {
            return Err(MoistureError::NotConverged(solution.report));
        }
        let candidate: Vec<f64> = solution.x.iter().map(|delta| reference + delta).collect();
        if let Some(known) = known {
            known.check_range(&candidate, reservoir, system.boundary, cx)?;
        }
        let mut moisture = Vec::with_capacity(candidate.len());
        let mut transport = Vec::with_capacity(candidate.len());
        let mut changes = Vec::with_capacity(candidate.len());
        let roundoff = 128.0 * f64::EPSILON;
        for (i, &a) in candidate.iter().enumerate() {
            poll(cx)?;
            if !a.is_finite() || a < span[0] - roundoff || a > span[1] + roundoff {
                return Err(MoistureError::Invalid(
                    "layered candidate violates the activity maximum principle",
                ));
            }
            let source = &state.sources[state.layer_for_cell[i]];
            let m = source.sorption.moisture(a)?;
            let point = source
                .transport
                .query()
                .with("moisture_content", m)
                .map_err(|_| {
                    MoistureError::Invalid("unrepresentable dry-basis moisture coordinate")
                })?;
            let endpoint = source.transport.at_query_point(source.card, &point)?;
            changes.push(
                state.dry_masses[i]
                    * (source.sorption.moisture_change(state.activity[i], a)?
                        + (source.sorption.moisture(state.activity[i])? - state.original[i])),
            );
            moisture.push(m);
            transport.push(endpoint);
        }
        let water_change = reservoir::water_sum(changes);
        let endpoint = assemble(&candidate)?;
        if let Some(known) = known {
            known.check_range(&candidate, reservoir, endpoint.boundary, cx)?;
        }
        let (_, relative) = endpoint.operator.defect(&solution.x, &endpoint.rhs, cx)?;
        let boundary = std::array::from_fn(|side| {
            let cell = if side == 0 { 0 } else { candidate.len() - 1 };
            let current =
                reservoir[side].map_or(0.0, |a| endpoint.boundary[side] * (a - candidate[cell]));
            current + known.map_or(0.0, |load| load.boundary[side])
        });
        let water_after = moisture
            .iter()
            .zip(state.dry_masses)
            .map(|(m, mass)| m * mass)
            .sum::<f64>();
        let received = boundary[0] + boundary[1];
        let balance = water_change - received;
        let stored_balance = water_after - state.water_mass_kg - received;
        let nonlinear = MoistureNonlinearReport {
            iterations: iteration,
            tridiagonal_solves: solves,
            relative_residual: relative,
            water_balance_residual_kg: balance,
        };
        if relative < controls.relative_tolerance
            && balance.is_finite()
            && balance.abs() <= step.water_balance_tolerance_kg
            && stored_balance.is_finite()
            && stored_balance.abs() <= step.water_balance_tolerance_kg
        {
            for (cell, (&m, &old)) in moisture.iter().zip(state.original).enumerate() {
                let increment = (m - old).abs();
                if increment > step.maximum_content_increment {
                    return Err(MoistureError::ContentIncrement {
                        cell,
                        increment,
                        maximum: step.maximum_content_increment,
                    });
                }
            }
            if !positive(state.dry_mass_kg + water_after) {
                return Err(MoistureError::Invalid(
                    "layered endpoint wet inventory is unrepresentable",
                ));
            }
            let first_load = capture_weight
                .map(|weight| {
                    StageLoad::capture(
                        &state, &endpoint, &candidate, reservoir, boundary, weight, cx,
                    )
                })
                .transpose()?;
            let report = MoistureLayeredReport {
                interval_s: [state.time_s, end],
                boundary_water_kg: boundary,
                surface_activities: endpoint.surfaces,
                interfaces: endpoint.interfaces,
                stored_water_kg: [state.water_mass_kg, water_after],
                water_balance_residual_kg: balance,
                diffusivity_sources: profiles
                    .iter()
                    .map(KirchhoffMoistureDiffusivity::identity)
                    .collect(),
                sorption_sources: state
                    .sources
                    .iter()
                    .map(|s| s.sorption.identity())
                    .collect(),
                transport_end: transport,
                linear: solution.report,
                nonlinear,
            };
            poll(cx)?;
            return Ok(LayeredCandidate {
                activity: candidate,
                moisture,
                report,
                first_load,
            });
        }
        if iteration == controls.maximum_iterations {
            return Err(MoistureError::NonlinearNotConverged(nonlinear));
        }
        if controls.relaxation == 1.0 {
            iterate = candidate;
            system = endpoint;
        } else {
            for (a, candidate) in iterate.iter_mut().zip(candidate) {
                *a = controls.relaxation.mul_add(candidate - *a, *a);
            }
            system = assemble(&iterate)?;
        }
    }
    unreachable!("validated positive layered iteration budget")
}
