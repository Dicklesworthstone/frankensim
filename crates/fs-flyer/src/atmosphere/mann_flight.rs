//! One coherent Mann history at the actual body/shaft/wake force stages.
//!
//! The original fs-atmo fixed-step owner advances by half a flight step.
//! Start and midpoint loads read their own accepted OU quadratures and
//! advection phase; every wing, hub, blade and drag point shares that history.
//! The endpoint becomes public only after the native mechanical solve succeeds.
//! No velocity average is substituted for a nonlinear force average.
//!
//! World coordinates are [flight-line forward, right, down]. Atmospheric
//! coordinates are [-world_x, world_y, ground_down - world_z]; the prescribed
//! log-law headwind points along negative world x. Rotation of the aircraft is
//! already owned by the native propulsion sampler. This neutral, finite-band
//! source adds neither resolved eddy dynamics nor unsteady airfoil memory.

use super::wind_refusal;
use crate::Refusal;
use crate::spine::rotors::propulsion::PropulsionInputs;
use crate::spine::rotors::propulsion::dynamic::first_harmonic::{
    HarmonicPropulsion, HarmonicState, HarmonicStep,
};
use crate::spine::rotors::propulsion::dynamic::{DynamicPropulsion, DynamicState, DynamicStep};
use fs_atmo::mann_field::fixed_step::{FixedStepMannCheckpoint, FixedStepMannField};
use fs_atmo::mann_field::wall::MannPlaneWall;
use fs_atmo::{FlatSiteLogLaw, FlowSample};

/// Explicit applicability of the underlying complete Mann spectrum.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum MannFlightBoundary {
    /// Homogeneous outer flow, admitted only at/above this positive height.
    /// The authored bound does not establish that real ground effects vanish.
    HomogeneousAbove { minimum_height_m: f64 },
    /// Harmonic no-penetration correction at the declared ground plane.
    /// This is a slip boundary, without a viscous layer or wall-generated eddies.
    GroundPlane,
}

/// Accepted atmospheric state; scratch stages are rebuilt on restoration.
/// Store this alongside the returned native mechanical/wake state.
#[derive(Clone, Debug, PartialEq)]
pub struct MannFlightWindCheckpoint {
    pub version: u32,
    pub field: FixedStepMannCheckpoint,
    pub mean: FlatSiteLogLaw,
    pub boundary: MannFlightBoundary,
    pub ground_down_m: f64,
}

/// Prepared start/midpoint/endpoint storage for the original force owners.
/// Successful wind advances and queries reuse numerical storage. The native
/// aerodynamic solves retain their existing allocation and work behavior.
#[derive(Clone, Debug)]
pub struct MannFlightWind {
    accepted: FixedStepMannField,
    midpoint: FixedStepMannField,
    endpoint: FixedStepMannField,
    mean: FlatSiteLogLaw,
    boundary: MannFlightBoundary,
    ground_down_m: f64,
}

fn invalid(code: &'static str, message: &str) -> Refusal {
    Refusal {
        code,
        message: message.into(),
        ranked_repairs: vec![
            "use the same admitted half-step wind clock and keep all material points in its domain"
                .into(),
        ],
    }
}

fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), Refusal> {
    if cancelled() {
        Err(invalid(
            "mann-flight-cancelled",
            "flight wind step cancelled",
        ))
    } else {
        Ok(())
    }
}

impl MannFlightWind {
    /// Use an explicitly selected fixed-step realization. Its dt is HALF the
    /// native mechanical dt. Its phase origin, amplitudes, modal lifetimes and
    /// future innovation namespace are retained; no legacy 120 Hz path changes.
    /// The budget conservatively counts all three owners' temporal capacities,
    /// including their shared immutable temporal arrays three times. It excludes
    /// the already-owned immutable spectral geometry and fixed-size metadata.
    pub fn prepare(
        field: FixedStepMannField,
        mean: FlatSiteLogLaw,
        boundary: MannFlightBoundary,
        ground_down_m: f64,
        maximum_workspace_bytes: usize,
    ) -> Result<Self, Refusal> {
        mean.admit().map_err(wind_refusal)?;
        let minimum = match boundary {
            MannFlightBoundary::HomogeneousAbove { minimum_height_m } => minimum_height_m,
            MannFlightBoundary::GroundPlane => 0.0,
        };
        if !ground_down_m.is_finite()
            || !minimum.is_finite()
            || minimum < 0.0
            || (matches!(boundary, MannFlightBoundary::HomogeneousAbove { .. }) && minimum == 0.0)
            || !(2.0 * field.dt_s()).is_finite()
        {
            return Err(invalid(
                "mann-flight-invalid",
                "invalid wind domain or flight timestep",
            ));
        }
        let required = field
            .workspace_bytes()
            .checked_mul(3)
            .ok_or_else(|| invalid("mann-flight-budget", "temporal workspace size overflows"))?;
        if required > maximum_workspace_bytes {
            return Err(invalid(
                "mann-flight-budget",
                "temporal workspace exceeds its byte budget",
            ));
        }
        let result = Self {
            midpoint: field.clone(),
            endpoint: field.clone(),
            accepted: field,
            mean,
            boundary,
            ground_down_m,
        };
        if result.workspace_bytes() > maximum_workspace_bytes {
            return Err(invalid(
                "mann-flight-budget",
                "allocated temporal capacities exceed the budget",
            ));
        }
        Ok(result)
    }

    /// Full mechanical step, comprising two sequential wind transitions.
    #[must_use]
    pub fn dt_s(&self) -> f64 {
        2.0 * self.accepted.dt_s()
    }

    /// Physical time of the accepted wind and corresponding mechanical state.
    #[must_use]
    pub fn time_s(&self) -> f64 {
        self.accepted.time_s()
    }

    /// Conservative numerical capacity charged by `prepare`.
    #[must_use]
    pub fn workspace_bytes(&self) -> usize {
        self.accepted.workspace_bytes()
            + self.midpoint.workspace_bytes()
            + self.endpoint.workspace_bytes()
    }

    /// Observe the accepted field in world coordinates without advancing time.
    pub fn sample(&self, position_world_m: [f64; 3]) -> Result<FlowSample, Refusal> {
        sample_world(
            &self.accepted,
            self.mean,
            self.boundary,
            self.ground_down_m,
            position_world_m,
        )
    }

    /// Checkpoint only the accepted source and its physical interpretation.
    #[must_use]
    pub fn checkpoint(&self) -> MannFlightWindCheckpoint {
        MannFlightWindCheckpoint {
            version: 1,
            field: self.accepted.checkpoint(),
            mean: self.mean,
            boundary: self.boundary,
            ground_down_m: self.ground_down_m,
        }
    }

    /// Restore the identical future trajectory, using the same total budget.
    pub fn resume(
        checkpoint: MannFlightWindCheckpoint,
        maximum_workspace_bytes: usize,
    ) -> Result<Self, Refusal> {
        if checkpoint.version != 1 {
            return Err(invalid(
                "mann-flight-checkpoint",
                "unsupported flight wind checkpoint version",
            ));
        }
        let field = FixedStepMannField::resume(checkpoint.field, maximum_workspace_bytes / 3)
            .map_err(wind_refusal)?;
        Self::prepare(
            field,
            checkpoint.mean,
            checkpoint.boundary,
            checkpoint.ground_down_m,
            maximum_workspace_bytes,
        )
    }

    /// Advance the actual six-DOF carrier, geared shafts and uniform wakes.
    /// Select `BladeIncident::WindField` or `WindAndWingField` on an oblique
    /// propulsion model to resolve this field at every native blade section;
    /// the adapter preserves the caller's blade-resolution model choice.
    ///
    /// The work budget counts all eight normal innovations per spectral mode
    /// across both half steps. A wind, force, wake or cancellation refusal
    /// publishes neither a mechanical result nor a changed accepted wind.
    /// Input-callback side effects remain caller-owned. Publish `result.state`
    /// and this owner's checkpoint together after success.
    #[allow(clippy::too_many_arguments)]
    pub fn step_dynamic<const N: usize, C>(
        &mut self,
        model: &DynamicPropulsion<'_, '_, N>,
        state: &DynamicState,
        angular_origin_world: [f64; 3],
        maximum_mode_transitions: usize,
        cancelled: impl FnMut() -> bool,
        inputs: C,
    ) -> Result<DynamicStep<N>, Refusal>
    where
        C: FnMut(f64, &DynamicState) -> Result<PropulsionInputs, Refusal>,
    {
        self.with_stages(maximum_mode_transitions, cancelled, |stages| {
            let result = model.step(
                state,
                stages.times[0],
                stages.dt_s,
                angular_origin_world,
                inputs,
                |time, point| stages.sample(time, point).map(|flow| flow.u),
            )?;
            // Admit the accepted CG before publishing either owner. The next
            // force stage will separately admit its full material-point layout.
            stages.sample(stages.times[2], result.state.motion.rigid.pos_m)?;
            Ok(result)
        })
    }

    /// The same coherent stage history with the actual first-harmonic wakes.
    /// Refusal, work and publication semantics are those of `step_dynamic`.
    #[allow(clippy::too_many_arguments)]
    pub fn step_harmonic<const N: usize, C>(
        &mut self,
        model: &HarmonicPropulsion<'_, '_, N>,
        state: &HarmonicState,
        angular_origin_world: [f64; 3],
        maximum_mode_transitions: usize,
        cancelled: impl FnMut() -> bool,
        inputs: C,
    ) -> Result<HarmonicStep<N>, Refusal>
    where
        C: FnMut(f64, &HarmonicState) -> Result<PropulsionInputs, Refusal>,
    {
        self.with_stages(maximum_mode_transitions, cancelled, |stages| {
            let result = model.step(
                state,
                stages.times[0],
                stages.dt_s,
                angular_origin_world,
                inputs,
                |time, point| stages.sample(time, point).map(|flow| flow.u),
            )?;
            stages.sample(stages.times[2], result.state.motion.rigid.pos_m)?;
            Ok(result)
        })
    }

    fn with_stages<T>(
        &mut self,
        maximum_mode_transitions: usize,
        mut cancelled: impl FnMut() -> bool,
        solve: impl FnOnce(&Stages<'_>) -> Result<T, Refusal>,
    ) -> Result<T, Refusal> {
        poll(&mut cancelled)?;
        let half_work = 4 * self.accepted.model().mode_count();
        if 2 * half_work > maximum_mode_transitions {
            return Err(invalid(
                "mann-flight-work",
                "two wind half steps exceed the work budget",
            ));
        }
        let last = self
            .accepted
            .step()
            .checked_add(2)
            .ok_or_else(|| invalid("mann-flight-clock", "flight wind clock overflows"))?;
        let dt_s = self.dt_s();
        let time = self.time_s();
        // These are the native mechanical integrator's stage labels. The
        // source evaluates phase at its own canonical integer half-step time;
        // reassociation of floating-point additions does not interpolate OU paths.
        let times = [time, time + 0.5 * dt_s, time + dt_s];
        if !times.iter().all(|t| t.is_finite()) || times[0] >= times[1] || times[1] >= times[2] {
            return Err(invalid(
                "mann-flight-clock",
                "flight stage times are not representable",
            ));
        }
        self.midpoint
            .copy_temporal_state_from(&self.accepted)
            .map_err(wind_refusal)?;
        self.midpoint
            .advance_to(last - 1, half_work, &mut cancelled)
            .map_err(wind_refusal)?;
        self.endpoint
            .copy_temporal_state_from(&self.midpoint)
            .map_err(wind_refusal)?;
        self.endpoint
            .advance_to(last, half_work, &mut cancelled)
            .map_err(wind_refusal)?;
        let result = solve(&Stages {
            fields: [&self.accepted, &self.midpoint, &self.endpoint],
            times,
            dt_s,
            mean: self.mean,
            boundary: self.boundary,
            ground_down_m: self.ground_down_m,
        })?;
        poll(&mut cancelled)?;
        core::mem::swap(&mut self.accepted, &mut self.endpoint);
        Ok(result)
    }
}

struct Stages<'a> {
    fields: [&'a FixedStepMannField; 3],
    times: [f64; 3],
    dt_s: f64,
    mean: FlatSiteLogLaw,
    boundary: MannFlightBoundary,
    ground_down_m: f64,
}

impl Stages<'_> {
    fn sample(&self, time_s: f64, point: [f64; 3]) -> Result<FlowSample, Refusal> {
        let stage = self
            .times
            .iter()
            .position(|&time| time_s == time)
            .ok_or_else(|| {
                invalid(
                    "mann-flight-stage-time",
                    "wind query is not at this step's start, midpoint or endpoint",
                )
            })?;
        sample_world(
            self.fields[stage],
            self.mean,
            self.boundary,
            self.ground_down_m,
            point,
        )
    }
}

fn sample_world(
    field: &FixedStepMannField,
    mean: FlatSiteLogLaw,
    boundary: MannFlightBoundary,
    ground_down_m: f64,
    [x, y, z]: [f64; 3],
) -> Result<FlowSample, Refusal> {
    let height = ground_down_m - z;
    let minimum = match boundary {
        MannFlightBoundary::HomogeneousAbove { minimum_height_m } => minimum_height_m,
        MannFlightBoundary::GroundPlane => 0.0,
    };
    if !height.is_finite() || height < minimum {
        return Err(invalid(
            "mann-flight-applicability",
            "material point lies outside the wind's admitted height domain",
        ));
    }
    let point = [-x, y, height];
    let mut flow = match boundary {
        MannFlightBoundary::HomogeneousAbove { .. } => field.sample(point),
        MannFlightBoundary::GroundPlane => {
            field.sample_wall_blocked(point, MannPlaneWall::new(0.0).map_err(wind_refusal)?)
        }
    }
    .map_err(wind_refusal)?;
    flow.u[0] += mean.speed(height);
    flow.grad[0][2] += mean.dspeed_dh(height);
    let sign = [-1.0, 1.0, -1.0];
    for i in 0..3 {
        flow.u[i] *= sign[i];
        for j in 0..3 {
            flow.grad[i][j] *= sign[i] * sign[j];
        }
    }
    if !flow
        .u
        .iter()
        .chain(flow.grad.iter().flatten())
        .all(|x| x.is_finite())
    {
        return Err(invalid(
            "mann-flight-invalid",
            "world wind or gradient is not representable",
        ));
    }
    Ok(flow)
}

#[cfg(test)]
mod tests;
