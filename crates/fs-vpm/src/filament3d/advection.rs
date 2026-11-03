//! Connected material-vertex wake evolution. Both the near-only and hybrid
//! near/far paths share the same finite-core midpoint stages and admission.

use super::{FilamentWake, Refusal, ShedRow};
use crate::farfield3d::WakeCoreEvolutionMode;
use crate::images3d::ImagePlane;

#[path = "advection/stage.rs"]
pub(crate) mod stage;
use stage::{PreparedSources, StageSettings, integrate_nodes};

/// Prescribed geometry and stationary slip boundary for one shedding interval.
#[derive(Clone, Copy, Debug, Default)]
pub struct WakeBoundary<'a> {
    /// Endpoint lifting-line nodes [m]. Each node follows a linear trajectory
    /// from its current position; `None` holds the line fixed. The new row is
    /// emitted at this endpoint. Free wake nodes are never pinned to the line.
    pub line_endpoint_m: Option<&'a [[f64; 3]]>,
    /// Stationary horizontal plane. Every physical node must remain strictly
    /// above it at both RK stages and the endpoint. Mirrored sources enforce
    /// zero induced normal velocity at the plane and never enter physical
    /// state or ledgers. The caller supplies a boundary-compatible incident
    /// field; images do not cancel an imposed normal flow.
    pub image_plane: Option<ImagePlane>,
}

/// Inputs to one self-induced shedding interval of length `1 / shed_hz`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SelfInducedShed {
    /// Uniform imposed fluid velocity [m/s], added to wake induction.
    pub ambient_velocity_mps: [f64; 3],
    /// Positive base radius [m] assigned to the row emitted at the endpoint.
    /// Existing rows retain their own base cores.
    pub emitted_core_radius_m: f64,
    /// Same age-dependent algebraic spreading law used by the far field.
    pub core_mode: WakeCoreEvolutionMode,
    /// Fixed positive number of midpoint RK2 substeps during this interval.
    pub substeps: usize,
    /// Maximum full-substep displacement divided by the smaller of the
    /// smallest active core radius and resolved active segment length, in (0, 1].
    /// Checked at both RK stages; this is a resolution guard, not an error bound.
    /// Edges born with zero length remain subject to the core guard and join
    /// the edge-length guard in the next shedding interval.
    pub max_displacement_fraction: f64,
    /// Upper bound on direct interaction work, including source assembly and
    /// target updates and final source admission:
    /// `(2 * substeps * (targets + 1) + 1) * (source_bound + 1)`.
    /// The source bound includes zero-strength edges and doubles for images,
    /// so admission is conservative.
    pub max_work: u64,
}

impl FilamentWake {
    /// Convect the complete near wake by finite-core self-induction plus uniform
    /// ambient flow, hold the lifting line fixed, and emit one endpoint row.
    ///
    /// This is the default-environment form of
    /// [`Self::shed_self_induced_with_environment`]. Attached far interfaces
    /// require the corresponding coupled [`crate::farfield3d::FarField`] method.
    ///
    /// # Errors
    /// Invalid geometry, time/core settings, attached far interface, exceeded
    /// work/displacement limits, non-finite numerics, or cancellation. Errors
    /// leave the original wake unchanged.
    pub fn shed_self_induced(
        &mut self,
        gamma: &[f64],
        config: SelfInducedShed,
        cancelled: impl FnMut() -> bool,
    ) -> Result<(), Refusal> {
        self.shed_self_induced_with_environment(
            gamma,
            config,
            WakeBoundary::default(),
            |_, _| Ok([0.0; 3]),
            cancelled,
        )
    }

    /// Evolve the held-strength material sheet with explicit midpoint RK2,
    /// prescribed line motion, an incident-flow variation, and optional images.
    ///
    /// `incident_variation(point, elapsed_seconds)` ADDS to the uniform ambient
    /// velocity. It is sampled at each actual initial/midpoint position and
    /// stage time within this interval. It must not include this wake's own
    /// induction or its images. Callback errors propagate without publication;
    /// side effects inside callbacks remain the caller's responsibility.
    ///
    /// Both stages rebuild the complete source geometry and elapsed-age cores.
    /// Existing circulations and base cores are retained; ages advance once.
    /// `gamma` is captured in one age-zero row at the endpoint lifting line.
    /// Second-order accuracy applies to the smooth held-strength problem, not
    /// a discontinuous change of boundary circulation. Linear node trajectories
    /// represent prescribed deformation/translation, not exact finite rotation.
    ///
    /// The prescribed line exchanges impulse with the wake. This direct
    /// O(nodes × sources × substeps) reference does not assert exact impulse or
    /// energy conservation, resolve viscous tube strain, or remesh automatically.
    /// Cancellation is checked before work, per substep/target, within every
    /// 256 source-target contributions, and before publication. Source assembly
    /// is a bounded linear pass between those checkpoints.
    ///
    /// # Errors
    /// The same refusals as [`Self::shed_self_induced`], plus invalid line motion,
    /// an incident-field error, or a stage/endpoint image-plane crossing.
    pub fn shed_self_induced_with_environment(
        &mut self,
        gamma: &[f64],
        config: SelfInducedShed,
        boundary: WakeBoundary<'_>,
        mut incident_variation: impl FnMut([f64; 3], f64) -> Result<[f64; 3], Refusal>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<(), Refusal> {
        checkpoint(&mut cancelled)?;
        if self.downstream_gamma.is_some() {
            return Err(refusal(
                "wake-advection-interface-attached",
                "advance an attached interface with the coupled near/far method",
            ));
        }
        let proposal = prepare_near(self, gamma, config)?;
        let mut state = proposal.near;
        let line = LineMotion::new(&state, boundary)?;
        let mut points: Vec<_> = state
            .rows
            .iter()
            .flat_map(|r| r.nodes.iter().copied())
            .collect();
        let settings = StageSettings {
            dt_s: proposal.dt_s,
            source_bound: state.rows.len() * (2 * state.cert.n_stations + 1)
                + state.cert.n_stations,
            image_plane: boundary.image_plane,
            line_step_m: line.step_distance(config.substeps),
        };
        integrate_nodes(
            &mut points,
            config,
            settings,
            |stage_points, fraction, sources| {
                scatter(&mut state, stage_points);
                line.apply(&mut state, fraction)?;
                append_near_sources(&state, proposal.growth_per_tick_m2, fraction, sources)
            },
            &mut incident_variation,
            &mut cancelled,
        )?;
        scatter(&mut state, &points);
        line.apply(&mut state, 1.0)?;
        for (row, original) in state.rows.iter_mut().zip(&self.rows) {
            row.age_ticks = original.age_ticks + 1.0;
        }
        let mut emitted = proposal.emitted;
        emitted.nodes.clone_from(&state.line_nodes);
        state.rows.push(emitted);
        checkpoint(&mut cancelled)?;
        *self = state;
        Ok(())
    }
}

fn scatter(near: &mut FilamentWake, points: &[[f64; 3]]) {
    for (p, q) in near.rows.iter_mut().flat_map(|r| &mut r.nodes).zip(points) {
        *p = *q;
    }
}

pub(crate) struct ShedProposal {
    pub near: FilamentWake,
    pub emitted: ShedRow,
    pub dt_s: f64,
    pub growth_per_tick_m2: f64,
}

/// Prepare a private endpoint emission but keep the old sheet at its initial
/// age. Valid attached interfaces are admitted here for the hybrid owner.
pub(crate) fn prepare_near(
    near: &FilamentWake,
    gamma: &[f64],
    config: SelfInducedShed,
) -> Result<ShedProposal, Refusal> {
    near.cert.admit()?;
    near.admit_geometry()?;
    config.core_mode.admit()?;
    let dt = 1.0 / near.cert.shed_hz;
    let h = dt / config.substeps as f64;
    let emitted_core2 = config.emitted_core_radius_m.powi(2);
    let growth = config.core_mode.core_growth_m2(dt);
    if config.substeps == 0
        || !h.is_finite()
        || h <= 0.0
        || !dt.is_finite()
        || !growth.is_finite()
        || config.ambient_velocity_mps.iter().any(|v| !v.is_finite())
        || !config.emitted_core_radius_m.is_finite()
        || config.emitted_core_radius_m <= 0.0
        || !emitted_core2.is_finite()
        || emitted_core2 <= 0.0
        || emitted_core2 > f64::MAX / 16.0
        || !config.max_displacement_fraction.is_finite()
        || !(0.0 < config.max_displacement_fraction && config.max_displacement_fraction <= 1.0)
    {
        return Err(refusal(
            "wake-advection-invalid",
            "invalid time, core or step settings",
        ));
    }
    if near.downstream_gamma.as_ref().is_some_and(|g| {
        near.rows.is_empty() || g.len() != near.cert.n_stations || g.iter().any(|v| !v.is_finite())
    }) {
        return Err(refusal(
            "wake-advection-topology-invalid",
            "invalid attached boundary strength",
        ));
    }
    for row in &near.rows {
        let final_core = row.core2_m2 + growth * (row.age_ticks + 1.0);
        if !row.core2_m2.is_finite()
            || row.core2_m2 < 0.0
            || !row.age_ticks.is_finite()
            || row.age_ticks < 0.0
            || !final_core.is_finite()
        {
            return Err(refusal(
                "wake-advection-core-invalid",
                "invalid retained age or core",
            ));
        }
    }
    let mut state = near.clone();
    state.shed(gamma, [0.0; 3])?;
    let mut emitted = state.rows.pop().expect("shed appended one row");
    emitted.core2_m2 = emitted_core2;
    if emitted.gamma.windows(2).any(|g| !(g[1] - g[0]).is_finite())
        || state.rows.last().is_some_and(|row| {
            row.gamma
                .iter()
                .zip(&emitted.gamma)
                .any(|(a, b)| !(a - b).is_finite())
        })
    {
        return Err(refusal(
            "wake-advection-circulation-invalid",
            "emission would overflow a shared-edge circulation",
        ));
    }
    for (row, original) in state.rows.iter_mut().zip(&near.rows) {
        row.age_ticks = original.age_ticks;
    }
    Ok(ShedProposal {
        near: state,
        emitted,
        dt_s: dt,
        growth_per_tick_m2: growth,
    })
}

pub(crate) fn append_near_sources(
    near: &FilamentWake,
    growth: f64,
    fraction: f64,
    sources: &mut PreparedSources,
) -> Result<(), Refusal> {
    let mut result = Ok(());
    near.visit_segments_with_core_growth(growth, |s| {
        if result.is_ok() {
            result = sources.push(
                s.start,
                s.end,
                s.circulation,
                s.core2_m2 + growth * fraction,
            );
        }
    });
    result
}

pub(crate) struct LineMotion<'a> {
    start: Vec<[f64; 3]>,
    end: Option<&'a [[f64; 3]]>,
    plane: Option<ImagePlane>,
}

impl<'a> LineMotion<'a> {
    pub(crate) fn new(near: &FilamentWake, boundary: WakeBoundary<'a>) -> Result<Self, Refusal> {
        if let Some(plane) = boundary.image_plane {
            plane.admit(near)?;
        }
        if boundary.line_endpoint_m.is_some_and(|end| {
            end.len() != near.line_nodes.len()
                || end.iter().zip(&near.line_nodes).any(|(q, p)| {
                    q.iter().any(|v| !v.is_finite())
                        || !distance(*p, *q).is_finite()
                        || boundary.image_plane.is_some_and(|plane| q[2] <= plane.z0_m)
                })
        }) {
            return Err(refusal(
                "wake-advection-line-invalid",
                "invalid endpoint lifting line",
            ));
        }
        Ok(Self {
            start: near.line_nodes.clone(),
            end: boundary.line_endpoint_m,
            plane: boundary.image_plane,
        })
    }

    pub(crate) fn step_distance(&self, substeps: usize) -> f64 {
        self.end.map_or(0.0, |end| {
            self.start
                .iter()
                .zip(end)
                .map(|(p, q)| distance(*p, *q))
                .fold(0.0, f64::max)
                / substeps as f64
        })
    }

    pub(crate) fn apply(&self, near: &mut FilamentWake, fraction: f64) -> Result<(), Refusal> {
        let Some(end) = self.end else { return Ok(()) };
        for ((p, start), end) in near.line_nodes.iter_mut().zip(&self.start).zip(end) {
            *p = if fraction == 0.0 {
                *start
            } else if fraction == 1.0 {
                *end
            } else {
                std::array::from_fn(|i| start[i] + fraction * (end[i] - start[i]))
            };
            if p.iter().any(|v| !v.is_finite())
                || self.plane.is_some_and(|plane| p[2] <= plane.z0_m)
            {
                return Err(refusal(
                    "wake-advection-line-invalid",
                    "unrepresentable line trajectory",
                ));
            }
        }
        Ok(())
    }
}

fn distance(p: [f64; 3], q: [f64; 3]) -> f64 {
    (q[0] - p[0]).hypot(q[1] - p[1]).hypot(q[2] - p[2])
}

pub(crate) fn checkpoint(cancelled: &mut impl FnMut() -> bool) -> Result<(), Refusal> {
    if cancelled() {
        Err(refusal(
            "wake-advection-cancelled",
            "self-induced shedding was cancelled",
        ))
    } else {
        Ok(())
    }
}

pub(crate) fn refusal(code: &'static str, message: &str) -> Refusal {
    Refusal {
        code,
        message: message.into(),
        ranked_repairs: vec![
            "use a complete finite-core wake and admitted integration settings".into(),
        ],
    }
}
