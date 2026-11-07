//! Native Hess–Smith/Kutta coupling of the existing clockwise-positive wake.
//!
//! The retained vortices keep their circulation. One newly released packet
//! acquires the bound-circulation change through an algebraic Kelvin solve at
//! each midpoint stage. The same solved body field convects that wake.

use super::{
    MAX_WAKE_STEPS, WakeSim, WakeStep, WakeVortex, local_configuration, validate_direct_wake_work,
};
use crate::BemError;
use crate::panel2d::{Airfoil2d, PanelState2d, PreparedPanel2d};
use fs_vpm::{VortexParticle, induced_velocity};

#[derive(Debug, Clone)]
pub(super) struct PanelCoupling {
    panels: PreparedPanel2d,
    state: Option<PanelState2d>,
    steady_basis: [f64; 2],
}

/// Explicit transport refinement within ONE existing shedding interval.
/// Time and positions retain the wake's unit-chord/unit-freestream convention.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PanelWakeTransportOptions {
    /// Positive number of equal midpoint substeps. This does not change the
    /// shedding interval, release position, or number of emitted packets.
    pub substeps: usize,
    /// Whole-interval work cap: `2 * substeps + 1` coupled boundary solves.
    /// Each makes two solves with the existing LU; the
    /// original panel and vortex ceilings bound the work within each stage.
    pub max_coupled_solves: usize,
    /// Optional positive fraction in (0, 1] limiting `h * |u|` at BOTH stages
    /// against the smaller of the velocity core and distance to the body.
    /// This resolves near-core and near-body motion; it is not a temporal
    /// error estimate or continuous collision certificate. None disables it.
    pub max_displacement_fraction: Option<f64>,
}

/// Work and resolution observations from one completely published interval.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PanelWakeTransportReport {
    /// Completed midpoint substeps; exactly one history row is published.
    pub substeps: usize,
    /// Coupled initial/midpoint solves plus the final endpoint solve.
    pub coupled_solves: usize,
    /// Largest `h * |u| / min(core, body clearance)` over sampled stages.
    /// This is reported even when the optional resolution guard is disabled.
    pub max_displacement_fraction: f64,
}

impl WakeSim {
    /// Prepare a wake whose circulation is determined by the native panel
    /// normal-flow equations, Kutta condition, and Kelvin conservation.
    ///
    /// Coordinates, time and vortex cores use the same local unit-chord and
    /// unit-freestream normalization as [`Self::new`]. The Hess–Smith matrix
    /// is factored once. Each step reuses that factorization for the initial,
    /// midpoint and accepted-endpoint wake geometry. The new packet's strength
    /// and the bound sheet are solved together; no circulation relaxation or
    /// pressure-derived `Cl/2` normalization enters this mode.
    ///
    /// Explicit midpoint transport uses the complete solved source/vortex
    /// panel field and the original VPM wake kernel. It is second order for
    /// the smooth algebraically constrained motion within a release interval;
    /// packet placement and discrete shedding retain their own refinement
    /// error. This is inviscid attached-flow evolution, without an unsteady
    /// pressure/load, separation, viscous wall layer or collision model.
    ///
    /// # Errors
    /// Original geometry/time/core admissions, failed panel factorization,
    /// bounded allocation or non-finite numerics. Future steps also refuse
    /// when a sampled physical vortex lies on or inside the section.
    pub fn new_panel_coupled(
        foil: &Airfoil2d,
        alpha: f64,
        dt: f64,
        core: f64,
    ) -> Result<Self, BemError> {
        let local_alpha = local_configuration(foil, alpha, dt, core)?;
        let (leading, tangent, chord) = foil.chord_frame();
        let mut nodes = reserved(foil.nodes().len(), "normalized wake section")?;
        for p in foil.nodes() {
            let r = [(p[0] - leading[0]) / chord, (p[1] - leading[1]) / chord];
            nodes.push([
                r[0] * tangent[0] + r[1] * tangent[1],
                -r[0] * tangent[1] + r[1] * tangent[0],
            ]);
        }
        // Bind the numerical trailing edge to the wake's exact reference.
        nodes[0] = [1.0, 0.0];
        let normalized = Airfoil2d::new(nodes)?;
        let panels = PreparedPanel2d::new(&normalized)?;
        let freestream = [
            fs_math::det::cos(local_alpha),
            fs_math::det::sin(local_alpha),
        ];
        let mut incident = reserved(panels.midpoints().len(), "steady panel incident field")?;
        incident.resize(panels.midpoints().len(), [1.0, 0.0]);
        let steady_x = panels.solve_incident(&incident)?.bound_circulation;
        incident.fill([0.0, 1.0]);
        let steady_y = panels.solve_incident(&incident)?.bound_circulation;
        let steady_gamma = steady_x * freestream[0] + steady_y * freestream[1];
        let state = Self {
            alpha: local_alpha,
            dt,
            core2: core * core,
            wake: Vec::new(),
            history: Vec::new(),
            steady_gamma,
            circulation_slope: 0.0,
            te: [1.0, 0.0],
            panel_coupling: Some(PanelCoupling {
                panels,
                state: None,
                steady_basis: [steady_x, steady_y],
            }),
        };
        state.validate_state()?;
        Ok(state)
    }

    /// Accepted panel solution in the same local frame as [`Self::wake`].
    /// Its boundary residuals include every accepted wake vortex. `None`
    /// denotes the original screening mode or no completed native step yet.
    #[must_use]
    pub fn panel_state(&self) -> Option<&PanelState2d> {
        self.panel_coupling.as_ref()?.state.as_ref()
    }

    /// Advance a native panel-coupled wake at a prescribed freestream angle
    /// relative to its fixed chord. The angle is held during this interval;
    /// previously shed circulation is retained and the new packet carries
    /// the response to the changed incidence. This supplies a gust/incidence
    /// schedule without refactorization; it does not rotate or deform the body.
    ///
    /// # Errors
    /// Requires [`Self::new_panel_coupled`] and a finite angle in [-pi/2, pi/2].
    /// The original wake, angle, panel state and history survive any refusal.
    pub fn step_panel_coupled_at_angle(&mut self, alpha_chord: f64) -> Result<(), BemError> {
        self.admit_panel_incidence(alpha_chord)?;
        self.advance_panels(alpha_chord, 1, None, &mut || false)
            .map(|_| ())
    }

    /// Refine wake transport without emitting extra vortices or changing the
    /// prescribed shedding schedule. One packet is released before all the
    /// substeps; its strength and the body circulation satisfy Kelvin at every
    /// stage with the PREVIOUS interval's bound circulation as their target.
    /// Older packet strengths remain unchanged. The accepted body solution is
    /// recomputed at the final geometry and one history row is published.
    ///
    /// The caller selects the substep count and optional displacement guard.
    /// This is fixed refinement of the same second-order constrained ODE,
    /// not adaptive error control or refinement of the shedding discretization.
    /// One substep without a guard retains the original step's arithmetic.
    /// Wake and midpoint buffers are reused; the existing bounded panel solves
    /// retain their internal coefficient/surface-vector allocations.
    ///
    /// Cancellation is checked before/after each coupled boundary solve,
    /// within source-copy/incident loops and before each wake velocity query.
    /// A dense panel solve itself is not interruptible; its size is bounded by
    /// the original panel ceiling. Refusal publishes no physical state.
    ///
    /// # Errors
    /// Original native-wake refusals, invalid refinement/guard settings,
    /// insufficient solve budget, excessive stage motion or cancellation.
    pub fn step_panel_coupled_with_transport(
        &mut self,
        alpha_chord: f64,
        options: PanelWakeTransportOptions,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PanelWakeTransportReport, BemError> {
        checkpoint(&mut cancelled)?;
        self.admit_panel_incidence(alpha_chord)?;
        if options.substeps == 0 || !(self.dt / options.substeps as f64 > 0.0) {
            return Err(BemError::InvalidScalar {
                name: "panel wake transport substeps",
                value: options.substeps as f64,
                requirement: "positive count with a representable positive substep",
            });
        }
        if let Some(fraction) = options.max_displacement_fraction {
            if !fraction.is_finite() || fraction <= 0.0 || fraction > 1.0 {
                return Err(BemError::InvalidScalar {
                    name: "panel wake displacement fraction",
                    value: fraction,
                    requirement: "finite fraction in (0, 1]",
                });
            }
        }
        let solves = options
            .substeps
            .checked_mul(2)
            .and_then(|n| n.checked_add(1));
        if solves.is_none_or(|n| n > options.max_coupled_solves) {
            return Err(BemError::WorkEnvelopeExceeded {
                operation: "panel wake coupled solves",
                requested: solves.unwrap_or(usize::MAX),
                max: options.max_coupled_solves,
            });
        }
        self.advance_panels(
            alpha_chord,
            options.substeps,
            options.max_displacement_fraction,
            &mut cancelled,
        )
    }

    fn admit_panel_incidence(&self, alpha_chord: f64) -> Result<(), BemError> {
        self.validate_state()?;
        if self.panel_coupling.is_none() {
            return Err(BemError::InvalidScalar {
                name: "panel-coupled wake mode",
                value: 0.0,
                requirement: "construct with WakeSim::new_panel_coupled",
            });
        }
        if !alpha_chord.is_finite() || alpha_chord.abs() > super::MAX_ABS_WAKE_ALPHA {
            return Err(BemError::InvalidScalar {
                name: "panel wake chord-relative angle",
                value: alpha_chord,
                requirement: "finite radians in [-pi/2, pi/2]",
            });
        }
        Ok(())
    }

    pub(super) fn step_panels(&mut self) -> Result<(), BemError> {
        self.advance_panels(self.alpha, 1, None, &mut || false)
            .map(|_| ())
    }

    fn advance_panels(
        &mut self,
        alpha: f64,
        substeps: usize,
        max_displacement_fraction: Option<f64>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<PanelWakeTransportReport, BemError> {
        checkpoint(cancelled)?;
        if self.history.len() >= MAX_WAKE_STEPS {
            return Err(BemError::WorkEnvelopeExceeded {
                operation: "wake steps",
                requested: self.history.len().saturating_add(1),
                max: MAX_WAKE_STEPS,
            });
        }
        let count = self
            .wake
            .len()
            .checked_add(1)
            .ok_or(BemError::WorkEnvelopeExceeded {
                operation: "panel-coupled wake vortices",
                requested: usize::MAX,
                max: super::MAX_DIRECT_WAKE_VORTICES,
            })?;
        validate_direct_wake_work(count)?;
        let coupling = self.panel_coupling.as_ref().expect("native mode dispatch");
        let mut candidate = reserved(count, "panel wake candidate")?;
        candidate.extend_from_slice(&self.wake);
        candidate.push(WakeVortex {
            pos: [self.te[0] + 0.3 * self.dt, self.te[1]],
            gamma: 0.0,
        });
        let mut midpoint = reserved(count, "panel wake midpoint")?;
        midpoint.extend_from_slice(&candidate);
        let mut workspace = CouplingWorkspace::new(count, coupling.panels.midpoints().len())?;
        let previous_bound = self.history.last().map_or(0.0, |row| row.bound);
        let freestream = [fs_math::det::cos(alpha), fs_math::det::sin(alpha)];
        let steady_gamma =
            coupling.steady_basis[0] * freestream[0] + coupling.steady_basis[1] * freestream[1];
        let core = fs_math::det::sqrt(self.core2);
        let h = self.dt / substeps as f64;
        let mut peak = 0.0_f64;
        let mut maximum_fraction = 0.0_f64;
        for step in 0..substeps {
            let first = workspace.solve(
                &coupling.panels,
                &mut candidate,
                freestream,
                core,
                previous_bound,
                cancelled,
            )?;
            for (i, (p, q)) in candidate.iter().zip(&mut midpoint).enumerate() {
                checkpoint(cancelled)?;
                let v = workspace.velocity(&coupling.panels, &first, p.pos, freestream, core)?;
                if step == 0 {
                    peak = peak.max(fs_math::det::sqrt(v[0].mul_add(v[0], v[1] * v[1])));
                }
                maximum_fraction = maximum_fraction.max(stage_fraction(
                    v,
                    h,
                    core,
                    workspace.clearance[i],
                    max_displacement_fraction,
                )?);
                q.pos = translated(p.pos, v, 0.5 * h)?;
            }
            let middle = workspace.solve(
                &coupling.panels,
                &mut midpoint,
                freestream,
                core,
                previous_bound,
                cancelled,
            )?;
            for (i, (p, q)) in candidate.iter_mut().zip(&midpoint).enumerate() {
                checkpoint(cancelled)?;
                let v = workspace.velocity(&coupling.panels, &middle, q.pos, freestream, core)?;
                maximum_fraction = maximum_fraction.max(stage_fraction(
                    v,
                    h,
                    core,
                    workspace.clearance[i],
                    max_displacement_fraction,
                )?);
                p.pos = translated(p.pos, v, h)?;
            }
        }
        // The endpoint solve publishes a body state satisfying the boundary
        // equations for the accepted geometry, rather than a stale midpoint.
        let endpoint = workspace.solve(
            &coupling.panels,
            &mut candidate,
            freestream,
            core,
            previous_bound,
            cancelled,
        )?;
        let time = (self.history.len() as f64 + 1.0) * self.dt;
        if !time.is_finite() || !peak.is_finite() || !steady_gamma.is_finite() {
            return Err(BemError::NonFiniteResult {
                operation: "panel wake endpoint",
            });
        }
        if candidate.last().is_some_and(|vortex| vortex.gamma == 0.0) {
            candidate.pop();
        }
        self.history
            .try_reserve(1)
            .map_err(|_| BemError::AllocationFailed {
                operation: "panel wake history",
            })?;
        let row = WakeStep {
            t: time,
            bound: endpoint.bound_circulation,
            vortices: candidate.len(),
            peak_speed: peak,
        };
        checkpoint(cancelled)?;
        // Every fallible stage is complete; no partial physical state escapes.
        self.alpha = alpha;
        self.steady_gamma = steady_gamma;
        self.wake = candidate;
        self.panel_coupling
            .as_mut()
            .expect("native mode dispatch")
            .state = Some(endpoint);
        self.history.push(row);
        Ok(PanelWakeTransportReport {
            substeps,
            coupled_solves: 2 * substeps + 1,
            max_displacement_fraction: maximum_fraction,
        })
    }
}

struct CouplingWorkspace {
    particles: Vec<VortexParticle>,
    incident: Vec<[f64; 2]>,
    newborn: Vec<[f64; 2]>,
    clearance: Vec<f64>,
}

impl CouplingWorkspace {
    fn new(particles: usize, panels: usize) -> Result<Self, BemError> {
        let mut values = reserved(particles, "panel wake particle scratch")?;
        values.resize(particles, VortexParticle::new([0.0; 2], 0.0));
        let mut incident = reserved(panels, "panel wake incident scratch")?;
        let mut newborn = reserved(panels, "panel wake newborn scratch")?;
        incident.resize(panels, [0.0; 2]);
        newborn.resize(panels, [0.0; 2]);
        let mut clearance = reserved(particles, "panel wake body clearance")?;
        clearance.resize(particles, 0.0);
        Ok(Self {
            particles: values,
            incident,
            newborn,
            clearance,
        })
    }

    fn solve(
        &mut self,
        panels: &PreparedPanel2d,
        wake: &mut [WakeVortex],
        freestream: [f64; 2],
        core: f64,
        previous_bound: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<PanelState2d, BemError> {
        for (i, (vortex, particle)) in wake.iter().zip(&mut self.particles).enumerate() {
            if i % 64 == 0 {
                checkpoint(cancelled)?;
            }
            self.clearance[i] = admit_exterior(panels.nodes(), vortex.pos)?;
            *particle = vortex.as_particle();
        }
        let retained = wake.len() - 1;
        let unit = WakeVortex {
            pos: wake[retained].pos,
            gamma: 1.0,
        }
        .as_particle();
        for (i, &point) in panels.midpoints().iter().enumerate() {
            if i % 64 == 0 {
                checkpoint(cancelled)?;
            }
            let induced = induced_velocity(&self.particles[..retained], point, core);
            self.incident[i] = [freestream[0] + induced[0], freestream[1] + induced[1]];
            self.newborn[i] = induced_velocity(core::slice::from_ref(&unit), point, core);
        }
        checkpoint(cancelled)?;
        let (state, shed) = panels.solve_kelvin(&self.incident, &self.newborn, previous_bound)?;
        checkpoint(cancelled)?;
        wake[retained].gamma = shed;
        self.particles[retained].circulation = -shed;
        Ok(state)
    }

    fn velocity(
        &self,
        panels: &PreparedPanel2d,
        state: &PanelState2d,
        point: [f64; 2],
        freestream: [f64; 2],
        core: f64,
    ) -> Result<[f64; 2], BemError> {
        let body = panels.body_velocity(state, point)?;
        let wake = induced_velocity(&self.particles, point, core);
        let velocity = [
            freestream[0] + body[0] + wake[0],
            freestream[1] + body[1] + wake[1],
        ];
        if velocity.iter().any(|x| !x.is_finite()) {
            return Err(BemError::NonFiniteResult {
                operation: "panel wake velocity",
            });
        }
        Ok(velocity)
    }
}

fn reserved<T>(count: usize, operation: &'static str) -> Result<Vec<T>, BemError> {
    let mut values = Vec::new();
    values
        .try_reserve_exact(count)
        .map_err(|_| BemError::AllocationFailed { operation })?;
    Ok(values)
}

fn translated(point: [f64; 2], velocity: [f64; 2], dt: f64) -> Result<[f64; 2], BemError> {
    let next = [point[0] + dt * velocity[0], point[1] + dt * velocity[1]];
    if next.iter().any(|x| !x.is_finite()) {
        return Err(BemError::NonFiniteResult {
            operation: "panel wake transport",
        });
    }
    Ok(next)
}

fn checkpoint(cancelled: &mut impl FnMut() -> bool) -> Result<(), BemError> {
    if cancelled() {
        Err(BemError::Cancelled {
            operation: "panel wake transport",
        })
    } else {
        Ok(())
    }
}

fn stage_fraction(
    velocity: [f64; 2],
    dt: f64,
    core: f64,
    clearance: f64,
    maximum: Option<f64>,
) -> Result<f64, BemError> {
    let ratio = fs_math::det::hypot(dt * velocity[0], dt * velocity[1]) / core.min(clearance);
    if !ratio.is_finite() {
        return Err(BemError::NonFiniteResult {
            operation: "panel wake stage displacement",
        });
    }
    if maximum.is_some_and(|maximum| ratio > maximum) {
        return Err(BemError::InvalidScalar {
            name: "panel wake stage displacement fraction",
            value: ratio,
            requirement: "no larger than the declared fraction; increase transport substeps",
        });
    }
    Ok(ratio)
}

// A closed polygon test in the native unit-chord frame. The short boundary
// tolerance rejects unresolved panel contact rather than sampling its jump.
// This admits sampled stages only; it is not a continuous collision detector.
// Return closest-body clearance from the SAME pass for the resolution guard.
fn admit_exterior(nodes: &[[f64; 2]], point: [f64; 2]) -> Result<f64, BemError> {
    if point.iter().any(|x| !x.is_finite()) {
        return Err(BemError::NonFiniteResult {
            operation: "panel wake position",
        });
    }
    let mut inside = false;
    let mut clearance2 = f64::INFINITY;
    for i in 0..nodes.len() {
        let a = nodes[i];
        let b = nodes[(i + 1) % nodes.len()];
        let edge = [b[0] - a[0], b[1] - a[1]];
        let r = [point[0] - a[0], point[1] - a[1]];
        let length2 = edge[0].mul_add(edge[0], edge[1] * edge[1]);
        let projection = ((r[0] * edge[0] + r[1] * edge[1]) / length2).clamp(0.0, 1.0);
        let dx = r[0] - projection * edge[0];
        let dy = r[1] - projection * edge[1];
        let distance2 = dx.mul_add(dx, dy * dy);
        clearance2 = clearance2.min(distance2);
        if distance2 <= 1e-24 {
            return Err(BemError::InvalidScalar {
                name: "panel wake exterior clearance",
                value: 0.0,
                requirement: "every sampled vortex must be outside the section and its 1e-12 chord contact layer",
            });
        }
        if (a[1] > point[1]) != (b[1] > point[1])
            && point[0] < a[0] + (point[1] - a[1]) * (b[0] - a[0]) / (b[1] - a[1])
        {
            inside = !inside;
        }
    }
    if inside {
        return Err(BemError::InvalidScalar {
            name: "panel wake exterior clearance",
            value: -1.0,
            requirement: "every sampled vortex must remain outside the airfoil",
        });
    }
    Ok(fs_math::det::sqrt(clearance2))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::panel2d::naca4_symmetric;

    #[test]
    fn refused_newborn_body_contact_leaves_the_native_state_unchanged() {
        let foil = naca4_symmetric(0.12, 24).unwrap();
        let mut wake = WakeSim::new_panel_coupled(&foil, 0.08, 0.02, 0.03).unwrap();
        wake.step().unwrap();
        // Inject a release-point failure into the real step transaction.
        wake.te = [0.5, 0.0];
        let before = wake.clone();
        assert!(wake.step().is_err());
        assert_eq!(wake.trace_json(1).unwrap(), before.trace_json(1).unwrap());
        assert_eq!(
            wake.panel_state().unwrap().sources,
            before.panel_state().unwrap().sources
        );
        for (a, b) in wake.wake().iter().zip(before.wake()) {
            assert_eq!(a.pos, b.pos);
            assert_eq!(a.gamma, b.gamma);
        }
    }
}
