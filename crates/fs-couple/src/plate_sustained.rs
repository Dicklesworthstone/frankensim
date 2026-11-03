//! Finite-interval frictionless contact between two original linear plates.
//!
//! The original sustained force law sees the endpoint gaps of actual native
//! generalized-alpha steps at rho_inf=1 (the midpoint setting). Contact forces
//! are held over each step; they are never inferred from an impact impulse.
//! Initially closed, stationary normal contact is required. Collision location,
//! rotating normals, friction and heat conversion belong to other owners.
//!
//! Work and Rayleigh dissipation are midpoint quadratures. The quadratic gap
//! interpolant is the declared discrete path, not an exact continuous orbit.
//! Both trajectories are published only after independent path, momentum,
//! support and work checks. Mass-proportional damping exchanges momentum with
//! the modeled environment and is reported explicitly.

use crate::plate_impact::PlateRigidImpactError;
use crate::plate_impact::multiple::{
    PlateImpactSite, PlateRigidImpactSet, PlateRigidImpactSetBudget, PlateRigidImpactSetControls,
};
use fs_contact::sustained::NaturalMapScales;
use fs_contact::sustained::solve::{self, SolveControls, SolveError, SolveReport};
use fs_math::dd::Dd;
use fs_plate::PlateMesh;
use fs_plate_transient::{PlateDynamics, PlateLoad};
use fs_sparse::Csr;
use fs_time::galpha::initialization::{
    InitialSolveError, initial_workspace_components, second_order_acceleration,
};
use fs_time::galpha::second_order_adjoint::trajectory::StructuralTrajectoryModel;
use fs_time::galpha::{
    ImplicitSolveConfig, OperatorGeneralizedAlpha, SecondOrderState, TimeSolveError,
};

/// Native spatial model and its original complete, uneliminated matrices.
/// Supports are stationary; the undeformed faces use one world +z normal.
pub struct PlateSustainedSurface<'a, L: ?Sized> {
    pub dynamics: PlateDynamics<'a, L>,
    pub mesh: &'a PlateMesh,
    pub full_mass: &'a Csr,
    pub full_stiffness: &'a Csr,
    pub reference_height_m: f64,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PlateSustainedSite {
    pub point_xy: [f64; 2],
    /// Conditioning of the original natural map, never a penalty stiffness.
    pub scales: NaturalMapScales,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PlateSustainedControls {
    /// Reuses original point/mass admission and dimensional balance tolerances.
    /// Its impact-law controls apply only to preparation, not the force law.
    pub mechanics: PlateRigidImpactSetControls,
    pub step: ImplicitSolveConfig,
    pub contact: SolveControls,
    pub gap_tolerance_m: f64,
    /// Absolute native midpoint kinematic residual, per generalized coordinate
    /// (metres for w, radians for the two slope coordinates).
    pub kinematic_tolerance: f64,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct PlateSustainedBudget {
    pub positive: PlateRigidImpactSetBudget,
    pub negative: PlateRigidImpactSetBudget,
    /// Sum of entries in both complete stiffness matrices.
    pub max_stiffness_nonzeros: usize,
    /// Live scalar storage including both preparation ceilings, force-response
    /// states, native time solves and audits; excludes borrowed inputs and the
    /// original contact solver's separately bounded workspace.
    pub max_workspace_components: usize,
    /// Maximum accepted history entries in each continuing native trajectory.
    pub max_history_entries: usize,
}

/// Complete generalized source samples, including fixed rows. Retained entries
/// must exactly match the original load callback at start, midpoint and end.
/// Between those samples only the native midpoint source quadrature is claimed.
#[derive(Clone, Copy, Debug)]
pub struct PlateSustainedLoads<'a> {
    pub start: &'a [f64],
    pub midpoint: &'a [f64],
    pub end: &'a [f64],
}
impl<'a> PlateSustainedLoads<'a> {
    #[must_use]
    pub fn held(full: &'a [f64]) -> Self {
        Self {
            start: full,
            midpoint: full,
            end: full,
        }
    }
}

#[derive(Debug)]
pub enum PlateSustainedSourceError {
    Invalid(&'static str),
    Budget(&'static str),
    Preparation(PlateRigidImpactError),
    Initialization(InitialSolveError),
    Step(TimeSolveError),
    Load(String),
    Balance {
        what: &'static str,
        residual: f64,
        allowance: f64,
    },
    Cancelled,
}
#[derive(Debug)]
pub enum PlateSustainedError {
    Source(PlateSustainedSourceError),
    Contact(SolveError<PlateSustainedSourceError>),
}
impl From<PlateSustainedSourceError> for PlateSustainedError {
    fn from(error: PlateSustainedSourceError) -> Self {
        Self::Source(error)
    }
}
impl core::fmt::Display for PlateSustainedError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "sustained plate contact refused: {self:?}")
    }
}
impl core::error::Error for PlateSustainedError {}

/// Resultants are [Pz, Lx, Ly] about the original world origin. The plate model
/// does not represent in-plane momentum or drilling angular momentum.
#[derive(Debug)]
pub struct PlateSustainedBalance {
    /// Complete stationary-support impulse acting on each fixed native row.
    pub support_impulses: Vec<(usize, f64)>,
    pub momentum_change: [f64; 3],
    /// Source impulse acting on the plate, using the native midpoint sample.
    pub external_impulse: [f64; 3],
    /// Resistive h*C*v_mid resultant. The modeled damping impulse acting ON
    /// the plate is its negative; the balance therefore adds this quantity.
    pub damping_impulse: [f64; 3],
    /// Complete stationary-support resultant acting on the plate.
    pub support_impulse: [f64; 3],
    pub kinetic_energy_change_j: f64,
    pub elastic_energy_change_j: f64,
    pub external_work_midpoint_j: f64,
    pub damping_midpoint_j: f64,
    pub contact_work_j: f64,
    pub energy_residual_j: f64,
}

#[derive(Debug)]
pub struct PlateSustainedReport {
    pub contact: SolveReport,
    pub gap_bounds_m: Vec<[f64; 2]>,
    pub relative_normal_speeds_after: Vec<f64>,
    pub positive: PlateSustainedBalance,
    pub negative: PlateSustainedBalance,
    /// Combined momentum + damping - external - support; internal contact
    /// cancels only after contributions from both actual native updates.
    pub momentum_balance_residual: [f64; 3],
    pub energy_residual_j: f64,
}

struct Surface<'a, L: ?Sized> {
    native: PlateSustainedSurface<'a, L>,
    points: PlateRigidImpactSet<'a, L>,
}

pub struct PlateSustainedPair<'a, P: ?Sized, N: ?Sized> {
    positive: Surface<'a, P>,
    negative: Surface<'a, N>,
    sites: Vec<PlateSustainedSite>,
    controls: PlateSustainedControls,
    budget: PlateSustainedBudget,
}

type SourceResult<T> = Result<T, PlateSustainedSourceError>;
use PlateSustainedSourceError as Source;

impl<'a, P: PlateLoad + ?Sized, N: PlateLoad + ?Sized> PlateSustainedPair<'a, P, N> {
    /// Prepare original point transpose rows and admit the complete matrices.
    /// The cached impact mobilities are not used as a finite-step response.
    pub fn prepare(
        positive: PlateSustainedSurface<'a, P>,
        negative: PlateSustainedSurface<'a, N>,
        sites: &[PlateSustainedSite],
        controls: PlateSustainedControls,
        budget: PlateSustainedBudget,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, PlateSustainedError> {
        poll(cancelled)?;
        if sites.is_empty()
            || sites.len() > solve::MAX_CONTACTS
            || sites.len() > controls.contact.maximum_contacts
        {
            return Err(Source::Budget("sustained contact count").into());
        }
        if !positive.reference_height_m.is_finite()
            || !negative.reference_height_m.is_finite()
            || !controls.gap_tolerance_m.is_finite()
            || controls.gap_tolerance_m < 0.0
            || !controls.kinematic_tolerance.is_finite()
            || controls.kinematic_tolerance < 0.0
            || controls.step.max_newton_iterations == 0
            || controls.step.max_newton_iterations > 64
            || controls.step.newton.linear_restart > 64
            || controls.step.newton.max_linear_cycles > 128
            || controls.step.newton.validate().is_err()
        {
            return Err(
                Source::Invalid("geometry, tolerances or bounded native step controls").into(),
            );
        }
        let stiffness = positive
            .full_stiffness
            .nnz()
            .checked_add(negative.full_stiffness.nnz());
        if stiffness.is_none_or(|n| n > budget.max_stiffness_nonzeros) {
            return Err(Source::Budget("complete stiffness entries").into());
        }
        let n = positive
            .dynamics
            .model()
            .free
            .checked_add(negative.dynamics.model().free);
        let full = positive
            .dynamics
            .model()
            .dof_map
            .len()
            .checked_add(negative.dynamics.model().dof_map.len());
        let scratch = (|| {
            let n = n?;
            let full = full?;
            let restart = controls.step.newton.linear_restart;
            let init = initial_workspace_components(n, 0, controls.mechanics.impulse.mass_solve)?;
            budget
                .positive
                .plate
                .max_workspace_components
                .checked_add(budget.negative.plate.max_workspace_components)?
                .checked_add(init.checked_mul(2)?)?
                .checked_add(n.checked_mul(128 + 4 * restart)?)?
                .checked_add(full.checked_mul(80)?)?
                .checked_add((restart + 1).checked_mul(restart)?.checked_mul(4)?)?
                .checked_add(controls.step.max_newton_iterations.checked_mul(256)?)?
                .checked_add(
                    (restart * controls.step.newton.max_linear_cycles + 1).checked_mul(16)?,
                )?
                .checked_add(sites.len().checked_mul(128)?)
        })();
        if scratch.is_none_or(|n| {
            n > budget.max_workspace_components
                || n > isize::MAX as usize / core::mem::size_of::<f64>()
        }) {
            return Err(Source::Budget("combined native workspace").into());
        }
        let impact_sites: Vec<_> = sites
            .iter()
            .map(|site| PlateImpactSite {
                point_xy: site.point_xy,
                activation_tolerance_m: controls.gap_tolerance_m,
                restitution: 0.0,
            })
            .collect();
        for site in sites {
            if [site.scales.normal_n_per_m, site.scales.tangent_n_per_m]
                .iter()
                .any(|value| !value.is_finite() || *value <= 0.0)
            {
                return Err(Source::Invalid("natural-map conditioning").into());
            }
        }
        let positive = Surface::prepare(
            positive,
            &impact_sites,
            controls,
            budget.positive,
            cancelled,
        )?;
        let negative = Surface::prepare(
            negative,
            &impact_sites,
            controls,
            budget.negative,
            cancelled,
        )?;
        Ok(Self {
            positive,
            negative,
            sites: sites.to_vec(),
            controls,
            budget,
        })
    }

    /// Advance closed contact over one interval. A separated endpoint is
    /// admissible only at zero force; locating any later activation is the
    /// caller's responsibility. No original state changes on refusal.
    #[allow(clippy::too_many_arguments)]
    pub fn step(
        &self,
        positive: &mut SecondOrderState,
        negative: &mut SecondOrderState,
        step_seconds: f64,
        positive_loads: PlateSustainedLoads<'_>,
        negative_loads: PlateSustainedLoads<'_>,
        initial_forces_n: &[f64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<PlateSustainedReport, PlateSustainedError> {
        poll(cancelled)?;
        let h = step_seconds;
        if !h.is_finite()
            || h <= 0.0
            || positive.t != negative.t
            || !(positive.t + h).is_finite()
            || positive.t + h <= positive.t
            || initial_forces_n.len() != self.sites.len()
            || initial_forces_n.iter().any(|x| !x.is_finite())
        {
            return Err(Source::Invalid(
                "common time, positive representable interval or force shape",
            )
            .into());
        }
        self.positive
            .admit_state(positive, self.budget.max_history_entries)?;
        self.negative
            .admit_state(negative, self.budget.max_history_entries)?;
        self.positive
            .admit_loads(positive.t, h, positive_loads, cancelled)?;
        self.negative
            .admit_loads(negative.t, h, negative_loads, cancelled)?;
        let mut before = Vec::with_capacity(self.sites.len());
        for i in 0..self.sites.len() {
            let (gap, speed) = self.gap_speed(i, positive, negative);
            check("initial closed gap", gap, self.controls.gap_tolerance_m)?;
            check(
                "initial zero normal speed (impacts require event owner)",
                speed,
                self.controls.mechanics.speed_tolerance_m_s,
            )?;
            before.push((gap, speed));
        }
        let scales: Vec<_> = self.sites.iter().map(|site| site.scales).collect();
        let friction = vec![0.0; self.sites.len()];
        let initial: Vec<_> = initial_forces_n.iter().map(|&x| [x, 0.0, 0.0]).collect();
        let mut candidate = None;
        let contact = solve::solve(
            &friction,
            &scales,
            &initial,
            self.controls.contact,
            &mut |forces, stop| {
                // The last callback is the law owner's mandatory fresh replay.
                // Keep only its two small native states, never clone histories.
                let p = self.positive.response(
                    positive,
                    h,
                    positive_loads,
                    forces,
                    1.0,
                    self.controls,
                    self.budget.max_workspace_components,
                    &mut || stop(),
                )?;
                let n = self.negative.response(
                    negative,
                    h,
                    negative_loads,
                    forces,
                    -1.0,
                    self.controls,
                    self.budget.max_workspace_components,
                    &mut || stop(),
                )?;
                let paths: Vec<_> = (0..self.sites.len())
                    .map(|i| [self.gap_speed(i, &p.state, &n.state).0, 0.0, 0.0])
                    .collect();
                candidate = Some((p, n, forces.to_vec()));
                Ok(paths)
            },
            cancelled,
        )
        .map_err(PlateSustainedError::Contact)?;
        let (mut p, mut n, replay_forces) =
            candidate.ok_or(Source::Invalid("missing final response"))?;
        if replay_forces != contact.forces_n {
            return Err(Source::Invalid("force-law final response mismatch").into());
        }
        let force_tol = self.controls.contact.force_residual_tolerance_n;
        let mut gap_bounds_m = Vec::with_capacity(self.sites.len());
        let mut relative_normal_speeds_after = Vec::with_capacity(self.sites.len());
        let mut normal_work = Dd::ZERO;
        for (i, &(g0, u0)) in before.iter().enumerate() {
            poll(cancelled)?;
            let (g1, u1) = self.gap_speed(i, &p.state, &n.state);
            let f = contact.forces_n[i];
            if f[0] < -force_tol {
                return Err(Source::Invalid("tensile normal contact force").into());
            }
            check("absent first tangential force", f[1], force_tol)?;
            check("absent second tangential force", f[2], force_tol)?;
            let bounds = quadratic_bounds(g0, h * u0, g1)?;
            if bounds[0] < -self.controls.gap_tolerance_m {
                return Err(Source::Invalid("penetrating midpoint gap path").into());
            }
            if f[0] > force_tol {
                check(
                    "loaded contact stays closed",
                    bounds[0].abs().max(bounds[1].abs()),
                    self.controls.gap_tolerance_m,
                )?;
                check(
                    "loaded contact endpoint normal speed",
                    u1,
                    self.controls.mechanics.speed_tolerance_m_s,
                )?;
            }
            normal_work = normal_work + dd(f[0]) * (dd(g1) - dd(g0));
            gap_bounds_m.push(bounds);
            relative_normal_speeds_after.push(u1);
        }
        let pa = self.positive.audit(
            positive,
            &p,
            h,
            positive_loads.midpoint,
            self.controls,
            cancelled,
        )?;
        let na = self.negative.audit(
            negative,
            &n,
            h,
            negative_loads.midpoint,
            self.controls,
            cancelled,
        )?;
        let mut momentum_balance_residual = [0.0; 3];
        for (i, out) in momentum_balance_residual.iter_mut().enumerate() {
            let value = pa.momentum[i] + na.momentum[i] + pa.damping[i] + na.damping[i]
                - pa.external[i]
                - na.external[i]
                - pa.support[i]
                - na.support[i];
            let allowance = if i == 0 {
                self.controls.mechanics.linear_impulse_tolerance_ns
            } else {
                self.controls.mechanics.angular_impulse_tolerance_nms
            };
            *out = value.to_f64();
            check(
                "combined world momentum with damping exchange",
                *out,
                allowance,
            )?;
        }
        let contact_work = pa.contact_work + na.contact_work;
        check(
            "reciprocal point work",
            (contact_work - normal_work).to_f64(),
            self.controls.mechanics.energy_tolerance_j,
        )?;
        check(
            "ideal sustained normal work",
            contact_work.to_f64(),
            self.controls.mechanics.energy_tolerance_j,
        )?;
        let energy_residual_j = (pa.energy_residual + na.energy_residual + contact_work).to_f64();
        check(
            "combined midpoint energy",
            energy_residual_j,
            self.controls.mechanics.energy_tolerance_j,
        )?;
        // Carry a physically consistent endpoint acceleration into the next
        // native step even when the external source changes during the interval.
        self.positive.restart(
            &mut p,
            positive_loads.end,
            self.controls,
            self.budget.max_workspace_components,
            cancelled,
        )?;
        self.negative.restart(
            &mut n,
            negative_loads.end,
            self.controls,
            self.budget.max_workspace_components,
            cancelled,
        )?;
        positive
            .history
            .try_reserve(1)
            .map_err(|_| Source::Budget("positive history allocation"))?;
        negative
            .history
            .try_reserve(1)
            .map_err(|_| Source::Budget("negative history allocation"))?;
        poll(cancelled)?;
        let report = PlateSustainedReport {
            contact,
            gap_bounds_m,
            relative_normal_speeds_after,
            positive: pa.report(),
            negative: na.report(),
            momentum_balance_residual,
            energy_residual_j,
        };
        publish(positive, p.state);
        publish(negative, n.state);
        Ok(report)
    }

    fn gap_speed(
        &self,
        i: usize,
        positive: &SecondOrderState,
        negative: &SecondOrderState,
    ) -> (f64, f64) {
        let p = self.positive.points.row(i).1;
        let n = self.negative.points.row(i).1;
        let gap = dd(self.positive.native.reference_height_m)
            - dd(self.negative.native.reference_height_m)
            + dot(p, &positive.q)
            - dot(n, &negative.q);
        let speed = dot(p, &positive.v) - dot(n, &negative.v);
        (gap.to_f64(), speed.to_f64())
    }
}

struct Candidate {
    state: SecondOrderState,
    contact_force: Vec<Dd>,
}

impl<'a, L: PlateLoad + ?Sized> Surface<'a, L> {
    fn prepare(
        native: PlateSustainedSurface<'a, L>,
        sites: &[PlateImpactSite],
        controls: PlateSustainedControls,
        budget: PlateRigidImpactSetBudget,
        cancelled: &mut impl FnMut() -> bool,
    ) -> SourceResult<Self> {
        let points = PlateRigidImpactSet::prepare(
            native.dynamics,
            native.mesh,
            native.full_mass,
            sites,
            controls.mechanics,
            budget,
            cancelled,
        )
        .map_err(Source::Preparation)?;
        admit_matrix(
            native.full_mass,
            &native.dynamics.model().m,
            &native.dynamics.model().dof_map,
            cancelled,
        )?;
        admit_matrix(
            native.full_stiffness,
            &native.dynamics.model().k,
            &native.dynamics.model().dof_map,
            cancelled,
        )?;
        Ok(Self { native, points })
    }

    fn admit_state(&self, state: &SecondOrderState, maximum_history: usize) -> SourceResult<()> {
        let n = self.native.dynamics.model().free;
        if !state.t.is_finite()
            || [&state.q, &state.v, &state.a]
                .iter()
                .any(|x| x.len() != n || x.iter().any(|v| !v.is_finite()))
        {
            return Err(Source::Invalid("finite native state dimensions"));
        }
        if state.steps.checked_add(1).is_none()
            || state
                .history
                .len()
                .checked_add(1)
                .is_none_or(|n| n > maximum_history)
        {
            return Err(Source::Budget("native accepted history or step counter"));
        }
        Ok(())
    }

    fn admit_loads(
        &self,
        t: f64,
        h: f64,
        loads: PlateSustainedLoads<'_>,
        cancelled: &mut impl FnMut() -> bool,
    ) -> SourceResult<()> {
        let map = &self.native.dynamics.model().dof_map;
        let mut sample = vec![f64::NAN; self.native.dynamics.model().free];
        for (time, full) in [
            (t, loads.start),
            (t + 0.5 * h, loads.midpoint),
            (t + h, loads.end),
        ] {
            poll(cancelled)?;
            if full.len() != map.len() || full.iter().any(|x| !x.is_finite()) {
                return Err(Source::Invalid(
                    "complete finite generalized source samples",
                ));
            }
            sample.fill(f64::NAN);
            self.native
                .dynamics
                .forcing(time, &mut sample)
                .map_err(Source::Load)?;
            for (i, &r) in map.iter().enumerate() {
                if let Some(r) = r {
                    if full[i] != sample[r] {
                        return Err(Source::Invalid("source callback/sample mismatch"));
                    }
                }
            }
        }
        poll(cancelled)
    }

    #[allow(clippy::too_many_arguments)]
    fn response(
        &self,
        before: &SecondOrderState,
        h: f64,
        loads: PlateSustainedLoads<'_>,
        forces: &[[f64; 3]],
        sign: f64,
        controls: PlateSustainedControls,
        workspace: usize,
        cancelled: &mut impl FnMut() -> bool,
    ) -> SourceResult<Candidate> {
        poll(cancelled)?;
        let mut contact_force = vec![Dd::ZERO; self.native.dynamics.model().dof_map.len()];
        for (i, force) in forces.iter().enumerate() {
            poll(cancelled)?;
            for (j, &value) in self.points.row(i).2.iter().enumerate() {
                if j % 256 == 0 {
                    poll(cancelled)?;
                }
                contact_force[j] = contact_force[j] + dd(value) * dd(force[0]) * dd(sign);
            }
        }
        let start = self.total_force(loads.start, &contact_force)?;
        let acceleration = second_order_acceleration(
            &self.native.dynamics,
            &before.q,
            &before.v,
            &start,
            &self.native.dynamics.mass_preconditioner(),
            controls.mechanics.impulse.mass_solve,
            workspace,
            cancelled,
        )
        .map_err(Source::Initialization)?;
        let mut state = SecondOrderState::new(before.t, &before.q, &before.v, &acceleration.value);
        state.steps = before.steps;
        let midpoint = self.total_force(loads.midpoint, &contact_force)?;
        OperatorGeneralizedAlpha::new(before.q.len(), h, 1.0, controls.step)
            .step_controlled(&mut state, &self.native.dynamics, &midpoint, cancelled)
            .map_err(Source::Step)?;
        Ok(Candidate {
            state,
            contact_force,
        })
    }

    fn total_force(&self, source: &[f64], contact: &[Dd]) -> SourceResult<Vec<f64>> {
        let mut force = vec![0.0; self.native.dynamics.model().free];
        for (i, &r) in self.native.dynamics.model().dof_map.iter().enumerate() {
            if let Some(r) = r {
                force[r] = (dd(source[i]) + contact[i]).to_f64();
                if !force[r].is_finite() {
                    return Err(Source::Invalid("unrepresentable native load"));
                }
            }
        }
        Ok(force)
    }

    fn restart(
        &self,
        candidate: &mut Candidate,
        end: &[f64],
        controls: PlateSustainedControls,
        workspace: usize,
        cancelled: &mut impl FnMut() -> bool,
    ) -> SourceResult<()> {
        let end = self.total_force(end, &candidate.contact_force)?;
        candidate.state.a = second_order_acceleration(
            &self.native.dynamics,
            &candidate.state.q,
            &candidate.state.v,
            &end,
            &self.native.dynamics.mass_preconditioner(),
            controls.mechanics.impulse.mass_solve,
            workspace,
            cancelled,
        )
        .map_err(Source::Initialization)?
        .value;
        Ok(())
    }

    fn audit(
        &self,
        before: &SecondOrderState,
        candidate: &Candidate,
        h: f64,
        source: &[f64],
        controls: PlateSustainedControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> SourceResult<Audit> {
        let map = &self.native.dynamics.model().dof_map;
        let after = &candidate.state;
        let mut dq = vec![Dd::ZERO; map.len()];
        let mut dv = dq.clone();
        let mut qm = dq.clone();
        let mut vm = dq.clone();
        for (i, &r) in map.iter().enumerate() {
            if i % 256 == 0 {
                poll(cancelled)?;
            }
            if let Some(r) = r {
                dq[i] = dd(after.q[r]) - dd(before.q[r]);
                dv[i] = dd(after.v[r]) - dd(before.v[r]);
                qm[i] = dd(0.5 * before.q[r]) + dd(0.5 * after.q[r]);
                vm[i] = dd(0.5 * before.v[r]) + dd(0.5 * after.v[r]);
                check(
                    "actual native midpoint kinematics",
                    (dq[i] - dd(h) * vm[i]).to_f64(),
                    controls.kinematic_tolerance,
                )?;
            }
        }
        let params = self.native.dynamics.parameters();
        let mut audit = Audit::new(map.len() - before.q.len());
        let mut elastic = [Dd::ZERO; 3];
        for (i, &r) in map.iter().enumerate() {
            if i % 256 == 0 {
                poll(cancelled)?;
            }
            let dp = dd(params.mass_scale) * row_product(self.native.full_mass, i, &dv);
            let km = dd(params.stiffness_scale) * row_product(self.native.full_stiffness, i, &qm);
            let dk = dd(params.stiffness_scale) * row_product(self.native.full_stiffness, i, &dq);
            let cv = dd(params.mass_scale)
                * dd(params.mass_damping_per_s)
                * row_product(self.native.full_mass, i, &vm)
                + dd(params.stiffness_scale)
                    * dd(params.stiffness_damping_s)
                    * row_product(self.native.full_stiffness, i, &vm);
            let ji = dd(h) * km;
            let jd = dd(h) * cv;
            let je = dd(h) * dd(source[i]);
            let jc = dd(h) * candidate.contact_force[i];
            let residual = dp + ji + jd - je - jc;
            if [dp, ji, jd, je, jc, residual]
                .iter()
                .any(|x| !x.to_f64().is_finite())
            {
                return Err(Source::Invalid("unrepresentable full-row midpoint balance"));
            }
            if r.is_some() {
                let scale = [dp, ji, jd, je, jc]
                    .iter()
                    .map(|x| x.to_f64().abs())
                    .fold(0.0, f64::max);
                check(
                    "actual native midpoint force balance",
                    residual.to_f64(),
                    controls.mechanics.impulse.absolute_impulse_tolerance
                        + controls.mechanics.impulse.relative_impulse_tolerance * scale,
                )?;
            } else {
                audit.support_impulses.push((i, residual.to_f64()));
                add_resultant(self.native.mesh, i, residual, &mut audit.support);
            }
            add_resultant(self.native.mesh, i, dp, &mut audit.momentum);
            add_resultant(self.native.mesh, i, je, &mut audit.external);
            add_resultant(self.native.mesh, i, jd, &mut audit.damping);
            add_resultant(self.native.mesh, i, ji, &mut elastic);
            audit.kinetic = audit.kinetic + vm[i] * dp;
            audit.elastic = audit.elastic + qm[i] * dk;
            audit.external_work = audit.external_work + dd(source[i]) * dq[i];
            audit.contact_work = audit.contact_work + candidate.contact_force[i] * dq[i];
            audit.damping_work = audit.damping_work + vm[i] * jd;
        }
        // Stiffness supplies no net rigid transverse force or moment in this
        // admitted step. Foundation springs require a distinct external owner.
        for (i, value) in elastic.iter().enumerate() {
            check(
                "internal elastic resultant",
                value.to_f64(),
                if i == 0 {
                    controls.mechanics.linear_impulse_tolerance_ns
                } else {
                    controls.mechanics.angular_impulse_tolerance_nms
                },
            )?;
        }
        if !audit.damping_work.to_f64().is_finite()
            || audit.damping_work.to_f64() < -controls.mechanics.energy_tolerance_j
        {
            return Err(Source::Invalid("nonpassive midpoint Rayleigh damping"));
        }
        audit.energy_residual = audit.kinetic + audit.elastic + audit.damping_work
            - audit.external_work
            - audit.contact_work;
        check(
            "individual native midpoint work",
            audit.energy_residual.to_f64(),
            controls.mechanics.impulse.work_tolerance,
        )?;
        Ok(audit)
    }
}

struct Audit {
    support_impulses: Vec<(usize, f64)>,
    momentum: [Dd; 3],
    external: [Dd; 3],
    damping: [Dd; 3],
    support: [Dd; 3],
    kinetic: Dd,
    elastic: Dd,
    external_work: Dd,
    damping_work: Dd,
    contact_work: Dd,
    energy_residual: Dd,
}
impl Audit {
    fn new(fixed: usize) -> Self {
        Self {
            support_impulses: Vec::with_capacity(fixed),
            momentum: [Dd::ZERO; 3],
            external: [Dd::ZERO; 3],
            damping: [Dd::ZERO; 3],
            support: [Dd::ZERO; 3],
            kinetic: Dd::ZERO,
            elastic: Dd::ZERO,
            external_work: Dd::ZERO,
            damping_work: Dd::ZERO,
            contact_work: Dd::ZERO,
            energy_residual: Dd::ZERO,
        }
    }
    fn report(self) -> PlateSustainedBalance {
        PlateSustainedBalance {
            support_impulses: self.support_impulses,
            momentum_change: self.momentum.map(Dd::to_f64),
            external_impulse: self.external.map(Dd::to_f64),
            damping_impulse: self.damping.map(Dd::to_f64),
            support_impulse: self.support.map(Dd::to_f64),
            kinetic_energy_change_j: self.kinetic.to_f64(),
            elastic_energy_change_j: self.elastic.to_f64(),
            external_work_midpoint_j: self.external_work.to_f64(),
            damping_midpoint_j: self.damping_work.to_f64(),
            contact_work_j: self.contact_work.to_f64(),
            energy_residual_j: self.energy_residual.to_f64(),
        }
    }
}

fn admit_matrix(
    full: &Csr,
    reduced: &Csr,
    map: &[Option<usize>],
    cancelled: &mut impl FnMut() -> bool,
) -> SourceResult<()> {
    if full.nrows() != map.len() || full.ncols() != map.len() {
        return Err(Source::Invalid("complete matrix dimensions"));
    }
    let mut inverse = vec![0; reduced.nrows()];
    for (i, &r) in map.iter().enumerate() {
        poll(cancelled)?;
        if let Some(r) = r {
            inverse[r] = i;
        }
        let (cols, values) = full.row(i);
        for (&j, &value) in cols.iter().zip(values) {
            if !value.is_finite() || entry(full, j, i) != value {
                return Err(Source::Invalid(
                    "finite symmetric original complete matrix required",
                ));
            }
            if let (Some(r), Some(c)) = (r, map[j]) {
                if entry(reduced, r, c) != value {
                    return Err(Source::Invalid("complete/reduced matrix mismatch"));
                }
            }
        }
    }
    for r in 0..reduced.nrows() {
        poll(cancelled)?;
        let (cols, values) = reduced.row(r);
        for (&c, &value) in cols.iter().zip(values) {
            if entry(full, inverse[r], inverse[c]) != value {
                return Err(Source::Invalid("complete matrix omits retained entry"));
            }
        }
    }
    Ok(())
}
fn entry(matrix: &Csr, row: usize, col: usize) -> f64 {
    let (columns, values) = matrix.row(row);
    columns.binary_search(&col).map_or(0.0, |i| values[i])
}
fn row_product(matrix: &Csr, row: usize, vector: &[Dd]) -> Dd {
    let (columns, values) = matrix.row(row);
    columns
        .iter()
        .zip(values)
        .fold(Dd::ZERO, |sum, (&col, &value)| {
            sum + dd(value) * vector[col]
        })
}
fn dot(row: &[f64], vector: &[f64]) -> Dd {
    row.iter()
        .zip(vector)
        .fold(Dd::ZERO, |sum, (&a, &b)| sum + dd(a) * dd(b))
}
fn dd(value: f64) -> Dd {
    Dd::from_f64(value)
}
fn add_resultant(mesh: &PlateMesh, row: usize, value: Dd, result: &mut [Dd; 3]) {
    let (x, y) = mesh.nodes[row / 3];
    match row % 3 {
        0 => {
            result[0] = result[0] + value;
            result[1] = result[1] + dd(y) * value;
            result[2] = result[2] - dd(x) * value;
        }
        1 => result[2] = result[2] - value,
        _ => result[1] = result[1] + value,
    }
}
fn quadratic_bounds(g0: f64, tangent: f64, g1: f64) -> SourceResult<[f64; 2]> {
    let mut lo = g0.min(g1);
    let mut hi = g0.max(g1);
    let curvature = dd(g1) - dd(g0) - dd(tangent);
    if curvature.to_f64() != 0.0 {
        let theta = (-dd(tangent) / (dd(2.0) * curvature)).to_f64();
        if theta > 0.0 && theta < 1.0 {
            let value = (dd(g0) + dd(theta) * (dd(tangent) + dd(theta) * curvature)).to_f64();
            if !value.is_finite() {
                return Err(Source::Invalid("unrepresentable gap interpolant"));
            }
            lo = lo.min(value);
            hi = hi.max(value);
        }
    }
    if !g0.is_finite() || !g1.is_finite() || !tangent.is_finite() || !curvature.to_f64().is_finite()
    {
        return Err(Source::Invalid("unrepresentable gap interpolant"));
    }
    Ok([lo, hi])
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> SourceResult<()> {
    if cancelled() {
        Err(Source::Cancelled)
    } else {
        Ok(())
    }
}
fn check(what: &'static str, residual: f64, allowance: f64) -> SourceResult<()> {
    if !residual.is_finite()
        || !allowance.is_finite()
        || allowance < 0.0
        || residual.abs() > allowance
    {
        Err(Source::Balance {
            what,
            residual,
            allowance,
        })
    } else {
        Ok(())
    }
}
fn publish(target: &mut SecondOrderState, mut candidate: SecondOrderState) {
    target.t = candidate.t;
    target.q = candidate.q;
    target.v = candidate.v;
    target.a = candidate.a;
    target.steps = candidate.steps;
    target.history.append(&mut candidate.history);
}
