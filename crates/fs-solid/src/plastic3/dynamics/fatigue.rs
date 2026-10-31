//! Source-bound strain life on the original physical-time inelastic solve.
//!
//! Each accepted interval calls `advance_backward_euler` with the complete
//! previous displacement, velocity and J2/Prager history. Every actual native
//! stress and total mechanical strain sample enters the same fixed-plane
//! observer used by the quasistatic J2 path. Initial velocity can therefore
//! produce fatigue with no applied force; load factors are not fatigue samples.
//!
//! The source-resolved reference density, constitutive coefficients, external
//! conditions, mesh and support DOFs stay fixed. Absolute endpoint load factors
//! scale the original prescribed displacement and dead-force patterns. This
//! adapter supplies no thermal eigenstrain, evolving-temperature law, arbitrary
//! imported state, nonproportional correction or fatigue-to-stiffness feedback.
//!
//! Backward Euler has numerical damping. Its original inertial energy receipt
//! remains distinct from physical plastic dissipation and Miner damage. The
//! caller must refine physical time and resolve all relevant turning points;
//! exact rainflow counting of accepted samples cannot recover missed extrema.

use std::sync::Arc;

use fs_blake3::ContentHash;
use fs_exec::Cx;
use fs_material::fatigue::strain_life::{StrainLifeCurve, StrainLifeSummary};
use fs_material::plastic::J2State;
use fs_material::plastic::source::ResolvedJ2Plasticity;
use fs_material::state_point::DENSITY_PROPERTY;
use fs_qty::{Density, QuantitySpec};

use super::InertialEnergy;
use crate::plastic3::fatigue::{
    J2FatigueConfig, J2FatigueError, J2FatigueGauge, J2FatigueObserver, J2FatigueSample, poll,
    with_path,
};
use crate::plastic3::{SmallStrainTetProblem, SmallStrainTetSolution};

/// Whole-run work bounds. The history budget includes the real initial sample
/// and every physical interval, regardless of how calls are divided into blocks.
#[derive(Debug, Clone, Copy)]
pub struct J2DynamicFatigueConfig {
    pub fatigue: J2FatigueConfig,
    pub maximum_steps: usize,
}

/// One actual physical endpoint. Time must strictly increase; the load factor
/// is absolute against the original support and dead-force patterns.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct J2DynamicFatigueTarget {
    pub time_s: f64,
    pub load_factor: f64,
}

/// Latest accepted physical interval and the ORIGINAL native inertia receipt.
/// Full native stresses, plastic state and reactions remain in `solution()`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct J2DynamicFatigueStep {
    pub index: usize,
    pub time_s: f64,
    pub dt_s: f64,
    pub load_factor: f64,
    pub energy: InertialEnergy,
}

/// Owning, resumable native motion and bounded strain-life observation.
/// Clones preserve the complete accepted mechanics and unresolved cycles.
/// Only the latest native solution/energy receipt is retained; source density
/// is shared immutably and no external stress/strain arrays can be injected.
#[derive(Debug, Clone)]
pub struct J2DynamicFatigueRun<'p, 'm> {
    problem: &'p SmallStrainTetProblem<'m, ResolvedJ2Plasticity>,
    solution: SmallStrainTetSolution<J2State>,
    velocity_mps: Vec<[f64; 3]>,
    reference_density_kg_m3: Arc<[f64]>,
    observer: J2FatigueObserver,
    time_s: f64,
    maximum_steps: usize,
    accepted_steps: usize,
    last_step: Option<J2DynamicFatigueStep>,
}

impl<'p, 'm> J2DynamicFatigueRun<'p, 'm> {
    /// Admit the real virgin zero-load native solution as sample zero, with an
    /// explicit finite initial nodal velocity and nonnegative time origin.
    /// Prescribed DOFs initially have zero displacement AND zero velocity.
    /// The template must have `load_path == [1.0]`; its support/force values are
    /// patterns for future absolute targets, not an ignored quasistatic path.
    /// Every original J2 material bundle must also contain a positive, static,
    /// source-resolved SI density at the same frozen external conditions.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        problem: &'p SmallStrainTetProblem<'m, ResolvedJ2Plasticity>,
        mesh_frame: ContentHash,
        gauges: Vec<J2FatigueGauge>,
        config: J2DynamicFatigueConfig,
        initial_time_s: f64,
        initial_velocity_mps: Vec<[f64; 3]>,
        cx: &Cx<'_>,
    ) -> Result<Self, J2FatigueError> {
        poll(cx)?;
        if problem.load_path != [1.0]
            || !initial_time_s.is_finite()
            || initial_time_s < 0.0
            || initial_velocity_mps.len() != problem.nodes_m.len()
            || initial_velocity_mps
                .iter()
                .flatten()
                .any(|v| !v.is_finite())
        {
            return Err(J2FatigueError::Invalid(
                "one unit endpoint pattern, nonnegative finite time and finite full nodal velocity required",
            ));
        }
        if config.maximum_steps == 0
            || u64::try_from(config.maximum_steps)
                .ok()
                .and_then(|n| n.checked_add(1))
                .is_none_or(|n| n > config.fatigue.history.maximum_samples)
        {
            return Err(J2FatigueError::Budget);
        }
        let mut observer = J2FatigueObserver::new(problem, mesh_frame, gauges, config.fatigue, cx)?;
        let reference_density_kg_m3 = source_density(problem, cx)?;
        // This original solve validates the complete mesh, supports, law field
        // and native budgets before any constrained velocity is indexed.
        let solution = with_path(problem, &[0.0]).solve(cx)?;
        if problem
            .prescribed_m
            .iter()
            .any(|&(dof, _)| initial_velocity_mps.as_flattened()[dof] != 0.0)
        {
            return Err(J2FatigueError::Invalid(
                "initial velocity must be zero on the original prescribed DOFs",
            ));
        }
        observer.observe(problem, &solution, cx)?;
        poll(cx)?;
        Ok(Self {
            problem,
            solution,
            velocity_mps: initial_velocity_mps,
            reference_density_kg_m3: reference_density_kg_m3.into(),
            observer,
            time_s: initial_time_s,
            maximum_steps: config.maximum_steps,
            accepted_steps: 0,
            last_step: None,
        })
    }

    pub fn mesh_frame(&self) -> ContentHash {
        self.observer.mesh_frame()
    }
    pub fn time_s(&self) -> f64 {
        self.time_s
    }
    pub fn accepted_steps(&self) -> usize {
        self.accepted_steps
    }
    /// Remaining work budget, not a requirement to advance before finishing.
    pub fn remaining_steps(&self) -> usize {
        self.maximum_steps - self.accepted_steps
    }
    pub fn solution(&self) -> &SmallStrainTetSolution<J2State> {
        &self.solution
    }
    pub fn velocity_mps(&self) -> &[[f64; 3]] {
        &self.velocity_mps
    }
    pub fn reference_density_kg_m3(&self) -> &[f64] {
        &self.reference_density_kg_m3
    }
    pub fn last_step(&self) -> Option<&J2DynamicFatigueStep> {
        self.last_step.as_ref()
    }
    pub fn summaries(&self) -> Vec<StrainLifeSummary> {
        self.observer.summaries()
    }
    pub fn samples(&self) -> Vec<J2FatigueSample> {
        self.observer.samples()
    }
    pub fn curve(&self, gauge: usize) -> Option<&StrainLifeCurve> {
        self.observer.curve(gauge)
    }
    /// Complete original constitutive and density source receipts.
    pub fn material(&self, gauge: usize) -> Option<&ResolvedJ2Plasticity> {
        self.observer.material(self.problem, gauge)
    }

    /// Accept one physical interval atomically across mechanics and all gauges.
    pub fn advance_to(
        &mut self,
        time_s: f64,
        load_factor: f64,
        cx: &Cx<'_>,
    ) -> Result<J2DynamicFatigueStep, J2FatigueError> {
        self.advance(
            &[J2DynamicFatigueTarget {
                time_s,
                load_factor,
            }],
            cx,
        )
    }

    /// Accept every endpoint in a bounded block as one transaction. Every
    /// interval calls the original native inertia/Newton owner and every actual
    /// endpoint enters the gauges. Cancellation, nonconvergence, source-domain
    /// or ANY gauge refusal rolls back the entire requested block, including
    /// physical time, velocity, constitutive history and open cycles.
    ///
    /// Returns the latest interval's original energy; earlier interval records
    /// are not retained. Use single-step calls when a full energy trace is needed.
    pub fn advance(
        &mut self,
        targets: &[J2DynamicFatigueTarget],
        cx: &Cx<'_>,
    ) -> Result<J2DynamicFatigueStep, J2FatigueError> {
        self.observer.check_sample(cx)?;
        if targets.is_empty() {
            return Err(J2FatigueError::Invalid(
                "declare at least one physical endpoint",
            ));
        }
        if targets.len() > self.remaining_steps() {
            return Err(J2FatigueError::Budget);
        }
        let mut previous_time = self.time_s;
        for (index, target) in targets.iter().enumerate() {
            if index % 256 == 0 {
                poll(cx)?;
            }
            if !target.time_s.is_finite()
                || target.time_s <= previous_time
                || !target.load_factor.is_finite()
                || !(target.time_s - previous_time).is_finite()
            {
                return Err(J2FatigueError::Invalid(
                    "physical endpoints need strictly increasing finite time and finite absolute load",
                ));
            }
            previous_time = target.time_s;
        }
        let mut solution = self.solution.clone();
        let mut velocity = self.velocity_mps.clone();
        let mut observer = self.observer.clone();
        let mut prescribed = self.problem.prescribed_m.to_vec();
        let mut forces = self.problem.nodal_forces_n.to_vec();
        let mut time_s = self.time_s;
        let mut last_step = None;
        for (offset, target) in targets.iter().enumerate() {
            observer.check_sample(cx)?;
            scale_pattern(
                self.problem.prescribed_m,
                &mut prescribed,
                target.load_factor,
                cx,
            )?;
            scale_pattern(
                self.problem.nodal_forces_n,
                &mut forces,
                target.load_factor,
                cx,
            )?;
            let endpoint = SmallStrainTetProblem {
                prescribed_m: &prescribed,
                nodal_forces_n: &forces,
                ..with_path(self.problem, &[1.0])
            };
            let dt_s = target.time_s - time_s;
            let motion = endpoint.advance_backward_euler(
                solution.displacement_m.as_flattened(),
                &velocity,
                &solution.element_states,
                &self.reference_density_kg_m3,
                dt_s,
                cx,
            )?;
            observer.observe(self.problem, &motion.solution, cx)?;
            last_step = Some(J2DynamicFatigueStep {
                index: self.accepted_steps + offset + 1,
                time_s: target.time_s,
                dt_s,
                load_factor: target.load_factor,
                energy: motion.energy,
            });
            solution = motion.solution;
            velocity = motion.velocity_mps;
            time_s = target.time_s;
        }
        let report = last_step.ok_or(J2FatigueError::Invalid("missing physical endpoint"))?;
        poll(cx)?;
        self.solution = solution;
        self.velocity_mps = velocity;
        self.observer = observer;
        self.time_s = time_s;
        self.accepted_steps += targets.len();
        self.last_step = Some(report);
        Ok(report)
    }

    /// Finalize all endpoint half cycles together; a failed finish changes
    /// nothing. Successful finish is idempotent and prohibits later motion.
    pub fn finish(&mut self, cx: &Cx<'_>) -> Result<Vec<StrainLifeSummary>, J2FatigueError> {
        self.observer.finish(cx)
    }
}

fn source_density(
    problem: &SmallStrainTetProblem<'_, ResolvedJ2Plasticity>,
    cx: &Cx<'_>,
) -> Result<Vec<f64>, J2FatigueError> {
    let mut densities = Vec::with_capacity(problem.tetrahedra.len());
    for element in 0..problem.tetrahedra.len() {
        if element % 256 == 0 {
            poll(cx)?;
        }
        let density = problem
            .laws
            .at(element)
            .material()
            .property(DENSITY_PROPERTY)
            .ok_or(J2FatigueError::Invalid(
                "every original J2 material needs one unambiguous source density",
            ))?;
        if density.requirement().quantity() != QuantitySpec::dimensional(Density::DIMS)
            || !density.is_static_at_fixed_external_conditions()
            || !density.value_si().is_finite()
            || density.value_si() <= 0.0
        {
            return Err(J2FatigueError::Invalid(
                "reference density must be positive SI mass density at frozen external conditions without mechanical/frequency dependence",
            ));
        }
        densities.push(density.value_si());
    }
    poll(cx)?;
    Ok(densities)
}

fn scale_pattern(
    original: &[(usize, f64)],
    target: &mut [(usize, f64)],
    factor: f64,
    cx: &Cx<'_>,
) -> Result<(), J2FatigueError> {
    for (index, (&(_, value), (_, endpoint))) in original.iter().zip(target).enumerate() {
        if index % 256 == 0 {
            poll(cx)?;
        }
        *endpoint = value * factor;
        if !endpoint.is_finite() {
            return Err(J2FatigueError::Invalid(
                "absolute endpoint load is not representable",
            ));
        }
    }
    Ok(())
}
