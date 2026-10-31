//! Bounded strain-life observation on the native source-bound J2 load path.
//!
//! Each accepted load factor uses the original `solve_from` owner and retains
//! its plastic strain and isotropic/kinematic history. Gauges observe signed
//! TOTAL MECHANICAL strain in a fixed reference plane, not accumulated plastic
//! strain or unsigned equivalent stress. `ResolvedJ2Plasticity` has no thermal
//! eigenstrain: its source temperature and coefficients stay fixed throughout.
//! Thermal-wrapped laws and changing-temperature paths need their own admitted
//! mechanical-strain transfer and cannot enter this fixed-state adapter.
//!
//! A Coffin-Manson-Basquin calibration is a caller declaration. Matching the
//! original material state, process, modulus and temperature prevents accidental
//! rebinding; it does not validate that calibration or predict component life.
//! Miner damage changes neither stiffness, heat nor failure topology. Native
//! small-strain/P1 and linear-hardening limits remain in force. Resolve every
//! physical turning point in `load_path`; no extrema between targets are inferred.

use std::sync::Arc;

use fs_blake3::ContentHash;
use fs_exec::Cx;
use fs_material::fatigue::strain_life::{
    StrainLifeCurve, StrainLifeMeanStress, StrainLifeMiner, StrainLifeSummary,
};
use fs_material::fatigue::{FatigueError, NormalStressDirection, RainflowLimits};
use fs_material::plastic::J2State;
use fs_material::plastic::source::ResolvedJ2Plasticity;

use super::{
    SmallStrainLawField, SmallStrainTetError, SmallStrainTetProblem, SmallStrainTetSolution,
};

/// Fixed-plane observation at one original tetrahedron, with explicit source
/// identity and total mechanical strain-life calibration.
#[derive(Debug, Clone)]
pub struct J2FatigueGauge {
    pub element: usize,
    pub normal_in_mesh_frame: [f64; 3],
    pub material_state_identity: ContentHash,
    pub curve: Arc<StrainLifeCurve>,
    pub mean_stress: StrainLifeMeanStress,
}

/// Logical storage/work limits; transactional proposals need an additional
/// copy of the latest native solution and these bounded gauge histories.
#[derive(Debug, Clone, Copy)]
pub struct J2FatigueConfig {
    pub history: RainflowLimits,
    pub maximum_gauges: usize,
    pub maximum_total_residue: usize,
}

#[derive(Debug, Clone, PartialEq)]
pub enum J2FatigueError {
    Invalid(&'static str),
    Budget,
    SourceMismatch { gauge: usize, what: &'static str },
    Gauge { gauge: usize, error: FatigueError },
    Mechanics(SmallStrainTetError),
    Cancelled,
    Finished,
    PathIncomplete,
}
impl core::fmt::Display for J2FatigueError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "cyclic J2 strain-life path refused: {self:?}")
    }
}
impl std::error::Error for J2FatigueError {}
impl From<SmallStrainTetError> for J2FatigueError {
    fn from(error: SmallStrainTetError) -> Self {
        match error {
            SmallStrainTetError::Cancelled => Self::Cancelled,
            error => Self::Mechanics(error),
        }
    }
}

/// Actual accepted plane observations. Stress is diagnostic; the declared
/// uncorrected strain-life model does not infer a mean-stress correction.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct J2FatigueSample {
    pub mechanical_strain: f64,
    pub normal_stress_pa: f64,
}

#[derive(Debug, Clone)]
struct GaugeState {
    element: usize,
    normal: NormalStressDirection,
    counter: StrainLifeMiner,
    sample: J2FatigueSample,
}

/// Resumable observation of one immutable native load path. `advance` commits
/// a requested block only after every mechanical increment and gauge accepts.
/// Cloning retains the complete current plastic and unresolved cycle histories.
/// Only the latest native increment record is retained, not the whole trace.
#[derive(Debug, Clone)]
pub struct J2FatiguePath<'p, 'm> {
    problem: &'p SmallStrainTetProblem<'m, ResolvedJ2Plasticity>,
    mesh_frame: ContentHash,
    solution: SmallStrainTetSolution<J2State>,
    gauges: Vec<GaugeState>,
    accepted_increments: usize,
    finished: bool,
}

impl<'p, 'm> J2FatiguePath<'p, 'm> {
    /// Solve the real virgin, zero-load state and accept it as the first
    /// observation. The declared load path then supplies absolute targets.
    /// Source identities, labels, fixed positive T and modulus must agree with
    /// the actual resolved constitutive field. No thermal strain is subtracted
    /// or guessed, and no pre-existing plastic/fatigue history is fabricated.
    pub fn new(
        problem: &'p SmallStrainTetProblem<'m, ResolvedJ2Plasticity>,
        mesh_frame: ContentHash,
        gauges: Vec<J2FatigueGauge>,
        config: J2FatigueConfig,
        cx: &Cx<'_>,
    ) -> Result<Self, J2FatigueError> {
        poll(cx)?;
        if mesh_frame == ContentHash([0; 32]) || problem.load_path.is_empty() {
            return Err(J2FatigueError::Invalid(
                "declare a nonzero mesh frame and finite nonempty absolute load path",
            ));
        }
        let samples = u64::try_from(problem.load_path.len())
            .ok()
            .and_then(|n| n.checked_add(1));
        let residue = gauges.len().checked_mul(config.history.maximum_residue);
        if gauges.is_empty()
            || gauges.len() > config.maximum_gauges
            || samples.is_none_or(|n| n > config.history.maximum_samples)
            || residue.is_none_or(|n| n > config.maximum_total_residue)
        {
            return Err(J2FatigueError::Budget);
        }
        for chunk in problem.load_path.chunks(256) {
            poll(cx)?;
            if chunk.iter().any(|x| !x.is_finite()) {
                return Err(J2FatigueError::Invalid("load factors must be finite"));
            }
        }
        problem.laws.validate_count(problem.tetrahedra.len())?;
        let mut prepared = Vec::with_capacity(gauges.len());
        for (index, gauge) in gauges.into_iter().enumerate() {
            poll(cx)?;
            if gauge.element >= problem.tetrahedra.len() {
                return Err(J2FatigueError::Invalid(
                    "gauge is outside the original tetrahedral field",
                ));
            }
            let law = problem.laws.at(gauge.element);
            let material = law.material();
            let scope = gauge.curve.scope();
            if gauge.material_state_identity != material.identity()
                || scope.material != material.material().chemistry
                || scope.process != material.material().process
            {
                return Err(J2FatigueError::SourceMismatch {
                    gauge: index,
                    what: "strain-life identity/material/process differs from the original J2 state",
                });
            }
            if material.state_coordinate("T") != Some(scope.temperature_k)
                || gauge.curve.parameters().youngs_modulus_pa != law.law().elastic.youngs
            {
                return Err(J2FatigueError::SourceMismatch {
                    gauge: index,
                    what: "strain-life temperature/modulus differs from the frozen J2 law",
                });
            }
            let normal =
                NormalStressDirection::new(gauge.normal_in_mesh_frame).map_err(|error| {
                    J2FatigueError::Gauge {
                        gauge: index,
                        error,
                    }
                })?;
            let counter = StrainLifeMiner::new(gauge.curve, gauge.mean_stress, config.history)
                .map_err(|error| J2FatigueError::Gauge {
                    gauge: index,
                    error,
                })?;
            prepared.push(GaugeState {
                element: gauge.element,
                normal,
                counter,
                sample: J2FatigueSample {
                    mechanical_strain: 0.0,
                    normal_stress_pa: 0.0,
                },
            });
        }
        let solution = with_path(problem, &[0.0]).solve(cx)?;
        observe(problem, &solution, &mut prepared, cx)?;
        poll(cx)?;
        Ok(Self {
            problem,
            mesh_frame,
            solution,
            gauges: prepared,
            accepted_increments: 0,
            finished: false,
        })
    }

    pub fn mesh_frame(&self) -> ContentHash {
        self.mesh_frame
    }
    pub fn accepted_increments(&self) -> usize {
        self.accepted_increments
    }
    pub fn remaining_increments(&self) -> usize {
        self.problem.load_path.len() - self.accepted_increments
    }
    pub fn solution(&self) -> &SmallStrainTetSolution<J2State> {
        &self.solution
    }
    pub fn summaries(&self) -> Vec<StrainLifeSummary> {
        self.gauges.iter().map(|g| g.counter.summary()).collect()
    }
    pub fn samples(&self) -> Vec<J2FatigueSample> {
        self.gauges.iter().map(|g| g.sample).collect()
    }
    pub fn curve(&self, gauge: usize) -> Option<&StrainLifeCurve> {
        self.gauges.get(gauge).map(|g| g.counter.curve())
    }
    /// Original complete property receipts and native constitutive metadata.
    pub fn material(&self, gauge: usize) -> Option<&ResolvedJ2Plasticity> {
        self.gauges
            .get(gauge)
            .map(|g| self.problem.laws.at(g.element))
    }

    /// Accept exactly `count` planned increments as one transaction. Every
    /// original load factor is observed; callback boundaries close no cycles.
    /// Nonconvergence, cancellation, strain/calibration refusal or any gauge
    /// failure preserves the complete prior mechanics and fatigue histories.
    pub fn advance(&mut self, count: usize, cx: &Cx<'_>) -> Result<(), J2FatigueError> {
        poll(cx)?;
        if self.finished {
            return Err(J2FatigueError::Finished);
        }
        if count == 0 || count > self.remaining_increments() {
            return Err(J2FatigueError::Invalid(
                "advance needs a nonzero count within the remaining load path",
            ));
        }
        let end = self.accepted_increments + count;
        let mut solution = self.solution.clone();
        let mut gauges = self.gauges.clone();
        for &factor in &self.problem.load_path[self.accepted_increments..end] {
            poll(cx)?;
            solution = with_path(self.problem, &[factor]).solve_from(
                solution.displacement_m.as_flattened(),
                &solution.element_states,
                cx,
            )?;
            observe(self.problem, &solution, &mut gauges, cx)?;
        }
        poll(cx)?;
        self.solution = solution;
        self.gauges = gauges;
        self.accepted_increments = end;
        Ok(())
    }

    /// Finish the entire declared observation once, after all planned targets.
    /// Open endpoint half cycles are finalized together for every gauge.
    pub fn finish(&mut self, cx: &Cx<'_>) -> Result<Vec<StrainLifeSummary>, J2FatigueError> {
        poll(cx)?;
        if self.remaining_increments() != 0 {
            return Err(J2FatigueError::PathIncomplete);
        }
        if self.finished {
            return Ok(self.summaries());
        }
        let mut gauges = self.gauges.clone();
        for (index, gauge) in gauges.iter_mut().enumerate() {
            gauge
                .counter
                .finish_with_cancel(&mut || cx.checkpoint().is_err())
                .map_err(|error| J2FatigueError::Gauge {
                    gauge: index,
                    error,
                })?;
        }
        poll(cx)?;
        self.gauges = gauges;
        self.finished = true;
        Ok(self.summaries())
    }
}

fn poll(cx: &Cx<'_>) -> Result<(), J2FatigueError> {
    cx.checkpoint().map_err(|_| J2FatigueError::Cancelled)
}

fn with_path<'a>(
    problem: &'a SmallStrainTetProblem<'_, ResolvedJ2Plasticity>,
    path: &'a [f64],
) -> SmallStrainTetProblem<'a, ResolvedJ2Plasticity> {
    SmallStrainTetProblem {
        nodes_m: problem.nodes_m,
        tetrahedra: problem.tetrahedra,
        laws: match &problem.laws {
            SmallStrainLawField::Uniform(law) => SmallStrainLawField::Uniform(law),
            SmallStrainLawField::PerElement(laws) => SmallStrainLawField::PerElement(laws),
        },
        prescribed_m: problem.prescribed_m,
        nodal_forces_n: problem.nodal_forces_n,
        load_path: path,
        budget: problem.budget,
        settings: problem.settings,
    }
}

fn observe(
    problem: &SmallStrainTetProblem<'_, ResolvedJ2Plasticity>,
    solution: &SmallStrainTetSolution<J2State>,
    gauges: &mut [GaugeState],
    cx: &Cx<'_>,
) -> Result<(), J2FatigueError> {
    // Enforce the original law's actual componentwise tensor-Voigt envelope at
    // EVERY element, including unobserved ones. The generic native solver
    // deliberately leaves this model-domain obligation to its owning caller.
    for (index, strain) in solution.element_strain.iter().enumerate() {
        if index % 256 == 0 {
            poll(cx)?;
        }
        let limit = problem.laws.at(index).law().elastic.strain_limit.min(0.05);
        let state = &solution.element_states[index];
        if strain.iter().any(|x| !x.is_finite() || x.abs() > limit)
            || solution.element_stress[index]
                .iter()
                .any(|x| !x.is_finite())
            || !state.alpha.is_finite()
            || state.alpha < 0.0
            || state.plastic_strain.iter().any(|x| !x.is_finite())
        {
            return Err(J2FatigueError::Invalid(
                "native response exceeds the admitted J2 small-strain/finite-state domain",
            ));
        }
    }
    for (index, gauge) in gauges.iter_mut().enumerate() {
        let strain = gauge
            .normal
            .project_strain(&solution.element_strain[gauge.element])
            .map_err(|error| J2FatigueError::Gauge {
                gauge: index,
                error,
            })?;
        let stress = gauge
            .normal
            .project(&solution.element_stress[gauge.element])
            .map_err(|error| J2FatigueError::Gauge {
                gauge: index,
                error,
            })?;
        gauge
            .counter
            .push_with_cancel(strain, &mut || cx.checkpoint().is_err())
            .map_err(|error| J2FatigueError::Gauge {
                gauge: index,
                error,
            })?;
        gauge.sample = J2FatigueSample {
            mechanical_strain: strain,
            normal_stress_pa: stress,
        };
    }
    Ok(())
}
