//! Scalar source response at the fields recovered by the original spatial solve.
//!
//! Each bounded iteration assembles the SAME FEEC conductivity operator, solves
//! its original voltage/current traces, and refreshes complete pinned sources at
//! actual |E| and |J|. Applied/evaluated conductivity, current and power must all
//! agree. Iteration error is never heat. This Picard closure need not converge
//! for every positive source curve and proves neither uniqueness nor stability.

use super::*;
use fs_material::conductor::{ConductorError, response::ScalarConductorResponse};
use std::sync::Arc;

pub const MAX_RESPONSE_ITERATIONS: usize = 128;

/// Constitutive gates in addition to the native algebraic/current/work gates.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ConductionResponseControls {
    pub maximum_iterations: usize,
    pub absolute_conductivity_tolerance_s_m: f64,
    pub relative_conductivity_tolerance: f64,
    pub current_density_tolerance_a_m2: f64,
}
impl ConductionResponseControls {
    fn validate(self) -> Result<(), ConductionResponseError> {
        if !(1..=MAX_RESPONSE_ITERATIONS).contains(&self.maximum_iterations)
            || !self.absolute_conductivity_tolerance_s_m.is_finite()
            || self.absolute_conductivity_tolerance_s_m < 0.0
            || !self.relative_conductivity_tolerance.is_finite()
            || !(0.0..1.0).contains(&self.relative_conductivity_tolerance)
            || (self.absolute_conductivity_tolerance_s_m == 0.0
                && self.relative_conductivity_tolerance == 0.0)
            || !positive(self.current_density_tolerance_a_m2)
        {
            return Err(ConductionResponseError::Invalid(
                "invalid bounded spatial response controls",
            ));
        }
        Ok(())
    }
}

#[derive(Debug)]
pub enum ConductionResponseError {
    Field(ConductionError),
    Material(ConductorError),
    Invalid(&'static str),
    NotConverged {
        iterations: usize,
        maximum_conductivity_defect_s_m: f64,
        maximum_current_density_defect_a_m2: f64,
        absolute_power_defect_w: f64,
    },
    Cancelled,
}
impl core::fmt::Display for ConductionResponseError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "spatial conductor response: {self:?}")
    }
}
impl std::error::Error for ConductionResponseError {}
impl From<ConductionError> for ConductionResponseError {
    fn from(error: ConductionError) -> Self {
        match error {
            ConductionError::Cancelled => Self::Cancelled,
            other => Self::Field(other),
        }
    }
}
impl From<ConductorError> for ConductionResponseError {
    fn from(error: ConductorError) -> Self {
        Self::Material(error)
    }
}

/// Same geometry and original sources; no independent circuit or thermal state.
/// A frozen linearization is private. Its actual-amplitude solve is repeated
/// inside the source closure; accepted nonlinear fields cannot be scaled to a
/// different excitation or obtained by silently freezing an operating point.
#[derive(Clone)]
pub struct PlanarConductionResponse {
    native: PlanarConduction,
    sources: Arc<[ScalarConductorResponse]>,
}

#[derive(Clone, Debug, PartialEq)]
pub struct ConductionResponseSample {
    pub triangle: usize,
    /// Complete original material evaluated at the actual reported E and J.
    pub conductor: ResolvedConductorState,
    /// Evaluated minus applied conductivity [S/m].
    pub conductivity_defect_s_m: f64,
    /// Native J minus the source-evaluated J(E) [A/m2], in mesh axes.
    pub current_density_defect_a_m2: [f64; 2],
    /// Native local E.J minus source-evaluated power over this cell [W].
    pub power_defect_w: f64,
}

/// Original native field observations and separate constitutive residuals.
/// Applied sources remain in `candidate.problem().conductors()`; evaluated
/// sources below agree to the declared gates, without post-solve substitution.
#[derive(Clone, Debug)]
pub struct ConductionResponseReport {
    pub voltage_v: Vec<f64>,
    pub electric_field_v_m: Vec<[f64; 2]>,
    pub current_density_a_m2: Vec<[f64; 2]>,
    pub cell_power_w: Vec<f64>,
    pub power_w: f64,
    pub nodal_boundary_current_a: Vec<Option<f64>>,
    pub nodal_applied_current_a: Vec<Option<f64>>,
    pub boundary_current_a: f64,
    pub maximum_free_current_residual_a: f64,
    pub power_balance_residual_w: f64,
    pub linear: Option<SolveReport>,
    pub evaluated: Vec<ConductionResponseSample>,
    /// Every native spatial solve, including discarded constitutive trials.
    pub field_solves: usize,
    /// Complete actual-field source queries, excluding initial binding admission.
    pub source_evaluations: usize,
}

pub struct ConductionResponseCandidate {
    problem: PlanarConductionResponse,
    report: ConductionResponseReport,
}
impl ConductionResponseCandidate {
    pub fn problem(&self) -> &PlanarConductionResponse {
        &self.problem
    }
    pub fn report(&self) -> &ConductionResponseReport {
        &self.report
    }
    pub fn into_parts(self) -> (PlanarConductionResponse, ConductionResponseReport) {
        (self.problem, self.report)
    }
}

#[derive(Clone, Copy)]
enum ResponseDrive<'a> {
    Traces {
        voltage: &'a [(usize, f64)],
        currents: &'a [ConductionCurrentTrace],
    },
    Current {
        positive: &'a [usize],
        negative: &'a [usize],
        current_a: f64,
    },
}

impl PlanarConductionResponse {
    /// Bind one source law per original triangle. Every source is a positive
    /// memoryless scalar conductor with explicitly owned magnitude coordinates.
    /// Scalar fixed-state and nonlinear regions may coexist. Ordinary tensor,
    /// field-storage and moving-conductor paths retain their own admission.
    pub fn new(
        mesh: TriComplex2,
        sources: Vec<ScalarConductorResponse>,
        limits: WeightedAssemblyLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, ConductionResponseError> {
        poll(cancelled)?;
        if sources.len() != mesh.faces().len() || sources.len() > limits.max_cells {
            return Err(ConductionResponseError::Invalid(
                "one bounded source response per triangle",
            ));
        }
        let mut conductors = Vec::with_capacity(sources.len());
        for source in &sources {
            poll(cancelled)?;
            source.admit(source.original())?;
            conductors.push(source.original().clone());
        }
        let native = PlanarConduction::assemble(mesh, conductors, limits, cancelled, |_| true)?;
        Ok(Self {
            native,
            sources: sources.into(),
        })
    }

    pub fn mesh(&self) -> &TriComplex2 {
        self.native.mesh()
    }
    pub fn conductors(&self) -> &[ResolvedConductorState] {
        self.native.conductors()
    }
    pub fn sources(&self) -> &[ScalarConductorResponse] {
        &self.sources
    }

    /// Solve the original voltage/current traces at this fixed external context.
    /// Prescribed current density uses actual native boundary loads, never a
    /// scaled unit solution. Equipotential total-current electrodes use
    /// [`Self::solve_current_candidate`] and its actual terminal constraint.
    #[allow(clippy::too_many_arguments)]
    pub fn solve_candidate(
        &self,
        voltage: &[(usize, f64)],
        currents: &[ConductionCurrentTrace],
        assembly: WeightedAssemblyLimits,
        field: ConductionControls,
        controls: ConductionResponseControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<ConductionResponseCandidate, ConductionResponseError> {
        self.solve_with_materials_candidate(
            self.conductors().to_vec(),
            voltage,
            currents,
            assembly,
            field,
            controls,
            cancelled,
        )
    }

    /// Close E/J at an explicitly refreshed original material context. An outer
    /// finite-heat owner supplies each proposed T; all other source roles/pins,
    /// geometry and basis stay original. Every rejected trial leaves this owner
    /// unchanged. Actual E.J and terminal work are authoritative; nonlinear
    /// coefficient/current/power defects remain separate numerical residuals.
    #[allow(clippy::too_many_arguments)]
    pub fn solve_with_materials_candidate(
        &self,
        seeds: Vec<ResolvedConductorState>,
        voltage: &[(usize, f64)],
        currents: &[ConductionCurrentTrace],
        assembly: WeightedAssemblyLimits,
        field_controls: ConductionControls,
        controls: ConductionResponseControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<ConductionResponseCandidate, ConductionResponseError> {
        self.solve_driven_with_materials_candidate(
            seeds,
            ResponseDrive::Traces { voltage, currents },
            assembly,
            field_controls,
            controls,
            cancelled,
        )
    }

    /// Close source response under an ideal equipotential total-current drive.
    /// Every frozen FEEC trial enforces the authored current on the original
    /// electrode geometry, then evaluates source data at its actual E and J.
    /// No unit-field source query or post-convergence amplitude scaling occurs.
    #[allow(clippy::too_many_arguments)]
    pub fn solve_current_candidate(
        &self,
        positive: &[usize],
        negative: &[usize],
        current_a: f64,
        assembly: WeightedAssemblyLimits,
        field_controls: ConductionControls,
        controls: ConductionResponseControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<ConductionResponseCandidate, ConductionResponseError> {
        self.solve_current_with_materials_candidate(
            self.conductors().to_vec(),
            positive,
            negative,
            current_a,
            assembly,
            field_controls,
            controls,
            cancelled,
        )
    }

    /// The same total-current closure at an explicitly refreshed original
    /// thermal/phase context. Source roles, pins and geometry remain original.
    #[allow(clippy::too_many_arguments)]
    pub fn solve_current_with_materials_candidate(
        &self,
        seeds: Vec<ResolvedConductorState>,
        positive: &[usize],
        negative: &[usize],
        current_a: f64,
        assembly: WeightedAssemblyLimits,
        field_controls: ConductionControls,
        controls: ConductionResponseControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<ConductionResponseCandidate, ConductionResponseError> {
        self.solve_driven_with_materials_candidate(
            seeds,
            ResponseDrive::Current {
                positive,
                negative,
                current_a,
            },
            assembly,
            field_controls,
            controls,
            cancelled,
        )
    }

    fn solve_driven_with_materials_candidate(
        &self,
        seeds: Vec<ResolvedConductorState>,
        drive: ResponseDrive<'_>,
        assembly: WeightedAssemblyLimits,
        field_controls: ConductionControls,
        controls: ConductionResponseControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<ConductionResponseCandidate, ConductionResponseError> {
        poll(cancelled)?;
        controls.validate()?;
        if seeds.len() != self.sources.len() || seeds.len() > assembly.max_cells {
            return Err(ConductionResponseError::Invalid(
                "complete bounded source response seeds required",
            ));
        }
        for (source, seed) in self.sources.iter().zip(&seeds) {
            poll(cancelled)?;
            source.admit(seed)?;
        }
        let mut applied = seeds.clone();
        let (mut maximum_sigma, mut maximum_current, mut absolute_power) =
            (0.0_f64, 0.0_f64, 0.0_f64);
        for iteration in 1..=controls.maximum_iterations {
            poll(cancelled)?;
            let native = PlanarConduction::assemble(
                self.mesh().clone(),
                applied,
                assembly,
                cancelled,
                |_| true,
            )?;
            let field =
                match drive {
                    ResponseDrive::Traces { voltage, currents } => native
                        .solve_with_current_traces(voltage, currents, field_controls, cancelled)?,
                    ResponseDrive::Current {
                        positive,
                        negative,
                        current_a,
                    } => native.solve_current(
                        positive,
                        negative,
                        current_a,
                        field_controls,
                        cancelled,
                    )?,
                };
            let mut evaluated = Vec::with_capacity(seeds.len());
            let mut next = Vec::with_capacity(seeds.len());
            let mut converged = true;
            maximum_sigma = 0.0;
            maximum_current = 0.0;
            absolute_power = 0.0;
            for (triangle, (source, seed)) in self.sources.iter().zip(&seeds).enumerate() {
                poll(cancelled)?;
                let e = field.electric_field_v_m()[triangle];
                let j = field.current_density_a_m2()[triangle];
                let conductor = source.at_fields(seed, e, j)?;
                let sigma = finite(conductor.resistivity_ohm_m().recip())?;
                let used = native.conductivity_s_m[triangle];
                let conductivity_defect_s_m = finite(sigma - used)?;
                let evaluated_j = conductor.conductivity().current_density_a_m2(e);
                let defect = [
                    finite(j[0] - evaluated_j[0])?,
                    finite(j[1] - evaluated_j[1])?,
                ];
                let current_defect = finite(defect[0].hypot(defect[1]))?;
                let volume =
                    self.mesh()
                        .face_measure(triangle)
                        .ok_or(ConductionResponseError::Invalid(
                            "native conductor cell volume missing",
                        ))?;
                let power_defect_w = finite(volume * e[0].mul_add(defect[0], e[1] * defect[1]))?;
                let tolerance = finite(
                    controls.absolute_conductivity_tolerance_s_m
                        + controls.relative_conductivity_tolerance * sigma.max(used),
                )?;
                maximum_sigma = maximum_sigma.max(conductivity_defect_s_m.abs());
                maximum_current = maximum_current.max(current_defect);
                absolute_power = finite(absolute_power + power_defect_w.abs())?;
                converged &= conductivity_defect_s_m.abs() <= tolerance
                    && current_defect <= controls.current_density_tolerance_a_m2;
                next.push(conductor.clone());
                evaluated.push(ConductionResponseSample {
                    triangle,
                    conductor,
                    conductivity_defect_s_m,
                    current_density_defect_a_m2: defect,
                    power_defect_w,
                });
            }
            poll(cancelled)?;
            if converged && absolute_power <= field_controls.power_tolerance_w {
                let report = ConductionResponseReport {
                    voltage_v: field.voltage_v().to_vec(),
                    electric_field_v_m: field.electric_field_v_m().to_vec(),
                    current_density_a_m2: field.current_density_a_m2().to_vec(),
                    cell_power_w: field.cell_power_w().to_vec(),
                    power_w: field.power_w(),
                    nodal_boundary_current_a: field.nodal_boundary_current_a().to_vec(),
                    nodal_applied_current_a: field.nodal_applied_current_a().to_vec(),
                    boundary_current_a: field.boundary_current_a(),
                    maximum_free_current_residual_a: field.maximum_free_current_residual_a(),
                    power_balance_residual_w: field.power_balance_residual_w(),
                    linear: field.report().cloned(),
                    evaluated,
                    field_solves: iteration,
                    source_evaluations: iteration.checked_mul(seeds.len()).ok_or(
                        ConductionResponseError::Invalid("source evaluation count overflows"),
                    )?,
                };
                return Ok(ConductionResponseCandidate {
                    problem: Self {
                        native,
                        sources: self.sources.clone(),
                    },
                    report,
                });
            }
            applied = next;
        }
        Err(ConductionResponseError::NotConverged {
            iterations: controls.maximum_iterations,
            maximum_conductivity_defect_s_m: maximum_sigma,
            maximum_current_density_defect_a_m2: maximum_current,
            absolute_power_defect_w: absolute_power,
        })
    }
}
