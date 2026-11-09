//! Pointwise gray radiation in the existing capacity/bulk Newton response.
//! Surface loads and tangents belong to the shared P1 radiation owner. Only
//! this adapter adds the storage coordinate, boundary elimination and time scale.

use std::cell::RefCell;

use fs_exec::Cx;
use fs_solver::{CsrOp, LinearOp};

use super::{BackwardEuler, EndpointEvaluation, NonlinearStepConfig, NonlinearStepSolution,
    StepConfig, finite, invalid, poll};
use super::nonlinear::Evaluation;
use crate::{ConductionError, ConductionProblem, PointwiseAmbientRadiation,
    PointwiseRadiationAssembly, PointwiseRadiationReport, ThermalInterfaces, assemble::DofMap};

/// Actual local emission and complete physical balance at one unchanged field.
#[derive(Debug, Clone)]
pub struct PointwiseRadiationEndpointEvaluation {
    /// Original convection/Neumann heat, local radiation, full free joule
    /// residual and prescribed reactions including boundary-node storage.
    pub physical: EndpointEvaluation,
    /// Integrated P1 surface emission with original emissivity/bath identities.
    pub radiation: PointwiseRadiationReport,
}

/// Accepted capacity endpoint with implicit spatially resolved gray emission.
#[derive(Debug, Clone)]
pub struct PointwiseRadiationStepSolution {
    /// The existing Newton response. Robin heat remains original convection;
    /// its energy residual additionally includes the radiation reported below.
    pub conduction: NonlinearStepSolution,
    /// Independently re-evaluated physical endpoint; no extra solve is run.
    pub physical: EndpointEvaluation,
    /// Local surface heat, separate from coolant heat and prescribed flux.
    pub radiation: PointwiseRadiationReport,
}

impl BackwardEuler<'_> {
    /// Advance capacity and k(T) with exact P1 ambient radiation in one Newton solve.
    ///
    /// Every residual and rejected line-search trial retains the SAME old
    /// temperature field. The radiation owner evaluates local T^4 and supplies
    /// its analytic tangent; no patch-mean secant or apparent enthalpy replaces
    /// either model. Radiation on fixed vertices enters the bath reaction once.
    /// Frozen contact operators and original Robin/Neumann heat are retained.
    ///
    /// Work limits and immutable-history semantics are those of
    /// [`Self::advance_nonlinear_prescribed`]. This method supplies no live-contact,
    /// material-emissivity, enclosure or trajectory-adjoint composition.
    ///
    /// # Errors
    /// Refuses invalid/mismatched input, unsupported source temperatures,
    /// exhausted nonlinear/Krylov work, cancellation or failed physical gates.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_prescribed_with_pointwise_ambient_radiation(
        &self, cx: &Cx<'_>, problem: ConductionProblem<'_>,
        interfaces: Option<&ThermalInterfaces>, old: &[f64], dt_s: f64,
        config: StepConfig, nonlinear: NonlinearStepConfig,
        radiation: &PointwiseAmbientRadiation<'_>,
    ) -> Result<PointwiseRadiationStepSolution, ConductionError> {
        let conduction = self.advance_nonlinear_with_boundary(cx, problem, interfaces,
            old, dt_s, config, nonlinear, true, Some(radiation))?;
        let audit = self.evaluate_prescribed_pointwise_ambient_endpoint(cx, problem,
            interfaces, old, &conduction.step.temperature, dt_s, config, radiation)?;
        if audit.physical.residual_norm_j > conduction.threshold_j
            || audit.physical.energy_residual_j.abs() > config.energy_tolerance_j {
            return Err(invalid("pointwise capacity endpoint failed its original physical residual or energy gate"));
        }
        Ok(PointwiseRadiationStepSolution {
            conduction, physical: audit.physical, radiation: audit.radiation,
        })
    }

    /// Audit an arbitrary candidate using local surface emission and actual storage.
    /// No nonlinear solve or history update occurs. Off-equilibrium candidates
    /// retain their defects, so an outer coupling owner can apply its own original
    /// residual target. `combined_boundary` is the unchanged convective boundary:
    /// pointwise radiation is an independent nodal load, never a coolant Robin row.
    ///
    /// # Errors
    /// Refuses invalid endpoint data, source-domain/identity mismatch or cancellation.
    #[allow(clippy::too_many_arguments)]
    pub fn evaluate_prescribed_pointwise_ambient_endpoint(
        &self, cx: &Cx<'_>, problem: ConductionProblem<'_>,
        interfaces: Option<&ThermalInterfaces>, old: &[f64], candidate: &[f64],
        dt_s: f64, config: StepConfig, radiation: &PointwiseAmbientRadiation<'_>,
    ) -> Result<PointwiseRadiationEndpointEvaluation, ConductionError> {
        let (physical, report) = self.evaluate_prescribed_endpoint_with_pointwise(
            cx, problem, interfaces, old, candidate, dt_s, config, None, Some(radiation))?;
        Ok(PointwiseRadiationEndpointEvaluation {
            physical,
            radiation: report.ok_or_else(|| invalid("pointwise endpoint lost its original radiation report"))?,
        })
    }

    #[allow(clippy::too_many_arguments)]
    pub(super) fn evaluate_nonlinear_endpoint(
        &self, cx: &Cx<'_>, problem: ConductionProblem<'_>,
        interfaces: Option<&ThermalInterfaces>, old: &[f64], temperature: &[f64],
        dt: f64, dofs: &DofMap, radiation: Option<&PointwiseAmbientRadiation<'_>>,
    ) -> Result<Evaluation, ConductionError> {
        let mut value = self.evaluate_endpoint(cx, problem, interfaces, old, temperature, dt, dofs)?;
        if let Some(radiation) = radiation {
            let (load, _) = radiation.load_and_report(cx, temperature)?;
            for (index, (&vertex, residual)) in dofs.free().iter().zip(&mut value.residual).enumerate() {
                if index % 512 == 0 { poll(cx, index)?; }
                *residual = finite(*residual + dt * load[vertex])?;
            }
            value.norm = finite(fs_solver::norm2(&value.residual))?;
        }
        Ok(value)
    }
}

/// Free-coordinate capacity/bulk matrix plus the original matrix-free surface
/// tangent. Fixed increments are zero; fixed reaction rows never enter Krylov.
pub(super) struct RadiativeDirection<'a, 'cx> {
    cx: &'a Cx<'cx>,
    base: CsrOp,
    radiation: Option<(&'a PointwiseRadiationAssembly, &'a DofMap, f64)>,
    scratch: RefCell<(Vec<f64>, Vec<f64>)>,
    failure: RefCell<Option<ConductionError>>,
}

impl<'a, 'cx> RadiativeDirection<'a, 'cx> {
    pub(super) fn new(cx: &'a Cx<'cx>, base: CsrOp,
        radiation: Option<(&'a PointwiseRadiationAssembly, &'a DofMap, f64)>) -> Self {
        let n = radiation.map_or(0, |(assembly, _, _)| assembly.outward_load_w.len());
        Self { cx, base, radiation, scratch: RefCell::new((vec![0.0; n], vec![0.0; n])),
            failure: RefCell::new(None) }
    }

    pub(super) fn check(&self) -> Result<(), ConductionError> {
        match self.failure.take() { Some(error) => Err(error), None => Ok(()) }
    }

    fn add_radiation(&self, x: &[f64], y: &mut [f64]) -> Result<(), ConductionError> {
        let Some((radiation, dofs, dt)) = self.radiation else { return Ok(()); };
        poll(self.cx, 0)?;
        let mut scratch = self.scratch.borrow_mut();
        let (direction, output) = &mut *scratch;
        direction.fill(0.0);
        for (index, &vertex) in dofs.free().iter().enumerate() {
            if index % 512 == 0 { poll(self.cx, index)?; }
            direction[vertex] = x[index];
        }
        radiation.tangent_apply(self.cx, direction, output)?;
        for (index, &vertex) in dofs.free().iter().enumerate() {
            if index % 512 == 0 { poll(self.cx, index)?; }
            y[index] = finite(y[index] + dt * output[vertex])?;
        }
        Ok(())
    }
}

impl LinearOp for RadiativeDirection<'_, '_> {
    fn n(&self) -> usize { self.base.n() }

    fn apply(&self, x: &[f64], y: &mut [f64]) {
        self.base.apply(x, y);
        if let Err(error) = self.add_radiation(x, y) {
            self.failure.replace(Some(error));
            y.fill(f64::NAN);
        }
    }

    fn apply_transpose(&self, x: &[f64], y: &mut [f64]) {
        self.base.apply_transpose(x, y);
        // Frozen-emissivity local emission is symmetric. The bulk k'(T)
        // matrix above still requires its genuine nonsymmetric transpose.
        if let Err(error) = self.add_radiation(x, y) {
            self.failure.replace(Some(error));
            y.fill(f64::NAN);
        }
    }
}
