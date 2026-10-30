//! Physical balance at a supplied enthalpy endpoint, without advancing time.
//!
//! An outer boundary iteration must check the complete nodal residual under
//! its updated physical law. Total energy alone can hide opposite nodal defects,
//! and equal temperatures can hide different enthalpies on a latent plateau.

use fs_exec::Cx;
use fs_solver::norm2;

use super::{
    ASSEMBLY_TILE, DofMap, EnthalpyBackwardEuler, EnthalpyError, EnthalpyStepConfig, StepContext,
    energy_balance, finite, poll, radiation::boundary_at_temperature, PointwiseBoundary,
};
use crate::{
    AmbientRadiationPatch, ConductionProblem, PointwiseAmbientRadiation,
    PointwiseEnclosureRadiation, RobinFlux, ThermalInterfaces,
};

/// Physical balances of an arbitrary candidate, with immutable physical history.
///
/// These measurements do not assert convergence. The outer coupling owner must
/// apply its fixed residual and energy tolerances before publishing a new state.
#[derive(Debug, Clone)]
pub struct EnthalpyEndpointEvaluation {
    /// Candidate temperatures resolved from each vertex's declared phase chart, K.
    pub temperature: Vec<f64>,
    /// Original assembled balance in canonical vertex order, in joules.
    /// Prescribed-temperature rows are zero after their exact enthalpy check;
    /// their physical heat reaction remains in `dirichlet_in_w`. Coupled owners
    /// may solve these rows directly without assembling a second heat operator.
    pub residual_j: Vec<f64>,
    /// Euclidean norm of the full nodal balance, including local radiation
    /// where selected: `m*(h-h_old) + dt*(A(T,f)*T-b+q_rad)`, in joules.
    pub residual_norm_j: f64,
    /// Number of unconstrained enthalpy unknowns included in the norm.
    pub free_dofs: usize,
    /// Sum of reference nodal mass times actual specific-enthalpy change, J.
    pub stored_energy_change_j: f64,
    /// Integrated endpoint volumetric generation, W.
    pub source_w: f64,
    /// Outward endpoint Neumann heat, W.
    pub neumann_out_w: f64,
    /// Inward prescribed-temperature reaction, including fixed-node storage, W.
    pub dirichlet_in_w: f64,
    /// Original boundary convection only, W outward.
    pub convective_out_w: f64,
    /// Actual outward radiative heat, W, under the selected spatial model.
    /// For a closed enclosure this is its retained net internal exchange defect.
    pub radiation_out_w: f64,
    /// Storage minus dt times complete physical net input, J.
    pub energy_residual_j: f64,
    /// Original convective traces, excluding radiation from reservoir exchange.
    pub convective_robin_fluxes: Vec<RobinFlux>,
}

impl EnthalpyBackwardEuler<'_, '_> {
    /// Evaluate a candidate using the same transport and source/boundary
    /// assemblers as the nonlinear enthalpy step, including phase conductivity,
    /// heterogeneous charts, thermal contact and optional ambient radiation.
    ///
    /// No nonlinear/linear iteration is performed and no tangent is required:
    /// constitutive corners and latent plateaus remain valid physical states.
    /// Neither the candidate nor its old history is modified. The boundary
    /// supplies the outer owner's current convection coefficients; radiation
    /// coefficients are evaluated at this candidate's actual patch temperature.
    /// An off-equilibrium candidate returns its actual residual and energy
    /// defect even when they exceed the supplied step's acceptance tolerances.
    ///
    /// # Errors
    /// Refuses invalid step controls, mismatched mesh/fields, ambiguous prescribed
    /// temperatures or a candidate that violates prescribed enthalpy,
    /// invalid phase/material/boundary states, nonfinite arithmetic or cancellation.
    #[allow(clippy::too_many_arguments)]
    pub fn evaluate_endpoint(
        &self,
        cx: &Cx<'_>,
        problem: ConductionProblem<'_>,
        interfaces: Option<&ThermalInterfaces>,
        old_specific_h: &[f64],
        candidate_specific_h: &[f64],
        dt_s: f64,
        config: EnthalpyStepConfig,
        patches: Option<&[AmbientRadiationPatch]>,
    ) -> Result<EnthalpyEndpointEvaluation, EnthalpyError> {
        self.evaluate_endpoint_radiative(
            cx, problem, interfaces, old_specific_h, candidate_specific_h,
            dt_s, config, patches, None,
        )
    }

    /// Measure a supplied endpoint with local P1 ambient radiation, using the
    /// same physical residual and prescribed-temperature reaction as the primal.
    /// Original Robin heat alone is exposed to convection/coolant consumers.
    ///
    /// This performs no Newton step, needs no constitutive tangent and does
    /// not apply an acceptance gate: the outer coupling owner must check both
    /// the returned nodal residual and energy defect against its fixed limits.
    /// Radiation must retain this exact mesh and current physical boundary.
    #[allow(clippy::too_many_arguments)]
    pub fn evaluate_endpoint_with_pointwise_ambient_radiation(
        &self,
        cx: &Cx<'_>,
        problem: ConductionProblem<'_>,
        interfaces: Option<&ThermalInterfaces>,
        old_specific_h: &[f64],
        candidate_specific_h: &[f64],
        dt_s: f64,
        config: EnthalpyStepConfig,
        radiation: &PointwiseAmbientRadiation<'_>,
    ) -> Result<EnthalpyEndpointEvaluation, EnthalpyError> {
        poll(cx, 0)?;
        if !radiation.matches(problem.mesh, problem.boundary) {
            return Err(EnthalpyError::InvalidInput(
                "pointwise radiation must retain the exact transport mesh and boundary",
            ));
        }
        self.evaluate_endpoint_radiative(
            cx, problem, interfaces, old_specific_h, candidate_specific_h,
            dt_s, config, None, Some(PointwiseBoundary::Ambient(radiation)),
        )
    }

    /// Measure a supplied endpoint with local emission and fresh reflected
    /// enclosure exchange. The same radiosity and enclosure-closure gates as
    /// the primal apply; net internal radiation stays separate from convection.
    /// No history is advanced, tangent requested or convergence asserted.
    ///
    /// A closed enclosure can have zero total heat while individual vertices
    /// remain far from equilibrium. Consumers must check the nodal residual,
    /// not only the summed energy balance, before publishing an outer iterate.
    #[allow(clippy::too_many_arguments)]
    pub fn evaluate_endpoint_with_pointwise_enclosure_radiation(
        &self,
        cx: &Cx<'_>,
        problem: ConductionProblem<'_>,
        interfaces: Option<&ThermalInterfaces>,
        old_specific_h: &[f64],
        candidate_specific_h: &[f64],
        dt_s: f64,
        config: EnthalpyStepConfig,
        radiation: &PointwiseEnclosureRadiation<'_>,
    ) -> Result<EnthalpyEndpointEvaluation, EnthalpyError> {
        poll(cx, 0)?;
        if !radiation.matches(problem.mesh, problem.boundary) {
            return Err(EnthalpyError::InvalidInput(
                "pointwise enclosure must retain the exact transport mesh and boundary",
            ));
        }
        self.evaluate_endpoint_radiative(
            cx, problem, interfaces, old_specific_h, candidate_specific_h,
            dt_s, config, None, Some(PointwiseBoundary::Enclosure(radiation)),
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn evaluate_endpoint_radiative(
        &self,
        cx: &Cx<'_>,
        problem: ConductionProblem<'_>,
        interfaces: Option<&ThermalInterfaces>,
        old_specific_h: &[f64],
        candidate_specific_h: &[f64],
        dt_s: f64,
        config: EnthalpyStepConfig,
        patches: Option<&[AmbientRadiationPatch]>,
        pointwise_radiation: Option<PointwiseBoundary<'_>>,
    ) -> Result<EnthalpyEndpointEvaluation, EnthalpyError> {
        self.admit_step(cx, problem, old_specific_h, dt_s, config)?;
        config.newton.validate().map_err(EnthalpyError::Newton)?;
        let n = self.masses.len();
        if candidate_specific_h.len() != n {
            return Err(EnthalpyError::InvalidInput(
                "candidate enthalpy must match the reference mass mesh",
            ));
        }
        for (vertex, &value) in old_specific_h.iter().enumerate() {
            if vertex % ASSEMBLY_TILE == 0 {
                poll(cx, vertex)?;
            }
            self.curve_for_vertex(vertex)
                .state_at_specific_enthalpy(value)
                .map_err(|source| EnthalpyError::Phase { vertex, source })?;
        }
        let context = StepContext {
            storage: self,
            cx,
            problem,
            interfaces,
            old: old_specific_h,
            dt: dt_s,
            pointwise_radiation,
        };
        context.check_prescribed_values(candidate_specific_h)?;
        let temperature = context.temperatures(candidate_specific_h)?;
        let combined_boundary = match patches {
            Some(patches) => Some(boundary_at_temperature(cx, problem, patches, &temperature)?),
            None => None,
        };
        let physical = StepContext {
            problem: ConductionProblem {
                boundary: combined_boundary.as_ref().unwrap_or(problem.boundary),
                ..problem
            },
            ..context
        };
        let phase_scales = physical.phase_scales(candidate_specific_h)?;
        let system = physical.assemble(&temperature, phase_scales.as_deref())?;
        let mut residual = vec![0.0; n];
        physical.assembled_residual(candidate_specific_h, &temperature, &system, &mut residual)?;
        let residual_norm_j = finite(norm2(&residual))?;
        let dofs = DofMap::new_allow_fully_prescribed(problem.boundary, n)?;
        let (energy, convective_robin_fluxes) = energy_balance(
            self.mesh,
            problem.boundary,
            problem.source,
            &system,
            &dofs,
            &temperature,
        );
        poll(cx, n)?;
        let source_w = finite(energy.source_w)?;
        let neumann_out_w = finite(energy.neumann_out_w)?;
        let convective_out_w = finite(energy.robin_out_w)?;
        let dirichlet_in_w =
            physical.prescribed_heat_input_w(candidate_specific_h, &temperature, &system)?;
        for (index, flux) in convective_robin_fluxes.iter().enumerate() {
            if index % ASSEMBLY_TILE == 0 {
                poll(cx, index)?;
            }
            for value in [
                flux.area_m2,
                flux.mean_htc_w_per_m2_k,
                flux.mean_wall_temperature_k,
                flux.mean_reference_temperature_k,
                flux.heat_rate_w,
            ] {
                finite(value)?;
            }
        }
        // A radiation-only patch starts from a named zero-flux boundary, so it
        // has NO original convection row. Obtain its area/mean from the actual
        // radiative secant boundary using the same trace integrator as the
        // primal. Never return that radiation row as a coolant/air heat port.
        // Existing convection-only/mixed-Robin calls retain their exact path.
        let needs_radiative_traces = patches.is_some_and(|patches| patches.iter().any(|patch|
            !convective_robin_fluxes.iter().any(|row| row.region == patch.region())));
        let radiative_traces = if needs_radiative_traces {
            let boundary = combined_boundary.as_ref().ok_or(EnthalpyError::InvalidInput(
                "radiation endpoint lost its physical boundary",
            ))?;
            let (_, traces) = energy_balance(
                self.mesh, boundary, problem.source, &system, &dofs, &temperature,
            );
            poll(cx, n)?;
            traces
        } else {
            Vec::new()
        };
        let mut radiation_out_w = match pointwise_radiation {
            Some(radiation) => {
                let (_, report) = radiation.load_and_report(cx, &temperature)?;
                finite(report.outward_heat_w())?
            }
            None => 0.0,
        };
        if let Some(patches) = patches {
            for (index, patch) in patches.iter().enumerate() {
                poll(cx, index)?;
                let flux = convective_robin_fluxes
                    .iter()
                    .chain(&radiative_traces)
                    .find(|row| row.region == patch.region())
                    .ok_or(EnthalpyError::InvalidInput(
                        "evaluated radiation patch has no physical trace",
                    ))?;
                radiation_out_w = finite(
                    radiation_out_w
                        + finite(
                            patch.heat_flux_w_m2(flux.mean_wall_temperature_k)? * flux.area_m2,
                        )?,
                )?;
            }
        }
        let mut stored_energy_change_j = 0.0;
        for (vertex, (&new, &old)) in candidate_specific_h.iter().zip(old_specific_h).enumerate() {
            if vertex % ASSEMBLY_TILE == 0 {
                poll(cx, vertex)?;
            }
            stored_energy_change_j =
                finite(stored_energy_change_j + self.masses[vertex] * (new - old))?;
        }
        let net_input_w =
            finite(source_w + dirichlet_in_w - neumann_out_w - convective_out_w - radiation_out_w)?;
        let energy_residual_j = finite(stored_energy_change_j - dt_s * net_input_w)?;
        poll(cx, n)?;
        Ok(EnthalpyEndpointEvaluation {
            temperature,
            residual_j: residual,
            residual_norm_j,
            free_dofs: dofs.n(),
            stored_energy_change_j,
            source_w,
            neumann_out_w,
            dirichlet_in_w,
            convective_out_w,
            radiation_out_w,
            energy_residual_j,
            convective_robin_fluxes,
        })
    }
}
