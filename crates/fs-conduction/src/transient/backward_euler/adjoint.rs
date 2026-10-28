//! Discrete endpoint derivatives for immutable-history backward Euler.
//!
//! R = C(T-new - T-old)/dt + A(T-new)T-new - b. Its endpoint Jacobian is
//! C/dt + J-steady, NOT a steady-state tangent and NOT J evaluated at T-old.
//! A coupled air adjoint can use the same RobinResponse, then carry its total
//! nodal-load multiplier backwards as (C/dt)^T lambda. Time grids, geometry,
//! prescribed temperatures, material laws and contact resistance stay fixed.

use super::*;
use crate::adjoint::robin::RobinResponse;
use crate::assemble::assemble_jacobian_with_optional_interfaces;

mod contact;

/// A residual- and energy-checked transient endpoint and its response operator.
/// Dereferencing exposes the same Robin tangent/adjoint operations as the
/// steady producer, without representing this state as a steady solution.
pub struct StepLinearization<'a> {
    primal: StepSolution,
    response: RobinResponse,
    capacity: &'a Csr,
    mesh: &'a ConductionMesh,
    dofs: DofMap,
    old: Vec<f64>,
    dt_s: f64,
    source_load_w: Vec<f64>,
}

impl std::ops::Deref for StepLinearization<'_> {
    type Target = RobinResponse;
    fn deref(&self) -> &Self::Target { &self.response }
}

impl BackwardEuler<'_> {
    /// Solve one real endpoint and prepare its discrete response operator.
    /// Use the same old field, dt, references and policies as the forward run
    /// when reconstructing an accepted endpoint during a reverse traversal.
    ///
    /// No physical history is mutated. The primal obeys the existing residual,
    /// energy, material-range and cancellation gates. The derivative has the
    /// caller's linear budget. A material kink retains the primal but refuses
    /// tangent/adjoint requests, as on the steady Robin path.
    ///
    /// # Errors
    /// All errors from advance/advance_nonlinear and Robin port admission.
    #[allow(clippy::too_many_arguments)]
    pub fn linearize_step<'a>(
        &'a self, cx: &Cx<'_>, problem: ConductionProblem<'_>,
        interfaces: Option<&ThermalInterfaces>, old: &[f64], dt_s: f64,
        config: StepConfig, nonlinear: Option<NonlinearStepConfig>, regions: &[&str],
    ) -> Result<StepLinearization<'a>, ConductionError> {
        let primal = match nonlinear {
            Some(policy) => self.advance_nonlinear(cx, problem, interfaces, old, dt_s, config, policy)?.step,
            None => self.advance(cx, problem, interfaces, old, dt_s, config)?,
        };
        self.linearize_primal(cx, problem, interfaces, old, dt_s, config, regions, primal)
    }

    /// Solve and differentiate an endpoint with prescribed temperatures that
    /// may differ from the old physical field. This is the derivative companion
    /// of `advance_prescribed` / `advance_nonlinear_prescribed`, not a steady
    /// solve and not an overwrite of the supplied history.
    ///
    /// The same endpoint Dirichlet lift, boundary-node storage reaction,
    /// nonlinear material Jacobian, matching contact and energy gates apply.
    /// History, source and capacity pullbacks retain their existing meanings;
    /// prescribed ENDPOINT temperatures are held fixed, not differentiated.
    /// Callers reconstructing a trajectory must compare the returned primal
    /// with their retained endpoint before using any derivatives.
    ///
    /// # Errors
    /// Refuses every invalid endpoint, material, port, work or cancellation
    /// state rejected by the corresponding prescribed forward producer. No
    /// endpoint or derivative is published after a failed forward solve.
    #[allow(clippy::too_many_arguments)]
    pub fn linearize_prescribed_step<'a>(
        &'a self, cx: &Cx<'_>, problem: ConductionProblem<'_>,
        interfaces: Option<&ThermalInterfaces>, old: &[f64], dt_s: f64,
        config: StepConfig, nonlinear: Option<NonlinearStepConfig>, regions: &[&str],
    ) -> Result<StepLinearization<'a>, ConductionError> {
        let primal = match nonlinear {
            Some(policy) => self.advance_nonlinear_prescribed(
                cx, problem, interfaces, old, dt_s, config, policy)?.step,
            None => self.advance_prescribed(cx, problem, interfaces, old, dt_s, config)?,
        };
        self.linearize_primal(cx, problem, interfaces, old, dt_s, config, regions, primal)
    }

    /// Bind the actual retained prescribed endpoint without another primal solve.
    ///
    /// The caller supplies the complete Robin law values evaluated at
    /// `candidate`, including any combined convection/radiation coefficients,
    /// and the fixed physical residual target accepted by its forward owner.
    /// This method independently reassembles C(T - T_old) + dt(K(T)T - b),
    /// requires its free-node norm to be at most `residual_tolerance_j`, and
    /// checks the unchanged energy budget. Prescribed endpoint temperatures,
    /// material ranges and matching contacts use the same production audit as
    /// the forward endpoint owner. No candidate or history value is changed.
    ///
    /// The returned response contains C/dt + J(T) at this exact field. For
    /// state-dependent boundary laws, call
    /// [`StepLinearization::pullback_with_robin_feedback`] with their actual
    /// coefficient and reference derivatives; the ordinary Robin pullback
    /// holds those boundary values fixed. Radiative reservoirs and convection
    /// controls must remain separately composed by their physical owner.
    ///
    /// No Krylov primal solve runs: the primal iteration count and its inner
    /// relative-residual field are zero. Those fields do not replace the
    /// independently checked absolute physical residual and energy gates.
    /// A finite zero physical target is retained exactly and requires an
    /// exactly zero free residual, including a fully prescribed endpoint.
    /// The returned owner borrows only this capacity engine, not the supplied
    /// candidate, boundary or law arrays.
    ///
    /// # Errors
    ///
    /// Refuses invalid fields, budgets, material/contact data, a retained
    /// endpoint exceeding either physical gate, invalid selected Robin ports,
    /// or cancellation. Material kinks refuse derivative requests, as for a
    /// newly solved endpoint. No hidden replacement solve or widened target
    /// is attempted.
    #[allow(clippy::too_many_arguments)]
    pub fn linearize_prescribed_endpoint<'a>(
        &'a self, cx: &Cx<'_>, problem: ConductionProblem<'_>,
        interfaces: Option<&ThermalInterfaces>, old: &[f64], candidate: &[f64],
        dt_s: f64, config: StepConfig, residual_tolerance_j: f64, regions: &[&str],
    ) -> Result<StepLinearization<'a>, ConductionError> {
        poll(cx, 0)?;
        if !residual_tolerance_j.is_finite() || residual_tolerance_j < 0.0 {
            return Err(invalid("retained endpoint needs a nonnegative finite physical residual target in joules"));
        }
        let evaluation = self.evaluate_prescribed_endpoint(
            cx, problem, interfaces, old, candidate, dt_s, config, None,
        )?;
        if evaluation.residual_norm_j > residual_tolerance_j {
            return Err(invalid(&format!(
                "retained endpoint residual {} J exceeds {residual_tolerance_j} J",
                evaluation.residual_norm_j,
            )));
        }
        if evaluation.energy_residual_j.abs() > config.energy_tolerance_j {
            return Err(invalid(&format!(
                "retained endpoint energy residual {} J exceeds {} J",
                evaluation.energy_residual_j, config.energy_tolerance_j,
            )));
        }
        let contact_fluxes = match interfaces {
            Some(interfaces) => interfaces.fluxes(candidate)?,
            None => Vec::new(),
        };
        let primal = StepSolution {
            temperature: candidate.to_vec(),
            robin_fluxes: evaluation.convective_robin_fluxes,
            contact_fluxes,
            source_w: evaluation.source_w,
            neumann_out_w: evaluation.neumann_out_w,
            robin_out_w: evaluation.convective_out_w,
            dirichlet_in_w: evaluation.dirichlet_in_w,
            stored_energy_change_j: evaluation.stored_energy_change_j,
            energy_residual_j: evaluation.energy_residual_j,
            relative_residual: 0.0,
            krylov_iterations: 0,
        };
        self.linearize_primal(cx, problem, interfaces, old, dt_s, config, regions, primal)
    }

    // All entry points share one endpoint-Jacobian and source projection.
    // Only an in-crate residual- and energy-checked state reaches this helper;
    // a caller-constructed StepSolution is never accepted.
    #[allow(clippy::too_many_arguments)]
    fn linearize_primal<'a>(
        &'a self, cx: &Cx<'_>, problem: ConductionProblem<'_>,
        interfaces: Option<&ThermalInterfaces>, old: &[f64], dt_s: f64,
        config: StepConfig, regions: &[&str], primal: StepSolution,
    ) -> Result<StepLinearization<'a>, ConductionError> {
        let dofs = DofMap::new(problem.boundary, self.mesh.vertex_count())?;
        let jacobian = assemble_jacobian_with_optional_interfaces(cx, self.mesh,
            problem.boundary, problem.material, &primal.temperature, interfaces, problem.element_materials)?;
        let inverse_dt = finite(1.0 / dt_s)?;
        let full = axpy_csr(&self.capacity, inverse_dt, &jacobian, 1.0);
        let (matrix, _) = reduce_matrix_and_lift(&full, &dofs);
        let response = RobinResponse::for_endpoint(cx, problem, &primal, matrix,
            dofs.clone(), config.linear, regions)?;
        // Same P1 source integration and accumulation order as assemble.rs;
        // isolate the body load rather than subtracting large Robin loads.
        let mut source_load_w = vec![0.0; self.mesh.vertex_count()];
        for element in 0..self.mesh.element_count() {
            if element % 512 == 0 { poll(cx, element)?; }
            let tet = self.mesh.complex().tets[element];
            let volume = self.mesh.element_volume(element);
            for a in 0..4 {
                let mut load = 0.0;
                for (b, &vertex) in tet.iter().enumerate() {
                    let weight = if a == b { volume / 10.0 } else { volume / 20.0 };
                    load = weight.mul_add(problem.source.at(vertex as usize), load);
                }
                let vertex = tet[a] as usize;
                source_load_w[vertex] = finite(source_load_w[vertex] + load)?;
            }
        }
        poll(cx, 0)?;
        Ok(StepLinearization { primal, response, capacity: &self.capacity, mesh: self.mesh, dofs,
            old: old.to_vec(), dt_s, source_load_w })
    }
}

impl StepLinearization<'_> {
    /// The actual residual- and energy-checked transient primal.
    #[must_use]
    pub const fn primal(&self) -> &StepSolution { &self.primal }

    /// Differentiate a nodal-temperature goal through complete mean-Robin feedback.
    ///
    /// All arrays follow the selected port order. `htc_wall_derivatives`
    /// supplies dh_i/dmean_i in W/(m² K²); `reference_wall_derivatives`
    /// supplies the local dref_i/dmean_i in K/K. The flat row-major m-by-m
    /// `reference_feedback` supplies ADDITIONAL dref_i/dmean_j entries.
    /// Diagonal entries add to the local reference slope: never supply the
    /// same derivative twice. Unselected boundary values remain fixed.
    ///
    /// The response uses the same consistent P1 face mass and area means as
    /// the Robin owner. Its full nonsymmetric C/dt + J(T) + feedback operator
    /// is transposed by the existing bounded FGMRES path. The returned physical
    /// multiplier is checked again after rescaling; no new primal solve runs.
    /// At most 64 ports and four nodal factor vectors per port must fit
    /// `max_feedback_entries`. This cap is not a total-memory estimate.
    ///
    /// Reference and log(h) entries are independent partials of the BOUND
    /// Robin rows. A combined convection/radiation law must map them back to
    /// its original physical controls. Feed `nodal_load` directly into the
    /// history, source and capacity pullbacks below: C/dt is already in the
    /// Jacobian, so multiplying those contractions by dt again is incorrect.
    ///
    /// # Errors
    ///
    /// Refuses malformed/nonfinite controls, a material kink, exhausted
    /// storage/work allowances, a failed complete transpose residual, or
    /// cancellation. A zero goal does not bypass law or budget admission.
    pub fn pullback_with_robin_feedback(
        &self, cx: &Cx<'_>, htc_wall_derivatives: &[f64],
        reference_wall_derivatives: &[f64], reference_feedback: &[f64],
        nodal_weights: &[f64], max_feedback_entries: usize,
    ) -> Result<crate::adjoint::robin::RobinGradient, ConductionError> {
        self.response.pullback_endpoint_with_mean_feedback(
            cx, htc_wall_derivatives, reference_wall_derivatives,
            reference_feedback, nodal_weights, max_feedback_entries,
        )
    }

    fn multiplier(&self, cx: &Cx<'_>, values: &[f64]) -> Result<Vec<f64>, ConductionError> {
        if values.len() != self.old.len() {
            return Err(ConductionError::FieldLength { field: "transient adjoint multiplier",
                expected: self.old.len(), found: values.len() });
        }
        let mut result = values.to_vec();
        for (i, &value) in values.iter().enumerate() {
            if i % 512 == 0 { poll(cx, i)?; }
            finite(value)?;
        }
        for &vertex in self.dofs.fixed() { result[vertex] = 0.0; }
        Ok(result)
    }

    /// Pull a TOTAL endpoint load multiplier back to the previous temperature.
    /// In a coupled problem pass the converged coupled nodal-load adjoint, not
    /// a frozen-reference solid adjoint. Prescribed-temperature controls are zero.
    pub fn previous_temperature_pullback(&self, cx: &Cx<'_>, nodal_load_adjoint: &[f64])
        -> Result<Vec<f64>, ConductionError> {
        let lambda = self.multiplier(cx, nodal_load_adjoint)?;
        let mut result = vec![0.0; lambda.len()];
        self.capacity.spmv(&lambda, &mut result);
        for (i, value) in result.iter_mut().enumerate() {
            if i % 512 == 0 { poll(cx, i)?; }
            *value = finite(*value / self.dt_s)?;
        }
        for &vertex in self.dofs.fixed() { result[vertex] = 0.0; }
        poll(cx, 0)?;
        Ok(result)
    }

    /// Pull a TOTAL endpoint load adjoint back to P1 nodal source densities.
    ///
    /// The returned vector is `M_source^T lambda`, with the same consistent
    /// tetrahedral integration `M_ab = V(1+delta_ab)/20` as the primal source.
    /// For a temperature objective its units are K per (W/m3), not K/W.
    /// Contract with a density derivative to obtain a component-watt or other
    /// source control. This requires no solve and no nonzero baseline source.
    /// The endpoint Jacobian already contains C/dt; multiplying by dt again
    /// would give an incorrect trajectory derivative.
    ///
    /// Prescribed equation multipliers are zeroed, but source-density entries
    /// on prescribed vertices are NOT: their P1 support can load neighboring
    /// free equations. This is distinct from changing prescribed temperatures.
    ///
    /// # Errors
    /// Refuses a malformed/nonfinite multiplier, arithmetic overflow, or
    /// cancellation. No partially accumulated vector is returned.
    pub fn source_density_pullback(&self, cx: &Cx<'_>, nodal_load_adjoint: &[f64])
        -> Result<Vec<f64>, ConductionError> {
        let lambda = self.multiplier(cx, nodal_load_adjoint)?;
        let mut result = vec![0.0; lambda.len()];
        for element in 0..self.mesh.element_count() {
            if element % 512 == 0 { poll(cx, element)?; }
            let tet = self.mesh.complex().tets[element];
            let volume = self.mesh.element_volume(element);
            for (b, &vertex) in tet.iter().enumerate() {
                let mut value = 0.0;
                for (a, &row) in tet.iter().enumerate() {
                    let mass = if a == b { volume / 10.0 } else { volume / 20.0 };
                    value = finite(mass.mul_add(lambda[row as usize], value))?;
                }
                let slot = &mut result[vertex as usize];
                *slot = finite(*slot + value)?;
            }
        }
        poll(cx, 0)?;
        Ok(result)
    }

    /// Derivative with respect to a multiplier on THIS endpoint's entire
    /// volumetric source at multiplier one. Zero source gives zero, not a
    /// division by its amplitude. This is not a Neumann or Robin control.
    pub fn source_multiplier_pullback(&self, cx: &Cx<'_>, nodal_load_adjoint: &[f64])
        -> Result<f64, ConductionError> {
        let lambda = self.multiplier(cx, nodal_load_adjoint)?;
        let mut value = 0.0;
        for (i, (&weight, &load)) in lambda.iter().zip(&self.source_load_w).enumerate() {
            if i % 512 == 0 { poll(cx, i)?; }
            value = finite(value + finite(weight * load)?)?;
        }
        Ok(value)
    }

    /// Derivative with respect to a common multiplier on the entire capacity
    /// matrix at multiplier one, holding the OLD state independently fixed.
    /// A trajectory sums this local contraction while propagating history.
    pub fn capacity_multiplier_pullback(&self, cx: &Cx<'_>, nodal_load_adjoint: &[f64])
        -> Result<f64, ConductionError> {
        let lambda = self.multiplier(cx, nodal_load_adjoint)?;
        let delta = self.primal.temperature.iter().zip(&self.old)
            .map(|(new, old)| finite(new - old)).collect::<Result<Vec<_>, _>>()?;
        let mut storage = vec![0.0; delta.len()];
        self.capacity.spmv(&delta, &mut storage);
        let mut value = 0.0;
        for (i, (&weight, &change)) in lambda.iter().zip(&storage).enumerate() {
            if i % 512 == 0 { poll(cx, i)?; }
            value = finite(value - finite(weight * finite(change / self.dt_s)?)?)?;
        }
        Ok(value)
    }
}