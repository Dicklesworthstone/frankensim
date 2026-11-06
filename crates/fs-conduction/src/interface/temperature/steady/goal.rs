//! Goal comparison with the original temperature-dependent interface law.

use super::{
    ContactSteadyLinearization, TemperatureContactEvaluation, TemperatureContacts, finite, invalid,
    poll, positive,
};
use crate::adjoint::DiscreteGoalComparison;
use crate::adjoint::robin::{RobinGradient, admit_linear, true_residual};
use crate::assemble::{assemble_operator_scaled_with_interfaces, reduce, residual};
use crate::{
    ConductionError, ConductionProblem, DofMap, LinearConfig, ScalarField, ThermalBc,
    ThermalBoundary,
};
use fs_exec::Cx;

impl TemperatureContacts<'_, '_> {
    /// Compare two fields using the complete source-backed steady residual.
    ///
    /// The reference is rechecked against the original absolute free-node and
    /// heat-balance targets in watts, independently of the dual's relative
    /// tolerance. Its complete bulk/contact tangent supplies one true transpose
    /// solve. Both residuals query the original selected contact law at their
    /// own temperatures, retaining independent A/B traces, matching face-mean
    /// or pointwise rules, and nonmatching common-refinement quadrature.
    ///
    /// Returned nodal contributions are lambda_i [R(T_ref)-R(T_approx)]_i.
    /// The observed remainder closes their sum to weights . (T_ref-T_approx);
    /// it includes the finite nonlinear remainder and cannot be discarded.
    /// These discrete diagnostics can guide refinement but are not continuum
    /// bounds. No extra primal solve or source-claim selection is performed.
    ///
    /// Geometry, source loads, nonthermal contact coordinates and prescribed
    /// boundary data are identical at both fields. Derived boundary feedback
    /// needs its own complete residual owner and is not supplied by this API.
    /// Approximate fields may cross source segments; a nonsmooth reference
    /// refuses a unique tangent. Phase-state-dependent contacts retain their
    /// temperature-only refusal. Inputs and work limits are checked before
    /// allocating the numerical comparison state; cancellation is propagated.
    #[allow(clippy::too_many_arguments)]
    pub fn compare_steady_goal(
        &self,
        cx: &Cx<'_>,
        problem: ConductionProblem<'_>,
        reference_temperature: &[f64],
        approximate_temperature: &[f64],
        linear: LinearConfig,
        residual_tolerance_w: f64,
        energy_tolerance_w: f64,
        full_nodal_weights: &[f64],
        max_face_pairs: usize,
    ) -> Result<DiscreteGoalComparison, ConductionError> {
        self.admit_goal_comparison(
            cx,
            problem,
            reference_temperature,
            approximate_temperature,
            linear,
            residual_tolerance_w,
            energy_tolerance_w,
            full_nodal_weights,
            max_face_pairs,
        )?;
        let endpoint = self.linearize_steady_endpoint(
            cx,
            problem,
            reference_temperature,
            linear,
            residual_tolerance_w,
            energy_tolerance_w,
            &[],
            max_face_pairs,
        )?;
        let approximate_contacts = self.evaluate(cx, problem.boundary, approximate_temperature)?;
        compare_goal_fields(
            cx,
            problem,
            problem.boundary,
            approximate_temperature,
            full_nodal_weights,
            &endpoint,
            &approximate_contacts,
            || endpoint.pullback(cx, full_nodal_weights, &[], &[]),
        )
    }

    /// Compare live-contact fields under actual area-mean-dependent Robin laws.
    ///
    /// In selected region order, `reference_points` contains `[h, r, dh/dmean,
    /// dr/dmean]` at the reference field, in W/(m² K), K, W/(m² K²), and K/K.
    /// `approximate_values` contains `[h, r]` at the approximate field. The
    /// row-major `reference_feedback[i*m+j]` is ADDITIONAL dr_i/dmean_j;
    /// its diagonal must not repeat the local reference slope. The boundary
    /// law owner supplies these actual values, domains and derivatives.
    ///
    /// Only unique existing uniform Robin regions are rebound. Both residuals
    /// use their own boundary values AND original contact-source queries. One
    /// complete transpose includes bulk/contact, local h/reference slopes and
    /// cross-region reference feedback, retaining consistent P1 face mass.
    /// Original geometry, fixed temperatures, other boundaries and loads are
    /// preserved. The reference passes the unchanged absolute watt gates even
    /// for a zero goal. No extra primal solve or frozen approximate law enters.
    ///
    /// The signed contributions and observed nonlinear remainder have the same
    /// discrete meaning as `compare_steady_goal`; no continuum bound follows.
    /// At most 64 selected regions and four nodal factor vectors per region
    /// must fit `max_feedback_entries`, separately from `max_face_pairs`.
    /// Malformed, duplicate, nonfinite, over-budget or cancelled input refuses
    /// before any boundary clone or numerical state allocation. This method
    /// supplies no pointwise or radiative boundary-feedback law.
    #[allow(clippy::too_many_arguments)]
    pub fn compare_steady_goal_with_robin_feedback(
        &self,
        cx: &Cx<'_>,
        problem: ConductionProblem<'_>,
        reference_temperature: &[f64],
        approximate_temperature: &[f64],
        linear: LinearConfig,
        residual_tolerance_w: f64,
        energy_tolerance_w: f64,
        regions: &[&str],
        reference_points: &[[f64; 4]],
        approximate_values: &[[f64; 2]],
        reference_feedback: &[f64],
        full_nodal_weights: &[f64],
        max_feedback_entries: usize,
        max_face_pairs: usize,
    ) -> Result<DiscreteGoalComparison, ConductionError> {
        self.admit_goal_comparison(
            cx,
            problem,
            reference_temperature,
            approximate_temperature,
            linear,
            residual_tolerance_w,
            energy_tolerance_w,
            full_nodal_weights,
            max_face_pairs,
        )?;
        let m = regions.len();
        let entries = problem
            .mesh
            .vertex_count()
            .checked_mul(4)
            .and_then(|n| n.checked_mul(m));
        if m > 64
            || reference_points.len() != m
            || approximate_values.len() != m
            || m.checked_mul(m) != Some(reference_feedback.len())
            || entries.is_none_or(|n| {
                n > max_feedback_entries
                    || n.checked_mul(8)
                        .is_none_or(|bytes| bytes > isize::MAX as usize)
            })
        {
            return Err(invalid(
                "<steady-contact-goal>",
                "Robin feedback needs matching law/map dimensions within 64 regions and its factor-entry allowance",
            ));
        }
        for (i, &value) in reference_feedback.iter().enumerate() {
            if i % 256 == 0 {
                poll(cx, i)?;
            }
            finite(value)?;
        }
        let mut reference_rows = [(0, 0.0, 0.0); 64];
        let mut approximate_rows = [(0, 0.0, 0.0); 64];
        let mut slopes = [0.0; 64];
        let mut reference_slopes = [0.0; 64];
        for (j, &name) in regions.iter().enumerate() {
            poll(cx, j)?;
            for &value in reference_points[j].iter().chain(&approximate_values[j]) {
                finite(value)?;
            }
            let index = problem
                .boundary
                .region_names()
                .iter()
                .position(|region| region == name)
                .ok_or_else(|| invalid("<steady-contact-goal>", "unknown Robin feedback region"))?;
            let [h, r, dh, dr] = reference_points[j];
            let [ha, ra] = approximate_values[j];
            if regions[..j].contains(&name)
                || h <= 0.0
                || ha <= 0.0
                || !matches!(
                    problem.boundary.conditions().get(index),
                    Some(ThermalBc::Robin {
                        htc: ScalarField::Uniform(_),
                        t_ref: ScalarField::Uniform(_),
                    })
                )
            {
                return Err(invalid(
                    "<steady-contact-goal>",
                    "unique existing uniform Robin regions and positive heat-transfer coefficients are required",
                ));
            }
            reference_rows[j] = (index, h, r);
            approximate_rows[j] = (index, ha, ra);
            slopes[j] = dh;
            reference_slopes[j] = dr;
        }
        let reference_boundary = problem
            .boundary
            .with_uniform_robin_replacements(&reference_rows[..m])?;
        let approximate_boundary = problem
            .boundary
            .with_uniform_robin_replacements(&approximate_rows[..m])?;
        let reference_problem = ConductionProblem {
            boundary: &reference_boundary,
            ..problem
        };
        let endpoint = self.linearize_steady_endpoint(
            cx,
            reference_problem,
            reference_temperature,
            linear,
            residual_tolerance_w,
            energy_tolerance_w,
            regions,
            max_face_pairs,
        )?;
        let approximate_contacts =
            self.evaluate(cx, &approximate_boundary, approximate_temperature)?;
        compare_goal_fields(
            cx,
            reference_problem,
            &approximate_boundary,
            approximate_temperature,
            full_nodal_weights,
            &endpoint,
            &approximate_contacts,
            || {
                endpoint.pullback_endpoint_with_mean_feedback(
                    cx,
                    &slopes[..m],
                    &reference_slopes[..m],
                    reference_feedback,
                    full_nodal_weights,
                    max_feedback_entries,
                )
            },
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn admit_goal_comparison(
        &self,
        cx: &Cx<'_>,
        problem: ConductionProblem<'_>,
        reference_temperature: &[f64],
        approximate_temperature: &[f64],
        linear: LinearConfig,
        residual_tolerance_w: f64,
        energy_tolerance_w: f64,
        full_nodal_weights: &[f64],
        max_face_pairs: usize,
    ) -> Result<(), ConductionError> {
        poll(cx, 0)?;
        self.check_mesh(problem.mesh)?;
        admit_linear(linear)?;
        if linear.restart == 0
            || !(residual_tolerance_w.is_finite() && residual_tolerance_w > 0.0)
            || !(energy_tolerance_w.is_finite() && energy_tolerance_w > 0.0)
        {
            return Err(invalid(
                "<steady-contact-goal>",
                "positive FGMRES restart and finite positive watt residual and heat targets are required",
            ));
        }
        let n = problem.mesh.vertex_count();
        for (field, values, temperature) in [
            (
                "steady contact reference temperature",
                reference_temperature,
                true,
            ),
            (
                "steady contact approximate temperature",
                approximate_temperature,
                true,
            ),
            ("steady contact goal weights", full_nodal_weights, false),
        ] {
            if values.len() != n {
                return Err(ConductionError::FieldLength {
                    field,
                    expected: n,
                    found: values.len(),
                });
            }
            for (i, &value) in values.iter().enumerate() {
                if i % 256 == 0 {
                    poll(cx, i)?;
                }
                if temperature {
                    positive(value)?;
                } else {
                    finite(value)?;
                }
            }
        }
        for &(vertex, prescribed) in problem.boundary.dirichlet() {
            poll(cx, vertex)?;
            if vertex >= n
                || reference_temperature[vertex].to_bits() != prescribed.to_bits()
                || approximate_temperature[vertex].to_bits() != prescribed.to_bits()
            {
                return Err(invalid(
                    "<steady-contact-goal>",
                    "both fields must retain the exact original prescribed temperatures",
                ));
            }
        }
        let count = self
            .variable_pairs
            .checked_add(self.nonmatching_triangles)
            .and_then(|base| {
                self.fixed.iter().try_fold(base, |count, surface| {
                    count.checked_add(surface.face_pairs().len())
                })
            })
            .ok_or_else(|| invalid("<steady-contact-goal>", "contact trace count overflows"))?;
        if count > max_face_pairs
            || count
                .checked_mul(432 * 24)
                .and_then(|bytes| {
                    n.checked_mul(16 * 8)
                        .and_then(|nodes| bytes.checked_add(nodes))
                })
                .is_none_or(|bytes| bytes > isize::MAX as usize)
        {
            return Err(invalid(
                "<steady-contact-goal>",
                "comparison exceeds its contact trace or addressable storage allowance",
            ));
        }
        Ok(())
    }
}

#[allow(clippy::too_many_arguments)]
fn compare_goal_fields(
    cx: &Cx<'_>,
    problem: ConductionProblem<'_>,
    approximate_boundary: &ThermalBoundary,
    approximate_temperature: &[f64],
    full_nodal_weights: &[f64],
    endpoint: &ContactSteadyLinearization<'_>,
    approximate_contacts: &TemperatureContactEvaluation,
    solve_dual: impl FnOnce() -> Result<RobinGradient, ConductionError>,
) -> Result<DiscreteGoalComparison, ConductionError> {
    let reference_temperature = endpoint.temperature();
    let n = problem.mesh.vertex_count();
    // Retain the endpoint's original source-law evaluation. The other
    // field gets its own fresh evaluation, never this frozen resistance.
    let dofs = DofMap::new(problem.boundary, n)?;
    let assemble = |boundary, temperature, interfaces| {
        assemble_operator_scaled_with_interfaces(
            cx,
            problem.mesh,
            boundary,
            problem.material,
            problem.source,
            temperature,
            None,
            Some(interfaces),
            problem.element_materials,
        )
    };
    let (reference_residual, primal_relative_residual) = {
        let system = assemble(
            problem.boundary,
            reference_temperature,
            endpoint.contact().physical_contacts().frozen_interfaces(),
        )?;
        let (matrix, rhs) = reduce(&system, &dofs);
        let relative = true_residual(&matrix, &dofs.gather(reference_temperature), &rhs)?;
        // The relative residual is diagnostic: physical admission above
        // uses the caller's watt tolerances, not the dual tolerance.
        (residual(&system, &dofs, reference_temperature), relative)
    };
    let approximate_residual = residual(
        &assemble(
            approximate_boundary,
            approximate_temperature,
            approximate_contacts.frozen_interfaces(),
        )?,
        &dofs,
        approximate_temperature,
    );
    let dual = solve_dual()?;
    let mut nodal_contributions = vec![0.0; n];
    let mut signed_residual_change = 0.0;
    for (index, &vertex) in dofs.free().iter().enumerate() {
        if index % 256 == 0 {
            poll(cx, index)?;
        }
        let difference = finite(reference_residual[index] - approximate_residual[index])?;
        let contribution = finite(dual.nodal_load[vertex] * difference)?;
        nodal_contributions[vertex] = contribution;
        signed_residual_change = finite(signed_residual_change + contribution)?;
    }
    let mut signed_goal_change = 0.0;
    for vertex in 0..n {
        if vertex % 256 == 0 {
            poll(cx, vertex)?;
        }
        let difference = finite(reference_temperature[vertex] - approximate_temperature[vertex])?;
        signed_goal_change =
            finite(signed_goal_change + finite(full_nodal_weights[vertex] * difference)?)?;
    }
    let linearization_remainder = finite(signed_goal_change - signed_residual_change)?;
    poll(cx, dual.iterations)?;
    Ok(DiscreteGoalComparison {
        signed_residual_change,
        linearization_remainder,
        signed_goal_change,
        nodal_contributions,
        primal_relative_residual,
        dual_relative_residual: dual.relative_residual,
        dual_iterations: dual.iterations,
        uses_nonlinear_jacobian: endpoint.uses_nonlinear_jacobian(),
    })
}
