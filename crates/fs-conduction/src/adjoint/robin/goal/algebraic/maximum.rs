//! A maximum is 1-Lipschitz in the infinity norm, even when its active node
//! changes. Reuse the checked stored-system inverse rather than relabeling
//! the adjoint of the currently hottest node as a maximum-error certificate.

use super::{
    ConductionError, ConductionProblem, Cx, LinearConfig, LinearGoalAnalysis,
    LinearGoalAnalysisConfig, LinearGoalAnalyzer, ThermalInterfaces, bounded_solve,
    invalid, map_enclosure, poll,
};

/// A region maximum and its conditional, stored-linear-system error bound.
/// No geometry, assembly, continuum, material or coupled-physics bound is
/// implied. Prescribed nodes remain exact declared values.
#[derive(Debug, Clone, PartialEq)]
pub struct LinearMaximumAnalysis {
    nominal_k: f64,
    interval_k: Option<[f64; 2]>,
    algebraic_half_width_k: Option<f64>,
    free_vertices: usize,
    analysis: LinearGoalAnalysis,
}

impl LinearMaximumAnalysis {
    /// Maximum of the supplied finite nodal field on the selected region.
    #[must_use]
    pub const fn nominal_k(&self) -> f64 { self.nominal_k }

    /// Outward enclosure of the exact stored system's regional maximum.
    /// `None` means no finite bound was established, not a zero error.
    #[must_use]
    pub const fn interval_k(&self) -> Option<[f64; 2]> { self.interval_k }

    /// Absolute algebraic error allowance, including residual-evaluation
    /// roundoff. The selected maximum may move to any selected free node.
    #[must_use]
    pub const fn algebraic_half_width_k(&self) -> Option<f64> {
        self.algebraic_half_width_k
    }

    /// Number of selected free vertices. Input selections are set-valued:
    /// duplicates are refused, while their order has no numerical effect.
    #[must_use]
    pub const fn free_vertices(&self) -> usize { self.free_vertices }

    /// Underlying inverse/residual evidence and actual preparation work.
    /// Its linear goal is not the regional maximum; only the inverse and
    /// full primal residual are used for this maximum enclosure.
    #[must_use]
    pub const fn linear_analysis(&self) -> &LinearGoalAnalysis { &self.analysis }

    /// A missing bound or invalid tolerance cannot admit a stopping decision.
    #[must_use]
    pub fn meets_absolute_tolerance(&self, tolerance_k: f64) -> bool {
        tolerance_k.is_finite() && tolerance_k > 0.0
            && self.algebraic_half_width_k.is_some_and(|bound| bound <= tolerance_k)
    }
}

impl<'m> LinearGoalAnalyzer<'m> {
    /// Prepare an operator for one or many regional maxima, with no nodal
    /// adjoint solve. The declared stability budget may improve an already
    /// finite but loose unscaled inverse bound. Every candidate is checked;
    /// only a strictly smaller certified bound replaces the existing one.
    /// This matters when interior row sums are tiny positive rounding values:
    /// accepting the first finite inverse would produce an unusably wide band.
    ///
    /// # Errors
    /// The same field, assembly, budget and cancellation refusals as `new`.
    pub fn new_for_maximum(
        cx: &Cx<'_>, problem: ConductionProblem<'m>, interfaces: Option<&ThermalInterfaces>,
        linear: LinearConfig, temperature: &[f64], config: LinearGoalAnalysisConfig,
    ) -> Result<Self, ConductionError> {
        poll(cx, 0)?;
        let n = problem.mesh.vertex_count();
        let mut weights = Vec::new();
        weights.try_reserve_exact(n)
            .map_err(|_| invalid("regional maximum weights allocation refused"))?;
        for index in 0..n {
            if index % 512 == 0 { poll(cx, index)?; }
            weights.push(0.0);
        }
        let mut analyzer = Self::new(
            cx, problem, interfaces, linear, temperature, &weights, config,
        )?;
        // `new` already spends this budget when the unscaled inverse is
        // absent. Never repeat that attempt or exceed the declared total.
        if config.max_stability_iterations > 0 && analyzer.stability_relative_residual.is_none() {
            let previous = analyzer.analyze(cx, temperature)?;
            let rhs = vec![1.0; analyzer.dofs().n()];
            let (scaling, residual, iterations) = bounded_solve(
                cx, &analyzer.response.matrix, &rhs,
                LinearConfig { tolerance: 0.01,
                    max_iterations: config.max_stability_iterations, restart: linear.restart },
                false,
            )?;
            analyzer.stability_iterations = iterations;
            analyzer.stability_relative_residual = Some(residual);
            if scaling.iter().all(|value| value.is_finite() && *value > 0.0) {
                let free_temperature = analyzer.dofs().gather(temperature);
                let checked = fs_solver::goal::enclose_goal_error(
                    &analyzer.response.matrix, &analyzer.rhs, &free_temperature,
                    &analyzer.weights, &analyzer.free_dual, Some(&scaling),
                    config.residual_limits, || cx.checkpoint().is_ok(),
                ).map_err(map_enclosure)?;
                if let Some(bound) = checked.inverse_infinity_upper()
                    && previous.enclosure.inverse_infinity_upper().is_none_or(|old| bound < old)
                {
                    analyzer.stability_scaling = Some(scaling);
                }
            }
        }
        poll(cx, analyzer.stability_iterations)?;
        Ok(analyzer)
    }

    /// Assess a regional nodal maximum on this analyzer's immutable system.
    /// No primal, dual or stability solve is repeated. The inverse estimate
    /// is checked by the existing outward evaluator on this exact matrix.
    ///
    /// For r = b - A T, ||T* - T||_inf <= ||A^-1||_inf ||r||_inf.
    /// Therefore |max(T*) - max(T)| has the same bound without assuming a
    /// unique or unchanged hottest node. Fixed vertices tighten the interval
    /// and an entirely prescribed selection has zero algebraic error.
    ///
    /// # Errors
    /// Refuses empty, repeated or out-of-range selections, invalid fields,
    /// changed prescribed values, material extrapolation and cancellation.
    /// Unavailable inverse bounds and unrepresentable endpoints stay `None`.
    pub fn analyze_maximum(
        &self, cx: &Cx<'_>, temperature: &[f64], vertices: &[usize],
    ) -> Result<LinearMaximumAnalysis, ConductionError> {
        poll(cx, 0)?;
        let n = self.problem.mesh.vertex_count();
        if vertices.is_empty() || vertices.len() > n {
            return Err(invalid("a regional maximum needs a nonempty vertex set within the mesh"));
        }
        // Selection storage is bounded by the already admitted mesh, never by
        // unchecked vertex identifiers. Validate before any field indexing.
        let mut selected = Vec::new();
        selected.try_reserve_exact(n)
            .map_err(|_| invalid("regional maximum selection allocation refused"))?;
        for index in 0..n {
            if index % 512 == 0 { poll(cx, index)?; }
            selected.push(false);
        }
        for (index, &vertex) in vertices.iter().enumerate() {
            if index % 512 == 0 { poll(cx, index)?; }
            if vertex >= n || selected[vertex] {
                return Err(invalid("regional maximum vertices must be distinct and in range"));
            }
            selected[vertex] = true;
        }
        let analysis = self.analyze(cx, temperature)?;
        let mut nominal = f64::NEG_INFINITY;
        let mut fixed_max = f64::NEG_INFINITY;
        let mut free_max = f64::NEG_INFINITY;
        let mut free_vertices = 0;
        // Canonical mesh order makes even tie/signed-zero behavior independent
        // of the caller's set ordering.
        for (vertex, included) in selected.into_iter().enumerate() {
            if vertex % 512 == 0 { poll(cx, vertex)?; }
            if !included { continue; }
            let value = if temperature[vertex] == 0.0 { 0.0 } else { temperature[vertex] };
            nominal = nominal.max(value);
            if self.dofs().slot_of(vertex).is_some() {
                free_max = free_max.max(value);
                free_vertices += 1;
            } else {
                fixed_max = fixed_max.max(value);
            }
        }
        let half_width = if free_vertices == 0 {
            Some(0.0)
        } else {
            analysis.enclosure.inverse_infinity_upper().and_then(|inverse| {
                let residual = analysis.enclosure.primal_residual_infinity_upper();
                if residual == 0.0 { return Some(0.0); }
                let bound = fs_math::next_up(inverse * residual);
                bound.is_finite().then_some(bound)
            })
        };
        let interval = half_width.and_then(|bound| {
            if free_vertices == 0 || bound == 0.0 { return Some([nominal, nominal]); }
            let lower = fixed_max.max(fs_math::next_down(free_max - bound));
            let upper = fixed_max.max(fs_math::next_up(free_max + bound));
            (lower.is_finite() && upper.is_finite()).then_some([lower, upper])
        });
        poll(cx, 0)?;
        Ok(LinearMaximumAnalysis {
            nominal_k: nominal, interval_k: interval,
            algebraic_half_width_k: interval.and(half_width), free_vertices, analysis,
        })
    }
}

/// Prepare and assess one linear regional maximum without a primal re-solve.
/// The zero linear goal avoids an unnecessary hottest-node adjoint. A bounded
/// positive-scaling proposal is still checked by outward arithmetic before
/// it supplies any inverse authority. Constant heterogeneous materials,
/// prescribed/Neumann/fixed Robin boundaries and matching contact use the
/// same production assembly as `LinearGoalAnalyzer`.
///
/// # Errors
/// The analyzer's field, assembly, budget, nonlinear-material and cancellation
/// refusals, plus the selection refusals of `analyze_maximum`.
#[allow(clippy::too_many_arguments)]
pub fn analyze_linear_maximum(
    cx: &Cx<'_>, problem: ConductionProblem<'_>, interfaces: Option<&ThermalInterfaces>,
    linear: LinearConfig, temperature: &[f64], vertices: &[usize],
    config: LinearGoalAnalysisConfig,
) -> Result<LinearMaximumAnalysis, ConductionError> {
    poll(cx, 0)?;
    let n = problem.mesh.vertex_count();
    if vertices.is_empty() || vertices.len() > n {
        return Err(invalid("a regional maximum needs a nonempty in-range vertex set"));
    }
    for (index, &vertex) in vertices.iter().enumerate() {
        if index % 512 == 0 { poll(cx, index)?; }
        if vertex >= n { return Err(invalid("regional maximum vertex is out of range")); }
    }
    let analyzer = LinearGoalAnalyzer::new_for_maximum(
        cx, problem, interfaces, linear, temperature, config,
    )?;
    analyzer.analyze_maximum(cx, temperature, vertices)
}
