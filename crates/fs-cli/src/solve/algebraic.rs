//! Algebraic error of the published region maximum. Only the complete fixed
//! linear solid model enters the outward analyzer; a frozen coupled or
//! nonlinear operator cannot certify the original equations.

use super::{
    BTreeMap, EvidenceWork, PropagatedTerm, QoiRegionTraceError, RungSolved, SolveRefusal,
    conduction_error, trace_qoi_region_vertices,
};

pub(super) fn maximum_term(
    cx: &fs_exec::Cx<'_>,
    solved: &RungSolved,
    region: &str,
    region_ids: &BTreeMap<String, u32>,
    memory_bytes: u64,
    work: EvidenceWork<'_>,
) -> Result<Option<PropagatedTerm>, SolveRefusal> {
    let cancelled = || {
        conduction_error(
            "cli-solve-cancelled",
            "the published-maximum algebraic analysis was cancelled",
            "rerun the solve",
        )
    };
    let gap = |reason: String| Ok(Some(PropagatedTerm::Unmeasured { reason }));
    let Some(data) = &solved.adjoint_data else {
        return gap("the published rung retained no final-state operator".to_string());
    };
    // The existing whole-model tolerance comparison remains available for
    // these models. Do not silently hold their feedback variables constant.
    if !data.air_paths.is_empty() || data.radiating_boundary.is_some() {
        return Ok(None);
    }
    for element in 0..solved.labels.len() {
        if element % 1024 == 0 {
            cx.checkpoint().map_err(|_| cancelled())?;
        }
        match data.materials.model_for(element) {
            Ok(model) if model.is_temperature_dependent() => return Ok(None),
            Ok(_) => {}
            Err(error) => return gap(format!("algebraic material binding refused: {error}")),
        }
    }
    let Some(&region_id) = region_ids.get(region) else {
        return gap("the temperature-maximum region has no element label".to_string());
    };
    let (vertices, _) = trace_qoi_region_vertices(
        &solved.labels,
        &solved.mesh.complex().tets,
        solved.mesh.vertex_count(),
        region_id,
        work,
    )
    .map_err(|error| match error {
        QoiRegionTraceError::Cancelled { .. } => cancelled(),
        QoiRegionTraceError::UnitCountOverflow => conduction_error(
            "cli-solve-conduction-algebraic",
            "the algebraic region trace exceeded its representable work count",
            "rerun with a smaller admitted mesh",
        ),
    })?;
    let problem = fs_conduction::ConductionProblem {
        mesh: &solved.mesh,
        boundary: &data.boundary,
        material: &data.fallback,
        element_materials: Some(&data.materials),
        source: &data.source,
    };
    // Explicit structural caps constrain the extra matrix/vector work. They
    // are not a claim about allocator peak memory, which the producer does
    // not measure. The existing Cx checkpoints bound cancellation latency.
    let config = fs_conduction::adjoint::LinearGoalAnalysisConfig {
        residual_limits: fs_solver::goal::GoalResidualLimits {
            max_rows: usize::try_from(memory_bytes / 256).unwrap_or(usize::MAX),
            max_nonzeros: usize::try_from(memory_bytes / 64).unwrap_or(usize::MAX),
        },
        max_stability_iterations: data.linear.max_iterations,
    };
    let analysis = match fs_conduction::adjoint::analyze_linear_maximum(
        cx,
        problem,
        data.interfaces.as_ref(),
        data.linear,
        &solved.solution.temperature,
        &vertices,
        config,
    ) {
        Ok(analysis) => analysis,
        Err(fs_conduction::ConductionError::Cancelled { .. }) => return Err(cancelled()),
        Err(error) => {
            return gap(format!(
                "published-maximum algebraic analysis refused: {error}"
            ));
        }
    };
    let linear_analysis = analysis.linear_analysis();
    let enclosure = &linear_analysis.enclosure;
    let Some(half_width_k) = analysis.algebraic_half_width_k() else {
        return gap(format!(
            "the stored linear operator has no finite verified maximum-error bound; \
             inverse status {:?}, primal residual upper {:e}; no tolerance \
             comparison is substituted for missing inverse evidence",
            enclosure.status(),
            enclosure.primal_residual_infinity_upper(),
        ));
    };
    Ok(Some(PropagatedTerm::Measured {
        half_width_k,
        method: "outward-linear-maximum-enclosure",
        detail: format!(
            "published mesh and temperature field; fixed linear conductivity, Robin \
             references and matching contact; {} region vertices ({} free), nominal \
             maximum {:e} K, enclosed maximum {:?} K; outward residual infinity \
             upper {:e}, verified inverse infinity upper {:?}, stability proposal \
             relative residual {:?}, stability iterations {}; \
             ||A^-1||_infinity * ||b-A*T||_infinity bounds the full region maximum, \
             including a change of hottest node and exact prescribed values; \
             stored floating-point system only, not assembly, discretization or physical \
             uncertainty; residual evaluation rounding is already included in this term; Estimated",
            vertices.len(),
            analysis.free_vertices(),
            analysis.nominal_k(),
            analysis.interval_k(),
            enclosure.primal_residual_infinity_upper(),
            enclosure.inverse_infinity_upper(),
            linear_analysis.stability_relative_residual,
            linear_analysis.stability_iterations,
        ),
        vertices: Vec::new(),
    }))
}
