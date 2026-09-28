//! Keep algebraic error small compared with the actual adaptive comparison.
//! Every changed field is prolonged and probed again before it can mark cells
//! or supply the published discretization estimate.

use super::{
    ADAPTIVE_ORDER_ONE_FACTOR, AdaptiveProbe, EvidenceWork, QoiRegionTraceError, RungSolved,
    SolveRefusal, adaptive_deadline, adaptive_probe, adaptive_prolongation, algebraic,
    canonical_f64, conduction_error, json_string, trace_qoi_region_vertices,
};

const MAX_CORRECTION_ROUNDS: usize = 4;

pub(super) struct BalancedPair {
    pub(super) probe: AdaptiveProbe,
    pub(super) receipt: String,
    /// Unsupported models retain their existing adaptive comparison. A covered
    /// pair with an unresolved algebraic budget cannot publish discretization.
    pub(super) may_use_comparison: bool,
}

fn number(value: Option<f64>) -> Result<String, SolveRefusal> {
    value.map_or_else(
        || Ok("null".to_string()),
        |value| {
            canonical_f64(value).ok_or_else(|| {
                conduction_error(
                    "cli-solve-conduction-adaptive-algebraic",
                    "adaptive algebraic-budget arithmetic is nonfinite",
                    "inspect the declared material and temperature scales",
                )
            })
        },
    )
}

fn sum_up(a: f64, b: f64) -> Option<f64> {
    let sum = if a == 0.0 {
        b
    } else if b == 0.0 {
        a
    } else {
        fs_math::next_up(a + b)
    };
    sum.is_finite().then_some(sum)
}

fn vertices(
    solved: &RungSolved,
    region: u32,
    work: EvidenceWork<'_>,
) -> Result<Vec<usize>, SolveRefusal> {
    trace_qoi_region_vertices(
        &solved.labels,
        &solved.mesh.complex().tets,
        solved.mesh.vertex_count(),
        region,
        work,
    )
    .map(|(vertices, _)| vertices)
    .map_err(|error| {
        conduction_error(
            if matches!(error, QoiRegionTraceError::Cancelled { .. }) {
                "cli-solve-cancelled"
            } else {
                "cli-solve-conduction-adaptive-algebraic"
            },
            "adaptive algebraic region tracing refused",
            "rerun with an admitted mesh and execution budget",
        )
    })
}

#[allow(clippy::too_many_arguments)]
pub(super) fn probe_pair(
    cx: &fs_exec::Cx<'_>,
    complex: &fs_mesh::LabeledTetComplex,
    split: &fs_mesh::TetRefinement,
    coarse: &mut RungSolved,
    fine: &mut RungSolved,
    region: u32,
    memory_bytes: u64,
    work: EvidenceWork<'_>,
    deadline: Option<(std::time::Instant, f64)>,
) -> Result<BalancedPair, SolveRefusal> {
    let probe_current = |coarse: &RungSolved, fine: &RungSolved| {
        adaptive_deadline(deadline)?;
        let approximate = adaptive_prolongation(cx, complex, coarse, split, fine)?;
        let probe = adaptive_probe(cx, coarse, fine, &approximate, region)?;
        adaptive_deadline(deadline)?;
        Ok::<_, SolveRefusal>(probe)
    };
    let covered = coarse.algebraic.linear_work.is_some() && fine.algebraic.linear_work.is_some();
    let mut probe = probe_current(coarse, fine)?;
    let mut rounds = 0;
    let mut field_updates = 0;
    let initial_coarse_iterations = coarse.algebraic.primal_iterations;
    let initial_fine_iterations = fine.algebraic.primal_iterations;
    let mut selected = None;
    let (status, discretization, allowance, coarse_bound, fine_bound, pair_bound) = loop {
        let discretization = ADAPTIVE_ORDER_ONE_FACTOR
            * probe
                .estimated_change_k
                .abs()
                .max(probe.measured_change_k.abs());
        // Division by the exact integer ten and downward rounding avoid
        // enlarging the allowance through binary64 evaluation of 0.1 * D.
        let allowance = fs_math::next_down(discretization / 10.0).max(0.0);
        let a = algebraic::maximum_bound(&coarse.algebraic);
        let b = algebraic::maximum_bound(&fine.algebraic);
        let pair = a.zip(b).and_then(|(a, b)| sum_up(a, b));
        let finish = |status| (status, discretization, allowance, a, b, pair);
        if !covered {
            break finish("unsupported-model");
        }
        if !discretization.is_finite() || !allowance.is_finite() {
            return Err(conduction_error(
                "cli-solve-conduction-adaptive-algebraic",
                "the observed discretization allowance is not finite",
                "inspect the declared temperature and material scales",
            ));
        }
        if pair.is_some_and(|bound| bound <= allowance) {
            break finish("balanced");
        }
        let Some((a, b)) = a.zip(b) else {
            break finish("bound-unavailable");
        };
        if rounds == MAX_CORRECTION_ROUNDS {
            break finish("correction-round-budget");
        }
        let target = fs_math::next_down(allowance / 2.0).max(0.0);
        if target == 0.0 {
            break finish("unallocated-discretization-scale");
        }
        if selected.is_none() {
            selected = Some((
                vertices(coarse, region, work)?,
                vertices(fine, region, work)?,
            ));
        }
        let (coarse_vertices, fine_vertices) = selected.as_ref().expect("selected pair region");
        adaptive_deadline(deadline)?;
        let changed_coarse = a > target
            && algebraic::retarget_linear_maximum(
                cx,
                coarse,
                coarse_vertices,
                target,
                memory_bytes,
            )?;
        adaptive_deadline(deadline)?;
        let changed_fine = b > target
            && algebraic::retarget_linear_maximum(cx, fine, fine_vertices, target, memory_bytes)?;
        rounds += 1;
        field_updates += usize::from(changed_coarse) + usize::from(changed_fine);
        adaptive_deadline(deadline)?;
        if !changed_coarse && !changed_fine {
            let a = algebraic::maximum_bound(&coarse.algebraic);
            let b = algebraic::maximum_bound(&fine.algebraic);
            let pair = a.zip(b).and_then(|(a, b)| sum_up(a, b));
            let status = if pair.is_some_and(|bound| bound <= allowance) {
                "balanced"
            } else {
                "correction-unresolved"
            };
            break (status, discretization, allowance, a, b, pair);
        }
        probe = probe_current(coarse, fine)?;
    };
    if covered {
        algebraic::finish_discretization_balance(coarse, discretization, allowance)?;
        algebraic::finish_discretization_balance(fine, discretization, allowance)?;
    }
    let additional_iterations = coarse
        .algebraic
        .primal_iterations
        .checked_sub(initial_coarse_iterations)
        .zip(
            fine.algebraic
                .primal_iterations
                .checked_sub(initial_fine_iterations),
        )
        .and_then(|(coarse, fine)| coarse.checked_add(fine))
        .ok_or_else(|| {
            conduction_error(
                "cli-solve-conduction-adaptive-algebraic",
                "adaptive pair correction work is outside its representable cumulative count",
                "report the inconsistent correction counters",
            )
        })?;
    let receipt = format!(
        "{{\"status\":{},\"correction_rounds\":{},\"correction_round_limit\":{},\"additional_correction_iterations\":{},\"field_updates\":{},\
         \"discretization_half_width_k\":{},\"pair_allowance_k\":{},\
         \"coarse_bound_k\":{},\"fine_bound_k\":{},\"pair_bound_k\":{},\
         \"coarse_correction_iterations\":{},\"fine_correction_iterations\":{},\
         \"rule\":\"coarse-plus-fine algebraic bound <= one tenth of the current observed discretization estimate\"}}",
        json_string(status),
        rounds,
        MAX_CORRECTION_ROUNDS,
        additional_iterations,
        field_updates,
        number(Some(discretization))?,
        number(Some(allowance))?,
        number(coarse_bound)?,
        number(fine_bound)?,
        number(pair_bound)?,
        coarse.algebraic.primal_iterations,
        fine.algebraic.primal_iterations,
    );
    Ok(BalancedPair {
        probe,
        receipt,
        may_use_comparison: !covered || status == "balanced",
    })
}
