//! Bounded summaries of the native pointwise affine stress envelopes.
use super::*;
use fs_topopt::sdf3::robust::envelope::affine_load_space::stress::{
    AffineStressEvaluation3, AffineStressLimit3, AffineStressWitness3,
};

pub(super) const SCOPE: &str = "Strength is the unweighted p-norm of native pointwise relaxed and physical von Mises envelopes on the eroded material. Every sampled point maximizes its own quadratic stress over the complete affine load ball, including common self-weight. Different points can have different adverse coefficients; the exported worst_displacement is the compliance witness. This conservative sampled aggregate has numerical oracle diagnostics and supplies no continuum or floating-point enclosure.";

pub(super) fn policy(limit: AffineStressLimit3, joint_sharpness: f64) -> String {
    format!(
        "{{\"kind\":\"affine-sampled-peak-bound\",\"stress_limit_pa\":{:.17e},\"relaxation_power\":{:.17e},\"aggregation_power\":{:.17e},\"max_stress_point_fields\":{},\"max_fields\":{},\"joint_sharpness\":{:.17e}}}",
        limit.limit,
        limit.options.relaxation_power,
        limit.options.aggregation_power,
        limit.options.max_points,
        limit.options.max_cases,
        joint_sharpness,
    )
}

fn witness(w: &AffineStressWitness3) -> String {
    let o = &w.oracle;
    format!(
        "{{\"cell\":{},\"position_m\":{:?},\"quadrature_weight_m3\":{:.17e},\"load_coefficients\":{:?},\"stress_layout\":\"xx, yy, zz, xy, yz, zx\",\"reference_stress_pa\":{:?},\"reference_von_mises_pa\":{:.17e},\"relaxed_von_mises_pa\":{:.17e},\"physical_von_mises_pa\":{:.17e},\"reference_envelope_pa\":{:.17e},\"relaxed_envelope_pa\":{:.17e},\"physical_envelope_pa\":{:.17e},\"relative_defect\":{:.17e},\"stress_normalization_pa\":{:.17e},\"oracle\":{{\"branch\":{:?},\"gradient_admissible\":{},\"normalized_value_without_offset\":{:.17e},\"normalized_dual_value_without_offset\":{:.17e},\"normalized_duality_gap\":{:.17e},\"relative_duality_gap\":{:.17e},\"relative_separation\":{:.17e},\"stationarity_residual\":{:.17e},\"feasibility_residual\":{:.17e},\"complementarity_residual\":{:.17e},\"eigenpair_residual\":{:.17e},\"iterations\":{}}}}}",
        w.cell,
        w.position,
        w.weight,
        w.coefficients(),
        w.reference_stress,
        w.reference_von_mises,
        w.relaxed_von_mises,
        w.physical_von_mises,
        w.worst_reference_von_mises,
        w.worst_relaxed_von_mises,
        w.worst_physical_von_mises,
        w.relative_defect,
        w.normalization,
        format!("{:?}", o.branch),
        o.gradient_admissible(),
        o.value,
        o.dual_value,
        o.duality_gap,
        o.relative_duality_gap,
        o.relative_separation,
        o.stationarity_residual,
        o.feasibility_residual,
        o.complementarity_residual,
        o.eigenpair_residual,
        o.iterations,
    )
}

fn governing(
    e: &AffineStressEvaluation3,
    measure: impl Fn(&AffineStressWitness3) -> f64,
) -> Result<&AffineStressWitness3> {
    let first = e.witnesses.first().ok_or_else(|| {
        fail(
            "cli-study-sdf3-strength-output",
            "retained strength has no sampled witness",
        )
    })?;
    // Keep the first retained point when several attain exactly the same value.
    Ok(e.witnesses.iter().skip(1).fold(first, |best, current| {
        if measure(current) > measure(best) {
            current
        } else {
            best
        }
    }))
}

pub(super) fn measurement(e: &AffineStressEvaluation3) -> Result<String> {
    let physical = governing(e, |w| w.worst_physical_von_mises)?;
    let relaxed = governing(e, |w| w.worst_relaxed_von_mises)?;
    let maximum = |measure: fn(&AffineStressWitness3) -> f64| {
        e.witnesses.iter().map(measure).fold(0.0_f64, f64::max)
    };
    let nonsmooth = e
        .witnesses
        .iter()
        .filter(|w| {
            (w.worst_reference_von_mises != 0.0 || w.reference_von_mises != 0.0)
                && !w.oracle.gradient_admissible()
        })
        .count();
    Ok(format!(
        "{{\"stress_limit_pa\":{:.17e},\"aggregate_pa\":{:.17e},\"constraint\":{:.17e},\"feasible\":{},\"gradient_admissible\":{},\"sampled_relaxed_max_pa\":{:.17e},\"sampled_physical_max_pa\":{:.17e},\"relaxed_envelope_max_pa\":{:.17e},\"physical_envelope_max_pa\":{:.17e},\"point_count\":{},\"point_fields\":{},\"adjoint_count\":{},\"worst_physical\":{},\"worst_relaxed\":{},\"oracle_diagnostics\":{{\"maximum_relative_defect\":{:.17e},\"maximum_relative_duality_gap\":{:.17e},\"maximum_stationarity_residual\":{:.17e},\"maximum_feasibility_residual\":{:.17e},\"maximum_complementarity_residual\":{:.17e},\"maximum_eigenpair_residual\":{:.17e},\"nonsmooth_positive_points\":{nonsmooth}}}}}",
        e.limit.limit,
        e.aggregate,
        e.constraint,
        e.physically_feasible(),
        e.gradient_admissible(),
        e.sampled_relaxed_max,
        e.sampled_physical_max,
        relaxed.worst_relaxed_von_mises,
        physical.worst_physical_von_mises,
        e.point_count,
        e.point_fields,
        e.adjoints.len(),
        witness(physical),
        witness(relaxed),
        maximum(|w| w.relative_defect),
        maximum(|w| w.oracle.relative_duality_gap),
        maximum(|w| w.oracle.stationarity_residual),
        maximum(|w| w.oracle.feasibility_residual),
        maximum(|w| w.oracle.complementarity_residual),
        maximum(|w| w.oracle.eigenpair_residual),
    ))
}

pub(super) fn fields(e: &AffineStressEvaluation3, cells: usize) -> Result<String> {
    if e.witnesses.len() != e.point_count || e.witnesses.is_empty() {
        return Err(fail(
            "cli-study-sdf3-strength-output",
            "incomplete retained stress samples",
        ));
    }
    let mut physical = vec![0.0_f64; cells];
    let mut relaxed = vec![0.0_f64; cells];
    let mut counts = vec![0_usize; cells];
    for w in &e.witnesses {
        if w.cell >= cells
            || [w.worst_physical_von_mises, w.worst_relaxed_von_mises]
                .iter()
                .any(|v| !v.is_finite() || *v < 0.0)
        {
            return Err(fail(
                "cli-study-sdf3-strength-output",
                "invalid retained cell stress envelope",
            ));
        }
        physical[w.cell] = physical[w.cell].max(w.worst_physical_von_mises);
        relaxed[w.cell] = relaxed[w.cell].max(w.worst_relaxed_von_mises);
        counts[w.cell] = counts[w.cell].checked_add(1).ok_or_else(|| {
            fail(
                "cli-study-sdf3-strength-output",
                "stress sample count overflow",
            )
        })?;
    }
    Ok(format!(
        "{{\"schema\":\"sdf3-affine-strength-fields.v1\",\"cell_physical_envelope_pa\":{physical:?},\"cell_relaxed_envelope_pa\":{relaxed:?},\"cell_sample_count\":{counts:?}}}"
    ))
}

pub(super) fn html(
    selected: &AffineStressEvaluation3,
    accepted: &AffineStressEvaluation3,
) -> String {
    format!(
        "<h2>Strength over the affine load family</h2><p>The declared stress limit is {:.8e} Pa. The selected design has sampled strength aggregate {:.8e} Pa (cap-feasible: {}); the last accepted design has aggregate {:.8e} Pa (cap-feasible: {}). The complete joint constraint also includes compliance and its declared conservative aggregation.</p><p>{SCOPE}</p>",
        selected.limit.limit,
        selected.aggregate,
        selected.physically_feasible(),
        accepted.aggregate,
        accepted.physically_feasible(),
    )
}
