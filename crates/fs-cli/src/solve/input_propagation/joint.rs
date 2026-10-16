//! Complete boundary/model-form vertex designs, including heated and cooled walls.
use super::{PropagatedTerm, SolveRefusal, account_joint_corner, conduction_error};

/// Two ambient ends crossed with two fan-pressure ends and two coefficient ends.
/// This bounds this additional design, not the complete solve invocation or RSS.
pub(in crate::solve) const MAX_CORNERS: usize = 8;
const MAX_BOUNDARY_POINTS: usize = 4;

fn bad(what: impl Into<String>) -> SolveRefusal {
    conduction_error("cli-solve-conduction-propagation", what,
        "retain the complete finite boundary/model-form corner design and its original physical inputs")
}

/// Cross EVERY isolated boundary point with BOTH coefficient extremes.
/// Stronger convection may increase heating or make a cold-side deviation larger;
/// selecting the hottest boundary and weaker coefficient is not a valid substitute.
///
/// Admission and deduplication precede physical work. Equal input pairs must have
/// equal retained results; a zero coefficient allowance reuses the boundary points.
/// The callback receives original input coordinates, not an updated previous corner.
/// It must propagate cancellation/work refusal as the outer Err. A physical refusal
/// is the inner Err and makes the whole interaction unknown, never a partial maximum.
///
/// Individual vertices stay in their original order. Every successful joint result
/// is retained in the term's detail, including the controlling largest ABSOLUTE
/// deviation. Only excess over the original individual widths is charged, once.
/// This is an Estimated vertex envelope under the existing response assumptions,
/// not a proof of interior monotonicity or a continuum/physical uncertainty bound.
#[allow(clippy::too_many_arguments)]
pub(in crate::solve) fn account_joint_design(
    boundary: PropagatedTerm,
    model_half_width_k: f64,
    nominal_k: f64,
    boundary_points: &[(f64, f64, f64)],
    coefficient_allowance: f64,
    max_corners: usize,
    mut evaluate: impl FnMut(String, f64, f64, f64)
        -> Result<Result<(String, f64), String>, SolveRefusal>,
) -> Result<(PropagatedTerm, f64), SolveRefusal> {
    let Some(boundary_width) = boundary.half_width() else {
        return Ok((boundary, f64::NAN));
    };
    if ![boundary_width, model_half_width_k].iter().all(|v| v.is_finite() && *v >= 0.0)
        || !(boundary_width + model_half_width_k).is_finite() || !nominal_k.is_finite()
        || !(coefficient_allowance.is_finite() && (0.0..1.0).contains(&coefficient_allowance))
        || boundary_points.is_empty() || boundary_points.len() > MAX_BOUNDARY_POINTS
    {
        return Err(bad("invalid individual evidence or boundary/model-form design"));
    }
    let mut points = Vec::new();
    points.try_reserve_exact(boundary_points.len())
        .map_err(|_| bad("boundary/model-form design allocation refused"))?;
    for &(temperature, factor, value) in boundary_points {
        if ![temperature, factor].iter().all(|v| v.is_finite() && *v > 0.0)
            || !value.is_finite() || !(value - nominal_k).is_finite()
        { return Err(bad("nonfinite or nonpositive boundary design coordinate")); }
        if let Some(&(_, _, retained)) = points.iter().find(|&&(t, f, _)| t == temperature && f == factor) {
            if retained != value {
                return Err(bad("identical boundary inputs carry inconsistent retained temperatures"));
            }
        } else {
            points.push((temperature, factor, value));
        }
    }
    let scales = [1.0 - coefficient_allowance, 1.0 + coefficient_allowance];
    if scales[0] == scales[1] {
        if model_half_width_k != 0.0 {
            return Err(bad("a zero coefficient allowance has a nonzero isolated model effect"));
        }
        // No model perturbation exists: these exact physical points already ran.
        // This is established zero interaction, unlike an unavailable model term.
        return Ok((boundary, 0.0));
    }
    let required = points.len().checked_mul(scales.len())
        .ok_or_else(|| bad("boundary/model-form corner count overflow"))?;
    let allowed = max_corners.min(MAX_CORNERS);
    if required > allowed {
        return account_joint_corner(boundary, model_half_width_k, nominal_k, Err(format!(
            "complete boundary/model-form design requires {required} physical corners, allowance {allowed}; no joint corner was evaluated"
        )));
    }
    let mut completed = Vec::new();
    completed.try_reserve_exact(required)
        .map_err(|_| bad("joint evidence allocation refused"))?;
    let mut worst = None::<(String, f64)>;
    for &(temperature, factor, _) in &points {
        for scale in scales {
            let label = format!("joint fluid {temperature} K, fan pressure x{factor}, card coefficient x{scale}");
            let result = evaluate(label.clone(), temperature, factor, scale)?;
            let (returned_label, value) = match result {
                Ok(row) => row,
                Err(reason) => return account_joint_corner(boundary, model_half_width_k, nominal_k,
                    Err(format!("{label}: {reason}; {}/{} joint corners completed; retained joint vertices: {}",
                        completed.len(), required, completed.join("; ")))),
            };
            if returned_label != label {
                return Err(bad("joint-corner result does not name the requested physical inputs"));
            }
            let deviation = (value - nominal_k).abs();
            if !value.is_finite() || !deviation.is_finite() {
                return Err(bad("nonfinite result in complete boundary/model-form design"));
            }
            completed.push(format!("{label} = {value} K"));
            if worst.as_ref().is_none_or(|(_, prior)| deviation > (prior - nominal_k).abs()) {
                worst = Some((label, value));
            }
        }
    }
    let worst = worst.ok_or_else(|| bad("complete boundary/model-form design produced no result"))?;
    let (mut term, excess) = account_joint_corner(boundary, model_half_width_k, nominal_k, Ok(worst))?;
    if let PropagatedTerm::Measured { detail, .. } = &mut term {
        *detail = format!("{detail}; full boundary/model-form design: {required}/{required} distinct joint corners; {}",
            completed.join("; "));
    }
    Ok((term, excess))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn measured(width: f64) -> PropagatedTerm {
        PropagatedTerm::Measured { half_width_k: width, method: "interval-vertex-resolve",
            detail: "isolated boundary evidence".into(),
            vertices: vec![("cold".into(), 290.0), ("hot".into(), 310.0)] }
    }
    fn points() -> [(f64, f64, f64); 2] { [(280.0, 1.0, 290.0), (320.0, 1.0, 310.0)] }

    #[test]
    fn stronger_exchange_can_control_both_cold_and_hot_deviations() {
        // Exact thermal balance: Q=100 W, fixed sink 290 K through 10 W/K,
        // fluid coefficient 10*s W/K. Nominal T=300 K for every s: isolated
        // model width is zero, but a changed ambient couples to that coefficient.
        let thermal = |ambient: f64, s: f64| (100.0 + 2900.0 + 10.0*s*ambient)/(10.0 + 10.0*s);
        let old_hottest_weak = (thermal(320.0, 0.5)-300.0).abs();
        assert!(old_hottest_weak < 10.0, "the old shortcut would charge zero");
        let mut calls = Vec::new();
        let (term, excess) = account_joint_design(measured(10.0), 0.0, 300.0,
            &points(), 0.5, 8, |label, t, f, s| {
                calls.push((t, f, s)); Ok(Ok((label, thermal(t, s))))
            }).unwrap();
        assert_eq!(calls, vec![(280.0,1.0,0.5), (280.0,1.0,1.5), (320.0,1.0,0.5), (320.0,1.0,1.5)]);
        assert_eq!(excess, 2.0);
        let PropagatedTerm::Measured { half_width_k, vertices, detail, .. } = term else { panic!("measured") };
        assert_eq!(half_width_k, 12.0);
        assert_eq!(vertices, vec![("cold".into(),290.0),("hot".into(),310.0)]);
        assert!(detail.contains("4/4 distinct joint corners"));
        assert!(detail.contains("card coefficient x1.5 = 288 K"));
        assert!(detail.contains("card coefficient x1.5 = 312 K"));
    }

    #[test]
    fn complete_cartesian_design_is_admitted_before_any_physical_work() {
        let rows = [(290.0,0.9,299.0),(290.0,1.1,298.0),(310.0,0.9,303.0),(310.0,1.1,302.0)];
        let mut calls = Vec::new();
        let (gap, excess) = account_joint_design(measured(3.0), 2.0, 300.0,
            &rows, 0.2, 7, |label,t,f,s| { calls.push((t,f,s)); Ok(Ok((label,300.0))) }).unwrap();
        assert!(calls.is_empty()); assert!(gap.half_width().is_none()); assert!(excess.is_nan());
        account_joint_design(measured(3.0), 2.0, 300.0, &rows, 0.2, 8,
            |label,t,f,s| { calls.push((t,f,s)); Ok(Ok((label,300.0))) }).unwrap();
        assert_eq!(calls.len(),8);
        for &(t,f,_) in &rows { for s in [0.8,1.2] { assert!(calls.contains(&(t,f,s))); } }
    }

    #[test]
    fn any_failed_corner_is_unknown_and_cancellation_is_not_swallowed() {
        let mut calls = 0;
        let (gap, excess) = account_joint_design(measured(10.0), 1.0, 300.0,
            &points(), 0.5, 8, |label,_,_,_| {
                calls += 1;
                if calls == 2 { Ok(Err("outside the material span".into())) }
                else { Ok(Ok((label,294.0))) }
            }).unwrap();
        assert_eq!(calls,2); assert!(excess.is_nan());
        let PropagatedTerm::Unmeasured { reason } = gap else { panic!("a subset is not a complete result") };
        assert!(reason.contains("outside the material span"));
        assert!(reason.contains("1/4 joint corners completed"));
        assert!(reason.contains("= 294 K"));
        for code in ["cli-solve-cancelled", "cli-solve-work-envelope"] {
            let error = account_joint_design(measured(10.0),1.0,300.0,&points(),0.5,8,
                |_,_,_,_| Err(conduction_error(code,"stop","retain prefix"))).err().unwrap();
            assert_eq!(error.code,code);
        }
    }

    #[test]
    fn duplicate_and_zero_width_coordinates_do_not_spend_duplicate_solves() {
        let points = [(300.0,1.0,300.0),(300.0,1.0,300.0)];
        let mut calls = 0;
        account_joint_design(measured(0.0),1.0,300.0,&points,0.2,2,
            |label,_,_,_| {calls+=1;Ok(Ok((label,300.0)))}).unwrap();
        assert_eq!(calls,2);
        let (term, excess) = account_joint_design(measured(0.0),0.0,300.0,&points,0.0,0,
            |_,_,_,_| {calls+=1;unreachable!("unchanged boundary result is already available")}).unwrap();
        assert_eq!(calls,2); assert_eq!(excess,0.0); assert_eq!(term.half_width(),Some(0.0));
        let inconsistent=[points[0],(300.0,1.0,301.0)];
        assert!(account_joint_design(measured(1.0),1.0,300.0,&inconsistent,0.2,8,
            |_,_,_,_| unreachable!("inconsistent evidence must refuse before solving")).is_err());
    }

    #[test]
    fn nonfinite_evidence_and_misaddressed_results_cannot_become_zero() {
        for bad_value in [f64::NAN,f64::INFINITY,f64::NEG_INFINITY] {
            assert!(account_joint_design(measured(10.0),1.0,300.0,&points(),0.5,8,
                |label,_,_,_| Ok(Ok((label,bad_value)))).is_err());
            assert!(account_joint_design(measured(10.0),bad_value,300.0,&points(),0.5,8,
                |_,_,_,_| unreachable!()).is_err());
        }
        for allowance in [-0.1,1.0,f64::NAN] {
            assert!(account_joint_design(measured(10.0),1.0,300.0,&points(),allowance,8,
                |_,_,_,_| unreachable!()).is_err());
        }
        assert!(account_joint_design(measured(10.0),1.0,300.0,&points(),0.5,8,
            |_,_,_,_| Ok(Ok(("a different physical corner".into(),300.0)))).is_err());
    }
}
