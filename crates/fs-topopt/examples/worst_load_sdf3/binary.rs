//! Convert the gray cap-constrained result into independently re-solved endpoints.
use super::*;
use fs_topopt::sdf3::robust::envelope::limits::{ComplianceCaps3, ComplianceLimitedOptions3,
    binary::{BinaryComplianceScanOptions3, BinaryComplianceScanStop3}};

#[allow(clippy::too_many_arguments)]
pub fn finish(study: &mut CutDensityStudy3, forces: &[&[f64]], rho: &[f64], caps: &ComplianceCaps3,
    design: ComplianceLimitedOptions3, keys: &[[usize; 3]], nodes: &[[f64; 3]], control: &mut SolveControl<'_>)
    -> Result<(), Box<dyn std::error::Error>> {
    let thresholds = [0.2, 0.35, 0.5, 0.65, 0.8];
    // Separate, explicit physical endpoint policy; never infer zero-void
    // behavior from the gray SIMP regularization parameter.
    let options = BinaryComplianceScanOptions3 { robustness: design.robustness,
        void_stiffness: 1e-3, max_load_cases: 2, max_candidates: 5,
        dilated_volume_limit: 1.0, smooth_tolerance: design.optimizer.tolerance };
    eprintln!("mode=binary-cap-endpoints; thresholds={thresholds:?}; endpoint_policy={options:?}; caps={:?}; continuum_certificate=false; geometry_extracted=false; gray_feasibility_is_binary_feasibility=false",
        caps.caps());
    let report = study.scan_binary_compliance_caps(rho, forces, caps, &thresholds, options, control)?;
    println!("candidate,threshold,eroded_volume,nominal_volume,dilated_volume,envelope,maximum_ratio,all_physical_caps,feasible");
    for c in &report.candidates {
        println!("{},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e},{},{}", c.index, c.threshold,
            c.eroded_volume, c.nominal_volume, c.dilated_volume, c.envelope, c.maximum_ratio,
            c.physical_caps_satisfied, c.feasible);
    }
    eprintln!("binary_scan_stop={:?}; binary_error={:?}; selected_index={:?}; cumulative_solve_work={:?}",
        report.stop, report.error, report.selected_index, report.work);
    let Some(best) = &report.best else {
        return Err("no completed binary endpoint met every cap and material limit; gray output is not a binary fallback".into());
    };
    eprintln!("exported_design=least-dilated-material-binary-endpoint; threshold={}; void_stiffness={}; independently_solved_binary_eroded_fields=true; zero_void_stiffness_claimed=false; manufactured_shape=false",
        best.options.threshold, best.options.void_stiffness);
    println!("load_case,binary_eroded_compliance,physical_cap,cap_satisfied");
    for (i, (&c, &cap)) in best.case_compliances.iter().zip(caps.caps()).enumerate() {
        println!("{i},{c:.17e},{cap:.17e},{}", c <= cap);
    }
    let scales = best.eroded_stiffness_scales();
    println!("cell_x,cell_y,cell_z,eroded_solid,nominal_solid,dilated_solid,eroded_stiffness");
    for (i, k) in keys.iter().enumerate() {
        println!("{},{},{},{},{},{},{:.17e}", k[0], k[1], k[2], u8::from(best.eroded.solid[i]),
            u8::from(best.nominal.solid[i]), u8::from(best.dilated.solid[i]), scales[i]);
    }
    println!("node_x,node_y,node_z,binary_uy_x,binary_uy_y,binary_uy_z,binary_uz_x,binary_uz_y,binary_uz_z");
    for (i, p) in nodes.iter().enumerate() {
        let uy = &best.displacements[0]; let uz = &best.displacements[1];
        println!("{:.17e},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e}",
            p[0], p[1], p[2], uy[3*i], uy[3*i+1], uy[3*i+2], uz[3*i], uz[3*i+1], uz[3*i+2]);
    }
    // A late stop retains a useful binary prefix, but is not a completed scan.
    if report.stop != BinaryComplianceScanStop3::ScheduleComplete {
        return Err("binary endpoint search stopped early; exported binary result is only the retained feasible prefix".into());
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    #[test]
    fn malformed_binary_mode_refuses_before_geometry_work() {
        assert!(super::super::limits::run_binary(&[]).is_err());
        let invalid: Vec<String> = ["NaN", "1500"].iter().map(|s| s.to_string()).collect();
        assert!(super::super::limits::run_binary(&invalid).is_err());
        let invalid: Vec<String> = ["1000", "1500", "not-an-iteration-count"].iter().map(|s| s.to_string()).collect();
        assert!(super::super::limits::run_binary(&invalid).is_err());
    }

    #[test]
    fn binary_mode_runs_the_restored_gray_and_binary_physics_path() {
        // Loose authored caps isolate routing/physical integration, not an
        // optimization-quality claim. Zero updates still runs the seed solve,
        // gradient gate, fresh recheck and all actual binary load families.
        let args: Vec<String> = ["100000000", "100000000", "0", "250000"]
            .iter().map(|s| s.to_string()).collect();
        super::super::limits::run_binary(&args).unwrap();
    }
}
