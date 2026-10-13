//! An endpoint failure is printed, not promoted to a usable manufactured shape.
use super::*;
use fs_topopt::sdf3::response::limits::binary::BinaryScanOptions3;

pub(super) fn run(study: &mut CutDensityStudy3<AdaptiveSolveSpace3>, source: &ResponseLimitEvaluation3,
    cases: &[ResponseLimitCase3<'_>], options: MinimumMaterialOptions3, control: &mut SolveControl<'_>)
    -> Result<(), Box<dyn std::error::Error>> {
    // Explicit finite schedule, not a proof of optimality over all binary masks.
    // The API accepts other schedules, caps and floors; this example retains
    // the source material floor and demands zero normalized limit relaxation.
    let thresholds = [0.2, 0.35, 0.5, 0.65, 0.8];
    let policy = BinaryScanOptions3 { void_stiffness: study.params().e_min,
        limits: options.limits, feasibility_tolerance: 0.0, ..Default::default() };
    let scan = study.scan_binary_response_limits(&source.rho, cases, &thresholds, policy, control)?;
    println!("binary_candidate,threshold,solid_volume_fraction,smooth_constraint,worst_original_residual,response_feasible,volume_feasible");
    for (i, row) in scan.candidates.iter().enumerate() {
        println!("{i},{:.17e},{:.17e},{:.17e},{:.17e},{},{}", row.threshold, row.volume_fraction,
            row.constraint, row.worst_residual, row.response_feasible, row.volume_feasible);
    }
    // Retained endpoint evidence is emitted even if a subsequent trial failed.
    if let Some(best) = &scan.best {
        println!("binary_cell_level,cell_i,cell_j,cell_k,source_raw,source_projected,solid,stiffness_scale");
        for (i, leaf) in study.operator().elasticity().leaves().iter().enumerate() {
            let [x,y,z] = leaf.index();
            let stiffness = if best.solid[i] { 1.0 } else { best.void_stiffness };
            println!("{},{x},{y},{z},{:.17e},{:.17e},{},{stiffness:.17e}", leaf.level(),
                source.rho[i], source.projected_rho[i], best.solid[i]);
        }
        println!("binary_case,limit,sense,value,bound,scale,original_residual");
        for (i, case) in cases.iter().enumerate() { for (j, limit) in case.limits.iter().enumerate() {
            println!("{i},{j},{:?},{:.17e},{:.17e},{:.17e},{:.17e}", limit.sense,
                best.measurements[i][j], limit.bound, limit.scale, best.residuals[i][j]);
        } }
        println!("binary_case,master_node,x,y,z,ux,uy,uz");
        for (case, field) in best.displacements.iter().enumerate() {
            for (node, (p,u)) in study.operator().elasticity().nodes().iter().zip(field.chunks_exact(3)).enumerate() {
                println!("{case},{node},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e}",
                    p[0],p[1],p[2],u[0],u[1],u[2]);
            }
        }
        eprintln!("binary_best_index={:?}; threshold={:.17e}; gray_volume={:.17e}; solid_volume={:.17e}",
            scan.best_index,best.threshold,best.source_volume_fraction,best.volume_fraction);
    }
    eprintln!("binary_completed={}/{}; void_stiffness={:.17e}; volume_cap={}; response_tolerance={}; has_feasible_endpoint={}; work={:?}; source_material_preserved=true; removed_cell_geometry=false; continuum_certified=false; manufacturability_certified=false",
        scan.candidates.len(),thresholds.len(),policy.void_stiffness,policy.volume_cap,
        policy.feasibility_tolerance,scan.best.is_some(),scan.work);
    scan.outcome?;
    if scan.best.is_none() {
        return Err("no tested binary threshold met every physical limit and the volume cap; gray feasibility does not transfer".into());
    }
    Ok(())
}
