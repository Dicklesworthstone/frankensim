//! Binary selection keeps the authored continuous load set and physical budget.
use super::*;
use fs_topopt::sdf3::robust::envelope::load_space::design::limits::binary::{
    BinaryLoadSpaceScanOptions3, BinaryLoadSpaceScanStop3,
};

#[allow(clippy::too_many_arguments)]
pub fn finish(study: &mut CutDensityStudy3, basis: &[&[f64]], rho: &[f64], cap: &LoadSpaceCap3,
    design: LoadSpaceLimitedOptions3, volume_limit: f64, keys: &[[usize;3]], nodes: &[[f64;3]],
    control: &mut SolveControl<'_>) -> Result<(), Box<dyn std::error::Error>> {
    // Binary stiffness is explicit, not inferred from gray SIMP regularization.
    // Enumerate every distinct three-field mask state rather than relying on a
    // coarse hand-authored threshold list which can skip feasible endpoints.
    let options = BinaryLoadSpaceScanOptions3 { robustness: design.robustness,
        void_stiffness: 1e-3, max_load_cases: 2, max_candidates: 64,
        dilated_volume_limit: volume_limit, smooth_tolerance: design.optimizer.tolerance };
    eprintln!("mode=binary-simultaneous-load; authored_cap={}; threshold_policy=all-distinct-three-field-states; endpoint_policy={options:?}; gray_feasibility_is_binary_feasibility=false; separate_columns_authorize_pass=false",
        cap.cap());
    let report = study.scan_all_binary_load_space_cap(rho, basis, cap, options, control)?;
    println!("candidate,threshold,eroded_volume,nominal_volume,dilated_volume,worst_compliance,physical_cap,spectral_maximum,smooth_maximum,load_set_feasible,feasible");
    for c in &report.candidates {
        println!("{},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e},{},{}",
            c.index, c.threshold, c.eroded_volume, c.nominal_volume, c.dilated_volume,
            c.worst_compliance, cap.cap(), c.maximum, c.smooth_maximum, c.load_set_feasible, c.feasible);
    }
    eprintln!("binary_scan_stop={:?}; binary_error={:?}; selected_index={:?}; cumulative_physical_work={:?}",
        report.stop, report.error, report.selected_index, report.work);
    let Some(best) = &report.best else {
        return Err("no completed binary candidate met the full load-set cap and material limit; gray output is not a binary fallback".into());
    };
    eprintln!("exported_source=least-dilated-binary-load-space-endpoint; threshold={}; void_stiffness={}; combined_load_feasible={}; witness_relative_defect={:.3e}; continuum_certificate=false; zero_void_stiffness=false; geometry_extracted=false",
        best.physical.options.threshold, best.physical.options.void_stiffness,
        best.feasible(options.smooth_tolerance), best.witness_relative_defect);
    println!("basis_index,binary_maximizing_coefficient,binary_basis_compliance_diagnostic");
    for (i, (&q, &c)) in best.worst_coefficients.iter().zip(&best.physical.case_compliances).enumerate() {
        println!("{i},{q:.17e},{c:.17e}");
    }
    let scales = best.physical.eroded_stiffness_scales();
    println!("cell_x,cell_y,cell_z,eroded_solid,nominal_solid,dilated_solid,binary_eroded_stiffness");
    for (i,k) in keys.iter().enumerate() {
        println!("{},{},{},{},{},{},{:.17e}", k[0],k[1],k[2], u8::from(best.physical.eroded.solid[i]),
            u8::from(best.physical.nominal.solid[i]), u8::from(best.physical.dilated.solid[i]), scales[i]);
    }
    println!("node_x,node_y,node_z,binary_worst_fx,binary_worst_fy,binary_worst_fz,binary_worst_ux,binary_worst_uy,binary_worst_uz");
    for (i,p) in nodes.iter().enumerate() {
        let f = &best.worst_force; let u = &best.worst_displacement;
        println!("{:.17e},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e}",
            p[0],p[1],p[2],f[3*i],f[3*i+1],f[3*i+2],u[3*i],u[3*i+1],u[3*i+2]);
    }
    if report.stop != BinaryLoadSpaceScanStop3::ScheduleComplete {
        return Err("binary load-space selection stopped; exported endpoint is the retained completed prefix".into());
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn binary_load_space_command_checks_inputs_and_runs_full_physics() {
        for args in [vec![], vec!["NaN"], vec!["1000","bad-count"],
            vec!["1000","0","250000","NaN"], vec!["1000","0","250000","1.1"],
            vec!["1000","0","250000","1","extra"]] {
            let args: Vec<String> = args.into_iter().map(String::from).collect();
            assert!(super::super::run_binary(&args).is_err());
        }
        // Loose authored cap isolates command wiring, not optimization quality.
        // Zero updates still restores, audits, rechecks and solves all candidates.
        let args = ["--binary-material","100000000","0","250000"].map(String::from);
        super::super::super::run(&args).unwrap();
    }

    #[test]
    fn no_binary_material_admission_is_an_error_not_a_gray_success() {
        let mut poll = |_| ControlFlow::Continue(());
        let mut q = QuadratureControl3::new(QuadratureOptions3 { depth: 1, ..Default::default() }, &mut poll).unwrap();
        let op = CutElasticity3::build(HexCell::try_new([0.0;3],[1.0;3]).unwrap(), [2;3], &Slab,
            &IsotropicElastic::new(1.0,0.3,1.0).unwrap(), &|p| p[0]==0.0,
            ElasticityOptions3::default(), &mut q).unwrap();
        let y = op.body_load(&|_| [0.0,-1.0,0.0], || ControlFlow::Continue(())).unwrap();
        let z = op.body_load(&|_| [0.0,0.0,-1.0], || ControlFlow::Continue(())).unwrap();
        let keys = op.cell_keys(); let nodes = op.nodes().to_vec();
        let mut study = CutDensityStudy3::new(op, 0.15, SimpParams::default());
        let raw = vec![1.0; study.cells()];
        let cap = LoadSpaceCap3::new(2,1e9,0.02).unwrap();
        let mut poll = |_| ControlFlow::Continue(());
        let mut control = SolveControl::new(SolveBudget::default(), &mut poll);
        // Every distinct threshold preserves full raw-one material; none fits
        // zero volume. The command must not reuse the feasible gray result instead.
        assert!(finish(&mut study, &[&y,&z], &raw, &cap, LoadSpaceLimitedOptions3::default(),
            0.0, &keys, &nodes, &mut control).is_err());
        assert!(control.work().linear_solves >= 3);
    }
}
