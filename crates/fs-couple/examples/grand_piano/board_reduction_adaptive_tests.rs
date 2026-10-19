use super::*;

#[test]
fn adaptive_bridge_accuracy_survives_tail_rotation_and_full_physical_damping() {
    use fs_math::c64::C64;
    use fs_la::eigen_complex::lu_complex;

    // Analytic diagonal physical pencil with unequal nodal masses. Source
    // vectors are exact mass-normalized eigenvectors, not a reduced-basis
    // fixture; force and motion are projected again through the FINAL phi.
    let frequency=[90.0,470.0,550.0,640.0,750.0,890.0,1050.0,1280.0];
    let mass=[1.0_f64,1.7,0.8,2.3,1.2,0.65,1.9,1.4];
    let n=frequency.len();
    let damping_ratio=0.15;
    let ports=vec![vec![0.7,1.0,-0.6,0.8,-0.4,0.3,-0.12,0.04],
        vec![1.04,1.31,-0.81,1.18,-0.59,0.51,-0.138,0.047]];
    let source=SliceReport {
        window:(0.0,(TAU*1400.0).powi(2)),below_low:0,below_high:n,expected:n,
        modes:frequency.iter().enumerate().map(|(i,hz)| {
            let lambda=(TAU*hz).powi(2);
            let mut phi=vec![0.0;n];phi[i]=1.0/mass[i].sqrt();
            ModePair {lambda,phi,residual:0.0,interval:(lambda,lambda)}
        }).collect(),
        stats:fs_modal::SliceStats {shift:0.0,factorizations:0,lanczos_iters:0,
            restarts:0,factor_nnz_l:0,factor_peak_bytes:0,pivots_delayed:0},
    };
    let options=ritz::RitzOptions {max_modes:7,keep_low_modes:1,
        sample_hz:vec![120.0],selection:ritz::Selection::Pivoted,transfer_tolerance:Some(0.002)};
    let reduced=project_modal(n,damping_ratio,&source,&ports,&options,
        |x,y|for i in 0..n {y[i]=mass[i]*x[i];}).unwrap();
    let r=reduced.modes.len();
    assert!(r>2 && r<n,"fixture must rotate a genuinely reduced tail");
    assert_eq!(reduced.modes[0].lambda,source.modes[0].lambda);
    assert_eq!(reduced.modes[0].phi,source.modes[0].phi);
    assert_eq!(reduced.modes[0].interval,source.modes[0].interval);
    let report=reduced.report.adaptive_response.as_ref().expect("adaptive report survives projection");
    assert_eq!(report.tolerance,options.transfer_tolerance.unwrap());
    assert_eq!(report.port_count,ports.len());
    assert!(report.enrichment_rounds>0);
    assert_eq!(report.evaluation_hz[0],0.0);
    for &hz in &frequency {
        assert!(report.evaluation_hz.iter().any(|point|(point-hz).abs()<1e-8));
    }
    let mut off_diagonal=0.0_f64;
    for i in 0..r {for j in 0..r {
        let mass_product:f64=(0..n).map(|a|
            mass[a]*reduced.modes[i].phi[a]*reduced.modes[j].phi[a]).sum();
        assert!((mass_product-if i==j {1.0}else{0.0}).abs()<1e-12);
        let expected_loss:f64=(0..n).map(|a|mass[a]*2.0*damping_ratio
            *source.modes[a].lambda.sqrt()*reduced.modes[i].phi[a]*reduced.modes[j].phi[a]).sum();
        assert!((reduced.physical_damping[i*r+j]-expected_loss).abs()
            <1e-10*expected_loss.abs().max(1.0),"final physical loss congruence ({i},{j})");
        if i!=j {off_diagonal=off_diagonal.max(reduced.physical_damping[i*r+j].abs());}
        assert_eq!(reduced.physical_damping[i*r+j],reduced.physical_damping[j*r+i]);
    }}
    assert!(off_diagonal>1e-5,"mixed-tail loss must retain off-diagonal terms");

    let mut maximum=0.0_f64;let mut maximum_transfer=0.0_f64;
    let mut diagonalized_error=0.0_f64;
    let mut damping_response_change=0.0_f64;
    let projected_ports:Vec<Vec<f64>>=ports.iter().map(|port| {
        let physical_force:Vec<_>=port.iter().zip(mass).map(|(g,m)|g*m.sqrt()).collect();
        reduced.modes.iter().map(|mode|mode.phi.iter()
            .zip(&physical_force).map(|(phi,f)|phi*f).sum::<f64>()).collect()
    }).collect();
    for (drive,port) in ports.iter().enumerate() {
        let force=&projected_ports[drive];
        for &hz in &report.evaluation_hz {
            let omega=TAU*hz;
            let exact=source.modes.iter().zip(port).fold(C64::ZERO,|sum,(mode,g)| {
                sum+C64::new(g*g,0.0)/C64::new(mode.lambda-omega*omega,
                    -omega*2.0*damping_ratio*mode.lambda.sqrt())
            });
            let mut full_compliance=C64::ZERO;
            for diagonal_only in [false,true] {
                let mut dynamic=vec![C64::ZERO;r*r];
                for i in 0..r {for j in 0..r {
                    let loss=if diagonal_only && i!=j {0.0}
                        else {reduced.physical_damping[i*r+j]};
                    dynamic[i*r+j]=C64::new(
                        if i==j {reduced.modes[i].lambda-omega*omega}else{0.0},
                        -omega*loss);
                }}
                let mut response:Vec<_>=force.iter().map(|f|C64::new(*f,0.0)).collect();
                lu_complex(&dynamic,r).unwrap().solve(&mut response);
                let compliance=force.iter().zip(&response)
                    .fold(C64::ZERO,|sum,(g,x)|sum+x.scale(*g));
                let error=(compliance-exact).abs()/exact.abs();
                if diagonal_only {
                    diagonalized_error=diagonalized_error.max(error);
                    damping_response_change=damping_response_change.max((compliance-full_compliance).abs()/exact.abs());
                } else {
                    maximum=maximum.max(error);full_compliance=compliance;
                    for (receive,receiver) in ports.iter().enumerate() {
                        let transfer=projected_ports[receive].iter().zip(&response)
                            .fold(C64::ZERO,|sum,(g,x)|sum+x.scale(*g));
                        let mut expected=C64::ZERO;let mut self_response=C64::ZERO;
                        for ((mode,a),b) in source.modes.iter().zip(port).zip(receiver) {
                            let divisor=C64::new(mode.lambda-omega*omega,
                                -omega*2.0*damping_ratio*mode.lambda.sqrt());
                            expected=expected+C64::from_re(a*b)/divisor;
                            self_response=self_response+C64::from_re(b*b)/divisor;
                        }
                        let scale=exact.abs().sqrt()*self_response.abs().sqrt();
                        maximum_transfer=maximum_transfer.max((transfer-expected).abs()/scale);
                    }
                }
            }
        }
    }
    eprintln!("adaptive final projection: source={n}, retained={r}, checked={}, rounds={}, self error={maximum:e}, complete transfer={maximum_transfer:e}, report={:e}, diagonal-only={diagonalized_error:e}, offdiag C={off_diagonal:e}, damping response change={damping_response_change:e}",
        report.evaluation_hz.len(),report.enrichment_rounds,report.max_normalized_transfer_error);
    assert!(damping_response_change>1e-5, "fixture must detect loss of physical cross-damping");
    assert!(maximum<=report.tolerance+1e-9,
        "final physical driving-point error {maximum} exceeds {}",report.tolerance);
    assert!((maximum-report.max_relative_driving_point_error).abs()<1e-8,
        "pre-rotation report must describe final K, C and physical ports");
    assert!(maximum_transfer<=report.tolerance+1e-9,
        "final physical transfer error {maximum_transfer} exceeds {}",report.tolerance);
    assert!((maximum_transfer-report.max_normalized_transfer_error).abs()<1e-8,
        "complete transfer report must describe the final nodal basis and full damping");
}

