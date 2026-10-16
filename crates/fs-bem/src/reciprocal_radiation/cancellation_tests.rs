use super::*;

fn budget()->Budget {
    Budget {maximum_unknowns:64,maximum_cache_bytes:1024*1024,
        maximum_work:20_000_000_000,angular_refinement:0}
}
fn medium()->Medium {Medium {density:1.0,sound_speed:1.0}}
fn bits(v:&[C64])->Vec<(u64,u64)> {
    v.iter().map(|c|(c.re.to_bits(),c.im.to_bits())).collect()
}
fn same(a:&[Solution<'_,'_>],b:&[Solution<'_,'_>]) {
    assert_eq!(a.len(),b.len());
    for (a,b) in a.iter().zip(b) {
        assert_eq!(bits(&a.pressure),bits(&b.pressure));
        assert_eq!(bits(&a.velocity),bits(&b.velocity));
        assert_eq!(bits(&a.coefficients),bits(&b.coefficients));
        assert_eq!(bits(&[a.surface_power]),bits(&[b.surface_power]));
        assert_eq!([a.outgoing_power_w,a.flux_balance,a.backward_residual,a.panels_per_wavelength].map(f64::to_bits),
            [b.outgoing_power_w,b.flux_balance,b.backward_residual,b.panels_per_wavelength].map(f64::to_bits));
    }
}

#[test]
fn full_reciprocal_batches_cancel_inside_lu_solves_and_flux_without_poisoning_geometry() {
    let surface=SpherePanels::icosphere(1.0,0).unwrap();
    let geometry=Prepared::new(&surface,1.0,budget(),||false).unwrap();
    let fields=[vec![C64::ONE;surface.areas().len()],
        surface.normals().iter().map(|n|C64::from_re(n[2])).collect()];
    let refs:Vec<_>=fields.iter().map(Vec::as_slice).collect();
    let before:Vec<_>=fields.iter().map(|v|bits(v)).collect();
    let mut assembly_polls=0;
    let (_,waves)=geometry.assemble(0.5,&mut || {assembly_polls+=1;false}).unwrap();
    let n=geometry.unknowns();
    let before_lu=1+assembly_polls+n; // entry, assembly, congruence rows
    let factor_polls=2+n*(n+1)/2; // entry, pivots/elimination rows, publication
    let after_lu=before_lu+factor_polls;
    let per_source=1+(2+2*n)+n+waves.len(); // entry, substitutions, residual, flux
    let mut total=0;
    let baseline=geometry.solve_batch(0.5,medium(),&refs,|| {total+=1;false}).unwrap();
    assert_eq!(total,after_lu+1+fields.len()*per_source+1);
    assert!(baseline[0].outgoing_power_w>0.0);
    for stop in [1,1+assembly_polls/2,before_lu+3,after_lu,
        after_lu+4,after_lu+1+1+(2+2*n)+1,
        after_lu+1+per_source+1,total-waves.len()/2,total] {
        let mut calls=0;
        let result=geometry.solve_batch(0.5,medium(),&refs,|| {calls+=1;calls==stop});
        assert!(matches!(result,Err(HelmholtzError::Cancelled)),"checkpoint {stop}");
        assert_eq!(calls,stop);
        assert_eq!(fields.iter().map(|v|bits(v)).collect::<Vec<_>>(),before);
        let retry=geometry.solve_batch(0.5,medium(),&refs,||false).unwrap();
        same(&baseline,&retry);
    }
}

#[test]
fn mechanical_velocity_ports_propagate_mid_factor_and_late_projection_cancellation() {
    let surface=SpherePanels::icosphere(1.0,0).unwrap();
    let geometry=Prepared::new(&surface,1.0,budget(),||false).unwrap();
    let rows=[vec![1.0;surface.areas().len()],vec![-0.4;surface.areas().len()]];
    let refs:Vec<_>=rows.iter().map(Vec::as_slice).collect();
    let mut assembly_polls=0;
    geometry.assemble(0.5,&mut || {assembly_polls+=1;false}).unwrap();
    let before_lu=1+rows.len()+1+assembly_polls+geometry.unknowns();
    let mut total=0;
    let baseline=geometry.solve_velocity_ports(0.5,medium(),&refs,|| {total+=1;false}).unwrap();
    for stop in [before_lu+3,total-1,total] {
        let mut calls=0;
        let result=geometry.solve_velocity_ports(0.5,medium(),&refs,|| {calls+=1;calls==stop});
        assert!(matches!(result,Err(HelmholtzError::Cancelled)));
        assert_eq!(calls,stop);
        let retry=geometry.solve_velocity_ports(0.5,medium(),&refs,||false).unwrap();
        assert_eq!(bits(baseline.impedance()),bits(retry.impedance()));
        same(baseline.solutions(),retry.solutions());
    }
}

#[test]
fn lu_cancellation_is_not_misreported_as_an_acoustically_singular_operator() {
    assert_eq!(lu_error(LuComplexError::Cancelled),HelmholtzError::Cancelled);
    assert_eq!(lu_error(LuComplexError::Singular{index:3}),HelmholtzError::Singular);
}
