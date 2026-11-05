use super::*;
use core::f64::consts::TAU;
fn air() -> Medium { Medium { density:1.2, sound_speed:343.0 } }
fn phase(x:f64) -> C64 { C64::new(x.cos(),x.sin()) }
fn sample(t:[[f64;3];3],shape:Vec<f64>) -> AreaSample {
    let a=[t[1][0]-t[0][0],t[1][1]-t[0][1]];
    let b=[t[2][0]-t[0][0],t[2][1]-t[0][1]];
    AreaSample {position_m:[t[0][0]+a[0]/3.+b[0]/3.,t[0][1]+a[1]/3.+b[1]/3.],
        area_m2:0.5*(a[0]*b[1]-a[1]*b[0]).abs(),mode_shape:shape}
}
fn options(source:&PreparedCoherentRadiation<'_>,receivers:usize) -> TriangleReceiverOptions {
    TriangleReceiverOptions { minimum_clearance_m:1e-6,
        maximum_row_bytes:source.triangle_receiver_row_bytes(receivers).unwrap(),
        integration:near_field::Options {relative_tolerance:1e-7,maximum_depth:14,maximum_kernel_evaluations:2_000_000}}
}
#[test]
fn g1_finite_cells_resolve_the_piston_near_field_and_reactive_velocity() {
    // The receiver is 5 mm above a 40 mm piston, much closer than a coarse
    // cell's diameter. Its finite axial velocity cannot be recovered by
    // replacing the triangles with monopoles at their centroids.
    let (radius,z,omega)=(0.04_f64,0.005_f64,TAU*700.);
    let v=C64::new(0.7,-0.3);
    let k=omega/air().sound_speed;
    let edge=(radius*radius+z*z).sqrt();
    let expected_p=(phase(k*z)-phase(k*edge))*v.scale(air().density*air().sound_speed);
    let expected_u=(phase(k*z)-phase(k*edge).scale(z/edge))*v;
    let expected_g=expected_u*C64::new(0.,air().density*omega);
    let mut previous=[f64::INFINITY;2];
    for segments in [16,32,64] {
        let triangles:Vec<_>=(0..segments).map(|i| {
            let a=TAU*i as f64/segments as f64;
            let b=TAU*(i+1) as f64/segments as f64;
            [[0.,0.,0.],[radius*a.cos(),radius*a.sin(),0.],[radius*b.cos(),radius*b.sin(),0.]]
        }).collect();
        let cells:Vec<_>=triangles.iter().map(|&t|sample(t,vec![1.])).collect();
        let surface=PreparedSurface::new(&cells).unwrap();
        let source=surface.prepare_coherent_radiation(omega,air(),&[v],||false).unwrap();
        let rows=source.prepare_triangle_receivers(&triangles,&[[0.,0.,z]],options(&source,1),||false).unwrap();
        let field=rows.evaluate(&source,||false).unwrap()[0];
        let error=[(field.pressure_pa-expected_p).abs()/expected_p.abs(),
            (field.particle_velocity_m_s[2]-expected_u).abs()/expected_u.abs()];
        for axis in 0..2 { assert!(error[axis]<0.3*previous[axis],"{segments}: {error:?} after {previous:?}"); }
        assert!((field.pressure_gradient_pa_m[2]-expected_g).abs()/expected_g.abs()<1.01*error[1]+1e-12);
        assert!(field.particle_velocity_m_s[0].abs()+field.particle_velocity_m_s[1].abs()<1e-10);
        assert!((rows.clearances_m()[0]-z).abs()<1e-14);
        assert!(rows.kernel_evaluations()>80*segments,"near-cell quadrature must actually refine");
        assert!(field.pressure_quadrature_error_estimate_pa.is_finite() && field.pressure_quadrature_error_estimate_pa>=0.);
        previous=error;
    }
    assert!(previous[0]<0.002 && previous[1]<0.001,"finest finite-piston errors {previous:?}");
}
#[test]
fn g1_pressure_null_retains_its_gradient_and_reuses_the_same_rows() {
    let right=[[0.02,-0.02,0.],[0.06,-0.02,0.],[0.04,0.02,0.]];
    let left=right.map(|p|[-p[0],p[1],p[2]]);
    let triangles=[right,left];
    let cells=[sample(right,vec![1.,0.]),sample(left,vec![0.,1.])];
    let surface=PreparedSurface::new(&cells).unwrap();
    let omega=TAU*500.;
    let v=C64::new(0.7,0.3);
    let one=surface.prepare_coherent_radiation(omega,air(),&[v,C64::ZERO],||false).unwrap();
    let both=surface.prepare_coherent_radiation(omega,air(),&[v,-v],||false).unwrap();
    let rows=one.prepare_triangle_receivers(&triangles,&[[0.,0.,0.025]],options(&one,1),||false).unwrap();
    let baseline=rows.evaluate(&one,||false).unwrap()[0];
    let null=rows.evaluate(&both,||false).unwrap()[0];
    assert!(null.pressure_pa.abs()<1e-13*baseline.pressure_pa.abs());
    assert!(null.pressure_gradient_pa_m[0].abs()>1.);
    assert!((null.pressure_gradient_pa_m[0]-baseline.pressure_gradient_pa_m[0].scale(2.)).abs()
        <1e-13*baseline.pressure_gradient_pa_m[0].abs());
    for axis in 0..3 {
        assert!((null.pressure_gradient_pa_m[axis]-null.particle_velocity_m_s[axis]*C64::new(0.,air().density*omega)).abs()<1e-10);
    }
    assert!(null.particle_velocity_m_s[1].abs()+null.particle_velocity_m_s[2].abs()<1e-12);
    let shifted=surface.prepare_coherent_radiation(omega,air(),&[v*phase(0.7),-v*phase(0.7)],||false).unwrap();
    let changed=rows.evaluate(&shifted,||false).unwrap()[0];
    assert!((changed.particle_velocity_m_s[0]-null.particle_velocity_m_s[0]*phase(0.7)).abs()<1e-13);
}
#[test]
fn g1_finite_receivers_converge_to_the_existing_outgoing_image() {
    let triangles=[
        [[-0.04,-0.03,0.],[0.04,-0.03,0.],[0.04,0.03,0.]],
        [[-0.04,-0.03,0.],[0.04,0.03,0.],[-0.04,0.03,0.]],
    ];
    let cells=triangles.map(|t|sample(t,vec![1.]));
    let surface=PreparedSurface::new(&cells).unwrap();
    let omega=TAU*500.;
    let source=surface.prepare_coherent_radiation(omega,air(),&[C64::new(0.7,-0.3)],||false).unwrap();
    let points=[[0.,0.,0.8],[0.,0.,3.2],[0.,0.,12.8]];
    let rows=source.prepare_triangle_receivers(&triangles,&points,options(&source,3),||false).unwrap();
    assert_eq!(rows.points_m(),&points);
    let fields=rows.evaluate(&source,||false).unwrap();
    let far=source.far_field(&[[0.,0.,1.]],2,||false).unwrap()[0];
    let mut previous=f64::INFINITY;
    for (point,field) in points.iter().zip(fields) {
        let r=point[2];
        let value=field.pressure_pa.scale(r)*phase(-omega/air().sound_speed*r);
        let error=(value-far).abs()/far.abs();
        assert!(error<0.3*previous,"r={r}: {error} after {previous}");
        previous=error;
    }
    assert!(previous<0.001);
}
#[test]
fn g2_finite_cell_receiver_keeps_the_original_dark_modal_residual() {
    let triangle=[[-0.02,-0.02,0.],[0.02,-0.02,0.],[0.,0.02,0.]];
    let epsilon=2.0_f64.powi(-36);
    let surface=PreparedSurface::new(&[sample(triangle,vec![1.,1.+epsilon])]).unwrap();
    let v=C64::new(0.7,-0.3);
    let one=surface.prepare_coherent_radiation(TAU*700.,air(),&[v,C64::ZERO],||false).unwrap();
    let dark=surface.prepare_coherent_radiation(TAU*700.,air(),&[v,-v],||false).unwrap();
    let rows=one.prepare_triangle_receivers(&[triangle],&[[0.001,0.003,0.01]],options(&one,1),||false).unwrap();
    let a=rows.evaluate(&one,||false).unwrap()[0];
    let b=rows.evaluate(&dark,||false).unwrap()[0];
    for (value,reference) in [b.pressure_pa].into_iter().chain(b.pressure_gradient_pa_m).chain(b.particle_velocity_m_s)
        .zip([a.pressure_pa].into_iter().chain(a.pressure_gradient_pa_m).chain(a.particle_velocity_m_s)) {
        let expected=reference.scale(-epsilon);
        assert!((value-expected).abs()<2e-14*expected.abs());
    }
}
#[test]
fn g0_finite_cells_refuse_geometry_clearance_budget_and_source_mismatch() {
    let triangle=[[-0.02,-0.02,0.],[0.02,-0.02,0.],[0.,0.02,0.]];
    let surface=PreparedSurface::new(&[sample(triangle,vec![1.])]).unwrap();
    let source=surface.prepare_coherent_radiation(TAU*700.,air(),&[C64::ONE],||false).unwrap();
    let points=[[0.,0.,0.01]];
    let settings=options(&source,1);
    let rows=source.prepare_triangle_receivers(&[triangle],&points,settings,||false).unwrap();
    let wrong_geometry=triangle.map(|p|[p[0]+0.001,p[1],p[2]]);
    assert!(source.prepare_triangle_receivers(&[wrong_geometry],&points,settings,||false).is_err());
    let mut wrong_area=triangle;wrong_area[1][0]+=0.001;
    assert!(source.prepare_triangle_receivers(&[wrong_area],&points,settings,||false).is_err());
    assert!(source.prepare_triangle_receivers(&[triangle],&[[0.,0.,0.]],settings,||false).is_err());
    assert!(source.prepare_triangle_receivers(&[triangle],&[[0.,0.,-0.01]],settings,||false).is_err());
    assert!(source.prepare_triangle_receivers(&[triangle],&points,
        TriangleReceiverOptions {maximum_row_bytes:settings.maximum_row_bytes-1,..settings},||false).is_err());
    assert!(source.prepare_triangle_receivers(&[triangle],&[points[0];2],
        TriangleReceiverOptions {maximum_row_bytes:source.triangle_receiver_row_bytes(2).unwrap(),
            integration:near_field::Options {maximum_kernel_evaluations:80,..settings.integration},..settings},||false).is_err());
    assert!(source.triangle_receiver_row_bytes(65).is_err());
    let other=PreparedSurface::new(&[sample(triangle,vec![1.])]).unwrap();
    let other=other.prepare_coherent_radiation(TAU*700.,air(),&[C64::ONE],||false).unwrap();
    assert!(rows.evaluate(&other,||false).is_err());
    let frequency=surface.prepare_coherent_radiation(TAU*701.,air(),&[C64::ONE],||false).unwrap();
    assert!(rows.evaluate(&frequency,||false).is_err());
    assert!(matches!(source.prepare_triangle_receivers(&[triangle],&points,settings,||true),Err(HelmholtzError::Cancelled)));
    assert!(matches!(rows.evaluate(&source,||true),Err(HelmholtzError::Cancelled)));
    let empty=source.prepare_triangle_receivers(&[triangle],&[],options(&source,0),||false).unwrap();
    assert_eq!(empty.kernel_evaluations(),0);
    assert!(empty.evaluate(&source,||false).unwrap().is_empty());
}
